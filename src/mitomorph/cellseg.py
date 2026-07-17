"""Cell-segmentation compute (Module A): nuclei + Cellpose bodies + territories.

CELL segmentation only — extracted out of the old ``cellqc`` stage. Per sample:
nuclei are Otsu/watershed-segmented, Cellpose instance bodies are detected from the
mito channel, and the two are reconciled into nucleus-seeded territories (one raster
label per cell, cell_id == nucleus_id == territory label).

This stage does NOT threshold mito: it READS the mito mask produced upstream by the
``mitoseg`` stage (``03_mitoseg/{group}/{filestem}/masks/mito.npy``) — the territory
watershed needs the mito signal to price ridge costs across mito-empty gaps.

Two flag-only tables are written (rows are NEVER dropped): a per-cell ``cellseg``
table, and its ``cellqc`` projection — the frozen Module A->B contract.
"""
from __future__ import annotations

import numpy as np

from .config import ExperimentConfig
from .io import read_manifest, resolve_samples, to_uint8, write_csv
from .segment import build_territories, cell_bodies, resolve_instances, segment_nuclei

# Per-cell cell-segmentation table (frozen columns): one row per territory.
CELLSEG_FIELDS = [
    "experiment", "group", "replicate", "sample", "px_um", "cell_id", "nucleus_id",
    "nucleus_area_um2", "territory_area_um2", "edge_distance_px",
    "border_touch", "footprint_confidence",
]

# Module A -> B contract: one row per cell, keyed on the territory raster label.
# nucleus_id == cell_id (each territory is labelled by its seeding nucleus).
# footprint_confidence: fraction of the territory covered by the Cellpose footprint
# (1.0 = fully body-backed, 0.0 = radius-disk fallback; NaN if Cellpose unavailable).
QC_FIELDS = [
    "experiment", "group", "replicate", "sample",
    "cell_id", "nucleus_id", "area_um2", "border_touch", "footprint_confidence",
]


def _border_cells(labels: np.ndarray) -> set[int]:
    """Cell ids whose BODY reaches the frame edge (true border_touch)."""
    edge = set(labels[0, :]) | set(labels[-1, :]) | set(labels[:, 0]) | set(labels[:, -1])
    return {int(c) for c in edge if c > 0}


def _cell_rows(
    cfg, group, replicate, filestem, px_um, nucleus_labels, territories, footprint,
) -> list[dict]:
    """One row per territory with the frozen cellseg columns (flag-only)."""
    rows = []
    px_area_um2 = px_um**2
    h, w = territories.shape
    border_ids = _border_cells(territories)

    for cell_id in range(1, int(nucleus_labels.max()) + 1):
        territory = territories == cell_id
        if not territory.any():
            continue
        # Crop to territory bbox: areas / edge distance are local to the cell.
        ys, xs = np.nonzero(territory)
        y0, y1 = int(ys.min()), int(ys.max()) + 1
        x0, x1 = int(xs.min()), int(xs.max()) + 1
        terr_c = territory[y0:y1, x0:x1]
        nuc_c = nucleus_labels[y0:y1, x0:x1] == cell_id
        if not nuc_c.any():
            continue

        territory_area_um2 = float(terr_c.sum()) * px_area_um2
        min_edge_distance = min(y0, x0, h - y1, w - x1)
        terr_px = int(terr_c.sum())
        if footprint is not None:
            fp_cover = int((footprint[y0:y1, x0:x1] & terr_c).sum())
            footprint_confidence = fp_cover / terr_px if terr_px else 0.0
        else:
            footprint_confidence = float("nan")

        rows.append({
            "experiment": cfg.experiment, "group": group, "replicate": replicate,
            "sample": filestem, "px_um": px_um, "cell_id": cell_id,
            "nucleus_id": cell_id,
            "nucleus_area_um2": float(nuc_c.sum()) * px_area_um2,
            "territory_area_um2": territory_area_um2,
            "edge_distance_px": min_edge_distance,
            "border_touch": cell_id in border_ids,
            "footprint_confidence": footprint_confidence,
        })
    return rows


def _qc_row(r: dict) -> dict:
    """Project a cellseg row onto the Module A->B QC contract (QC_FIELDS)."""
    return {
        "experiment": r["experiment"], "group": r["group"],
        "replicate": r["replicate"], "sample": r["sample"],
        "cell_id": r["cell_id"], "nucleus_id": r["nucleus_id"],
        "area_um2": r["territory_area_um2"],
        "border_touch": r["border_touch"],
        "footprint_confidence": r["footprint_confidence"],
    }


def compute_one(cfg: ExperimentConfig, filestem: str, manifest_rows: dict) -> None:
    """Run cell segmentation (Module A) for one sample."""
    manifest_row = manifest_rows[(cfg.experiment, f"{filestem}.czi")]
    px_um = float(manifest_row["px_x_um"])
    if int(manifest_row["size_z"]) == 1:
        print(f"WARNING: skipping {cfg.experiment}/{filestem} (size_z == 1; cannot MIP)")
        return

    group, replicate = cfg.parse_sample(filestem)
    czi = cfg.raw / f"{filestem}.czi"
    mip_path = cfg.results / "02_extract" / group / filestem / "figures" / f"{filestem}_mito_mip.tif"
    if not mip_path.exists():
        print(f"ERROR: missing mito MIP: {mip_path}. Run `mito extract` first.")
        return
    if not czi.exists():
        print(f"ERROR: missing raw CZI: {czi}")
        return

    mito_mask_path = cfg.results / "03_mitoseg" / group / filestem / "masks" / "mito.npy"
    if not mito_mask_path.exists():
        raise ValueError(f"missing mito mask: {mito_mask_path}. Run `mito mitoseg` first.")

    import tifffile

    from .io import channel_mip

    mito_mip = tifffile.imread(mip_path)
    nucleus_mip = to_uint8(channel_mip(czi, cfg.nucleus_channel_index)[0])
    mito_mask = np.load(mito_mask_path).astype(bool)
    nucleus_labels, _ = segment_nuclei(nucleus_mip, cfg, px_um)

    # Module A: trust Cellpose instance bodies, reconcile with nucleus seeds
    # (1-nucleus bodies kept intact; multi-nucleus bodies split). Fall back to the
    # legacy Voronoi only if Cellpose is unavailable in this env.
    try:
        bodies = cell_bodies(mito_mip, nucleus_mip, cfg, px_um)
        territories = resolve_instances(nucleus_labels, bodies, mito_mask, cfg, px_um)
        footprint = bodies > 0  # for footprint_confidence (fraction body-backed)
    except ImportError:
        print(f"WARNING: cellpose unavailable; legacy Voronoi territories for {filestem}")
        footprint = None
        territories = build_territories(nucleus_labels, cfg, px_um)

    sample_dir = cfg.results / "04_cellseg" / group / filestem
    (sample_dir / "masks").mkdir(parents=True, exist_ok=True)
    np.save(sample_dir / "masks" / "nuclei.npy", nucleus_labels)
    np.save(sample_dir / "masks" / "territories.npy", territories)

    rows = _cell_rows(
        cfg, group, replicate, filestem, px_um, nucleus_labels, territories, footprint,
    )

    table_dir = sample_dir / "tables"
    write_csv(table_dir / f"{filestem}_cellseg.csv", rows, CELLSEG_FIELDS)
    write_csv(table_dir / f"{filestem}_cellqc.csv", [_qc_row(r) for r in rows], QC_FIELDS)
    print(f"cellseg {cfg.experiment}/{filestem}: {len(rows)} cells", flush=True)


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Batch cell segmentation over requested samples (or all raw CZI)."""
    manifest_rows = read_manifest(cfg)
    for filestem, _czi in resolve_samples(cfg, samples):
        if filestem in cfg.skip_samples:
            print(f"WARNING: skipping {cfg.experiment}/{filestem} (config skip_samples)")
            continue
        try:
            compute_one(cfg, filestem, manifest_rows)
        except (KeyError, ValueError) as exc:
            print(f"ERROR: {cfg.experiment}/{filestem}: {exc}")
