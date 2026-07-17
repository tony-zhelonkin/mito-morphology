"""Cell-QC compute: nucleus-seeded territories + per-cell mitochondrial metrics.

QC is objective, condition-blind, and flag-only (rows are never dropped):
  - border_touch: the true cell BODY reaches the frame edge (not an inflated bbox)
  - out_of_focus: focus_score below a per-image relative threshold (DAPI-based)
  - empty:        mito area below mito_min_um2

Territories come from the footprint-bounded, nucleus-seeded watershed
(``build_territories_v2``); merges/over-splits are resolved upstream, so there is
no "exactly one nucleus" post-hoc filtering here.

Focus is flagged relative to each image's own median, never by a cross-experiment
constant, because focus_score is not comparable across imaging batches.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from skimage import morphology, segmentation

from .config import ExperimentConfig
from .io import read_manifest, resolve_samples, to_uint8, write_csv
from .metrics import compute_skeleton_metrics, focus_score, mito_shape_metrics
from .segment import (
    build_territories,
    cell_bodies,
    pipeline_li,
    resolve_instances,
    segment_nuclei,
)

CELL_FIELDS = [
    "experiment", "group", "replicate", "sample", "px_um", "cell_id", "nucleus_id",
    "nucleus_area_um2", "territory_area_um2", "mito_area_um2", "mito_area_fraction",
    "skeleton_length_um", "branch_count", "mito_solidity", "mito_form_factor",
    "mito_border_frac", "focus_score", "edge_distance_px", "footprint_confidence",
    "border_touch", "out_of_focus", "empty", "qc_pass",
]

# Module A -> B contract: one row per cell, keyed on the territory raster label.
# nucleus_id == cell_id (each territory is labelled by its seeding nucleus).
# footprint_confidence: fraction of the territory covered by the Cellpose footprint
# (1.0 = fully body-backed, 0.0 = radius-disk fallback; NaN if Cellpose unavailable).
QC_FIELDS = [
    "experiment", "group", "replicate", "sample",
    "cell_id", "nucleus_id", "area_um2", "border_touch", "footprint_confidence",
]

METRIC_FIELDS = [
    "nucleus_area_um2", "territory_area_um2", "mito_area_um2", "mito_area_fraction",
    "skeleton_length_um", "branch_count", "mito_solidity", "mito_form_factor", "focus_score",
]

SUMMARY_FIELDS = [
    "experiment", "group", "replicate", "sample", "px_um", "li_threshold",
    "focus_threshold", "cell_N", "n_qc_pass", "n_objective_pass",
    "n_border", "n_out_of_focus", "n_empty",
] + [f"median_{f}" for f in METRIC_FIELDS] + [f"objective_median_{f}" for f in METRIC_FIELDS]

SKIPPED_FIELDS = ["experiment", "sample", "px_um", "size_z", "skipped", "skip_reason"]


def _border_cells(labels: np.ndarray) -> set[int]:
    """Cell ids whose BODY reaches the frame edge (true border_touch)."""
    edge = set(labels[0, :]) | set(labels[-1, :]) | set(labels[:, 0]) | set(labels[:, -1])
    return {int(c) for c in edge if c > 0}


def _cell_metrics(
    cfg, group, replicate, filestem, px_um,
    nucleus_mip, mito_mask, nucleus_labels, territories, footprint,
) -> list[dict]:
    """One metrics row per territory (focus flag applied later, per-image)."""
    rows = []
    lap = ndimage.laplace(nucleus_mip.astype(float))
    territory_boundaries = segmentation.find_boundaries(territories, mode="inner") & (
        territories > 0
    )
    px_area_um2 = px_um**2
    h, w = territories.shape
    border_ids = _border_cells(territories)

    for cell_id in range(1, int(nucleus_labels.max()) + 1):
        territory = territories == cell_id
        if not territory.any():
            continue
        # Crop to territory bbox: skeletonize / EDT are local, so cropping is exact
        # but keeps per-cell cost cell-sized instead of full-frame.
        ys, xs = np.nonzero(territory)
        y0, y1 = int(ys.min()), int(ys.max()) + 1
        x0, x1 = int(xs.min()), int(xs.max()) + 1
        terr_c = territory[y0:y1, x0:x1]
        nuc_c = nucleus_labels[y0:y1, x0:x1] == cell_id
        if not nuc_c.any():
            continue
        mito_c = mito_mask[y0:y1, x0:x1] & terr_c

        skeleton = morphology.skeletonize(mito_c)
        skel_metrics = compute_skeleton_metrics(skeleton, mito_c, px_um)
        mito_area_um2 = float(mito_c.sum()) * px_area_um2
        territory_area_um2 = float(terr_c.sum()) * px_area_um2
        mito_solidity, mito_form_factor = mito_shape_metrics(mito_c)

        score = focus_score(nucleus_mip[y0:y1, x0:x1], lap[y0:y1, x0:x1], terr_c)
        mito_px = int(mito_c.sum())
        mito_border_frac = (
            float((mito_c & territory_boundaries[y0:y1, x0:x1]).sum()) / float(mito_px)
            if mito_px > 0
            else 0.0
        )
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
            "mito_area_um2": mito_area_um2,
            "mito_area_fraction": mito_area_um2 / territory_area_um2
            if territory_area_um2 > 0 else 0.0,
            "skeleton_length_um": skel_metrics["total_length_um"],
            "branch_count": skel_metrics["branch_count"],
            "mito_solidity": mito_solidity,
            "mito_form_factor": mito_form_factor,
            "mito_border_frac": mito_border_frac,
            "focus_score": score,
            "edge_distance_px": min_edge_distance,
            "footprint_confidence": footprint_confidence,
            "border_touch": cell_id in border_ids,
            "empty": mito_area_um2 < cfg.mito_min_um2,
        })
    return rows


def _apply_focus_flag(rows: list[dict], cfg: ExperimentConfig) -> float:
    """Flag out_of_focus per-image: below focus_rel_factor * median (real cells)."""
    real = [r["focus_score"] for r in rows if not r["border_touch"] and not r["empty"]]
    median = float(np.median(real)) if real else 0.0
    threshold = cfg.focus_rel_factor * median
    for r in rows:
        r["out_of_focus"] = r["focus_score"] < threshold
        r["qc_pass"] = not (r["border_touch"] or r["out_of_focus"] or r["empty"])
    return threshold


def _image_summary(cfg, group, replicate, filestem, px_um, li_threshold, focus_threshold, rows):
    """Summarize cell counts, QC counts, and metric medians."""
    pass_rows = [r for r in rows if r["qc_pass"]]
    objective_rows = [r for r in rows if not r["border_touch"]]  # empty is content, not quality
    summary: dict[str, float | int | str] = {
        "experiment": cfg.experiment, "group": group, "replicate": replicate,
        "sample": filestem, "px_um": px_um, "li_threshold": li_threshold,
        "focus_threshold": focus_threshold, "cell_N": len(rows),
        "n_qc_pass": len(pass_rows), "n_objective_pass": len(objective_rows),
        "n_border": sum(bool(r["border_touch"]) for r in rows),
        "n_out_of_focus": sum(bool(r["out_of_focus"]) for r in rows),
        "n_empty": sum(bool(r["empty"]) for r in rows),
    }
    for field in METRIC_FIELDS:
        values = [float(r[field]) for r in pass_rows]
        summary[f"median_{field}"] = float(np.median(values)) if values else np.nan
        obj = [float(r[field]) for r in objective_rows]
        summary[f"objective_median_{field}"] = float(np.median(obj)) if obj else np.nan
    return summary


def _qc_row(r: dict) -> dict:
    """Project a metrics row onto the Module A->B QC contract (QC_FIELDS)."""
    return {
        "experiment": r["experiment"], "group": r["group"],
        "replicate": r["replicate"], "sample": r["sample"],
        "cell_id": r["cell_id"], "nucleus_id": r["nucleus_id"],
        "area_um2": r["territory_area_um2"],
        "border_touch": r["border_touch"],
        "footprint_confidence": r["footprint_confidence"],
    }


def _write_skipped(cfg, filestem, manifest_row, reason):
    row = {
        "experiment": cfg.experiment, "sample": filestem,
        "px_um": float(manifest_row["px_x_um"]), "size_z": int(manifest_row["size_z"]),
        "skipped": True, "skip_reason": reason,
    }
    out = cfg.results / "04_cellqc" / "_skipped" / cfg.experiment / f"{filestem}_skipped.csv"
    write_csv(out, [row], SKIPPED_FIELDS)


def compute_one(cfg: ExperimentConfig, filestem: str, manifest_rows: dict) -> None:
    """Run cell-QC computation for one sample."""
    manifest_row = manifest_rows[(cfg.experiment, f"{filestem}.czi")]
    px_um = float(manifest_row["px_x_um"])
    if int(manifest_row["size_z"]) == 1:
        reason = "size_z == 1; cannot MIP"
        _write_skipped(cfg, filestem, manifest_row, reason)
        print(f"WARNING: skipping {cfg.experiment}/{filestem} ({reason})")
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

    import tifffile
    from .io import channel_mip

    mito_mip = tifffile.imread(mip_path)
    nucleus_mip = to_uint8(channel_mip(czi, cfg.nucleus_channel_index)[0])
    mito_mask, li_threshold = pipeline_li(mito_mip, px_um, cfg)
    mito_mask = mito_mask.astype(bool)
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

    sample_dir = cfg.results / "04_cellqc" / group / filestem
    (sample_dir / "masks").mkdir(parents=True, exist_ok=True)
    np.save(sample_dir / "masks" / "nuclei.npy", nucleus_labels)
    np.save(sample_dir / "masks" / "territories.npy", territories)

    rows = _cell_metrics(
        cfg, group, replicate, filestem, px_um,
        nucleus_mip, mito_mask, nucleus_labels, territories, footprint,
    )
    focus_threshold = _apply_focus_flag(rows, cfg)
    summary = _image_summary(
        cfg, group, replicate, filestem, px_um, li_threshold, focus_threshold, rows
    )

    table_dir = sample_dir / "tables"
    write_csv(table_dir / f"{filestem}_cells.csv", rows, CELL_FIELDS)
    write_csv(table_dir / f"{filestem}_cellqc.csv", [_qc_row(r) for r in rows], QC_FIELDS)
    write_csv(table_dir / f"{filestem}_image_summary.csv", [summary], SUMMARY_FIELDS)
    print(f"cellqc {cfg.experiment}/{filestem}: {len(rows)} cells, "
          f"{summary['n_qc_pass']} qc_pass", flush=True)


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Batch cell-QC over requested samples (or all raw CZI)."""
    manifest_rows = read_manifest(cfg)
    for filestem, _czi in resolve_samples(cfg, samples):
        if filestem in cfg.skip_samples:
            print(f"WARNING: skipping {cfg.experiment}/{filestem} (config skip_samples)")
            continue
        try:
            compute_one(cfg, filestem, manifest_rows)
        except (KeyError, ValueError) as exc:
            print(f"ERROR: {cfg.experiment}/{filestem}: {exc}")
