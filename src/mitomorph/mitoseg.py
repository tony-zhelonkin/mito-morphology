"""Mitochondria segmentation as its own compute stage.

Reads the mito max-Z MIP from ``02_extract``, thresholds it with the shared
Li pipeline (pure skimage, no GPU), and saves the boolean mito mask plus
image-level and per-object shape tables under ``03_mitoseg``.

Single-Z images are skipped (no MIP), mirroring cellqc's skip behavior.
"""
from __future__ import annotations

import numpy as np
import tifffile

from .config import ExperimentConfig
from .io import read_manifest, resolve_samples, write_csv
from .metrics import mito_object_props
from skimage import measure

IMAGE_FIELDS = [
    "experiment", "group", "replicate", "sample", "px_um", "li_threshold",
    "mito_area_um2", "mito_area_fraction_frame", "n_objects",
]

OBJECT_FIELDS = [
    "experiment", "group", "replicate", "sample", "object_id",
    "area_um2", "perimeter_um", "solidity", "eccentricity",
    "major_axis_um", "minor_axis_um", "form_factor", "centroid_y", "centroid_x",
]

SKIPPED_FIELDS = ["experiment", "sample", "px_um", "size_z", "skipped", "skip_reason"]


def _write_skipped(cfg, filestem, manifest_row, reason):
    row = {
        "experiment": cfg.experiment, "sample": filestem,
        "px_um": float(manifest_row["px_x_um"]), "size_z": int(manifest_row["size_z"]),
        "skipped": True, "skip_reason": reason,
    }
    out = cfg.results / "03_mitoseg" / "_skipped" / cfg.experiment / f"{filestem}_skipped.csv"
    write_csv(out, [row], SKIPPED_FIELDS)


def compute_one(cfg: ExperimentConfig, filestem: str, manifest_rows: dict) -> None:
    """Segment mitochondria for one sample and write mask + tables."""
    manifest_row = manifest_rows[(cfg.experiment, f"{filestem}.czi")]
    px_um = float(manifest_row["px_x_um"])
    if int(manifest_row["size_z"]) == 1:
        reason = "size_z == 1; cannot MIP"
        _write_skipped(cfg, filestem, manifest_row, reason)
        print(f"WARNING: skipping {cfg.experiment}/{filestem} ({reason})")
        return

    group, replicate = cfg.parse_sample(filestem)
    mip_path = cfg.results / "02_extract" / group / filestem / "figures" / f"{filestem}_mito_mip.tif"
    if not mip_path.exists():
        print(f"ERROR: missing mito MIP: {mip_path}. Run `mito extract` first.")
        return

    from .segment import pipeline_li

    mip = tifffile.imread(mip_path)
    mito_mask, li_threshold = pipeline_li(mip, px_um, cfg)
    mito_mask = mito_mask.astype(bool)

    sample_dir = cfg.results / "03_mitoseg" / group / filestem
    (sample_dir / "masks").mkdir(parents=True, exist_ok=True)
    np.save(sample_dir / "masks" / "mito.npy", mito_mask)

    object_rows = mito_object_props(mito_mask, px_um)
    for r in object_rows:
        r.update({
            "experiment": cfg.experiment, "group": group,
            "replicate": replicate, "sample": filestem,
        })
    n_objects = int(measure.label(mito_mask, background=0).max())
    mito_area_um2 = float(mito_mask.sum()) * px_um**2

    image_row = {
        "experiment": cfg.experiment, "group": group, "replicate": replicate,
        "sample": filestem, "px_um": px_um, "li_threshold": li_threshold,
        "mito_area_um2": mito_area_um2,
        "mito_area_fraction_frame": float(mito_mask.mean()),
        "n_objects": n_objects,
    }

    table_dir = sample_dir / "tables"
    write_csv(table_dir / f"{filestem}_mito_image.csv", [image_row], IMAGE_FIELDS)
    write_csv(table_dir / f"{filestem}_mito_objects.csv", object_rows, OBJECT_FIELDS)
    print(f"mitoseg {cfg.experiment}/{filestem}: {n_objects} objects, "
          f"mito {mito_area_um2:.1f} um2", flush=True)


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Batch mitochondria segmentation over requested samples (or all raw CZI)."""
    manifest_rows = read_manifest(cfg)
    for filestem, _czi in resolve_samples(cfg, samples):
        if filestem in cfg.skip_samples:
            print(f"WARNING: skipping {cfg.experiment}/{filestem} (config skip_samples)")
            continue
        try:
            compute_one(cfg, filestem, manifest_rows)
        except (KeyError, ValueError) as exc:
            print(f"ERROR: {cfg.experiment}/{filestem}: {exc}")
