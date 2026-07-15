"""Extract the mitochondrial channel from raw CZI and save a 2D max-Z MIP."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from .config import ExperimentConfig
from .io import channel_mip, read_manifest, resolve_samples, write_csv

FIELDS = [
    "experiment", "group", "replicate", "sample",
    "mito_channel_index", "mito_channel_name",
    "size_x", "size_y", "size_z", "px_x_um",
    "mip_min", "mip_max", "mip_mean",
]


def extract_one(cfg: ExperimentConfig, filestem: str, czi: Path, manifest_rows: dict) -> None:
    """Extract and write one mitochondrial max-Z projection plus metadata."""
    group, replicate = cfg.parse_sample(filestem)
    manifest_row = manifest_rows[(cfg.experiment, czi.name)]
    idx = cfg.mito_channel_index

    print(f"extracting {cfg.experiment}/{czi.name} ...", flush=True)
    mip, channel_names = channel_mip(czi, idx)
    if mip.dtype != np.uint8:
        mip = np.clip(mip, 0, 255).astype(np.uint8)
    mito_channel_name = channel_names[idx] if idx < len(channel_names) else ""

    sample_dir = cfg.results / "02_extract" / group / filestem
    fig_dir = sample_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(fig_dir / f"{filestem}_mito_mip.tif", mip)

    row = {
        "experiment": cfg.experiment, "group": group, "replicate": replicate,
        "sample": filestem, "mito_channel_index": idx, "mito_channel_name": mito_channel_name,
        "size_x": manifest_row["size_x"], "size_y": manifest_row["size_y"],
        "size_z": manifest_row["size_z"], "px_x_um": manifest_row["px_x_um"],
        "mip_min": int(mip.min()), "mip_max": int(mip.max()), "mip_mean": float(mip.mean()),
    }
    write_csv(sample_dir / "tables" / f"{filestem}_extract.csv", [row], FIELDS)


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Batch mito extraction over requested samples (or all raw CZI)."""
    manifest_rows = read_manifest(cfg)
    for filestem, czi in resolve_samples(cfg, samples):
        if filestem in cfg.skip_samples:
            print(f"WARNING: skipping {cfg.experiment}/{filestem} (config skip_samples)")
            continue
        if not czi.exists():
            print(f"ERROR: missing raw CZI: {czi}")
            continue
        try:
            extract_one(cfg, filestem, czi, manifest_rows)
        except (KeyError, ValueError) as exc:
            print(f"ERROR: {cfg.experiment}/{filestem}: {exc}")
