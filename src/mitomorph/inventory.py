"""Probe every raw CZI and record dimensions / dtype / channels / pixel sizes."""
from __future__ import annotations

from pathlib import Path

from bioio import BioImage

from .config import ExperimentConfig
from .io import write_csv

FIELDS = [
    "experiment", "file", "size_x", "size_y", "size_z", "size_c", "size_t",
    "dtype", "px_x_um", "px_y_um", "px_z_um", "channel_names",
]


def _probe(path: Path, experiment: str) -> dict:
    img = BioImage(path)
    d = img.dims
    px = img.physical_pixel_sizes  # (Z, Y, X) micrometres
    try:
        ch = ",".join(img.channel_names)
    except Exception:
        ch = ""
    return {
        "experiment": experiment, "file": path.name,
        "size_x": d.X, "size_y": d.Y, "size_z": d.Z, "size_c": d.C, "size_t": d.T,
        "dtype": str(img.dtype),
        "px_x_um": px.X, "px_y_um": px.Y, "px_z_um": px.Z,
        "channel_names": ch,
    }


def run(cfg: ExperimentConfig) -> None:
    """Write the CZI manifest for this experiment."""
    rows = []
    for czi in sorted(cfg.raw.glob("*.czi")):
        print(f"probing {cfg.experiment}/{czi.name} ...", flush=True)
        try:
            rows.append(_probe(czi, cfg.experiment))
        except Exception as e:
            print(f"  !! FAILED: {e}")
            rows.append({**{k: "" for k in FIELDS},
                         "experiment": cfg.experiment, "file": czi.name,
                         "dtype": f"ERROR: {e}"})
    write_csv(cfg.manifest, rows, FIELDS)
    print(f"\nwrote {len(rows)} rows -> {cfg.manifest}")
