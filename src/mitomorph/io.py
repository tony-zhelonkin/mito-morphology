"""Shared I/O: manifest reading, sample resolution, CZI channel projection."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from bioio import BioImage

from .config import ExperimentConfig


def read_manifest(cfg: ExperimentConfig) -> dict[tuple[str, str], dict[str, str]]:
    """Read per-file image metadata; this is the source of pixel sizes."""
    with cfg.manifest.open(newline="") as fh:
        rows = csv.DictReader(fh)
        return {(r["experiment"], r["file"]): r for r in rows}


def resolve_samples(cfg: ExperimentConfig, samples: list[str] | None) -> list[tuple[str, Path]]:
    """Resolve --samples filestems (or all raw CZI) to (filestem, czi path)."""
    if samples:
        out = []
        for item in samples:
            filestem = Path(item.split("/", 1)[-1]).stem
            out.append((filestem, cfg.raw / f"{filestem}.czi"))
        return out
    return [(czi.stem, czi) for czi in sorted(cfg.raw.glob("*.czi"))]


def to_uint8(img: np.ndarray) -> np.ndarray:
    """Scale an image to uint8 using its current numeric range."""
    if img.dtype == np.uint8:
        return img
    f = img.astype(np.float32)
    if f.max() > f.min():
        f = (f - f.min()) / (f.max() - f.min()) * 255.0
    return np.clip(f, 0, 255).astype(np.uint8)


def channel_mip(czi: Path, channel_index: int) -> tuple[np.ndarray, list[str]]:
    """Return a raw max-Z projection for one channel plus channel names."""
    img = BioImage(czi)
    data = img.get_image_data("CZYX")
    mip = np.max(data[channel_index], axis=0)
    try:
        names = list(img.channel_names)
    except Exception:
        names = []
    return mip, names


def radius_px(radius_um: float, px_um: float) -> int:
    """Convert a physical radius to at least one pixel."""
    return max(1, int(radius_um / px_um))


def area_px(area_um2: float, px_um: float) -> int:
    """Convert a physical area to at least one pixel."""
    return max(1, int(math.ceil(area_um2 / (px_um**2))))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    """Write rows to CSV with a fixed header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
