"""Shared viz helpers for the per-stage overlay modules.

Split out of the original monolithic ``viz.py`` so each stage's viz module
(mitoseg/cellseg/quantify) reuses the same normalization/overlay primitives.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from skimage import segmentation


def is_true(value: str) -> bool:
    return value.lower() in {"true", "1", "yes"}


def read_cells(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def normalize_gray(img: np.ndarray) -> np.ndarray:
    f = img.astype(float)
    lo, hi = np.percentile(f, (0.5, 99.5))
    if hi > lo:
        return np.clip((f - lo) / (hi - lo), 0, 1)
    if f.max() > f.min():
        return (f - f.min()) / (f.max() - f.min())
    return np.zeros_like(f)


def centroids(labels: np.ndarray) -> dict[int, tuple[float, float]]:
    out = {}
    for cell_id in range(1, int(labels.max()) + 1):
        coords = np.argwhere(labels == cell_id)
        if coords.size:
            yx = coords.mean(axis=0)
            out[cell_id] = (float(yx[0]), float(yx[1]))
    return out


def overlay_rgb(mip: np.ndarray, nuclei: np.ndarray, territories: np.ndarray) -> np.ndarray:
    gray = normalize_gray(mip)
    rgb = np.dstack([gray, gray, gray])
    territory_edges = segmentation.find_boundaries(territories, mode="inner") & (territories > 0)
    nucleus_edges = segmentation.find_boundaries(nuclei, mode="inner") & (nuclei > 0)
    rgb[territory_edges] = [0.0, 0.85, 1.0]
    rgb[nucleus_edges] = [1.0, 0.85, 0.0]
    return rgb


def qc_reasons(row: dict[str, str]) -> list[str]:
    reasons = []
    if is_true(row["border_touch"]):
        reasons.append("border")
    if is_true(row["out_of_focus"]):
        reasons.append("focus")
    if is_true(row["empty"]):
        reasons.append("empty")
    return reasons
