"""Cell-QC visualization: per-sample territory overlays and an experiment montage.

Every table's source figure is saved next to it (compute/viz pairing).
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from skimage import segmentation  # noqa: E402

from .config import ExperimentConfig  # noqa: E402
from .io import read_manifest, resolve_samples  # noqa: E402


def _is_true(value: str) -> bool:
    return value.lower() in {"true", "1", "yes"}


def _read_cells(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _normalize_gray(img: np.ndarray) -> np.ndarray:
    f = img.astype(float)
    lo, hi = np.percentile(f, (0.5, 99.5))
    if hi > lo:
        return np.clip((f - lo) / (hi - lo), 0, 1)
    if f.max() > f.min():
        return (f - f.min()) / (f.max() - f.min())
    return np.zeros_like(f)


def _centroids(labels: np.ndarray) -> dict[int, tuple[float, float]]:
    out = {}
    for cell_id in range(1, int(labels.max()) + 1):
        coords = np.argwhere(labels == cell_id)
        if coords.size:
            yx = coords.mean(axis=0)
            out[cell_id] = (float(yx[0]), float(yx[1]))
    return out


def _qc_reasons(row: dict[str, str]) -> list[str]:
    reasons = []
    if _is_true(row["border_touch"]):
        reasons.append("border")
    if _is_true(row["out_of_focus"]):
        reasons.append("focus")
    if _is_true(row["empty"]):
        reasons.append("empty")
    return reasons


def _overlay_rgb(mip: np.ndarray, nuclei: np.ndarray, territories: np.ndarray) -> np.ndarray:
    gray = _normalize_gray(mip)
    rgb = np.dstack([gray, gray, gray])
    territory_edges = segmentation.find_boundaries(territories, mode="inner") & (territories > 0)
    nucleus_edges = segmentation.find_boundaries(nuclei, mode="inner") & (nuclei > 0)
    rgb[territory_edges] = [0.0, 0.85, 1.0]
    rgb[nucleus_edges] = [1.0, 0.85, 0.0]
    return rgb


def _add_panel(ax, mip, nuclei, territories, rows, title) -> None:
    ax.imshow(_overlay_rgb(mip, nuclei, territories))
    centroids = _centroids(nuclei)
    flagged = []
    for row in rows:
        cell_id = int(row["cell_id"])
        centroid = centroids.get(cell_id)
        if centroid is None:
            continue
        y, x = centroid
        ax.text(x, y, str(cell_id), color="white", fontsize=5, ha="center", va="center",
                bbox={"facecolor": "black", "alpha": 0.55, "pad": 0.5, "linewidth": 0})
        reasons = _qc_reasons(row)
        if reasons:
            flagged.append(cell_id)
            ax.text(x, y + 16, "X " + ",".join(reasons), color="tab:red", fontsize=5,
                    ha="center", va="center",
                    bbox={"facecolor": "white", "alpha": 0.75, "pad": 0.5, "linewidth": 0})
    for cell_id in flagged:
        ax.contour(territories == cell_id, levels=[0.5], colors=["tab:red"], linewidths=0.5)
    ax.set_title(title, fontsize=9)
    ax.set_axis_off()


def _load(cfg: ExperimentConfig, filestem: str):
    group, _ = cfg.parse_sample(filestem)
    mip = cfg.results / "02_extract" / group / filestem / "figures" / f"{filestem}_mito_mip.tif"
    sample = cfg.results / "04_cellqc" / group / filestem
    table = sample / "tables" / f"{filestem}_cells.csv"
    nuclei = sample / "masks" / "nuclei.npy"
    terr = sample / "masks" / "territories.npy"
    for p in (mip, table, nuclei, terr):
        if not p.exists():
            print(f"ERROR: missing {p}. Run `mito extract` and `mito cellqc` first.")
            return None
    return (group, tifffile.imread(mip), np.load(nuclei), np.load(terr), _read_cells(table))


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Render per-sample overlays and a combined montage."""
    manifest_rows = read_manifest(cfg)
    panels = []
    for filestem, _czi in resolve_samples(cfg, samples):
        if filestem in cfg.skip_samples:
            continue
        manifest_row = manifest_rows.get((cfg.experiment, f"{filestem}.czi"))
        if manifest_row is None or int(manifest_row["size_z"]) == 1:
            continue
        loaded = _load(cfg, filestem)
        if loaded is None:
            continue
        group, mip, nuclei, territories, rows = loaded
        fig_dir = cfg.results / "04_cellqc" / group / filestem / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
        _add_panel(ax, mip, nuclei, territories, rows, filestem)
        fig.savefig(fig_dir / f"{filestem}_cellqc_overlay.png", dpi=180)
        plt.close(fig)
        panels.append((filestem, mip, nuclei, territories, rows))

    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5.2),
                             constrained_layout=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (filestem, mip, nuclei, territories, rows) in zip(axes, panels):
        _add_panel(ax, mip, nuclei, territories, rows, filestem)
    out = cfg.results / "04_cellqc" / "04_cellqc_montage.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"wrote {len(panels)} overlays + montage -> {out}")
