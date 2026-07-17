"""Mito-segmentation visualization: per-sample mito-mask overlays and a montage.

Reads the shared MIP (02_extract) and the mito mask written by `mito mitoseg`;
overlays the mask boundaries on the normalized MIP (compute/viz pairing).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from skimage import segmentation  # noqa: E402

from .config import ExperimentConfig  # noqa: E402
from .io import read_manifest, resolve_samples  # noqa: E402
from .viz_common import normalize_gray  # noqa: E402


def _overlay_rgb(mip: np.ndarray, mito: np.ndarray) -> np.ndarray:
    gray = normalize_gray(mip)
    rgb = np.dstack([gray, gray, gray])
    mask = mito > 0
    # Translucent green fill for the mito mask, bright green boundaries on top.
    rgb[mask] = 0.6 * rgb[mask] + 0.4 * np.array([0.0, 0.9, 0.3])
    edges = segmentation.find_boundaries(mask, mode="inner") & mask
    rgb[edges] = [0.0, 1.0, 0.4]
    return rgb


def _add_panel(ax, mip, mito, title) -> None:
    ax.imshow(_overlay_rgb(mip, mito))
    ax.set_title(title, fontsize=9)
    ax.set_axis_off()


def _load(cfg: ExperimentConfig, filestem: str):
    group, _ = cfg.parse_sample(filestem)
    mip = cfg.results / "02_extract" / group / filestem / "figures" / f"{filestem}_mito_mip.tif"
    mito = cfg.results / "03_mitoseg" / group / filestem / "masks" / "mito.npy"
    for p in (mip, mito):
        if not p.exists():
            print(f"ERROR: missing {p}. Run `mito extract` and `mito mitoseg` first.")
            return None
    return (group, tifffile.imread(mip), np.load(mito))


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Render per-sample mito-mask overlays and a combined montage."""
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
        group, mip, mito = loaded
        fig_dir = cfg.results / "03_mitoseg" / group / filestem / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
        _add_panel(ax, mip, mito, filestem)
        fig.savefig(fig_dir / f"{filestem}_mito_mask_overlay.png", dpi=180)
        plt.close(fig)
        panels.append((filestem, mip, mito))

    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5.2),
                             constrained_layout=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (filestem, mip, mito) in zip(axes, panels):
        _add_panel(ax, mip, mito, filestem)
    out = cfg.results / "03_mitoseg" / "03_mitoseg_montage.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"wrote {len(panels)} overlays + montage -> {out}")
