"""Cell-segmentation visualization: per-sample nuclei/territory edge overlays.

Reads the shared MIP (02_extract) and the nuclei/territory masks written by
`mito cellseg`; draws territory + nucleus boundaries with nucleus centroids,
WITHOUT any mito or QC-flag annotations (compute/viz pairing).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import tifffile  # noqa: E402

from .config import ExperimentConfig  # noqa: E402
from .io import read_manifest, resolve_samples  # noqa: E402
from .viz_common import centroids, overlay_rgb  # noqa: E402


def _add_panel(ax, mip, nuclei, territories, title) -> None:
    ax.imshow(overlay_rgb(mip, nuclei, territories))
    for cell_id, (y, x) in centroids(nuclei).items():
        ax.text(x, y, str(cell_id), color="white", fontsize=5, ha="center", va="center",
                bbox={"facecolor": "black", "alpha": 0.55, "pad": 0.5, "linewidth": 0})
    ax.set_title(title, fontsize=9)
    ax.set_axis_off()


def _load(cfg: ExperimentConfig, filestem: str):
    group, _ = cfg.parse_sample(filestem)
    mip = cfg.results / "02_extract" / group / filestem / "figures" / f"{filestem}_mito_mip.tif"
    sample = cfg.results / "04_cellseg" / group / filestem
    nuclei = sample / "masks" / "nuclei.npy"
    terr = sample / "masks" / "territories.npy"
    for p in (mip, nuclei, terr):
        if not p.exists():
            print(f"ERROR: missing {p}. Run `mito extract` and `mito cellseg` first.")
            return None
    return (group, tifffile.imread(mip), np.load(nuclei), np.load(terr))


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Render per-sample cell-segmentation overlays and a combined montage."""
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
        group, mip, nuclei, territories = loaded
        fig_dir = cfg.results / "04_cellseg" / group / filestem / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
        _add_panel(ax, mip, nuclei, territories, filestem)
        fig.savefig(fig_dir / f"{filestem}_cellseg_overlay.png", dpi=180)
        plt.close(fig)
        panels.append((filestem, mip, nuclei, territories))

    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5.2),
                             constrained_layout=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (filestem, mip, nuclei, territories) in zip(axes, panels):
        _add_panel(ax, mip, nuclei, territories, filestem)
    out = cfg.results / "04_cellseg" / "04_cellseg_montage.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"wrote {len(panels)} overlays + montage -> {out}")
