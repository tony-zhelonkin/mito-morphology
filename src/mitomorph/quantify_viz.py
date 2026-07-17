"""Quantification QC visualization: per-sample flag overlays and a montage.

Reads the shared MIP (02_extract), the territory masks (04_cellseg) and the
per-cell table written by `mito quantify`; marks flagged cells with a red
contour and short `border/focus/empty` labels (compute/viz pairing).
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
from .viz_common import centroids, normalize_gray, qc_reasons, read_cells  # noqa: E402


def _add_panel(ax, mip, territories, rows, title) -> None:
    gray = normalize_gray(mip)
    ax.imshow(np.dstack([gray, gray, gray]))
    cents = centroids(territories)
    flagged = []
    for row in rows:
        cell_id = int(row["cell_id"])
        reasons = qc_reasons(row)
        if not reasons:
            continue
        flagged.append(cell_id)
        centroid = cents.get(cell_id)
        if centroid is None:
            continue
        y, x = centroid
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
    terr = cfg.results / "04_cellseg" / group / filestem / "masks" / "territories.npy"
    table = cfg.results / "05_quantify" / group / filestem / "tables" / f"{filestem}_cells.csv"
    for p in (mip, terr, table):
        if not p.exists():
            print(f"ERROR: missing {p}. Run `mito extract`, `mito cellseg` and `mito quantify` first.")
            return None
    return (group, tifffile.imread(mip), np.load(terr), read_cells(table))


def run(cfg: ExperimentConfig, samples: list[str] | None) -> None:
    """Render per-sample quantify flag overlays and a combined montage."""
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
        group, mip, territories, rows = loaded
        fig_dir = cfg.results / "05_quantify" / group / filestem / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
        _add_panel(ax, mip, territories, rows, filestem)
        fig.savefig(fig_dir / f"{filestem}_quantify_overlay.png", dpi=180)
        plt.close(fig)
        panels.append((filestem, mip, territories, rows))

    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5.2),
                             constrained_layout=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (filestem, mip, territories, rows) in zip(axes, panels):
        _add_panel(ax, mip, territories, rows, filestem)
    out = cfg.results / "05_quantify" / "05_quantify_montage.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"wrote {len(panels)} overlays + montage -> {out}")
