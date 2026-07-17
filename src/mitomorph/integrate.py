"""Integrate compute: experiment-level master tables from all quantify outputs.

Reshapes only — no metrics are recomputed. Concatenates every per-sample
`05_quantify/**/tables/<stem>_cells.csv` and `<stem>_image_summary.csv` across
ALL groups/samples of the experiment into two master tables under
`03_results/master/`, decorated with clean design columns derived from the
existing `group` column (see `_design_columns`). This is the input for
SuperPlots (cell-level) and replicate-level stats (image-level).

When ``samples`` is a non-empty list, both master tables are filtered to rows
whose ``sample`` column is in ``samples`` (so `mito all --samples X` integrates
only X, matching what `quantify` just processed, instead of folding in stale
on-disk outputs). When ``samples`` is None/empty, every quantify output on disk
for the experiment is integrated.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .config import ExperimentConfig
from .provenance import _mitomorph_version, git_provenance

MASTER_DIR = "master"
KNOWN_EXPERIMENTS = ("20260518_TBK1i", "20251029_glutamine")


def _design_columns(experiment: str, group: str) -> dict:
    """Derive clean design columns from the ``group`` label, keyed on experiment id.

    Per-row: an unparseable group yields a NaN factor for THAT row only (never
    drops the column for the whole experiment). Unknown experiment ids get no
    factor columns here (the caller still adds `condition`/`replicate`).
    """
    if experiment == "20260518_TBK1i":
        # groups: DMSO, ISD90_Only, TBK1i_Only, TBK1i_ISD90
        return {
            "tbk1i": int(group.startswith("TBK1i")),
            "isd90": int("ISD90" in group),
        }
    if experiment == "20251029_glutamine":
        # groups like 0mM_GLN, 0.2mM_GLN, 0.6mM_GLN, 2mM_GLN
        m = re.match(r"([0-9.]+)mM", group)
        # NaN fallback: one bad group NaNs only its own rows, not the column.
        return {"gln_mM": float(m.group(1)) if m else np.nan}
    return {}


def _decorate(df: pd.DataFrame, experiment: str) -> pd.DataFrame:
    """Add `condition`/`replicate` (+ experiment-specific factors) to a master table."""
    df = df.copy()
    df["condition"] = df["group"]
    # `replicate` already exists in the source tables; keep as-is.
    extra_cols = df["group"].apply(lambda g: _design_columns(experiment, g)).apply(pd.Series)
    if not extra_cols.empty:
        df = pd.concat([df, extra_cols], axis=1)
        # A KNOWN experiment must reliably get its factor: warn loudly (naming the
        # offending groups) if any factor value came out NaN — never silently drop.
        if experiment in KNOWN_EXPERIMENTS:
            for col in extra_cols.columns:
                bad = df.loc[df[col].isna(), "group"].unique().tolist()
                if bad:
                    print(f"WARNING: integrate: {experiment}: could not derive '{col}' for "
                          f"group(s) {bad} (left NaN)")
    elif experiment not in KNOWN_EXPERIMENTS:
        print(f"WARNING: integrate: unknown experiment id {experiment!r}; "
              f"only adding condition/replicate (no tbk1i/isd90/gln_mM columns)")
    return df


def _collect(cfg: ExperimentConfig, suffix: str) -> tuple[pd.DataFrame, list[Path]]:
    """Concatenate every `<stem>{suffix}` table under 05_quantify/**/tables/."""
    paths = sorted((cfg.results / "05_quantify").glob(f"*/*/tables/*{suffix}"))
    frames = [pd.read_csv(p) for p in paths]
    if not frames:
        return pd.DataFrame(), paths
    return pd.concat(frames, ignore_index=True), paths


def run(cfg: ExperimentConfig, samples: list[str] | None = None) -> None:
    """Build experiment-level master_cells.csv / master_image.csv / master_manifest.yaml.

    When ``samples`` is a non-empty list, both master tables are filtered to rows
    whose ``sample`` column is in ``samples`` (so `mito all --samples X` integrates
    only X, not stale on-disk outputs). When None/empty, all quantify outputs on
    disk for the experiment are integrated.
    """
    cells_df, cells_paths = _collect(cfg, "_cells.csv")
    image_df, image_paths = _collect(cfg, "_image_summary.csv")

    if cells_df.empty or image_df.empty:
        print(f"ERROR: no quantify outputs found under {cfg.results / '05_quantify'}. "
              f"Run `mito quantify` first.")
        return

    if samples:
        wanted = {Path(s.split("/", 1)[-1]).stem for s in samples}
        cells_df = cells_df[cells_df["sample"].isin(wanted)].reset_index(drop=True)
        image_df = image_df[image_df["sample"].isin(wanted)].reset_index(drop=True)

    cells_df = _decorate(cells_df, cfg.experiment)
    image_df = _decorate(image_df, cfg.experiment)

    out_dir = cfg.results / MASTER_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    cells_out = out_dir / "master_cells.csv"
    image_out = out_dir / "master_image.csv"
    cells_df.to_csv(cells_out, index=False)
    image_df.to_csv(image_out, index=False)

    _write_manifest(cfg, cells_paths, image_paths, len(cells_df), len(image_df), out_dir)
    print(f"integrate {cfg.experiment}: {len(image_df)} samples, {len(cells_df)} cells "
          f"-> {cells_out}, {image_out}", flush=True)


def _write_manifest(
    cfg: ExperimentConfig,
    cells_paths: list[Path],
    image_paths: list[Path],
    n_cells: int,
    n_samples: int,
    out_dir: Path,
) -> None:
    """Provenance for the master tables: mitomorph version/sha, source stage, sources."""
    import yaml

    prov = git_provenance()
    manifest = {
        "experiment": cfg.experiment,
        "config_version": cfg.config_version,
        "source_stage": "05_quantify",
        "n_samples": n_samples,
        "n_cells": n_cells,
        "mitomorph": {
            "version": _mitomorph_version(),
            "git_sha": prov["git_sha"],
            "git_dirty": prov["git_dirty"],
            "git_source": prov["source"],
        },
        "source_files": {
            "cells": [str(p.relative_to(cfg.root)) for p in cells_paths],
            "image_summary": [str(p.relative_to(cfg.root)) for p in image_paths],
        },
    }
    out = out_dir / "master_manifest.yaml"
    out.write_text(yaml.safe_dump(manifest, sort_keys=False))
