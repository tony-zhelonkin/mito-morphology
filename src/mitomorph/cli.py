"""`mito` command-line entry point for the shared pipeline.

Runs against the experiment.toml in the current directory (or --root):
  mito inventory
  mito to-zarr      [--samples ...]
  mito migrate-labels [--samples ...]  # .npy masks -> NGFF labels/ groups (opt-in)
  mito extract      [--samples ...]
  mito mitoseg      [--samples ...]
  mito mitoseg-viz  [--samples ...]
  mito cellseg      [--samples ...]
  mito cellseg-viz  [--samples ...]
  mito quantify     [--samples ...]
  mito quantify-viz [--samples ...]
  mito integrate    [--samples ...]
  mito all          [--samples ...]   # inventory -> extract -> mitoseg -> cellseg
                                       # -> quantify -> integrate -> *-viz
"""
from __future__ import annotations

import argparse

from . import (
    cellseg,
    cellseg_viz,
    extract,
    integrate,
    inventory,
    mitoseg,
    mitoseg_viz,
    omezarr,
    quantify,
    quantify_viz,
)
from .config import load_config
from .provenance import write_run_manifest


def main() -> None:
    p = argparse.ArgumentParser(prog="mito", description="Mitochondrial morphology pipeline.")
    p.add_argument("--root", default=".", help="Experiment repo root (has experiment.toml).")
    sub = p.add_subparsers(dest="command", required=True)
    for name in (
        "inventory",
        "to-zarr",
        "migrate-labels",
        "extract",
        "mitoseg",
        "mitoseg-viz",
        "cellseg",
        "cellseg-viz",
        "quantify",
        "quantify-viz",
        "integrate",
        "all",
    ):
        sp = sub.add_parser(name)
        if name != "inventory":
            sp.add_argument("--samples", nargs="*", help="Optional filestems to process.")

    args = p.parse_args()
    cfg = load_config(args.root)
    samples = getattr(args, "samples", None)
    # Stamp provenance for every run that produces analysis outputs (skip pure
    # inventory / zarr-conversion and viz-only stages, which record their own
    # state or none).
    if args.command in ("mitoseg", "cellseg", "quantify", "integrate", "all", "migrate-labels"):
        write_run_manifest(cfg, f"mito {args.command}", samples)

    if args.command == "inventory":
        inventory.run(cfg)
    elif args.command == "to-zarr":
        omezarr.run(cfg, samples)
    elif args.command == "migrate-labels":
        omezarr.run_migrate_labels(cfg, samples)
    elif args.command == "extract":
        extract.run(cfg, samples)
    elif args.command == "mitoseg":
        mitoseg.run(cfg, samples)
    elif args.command == "mitoseg-viz":
        mitoseg_viz.run(cfg, samples)
    elif args.command == "cellseg":
        cellseg.run(cfg, samples)
    elif args.command == "cellseg-viz":
        cellseg_viz.run(cfg, samples)
    elif args.command == "quantify":
        quantify.run(cfg, samples)
    elif args.command == "quantify-viz":
        quantify_viz.run(cfg, samples)
    elif args.command == "integrate":
        integrate.run(cfg, samples)
    elif args.command == "all":
        inventory.run(cfg)
        extract.run(cfg, samples)
        mitoseg.run(cfg, samples)
        cellseg.run(cfg, samples)
        quantify.run(cfg, samples)
        integrate.run(cfg, samples)
        mitoseg_viz.run(cfg, samples)
        cellseg_viz.run(cfg, samples)
        quantify_viz.run(cfg, samples)


if __name__ == "__main__":
    main()
