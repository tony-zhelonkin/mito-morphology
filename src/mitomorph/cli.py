"""`mito` command-line entry point for the shared pipeline.

Runs against the experiment.toml in the current directory (or --root):
  mito inventory
  mito extract   [--samples ...]
  mito cellqc    [--samples ...]
  mito cellqc-viz [--samples ...]
  mito all       [--samples ...]   # inventory -> extract -> cellqc -> cellqc-viz
"""
from __future__ import annotations

import argparse

from . import cellqc, extract, inventory, viz
from .config import load_config


def main() -> None:
    p = argparse.ArgumentParser(prog="mito", description="Mitochondrial morphology pipeline.")
    p.add_argument("--root", default=".", help="Experiment repo root (has experiment.toml).")
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("inventory", "extract", "cellqc", "cellqc-viz", "all"):
        sp = sub.add_parser(name)
        if name != "inventory":
            sp.add_argument("--samples", nargs="*", help="Optional filestems to process.")

    args = p.parse_args()
    cfg = load_config(args.root)
    samples = getattr(args, "samples", None)

    if args.command == "inventory":
        inventory.run(cfg)
    elif args.command == "extract":
        extract.run(cfg, samples)
    elif args.command == "cellqc":
        cellqc.run(cfg, samples)
    elif args.command == "cellqc-viz":
        viz.run(cfg, samples)
    elif args.command == "all":
        inventory.run(cfg)
        extract.run(cfg, samples)
        cellqc.run(cfg, samples)
        viz.run(cfg, samples)


if __name__ == "__main__":
    main()
