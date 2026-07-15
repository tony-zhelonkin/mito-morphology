# mito-morphology

Shared compute for mitochondrial network morphology in ISD90-activated BMDCs.

Installable package (`mitomorph`) used by each experiment repo. Everything
experiment-specific lives in that repo's `experiment.toml`; the code here is generic.

## Pipeline

```
mito inventory       # probe raw CZI -> 03_results/00_inventory/tables/czi_manifest.csv
mito extract         # mito channel max-Z MIP -> 03_results/02_extract/
mito cellqc          # nucleus-seeded territories + per-cell metrics -> 03_results/04_cellqc/
mito cellqc-viz      # territory overlays + montage (source figure for each table)
mito all             # the four steps above
```

Run from an experiment repo root (containing `experiment.toml`), or pass `--root`.

## Design

- **Compute ≠ viz**: compute writes tables/masks; viz reads them and saves the figure each table came from.
- **QC is flag-only**: `border_touch`, `out_of_focus`, `empty` are recorded; rows are never dropped.
- **Focus is per-image relative**: `out_of_focus` uses each image's own median, never a
  cross-experiment constant (focus_score is not comparable across imaging batches).
- **Mito mask = Li** (locked bake-off winner); per-cell via nucleus-seeded, mito-independent territories.
