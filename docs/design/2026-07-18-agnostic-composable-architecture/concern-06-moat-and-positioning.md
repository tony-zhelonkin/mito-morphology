# Concern 06 — Moat and Positioning

**Question:** What is mitomorph's *moat*? What justifies building it rather than
adopting Nellie / cellpose / empanada / napari-mito-hcs / CellProfiler?

**Method:** Read mitomorph source (`src/mitomorph/*.py`), read napari-mito-hcs
source (`init-draft/01_modules/napari-mito-hcs/src/`), read the experiments' R
stats (`04_stats/{contrasts,run_stats}.R`), and checked current (2025-2026)
capabilities of the competitors via web sources (cited at end).

**Verdict up front:** the moat is **not** segmentation, **not** per-cell
morphometry, and **not** "composable + napari-observable" (napari-mito-hcs
already ships all of that). The one genuine, defensible gap is the **design-aware
replicate-level inference layer** — factorial/dose contrasts on per-image medians
with transparent per-cell SuperPlots — that closes the loop from image to
publishable statistic. Everything upstream of that is an integration convenience
that assembles existing pieces; nobody has assembled *exactly* this loop, but
they could. The stats layer is the part no imaging tool does at all.

---

## What mitomorph actually does (from source)

End-to-end: CZI → MIP-collapse to 2D → **Li**-threshold mito mask
(`segment.pipeline_li`, physical min-object floor scaled by `px_um`) →
Cellpose `cyto3` whole-cell bodies on [enhanced-mito, nucleus] + Otsu/watershed
nucleus seeds → **nucleus-reconciled** cell territories (`resolve_instances`:
trust Cellpose boundary for 1-nucleus bodies, geodesic-watershed split only
genuine multi-nucleus merges, radius-disk fallback for missed nuclei) → per-CELL
mito morphometry (`metrics.py`: **skeleton network** — length, branch/junction/
endpoint counts, mean branch length/diameter; **shape** — area-weighted solidity
and form factor) → objective, condition-blind QC gating (border_touch,
out_of_focus relative to per-image DAPI median, empty; **flag-only, never
drops**) → per-IMAGE median aggregation (`_image_summary`, both `median_` on
qc_pass and `objective_median_` on non-border) → per-CONDITION replicate-level
contrasts in R (`emmeans`, per-metric `lm()`, explicit factorial/dose contrast
vectors, primary/secondary families) → SuperPlots. Provenance: git-sha manifest
+ pixi lockfile; masks are `.npy`, every table has a paired figure/MIP source.

## The closest competitor: napari-mito-hcs (Denali, SLAS Discovery 2025)

This is the honest threat, and it is much closer than cellpose/empanada. Its
`MitoHCSPipeline` runs, as **separate, config-driven, napari-observable steps**:
`segment_nuclei` → `segment_cells` (nucleus-guided splitting of touching cells)
→ `segment_mitochondria` (**assigns mito to a parent cell** via most-common
parent label, `stats.calc_parent_label`) → `calc_features` (shape-index texture:
spot/hole/ridge/valley/saddle) → `save_stats` (per-object geometry/intensity/
texture) → `calc_summary_stats` (per-FOV weighted means, SpotRidgeRatio,
AspectRatio). It has a napari GUI, TOML config save/load from an interactive
session, and a batch CLI.

So napari-mito-hcs **already** does: nuclei→cell→mito→per-cell assignment→per-cell
features→per-FOV summary, composably, observably, config-reproducibly. The
"declarative + composable + GUI-observable synthesis" framing is therefore **not
unique to mitomorph** — Denali published it in 2025.

What napari-mito-hcs does **not** do:
- Segmentation is **threshold-based** (no learned/adaptive cell boundaries; no
  Cellpose in the loop) — brittle on touching/faint cells vs mitomorph's
  Cellpose + nucleus reconciliation.
- Features are **texture (shape-index)**, not mito **network** metrics — no
  skeleton length/branch/junction graph.
- Aggregation stops at a **per-FOV summary table**. There is **no experimental
  design, no replicate-level inference, no contrasts, no SuperPlots.** It hands
  you a spreadsheet; you do the stats yourself, ad hoc.

## CellProfiler — the strongest "it already exists" argument

CellProfiler genuinely covers more of the *upstream* loop than people assume:
- **RelateObjects** gives nucleus→cell→organelle parent-child aggregation
  (per-parent means) — i.e. per-cell assignment.
- **MeasureObjectSkeleton** gives per-seed skeleton metrics: NumberTrunks,
  NumberNonTrunkBranches, NumberBranchEnds, TotalObjectSkeletonLength — i.e.
  per-cell mito **network** metrics, nucleus-seeded.
- **Metadata/Groups** modules do well/condition grouping; a **RunCellpose**
  plugin exists for adaptive seg.

So CellProfiler can, in principle, produce nucleus-seeded per-cell mito skeleton
morphometry grouped by condition. What it does **not** provide:
- **Design-aware inferential statistics.** CellProfiler *exports measurements*;
  it does not fit models or test factorial/dose contrasts. Grouping ≠ inference.
  Replicate-level contrasts + pseudoreplication-correct aggregation are entirely
  downstream and manual.
- **napari** step-observability (its own GUI, not the napari ecosystem the lab
  standardizes on).
- **Modern scripted, env-locked reproducibility** — a `.cppipe` is versioned but
  not a lockfile-pinned, git-sha-stamped, scriptable pipeline.

## Capability matrix

Legend: ✔ native · ~ partial / possible-with-effort / plugin · ✘ absent

| Capability | mitomorph | napari-mito-hcs | CellProfiler | Nellie | cellpose(-SAM) | empanada/MitoNet |
|---|---|---|---|---|---|---|
| Adaptive (learned) mito/cell seg | ✔ Cellpose+Li | ✘ threshold | ~ (RunCellpose plugin) | ✔ adaptive | ✔ cells only | ✔ EM mito only |
| Cell assignment (nucleus-seeded) | ✔ | ✔ | ✔ RelateObjects | ✘ | ✘ | ✘ |
| Mito **network** metrics (skeleton graph) | ✔ | ✘ (texture) | ~ MeasureObjectSkeleton | ✔ (best-in-class) | ✘ | ✘ |
| Design/condition-aware **stats** (contrasts, replicate-level) | ✔ R/emmeans | ✘ | ✘ (export only) | ✘ | ✘ | ✘ |
| Composable / reorderable steps | ✔ | ✔ | ✔ | ~ | n/a | n/a |
| napari step-observability | ~ (vision) | ✔ | ✘ | ✔ | ~ | ~ |
| Provenance + env-lock (pixi/git-sha) | ✔ | ~ (TOML cfg) | ~ (.cppipe) | ✘ | ✘ | ✘ |
| 2D-MIP fit | ✔ | ✔ (2D) | ✔ | ✔ | ✔ | ✘ (3D EM) |
| Handles anisotropy (CZI→MIP) | ✔ | ✘ (pre-made 2D) | ~ | ✔ (metadata-aware) | ~ | ~ |
| Transparent per-cell SuperPlots | ✔ | ✘ | ✘ | ✘ | ✘ | ✘ |

Two columns are decisive. The **"design-aware stats"** and **"SuperPlots"** rows
are ✔ for mitomorph and ✘ for *every* imaging tool. Every other row is matched by
at least one competitor.

## Blunt verdict on the moat

The honest classification is **(iii) with a hard (i) core**:

- The "composable + GUI-observable synthesis of existing pieces" story
  (the composable-synthesis framing) is **already realized by napari-mito-hcs** and largely by
  CellProfiler. As a *whole-loop* claim it is **(ii) integration convenience** —
  no single tool does *this exact* 2D-MIP→Cellpose-cell→per-cell-skeleton loop,
  but that is "nobody has assembled it," not "it cannot be assembled."
- The **genuine capability gap (i)** is narrow and specific: **no imaging tool
  carries per-cell morphometry through to design-aware, replicate-level
  statistical inference with pseudoreplication-correct aggregation and
  transparent per-cell SuperPlots.** napari-mito-hcs stops at a per-FOV
  spreadsheet; CellProfiler stops at exported measurements; Nellie/cellpose/
  empanada stop at features/masks. The image→inference bridge, made design-aware
  (2×2 factorial contrasts; dose sensitivity), is the part that exists nowhere.

Secondary, defensible-but-not-unique advantages: Cellpose+nucleus-reconciliation
seg is more robust than HCS's thresholds; skeleton-network metrics beat HCS's
texture for *network* phenotypes (though Nellie's are richer); pixi/git-sha
provenance exceeds all of them. None of these alone justifies a new tool — each
is available elsewhere.

## The one-sentence moat

> **mitomorph's justification is not its segmentation or its per-cell
> morphometry — both exist in napari-mito-hcs, CellProfiler, and Nellie — but
> that it is the only tool that carries nucleus-assigned per-cell mitochondrial
> network morphometry all the way through to design-aware, replicate-level
> statistical contrasts and transparent per-cell SuperPlots, closing the
> image→publishable-inference loop that every existing imaging tool leaves to
> ad-hoc downstream spreadsheet work.**

If that stats/design layer were stripped out, mitomorph would be a
(nicely-provenance'd) convenience re-assembly of napari-mito-hcs + Cellpose and
the honest recommendation would be to adopt napari-mito-hcs and bolt on Cellpose.
The stats/design layer is what makes it more than a wrapper.

## Sources
- Nellie — Nat Methods 2025: https://www.nature.com/articles/s41592-025-02612-7 ; repo https://github.com/aelefebv/nellie
- napari-mito-hcs — Chin et al., SLAS Discovery 2025: https://doi.org/10.1016/j.slasd.2025.100208 (source read locally under `init-draft/01_modules/napari-mito-hcs/`)
- CellProfiler MeasureObjectSkeleton / RelateObjects / Groups: https://cellprofiler-manual.s3.amazonaws.com/CellProfiler-4.2.8/modules/measurement.html
- Segment Anything for Microscopy (cellpose-SAM context): https://www.nature.com/articles/s41592-024-02580-4
