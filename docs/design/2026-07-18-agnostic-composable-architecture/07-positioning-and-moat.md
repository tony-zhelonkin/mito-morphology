# Positioning & moat — why mitomorph exists, and where Nellie fits

**Date:** 2026-07-18 · **Status:** positioning statement (grounded in evidence).
**Evidence:** `concern-06-moat-and-positioning.md` (mitomorph + napari-mito-hcs source read;
CellProfiler/Nellie/cellpose/empanada current capabilities), `concern-05` (Nellie anisotropy),
`concern-02` (prior-art). This is the blunt, honest answer to "what justifies building this?"

---

## 1. The blunt verdict (revised — two differentiators, and statistics stays OUT)

After reading the competitors' source, the current (2025-26) field, AND the project's own
measurement rationale (`TBK1i.../03_results/05_quantify/README.md`, umbrella
`METHODS.md`, `docs/llm-research/2026-07-18-mito-metrics-provenance/`), the moat is
**narrow but real**, and it is **not** where intuition first pointed. Two things it
is *not*:

- It is **not** the segmentation scaffold (Cellpose/Li assembled from cellpose + skimage).
- It is **not** "composable + intent-driven + napari-observable synthesis" — that
  exact framing was **published by Denali as napari-mito-hcs (SLAS Discovery 2025)**.
  The vision is *good* but, as a bare workflow differentiator, already taken.

And — a correction to the earlier draft of this doc — the moat is **not "the package
does your statistics" either.** Statistical *modelling* is inherently per-design and
**cannot be meaningfully abstracted** into a general library (§4b); it correctly lives
in each repo's `04_stats/`, not in `mitomorph`. What the package owns is one thin,
design-*independent* slice of that pipeline (the pseudoreplication-correct,
SuperPlot-ready substrate), not the models.

The two **genuine, defensible differentiators** are:

> **(1) The measurement itself.** mitomorph quantifies the mitochondrial **network
> phenotype** — MiNA-lineage skeleton metrics (length, branch/junction/endpoint) +
> Koopman-lineage shape (solidity, form factor) — that the fusion↔fission biology
> requires. napari-mito-hcs measures **shape-index *texture*** (spot/ridge/valley),
> a screening descriptor, **not network topology** — it cannot produce the project's
> primary readouts. See §3b.
>
> **(2) The analysis-ready handoff.** mitomorph delivers **nucleus-assigned per-cell
> network morphometry, collapsed to the biological replicate (pseudoreplication done
> right), factor-decorated and SuperPlot-ready** — the correct input for whatever
> bespoke contrasts the design demands. No imaging tool does this collapse or ships
> SuperPlots; all stop at a per-object/per-FOV spreadsheet. The *modelling* on top is
> per-repo, by design (§4b).

## 2. The capability matrix (condensed)

Two rows are ✔ for mitomorph and ✘ for **every** imaging tool. Every other row is
matched by at least one competitor.

| Capability | mitomorph | napari-mito-hcs | CellProfiler | Nellie | cellpose | empanada |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| Adaptive mito/cell seg | ✔ | ✘ thresh | ~ plugin | ✔ | ✔ cells | ✔ EM |
| Cell assignment (nucleus-seeded) | ✔ | ✔ | ✔ | ✘ | ✘ | ✘ |
| Mito **network** (skeleton graph) metrics | **✔** | **✘ texture** | ~ | ✔ best | ✘ | ✘ |
| **Pseudoreplication-correct replicate collapse (substrate)** | **✔** | **✘** | **✘** | **✘** | **✘** | **✘** |
| **Transparent per-cell SuperPlots** | **✔** | **✘** | **✘** | **✘** | **✘** | **✘** |
| Composable / reorderable steps | ✔ | ✔ | ✔ | ~ | – | – |
| napari step-observability | ~ vision | ✔ | ✘ | ✔ | ~ | ~ |
| Provenance + env-lock | ✔ | ~ | ~ | ✘ | ✘ | ✘ |
| 2D-MIP / anisotropy fit | ✔ | ✔ | ✔ | ✔* | ~ | ✘ |
| Design-specific *modelling* (contrasts) | per-repo | – | – | – | – | – |

*Nellie is metadata/anisotropy-aware but sub-Nyquist-Z still fragments it — see `05`.

**The decisive rows: network metrics (the *right phenotype*, vs HCS's texture), the
replicate-collapsed substrate, and SuperPlots.** Design-specific *modelling* is
deliberately **per-repo** (last row) — not a package capability, by design (§4b).

## 3. The competitors, precisely

- **napari-mito-hcs (the closest threat).** Does nuclei→cell→mito→per-cell
  assignment→per-cell features→per-FOV summary as config-driven napari steps with a
  batch CLI. Stops at a **per-FOV spreadsheet** — no experimental design, no
  replicate-level inference, no SuperPlots. Segmentation is threshold-based; features
  are shape-index *texture*, not skeleton *network*.
- **CellProfiler (the strongest "it already exists").** `RelateObjects`
  (nucleus→cell→organelle), `MeasureObjectSkeleton` (per-cell branch/trunk/endpoint
  network metrics), `Groups` (condition grouping), `RunCellpose` plugin. It can, in
  principle, produce nucleus-seeded per-cell mito skeleton metrics grouped by
  condition — but it **exports measurements; it does not fit models or test
  factorial/dose contrasts**. Grouping ≠ inference. No napari; `.cppipe` is versioned
  but not lockfile-pinned/git-sha reproducible.
- **Nellie.** Best-in-class mito network features, napari, metadata-aware — but **no
  cell assignment, no design stats.** And on our anisotropic sub-Nyquist-Z data it
  fragments (`05`), so it is not even the right *segmenter* here.
- **cellpose / empanada:** masks only.

## 3b. The measurement differentiator — network topology vs texture

This is the deciding gap, and it comes straight from the project's own measurement
rationale (`05_quantify/README.md`, `METHODS.md`, metrics-provenance research).

**The biological question** (METHODS.md, verbatim): *how innate-immune activation
(cGAS–STING–TBK1) reshapes the mitochondrial **network*** in ISD90-activated BMDCs.
The phenotype is **network remodeling: fusion ↔ fission / fragmentation.** That
dictates the metric family:

- **Skeleton / network family — MiNA lineage** (Valente 2017; Fiji `AnalyzeSkeleton`):
  `skeleton_length_um` (network extent), `branch_count` (complexity/fragmentation),
  junctions, endpoints, mean branch length/diameter. These map *directly and
  interpretably* onto fusion/fission (fused = long skeleton, few branches,
  filamentous; fragmented = short, many disconnected pieces, round).
- **Object-shape family — Koopman/MitoAnalyzer lineage** (Koopman 2005/6; Dagda 2009):
  `mito_area(_fraction)`, `solidity` (compact↔ramified), `form_factor` (round↔
  filamentous). *(Documented caveat: mitomorph's `form_factor` is mathematically
  circularity — bounds correct, name inverted from convention; report as circularity.)*

These are the **field-standard readouts of mitochondrial dynamics** and are what a
mitochondrial-morphology reviewer expects.

**What napari-mito-hcs measures instead:** shape-index **texture** — spot / hole /
ridge / valley / saddle proportions and ratios (SpotRidgeRatio, etc.). That descriptor
was built for high-content **screening** (reduce each cell to a texture fingerprint to
rank thousands of wells). It *correlates* with fragmentation but is **not** network
topology: no branch count, no junctions, no skeleton length. It **cannot produce
the project's primary readouts** without bolting the entire skeleton-network family on top.

**Verdict:** the "~" this doc's matrix gave HCS on network metrics is really a **✘ for
the project's question.** HCS looks similar at the *workflow* level (nuclei→cell→mito→per-
cell→batch, composable, napari) but measures a *different phenotype* at the
*measurement* level. This is a real, biology-driven capability gap, not cosmetic.

## 4. Honest classification

**(iii) declarative/integration synthesis with a hard (i) core — now resting on the
measurement (§3b), not on owning the statistics.**

- The whole-loop "composable synthesis of existing pieces" claim is **integration
  convenience** — no single tool does *this exact* 2D-MIP→Cellpose-cell→per-cell
  skeleton loop, but that is "nobody assembled it," not "it can't be assembled."
- The **genuine, defensible core** is now (1) the **network-phenotype measurement**
  HCS doesn't compute (§3b), and (2) the **analysis-ready handoff** — nucleus-assigned
  per-cell metrics collapsed to the biological replicate (pseudoreplication done
  right), factor-decorated, SuperPlot-ready. Both exist nowhere else in the
  imaging-tool landscape (tools stop at per-object/per-FOV spreadsheets).

**Falsifier, restated:** strip *both* the network-metric family *and* the
replicate-collapsed SuperPlot substrate, and mitomorph is a convenience re-assembly of
HCS + Cellpose. The network metrics + the analysis-ready collapse are what make it more
than a wrapper. *(Note: the earlier draft located the moat in "the statistics"; §4b
corrects that — the reusable slice is the substrate, not the model.)*

## 4b. Does condition-aware statistics belong in the package? — No; the substrate does

The challenge raised (2026-07-18): statistical *modelling* has so many design-specific
degrees of freedom that it can't be abstractably or maintainably owned by a general
library. **This is correct, and the current architecture already honors it.**

**Evidence it isn't abstractable:** the two experiments need *entirely different*
modelling stacks — TBK1i a 2×2 factorial `makeContrasts`+`emmeans`+BY-FDR; glutamine a
six-model dose-response bake-off (linear / log-dose / ordered-poly / each-vs-0 /
Jonckheere / drc). The model *is* the design; the design is per-project. A general
"stats module" would be a leaky god-object or a toy. Accordingly, the R contrasts
**already live in each repo's `04_stats/`, not in `mitomorph`** — same seam as the
experiment-design metadata (`DesignFn` / per-repo `design.py`, `concern-04 §5`).

**But one thin, design-*independent* slice IS abstractable — and it's the slice
biologists most often get wrong. Keep it in the package:**

1. **Pseudoreplication collapse** — cell → per-image replicate median with QC gating
   (`objective_median`). The canonical failure mode (treating N cells as N) and the
   whole point of SuperPlots/Lord 2020. Design-independent: the same correct move for a
   factorial or a dose series.
2. **The analysis-ready tidy contract** — `master_cells` (per-cell, for transparent
   dots) + `master_image` (per-replicate, factor-decorated). The correct input any
   downstream model needs.
3. **SuperPlot rendering** — transparent per-cell dots + replicate-median markers.
   Design-independent viz.
4. **Provenance / env-lock** around all of the above.

**The line:** the package owns *"produce the correct, QC-gated, pseudoreplication-
collapsed, factor-decorated, SuperPlot-ready tables (and the SuperPlot)"*; it **hands
off** to whatever bespoke per-repo model the design demands. It owns the
**analysis-ready substrate, never the model.** In pipeline terms: `integrate` +
SuperPlot viz stay in the package; `contrasts.R`/`run_stats.R`/`dose_sensitivity.R`
stay per-repo (as they already are).

## 5. Where Nellie fits — as an optional segmentation back-end

Reconciling this with `00-SYNTHESIS §2` (which reads too pro-Nellie): Nellie is a
**candidate registered segmentation/feature back-end behind the `MITO_METHODS`
seam — not the default**, for three evidence-based reasons:

1. **Wrong for our data now.** On 9:1 / 14:1 anisotropic sub-Nyquist-Z stacks it
   fragments (`05`); tuned classical 2D-MIP thresholding wins on the data we have.
2. **Doesn't touch the moat.** Nellie stops at features; mitomorph's value is the
   *right network phenotype* + the analysis-ready replicate-collapsed substrate —
   orthogonal to which segmenter feeds it.
3. **Interop, not dependency.** Borrow Nellie's *data model* (metadata/scale-aware,
   hierarchical, provenance-carrying) and offer it as a plug-in method for future
   **isotropic** data (finer z-step / Airyscan / deconvolution), where its 3D
   machinery would finally pay off. Until then it is a documented option, not the road.

The registry architecture (`concern-04 §2`) makes this a clean, non-committal
statement: Li is the realized default; cellpose/Otsu/Nellie/MitoNet are slots. The
package's identity is the **observable, provenance-locked, design-aware
quantification pipeline** — segmentation is swappable underneath it.

## 6. One-line elevator pitch (for grants/README)

> *mitomorph turns fluorescence z-stacks into **analysis-ready mitochondrial-network
> measurements**: nucleus-assigned per-cell skeleton + shape metrics on a 2D-MIP
> substrate, collapsed to the biological replicate (pseudoreplication done right),
> factor-decorated and SuperPlot-ready — composable, napari-observable, with swappable
> segmentation back-ends (mito/cell, per channel availability) and expert GUI mask
> overrides logged transparently. Reproducible and provenance-locked end to end.
> It hands off to whatever bespoke statistical model the design demands; the
> measurement and the honest, replicate-level handoff are the point.*
