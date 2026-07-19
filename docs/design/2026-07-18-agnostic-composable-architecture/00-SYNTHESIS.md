# mitomorph — Agnostic, Composable, Observable Architecture (synthesis)

**Date:** 2026-07-18
**Author:** Claude (Opus 4.8), orchestrator — synthesis of a 4-agent parallel exploration.
**Status:** EXPLORATION ONLY. 
**No source code was changed.** 
This is a design
substrate for discussion before any refactor is authorized.
**Concern docs synthesized here:**
- `concern-01-blast-radius-coupling.md` — where the couplings physically live (file:line).
- `concern-02-prior-art-init-draft.md` — why the re-implementation happened (13-tool survey + the MVP verdict).
- `concern-03-landscape-web-research.md` — external evidence (agy) on the field + modern design patterns.
- `concern-04-target-architecture-seams.md` — the proposed boundaries & interfaces.
- `concern-05-nellie-anisotropy-evidence.md` — is the Nellie over-seg really anisotropy? (evidence-graded).
- `concern-06-moat-and-positioning.md` — what justifies mitomorph vs Nellie/HCS/CellProfiler.
- `05-data-and-2d-mip-rationale.md` — the data, its limits, why 2D-MIP not 3D (readable reference).
- `07-positioning-and-moat.md` — the grounded moat statement + Nellie-as-optional-backend.
- External artifacts: `mitoISDscopy/docs/llm-research/2026-07-18-mito-package-landscape/{brief_1,report_1}`.

---

## 1. The Vision
`mitomorph` should stop being *"the compute engine for these two specific experiments read
from CZI in this folder tree"* and become a **generic library** with clean seams:

1. **Experiment-agnostic** — no experiment ids, no group-label grammar in shared code.
2. **Vendor-agnostic** — default input is open **OME-NGFF (OME-Zarr)**; CZI/Zeiss
   conversion is outsourced to the user / per-repo ("CZI just happened to be the Zeiss output").
3. **Composable** — mito-mask *methods* (Li today; Otsu/ML/Nellie later), cell-mask
   *recipes* (mito+nucleus today; membrane+nucleus later), and *metrics* are all pluggable;
   the user runs operations in any order/params the science requires.
4. **Observable** — every single operation is openable in napari; the *same* op stack runs
   headless (CLI/Python) or interactively, producing a clean step-by-step "sequence of events"
   figure for papers/grants. Each non-obvious metric (form_factor, solidity, branch_count,
   mean_diameter_um) carries a napari visual of *what it measured*.
5. **Idempotent & layout-blind** — the user says where outputs land; the library does not
   assume a repo structure.
6. **Expert-in-the-loop** — a hand-drawn mask in napari may override the mask it overlaps,
   but the composed op stack + params + the override event are **logged** (provenance);
   the original mask is never destroyed.
7. **2D now, 3D deferred but not designed out.**

---

## 2. Why this is the right move (prior-art + landscape convergence)

Both the local prior-art read (concern-02) and the external web research (concern-03)
land on the **same** conclusion, independently:

- The classic morphometry tools conceptually were — **MiNA** and
  **Mitochondria Analyzer** — both are **Fiji/ImageJ-macro-locked and dormant (~2021-22)**:
  no Python API, GUI-coupled (breaks headless), no provenance, macro language can't hold
  graph data or import the scientific-Python stack. That platform trap *is* the reason for
  a re-implementation. Their **ideas** (skeleton-network metrics; adaptive local thresholding;
  threshold-optimize; per-object shape) are worth keeping; their **substrate** is barely usable.
- **Nellie** (Python + napari, *Nature Methods* 2025) is the standout *feature* engine
  (metadata-adaptive, scale-aware, hierarchical pixel→branch→organelle→cell). Borrow its
  **data model**, but it is an **optional segmentation back-end, not the default** — on our
  anisotropic sub-Nyquist-Z data it fragments badly (see `05` / `concern-05`), and it does not
  touch the moat (it stops at features). Interop target, not the road.
- **cellpose/-SAM** and **empanada/MitoNet** are healthy but **segmenters only** — they slot
  in as *pluggable back-ends*, not as the package identity.
- The best architectural template in the field is **napari-mito-hcs** (`Configurable`→TOML,
  tune-in-napari→CLI-batch, clean module seams, real tests) — except it's rigidly ordered and
  2D-plate-only, threshold-seg, and stops at a per-FOV spreadsheet. **empanada's core/plugin
  split** is the role model for headless-library + interactive-plugin from one core.
- The scope-fixing negative result from `init-draft`: on anisotropic z-stacks the **3D route
  (Nellie) fragments structurally** (data-limited, not a tool bug — `05`), while **classical
  2D-MIP thresholding over-segmented only at first and was *tunable* back to scale** (a
  parameter fix, v2). Plus the `run()`-swallows-params Nellie gotcha → *exposed, tunable stages*
  matter. So the package's identity is **not** "our segmenter" — it is **the observable,
  provenance-logged, composable, design-aware quantification pipeline on an OME-Zarr
  substrate**, with segmentation swappable underneath (see the moat, §2b).

**External design patterns adopted (concern-03 §3):** API-first / GUI-optional (never import
napari in compute — already our house rule); OME-Zarr substrate with masks in `/labels`
(nearest-neighbor pyramids, `image-label` colors + `source` linkage, µm in axes +
`coordinateTransformations`); step-wise inspectable layers; HITL correction with hashed,
event-logged provenance; DAG/recipe re-execution; pluggable segmentation; a per-run
provenance file. These corroborate the direction rather than redirect it.

### 2b. The moat — two differentiators, and statistics stays OUT (grounded; full detail `07`)

Read honestly against the field, the moat is **narrow but real**, and it is *not*
where intuition points. Segmentation, per-cell morphometry, and even "composable +
napari-observable synthesis" are **already shipped** — the last by **napari-mito-hcs
(Denali, SLAS Discovery 2025)**. So the differentiators are two, both evidence-backed:

1. **The measurement.** mitomorph quantifies the mitochondrial **network phenotype**
   (MiNA skeleton: length/branch/junction; Koopman shape: solidity/form-factor) that
   the fusion↔fission biology requires. **napari-mito-hcs measures shape-index
   *texture*** (spot/ridge/valley) — a screening descriptor, **not network topology**;
   it cannot produce the primary readouts. *(Source: `05_quantify/README.md`,
   `METHODS.md`, `docs/llm-research/2026-07-18-mito-metrics-provenance/`.)*
2. **The analysis-ready handoff.** Nucleus-assigned per-cell metrics **collapsed to the
   biological replicate (pseudoreplication done right)**, factor-decorated and
   SuperPlot-ready. No imaging tool does this collapse or ships SuperPlots.

**Statistics-the-*modelling* is deliberately NOT in the package.** It is inherently
per-design (TBK1i's 2×2 `makeContrasts` vs glutamine's 6-model dose bake-off — two
bespoke stacks), un-abstractable, and already lives in each repo's `04_stats/` — the
same seam as experiment-design metadata. The package owns only the design-*independent*
substrate: pseudoreplication collapse + tidy `master_cells`/`master_image` contract +
SuperPlot rendering + provenance. It owns the **analysis-ready substrate, never the
model** (`07 §4b`).

Falsifier (kept explicit): strip *both* the network-metric family *and* the
replicate-collapsed SuperPlot substrate and mitomorph is a convenience re-assembly of
napari-mito-hcs + Cellpose. Those two — plus swappable seg back-ends and logged expert
overrides — are what make it more than a wrapper.

---

## 3. The blast radius (what actually has to move)

Concern-01 is the good news: **the coupling is concentrated and cleanly extractable.**
`segment.py`, `metrics.py`, `viz_common.py`, and most of `omezarr.py`'s primitives are
**already library-grade** (arrays in, arrays/rows out). Everything blocking agnosticism sits
in four pockets:

| Pocket | Where | Severity | The seam |
|---|---|---|---|
| **Experiment design logic** | `integrate.py:28` `KNOWN_EXPERIMENTS`; `:31-49` `_design_columns` (tbk1i/isd90/gln_mM string-parsing); `:52-71` id-keyed warning | **BLOCKS** | integrate becomes reshape-only + a user-supplied `DesignFn`; tokens already appear *only* here + in per-repo `experiment.toml` |
| **Vendor/CZI in core** | `io.py:9,42-51` (bioio `channel_mip` — the read primitive every stage funnels through); `inventory.py`, `extract.py`, `omezarr.czi_to_zarr`; literal `.czi` threaded through ~10 files as manifest keys/paths | **BLOCKS** | `ImageReader`/`Substrate` protocol; core consumes NGFF only; CZI→NGFF moves to `adapters/czi.py` (optional extra) |
| **Fixed on-disk layout** | `NN_stage/<group>/<sample>` string literals in every `run`/`compute_one`; `config.py:66-75` roots; `integrate._collect` hardcodes `*/*/tables/` glob depth | **HIGH** | user-specified out-dir; layout resolver with the current tree as the *default*; compute returns arrays/rows, orchestration writes |
| **Config shape** | `config.py:104-110` walks parents for the umbrella `run-config.yaml`; `experiment.toml`-in-cwd is the *only* entry; two-channel-index identity (`mito_/nucleus_channel_index`) | **HIGH** | accept explicit config path/dict; named **channel-role map** (arbitrary channels); optional sample parsing |

Lower-severity: `.npy` mask filenames as the A→B contract (→ NGFF `labels/`, already modeled);
mandatory `parse_sample` regex (→ optional); NGFF version pinned in two places (→ single source);
`size_z==1`/MIP assumption (→ make MIP-collapse just the first op).

The main *structural* entanglement to unpick: the compute stages
(`mitoseg`/`cellseg`/`quantify`) have **pure kernels** but their `compute_one` wrappers mix in
manifest reads, fixed paths, `mkdir`, and `write_csv`. Splitting kernel (pure) from
orchestration (I/O) is the recurring move.

---

## 4. Target architecture (the seams)

Full detail in concern-04. The shape:

**Four rings, dependencies inward only.**
```
per-repo glue (experiment.toml · design.py · CZI→NGFF choice · out-dir · recipe+params)
        │ depends on
orchestration (RecipeSpec runner · provenance log · out-dir)  ◄──  napari-plugin (same Ops, interactive)
        │                         both depend on                        │
        └──────────────►  CORE (pure)  ◄──────────────────────────────┘
                          registries: mito-methods · cell-recipes · metrics
                          Operation protocol · OpResult(array+view) · Substrate(NGFF) contract
                          segment.py / metrics.py — unchanged math, no I/O
                                    ▲ depends on
                          io/adapters (vendor→NGFF), OUTSIDE core  —  delete it, core still runs on any NGFF store
```

**Three registries, one pattern** (decorator-registry keyed by name; core enumerates nothing):
- **`MITO_METHODS`** — `(mip, px_um, params) -> bool mask`; Li is the first registrant.
- **`CELL_RECIPES`** — declare **required channel roles** (`("mito","nucleus")` vs
  `("membrane","nucleus")`); per-repo config maps role→channel index, so the library is blind
  to a scope's layout.
- **`METRICS`** — each metric is **self-describing** (name, level `per_object|per_cell|per_image`,
  units, a transparent `describe` of *what & how*) + an optional **`view`** returning napari
  layers that illustrate the measurement (form_factor draws the object + its equal-area circle,
  etc.). `metrics=["form_factor"]` runs one; "all" = `METRICS.names()`. The registry *is* the
  transparent catalogue; a user registers their own metric from their repo.

**Composable operation model — identical headless & interactive.** The pipeline is a
**stack of typed `Operation`s** (not the `cli.py` if/elif ladder). Each reads/writes named
substrate layers and returns `OpResult{arrays, rows, view}`. A **`RecipeSpec`** (ordered ops +
params, a YAML) is the reproducible object — *that* is how "any order, any params" becomes
first-class instead of CLI-arg order. Headless runner and napari plugin consume the **same**
`RecipeSpec`; each op's `ViewSpec` is one captioned frame of the paper-figure "sequence of
events". A GUI exploration session **serializes back** to a `RecipeSpec` + override events —
the interactive session *is* a reproducible recipe.

**OME-NGFF is the working substrate (not a parallel export).** Image pyramid read-only under
`0/`; ops write `labels/` groups each carrying `mitomorph_provenance` (method/params/version/
git_sha); `.zattrs.mitomorph_log` is the on-store op stack. **Expert override:** a hand-drawn
`labels/cells_expert` is written; the runner detects overlap (IoU) with its target and records
an **override event** `{override_of, iou, n_objects_changed, drawn_by, ts}` into both the group's
provenance and the op log — downstream ops consume the expert mask where it overlaps; the
original is immutable. This *extends* the `omezarr.py` label discipline already in place.

**Experiment metadata moves OUT.** `integrate` → reshape-only + a `DesignFn` callback; a per-repo
`design.py` supplies `group→factors` and the sample-name parsing that used to be `parse_sample`;
`experiment.toml` gains a `design = "design:design"` pointer that orchestration imports. The
library never enumerates experiments.

---

## 5. Phased migration (mostly non-breaking; nothing done yet)

From concern-04 §6 — each phase is independently reviewable, the two live experiment repos keep
working until they opt in:

1. **Wrap existing fns as registrants** (non-breaking) — `pipeline_li`, `build_territories_v2`,
   `metrics.py` become `MitoMethod`/`CellRecipe`/`Metric`; stages call the registry by name;
   identical outputs, seams now exist.
2. **Invert the substrate** (non-breaking read path) — ops read MIP + roles from the `.zarr`
   store (already produced by `to-zarr`) instead of re-reading CZI; `bioio` isolated into
   `adapters/czi.py`, the only importer.
3. **`RecipeSpec` + runner** (additive) — `mito run recipe.yaml` beside current subcommands;
   `all` becomes a default recipe; op log to `.zattrs`.
4. **De-hardcode integrate** (mildly breaking for the two repos) — `DesignFn` + per-repo
   `design.py`; add `design.py` to both repos so master tables are unchanged.
5. **napari plugin over the same ops** (additive) — renders each `ViewSpec`, captures expert
   overrides; the matplotlib `*_viz.py` become the headless `ViewSpec` renderer, then retire.
6. **Out-dir decoupling** (breaking, last) — user-specified `--out`; current tree stays the default.

3D stays deferrable: `_label_axes` already handles 2D/3D, `px_um` generalizes to a scale vector,
ops are ndim-generic, and MIP-collapse is just the first op (`op: mip_collapse`) a 3D recipe omits.

---

## 6. Open decisions 

1. **Scope & appetite** — is this the "scanpy-for-mito" full redesign (registries + RecipeSpec +
   napari plugin), or do we land only the *cheap, high-value* de-couplings first
   (Phases 1–2 + 4: registries, CZI-into-adapter, `DesignFn`) and defer the plugin/recipe engine?
   Recommendation: **Phases 1, 2, 4 first** — they remove all four BLOCKS/HIGH couplings with
   near-zero behavior change and no new subsystem; the RecipeSpec engine + napari plugin (3, 5, 6)
   are a larger, separately-scoped effort.
2. **Nellie: interop vs reimplement** — adopt Nellie as a registered segmentation/metrics back-end
   (call it), or keep our own classical kernels and merely match its data model? Recommendation:
   **pluggable back-end** — our value is the observable pipeline, not the segmenter.
3. **Container object** — the `init-draft/docs/Design/02_architecture_decisions.md` "MitoData vs
   ImageContainer vs dict" question is still open. Recommendation: let the **`.zarr` store +
   `Substrate` protocol** *be* the container (avoids inventing an in-memory type); revisit only if
   ergonomics demand it.
4. **`04_stats/` location (repo-org, your point #1)** — you want statistics under `02_analysis/`
   rather than as a top-level `04_stats/` sibling. This is **analysis-repo organization, not
   mitomorph code**, so it's out of *this* exploration; noted here so it isn't lost. When you're
   ready I can reorganize both experiment repos (`02_analysis/stats/…` with the R+viz split intact,
   `git mv` to preserve history) — say the word.
5. **Channel-role vocabulary** — fix a small controlled set now (`mito`, `nucleus`, `membrane`,
   `cyto`?) so recipes and per-repo configs agree.



