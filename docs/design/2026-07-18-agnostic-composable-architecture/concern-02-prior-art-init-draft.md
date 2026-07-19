# Concern 02 — Prior Art in `init-draft`: Why `mitomorph` Was Re-implemented

> **Scope**: A historian's read of the reference packages collected under
> `init-draft/01_modules/`, the MVP he prototyped in `init-draft/02_analysis`,
> and the design docs in `init-draft/docs/Design/`. Read-only survey; no code run.
> **Question**: What pain point in each prior tool pushed the project to write its own
> package, and which design DNA is worth carrying into a generic
> composable + observable (napari-integrated) redesign.
>
> **Source of truth**: This doc condenses and re-frames
> `init-draft/docs/Design/01_field_landscape.md` (13-tool survey) and
> `02_architecture_decisions.md` (the "scanpy-for-mito" design space), read
> against what the MVP scripts actually computed.

---

## 0. The two things actually built in `init-draft`

Before judging the reference tools, it matters what was *proved* with them:

- **`02_analysis/01_raw_data_qc.py`** — a raw-data QA/QC pass on the ISD903
  confocal z-stack (per-slice intensity stats, z-montage, histograms, SNR
  estimate, saturation check). This is the "physical check first" principle
  encoded: characterize the signal before segmenting it.
- **`02_analysis/02_mito_threshold_explore.py`** — a **Python re-implementation
  of the MiNA + Mitochondria Analyzer pipeline** as ~25 composable
  `uint8 → uint8` step functions (contrast stretch, top-hat background
  subtract, denoise, CLAHE, gamma, Frangi/Sato/Meijering tubeness; local
  thresholds mean/Gaussian/Sauvola/Niblack/Phansalkar; cleanup, gap-reconnect,
  aggressive despeckle; Zhang-Suen skeletonize; skeleton metrics +
  regionprops). Four named pipelines (A MitoAnalyzer-replica, B Sauvola,
  C Frangi, D Phansalkar) plus C-value and Frangi-σ parameter sweeps.
- A later **Nellie MVP** (`02_analysis/nellie_mcp_log/`, `nellie_EDA/`) driving
  Nellie via napari-MCP on the same C2 channel, with a perinuclear/peripheral
  two-zone segmentation experiment and a full reproducibility recipe.

**What it measured** (MiNA feature set, ported): branch / junction / endpoint
counts, total & mean branch length, EDT-based mean diameter, and per-object
regionprops (area, perimeter, aspect ratio, form factor). This is essentially
MiNA's 9-parameter network readout + MitoAnalyzer's per-mitochondrion shape
descriptors, in Python.

**What went wrong (the MVP verdict)** — from
`03_results/03_threshold_tuned/README.md`: v1 pipelines **over-segmented**;
v2 re-tuning (proper multiplicative Sauvola/Phansalkar `k`, background-floor
masking, larger `min_size`) helped but the honest conclusion was that
**2D MIP + classical thresholding is a dead end for this data** — MIP discards
axial information and manual per-cell ROI selection does not scale. The README's
own recommendation is to move to **3D DL segmentation (Nellie / empanada /
MitoSegNet)** with manual ground-truth validation (Dice > 0.80). That verdict is
the hinge on which `mitomorph` turned.

---

## 1. Reference packages — one subsection each

Each entry: platform · what it computes · design approach · maintenance signal
(submodule's own last upstream commit) · **the pain point for the project**.

### 1.1 MiNA (Mitochondrial Network Analysis)
- **Platform**: Jython / Java, runs inside Fiji/ImageJ. GUI-only.
- **Computes**: Skeleton-based network morphology — "individuals" vs "networks",
  branch length stats, network size, footprint (9 parameters). Chains Fiji
  plugins (CLAHE, unsharp, median, auto-threshold, Skeletonize, AnalyzeSkeleton,
  Ridge Detection).
- **Design**: Macro that delegates to a small `src/mina/` library
  (`filters`, `statistics`, `tables`, `mina_view`). Metadata smuggled through
  free-text ImageJ comment field as `key=value` pairs. Pre/post-processor hooks.
- **Maintenance**: local checkout dated 2026-03; upstream is the classic 2017
  Stuart-Lab tool, known compatibility breakage on modern Fiji.
- **Pain point**: Fiji-locked, **no Python API**, 2D-only for real use, manual
  single-cell ROI selection (does not scale), no per-mitochondrion output, no
  saved parameters/provenance. This is *the* tool wanted in Python — its
  absence in Python is gap 4.1 of the landscape doc.

### 1.2 Mitochondria Analyzer (MitoAnalyzer)
- **Platform**: ImageJ macro (`.ijm`) + Java entry class. GUI-driven.
- **Computes**: 2D/3D/4D morphology; **adaptive local thresholding** (the key
  idea for uneven illumination); per-cell *and* per-mitochondrion shape
  descriptors (aspect ratio, form factor, solidity); morpho-functional dual-probe
  mode; threshold-optimize commands.
- **Design**: A pile of independent `.ijm` scripts (`2DAnalysis`, `3DAnalysis`,
  `2DThresholdOptimize`, …) registered via a Java plugin. Depends on MorphoLibJ,
  3D ImageJ Suite, and a platform-specific OpenCV adaptive-threshold plugin.
- **Maintenance**: upstream last commit **2022-12** — effectively frozen.
- **Pain point**: Same platform trap as MiNA — no Python API, no programmatic
  data model, macro language is hard to compose/extend, OpenCV dependency splits
  builds per-OS. Its *ideas* (adaptive threshold, threshold-optimize,
  per-object shape) are excellent; its *substrate* is unusable in a Python DS
  stack. Pipeline A of the MVP is a direct replica of its despeckle + outlier
  cleanup.

### 1.3 napari-mito-hcs (Denali Therapeutics)
- **Platform**: Python, napari plugin **+ CLI** (`mito-hcs-batch`).
- **Computes**: Full HCS pipeline — nuclei → cell → mito segmentation, shape-index
  texture features (spot/ridge), per-cell stats.
- **Design**: The **cleanest architecture in the field** for this project's purposes.
  Every configurable class inherits `Configurable` → serializes to/from **TOML**;
  factory defaults (`SegmentationPipeline.load_default('nuclei')`); clean split
  of `segmentation` / `feature` / `pipeline` / `stats` / `finder` / `widget`;
  real test suite. Interactive tuning in napari → export TOML → batch via CLI.
- **Maintenance**: upstream 2025-03; actively engineered, well-tested.
- **Pain point**: 2D plate-imaging only (no 3D confocal z-stacks); segmentation
  is **simple global threshold** (no adaptive, no DL); the pipeline is **rigid**
  (nuclei→cell→mito→features→stats, hard to reorder). Great *pattern*, wrong
  *domain fit* and too inflexible to be reused directly.

### 1.4 MitoClass (LARIS / Univ. Angers)
- **Platform**: Python, napari plugin, Keras/TensorFlow.
- **Computes**: **Classification, not segmentation** — labels 2D MIP patches as
  connected / fragmented / intermediate; heatmap overlay + Plotly class-proportion
  plot.
- **Design**: Patch-wise inference (`_pretreat`/`_processor`/`_utils`/`_widget`),
  pixelwise best-score aggregation. Ships a `.h5` model with no training recipe.
- **Maintenance**: upstream 2025-09.
- **Pain point**: Gives **no measurements** — no per-mito morphology, no network
  metrics, no masks. A coarse 3-class label is not a quantitative phenotype. Also
  2D-only and an opaque pre-trained model. Wrong output type for a morphometry
  pipeline.

### 1.5 empanada / empanada-napari (MitoNet)
- **Platform**: Python library (`empanada`) + napari plugin. Panoptic-DeepLab.
- **Computes**: Instance segmentation of mitochondria in **electron microscopy**;
  2D/3D inference, ortho-plane consensus, proofreading, few-shot fine-tuning.
- **Design**: **Exemplary core/plugin separation** — library runs headless on
  HPC; plugin adds interactive proofreading; model-registry pattern
  (MitoNet/NucleoNet/DropNet).
- **Maintenance**: upstream 2026-03; healthy, generalist model.
- **Pain point**: **EM-only** — does not transfer to fluorescence confocal;
  struggles with tightly-packed mito / MOAS; 2D-to-3D stacking adds split/merge
  errors. The architecture is a role model; the model itself is off-modality.

### 1.6 mitonet-seg (Hoogenboom Group)
- **Platform**: Python; thin wrapper script around `empanada`.
- **Computes**: MitoNet inference on volume-EM stored in WebKnossos.
- **Design**: Single `mitonet-inference.py` with hard-coded params + WebKnossos I/O.
- **Maintenance**: upstream 2024-11.
- **Pain point**: Not a library — a one-off integration script, no abstraction,
  no config. Only value here is as *proof* that empanada's core/plugin split
  enables reuse in a third context.

### 1.7 Nellie (Calico, *Nature Methods* 2025)
- **Platform**: Python, napari plugin, `nellie.plugins` entry point.
- **Computes**: Automated, **organelle-agnostic** segmentation + tracking +
  **hierarchical feature extraction** at voxel / node / branch / organelle
  levels; 2D & 3D, static & time-lapse.
- **Design**: **Metadata-adaptive multiscale Frangi** (auto-picks σ from voxel
  size — kills manual tuning); hierarchical decomposition with adjacency maps /
  k-d trees; mocap-marker tracking + flow interpolation (robust to fission/fusion).
- **Maintenance**: upstream 2026-01; state-of-the-art, active.
- **Pain point**: **Fully automated — no interactive parameter tuning**, which is
  exactly wrong for exploratory validation on a new dataset; feature space is
  huge and un-curated; classical (non-DL) core can underperform on specific
  organelles; API gotcha the MVP hit head-on — `nellie.run()` does **not**
  forward tuning params, you must call pipeline stages directly. Strongest ideas
  in the field, but a closed box that resists the "tune-then-batch" workflow.

### 1.8 MitAnZ-projections (Augustine et al.)
- **Platform**: Python script launching napari + napari-sam (SAM).
- **Computes**: Interactive SAM-assisted segmentation of cells & mito from
  z-projections; regionprops + custom shape + skeleton (branch/junction/length).
- **Design**: Single `fluorescence_analyzer_enhanced.py`; user draws
  `cell_masks`/`object_masks` layers, presses Enter, features dump to CSV.
- **Maintenance**: upstream 2025-05.
- **Pain point**: **Fully manual per-object segmentation** — does not batch;
  exact layer-name magic strings; no config. Good feature set, no automation, no
  provenance.

### 1.9 Trace-Ridges (Arafat et al.)
- **Platform**: **MATLAB** function.
- **Computes**: Automatic ridge/fibre tracing (watershed + edge detection);
  orientation, curvature, gap area, circularity, eccentricity, aspect ratio.
- **Design**: `Trace_Ridges(image, canny_size, gap_threshold)` → 3 structs.
- **Maintenance**: upstream 2025-06.
- **Pain point**: **MATLAB-only** (invisible to a Python stack), built for ECM
  fibres not organelles, no GUI, no batch framework. Interesting tracing metrics,
  wrong ecosystem.

### 1.10 3D-Mito (Augustine et al.)
- **Platform**: Python **Jupyter notebook**.
- **Computes**: Downstream stats on an existing measurements table —
  descriptives, box plots, Spearman/Pearson + Bonferroni, PCA, K-means (elbow),
  Kruskal-Wallis + Dunn's post-hoc heatmaps.
- **Design**: Edit `file_path`, run cells top-to-bottom.
- **Maintenance**: upstream 2025-05.
- **Pain point**: **Jupyter** (violates the "no Jupyter" rule); no image
  integration; hard-coded analysis, no standardized ingestion. But it correctly
  identifies the *last* pipeline stage — table → hypothesis test — which
  `mitomorph`'s R-in-Docker stats step now owns.

### 1.11 napari-skimage (Witz) — architectural reference
- **Platform**: Python, napari plugin (magicgui).
- **Computes**: Nothing mito-specific — wraps scikit-image filter / threshold /
  morphology / restoration / detection / label / regionprops as per-category
  widgets.
- **Design**: One widget module per operation category; magicgui auto-generates
  UI from function signatures; deliberately "simple for beginners".
- **Maintenance**: upstream 2025-09.
- **Pain point**: **No pipeline composition, no batch, no config, no
  multi-channel** — atomic operations only. Confirms the widget-per-step idiom is
  cheap, but a pile of widgets is not a pipeline.

---

## 2. Synthesis

### 2.1 What the MVP replicated
The `02_mito_threshold_explore.py` MVP is a faithful **Python port of
MiNA + MitoAnalyzer**: MiNA's skeleton-network paradigm (branch/junction counts,
network decomposition, footprint) and MitoAnalyzer's adaptive local thresholding
+ per-mitochondrion shape descriptors + aggressive despeckle cleanup, refactored
into composable `uint8→uint8` step functions with named pipelines and parameter
sweeps. The QC script encodes the "characterize signal before segmenting"
discipline. The Nellie MVP tested the automated 3D alternative on the same data.

### 2.2 The concrete gap that justified a fresh package
No single existing tool offered **all** of: (1) classical MiNA/MitoAnalyzer
morphometry, (2) *in Python*, (3) **composable** step functions you can reorder,
(4) **interactive tune-then-batch** parameter selection, (5) a shared,
serializable **container with provenance**, and (6) a path to **3D / DL**
segmentation. The field splits cleanly into "right features, wrong platform"
(MiNA, MitoAnalyzer, Trace-Ridges), "right platform, wrong output" (MitoClass),
"right architecture, wrong domain" (napari-mito-hcs, empanada), and "right ideas,
no tuning" (Nellie). The MVP itself demonstrated the deeper gap: **2D-MIP
classical thresholding over-segments and doesn't scale**, so the redesign needed
a substrate that carries 3D data and swappable (classical ↔ DL) segmentation
behind one composable, observable API — the "scanpy-for-mito" target set out in
`02_architecture_decisions.md`.

### 2.3 Design DNA — keep vs discard

| Source | Idea | Verdict | Why |
|--------|------|---------|-----|
| MiNA | Skeleton network metrics (branch/junction/endpoint, network vs individual, footprint) | **Keep** | The canonical, validated fluorescence readout; already ported in MVP |
| MiNA | Metadata as free-text `key=value` in an image comment | Discard | Fragile; replace with structured container / `experiment.toml` |
| MitoAnalyzer | Adaptive local thresholding (Sauvola/Phansalkar) for uneven illumination | **Keep** | Core to fluorescence segmentation quality |
| MitoAnalyzer | Threshold-Optimize (systematic param search) | **Keep** | Becomes the interactive napari tune-then-batch loop |
| MitoAnalyzer | Per-mitochondrion shape descriptors (aspect ratio, form factor, solidity) | **Keep** | Complements network metrics at object level |
| MitoAnalyzer / MiNA | ImageJ-macro substrate, GUI-only, per-OS OpenCV builds | Discard | The reason for the whole re-implementation |
| napari-mito-hcs | `Configurable` → **TOML** serialization; factory defaults; tune-in-napari → CLI batch | **Keep (central)** | Exactly the target reproducibility + tune-then-batch model |
| napari-mito-hcs | Clean module seam (segment / feature / stats / finder / widget) + real tests | **Keep** | Template for `mitomorph` module boundaries |
| napari-mito-hcs | Rigid fixed pipeline order | Discard | Need reorderable composable steps instead |
| empanada | **Core/plugin separation** — headless library + interactive plugin; model registry | **Keep (central)** | Enables HPC batch + napari proofreading from one core; swappable models |
| empanada | EM-only MitoNet weights | Discard (as model) | Off-modality for confocal fluorescence |
| Nellie | **Metadata-adaptive multiscale Frangi** (auto-σ from voxel size) | **Keep** | Removes the worst manual-tuning pain |
| Nellie | **Hierarchical decomposition** (voxel→node→branch→organelle) as the data model | **Keep** | The strongest observability idea; maps to biological scales |
| Nellie | mocap + flow tracking (fission/fusion-robust) | Keep (later) | For time-lapse; not needed for static v1 |
| Nellie | Fully automated, no tuning; `run()` swallows params | Discard | Must expose stages + interactive tuning |
| napari-skimage | magicgui widget-per-step idiom | Keep (thin) | Cheap way to surface steps; but wrap in a real pipeline |
| 3D-Mito | Rigorous downstream stats (non-parametric + multiple-testing correction, PCA/clustering) | **Keep** | The final table→hypothesis stage (now R-in-Docker) |
| 3D-Mito / MitoClass | Jupyter notebook / opaque `.h5` model | Discard | Violates env-lock + provenance rules |
| MitAnZ / SAM | Interactive SAM-assisted segmentation | Keep (optional) | Useful HITL seeding path; must be batchable + provenanced |
| squidpy (ref) | Container-centric, submodule API (`.pp/.tl/.pl`), in-place mutation + `.uns` history | **Keep (central)** | The overall API skeleton being targeted |

---

## 3. Loose ends worth noting
- The submodule "last-commit" dates above are the upstream repos' own histories
  (MitoAnalyzer 2022 is genuinely frozen; empanada/Nellie/MiNA are current). They
  are maintenance signals, not local checkout dates.
- The MVP's `PIXEL_SIZE_UM` was hard-coded to `0.1` while CZI metadata says
  `0.0852` µm/px — a provenance bug the container/`experiment.toml` design must
  prevent by carrying pixel size from acquisition metadata.
- `docs/Design/02_architecture_decisions.md` already lays out the three API
  shapes (functional scanpy-style / OOP pipeline / hybrid) and container options
  (AnnData-like `MitoData` / squidpy `ImageContainer` / dict). This concern doc
  is the "why", that doc is the "what to build".
