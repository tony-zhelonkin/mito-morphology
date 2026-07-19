# Concern 03 — Landscape & external design evidence (web research)

Synthesis of external web research (agy, 2026-07-18) on the mito-morphology tool
landscape and modern composable-pipeline design. Source artifacts:
`mitoISDscopy/docs/llm-research/2026-07-18-mito-package-landscape/`
(`brief_1` + verbatim `report_1`). Citations below are agy-sourced leads and
should be spot-checked before manuscript use.

## 1. Landscape at a glance

| Tool | Platform | 2025-26 status | Dim | Measures | Composable? | Mask/output | Provenance |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| **MiNA** | Fiji / ImageJ macro | Dormant (~2021) | 2D/3D | footprint, branch length, network branches, complexity | Low — GUI macro | CSV, binary skeleton | None |
| **Mitochondria Analyzer** | Fiji macro + Java | Dormant (~2021) | 2D/3D | count, perimeter, form factor, aspect ratio, junctions | Medium — GUI batch | CSV, ROI, mask TIFF | Manual param files |
| **MitoGraph** | C++ CLI | Dormant (~2024) | 3D | network length, width, branching, volume, graphs | High CLI, no Py API | VTK, CSV, TIFF | CLI flags only |
| **MitoHacker** | Python notebooks → cloud | Archived/proprietary | 2D | footprint, aspect ratio, circularity | Low (OSS basic) | CSV, overlays | Notebook history |
| **Nellie** | Python + napari | **Active** (Nat Methods 2025) | 2D/3D/+t | hierarchical pixel→branch→organelle→cell, tortuosity, tracking | **High** — API + GUI | numpy, CSV, Zarr | **High** — reads scale/metadata |
| **MitoTNT** | Python | Dormant (2023) | 4D | tracking, fission/fusion, velocity, MSD | High, batchable | CSV tracks, gltf | Config files |
| **MicroP** | MATLAB | Abandoned (2011) | 2D | morphological subtypes | Low (license) | MAT, CSV | None |
| **mitometer** | MATLAB | Abandoned (2021, → Nellie) | 2D/3D+t | size, velocity, fission/fusion | Low (license) | MAT, CSV | None |
| **cellpose / -SAM** | Python/PyTorch + napari | **Active** | 2D/3D | generalist masks (no morphometry) | **High** — API/CLI/GUI | labels, numpy, OME-Zarr | model versions/configs |
| **empanada / MitoNet** | Python/PyTorch + napari | **Active** | 2D/3D | panoptic mito masks (EM-trained) | **High** | labels, OME-Zarr, RLE | model cards |

Takeaways:
- The classic **morphometry** tools (MiNA, Mitochondria Analyzer) are dormant and
  Fiji-locked. The classic **MATLAB** tools (MicroP, mitometer) are abandoned.
- The only actively maintained tool that does *both* modern segmentation *and*
  hierarchical morphometry with metadata awareness is **Nellie** (Python + napari).
  It is the closest existing thing to our target and the strongest prior art /
  interop target.
- **cellpose** and **empanada/MitoNet** are actively maintained but are
  segmenters only — they produce masks, not network metrics. They slot in as
  pluggable segmentation back-ends, not as the whole pipeline.

## 2. Why the Fiji/ImageJ-Java mito tools died (and why we leave)

1. **Macro language (IJM) is a dead end** — typeless, no maps/trees for graph
   representations, global-variable scoping, no `try/catch`, cannot import
   scikit-image / PyTorch / pandas.
2. **Java/ImageJ2 lock-in** — driving plugins from Python needs PyImageJ, which
   boots a JVM (memory overhead, `-Xmx` tuning) and pays JNI array-copy costs to
   move masks across the Python/Java boundary.
3. **GUI-coupling defeats headless use** — commands depend on the active window
   (`selectWindow`, `Duplicate...`); `--headless` runs throw null-pointer errors
   because window containers are never instantiated. Processing logic is
   interleaved with `GenericDialog` prompts, so routines can't be isolated.
4. **No provenance / reproducibility** — macros log no version, params, or array
   shapes; binarization drifts with Java version and Fiji update-site state.
5. **Dependency rot + academic abandonment** — plugins pinned to Java 8; update
   sites introduce breaking changes; original authors move on (mitometer's author
   built Nellie to replace it; MitoHacker went closed-cloud).

## 3. Design patterns worth adopting

1. **API-first, GUI-optional.** Core segmentation/skeletonization/quantification
   are pure functions over numpy arrays / Zarr groups returning plain Python types.
   **Do not import napari/magicgui in core modules** — so the engine runs headless
   on HPC / in Docker. napari is a viewer layered on top, not a dependency of compute.
   (Matches our AGENTS.md rule: compute writes tables; viz reads tables.)
2. **OME-Zarr (OME-NGFF) as the substrate.** No custom folders / proprietary
   formats for intermediates. Chunked lazy loading, cloud-native, physical scaling
   preserved. Masks live in the `/labels` group.
3. **Step-wise inspectable layers.** Every stage (raw → denoised → mask →
   skeleton → graph nodes) is a discrete overlayable napari layer, so an expert
   can toggle/opacity-check that a skeleton traces the real signal. This is exactly
   the "prove thresholding visually before batching" MVP rule we already hold.
4. **Human-in-the-loop mask correction with provenance.** napari Labels layer for
   paint/erase/split; listen to `layer.events.data`, record who/what/when, hash the
   corrected mask for auditability.
5. **DAG workflow + caching.** Model raw→preprocess→segment→skeleton→features as a
   DAG (networkx / napari-workflows) so changing an upstream param re-runs only
   downstream steps.
6. **Pluggable segmentation.** Unified interface across classic thresholds
   (Li/Otsu), generalist DL (cellpose), and specialist models (MitoNet) — do not
   hardcode one segmenter.
7. **Autogenerated provenance file per run** — JSON/YAML with pipeline version, git
   commit, env package list, parameters, raw-data hashes, and manual-edit records.

Exemplars: **Nellie** (metadata-driven, parameter-free), **devbio-napari /
napari-workflows** (DAG re-execution), **spotiflow** (Python-first CLI/API with a
thin napari wrapper — the separation model we want).

## 4. OME-NGFF `labels/` best practice (mask storage)

- Masks go in a `labels/` group at the image root; `labels/.zattrs` holds a
  registry `{"labels": ["mitochondria_mask"]}` so viewers can discover them.
- Each mask group carries `image-label` metadata: `colors`
  (`label-value`→`rgba` uint8, background label 0 = `[0,0,0,0]` transparent),
  optional per-label `properties`, and `source: {"image": "../../"}` linking the
  mask back to its parent image.
- Mask **multiscale pyramids must use nearest-neighbor downsampling** (never
  averaging/Gaussian) or you invent non-existent label IDs.
- `axes` declare `type: "space"` + `unit: "micrometer"`; voxel size lives in the
  `coordinateTransformations` `scale` vector, in axis order. This is how physical
  units (our pixel size) ride with the data.

## 5. Implications for `mitomorph`

- **Anchor on Nellie's model** (Python + napari + metadata awareness) as prior art
  and interop target; don't reinvent hierarchical skeleton metrics if Nellie's
  outputs suffice.
- **Adopt OME-Zarr with `/labels`** as the module seam already sketched in the
  cell-seg architecture plan — this concern doc gives the external justification.
- **Keep compute headless and napari-free**, viz as a separate layer — already our
  house rule, now externally corroborated.
- **Segmentation is a pluggable back-end**, not the package's identity; the
  package's value is the observable, provenance-logged quantification pipeline.
