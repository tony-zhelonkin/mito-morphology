# Concern 01 — Blast-radius coupling audit

**Scope:** every file in `src/mitomorph/` (18 modules, ~2360 LOC).
**Goal:** inventory all coupling that blocks `mitomorph` from becoming a generic,
experiment-agnostic, vendor-agnostic, composable library defaulting to OME-NGFF.
**Method:** full read of each module + targeted grep for experiment ids, group
string-parsing, vendor tokens (`.czi`, `bioio`, `BioImage`, `CZYX`), and hardcoded
on-disk layout.

This is a **coupling + seam-location** inventory only. It does not propose the new
architecture (owned by another concern).

---

## Severity-ranked summary

| # | Finding | File:line | Severity | Seam (what user should supply) |
|---|---------|-----------|----------|--------------------------------|
| 1 | `KNOWN_EXPERIMENTS` hardcoded experiment ids | integrate.py:28 | **BLOCKS** | drop; no library-known experiment set |
| 2 | `_design_columns()` hardcodes tbk1i/isd90/gln_mM group-name parsing | integrate.py:31-49 | **BLOCKS** | user-supplied design callback / factor spec in config |
| 3 | `_decorate()` warns keyed on `KNOWN_EXPERIMENTS` | integrate.py:52-71 | **BLOCKS** | generic "unparsed factor" warning, not id-keyed |
| 4 | CZI reading wired directly into core compute/io (not behind a reader seam) | io.py:9,42-51; inventory.py:6-31; extract.py:27; omezarr.py:80-95 | **BLOCKS** | `ImageReader` protocol; core consumes OME-NGFF arrays |
| 5 | `.czi` suffix hardcoded across pipeline (manifest keys, path build, glob) | io.py:27,29; mitoseg.py:45; cellseg.py:117,124; quantify.py:187,196; omezarr.py:297,303; *_viz.py | **BLOCKS** | reader owns extension/format; keys by sample id, not filename |
| 6 | Fixed stage-numbered dir layout (`00_..05_`, `master/`) baked into every module | config.py:66-75; integrate.py:27,76,107; extract.py:32; mitoseg.py:39,54,65; cellseg.py:125,133,158; quantify.py:181,197-244; omezarr.py:301,385; all *_viz.py | **HIGH** | user-supplied output/layout resolver; library takes explicit paths |
| 7 | `run-config.yaml` discovery walks up to umbrella superproject | config.py:104-110,143-156 | **HIGH** | explicit config path/dict; no parent-walk to a specific repo shape |
| 8 | `experiment.toml`-in-cwd is the only config entry (`load_config`) | config.py:127-168; cli.py:61 | **HIGH** | accept a config object/dict; toml is one adapter, not the API |
| 9 | Two-channel-index identity model (`mito_channel_index`/`nucleus_channel_index`) | config.py:33-34,163-164; io.py:42-51; extract.py:24; quantify.py:225; cellseg.py:142 | **HIGH** | named-channel/role map, arbitrary channels |
| 10 | Fixed mask filenames + `.npy` substrate as the A→B contract | quantify.py:207-209; cellseg.py:160-161; mitoseg.py:67; omezarr.py:185-190 | MEDIUM | NGFF `labels/` as the contract (omezarr already models it) |
| 11 | `parse_sample` (group/replicate regex) mandatory for every stage | config.py:77-82; used in all stages | MEDIUM | optional; a sample is an id + free-form metadata dict |
| 12 | `schema_version`/`config_version` assert tied to umbrella run-config | config.py:85,148-154 | MEDIUM | reasonable per-config; decouple from mandatory umbrella file |
| 13 | Compute stages own file I/O + directory creation (not pure functions) | mitoseg.py, cellseg.py, quantify.py `compute_one` | MEDIUM | pure compute already exists in segment.py/metrics.py; orchestration should not be in the package core |
| 14 | `size_z==1` skip + MIP assumption (2D-from-3D-stack workflow) | mitoseg.py:47; cellseg.py:119; quantify.py:189; *_viz.py | LOW-MED | this specific imaging workflow; make a composable step |
| 15 | provenance keys off dist name `mito-morphology` + walks git of package dir | provenance.py:19,34,50-65,80-95 | LOW | fine for this package; only note if repackaged |
| 16 | NGFF version pinned two places (config default + module const) | config.py:59; omezarr.py:36 | COSMETIC | single source |

---

## 1. Experiment-specific hardcoding

The worst offenders all live in `integrate.py`.

- **integrate.py:28** — `KNOWN_EXPERIMENTS = ("20260518_TBK1i", "20251029_glutamine")`.
  A shared library literally naming this project's two experiments. Used at
  :62 and :68 to decide whether to warn. **Leak:** the package knows the caller's
  experiments. **Seam:** delete. The library should have no notion of a "known"
  experiment set; warning logic (#3) should be generic.

- **integrate.py:31-49** — `_design_columns(experiment, group)`. Hardcodes:
  - `experiment == "20260518_TBK1i"` → `tbk1i = int(group.startswith("TBK1i"))`,
    `isd90 = int("ISD90" in group)` (:38-43).
  - `experiment == "20251029_glutamine"` → `gln_mM = float(re.match(r"([0-9.]+)mM", group))` (:44-48).
  This is the single most experiment-specific function in the package: it encodes
  the factorial design of two named experiments and the string grammar of their
  group labels. **Leak:** design semantics (what factors exist, how to derive them
  from a group label) are project knowledge, not library knowledge. **Seam:** the
  user supplies design decoration — either (a) a callback `group -> dict[factor, value]`
  passed into `integrate.run`, (b) a declarative factor spec in per-repo config
  (regex/patterns → column), or (c) integrate emits only `condition`/`replicate`
  and per-repo code adds factor columns downstream. The library should never
  `if experiment == "<id>"`.

- **integrate.py:52-71** — `_decorate()`. The `condition = group` / keep-`replicate`
  logic is generic and fine; the branch at :62-70 keyed on `KNOWN_EXPERIMENTS` is
  the leak. **Seam:** warn generically when a user-supplied decorator yields NaN;
  no id-keyed special-casing.

No other module hardcodes experiment ids or condition strings — the tokens
`TBK1i`/`ISD90`/`GLN`/`DMSO`/`mM` appear **only** in integrate.py and in the
per-repo `experiment.toml` files (correct location). The design coupling is
therefore concentrated and extractable.

---

## 2. Vendor-format (CZI / Zeiss) coupling

CZI reading is **not** isolated behind a seam — it is imported and called directly
in four core modules, and the `.czi` string is threaded through the whole pipeline.

**Where CZI reading lives:**
- **io.py:9,42-51** — `from bioio import BioImage`; `channel_mip(czi, idx)` opens the
  CZI, requests `get_image_data("CZYX")`, max-projects one channel. This is the core
  read primitive every stage funnels through.
- **inventory.py:6-31** — `BioImage(path)`, reads `.dims`, `.physical_pixel_sizes`,
  `.channel_names`; produces `czi_manifest.csv` (the pixel-size source of truth).
- **extract.py:27** — calls `channel_mip` to write the MIP tiff.
- **omezarr.py:80-95** — `czi_to_zarr()`: `BioImage(czi_path).get_image_data("CZYX")`
  → OME-NGFF pyramid.

**What it produces:** a `(C,Z,Y,X)` numpy array + channel names + physical pixel
size, consumed downstream as (a) the manifest CSV, (b) per-channel MIP tiffs, (c)
OME-NGFF stores.

**The `.czi` string leak** (format assumed in path/key construction, not just in the
reader): io.py:27,29; inventory.py:37; mitoseg.py:45; cellseg.py:117,124;
quantify.py:187,196; omezarr.py:297,303,385; mitoseg_viz.py:58; cellseg_viz.py:52;
quantify_viz.py:66. Manifest keys are `(experiment, f"{filestem}.czi")` and raw
paths are `cfg.raw / f"{filestem}.czi"` everywhere.

**Where the seam SHOULD be:** `omezarr.czi_to_zarr` is already almost the right
shape — it is the one place that turns vendor bytes into the open substrate. The
target design is:
- **Vendor→NGFF conversion is outsourced** to the user / per-repo (bioio is one
  adapter; `czi_to_zarr` becomes an optional convenience, not core).
- **mitomorph core consumes OME-NGFF only** — `inventory`/`extract`/`mitoseg`/
  `cellseg`/`quantify` should read arrays + pixel size + channel metadata from a
  zarr store (or an in-memory array), never from `BioImage`. That removes bioio
  from io.py and inventory.py entirely and eliminates the `.czi` string from every
  path/key.
- An **`ImageReader` protocol** (`read_array`, `pixel_size`, `channel_names`) lets a
  user plug any vendor format without the package importing it.

---

## 3. Analysis-repo structure assumptions

The package assumes one fixed on-disk layout everywhere. `ExperimentConfig`
properties fix the roots (config.py:66-75):
- `raw = root/00_data/raw/<experiment>`
- `results = root/03_results`
- `manifest = results/00_inventory/tables/czi_manifest.csv`

Stage dirs are hardcoded string literals in each module:
- `02_extract/<group>/<filestem>/{figures,tables}` — extract.py:32.
- `03_mitoseg/<group>/<filestem>/{masks,tables,figures}` + `_skipped/<exp>` +
  `03_mitoseg_montage.png` — mitoseg.py:39,54,65; mitoseg_viz.py:43,65,81.
- `04_cellseg/<group>/<filestem>/{masks,tables,figures,annotations}` —
  cellseg.py:125,133,158; cellseg_viz.py:35-37; omezarr.py:187-189,334.
- `05_quantify/<group>/<filestem>/tables` + `_skipped` — quantify.py:181,197-242.
- `03_results/master/` + `05_quantify/*/*/tables/*{suffix}` glob — integrate.py:27,76,107.
- `00_ome_zarr/<group>/<filestem>.zarr` — omezarr.py:301,385.
- `run_manifest/` — provenance.py:132.

**Flag:** every stage both computes AND owns where its output lands, in a naming
convention (`NN_stage/<group>/<sample>`) specific to these repos. `integrate._collect`
even hardcodes the `*/*/tables/` glob depth (integrate.py:76), coupling it to the
`group/sample` nesting. **Seam:** output location must be user-specified. A layout
resolver (or explicit path args to each `run`) should map (stage, sample) → path;
the default resolver reproduces the current convention for these repos, but the
core functions take paths. Compute functions themselves should return arrays/rows,
not write to a fixed tree.

---

## 4. Config coupling

`config.py` resolves two files:
- **`experiment.toml` in cwd/root** (config.py:127-168) — the only entry point;
  `cli.py:61` calls `load_config(args.root)` and `--root` defaults to `.`. Requires
  keys `experiment`, `sample_regex`, `mito_channel_index`, `nucleus_channel_index`.
  **Reasonable:** per-repo identity separated from knobs. **Leak:** it is the *only*
  way to construct a config; a library user cannot pass a dict/object. Also the
  identity model itself (regex + exactly two channel indices) is experiment-shaped
  (#9, #11).
- **`run-config.yaml` by walking up parents** (config.py:104-110) — `_find_run_config`
  walks `root` + all parents looking for `config/run-config.yaml`. This hardwires the
  **umbrella superproject layout**: it assumes the experiment repo is nested under a
  parent that holds `config/run-config.yaml`. **Leak:** a standalone repo or a napari
  plugin user has no such umbrella. **Seam:** accept an explicit config path or an
  already-merged dict; parent-walk becomes an optional convenience of a CLI adapter,
  not core behavior.
- **`schema_version` assert** (config.py:85,148-154) — raises if the run-config's
  `schema_version != 1`. Reasonable as a config-contract guard; only the fact that
  it is triggered by the umbrella-walked file couples it (decouple with #7).
- **`config_version` / `resolved_knobs`** — provenance plumbing, generic and fine.

The knob set (`_KNOB_FIELDS`, config.py:88-101) is domain-specific (cellpose,
watershed, Li floor) but that is legitimate algorithm parameterization, not a
genericity leak — it just needs to travel with the compute functions rather than a
monolithic `ExperimentConfig`.

---

## 5. Compute / viz / io seams as they exist today

- **Pure compute (no I/O, no layout):** `metrics.py` (skeleton/shape/focus — clean,
  fully generic), `segment.py` (Li mask, nuclei, Cellpose bodies, territory
  watershed — takes arrays + `cfg` + `px_um`, returns arrays; the reusable core).
  These are the parts already close to library-grade.
- **Compute + orchestration mixed:** `mitoseg.py`, `cellseg.py`, `quantify.py`,
  `inventory.py`, `extract.py`. Each has a pure kernel but their `compute_one`
  functions read the manifest, build fixed paths, load `.npy`/tiff, `mkdir`, and
  `write_csv`. So "compute" modules **do** import io paths and assume directory
  structure — e.g. quantify.py:205-231 hardcodes six upstream paths across two stage
  dirs; `_collect` in integrate.py assumes the glob depth. This is the main
  compute↔io entanglement.
- **Viz:** `viz_common.py` (pure helpers, generic), `mitoseg_viz.py`,
  `cellseg_viz.py`, `quantify_viz.py` — read the same fixed tree via `_load`, write
  PNGs into the stage dirs. Correctly separated from compute (compute writes tables,
  viz reads them) per the project rule, but each hardcodes the layout (#6).
- **io / orchestration:** `io.py` (manifest read, sample resolution, CZI projection —
  carries the vendor + `.czi` coupling), `cli.py` (argparse dispatch; also decides
  which stages stamp provenance), `config.py`, `provenance.py`, `omezarr.py`
  (the one module already modelling the NGFF substrate + label-version contract —
  the intended future core).

**Coupling summary:** the pure kernels (`segment.py`, `metrics.py`, `viz_common.py`,
and most of `omezarr.py`'s read/write primitives) are genericity-ready. Everything
that blocks agnosticism is concentrated in (a) `integrate.py`'s design logic,
(b) `io.py`/`inventory.py`'s direct CZI dependence, (c) the fixed `NN_stage/<group>/
<sample>` layout threaded through every `run`/`compute_one`, and (d) `config.py`'s
umbrella-walk + toml-only + two-channel identity model.
