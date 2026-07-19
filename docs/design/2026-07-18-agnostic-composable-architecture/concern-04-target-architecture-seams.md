# Concern 04 — Target Architecture: Seams for an Agnostic, Composable `mitomorph`

Status: design proposal (seams, not a rewrite spec). Date: 2026-07-18.

This document proposes the *boundaries and interfaces* for a redesign in which
`mitomorph` becomes agnostic to experiment, repo layout, and vendor format;
composable across mask-creation methods, cell recipes, and metrics; and
observable (every operation openable in napari). It builds on what already
exists — the OME-NGFF substrate (`omezarr.py`), the physical-units convention
(`px_um` everywhere), and the immutable-label / provenance discipline — and
proposes where to cut so those good parts survive while the couplings go.

## Where the couplings live today (baseline)

| Coupling | Location | Target |
|---|---|---|
| CZI is read directly mid-pipeline | `io.channel_mip(czi, idx)` called in `extract`/`cellseg`/`quantify`; `cfg.raw` = CZI dir; manifest = `czi_manifest.csv` | Core reads **NGFF only**; CZI conversion is an out-of-core adapter |
| Fixed on-disk layout | `02_extract/…03_mitoseg…04_cellseg…05_quantify…/master` hardwired as string paths in every stage | Substrate is the `.zarr` store; sidecar tables go to a **user-specified out dir** |
| Fixed linear CLI | `cli.py` if/elif ladder; `all` = fixed order | **Composable op stack** (recipe), CLI runs a recipe |
| Only one mito method | `segment.pipeline_li` hardwired in `mitoseg.compute_one` | **Registry** of mito-mask methods |
| Only one cell recipe | `cellseg.compute_one` hardwires nuclei+Cellpose | **Registry** of cell-mask recipes |
| Metrics hardwired | `quantify.CELL_FIELDS` + inline computation in `_cell_metrics` | **Registry** of self-describing metrics |
| Experiment identity in library | `integrate.KNOWN_EXPERIMENTS`, `_design_columns` parse `TBK1i`/`gln_mM`; `config.parse_sample` regex | Moves **out to the repo** (design callback) |
| Viz = separate matplotlib | `*_viz.py` render PNGs from tables after the fact | Every op emits a **napari-renderable view**; same op headless or in GUI |

---

## 1. Module boundary map

Four rings, dependencies point **inward only**. Core never imports outward.

```
┌─────────────────────────────────────────────────────────────────────┐
│ per-repo experiment glue  (lives in each experiment repo, NOT library)│
│   experiment.toml · design.py (group→factors) · CZI→NGFF conversion   │
│   choice of recipe + params · out-dir location                        │
└───────────────┬───────────────────────────────────────────────────────┘
                │ depends on
┌───────────────▼───────────────┐   ┌──────────────────────────────────┐
│ orchestration (recipe runner)  │   │ napari-plugin layer               │
│   RecipeSpec → run headless    │◄──┤   same Operations, run interactively│
│   provenance log · out-dir mgmt│   │   step views · expert-override capture│
└───────────────┬────────────────┘   └───────────────┬──────────────────┘
                │            both depend on           │
        ┌───────▼─────────────────────────────────────▼───────┐
        │ CORE (pure)                                          │
        │   registries: mito-methods · cell-recipes · metrics  │
        │   Operation protocol · OpResult (array + view spec)  │
        │   NGFF substrate contract (Substrate: read/write     │
        │   images + labels + attrs)                           │
        │   pure array ops (segment.py, metrics.py — unchanged  │
        │   math, no I/O)                                       │
        └───────▲──────────────────────────────────────────────┘
                │ depends on (adapter implements Substrate write)
        ┌───────┴────────────────────────────────────┐
        │ io/adapters (vendor → NGFF), OUTSIDE core   │
        │   czi_adapter.py (bioio) · generic tiff/nd2 │
        │   NOT imported by core; produces the .zarr   │
        └──────────────────────────────────────────────┘
```

Key rules:
- **Core is NGFF-in / labels-out and pure.** It imports numpy/skimage/zarr, never
  `bioio`, never an experiment id, never a hardcoded results path. `segment.py`
  and `metrics.py` are *already* pure array functions — they stay, they just get
  wrapped as registered Operations.
- **Adapters are outside core.** `czi_to_zarr` (today inside `omezarr.py`) moves
  to `mitomorph.adapters.czi` (or a separate `mitomorph-czi` extra). Core depends
  on the `Substrate` protocol; the adapter is one producer of a conforming store.
  "CZI just happened to be what the Zeiss scope produced" is enforced structurally:
  delete the adapter and core still runs on any NGFF store.
- **Orchestration** owns the out-dir (idempotent, unaware of repo organization —
  the user passes `--out /wherever`), the recipe runner, and the provenance log.
- **napari-plugin** re-runs the *identical* Operation objects; it is the only place
  that knows about the GUI. Core emits view specs; the plugin renders them, headless
  emits them as PNG through the same spec.

---

## 2. Registry / plugin pattern for (A) mito methods, (B) cell recipes, (C) metrics

One shared idea: a **decorator-registry keyed by name**, each entry a callable with
a typed signature and a self-describing spec. Third parties (or a per-repo module)
register by importing and decorating; nothing in core enumerates the members.

```python
# core/registry.py — one generic registry, three instances
from typing import Callable, Generic, TypeVar
T = TypeVar("T")

class Registry(Generic[T]):
    def __init__(self, kind: str): self.kind, self._items = kind, {}
    def register(self, name: str) -> Callable[[T], T]:
        def deco(obj: T) -> T:
            if name in self._items:
                raise ValueError(f"{self.kind} {name!r} already registered")
            self._items[name] = obj; return obj
        return deco
    def get(self, name: str) -> T: return self._items[name]
    def names(self) -> list[str]: return sorted(self._items)

MITO_METHODS  = Registry["MitoMethod"]("mito-method")
CELL_RECIPES  = Registry["CellRecipe"]("cell-recipe")
METRICS       = Registry["Metric"]("metric")
```

### (A) Mito-mask method

A method turns one channel MIP + pixel size + params into a boolean mask, and
declares its params + view. Li is the first registrant; Otsu / ML / Nellie slot in
beside it with no core change.

```python
@dataclass(frozen=True)
class MitoMethod:
    name: str
    fn: Callable[[np.ndarray, float, Mapping], np.ndarray]  # (mip, px_um, params) -> bool mask
    param_schema: Mapping[str, ParamSpec]                   # name -> {type, default, doc, units}
    describe: str                                           # one-line "what it does"

@MITO_METHODS.register("li")
def _li(mip, px_um, params):                # wraps existing segment.pipeline_li
    ...
# registered as MitoMethod(name="li", fn=_li, param_schema={"min_object_um2": ...}, ...)
```

### (B) Cell-mask recipe

A recipe is *composed from named channel inputs*, so "mito + nucleus" and
"membrane + nucleus" are two registrants differing only in which channels they
request. The recipe declares its **required channel roles** — the per-repo config
maps roles → channel indices, keeping the library blind to a given scope's layout.

```python
@dataclass(frozen=True)
class CellRecipe:
    name: str
    required_roles: tuple[str, ...]         # e.g. ("mito","nucleus") or ("membrane","nucleus")
    fn: Callable[[Mapping[str, np.ndarray], float, Mapping], np.ndarray]  # roles->MIPs -> int labels
    param_schema: Mapping[str, ParamSpec]
    describe: str

@CELL_RECIPES.register("cellpose_mito_nucleus")   # today's build_territories_v2
@CELL_RECIPES.register("membrane_nucleus")        # future: watershed on membrane
```

### (C) Metric — self-describing + napari-visual contract

Every metric has a "parking spot": a name, a transparent definition of *what and
how*, its output dtype/units, the mask level it consumes (`per_object` |
`per_cell` | `per_image`), and — crucially — an optional `view` that returns the
napari layers illustrating what it measured. Non-obvious metrics
(`form_factor`, `solidity`, `branch_count`, `mean_diameter_um`) implement `view`;
trivial ones (`area_um2`) may skip it.

```python
@dataclass(frozen=True)
class Metric:
    name: str
    level: str                              # "per_object" | "per_cell" | "per_image"
    units: str                              # "um2", "count", "dimensionless", ...
    describe: str                           # transparent definition of the computation
    fn: Callable[[MetricInput], float]      # the number
    view: Callable[[MetricInput], list["LayerSpec"]] | None = None  # napari visual of the metric

@METRICS.register  # form_factor: 4πA/P² — view draws the object outline + its
                   # equal-area circle so "how round" is visually obvious.
def form_factor() -> Metric: ...
```

A user who needs a metric mitomorph lacks writes a `Metric` and registers it from
their repo; a user who wants only one metric asks for `metrics=["form_factor"]`;
"all" = `METRICS.names()`. The registry *is* the transparent catalogue.

---

## 3. Composable operation model (identical headless & interactive)

The pipeline is a **stack of typed Operations**, not the `cli.py` if/elif ladder.
An `Operation` is a pure-ish transform over the substrate: it reads named layers,
writes named layers and/or table rows, and returns an `OpResult` carrying the
array(s) **and** a `view` spec. The *same* `Operation` objects run in a headless
runner and inside napari — the only difference is who consumes the `view`.

```python
class Operation(Protocol):
    name: str
    params: Mapping[str, object]
    def inputs(self) -> tuple[str, ...]: ...           # substrate layer/table names read
    def outputs(self) -> tuple[str, ...]: ...          # names written
    def run(self, ctx: "OpContext") -> "OpResult": ...

@dataclass
class OpResult:
    arrays: dict[str, np.ndarray]          # e.g. {"labels/mito": mask}
    rows:   dict[str, list[dict]]          # table fragments (metric outputs)
    view:   "ViewSpec"                     # napari layers + a title = one "sequence" frame
    params: Mapping[str, object]           # frozen, for the log

@dataclass
class ViewSpec:                            # renders the SAME in napari or as a PNG
    layers: list[LayerSpec]                # image/labels/shapes with blending + colormap
    title: str                            # step caption for the paper figure
```

A **RecipeSpec** is the serializable stack (order + params) — this is what makes
"different parameters, different order, as the science requires" a first-class,
reproducible object rather than CLI argument order:

```yaml
# recipe.yaml — headless and GUI both consume this
substrate: {store: "path/or/-", channels: {mito: 0, nucleus: 1}}
steps:
  - op: mito_mask     {method: li, params: {min_object_um2: 0.14}}
  - op: cell_mask     {recipe: cellpose_mito_nucleus, params: {cellprob_threshold: -2.0}}
  - op: metrics       {names: [area_um2, form_factor, branch_count], level: per_cell}
```

- **Headless**: `run_recipe(spec, substrate, out)` executes each op, writes arrays to
  `labels/`, appends rows to tables in the out-dir, appends each `view` frame to a
  `steps/` figure sequence, and appends to the provenance log (§4).
- **Interactive**: the napari plugin loads the same `RecipeSpec`, and for each op
  shows its `ViewSpec` live; the scientist can re-run a single op with tweaked
  params, reorder, or inject a hand-drawn override. Their session serializes back
  to the identical `RecipeSpec` + override events — so the exploratory GUI session
  *is* a reproducible recipe.

The **step-sequence provenance** (the ordered list of `(op.name, params, view.title,
output layers, git/version)`) is exactly the "sequence of events" a paper/grant
figure needs: each frame is one op's `ViewSpec` rendered, captioned by its
`describe`.

---

## 4. OME-NGFF contract (substrate + labels + provenance + expert override)

The `.zarr` store *is* the working substrate — not a parallel export as
`migrate-labels` treats it today. Ops read images and labels from it and write
labels back into it.

```
sample.zarr/
  0/ … 3/                      image multiscale pyramid  (from an adapter; core never writes 0/)
  .zattrs: omero, multiscales, mitomorph_channels: {mito: 0, nucleus: 1}  ← role map
  labels/
    mito/            .zattrs.mitomorph_provenance = {method, params, version, git_sha}
    cells/           .zattrs.mitomorph_provenance = {recipe, params, ...}
    cells_expert/    .zattrs.mitomorph_provenance = {override_of: "cells", event: {...}}
    .zattrs.labels = [...]     (deduped, order-preserving)
  .zattrs.mitomorph_log = [ {op, params, inputs, outputs, ts, git_sha}, ... ]  ← the op stack
```

- **Image substrate**: written once by an adapter; core treats `0/` as read-only.
  `Substrate` protocol = `read_image(role) -> (mip, px_um)`, `read_labels(name)`,
  `write_labels(name, arr, provenance)`, `append_log(entry)`. Reuses the existing
  `_pyramid_scale_transforms` so labels align with the image at every level.
- **Labels/ groups**: keep today's discipline — `write_labels` immutable-versioned
  for the analytic A→B seam, `write_named_labels` idempotent for fixed roles. The op
  runner writes each op's output as a label group carrying its params in
  `mitomorph_provenance`.
- **Op log in attrs**: `.zattrs.mitomorph_log` is the on-substrate twin of the recipe
  stack — the store is self-describing about how every label was made, independent of
  the out-dir. (Complements `provenance.write_run_manifest`, which stays as the run-level record.)
- **Expert override** (in napari): a hand-drawn label group `cells_expert` is written;
  the runner **detects overlap** with the mask it targets (`cells`), and rather than
  silently clobbering, records an *override event* — `{override_of, iou, n_objects_changed,
  drawn_by, ts}` — into both the new group's provenance and the op log, then downstream
  ops consume `cells_expert` where it overlaps. The composed stack + params + the
  override are all logged; the original `cells` is never destroyed (immutable-version
  discipline already in place).

---

## 5. Moving experiment/design metadata OUT of the library

Today `integrate.py` hardcodes `KNOWN_EXPERIMENTS` and parses `TBK1i`/`gln_mM`; and
`config.parse_sample` bakes a group/replicate regex into the shared library. Both
move to the repo. Core integrate becomes **reshape-only** and takes a *callback*:

```python
# core/integrate.py — generic; no experiment id anywhere
def integrate(tables: Iterable[Path], design: "DesignFn | None") -> pd.DataFrame:
    df = pd.concat(map(pd.read_csv, tables), ignore_index=True)
    if design is not None:
        df = df.join(df.apply(design, axis=1, result_type="expand"))   # add factor columns
    return df

DesignFn = Callable[[pd.Series], dict]   # per-row: sample metadata -> design factors
```

The per-repo `design.py` supplies the callback (and the sample-name parsing that used
to be `parse_sample`):

```python
# experiment repo: design.py  (TBK1i example — lives WITH the data, not in the library)
def design(row):
    g = row["sample"]                       # repo owns naming; library stays blind
    return {"tbk1i": int("TBK1i" in g), "isd90": int("ISD90" in g)}
```

`experiment.toml` gains a pointer (`design = "design:design"`) that orchestration
imports; the library never enumerates experiments and never warns about "unknown
experiment id". Channel *roles* (mito/nucleus/membrane) likewise move into the
per-repo config as a role→index map, replacing `mito_channel_index` /
`nucleus_channel_index` constants — so a scope with a membrane channel just adds a
role.

---

## 6. Phased migration sketch (high level, mostly non-breaking)

1. **Extract pure ops behind registries (non-breaking).** Wrap existing
   `segment.pipeline_li`, `build_territories_v2`, and the `metrics.py` functions as
   `MitoMethod`/`CellRecipe`/`Metric` registrants. Current stages call the registry
   with the single registered name — identical outputs, seams now exist.
2. **Invert the substrate (non-breaking read path).** Make ops read the MIP + roles
   from the `.zarr` store (already produced by `to-zarr`) instead of re-reading CZI
   via `channel_mip`. CZI conversion becomes the only place `bioio` is imported;
   move it to `adapters/czi.py`. Sidecar tables still land in the current layout.
3. **Introduce `RecipeSpec` + runner (additive).** New `mito run recipe.yaml` beside
   the existing subcommands; `all` becomes a built-in default recipe. Op log written
   to `.zattrs`.
4. **De-hardcode integrate (mildly breaking for the two repos).** Replace
   `KNOWN_EXPERIMENTS`/`_design_columns` with the `DesignFn` callback + per-repo
   `design.py`. Add `design.py` to both existing experiment repos so their master
   tables are unchanged.
5. **napari plugin over the same ops (additive).** Plugin runs `RecipeSpec` step by
   step, renders each `ViewSpec`, captures expert overrides into the op log. The old
   matplotlib `*_viz.py` become the headless renderer of the same `ViewSpec`s, then
   can be retired.
6. **Out-dir decoupling (breaking, last).** Replace hardwired `02_extract…05_quantify`
   strings with a user-specified `--out`; keep the current tree as the default so the
   two repos keep working until they opt in.

3D is deferred but not designed out: `_label_axes` already handles 2D/3D, `px_um`
generalizes to a scale vector, and Operations are written over ndim-generic arrays;
the MIP-collapse is itself just the first Operation in the stack (`op: mip_collapse`),
so a future 3D recipe simply omits it.
