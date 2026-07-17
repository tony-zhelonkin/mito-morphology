"""Per-experiment configuration.

Two-file split, resolved into one `ExperimentConfig`:

- ``experiment.toml`` (per repo) = experiment IDENTITY: sample-name regex, channel
  indices, skip list. Intrinsically per-experiment, lives next to the data.
- ``config/run-config.yaml`` (centralized in the umbrella superproject, resolved by
  walking up from the experiment root) = shared pipeline KNOBS + sparse per-experiment
  overrides + provenance/changelog. This is the reproducible run contract.

Precedence (low -> high): dataclass defaults -> experiment.toml ``[qc]``/``[ome_zarr]``
-> run-config.yaml ``shared`` -> run-config.yaml ``experiments.<experiment>``.
The YAML is authoritative for knobs; the toml supplies identity (and remains a
fallback for the qc block so a repo works even without the umbrella present).
"""
from __future__ import annotations

import re
import tomllib
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class ExperimentConfig:
    """Resolved configuration for one experiment repo."""

    root: Path
    experiment: str
    sample_regex: str
    mito_channel_index: int
    nucleus_channel_index: int
    skip_samples: list[str] = field(default_factory=list)

    # QC / tuning (physical units; per-file pixel size from the manifest).
    min_nucleus_area_um2: float = 15.0
    nucleus_seed_min_distance_um: float = 3.0
    max_cell_radius_um: float = 12.0
    mito_min_um2: float = 1.0
    # Focus is flagged per-image (relative), never by a cross-experiment constant:
    # a cell is out_of_focus if its focus_score < focus_rel_factor * image median.
    focus_rel_factor: float = 0.33

    # --- Module A: cell footprint (Cellpose) + instance resolution (watershed) ---
    cellpose_model_type: str = "cyto3"
    cellprob_threshold: float = -2.0   # lower -> grow masks into dim edge mito
    flow_threshold: float = 0.7        # higher -> keep lower-quality masks
    mito_gamma: float = 0.5            # <1 brightens dim peripheral mito for recall
    elevation_lambda: float = 50.0     # watershed ridge cost across mito-empty gaps

    # --- Module B: mito mask (Li) --- physical noise floor (was a fixed 32 px;
    # scaled by px_um so the smallest retained object is the same PHYSICAL size
    # across batches with different pixel sizes).
    mito_min_object_um2: float = 0.14

    # OME-Zarr substrate (NGFF); Zarr v2 for broadest tool support.
    ngff_version: str = "0.4"

    # Provenance: which run-config produced this (stamped into the run manifest).
    schema_version: int = 1
    config_version: str = "uncontracted"  # set when a run-config.yaml is resolved

    @property
    def raw(self) -> Path:
        return self.root / "00_data" / "raw" / self.experiment

    @property
    def results(self) -> Path:
        return self.root / "03_results"

    @property
    def manifest(self) -> Path:
        return self.results / "00_inventory" / "tables" / "czi_manifest.csv"

    def parse_sample(self, filestem: str) -> tuple[str, int]:
        """Return (treatment group, replicate) from a CZI filename stem."""
        m = re.match(self.sample_regex, filestem)
        if not m:
            raise ValueError(f"could not parse sample name: {self.experiment}/{filestem}")
        return m.group("group"), int(m.group("replicate"))


EXPECTED_SCHEMA_VERSION = 1

# Knob fields that a run-config.yaml may set (identity fields stay in the toml).
_KNOB_FIELDS = {
    "min_nucleus_area_um2", "nucleus_seed_min_distance_um", "max_cell_radius_um",
    "mito_min_um2", "focus_rel_factor", "cellpose_model_type", "cellprob_threshold",
    "flow_threshold", "mito_gamma", "elevation_lambda", "mito_min_object_um2",
    "ngff_version",
}


def _find_run_config(root: Path) -> Path | None:
    """Walk up from an experiment root to the umbrella's config/run-config.yaml."""
    for d in (root, *root.parents):
        p = d / "config" / "run-config.yaml"
        if p.is_file():
            return p
    return None


def _flatten_knobs(section: dict) -> dict:
    """Flatten a run-config.yaml block's nested groups into ExperimentConfig fields.

    Groups (mito/cellpose/watershed/qc/ome_zarr) are cosmetic nesting; we pull every
    recognised knob leaf up to the flat field name the dataclass uses.
    """
    out: dict = {}
    for v in section.values():
        if isinstance(v, dict):
            out.update({k: val for k, val in v.items() if k in _KNOB_FIELDS})
    out.update({k: val for k, val in section.items() if k in _KNOB_FIELDS})
    return out


def load_config(root: Path | str = ".") -> ExperimentConfig:
    """Resolve experiment.toml (identity) + run-config.yaml (knobs) into a config."""
    root = Path(root).resolve()
    with (root / "experiment.toml").open("rb") as fh:
        data = tomllib.load(fh)

    # Start from toml [qc]/[ome_zarr] (fallback / backward compatible).
    qc = data.get("qc", {})
    ome = data.get("ome_zarr", {})
    knobs: dict = {}
    knobs.update({k: v for k, v in qc.items() if k in _KNOB_FIELDS})
    if "ngff_version" in ome:
        knobs["ngff_version"] = ome["ngff_version"]

    # Overlay the centralized run-config.yaml (shared -> per-experiment), authoritative.
    config_version = "uncontracted"
    rc_path = _find_run_config(root)
    if rc_path is not None:
        import yaml

        rc = yaml.safe_load(rc_path.read_text()) or {}
        sv = int(rc.get("schema_version", EXPECTED_SCHEMA_VERSION))
        if sv != EXPECTED_SCHEMA_VERSION:
            raise ValueError(
                f"run-config.yaml schema_version={sv} != expected {EXPECTED_SCHEMA_VERSION} "
                f"({rc_path}); update mitomorph or the config."
            )
        config_version = str(rc.get("config_version", "unversioned"))
        knobs.update(_flatten_knobs(rc.get("shared", {})))
        knobs.update(_flatten_knobs(rc.get("experiments", {}).get(data["experiment"], {}) or {}))

    typed = {k: _coerce(k, v) for k, v in knobs.items()}
    return ExperimentConfig(
        root=root,
        experiment=data["experiment"],
        sample_regex=data["sample_regex"],
        mito_channel_index=int(data["mito_channel_index"]),
        nucleus_channel_index=int(data["nucleus_channel_index"]),
        skip_samples=list(data.get("skip_samples", [])),
        config_version=config_version,
        **typed,
    )


def _coerce(key: str, value):
    """Coerce a knob value to the dataclass field type (str stays str, rest float)."""
    if key in ("cellpose_model_type", "ngff_version"):
        return str(value)
    return float(value)


def resolved_knobs(cfg: ExperimentConfig) -> dict:
    """The flat knob dict actually used — for the run provenance manifest."""
    d = asdict(cfg)
    d.pop("root", None)
    return d
