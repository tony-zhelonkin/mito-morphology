"""Per-experiment configuration loaded from an experiment.toml at the repo root.

Everything experiment-specific lives here (sample-name regex, channel indices,
QC tuning constants), so the compute code stays generic and shared.
"""
from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass, field
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

    # OME-Zarr substrate (NGFF); Zarr v2 for broadest tool support.
    ngff_version: str = "0.4"

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


def load_config(root: Path | str = ".") -> ExperimentConfig:
    """Load experiment.toml from a repo root into an ExperimentConfig."""
    root = Path(root).resolve()
    path = root / "experiment.toml"
    with path.open("rb") as fh:
        data = tomllib.load(fh)
    qc = data.get("qc", {})
    ome = data.get("ome_zarr", {})
    return ExperimentConfig(
        root=root,
        experiment=data["experiment"],
        sample_regex=data["sample_regex"],
        mito_channel_index=int(data["mito_channel_index"]),
        nucleus_channel_index=int(data["nucleus_channel_index"]),
        skip_samples=list(data.get("skip_samples", [])),
        min_nucleus_area_um2=float(qc.get("min_nucleus_area_um2", 15.0)),
        nucleus_seed_min_distance_um=float(qc.get("nucleus_seed_min_distance_um", 3.0)),
        max_cell_radius_um=float(qc.get("max_cell_radius_um", 12.0)),
        mito_min_um2=float(qc.get("mito_min_um2", 1.0)),
        focus_rel_factor=float(qc.get("focus_rel_factor", 0.33)),
        ngff_version=str(ome.get("ngff_version", "0.4")),
    )
