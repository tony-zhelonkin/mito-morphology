"""Run provenance: an append-only manifest stamped next to results each run.

The run-level twin of the OME-Zarr label `.zattrs`: every invocation records the
exact resolved parameters, the mitomorph version + git commit (and dirty flag),
and a timestamp — so any `master_table.csv` row is traceable to (config_version +
mitomorph SHA + timestamp). One file per run; never overwritten.
"""
from __future__ import annotations

import subprocess
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from .config import ExperimentConfig, resolved_knobs


def _mitomorph_version() -> str:
    try:
        return metadata.version("mito-morphology")
    except Exception:
        return "unknown"


def _git_state() -> tuple[str, bool]:
    """(short SHA, dirty?) of the mitomorph source tree, best-effort."""
    pkg_dir = Path(__file__).resolve().parent
    try:
        sha = subprocess.run(
            ["git", "-C", str(pkg_dir), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "-C", str(pkg_dir), "status", "--porcelain"],
            capture_output=True, text=True, check=True,
        ).stdout.strip())
        return sha, dirty
    except Exception:
        # Installed from a pinned git rev / wheel: no working tree to inspect.
        return "n/a", False


def write_run_manifest(
    cfg: ExperimentConfig, command: str, samples: list[str] | None
) -> Path:
    """Write an append-only run manifest under 03_results/run_manifest/."""
    import yaml

    sha, dirty = _git_state()
    stamp = datetime.now(timezone.utc)
    manifest = {
        "run": {
            "timestamp_utc": stamp.isoformat(timespec="seconds"),
            "command": command,
            "experiment": cfg.experiment,
            "samples": list(samples) if samples else "all",
        },
        "config": {
            "schema_version": cfg.schema_version,
            "config_version": cfg.config_version,
            "resolved_knobs": resolved_knobs(cfg),
        },
        "mitomorph": {
            "version": _mitomorph_version(),
            "git_sha": sha,
            "git_dirty": dirty,
        },
    }
    out_dir = cfg.results / "run_manifest"
    out_dir.mkdir(parents=True, exist_ok=True)
    # filename-safe timestamp; one file per run (never overwritten).
    fname = f"{stamp.strftime('%Y%m%dT%H%M%SZ')}_{cfg.config_version}_{command.replace(' ', '-')}.yaml"
    out = out_dir / fname
    out.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(f"run manifest: {out}", flush=True)
    return out
