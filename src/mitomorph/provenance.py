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


def _pep610_commit() -> str | None:
    """Full commit SHA from PEP 610 direct_url.json, if mitomorph was pip-installed
    from a git URL (the normal production path in a pixi-pinned experiment env)."""
    try:
        raw = metadata.distribution("mito-morphology").read_text("direct_url.json")
    except Exception:
        return None
    if not raw:
        return None
    try:
        import json
        data = json.loads(raw)
        return data.get("vcs_info", {}).get("commit_id")
    except Exception:
        return None


def _looks_like_mitomorph_checkout(pkg_dir: Path) -> bool:
    """Guard against `git -C pkg_dir` walking up into an unrelated repo (e.g. an
    experiment repo that happens to contain the pip-installed package nested inside
    its working tree, under .pixi/envs/.../site-packages/mitomorph/)."""
    parts = pkg_dir.parts
    if "site-packages" in parts or ".pixi" in parts:
        return False
    try:
        toplevel = Path(subprocess.run(
            ["git", "-C", str(pkg_dir), "rev-parse", "--show-toplevel"],
            capture_output=True, text=True, check=True,
        ).stdout.strip())
    except Exception:
        return False
    for name in ("pyproject.toml", "setup.py", "setup.cfg"):
        f = toplevel / name
        if f.exists():
            try:
                text = f.read_text(errors="ignore")
            except Exception:
                continue
            if "mito-morphology" in text or "mitomorph" in text:
                return True
    return False


def git_provenance() -> dict:
    """Resolve mitomorph provenance in priority order:
    1. PEP 610 direct_url.json (pinned pip-from-git install) — the normal production path.
    2. A genuine mitomorph git checkout (editable/dev install), guarded so we never
       fall back to an unrelated repo that happens to contain the installed package.
    3. Unknown.
    Returns {"git_sha": str, "git_dirty": bool, "source": "pep610"|"git"|"unknown"}.
    """
    commit = _pep610_commit()
    if commit:
        return {"git_sha": commit[:7], "git_dirty": False, "source": "pep610"}

    pkg_dir = Path(__file__).resolve().parent
    if _looks_like_mitomorph_checkout(pkg_dir):
        try:
            sha = subprocess.run(
                ["git", "-C", str(pkg_dir), "rev-parse", "--short", "HEAD"],
                capture_output=True, text=True, check=True,
            ).stdout.strip()
            dirty = bool(subprocess.run(
                ["git", "-C", str(pkg_dir), "status", "--porcelain"],
                capture_output=True, text=True, check=True,
            ).stdout.strip())
            return {"git_sha": sha, "git_dirty": dirty, "source": "git"}
        except Exception:
            pass

    return {"git_sha": "n/a", "git_dirty": False, "source": "unknown"}


def _git_state() -> tuple[str, bool]:
    """(short SHA, dirty?) of the mitomorph source, best-effort. Kept for backward
    compatibility; prefer `git_provenance()` for the full picture (incl. source)."""
    prov = git_provenance()
    return prov["git_sha"], prov["git_dirty"]


def write_run_manifest(
    cfg: ExperimentConfig, command: str, samples: list[str] | None
) -> Path:
    """Write an append-only run manifest under 03_results/run_manifest/."""
    import yaml

    prov = git_provenance()
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
            "git_sha": prov["git_sha"],
            "git_dirty": prov["git_dirty"],
            "git_source": prov["source"],
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
