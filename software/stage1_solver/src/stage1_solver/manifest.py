"""Run manifest writing for the freeze-boundary discipline."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
from typing import Any

from .backend import library_versions
from .config import stable_hash


def git_revision() -> dict[str, str]:
    def run(args: list[str]) -> str:
        try:
            return subprocess.check_output(args, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unavailable"

    return {
        "head": run(["git", "rev-parse", "HEAD"]),
        "status_short": run(["git", "status", "--short"]),
    }


def write_manifest(
    *,
    run_root: str,
    benchmark_name: str,
    grid_name: str,
    config: dict[str, Any],
    mesh: dict[str, Any],
    results: dict[str, Any],
) -> Path:
    run_dir = Path(run_root) / benchmark_name / grid_name
    run_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "benchmark": benchmark_name,
        "grid": grid_name,
        "dtype": config["backend"]["dtype"],
        "device": config["backend"]["device"],
        "solver_controls": config["newton"],
        "mesh": mesh,
        "config_hash": stable_hash(config),
        "git": git_revision(),
        "library_versions": library_versions(),
        "results": results,
    }
    path = run_dir / "manifest.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path
