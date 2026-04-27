#!/usr/bin/env python3
"""Run the reader-site math tests, generator, and validator."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ATLAS_SCRIPTS = ROOT / "atlas" / "scripts"


PIPELINE = (
    ("reader math fixtures", ATLAS_SCRIPTS / "test_reader_math_normalization.py"),
    ("reader site generator", ATLAS_SCRIPTS / "generate_reader_site.py"),
    ("reader site validator", ATLAS_SCRIPTS / "validate_reader_site.py"),
)


def run_step(label: str, script: Path) -> int:
    print(f"==> {label}", flush=True)
    completed = subprocess.run([sys.executable, str(script)], cwd=ROOT, check=False)
    if completed.returncode:
        print(f"FAILED: {label}", file=sys.stderr)
    return completed.returncode


def main() -> None:
    for label, script in PIPELINE:
        returncode = run_step(label, script)
        if returncode:
            sys.exit(returncode)
    print("OK: reader site pipeline completed", flush=True)


if __name__ == "__main__":
    main()
