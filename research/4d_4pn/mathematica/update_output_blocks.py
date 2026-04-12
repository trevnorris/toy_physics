#!/usr/bin/env python3
"""Refresh trailing Mathematica Output: blocks for the current 4d_4pn WL manifest."""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess


ROOT = pathlib.Path(__file__).resolve().parent
MANIFEST = ROOT / "wl_notes.txt"
OUTPUT_BLOCK_RE = re.compile(r'\n\(\*"\nOutput:\n.*?\n"\*\)\s*\Z', re.S)


def strip_warnings(stdout: str) -> str:
    lines = []
    for line in stdout.splitlines():
        if "OMP: Warning #179" in line:
            continue
        if "Function Can't set size of SHM failed" in line:
            continue
        lines.append(line)
    return "\n".join(lines).rstrip() + "\n"


def load_manifest() -> list[pathlib.Path]:
    files: list[pathlib.Path] = []
    for raw in MANIFEST.read_text().splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        files.append(ROOT / raw)
    return files


def refresh_file(path: pathlib.Path, dry_run: bool = False) -> None:
    proc = subprocess.run(
        ["math", "-script", str(path)],
        capture_output=True,
        text=True,
        check=True,
        cwd=ROOT.parent.parent.parent,
    )
    stdout = strip_warnings(proc.stdout)
    text = path.read_text()
    text = OUTPUT_BLOCK_RE.sub("", text).rstrip() + "\n\n"
    text += '(*"\nOutput:\n'
    text += stdout
    if not stdout.endswith("\n"):
        text += "\n"
    text += '"*)\n'
    if dry_run:
        print(f"would update {path.name}")
    else:
        path.write_text(text)
        print(f"updated {path.name}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    for path in load_manifest():
        refresh_file(path, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
