#!/usr/bin/env python3
"""
Filesystem scanner for an audit unit's artifacts.

Usage:
  scan.py <config-path> <unit-number-as-int>

Prints a YAML document on stdout describing what exists for that unit:

  stage: 58
  stage_str: "058"
  sympy:              {path: "scripts/...", exists: true, mtime: 1736969231}
  mathematica:        {path: "mathematica/...", exists: true, mtime: 1736969231}
  sympy_output:       {path: "scripts/output/...", exists: true, mtime: 1736969500}
  mathematica_output: {path: "mathematica/output/...", exists: true, mtime: 1736969500}

Paths are stored relative to the project root (the directory containing
.redteam-config.yaml) so the manifest is portable across checkout locations.
Callers resolve to absolute by joining with the project root at use-site.

Each field is null/{exists: false} when no match is found.
Multiple matches collapse to the first (lexicographic).

Note: tex and notes files are NOT scanned. The red-team is scripts-only;
doc alignment is handled manually, out-of-band.
"""
import os
import sys
import glob
import yaml
from pathlib import Path

def load_cfg(path):
    with open(path) as f:
        return yaml.safe_load(f)

def first_match(project_root, patterns, stage_str):
    """Return (relative_path, absolute_path) for the first match, or (None, None)."""
    for pat in patterns or []:
        expanded = pat.format(N=stage_str)
        full = str(project_root / expanded)
        matches = sorted(glob.glob(full))
        if matches:
            abs_path = matches[0]
            try:
                rel = str(Path(abs_path).resolve().relative_to(project_root.resolve()))
            except ValueError:
                # Match is outside the project root (unusual); keep absolute.
                rel = abs_path
            return rel, abs_path
    return None, None

def file_info(rel_path, abs_path):
    if rel_path is None:
        return {"path": None, "exists": False}
    try:
        st = os.stat(abs_path)
        return {"path": rel_path, "exists": True, "mtime": int(st.st_mtime)}
    except OSError:
        return {"path": rel_path, "exists": False}

def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    cfg_path = Path(sys.argv[1]).resolve()
    project_root = cfg_path.parent
    stage = int(sys.argv[2])
    cfg = load_cfg(cfg_path)
    pad = int(cfg.get("stage_pad", 3))
    stage_str = f"{stage:0{pad}d}"
    paths_cfg = cfg.get("paths", {}) or {}

    out = {
        "stage": stage,
        "stage_str": stage_str,
        "sympy": file_info(*first_match(project_root, paths_cfg.get("sympy"), stage_str)),
        "mathematica": file_info(*first_match(project_root, paths_cfg.get("mathematica"), stage_str)),
        "sympy_output": file_info(*first_match(project_root, paths_cfg.get("sympy_output"), stage_str)),
        "mathematica_output": file_info(*first_match(project_root, paths_cfg.get("mathematica_output"), stage_str)),
    }

    yaml.safe_dump(out, sys.stdout, default_flow_style=False, sort_keys=False)

if __name__ == "__main__":
    main()
