#!/usr/bin/env python3
"""Run Path-A C0f2 chunked timing rerun."""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from stage1_solver import patha_c0f2_timing_rerun as c0f2


def main() -> int:
    return c0f2.main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
