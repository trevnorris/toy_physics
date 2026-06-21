#!/usr/bin/env python3
"""Run the Path-A C0g B-3 deep-point characterization battery."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from stage1_solver import patha_c0_conditioning_spike as c0


def main() -> int:
    return c0.c0g_b3_deeppoint_main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
