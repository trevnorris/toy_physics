#!/usr/bin/env python3
"""Aggregate Path-A C0g staged Steps 0-3 and write the partial report."""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from stage1_solver import patha_c0_conditioning_spike as c0


def main() -> int:
    return c0.c0g_aggregate_main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
