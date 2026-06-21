#!/usr/bin/env python3
"""Run Path-A C0f globalization probe."""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from stage1_solver import patha_c0f_globalization_probe as c0f


def main() -> int:
    return c0f.main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
