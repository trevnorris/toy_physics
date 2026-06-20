#!/usr/bin/env python3
"""Run Path-A C0c null-mode identification from saved C0b artifacts."""

from __future__ import annotations

import sys

from stage1_solver import patha_c0_conditioning_spike as c0


def main() -> int:
    return c0.c0c_main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
