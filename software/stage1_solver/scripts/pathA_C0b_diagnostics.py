#!/usr/bin/env python3
"""Run only the bounded Path-A C0b linear/admissibility diagnostic phase."""

from __future__ import annotations

import sys

from stage1_solver import patha_c0_conditioning_spike as c0


def main() -> int:
    return c0.main(["--phase", "diagnostics", *sys.argv[1:]])


if __name__ == "__main__":
    raise SystemExit(main())
