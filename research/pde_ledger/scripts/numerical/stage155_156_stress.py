#!/usr/bin/env python3
"""Retired Stage 155/156 plain harness wrapper.

The exploratory plain harness was materially weaker than the fixed-point harness
and had drifted out of maintenance. Keep the entry point for discoverability,
but delegate all work to the authoritative fixed-point stress check.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def load_fixedpoint_module():
    target = Path(__file__).with_name("stage155_156_fixedpoint_stress.py")
    spec = importlib.util.spec_from_file_location("stage155_156_fixedpoint_stress", target)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {target}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, target


def main(argv: list[str]) -> int:
    module, target = load_fixedpoint_module()
    print(
        "stage155_156_stress.py is retired; delegating to "
        "stage155_156_fixedpoint_stress.py"
    )
    forwarded_argv = [str(target), *argv[1:]]
    return int(module.main(forwarded_argv))


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
