#!/usr/bin/env python3
"""Adapter audit for Stage 004.018 reduced one-port normal form."""

from __future__ import annotations

from pathlib import Path
import runpy


if __name__ == "__main__":
    target = Path(__file__).with_name("moving_throat_pde_stage004_maxwell_mixed_sympy_audit.py")
    runpy.run_path(str(target), run_name="__main__")
