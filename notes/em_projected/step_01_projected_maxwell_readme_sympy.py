#!/usr/bin/env python3
"""Bundle-level audit for step_01_projected_maxwell_readme.md.

The step-01 file is an index note for the first projected-Maxwell bundle rather
than a single derivation. This audit verifies the core formulas summarized
there and checks that the three underlying derivation scripts are present.
"""
from __future__ import annotations

from pathlib import Path

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.simplify(expr)
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def main() -> None:
    root = Path(__file__).resolve().parent
    expected = [
        "step_02_projected_maxwell_covariant_sympy.py",
        "step_03_projected_maxwell_vector_sympy.py",
        "step_04_projection_reduction_comparison_sympy.py",
    ]
    missing = [name for name in expected if not (root / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing bundle scripts: {missing}")

    w = sp.symbols("w", real=True)
    W = sp.Function("W")(w)
    Q = sp.Function("Q")(w)

    # Projection by parts identity:
    # int W Q_w dw = [WQ] - int W_w Q dw.
    by_parts_density = sp.diff(W * Q, w) - W * sp.diff(Q, w) - sp.diff(W, w) * Q
    assert_zero("projection integration-by-parts density", by_parts_density)

    t, x, y, z = sp.symbols("t x y z", real=True)
    E1 = sp.Function("E1")(t, x, y, z)
    E2 = sp.Function("E2")(t, x, y, z)
    E3 = sp.Function("E3")(t, x, y, z)
    B1 = sp.Function("B1")(t, x, y, z)
    B2 = sp.Function("B2")(t, x, y, z)
    B3 = sp.Function("B3")(t, x, y, z)

    F23, F31, F12 = B1, B2, B3
    F10, F20, F30 = E1, E2, E3
    F01, F02, F03 = -E1, -E2, -E3

    faraday = [
        sp.diff(F23, t) + sp.diff(F30, y) + sp.diff(F02, z)
        - (sp.diff(B1, t) + sp.diff(E3, y) - sp.diff(E2, z)),
        sp.diff(F31, t) + sp.diff(F10, z) + sp.diff(F03, x)
        - (sp.diff(B2, t) + sp.diff(E1, z) - sp.diff(E3, x)),
        sp.diff(F12, t) + sp.diff(F20, x) + sp.diff(F01, y)
        - (sp.diff(B3, t) + sp.diff(E2, x) - sp.diff(E1, y)),
    ]
    for i, residue in enumerate(faraday, start=1):
        assert_zero(f"Faraday component {i}", residue)

    lam, mu0 = sp.symbols("lambda mu0", positive=True, nonzero=True)
    Z = sp.exp(-w**2 / lam**2)
    Z_int = sp.simplify(sp.integrate(Z, (w, -sp.oo, sp.oo)))
    Z2_int = sp.simplify(sp.integrate(Z**2, (w, -sp.oo, sp.oo)))
    W_match = sp.simplify(Z / Z_int)
    I_WZ = sp.simplify(sp.integrate(W_match * Z, (w, -sp.oo, sp.oo)))
    I_WS_delta = sp.simplify(W_match.subs(w, 0))
    mu_proj_delta = sp.simplify(mu0 * I_WS_delta / I_WZ)
    mu_red = sp.simplify(mu0 / Z_int)

    assert_zero("Gaussian Z_int", Z_int - sp.sqrt(sp.pi) * lam)
    assert_zero("Gaussian Z2_int", Z2_int - sp.sqrt(2 * sp.pi) * lam / 2)
    assert_zero("matched-kernel I_WZ", I_WZ - sp.sqrt(2) / 2)
    assert_zero("delta-source projection/reduction ratio", mu_proj_delta / mu_red - sp.sqrt(2))

    print("STEP 01 PROJECTED MAXWELL README AUDIT")
    print("Checked script inventory, projection identity, vector Bianchi signs, and Gaussian coupling mismatch.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
