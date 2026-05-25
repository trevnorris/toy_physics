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
        "moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py",
        "moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py",
        "moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py",
    ]
    missing = [name for name in expected if not (root / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing bundle scripts: {missing}")

    w = sp.symbols("w", real=True)

    # Projection by parts identity:
    # Verify integration-by-parts at density level using concrete decaying
    # test functions so the boundary term [W Q] vanishes at +/- infinity:
    #   int W Q_w dw + int W_w Q dw = 0.
    lam_ibp = sp.Symbol("lam_ibp", positive=True)
    W_ex = sp.exp(-w**2 / lam_ibp**2)
    Q_ex = w * sp.exp(-w**2 / lam_ibp**2)
    ibp_lhs = sp.integrate(W_ex * sp.diff(Q_ex, w), (w, -sp.oo, sp.oo))
    ibp_rhs = -sp.integrate(sp.diff(W_ex, w) * Q_ex, (w, -sp.oo, sp.oo))
    assert_zero("projection integration-by-parts (decaying test functions)",
                sp.simplify(ibp_lhs - ibp_rhs))

    t, x, y, z = sp.symbols("t x y z", real=True)
    coords = (t, x, y, z)
    A0 = sp.Function("A0")(t, x, y, z)
    A1 = sp.Function("A1")(t, x, y, z)
    A2 = sp.Function("A2")(t, x, y, z)
    A3 = sp.Function("A3")(t, x, y, z)
    A_components = (A0, A1, A2, A3)

    # F_{mu nu} = d_mu A_nu - d_nu A_mu (antisymmetric by construction).
    def F(mu, nu):
        return (sp.diff(A_components[nu], coords[mu])
                - sp.diff(A_components[mu], coords[nu]))

    # Cyclic Bianchi: d_[alpha F_{beta gamma]} = 0 for F = dA, by Schwarz
    # symmetry of mixed partials on smooth A on real coords. The check is
    # non-tautological: a sign error in the F definition produces a nonzero
    # residual.
    for (alpha, beta, gamma) in [(0, 2, 3), (0, 3, 1), (0, 1, 2)]:
        cyc = (sp.diff(F(beta, gamma), coords[alpha])
               + sp.diff(F(gamma, alpha), coords[beta])
               + sp.diff(F(alpha, beta), coords[gamma]))
        assert_zero(
            f"cyclic Bianchi (alpha={alpha}, beta={beta}, gamma={gamma})", cyc)

    # Specialize via E_i = -F_{0i} and B_1 = F_{23}, B_2 = F_{31}, B_3 = F_{12},
    # then verify the three components of dB/dt + curl(E) = 0 reduce to
    # cyclic Bianchi (and hence vanish). A sign error in the E,B<->F map
    # produces a nonzero residual.
    E_from_A = (-F(0, 1), -F(0, 2), -F(0, 3))
    B_from_A = (F(2, 3), F(3, 1), F(1, 2))
    mf1 = (sp.diff(B_from_A[0], t)
           + sp.diff(E_from_A[2], y) - sp.diff(E_from_A[1], z))
    mf2 = (sp.diff(B_from_A[1], t)
           + sp.diff(E_from_A[0], z) - sp.diff(E_from_A[2], x))
    mf3 = (sp.diff(B_from_A[2], t)
           + sp.diff(E_from_A[1], x) - sp.diff(E_from_A[0], y))
    assert_zero("Maxwell-Faraday component 1 from A", mf1)
    assert_zero("Maxwell-Faraday component 2 from A", mf2)
    assert_zero("Maxwell-Faraday component 3 from A", mf3)

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
    print("Checked script inventory, projection identity, cyclic Bianchi from F=dA and Maxwell-Faraday reduction, and Gaussian coupling mismatch.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
