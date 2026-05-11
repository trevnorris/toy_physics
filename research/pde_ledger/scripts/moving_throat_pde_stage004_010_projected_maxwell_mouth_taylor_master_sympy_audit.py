#!/usr/bin/env python3
"""Master-note audit for step_11_projected_maxwell_mouth_taylor_master_notes.md."""
from __future__ import annotations

import sympy as sp


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def main() -> None:
    u, ell = sp.symbols("u ell", positive=True)
    X0, X1, X2 = sp.symbols("X0 X1 X2")
    W = sp.exp(-u)
    X = X0 + ell * u * X1 + ell**2 * u**2 * X2 / 2
    Xproj = sp.integrate(W * X, (u, 0, sp.oo))
    assert_zero("one-sided Taylor first moment", sp.series(Xproj, ell, 0, 2).removeO() - (X0 + ell * X1))
    assert_nonzero(
        "mutated one-sided Taylor first moment should fail",
        sp.series(Xproj, ell, 0, 2).removeO() - (X0 + 2 * ell * X1),
    )
    W2 = u * sp.exp(-u)
    mu1_W2 = sp.integrate(u * W2, (u, 0, sp.oo))
    Xproj_W2 = sp.integrate(W2 * X, (u, 0, sp.oo))
    assert_zero("one-sided Taylor W2 first moment", mu1_W2 - 2)
    assert_zero("one-sided Taylor arbitrary first moment", sp.series(Xproj_W2, ell, 0, 2).removeO() - (X0 + ell * mu1_W2 * X1))
    assert_nonzero(
        "mutated arbitrary first moment should fail",
        sp.series(Xproj_W2, ell, 0, 2).removeO() - (X0 + ell * X1),
    )

    Q, S2, Hport, Delta, P, Gw = sp.symbols("Q S2 H Delta P Gw", nonzero=True)
    D0, D2, D4, N0, Ptarget = sp.symbols("D0 D2 D4 N0 Ptarget", nonzero=True)
    mu1 = sp.symbols("mu1", nonzero=True)

    q1, s1, h1, d1, p1, g1 = sp.symbols("q1 s1 h1 d1 p1 g1")
    Qx, Sx, Hx, Dx, Px, Gx = sp.symbols("Qx S2x Hx Deltax Px Gwx")
    subs_der = {q1: mu1 * Qx, s1: mu1 * Sx, h1: mu1 * Hx, d1: mu1 * Dx, p1: mu1 * Px, g1: mu1 * Gx}

    z0 = (Delta * q1 - Q * d1) / Delta**2
    z2 = (-Delta**2 * h1 + Delta * (Hport * d1 + Q * s1 + S2 * q1) - 2 * Q * S2 * d1) / Delta**3
    z4 = (
        -Delta**2 * Hport * s1
        - Delta**2 * S2 * h1
        - Delta**2 * q1
        + 2 * Delta * Hport * S2 * d1
        + 2 * Delta * Q * S2 * s1
        + 2 * Delta * Q * d1
        + Delta * S2**2 * q1
        - 3 * Q * S2**2 * d1
    ) / Delta**4

    n0 = 2 * P * (Delta * p1 - P * d1) / Delta**3
    n2 = -(
        2 * Delta**2 * (Gw * p1 + P * g1)
        - 2 * Delta * P * (2 * Gw * d1 + P * s1 + 2 * S2 * p1)
        + 6 * P**2 * S2 * d1
    ) / Delta**4
    n4 = 2 * (
        Delta**3 * Gw * g1
        - Delta**2 * Gw**2 * d1
        - 2 * Delta**2 * Gw * P * s1
        - 2 * Delta**2 * Gw * S2 * p1
        - 2 * Delta**2 * P * S2 * g1
        - 2 * Delta**2 * P * p1
        + 6 * Delta * Gw * P * S2 * d1
        + 3 * Delta * P**2 * S2 * s1
        + 3 * Delta * P**2 * d1
        + 3 * Delta * P * S2**2 * p1
        - 6 * P**2 * S2**2 * d1
    ) / Delta**5

    Xi = sp.simplify((2 * p1 / P - 2 * d1 / Delta + q1 / (D0 * Delta) - Q * d1 / (D0 * Delta**2)).subs(subs_der) / mu1)
    K1 = sp.simplify((-(z2 + z0 / sp.Integer(9))).subs(subs_der) / mu1)
    H_even = sp.simplify(((-z4 + sp.Rational(2, 3) * z2 - z0 / sp.Integer(27))).subs(subs_der) / mu1)

    for sym in (Sx, Hx, Gx):
        assert_zero(f"Xi independence from {sym}", sp.diff(Xi, sym))
    for sym in (Px, Gx):
        assert_zero(f"even gates independence from {sym}", sp.diff(K1, sym))
        assert_zero(f"even gates independence from {sym}", sp.diff(H_even, sym))
    assert_zero("dXi/dPprime", sp.diff(Xi, Px) - 2 / P)

    deltaP2 = sp.simplify((D0**2 * n2 - 2 * D0 * D2 * n0 + 2 * D0 * N0 * z2 - 2 * D2 * N0 * z0) / D0**3)
    deltaP4 = sp.simplify(
        (
            D0**3 * n4
            - 2 * D0**2 * D2 * n2
            - 2 * D0**2 * D4 * n0
            + 2 * D0**2 * N0 * z4
            + 3 * D0 * D2**2 * n0
            - 2 * D0 * D2 * N0 * z2
            - 2 * D0 * D4 * N0 * z0
            + 2 * D2**2 * N0 * z0
        )
        / D0**4
    )
    deltaP2_der = sp.simplify(deltaP2.subs(subs_der) / mu1)
    deltaP4_der = sp.simplify(deltaP4.subs(subs_der) / mu1)
    assert_zero("d(delta P2)/dGprime", sp.diff(deltaP2_der, Gx) + 2 * P / (D0 * Delta**2))
    if sp.simplify(sp.diff(deltaP4_der, Gx)) == 0:
        raise AssertionError("delta P4 should depend on G_W prime")

    qd_only = sp.solve(
        [sp.Eq(K1.subs({Sx: 0, Hx: 0}), 0), sp.Eq(H_even.subs({Sx: 0, Hx: 0}), 0)],
        [Qx, Dx],
        dict=True,
    )
    sh_only = sp.solve(
        [sp.Eq(K1.subs({Qx: 0, Dx: 0}), 0), sp.Eq(H_even.subs({Qx: 0, Dx: 0}), 0)],
        [Sx, Hx],
        dict=True,
    )
    qd_matrix = sp.Matrix([
        [sp.diff(K1.subs({Sx: 0, Hx: 0}), Qx), sp.diff(K1.subs({Sx: 0, Hx: 0}), Dx)],
        [sp.diff(H_even.subs({Sx: 0, Hx: 0}), Qx), sp.diff(H_even.subs({Sx: 0, Hx: 0}), Dx)],
    ])
    sh_matrix = sp.Matrix([
        [sp.diff(K1.subs({Qx: 0, Dx: 0}), Sx), sp.diff(K1.subs({Qx: 0, Dx: 0}), Hx)],
        [sp.diff(H_even.subs({Qx: 0, Dx: 0}), Sx), sp.diff(H_even.subs({Qx: 0, Dx: 0}), Hx)],
    ])
    assert_nonzero("Xi should depend on Pprime", sp.diff(Xi, Px))
    assert_nonzero("deltaP4 should depend on G_W prime", sp.diff(deltaP4_der, Gx))
    assert_nonzero("source/denominator sieve determinant", qd_matrix.det())
    assert_nonzero("spectral sieve determinant", sh_matrix.det())
    if qd_only != [{Qx: 0, Dx: 0}]:
        raise AssertionError(f"Unexpected pure source/denominator solve: {qd_only}")
    if sh_only != [{Sx: 0, Hx: 0}]:
        raise AssertionError(f"Unexpected pure spectral solve: {sh_only}")

    print("STEP 11 PROJECTED MAXWELL MOUTH-TAYLOR MASTER AUDIT")
    print("Checked one-sided Taylor projection, bottleneck dependencies, G_W transport entry, and mechanism sieve.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
