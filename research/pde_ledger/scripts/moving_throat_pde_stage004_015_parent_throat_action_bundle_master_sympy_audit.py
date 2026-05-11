#!/usr/bin/env python3
"""Master-note audit for step_16_parent_throat_action_bundle_master_notes.md."""
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
    KSigma, MSigma = sp.symbols("KSigma MSigma", nonzero=True)
    B0, B2, B4, Z0, Z2, Z4 = sp.symbols("B0 B2 B4 Z0 Z2 Z4", nonzero=True)
    N0, Ptarget = sp.symbols("N0 Ptarget", nonzero=True)

    D0 = KSigma - B0 - Z0
    D2 = -(MSigma + B2 + Z2)
    D4 = -(B4 + Z4)
    u2 = -D2 / D0
    u4 = (D2**2 - D0 * D4) / D0**2
    assert_zero(
        "one-pole numerator equivalence",
        (u4 - 4*u2**2) - (D0*(B4 + Z4) - 3*(MSigma + B2 + Z2)**2)/D0**2,
    )

    K_from_one_pole = B0 + Z0 + 3 * (MSigma + B2 + Z2) ** 2 / (B4 + Z4)
    K_from_norm = B0 + Z0 + N0 / Ptarget
    assert_zero("one-pole K solve", (u4 - 4 * u2**2).subs(KSigma, K_from_one_pole))
    assert_zero("normalization K solve", (N0 / D0).subs(KSigma, K_from_norm) - Ptarget)
    assert_nonzero("mutated one-pole K solve should fail", (u4 - 3 * u2**2).subs(KSigma, K_from_one_pole))
    assert_nonzero("mutated normalization K solve should fail", (N0 / D0).subs(KSigma, K_from_norm) - 2 * Ptarget)
    compatibility = N0 / Ptarget - 3 * (MSigma + B2 + Z2) ** 2 / (B4 + Z4)
    assert_zero("compatibility equality", sp.simplify(K_from_norm - K_from_one_pole) - compatibility)

    dKSigma, dMSigma = sp.symbols("dKSigma dMSigma")
    B01, B21, B41 = sp.symbols("B01 B21 B41")
    Z01, Z21, Z41 = sp.symbols("Z01 Z21 Z41")
    N01 = sp.symbols("N01")

    D01 = dKSigma - B01 - Z01
    D21 = -(dMSigma + B21 + Z21)
    D41 = -(B41 + Z41)
    K1 = D21 + D01 / 9
    H_even = D41 - sp.Rational(2, 3) * D21 - D01 / 27
    Xi1 = N01 / N0 - D01 / D0
    coeff_matrix = sp.Matrix([
        [sp.diff(K1, dKSigma), sp.diff(K1, dMSigma)],
        [sp.diff(H_even, dKSigma), sp.diff(H_even, dMSigma)],
    ])
    assert_zero("even-gate solve determinant", coeff_matrix.det() - sp.Rational(1, 27))

    sol = sp.solve([sp.Eq(K1, 0), sp.Eq(H_even, 0)], [dKSigma, dMSigma], dict=True)[0]
    expected_dK = B01 + Z01 + 27 * (B41 + Z41)
    expected_dM = -(B21 + Z21) + 3 * (B41 + Z41)
    assert_zero("wall stiffness slope solve", sol[dKSigma] - expected_dK)
    assert_zero("wall inertia slope solve", sol[dMSigma] - expected_dM)
    assert_nonzero("mutated even-gate determinant should fail", coeff_matrix.det() + sp.Rational(1, 27))

    D01_comp = sp.simplify(D01.subs(sol))
    assert_zero("compensated D01", D01_comp - 27 * (B41 + Z41))
    assert_zero("compensated K1", K1.subs(sol))
    assert_zero("compensated H_even", H_even.subs(sol))
    assert_zero(
        "residual Xi1",
        Xi1.subs(sol) - (N01 / N0 - 27 * (B41 + Z41) / (KSigma - B0 - Z0)),
    )
    assert_nonzero(
        "mutated residual Xi1 should fail",
        Xi1.subs(sol) - (N01 / N0 + 27 * (B41 + Z41) / (KSigma - B0 - Z0)),
    )

    w = sp.symbols("w", real=True)
    beta = sp.exp(-w**2 / 2)
    MSigma_example = sp.integrate(beta**2, (w, -sp.oo, sp.oo))
    KSigma_example = sp.integrate(sp.diff(beta, w)**2 + beta**2, (w, -sp.oo, sp.oo))
    assert_zero("example wall inertia integral", MSigma_example - sp.sqrt(sp.pi))
    assert_zero("example wall stiffness integral", KSigma_example - 3*sp.sqrt(sp.pi)/2)

    print("STEP 16 PARENT THROAT ACTION BUNDLE MASTER AUDIT")
    print("Checked isotropic compatibility, exact weak-axisymmetric wall-slope solve, and residual Xi1.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
