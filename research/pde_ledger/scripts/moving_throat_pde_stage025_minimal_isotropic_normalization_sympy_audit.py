#!/usr/bin/env python3
"""
moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py

SymPy audit for Stage 25 of the moving-throat PDE program.

Scope
-----
This script verifies the minimal isotropic single-support/single-port closure that
follows after the Stage-024 angular theorem:

  • exact zero-frequency coefficients B0, Z0, N0,
  • exact closed formula for P0 = N0 / D0,
  • exact target equation on the minimal isotropic branch,
  • stability-domain positivity structure,
  • and exact monotonic derivatives with respect to bare stiffness K and
    support-softening X = C^2/varpi^2.

This is still a reduced-sector theorem, but it turns the remaining normalization
problem into one explicit scalar equation in radial/axial overlap amplitudes.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

K, varpi = sp.symbols("K varpi", positive=True, real=True)
C = sp.symbols("C", real=True)
OmegaU, OmegaW = sp.symbols("OmegaU OmegaW", positive=True, real=True)
GU, GW, R = sp.symbols("GU GW R", real=True)
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
mhat = sp.symbols("mhat", positive=True, real=True)

X = sp.symbols("X", nonnegative=True, real=True)

SAMPLE_POINT = {K: sp.Integer(2), varpi: sp.Integer(1), C: sp.Integer(1),
                OmegaU: sp.Integer(2), OmegaW: sp.Integer(2), R: sp.Integer(1),
                GU: sp.Integer(1), GW: sp.Integer(1)}


# ---------------------------------------------------------------------------
# I. Minimal isotropic zero-frequency coefficients
# ---------------------------------------------------------------------------

def zero_frequency_coefficients() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION I — MINIMAL ISOTROPIC ZERO-FREQUENCY COEFFICIENTS")

    Delta = sp.simplify(OmegaU**2 * OmegaW**2 - R**2)
    Q = sp.simplify(GU**2 * OmegaW**2 + 2 * GU * GW * R + GW**2 * OmegaU**2)
    P = sp.simplify(OmegaU**2 * GW + R * GU)

    B0 = sp.simplify(C**2 / varpi**2)
    Z0 = sp.simplify(Q / Delta)
    N0 = sp.simplify(P**2 / Delta**2)
    D0 = sp.simplify(K - B0 - Z0)

    print("Delta =", Delta)
    print("Q     =", Q)
    print("P     =", P)
    print("B0    =", B0)
    print("Z0    =", Z0)
    print("N0    =", N0)
    print("D0    =", D0)

    return Delta, Q, P, B0, Z0, N0, D0


# ---------------------------------------------------------------------------
# II. Exact normalization formula
# ---------------------------------------------------------------------------

def normalization_formula() -> None:
    banner("SECTION II — EXACT MINIMAL ISOTROPIC NORMALIZATION FORMULA")

    Delta, Q, P, B0, Z0, N0, D0 = zero_frequency_coefficients()
    P0 = sp.simplify(N0 / D0)
    P0_compact = sp.simplify(P**2 / (Delta * (K * Delta - Delta * C**2 / varpi**2 - Q)))

    subbanner("II.1 — P0 in raw and compact form")
    print("P0 raw    =", P0)
    print("P0 compact=", P0_compact)
    expect_zero("P0 - P0_compact", P0 - P0_compact)

    p0_value_raw = sp.nsimplify(P0.subs(SAMPLE_POINT))
    p0_value_compact = sp.nsimplify(P0_compact.subs(SAMPLE_POINT))
    print(f"P0 raw at sample     = {p0_value_raw}")
    print(f"P0 compact at sample = {p0_value_compact}")
    if p0_value_raw != sp.Rational(1, 3):
        raise AssertionError(f"P0 raw at sample point != 1/3: got {p0_value_raw}")
    if p0_value_compact != sp.Rational(1, 3):
        raise AssertionError(f"P0 compact at sample point != 1/3: got {p0_value_compact}")

    # target coefficient 54/5 is carried forward from
    # scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:321-342,
    # where Gamma5_port = a^5/(27*c_s^5) and gamma_GR = 2*G/(5*c^5)
    # imply mhat^2*P0 = gamma_GR/Gamma5_port = 54*G*c_s^5/(5*a^5*c^5).
    target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
    equation_residual = sp.simplify(mhat**2 * P0_compact - target)

    subbanner("II.2 — Exact target equation")
    print("Target residual =")
    sp.pprint(equation_residual)
    # Solvability check: mhat^2 = target / P0_compact must be positive on the
    # stability branch (Delta > 0, D0 > 0, so P0_compact > 0; and target > 0
    # since G, c_s, a, c > 0). On the sample point this is a definite rational.
    mhat_sq = sp.simplify(target / P0_compact)
    sample_with_units = dict(SAMPLE_POINT)
    sample_with_units.update({G: sp.Integer(1), c_s: sp.Integer(1), a: sp.Integer(1), c: sp.Integer(1)})
    mhat_sq_at_sample = sp.nsimplify(mhat_sq.subs(sample_with_units))
    print(f"mhat^2 on sample = {mhat_sq_at_sample}")
    if mhat_sq_at_sample <= 0:
        raise AssertionError(f"mhat^2 on sample is not positive: {mhat_sq_at_sample}")

    subbanner("II.3 — P = 0 corollary (paper Checks item 4)")
    # Forcing GW = -R*GU/OmegaU^2 gives P = OmegaU^2*GW + R*GU = 0 symbolically.
    P_zero_sub = {GW: -R * GU / OmegaU**2}
    N0_at_Pzero = sp.simplify(N0.subs(P_zero_sub))
    print(f"N0 at P=0 = {N0_at_Pzero}")
    expect_zero("N0 vanishes when P=0", N0_at_Pzero)
    residual_at_Pzero = sp.simplify((mhat**2 * P0_compact - target).subs(P_zero_sub))
    print(f"(mhat^2*P0 - target) at P=0 = {residual_at_Pzero}")
    expect_zero("(mhat^2*P0 - target) at P=0 equals -target", residual_at_Pzero + target)


# ---------------------------------------------------------------------------
# III. Stability and positivity structure
# ---------------------------------------------------------------------------

def stability_and_positivity() -> None:
    banner("SECTION III — STABILITY AND POSITIVITY STRUCTURE")

    Delta, Q, P, B0, Z0, N0, D0 = zero_frequency_coefficients()
    compact_denom = sp.simplify(Delta * D0)
    expect_zero("Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", compact_denom - (K * Delta - Delta * C**2 / varpi**2 - Q))
    P_raw = OmegaU**2 * GW + R * GU
    Delta_raw = OmegaU**2 * OmegaW**2 - R**2
    expect_zero("N0 reconstructed from raw symbols", N0 - P_raw**2 / Delta_raw**2)

    delta_value = sp.nsimplify(Delta.subs(SAMPLE_POINT))
    d0_value = sp.nsimplify(D0.subs(SAMPLE_POINT))
    p0_pos = sp.nsimplify((N0 / D0).subs(SAMPLE_POINT))
    print(f"Delta on sample = {delta_value}")
    print(f"D0    on sample = {d0_value}")
    print(f"P0    on sample = {p0_pos}")
    if delta_value <= 0:
        raise AssertionError(f"Delta on sample is not positive: {delta_value}")
    if d0_value <= 0:
        raise AssertionError(f"D0 on sample is not positive: {d0_value}")
    if p0_pos <= 0:
        raise AssertionError(f"P0 on sample is not positive: {p0_pos}")

    print("If Delta > 0 and D0 > 0, then P0 > 0 whenever P != 0.")


# ---------------------------------------------------------------------------
# IV. Exact monotonic derivatives on the stable branch
# ---------------------------------------------------------------------------

def monotonic_derivatives() -> None:
    banner("SECTION IV — EXACT MONOTONIC DERIVATIVES")

    Delta, Q, P, B0, Z0, N0, D0 = zero_frequency_coefficients()
    P0 = sp.simplify(N0 / (K - X - Q / Delta))

    dP0_dK = sp.simplify(sp.diff(P0, K))
    dP0_dX = sp.simplify(sp.diff(P0, X))

    subbanner("IV.1 — Derivatives with respect to K and X = C^2/varpi^2")
    print("dP0/dK =", dP0_dK)
    print("dP0/dX =", dP0_dX)

    expect_zero("dP0/dK + N0/(K - X - Q/Delta)^2", dP0_dK + N0 / (K - X - Q / Delta)**2)
    expect_zero("dP0/dX - N0/(K - X - Q/Delta)^2", dP0_dX - N0 / (K - X - Q / Delta)**2)
    expect_zero("dP0/dX + dP0/dK", dP0_dX + dP0_dK)

    sample_iv = dict(SAMPLE_POINT)
    sample_iv[X] = sp.Integer(1)  # X = C^2/varpi^2 = 1 at sample
    dP0_dK_val = sp.nsimplify(dP0_dK.subs(sample_iv))
    dP0_dX_val = sp.nsimplify(dP0_dX.subs(sample_iv))
    print(f"dP0/dK on sample = {dP0_dK_val}")
    print(f"dP0/dX on sample = {dP0_dX_val}")
    if dP0_dK_val >= 0:
        raise AssertionError(f"dP0/dK on sample is not negative: {dP0_dK_val}")
    if dP0_dX_val <= 0:
        raise AssertionError(f"dP0/dX on sample is not positive: {dP0_dX_val}")
    if dP0_dK_val + dP0_dX_val != 0:
        raise AssertionError(f"dP0/dK + dP0/dX on sample is not zero: {dP0_dK_val + dP0_dX_val}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    normalization_formula()
    stability_and_positivity()
    monotonic_derivatives()
