#!/usr/bin/env python3
"""
g2_step51_canonical_outgoing_notuning.py

Step 51 of the g-2 rebuild.

What this script does
---------------------
1. Uses the exact outgoing spherical Hankel l=2 DtN branch to derive the
   canonical compact outgoing fingerprint
       Yhat_2^out = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + ...
   directly.
2. Matches that explicit branch to the minimal retarded grouped-P2 module and
   proves that the canonical outgoing normalization factor is fixed exactly to
       chi_Q = 1.
3. Combines this with the natural source-map factorization
       mhat_0^2 chi_Q N_Q = 1,
   so that on the strict point-particle source-map branch one gets the no-tuning
   closure
       N_Q = 1.
4. Shows that pure scale deformations and pure scale+argument deformations cannot
   move chi_Q once the even moments are kept canonical.

Interpretation
--------------
There is now a clear distinction between two questions:
  - the canonical passive/outgoing branch itself is fixed with no tuning;
  - any nonzero electron anomaly must therefore come from a genuine isotropic
    throat-core DtN deformation away from that canonical branch.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STEP 51A — EXACT OUTGOING l=2 DtN FINGERPRINT")

z = sp.symbols("z", real=True)
# Exact closed form for the spherical Hankel l=2 mode.
h2 = sp.exp(sp.I * z) * (z**2 + 3*sp.I*z - 3) / z**3
Lambda_out = sp.simplify(z * sp.diff(sp.log(h2), z))
Lambda_series = sp.expand(sp.series(Lambda_out, z, 0, 8).removeO())
Yhat_out = sp.simplify(-sp.Integer(3) / Lambda_out)
Yhat_series = sp.expand(sp.series(Yhat_out, z, 0, 8).removeO())
Yhat_series_6 = sp.expand(sp.series(Yhat_out, z, 0, 6).removeO())

print("Lambda_2^out(z) =", Lambda_out)
print("Lambda_2^out series =", Lambda_series)
print("Yhat_2^out(z) =", Yhat_out)
print("Yhat_2^out series =", Yhat_series)

expect_zero(
    "canonical z^2 coefficient",
    sp.simplify(Yhat_series.coeff(z, 2) - sp.Rational(1, 9)),
)
expect_zero(
    "canonical z^4 coefficient",
    sp.simplify(Yhat_series.coeff(z, 4) - sp.Rational(4, 81)),
)
expect_zero(
    "canonical odd coefficient",
    sp.simplify(sp.im(Yhat_series.coeff(z, 5)) - sp.Rational(1, 27)),
)

print("\nReading:")
print("  The explicit outgoing spherical l=2 DtN branch fixes the canonical odd")
print("  quadrupole coefficient directly. It is not a free reduced parameter.")


banner("STEP 51B — MATCH TO THE MINIMAL RETARDED GROUPED-P2 MODULE")

omega, a, c_s, OmegaQ, chi_Q = sp.symbols("omega a c_s Omega_Q chi_Q", positive=True, real=True)

sigma_can = sp.simplify(9 / (8 * OmegaQ**5))
sigma_can_geom = sp.simplify(sigma_can.subs(OmegaQ, 3 * c_s / (2 * a)))
Yret = sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ**2 - sp.I * chi_Q * sigma_can * omega**5)
Yret_series = sp.expand(sp.series(Yret, omega, 0, 6).removeO())
Yret_geom = sp.expand(Yret_series.subs(OmegaQ, 3 * c_s / (2 * a)))

print("sigma_can(OmegaQ) =", sigma_can)
print("sigma_can(a,c_s)  =", sigma_can_geom)
print("Yret series       =", Yret_series)
print("Yret geom series  =", Yret_geom)

expect_zero(
    "sigma_can geometry rewrite",
    sp.simplify(sigma_can_geom - 4 * a**5 / (27 * c_s**5)),
)
expect_zero(
    "retarded z^2 coefficient",
    sp.simplify(Yret_geom.coeff(omega, 2) - a**2 / (9 * c_s**2)),
)
expect_zero(
    "retarded z^4 coefficient",
    sp.simplify(Yret_geom.coeff(omega, 4) - 4 * a**4 / (81 * c_s**4)),
)

chi_sol = sp.solve(
    sp.Eq(Yret_geom.coeff(omega, 5), sp.I * a**5 / (27 * c_s**5)),
    chi_Q,
)[0]
print("chi_Q from exact DtN match =", chi_sol)
expect_zero("canonical outgoing normalization", sp.simplify(chi_sol - 1))

print("\nReading:")
print("  The explicit compact outgoing DtN branch fixes the minimal retarded")
print("  grouped-P2 normalization to chi_Q = 1 exactly.")


banner("STEP 51C — NATURAL SOURCE-MAP NO-TUNING CLOSURE")

mhat0, NQ = sp.symbols("mhat_0 N_Q", positive=True, real=True)
factorization = sp.Eq(mhat0**2 * chi_Q * NQ, 1)
NQ_exact = sp.solve(factorization.subs(chi_Q, 1), NQ)[0]
NQ_point = sp.simplify(NQ_exact.subs(mhat0, 1))

print("factorization =", factorization)
print("N_Q exact when chi_Q=1 =", NQ_exact)
print("N_Q on strict source-map branch =", NQ_point)
expect_zero("no-tuning closure", sp.simplify(NQ_point - 1))

print("\nReading:")
print("  On the natural point-particle source-map branch, the canonical outgoing")
print("  DtN model closes the normalization stack with no further parameter fixing:")
print("      chi_Q = 1  ->  N_Q = 1.")


banner("STEP 51D — ROBUSTNESS: PURE SCALE AND PURE ARGUMENT RESCALINGS CANNOT SHIFT THE CANONICAL VALUE")

S, beta = sp.symbols("S beta", positive=True, real=True)
Y_scale = sp.simplify(S * Lambda_out)
Y_scale_norm = sp.simplify((-3 * S) / Y_scale)
Y_scale_series = sp.expand(sp.series(Y_scale_norm, z, 0, 6).removeO())

print("Y_scale_norm series =", Y_scale_series)
expect_zero("pure scale invariance", sp.simplify(Y_scale_series - Yhat_series_6))

Y_beta = sp.expand(sp.series((-3 * S) / (S * Lambda_out.subs(z, beta * z)), z, 0, 6).removeO())
print("Y_beta series =", Y_beta)
expect_zero("beta z^2 coefficient", sp.simplify(Y_beta.coeff(z, 2) - beta**2 / 9))
expect_zero("beta z^4 coefficient", sp.simplify(Y_beta.coeff(z, 4) - 4 * beta**4 / 81))

print("\nReading:")
print("  Pure scale deformations are invisible to the normalized outgoing branch,")
print("  and pure argument rescalings move the even moments unless beta = 1.")
print("  So once the canonical even fingerprint is frozen, those classes do not")
print("  generate a nontrivial chi_Q.")

print("\nFinal reading of Step 51:")
print("  - the canonical compact outgoing l=2 DtN branch is a genuine no-tuning")
print("    closure, not an unfixed phenomenological coefficient, and")
print("  - any nonzero electron anomaly must therefore come from a genuine isotropic")
print("    throat-core DtN deformation away from that canonical branch.")
