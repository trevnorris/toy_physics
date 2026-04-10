#!/usr/bin/env python3
"""
g2_step53_canonical_outgoing_dtn_closure.py

Step 53 of the g-2 chain.

What this script does
---------------------
1. Re-derives the canonical compact outgoing l=2 DtN fingerprint in the z=a*omega/c_s variable.
2. Matches that fingerprint to the minimal isotropic grouped-P2 retarded module.
3. Fixes the outgoing-normalization factor chi_Q exactly to 1 on the canonical DtN branch.
4. Inserts the natural source-map relation mhat_0 -> 1 and shows that the conservative
   quadrupole normalization defect collapses to N_Q = 1 on that branch.
"""

from __future__ import annotations

import sympy as sp

I = sp.I


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


banner("STEP 53 — CANONICAL OUTGOING DTN CLOSURE")

# ---------------------------------------------------------------------------
# 1. Exact compact outgoing l=2 DtN fingerprint in the low-frequency variable z
# ---------------------------------------------------------------------------

z = sp.symbols("z", real=True)
Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)
Y_out = sp.simplify(-3 / Lambda_out)
Y_out_series = sp.expand(sp.series(Y_out, z, 0, 6).removeO())

print("Lambda_out(z) =", Lambda_out)
print("Y_out(z)      =", Y_out_series)

expected_Y_out = 1 + z**2/sp.Integer(9) + 4*z**4/sp.Integer(81) + I*z**5/sp.Integer(27)
expect_zero("Y_out - expected", Y_out_series - expected_Y_out)

# ---------------------------------------------------------------------------
# 2. Minimal isotropic grouped-P2 retarded module and exact chi_Q matching
# ---------------------------------------------------------------------------

chi_Q = sp.symbols("chi_Q", real=True)

# With Omega_Q = 3 c_s / (2 a), one has omega^2/Omega_Q^2 = 4 z^2 / 9 and
# sigma_Q^can * omega^5 = 4 z^5 / 27.
Y_ret = sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - 4*z**2/sp.Integer(9) - I*chi_Q*4*z**5/sp.Integer(27))
Y_ret_series = sp.expand(sp.series(Y_ret, z, 0, 6).removeO())

print("Y_ret(z) =", Y_ret_series)

coeff_chi = sp.simplify(sp.expand(Y_ret_series).coeff(z, 5) / I)
print("odd z^5 coefficient in Y_ret =", coeff_chi)

chi_sol = sp.solve(sp.Eq(coeff_chi, sp.Rational(1, 27)), chi_Q)[0]
print("chi_Q =", chi_sol)
if chi_sol != 1:
    raise AssertionError("Canonical outgoing DtN branch did not fix chi_Q = 1.")

expect_zero("Y_ret(chi_Q=1) - Y_out", Y_ret_series.subs(chi_Q, chi_sol) - Y_out_series)

# ---------------------------------------------------------------------------
# 3. Natural source-map branch and conservative normalization defect
# ---------------------------------------------------------------------------

mhat0, N_Q = sp.symbols("mhat0 N_Q", positive=True, real=True)

# Reduced normalization stack from the moving-throat bridge:
# mhat0^2 * chi_Q * N_Q = 1.
NQ_sol = sp.solve(sp.Eq(mhat0**2 * chi_sol * N_Q, 1), N_Q)[0]
print("N_Q(mhat0) =", NQ_sol)

# Strict point-particle natural source-map limit.
NQ_pp = sp.simplify(NQ_sol.subs(mhat0, 1))
print("N_Q(point-particle) =", NQ_pp)
if NQ_pp != 1:
    raise AssertionError("Natural source-map point-particle limit did not give N_Q = 1.")

# ---------------------------------------------------------------------------
# 4. Canonical invariant coefficients on the no-tuning outgoing branch
# ---------------------------------------------------------------------------

G, c, a, c_s = sp.symbols("G c a c_s", positive=True, real=True)

K0bar = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
K2bar = sp.simplify(6 * G * c_s**3 / (5 * a**3 * c**5))
K4bar = sp.simplify(8 * G * c_s / (15 * a * c**5))
Gamma5bar = sp.simplify(2 * G / (5 * c**5))

print("K0bar     =", K0bar)
print("K2bar     =", K2bar)
print("K4bar     =", K4bar)
print("Gamma5bar =", Gamma5bar)

banner("FINAL LEDGER")
print("Canonical compact outgoing l=2 DtN branch fixes chi_Q = 1 exactly.")
print("On the natural point-particle source-map branch mhat_0 -> 1,")
print("the reduced normalization stack collapses to N_Q = 1.")
print("So the canonical outgoing grouped-P2 branch is a true no-tuning closure.")
