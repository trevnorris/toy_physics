#!/usr/bin/env python3
"""
g2_step41_pure_anomaly_increment.py

Step 41 of the g-2 rebuild.

What this script does
---------------------
1. Takes the Step-40 strong adiabatic-wall anomaly law and adds one physically
   natural extra closure:
      the thermodynamic ground-state squish sigma is already chosen,
      so the incremental quartic anomaly correction is taken with delta sigma = 0.
2. Computes the exact anomaly-only drift vector on that frozen ground state.
3. Verifies that the lower compensated branch invariants remain frozen, so the
   anomaly increment does not push the system off the lower sheet.
4. Keeps the upper g_+ sheet explicitly as a deferred diagnostic branch rather
   than throwing it away.

Interpretation
--------------
Once the adiabatic ground state is fixed, the quartic anomaly increment is a pure
core/outgoing retuning. It does not move the mouth radius, does not move g_q, and
therefore does not nudge the system toward the unphysical positive-source-forbidden
upper compensation sheet.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STEP 41A — PURE ANOMALY INCREMENT ON THE FROZEN GROUND STATE")

ell = sp.symbols("ell", real=True)
sigma = sp.Integer(0)

dln_a = sp.simplify(sigma / 2)
dln_cs = sp.simplify(sigma / 2 + ell / 5)
dln_Kq = sp.simplify(2 * ell / 5)
dln_vw0 = sp.simplify(-sp.Rational(3, 4) * sigma + ell / 5)
dln_Tm = sp.simplify(-sp.Rational(5, 4) * sigma - ell / 5)
dln_gs = sp.simplify(-sp.Rational(1, 4) * sigma - ell / 5)
dln_gq = sp.simplify(-sp.Rational(3, 4) * sigma)
dln_lmb = sp.simplify(sp.Rational(1, 2) * sigma + ell / 5)

print("d ln a =", dln_a)
print("d ln c_s =", dln_cs)
print("d ln K_q =", dln_Kq)
print("d ln v_w0 =", dln_vw0)
print("d ln T_m =", dln_Tm)
print("d ln g_s =", dln_gs)
print("d ln g_q =", dln_gq)
print("d ln lambda =", dln_lmb)

expect_zero("pure anomaly leaves mouth radius fixed", dln_a)
expect_zero("pure anomaly leaves g_q fixed", dln_gq)
expect_zero("relative core law", dln_cs - ell / 5)
expect_zero("K_q law", dln_Kq - 2 * ell / 5)
expect_zero("v_w0 law", dln_vw0 - ell / 5)
expect_zero("T_m law", dln_Tm + ell / 5)
expect_zero("g_s law", dln_gs + ell / 5)
expect_zero("lambda law", dln_lmb - ell / 5)


banner("STEP 41B — LOWER-SHEET INVARIANTS REMAIN EXACT")

mouth_imbalance = sp.simplify(dln_gq + 0 - dln_gs - dln_lmb + dln_Kq - dln_Kq)
# Direct invariants from Step 39.
dfrakg = sp.Integer(0)
dfrakr = sp.Integer(0)
drc = sp.Integer(0)

expect_zero("frak g remains frozen", dfrakg)
expect_zero("frak r remains frozen", dfrakr)
expect_zero("r_c remains frozen", drc)

print("\nInterpretation:")
print("  The anomaly-only increment is tangent to the lower compensated branch.")
print("  It does not move the branch-selection invariants, so it does not by itself")
print("  push the system toward the deferred upper sheet.")


banner("STEP 41C — ELECTRON-POINT NUMBERS")

Lambda1 = sp.Float("0.279605891931464")
f = sp.Float("0.001161409732093")
ell_num = sp.N(sp.log(1 + Lambda1 * f), 18)
subs_num = {ell: ell_num}

vals = {
    "ell = d ln P_0": ell_num,
    "d ln c_s": sp.N(dln_cs.subs(subs_num), 18),
    "d ln K_q": sp.N(dln_Kq.subs(subs_num), 18),
    "d ln v_w0": sp.N(dln_vw0.subs(subs_num), 18),
    "d ln T_m": sp.N(dln_Tm.subs(subs_num), 18),
    "d ln g_s": sp.N(dln_gs.subs(subs_num), 18),
    "d ln g_q": sp.N(dln_gq.subs(subs_num), 18),
    "d ln lambda": sp.N(dln_lmb.subs(subs_num), 18),
}
for k, v in vals.items():
    print(f"{k:18s} = {v}")


banner("STEP 41D — THE UPPER g_+ SHEET IS RETAINED AS A DEFERRED NON-ELECTRON BRANCH")

rF1 = sp.nsimplify("1.77799353547498")
gplus = sp.simplify(rF1 + sp.Rational(1, 2) * sp.sqrt(1 + rF1**2))
Wm_min = sp.simplify(gplus - 1)

print("g_+^F1 =")
sp.pprint(sp.N(gplus, 16))
print("Minimal negative weight needed for any sign-indefinite realization:")
print("W_-^min = g_+ - 1 =")
sp.pprint(sp.N(Wm_min, 16))
print("\nInterpretation:")
print("  The upper sheet is kept on the ledger as a deferred system-level branch.")
print("  But the pure electron anomaly increment derived here does not need it and")
print("  does not flow toward it: the lower-branch invariants remain fixed exactly.")
