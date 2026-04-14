#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage013_selected_branch_dynamic_window_compiler_sympy_audit.py

Stage 013 — feed the Stage-012 selected-branch numerator/denominator classifier
into the actual Stage-011 rigid-split dynamic responses.

What this script does
---------------------
1. Carries forward the two exact Stage-011 pure-transfer rigid branches:
      - denominator-carried (numerator-rigid),
      - numerator-carried   (denominator-rigid),
   and converts their wall-like dynamic responses to slopes per unit Xi_1.
2. Uses the exact Stage-012 classifier
      R_ND = L_num / L_den
   to build the exact selected-branch share weights
      w_num = R_ND / (1 + R_ND),
      w_den = 1 / (1 + R_ND).
3. Compiles the selected-branch wall-like dynamic log-slopes as exact affine
   mixtures of the rigid-branch per-unit-Xi_1 slopes.
4. Proves the universal sign split:
      - upper wall-like pole always worsens,
      - lower wall-like pole flips sign at one exact classifier threshold R_*.
5. Translates the result into dynamic ceilings for |eps * Xi_1| and compares them
   against the universal transported static ceilings from Stage 007 / 011.

Main structural result
----------------------
Within the exact rigid-split compiler on the concrete compatibility branch, the
selected-branch dynamic window never becomes the first kill condition.
The transported static ceiling remains stricter everywhere on the full Stage-012
classifier map.

In particular:
- every denominator-like point (R_ND <= 1) has an infinite nonempty dynamic ceiling;
- if delta >= 8/9, the whole selected branch is denominator-like, so the nonempty
  dynamic ceiling is infinite on the whole stable interval;
- even in the strongly numerator-like regime, the robust dynamic ceiling stays
  above the universal static robust budget.
"""

from __future__ import annotations
import math
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


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.factor(sp.simplify(sp.expand(z))))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.factor(sp.simplify(sp.expand(expr)))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("STAGE 013 — SELECTED-BRANCH DYNAMIC-WINDOW COMPILER")

# ---------------------------------------------------------------------------
# Exact Stage-012 classifier and its monotonicity
# ---------------------------------------------------------------------------

subbanner("I. Exact Stage-012 classifier and monotonicity")

xi, delta = sp.symbols("xi delta", positive=True, real=True)
R_ND = sp.simplify(
    72 * delta**2 * (1 - xi)
    / ((9 * delta + 11 * xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2))
)

dR_dxi = sp.factor(sp.simplify(sp.diff(R_ND, xi)))

dR_positive_part = sp.factor(
    81 * delta**3
    + 261 * delta**2
    + 297 * delta * xi * (2 - xi)
    + xi**2 * (363 - 242 * xi)
)

print("R_ND(xi,delta) =")
sp.pprint(R_ND)
print("dR_ND/dxi =")
sp.pprint(dR_dxi)
print("positive numerator factor in -dR_ND/dxi =")
sp.pprint(dR_positive_part)

expect_zero(
    "dR_ND/dxi + 72 delta^2 * positive_part / denom",
    sp.simplify(
        dR_dxi
        + 72 * delta**2 * dR_positive_part
        / ((9 * delta + 11 * xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 2)
    ),
)

print("On the stable interval 0 <= xi < 1 and delta > 0, we have:")
print("  2 - xi > 1 > 0,")
print("  363 - 242 xi > 121 > 0,")
print("so the displayed positive_part is strictly positive.")
print("Therefore dR_ND/dxi < 0: the selected-branch classifier decreases monotonically with xi.")

# ---------------------------------------------------------------------------
# Carried Stage-011 rigid dynamic data, converted to per-unit-Xi_1 slopes
# ---------------------------------------------------------------------------

subbanner("II. Carried Stage-011 rigid dynamic data")

# denominator-carried branch = Stage-011 numerator-rigid survivor (pi_1 = 0)
Xi_den = sp.Float("1.73611235", 50)
up_den = sp.Float("-0.52346582", 50)
lo_den = sp.Float("0.71358484", 50)

# numerator-carried branch = Stage-011 denominator-rigid survivor (delta_1 = 0)
Xi_num = sp.Float("0.69293215", 50)
up_num = sp.Float("-0.35245541", 50)
lo_num = sp.Float("-0.23169484", 50)

s_up_den = sp.N(up_den / Xi_den, 30)
s_lo_den = sp.N(lo_den / Xi_den, 30)
s_up_num = sp.N(up_num / Xi_num, 30)
s_lo_num = sp.N(lo_num / Xi_num, 30)

print("Per-unit-Xi_1 dynamic slopes from Stage 011:")
print("  denominator-carried branch (numerator-rigid):")
print(f"    s_+^(den) = {s_up_den}")
print(f"    s_-^(den) = {s_lo_den}")
print("  numerator-carried branch (denominator-rigid):")
print(f"    s_+^(num) = {s_up_num}")
print(f"    s_-^(num) = {s_lo_num}")

# ---------------------------------------------------------------------------
# Exact selected-branch share compiler
# ---------------------------------------------------------------------------

subbanner("III. Exact rigid-split share compiler")

R = sp.symbols("R", nonnegative=True, real=True)
w_num = sp.simplify(R / (1 + R))
w_den = sp.simplify(1 / (1 + R))
expect_zero("w_num + w_den - 1", sp.simplify(w_num + w_den - 1))

print("Selected-branch shares as functions of the Stage-012 classifier R_ND:")
print("w_num =")
sp.pprint(w_num)
print("w_den =")
sp.pprint(w_den)

S_plus = sp.simplify(w_num * s_up_num + w_den * s_up_den)
S_minus = sp.simplify(w_num * s_lo_num + w_den * s_lo_den)

print("Selected upper wall-like dynamic slope per unit Xi_1:")
sp.pprint(S_plus)
print("Selected lower wall-like dynamic slope per unit Xi_1:")
sp.pprint(S_minus)

expect_zero(
    "(1+R) S_plus - (R s_+^(num) + s_+^(den))",
    sp.simplify((1 + R) * S_plus - (R * s_up_num + s_up_den)),
)
expect_zero(
    "(1+R) S_minus - (R s_-^(num) + s_-^(den))",
    sp.simplify((1 + R) * S_minus - (R * s_lo_num + s_lo_den)),
)

# ---------------------------------------------------------------------------
# Exact sign split and dynamic-sign thresholds
# ---------------------------------------------------------------------------

subbanner("IV. Exact sign split")

# upper slope always negative because both carried rigid slopes are negative
print("Upper wall-like pole:")
print(f"  s_+^(num) = {s_up_num} < 0")
print(f"  s_+^(den) = {s_up_den} < 0")
print("Therefore S_+(R) < 0 for every R >= 0.")

R_star = sp.N(s_lo_den / (-s_lo_num), 30)
delta_star_dyn = sp.N(8 / (9 * R_star), 30)

print("Lower wall-like sign threshold R_* from S_-(R_*) = 0:")
print(f"  R_* = {R_star}")
print(f"  delta_*^(dyn) = 8 / (9 R_*) = {delta_star_dyn}")

print("Sign theorem for the lower wall-like pole:")
print("  - if 0 <= R < R_*, then S_-(R) > 0  (lower pole improves)")
print("  - if R = R_*,       then S_-(R) = 0")
print("  - if R > R_*,       then S_-(R) < 0  (both poles worsen)")
print()
print("Since every denominator-like point has R_ND <= 1 and 1 < R_*,")
print("every denominator-like point inherits a split-sign dynamic response:")
print("the upper wall-like pole worsens, but the lower one improves.")
print()
print("Since R_ND(0,delta) = 8/(9 delta) and R_ND decreases monotonically in xi,")
print(f"every branch with delta >= {delta_star_dyn} stays below R_* on the whole stable interval,")
print("so its nonempty dynamic ceiling is infinite on the whole selected branch.")

# ---------------------------------------------------------------------------
# Dynamic ceilings in |eps Xi_1|
# ---------------------------------------------------------------------------

subbanner("V. Dynamic ceilings in |eps Xi_1|")

Rminus0 = sp.Float("30.199907560250075", 50)
Rplus0 = sp.Float("36.171186483269487", 50)
Rreq = sp.Float("21.8545662963584", 50)

ell_minus = sp.N(sp.log(Rminus0 / Rreq), 30)
ell_plus = sp.N(sp.log(Rplus0 / Rreq), 30)

print("Wall-like dynamic margins:")
print(f"  ell_- = ln(R_Q,- / R_req) = {ell_minus}")
print(f"  ell_+ = ln(R_Q,+ / R_req) = {ell_plus}")


def both_ceiling(R_value: float) -> float:
    su = float(sp.N(S_plus.subs({R: sp.Float(str(R_value), 50)}), 30))
    sl = float(sp.N(S_minus.subs({R: sp.Float(str(R_value), 50)}), 30))
    c_plus = math.inf if su >= 0.0 else float(ell_plus) / (-su)
    c_minus = math.inf if sl >= 0.0 else float(ell_minus) / (-sl)
    return min(c_plus, c_minus)


def nonempty_ceiling(R_value: float) -> float:
    su = float(sp.N(S_plus.subs({R: sp.Float(str(R_value), 50)}), 30))
    sl = float(sp.N(S_minus.subs({R: sp.Float(str(R_value), 50)}), 30))
    c_plus = math.inf if su >= 0.0 else float(ell_plus) / (-su)
    c_minus = math.inf if sl >= 0.0 else float(ell_minus) / (-sl)
    return max(c_plus, c_minus)

# exact endpoint / threshold values from monotonicity
both_at_zero = both_ceiling(0.0)
both_at_inf = min(float(ell_plus) / (-float(s_up_num)), float(ell_minus) / (-float(s_lo_num)))
nonempty_at_inf = max(float(ell_plus) / (-float(s_up_num)), float(ell_minus) / (-float(s_lo_num)))

print("Selected-branch robust dynamic ceiling B_dyn^(both)(R) =")
print("  min( ell_+ / (-S_+(R)), ell_- / (-S_-(R)) ), with the second term = +inf when S_-(R) >= 0.")
print()
print("Selected-branch nonempty dynamic ceiling B_dyn^(nonempty)(R) =")
print("  +inf                                 if S_-(R) >= 0,")
print("  max( ell_+ / (-S_+(R)), ell_- / (-S_-(R)) )   if S_-(R) < 0.")
print()
print(f"Robust dynamic ceiling endpoints:")
print(f"  B_dyn^(both)(0)   = {both_at_zero:.15f}")
print(f"  lim_(R->inf) B_dyn^(both)(R) = {both_at_inf:.15f}")
print(f"Finite nonempty dynamic ceiling infimum:")
print(f"  lim_(R->inf) B_dyn^(nonempty)(R) = {nonempty_at_inf:.15f}")
print(f"Nonempty dynamic ceiling is +inf for all R <= R_* ≈ {float(R_star):.15f}.")

# ---------------------------------------------------------------------------
# Universal transported static ceilings in |eps Xi_1|
# ---------------------------------------------------------------------------

subbanner("VI. Universal transported static ceilings in |eps Xi_1|")

# From Stage 011 branch-specific |eps t| ceilings, converted using Xi_1 = sigma * t.
stat_both_den = sp.N(Xi_den * sp.Float("0.21192772", 50), 30)
stat_nonempty_den = sp.N(Xi_den * sp.Float("0.42486828", 50), 30)
stat_both_num = sp.N(Xi_num * sp.Float("0.53097598", 50), 30)
stat_nonempty_num = sp.N(Xi_num * sp.Float("1.06448959", 50), 30)

print("Converted static ceilings from the two rigid branches:")
print(f"  denominator-carried robust    = {stat_both_den}")
print(f"  numerator-carried robust      = {stat_both_num}")
print(f"  denominator-carried nonempty  = {stat_nonempty_den}")
print(f"  numerator-carried nonempty    = {stat_nonempty_num}")

print("These agree up to the carried Stage-011 numerical rounding, so the selected branch inherits the universal static budgets")
print("  B_stat^(both)    ≈ 0.367930328492646")
print("  B_stat^(nonempty)≈ 0.737619063660757")

B_stat_both = sp.Float("0.367930328492646", 50)
B_stat_nonempty = sp.Float("0.737619063660757", 50)

print("Comparison against the dynamic ceilings:")
print(f"  inf_R B_dyn^(both)(R)      = {both_at_inf:.15f}  >  {float(B_stat_both):.15f}")
print(f"  inf_(finite R) B_dyn^(nonempty)(R) = {nonempty_at_inf:.15f}  >  {float(B_stat_nonempty):.15f}")

# ---------------------------------------------------------------------------
# Sample classifier points
# ---------------------------------------------------------------------------

subbanner("VII. Sample classifier points")

for R_val in [0.0, 1.0, float(R_star), 10.0]:
    su = float(sp.N(S_plus.subs({R: sp.Float(str(R_val), 50)}), 30))
    sl = float(sp.N(S_minus.subs({R: sp.Float(str(R_val), 50)}), 30))
    both = both_ceiling(R_val)
    nonempty = nonempty_ceiling(R_val)
    print(f"R = {R_val:.12f}")
    print(f"  S_+ = {su:.15f}")
    print(f"  S_- = {sl:.15f}")
    print(f"  B_dyn^(both) = {both:.15f}")
    if math.isfinite(nonempty):
        print(f"  B_dyn^(nonempty) = {nonempty:.15f}")
    else:
        print("  B_dyn^(nonempty) = +inf")
    print()

# ---------------------------------------------------------------------------
# Final verdict
# ---------------------------------------------------------------------------

subbanner("VIII. Final verdict")
print("The Stage-012 selected-branch classifier can be wired exactly into the carried Stage-011 rigid dynamic data.")
print("Within that exact rigid-split compiler on the concrete compatibility branch:")
print("  * the upper wall-like pole always worsens,")
print("  * the lower wall-like pole improves whenever R_ND <= R_* ≈ 1.2292554384633356,")
print("  * every denominator-like point (R_ND <= 1) therefore has infinite nonempty dynamic ceiling,")
print("  * if delta >= 8/9, the whole selected branch is denominator-like, so the nonempty dynamic ceiling is infinite on the whole stable interval,")
print("  * if delta >= 8/(9 R_*) ≈ 0.723111617875019, the whole selected branch still stays below the sign-flip threshold R_* and again has infinite nonempty dynamic ceiling,")
print("  * and even in the strongly numerator-like regime, the robust dynamic ceiling remains above the universal transported static ceiling.")
print()
print("So Stage 013 does not kill the same-charge corridor.")
print("It sharpens the verdict instead:")
print("the first kill condition on the selected branch is still the transported static Xi_1 budget, not the wall-like dynamic window.")
