#!/usr/bin/env python3
"""
moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py

SymPy-backed audit for the exact lower-compensated Family-1 branch drift laws.

Checks:
1. D/N law implies d ln L_W = d ln a at fixed r_*.
2. Fixed-r and fixed-g conditions solve for d ln v_{w0} and d ln T_m.
3. Product/ratio transport laws collapse exactly.
4. n=5 wall EOS gives d ln c_{s,w} = 2 d ln rho_w.
5. The Stage 164 off-family channels vanish identically after substitution.
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

banner("STAGE 165 — EXACT LOWER-BRANCH DRIFT LAWS")

dZ, drho, dcsw, dcs, dT, dv, da, dLW = sp.symbols(
    "dZ drho dcsw dcs dT dv da dLW", real=True
)

r = sp.symbols("r", positive=True, real=True)
a, LW = sp.symbols("a LW", positive=True, real=True)

# 1. D/N law
LW_law = sp.pi * a * sp.sqrt(1 + r**2) / (2 * sp.sqrt(3))
dlog_LW = sp.simplify(a * sp.diff(sp.log(LW_law), a))
print("D/N law: L_W =", LW_law)
print("At fixed r_* : d ln L_W = d ln a")
expect_zero("d ln L_W - d ln a at fixed r_*", dlog_LW - 1)

# 2. Fixed-r and fixed-g channel conditions from Stage 164
eq_r = sp.Eq(dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW, 0)
eq_g = sp.Eq(dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW, 0)
sol = sp.solve([eq_r, eq_g], [dv, dT], dict=True)[0]
dv_sol = sp.simplify(sol[dv])
dT_sol = sp.simplify(sol[dT])

print("\nExact drift laws before using d ln L_W = d ln a:")
print("d ln v_w0 =", dv_sol)
print("d ln T_m  =", dT_sol)

dv_DN = sp.simplify(dv_sol.subs(dLW, da))
dT_DN = sp.simplify(dT_sol.subs(dLW, da))

print("\nAfter using d ln L_W = d ln a:")
print("d ln v_w0 =", dv_DN)
print("d ln T_m  =", dT_DN)

# 3. Product/ratio laws
ratio_law = sp.simplify(dv_DN - dT_DN)
prod_law = sp.simplify(dv_DN + dT_DN)
expect_zero("ratio law - (2 d ln c_s - d ln a)", ratio_law - (2*dcs - da))
expect_zero("product law - (dZ + 3 dcsw - drho - 4 da)", prod_law - (dZ + 3*dcsw - drho - 4*da))

# 4. n=5 wall EOS
dv_n5 = sp.simplify(dv_DN.subs(dcsw, 2*drho))
dT_n5 = sp.simplify(dT_DN.subs(dcsw, 2*drho))
print("\nAfter n=5 wall EOS (d ln c_s,w = 2 d ln rho_w):")
print("d ln v_w0 =", dv_n5)
print("d ln T_m  =", dT_n5)
expect_zero(
    "n=5 product law - (dZ + 5 drho - 4 da)",
    sp.simplify((dv_n5 + dT_n5) - (dZ + 5*drho - 4*da)),
)

# 5. Stage 164 off-family channels vanish identically
# NOTE: channel_g/channel_r are the LHS of eq_g/eq_r evaluated at the solved (dv,dT).
# They vanish by construction (substituting a linear solution into its own system),
# so this is a consistency print of the solver, not an independent verification of
# delta_perp = 0. Reported, not asserted.
channel_g = sp.simplify((dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW).subs({dv: dv_sol, dT: dT_sol}))
channel_r = sp.simplify((dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW).subs({dv: dv_sol, dT: dT_sol}))
print("fixed-g channel (solver consistency) =", channel_g)
print("fixed-r channel (solver consistency) =", channel_r)

print("\nCarry-forward formulas:")
print("  d ln L_W = d ln a")
print("  d ln v_w0 = 1/2 d ln(Z_q/rho_w) + 3/2 d ln c_s,w + d ln c_s - d ln a - 3/2 d ln L_W")
print("  d ln T_m  = 1/2 d ln(Z_q/rho_w) + 3/2 d ln c_s,w - d ln c_s - d ln a - 1/2 d ln L_W")
print("  d ln(v_w0/T_m) = 2 d ln c_s - d ln a")
print("  d ln(v_w0 T_m) = d ln Z_q + 3 d ln c_s,w - d ln rho_w - 4 d ln a")
print("  n=5 wall EOS: d ln c_s,w = 2 d ln rho_w")

banner("LOWER-BRANCH NUMERIC PREFACTORS")

def expect_close(name: str, value: sp.Expr, target: sp.Expr, tol: float = 1e-12) -> None:
    value_num = sp.N(value, 30)
    diff = abs(float(sp.N(value - target, 30)))
    print(f"{name} = {value_num}")
    print(f"{name} diff = {diff}")
    if not diff < tol:
        raise AssertionError(f"{name} mismatch: diff={diff}")

rstar = sp.sqrt(4107 - 100*sp.pi**2) / (10*sp.pi)
gstar = sp.Float("0.758035078944663", 30)

Tm_pref = 3*sp.sqrt(10)*sp.Integer(3)**sp.Rational(3, 4) / (5*sp.pi*gstar*(1 + rstar**2)**sp.Rational(1, 4))
v_pref = 9*sp.sqrt(10)*sp.Integer(3)**sp.Rational(1, 4)*rstar / (20*(1 + rstar**2)**sp.Rational(3, 4))
ratio_pref = sp.sqrt(3)*sp.pi*gstar*rstar / (4*sp.sqrt(1 + rstar**2))
prod_pref = 81*rstar / (10*sp.pi*gstar*(1 + rstar**2))

expect_close("Tm_pref", Tm_pref, sp.Float("1.2715890393387603", 30))
expect_close("v_pref", v_pref, sp.Float("1.1428896163056477", 30))
expect_close("ratio_pref", ratio_pref, sp.Float("0.8987885086678338", 30))
expect_close("prod_pref", prod_pref, sp.Float("1.4532859092683434", 30))
