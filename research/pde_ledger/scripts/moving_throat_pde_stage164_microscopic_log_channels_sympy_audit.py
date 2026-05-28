#!/usr/bin/env python3
"""
moving_throat_pde_stage164_microscopic_log_channels_sympy_audit.py

SymPy-backed audit for Stage 164.

Checks:
1. Exact identity of the Stage 163 logarithmic channels with d ln(g/r) and -d ln r_c.
2. Exact general Stage 118 channel products.
3. Uniform-overlap + D/N simplification.
4. Healing-locked branch product formulas.
5. Explicit delta_perp formula and tangency law on the healing-locked branch.
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

banner("STAGE 164 — MICROSCOPIC LOG-IMBALANCE CHANNELS")

# ------------------------------------------------------------------
# 1. Exact parent-ratio identities
# ------------------------------------------------------------------
Ks, Kq, lam, gs, gq = sp.symbols("K_s K_q lam g_s g_q", positive=True, real=True)

r = lam / sp.sqrt(Ks * Kq)
g = gq * sp.sqrt(Ks) / (gs * sp.sqrt(Kq))
rc = lam**2 / (Ks * Kq)

expect_zero("g/r - (g_q K_s)/(g_s lam)", g/r - gq*Ks/(gs*lam))
expect_zero("1/r_c - (K_s K_q)/lam^2", 1/rc - Ks*Kq/lam**2)
print("Therefore:")
print("  dln[(g_q K_s)/(g_s lam)] = dln(g/r)")
print("  dln[(K_s K_q)/lam^2] = - dln(r_c)")

# ------------------------------------------------------------------
# 2. General Stage 118 formulas
# ------------------------------------------------------------------
banner("General Stage 118 formulas")

Tm, Js, Zq, mu0, Lw, qstar, vw0, Isq, cs = sp.symbols(
    "T_m J_s Z_q mu_0 L_W qstar v_w0 I_sq c_s", positive=True, real=True
)

gs_gen = Tm * Js
gq_gen = Zq * sp.pi / (sp.sqrt(2) * mu0 * Lw**sp.Rational(3, 2))
Kq_gen = Zq * sp.pi**2 * cs**2 / (4 * mu0 * Lw**2)
lam_gen = -qstar * vw0 * Isq

first_prod_general = sp.simplify(gq_gen * Ks / (gs_gen * lam_gen))
second_prod_general = sp.simplify(Ks * Kq_gen / lam_gen**2)

print("general first product =", first_prod_general)
print("general second product =", second_prod_general)

# ------------------------------------------------------------------
# 3. Uniform-overlap + D/N simplification
# ------------------------------------------------------------------
banner("Uniform-overlap + D/N simplification")

a, ell = sp.symbols("a ell", positive=True, real=True)
Js_closure = 4 * sp.pi * a**2 * ell / 3
Iq = 2 * sp.sqrt(2 * Lw) / sp.pi
Isq_closure = Js_closure * Iq

first_prod_uniform = sp.simplify(first_prod_general.subs({Js: Js_closure, Isq: Isq_closure}))
second_prod_uniform = sp.simplify(second_prod_general.subs({Isq: Isq_closure}))

print("uniform-overlap first product =", first_prod_uniform)
print("uniform-overlap second product =", second_prod_uniform)

# ------------------------------------------------------------------
# 4. Healing-locked shell simplification
# ------------------------------------------------------------------
banner("Healing-locked shell branch")

mpsi, hbar, rho_w, csw = sp.symbols("mpsi hbar rho_w c_sw", positive=True, real=True)
ell_lock = hbar / (2 * mpsi * csw)
Ks_lock = sp.simplify(3 * sp.pi * a**2 * hbar**2 / (5 * mpsi * rho_w * ell_lock))

first_prod_heal = sp.simplify(first_prod_uniform.subs({Ks: Ks_lock, ell: ell_lock}))
second_prod_heal = sp.simplify(second_prod_uniform.subs({Ks: Ks_lock, ell: ell_lock}))

print("healing-locked first product =", first_prod_heal)
print("healing-locked second product =", second_prod_heal)

first_expected = -sp.Rational(27,40) * sp.pi * mpsi**2 * Zq * csw**3 / (hbar * mu0 * qstar * rho_w * Tm * vw0 * a**2 * Lw**2)
second_expected = sp.Rational(27,320) * sp.pi**3 * mpsi**2 * Zq * cs**2 * csw**3 / (hbar * mu0 * qstar**2 * rho_w * vw0**2 * a**2 * Lw**3)

expect_zero("healing first product exact formula", first_prod_heal - first_expected)
expect_zero("healing second product exact formula", second_prod_heal - second_expected)

# ------------------------------------------------------------------
# 5. Explicit delta_perp and tangency law
# ------------------------------------------------------------------
banner("delta_perp on the healing-locked branch")

dlnZ, dlnrho, dlncsw, dlncs, dlnTm, dlnv, dlna, dlnLw = sp.symbols(
    "dlnZ dlnrho dlncsw dlncs dlnTm dlnv dlna dlnLw", real=True
)
rstar = sp.symbols("rstar", positive=True, real=True)
gstar = sp.symbols("gstar", positive=True, real=True)
s = sp.sqrt(1 + rstar**2)
b = sp.Rational(1,4) / s

first_heal = dlnZ + 3*dlncsw - dlnrho - dlnTm - dlnv - 2*dlna - 2*dlnLw
second_heal = dlnZ + 2*dlncs + 3*dlncsw - dlnrho - 2*dlnv - 2*dlna - 3*dlnLw

delta_perp = sp.expand(gstar*first_heal + b*second_heal)
delta_perp_expected = sp.expand(
    (gstar + b)*(dlnZ - dlnrho)
    + 3*(gstar + b)*dlncsw
    + 2*b*dlncs
    - gstar*dlnTm
    - (gstar + 2*b)*dlnv
    - 2*(gstar + b)*dlna
    - (2*gstar + 3*b)*dlnLw
)
expect_zero("delta_perp compressed form", delta_perp - delta_perp_expected)
print("delta_perp =", sp.simplify(delta_perp_expected))

dlnTm_sol = sp.simplify(sp.solve(sp.Eq(delta_perp_expected, 0), dlnTm)[0])
print("tangency-law dln T_m =", dlnTm_sol)

# ------------------------------------------------------------------
# 6. Family-1 numerical coefficients
# ------------------------------------------------------------------
banner("Family-1 numerical coefficients")

g_num = sp.Float("0.758035078944663", 30)
r_num = sp.Float("1.77799353547498", 30)
s_num = sp.sqrt(1 + r_num**2)
b_num = sp.N(sp.Rational(1,4) / s_num, 20)
A_num = sp.N(g_num + b_num, 20)
B_num = sp.N(2*b_num, 20)

print("g_* =", g_num)
print("1/(4 sqrt(1+r_*^2)) =", b_num)
print("A_* = g_* + 1/(4 sqrt(...)) =", A_num)
print("B_* = 1/(2 sqrt(1+r_*^2)) =", B_num)
print("coeff[ln(Z_q/rho_w)] =", A_num)
print("coeff[ln(c_sw)]      =", sp.N(3*A_num, 20))
print("coeff[ln(c_s)]       =", B_num)
print("coeff[ln(T_m)]       =", sp.N(-g_num, 20))
print("coeff[ln(v_w0)]      =", sp.N(-(g_num + B_num), 20))
print("coeff[ln(a)]         =", sp.N(-2*A_num, 20))
print("coeff[ln(L_W)]       =", sp.N(-(2*g_num + 3*b_num), 20))

print("\nCarry-forward formulas:")
print("  dln[(g_q K_s)/(g_s lam)] = dln(g/r)")
print("  dln[(K_s K_q)/lam^2] = - dln(r_c)")
print("  On the healing-locked branch:")
print("    dln[(g_q K_s)/(g_s lam)] = dln Z_q + 3 dln c_sw - dln rho_w - dln T_m - dln v_w0 - 2 dln a - 2 dln L_W")
print("    dln[(K_s K_q)/lam^2] = dln Z_q + 2 dln c_s + 3 dln c_sw - dln rho_w - 2 dln v_w0 - 2 dln a - 3 dln L_W")
print("  delta_perp is the weighted sum of those two channels.")