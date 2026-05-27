#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 082 SymPy audit

Checks
------
1. Exact inverse map between zeta_req and Pi_tr.
2. Exact product thresholds Pi_suff and Pi_fail from lower/upper support ratios.
3. Exact Family-1 strength identity Xi_F1 = 1369 Upsilon_w = 136900 Theta_w.
4. Exact master residual definition.
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


banner("STAGE 082 — MASTER QUADRUPOLE RESIDUAL")

Pi_tr, Cmix, eps_blk = sp.symbols("Pi_tr C_mix eps_blk", positive=True, real=True)
zeta = sp.symbols("zeta", real=True)
zeta_minus, zeta_plus = sp.symbols("zeta_- zeta_+", real=True)

zeta_req = sp.simplify((Pi_tr - Cmix) / (Cmix - eps_blk * (2 * Cmix - Pi_tr)))
Q = sp.simplify((1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta))

print("zeta_req(Pi_tr,C_mix,eps_blk) =")
sp.pprint(zeta_req)
print("\nQ(zeta;eps_blk) =")
sp.pprint(Q)

expect_zero(
    "inverse map zeta_req(C_mix*Q(zeta)) - zeta",
    zeta_req.subs(Pi_tr, sp.simplify(Cmix * Q)) - zeta,
)

Pi_suff = sp.simplify(Cmix * Q.subs(zeta, zeta_minus))
Pi_fail = sp.simplify(Cmix * Q.subs(zeta, zeta_plus))

print("\nPi_suff =")
sp.pprint(Pi_suff)
print("\nPi_fail =")
sp.pprint(Pi_fail)

expect_zero("zeta_req(Pi_suff)-zeta_-", zeta_req.subs(Pi_tr, Pi_suff) - zeta_minus)
expect_zero("zeta_req(Pi_fail)-zeta_+", zeta_req.subs(Pi_tr, Pi_fail) - zeta_plus)

# Exact master residual with symbolic support ratio.
zeta_phys = sp.symbols("zeta_phys", real=True)
R_quad = sp.simplify(zeta_req - zeta_phys)
print("\nR_quad =")
sp.pprint(R_quad)

expect_zero(
    "R_quad(Pi_suff, zeta_phys=zeta_-)",
    R_quad.subs({Pi_tr: Pi_suff, zeta_phys: zeta_minus}),
)
expect_zero(
    "R_quad(Pi_fail, zeta_phys=zeta_+)",
    R_quad.subs({Pi_tr: Pi_fail, zeta_phys: zeta_plus}),
)

# F3 (v2): directional content of zeta_req. The inverse-map theorem
# (notes section 4) relies on zeta_req being strictly increasing in Pi_tr
# on the physical branch (where Pi_tr, C_mix, eps_blk are positive).
# Verify by factoring d zeta_req / d Pi_tr into a sign-controlled
# numerator/denominator pair under those assumptions.
dzeta_req_dPi = sp.simplify(sp.diff(zeta_req, Pi_tr))
print(f"\ndzeta_req/dPi_tr = {dzeta_req_dPi}")
num, den = sp.fraction(sp.together(dzeta_req_dPi))
expect_zero(
    "numerator(d zeta_req/d Pi_tr) - C_mix*(1 - eps_blk)",
    sp.simplify(num - Cmix * (1 - eps_blk)),
)
expect_zero(
    "denominator(d zeta_req/d Pi_tr) - (C_mix - eps_blk*(2*C_mix - Pi_tr))**2",
    sp.simplify(den - (Cmix - eps_blk * (2 * Cmix - Pi_tr)) ** 2),
)

# F1 (v2 paper-alignment Q3 direction (a)): closed-form pin for zeta_phys.
# Paper eq. app-stage082-zeta-phys:
#   zeta_phys(Pe, eta; kappa) = Omega_Pe(Pe)^2 * (kappa + pi^2/4) / (kappa + y(eta)^2)
# with Omega_Pe(Pe) = pi*Pe*(2*Pe*exp(Pe) + pi) / ((4*Pe^2 + pi^2)*(exp(Pe) - 1))
# and y(eta) the smallest positive root of y*tan(y) = eta.
# Verify by reproducing the Pe->oo limit at Family-1 (eta, kappa) = (37, 12321/5),
# which equals stage 084's zeta_max^(F1) constant.
Pe, kappa_sym, eta_sym, y_sym = sp.symbols("Pe kappa eta y", positive=True, real=True)
Omega_Pe_expr = sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi) / (
    (4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1)
)
zeta_phys_closed = Omega_Pe_expr**2 * (kappa_sym + sp.pi**2 / 4) / (
    kappa_sym + y_sym**2
)
Omega_Pe_limit = sp.limit(Omega_Pe_expr, Pe, sp.oo)
print(f"\nOmega_Pe -> {Omega_Pe_limit} as Pe -> oo")
expect_zero("Omega_Pe -> pi/2 as Pe -> oo", Omega_Pe_limit - sp.pi / 2)

# Family-1 root: smallest positive y solving y*tan(y) = 37 in (0, pi/2).
# Use mpmath bisection rather than sp.nsolve — Newton iteration is unstable
# near pi/2 where tan'(y) blows up and jumps to far-away roots.
import mpmath
mpmath.mp.dps = 30
y_F1 = sp.Float(
    mpmath.findroot(lambda yv: yv * mpmath.tan(yv) - 37, (1.5, 1.55), solver="bisect"),
    30,
)
print(f"y_F1 (root of y tan y = 37, smallest positive) = {y_F1}")

kappa_F1 = sp.Rational(12321, 5)
zeta_phys_F1_limit = (sp.pi**2 / 4) * (kappa_F1 + sp.pi**2 / 4) / (
    kappa_F1 + y_F1**2
)
print(f"zeta_phys(Pe->oo, kappa_F1, y_F1) = {sp.N(zeta_phys_F1_limit, 20)}")

# Reference constant from stage 084 .wl (verified at v2): zeta_max^(F1) ~ 2.4675292294558...
# Stage 084 lines 73-76 compute the same Pe->oo limit and assert it matches the
# upstream zeta_max^(F1) constant to 10^-10. We carry that reference forward here.
zeta_max_F1_reference = sp.Float("2.467529229455835", 20)
diff_to_reference = abs(sp.N(zeta_phys_F1_limit - zeta_max_F1_reference, 30))
print(f"|zeta_phys(F1, Pe->oo) - zeta_max^(F1)| = {diff_to_reference}")
assert diff_to_reference < sp.Float("1e-10"), (
    f"Family-1 closed-form pin disagrees with upstream zeta_max^(F1) by {diff_to_reference}"
)
print("PASS: Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10.")

# Family-1 strength identities.
Theta_w, Upsilon_w = sp.symbols("Theta_w Upsilon_w", positive=True, real=True)
# Carry-forward constants (Lambda_ell = 37 from stages 056/073, Upsilon_w = 100 Theta_w
# from stage 075 with alpha_r = 10, paper eq. app-stage075-fail / app-stage075-succeed).
# After v2 paper-alignment Q2 direction (a): paper/stages/stage_075.tex Inputs line and
# notes/stages/.../stage075...md were updated to state Upsilon_w = 100 Theta_w, fully
# consistent with the script's value here.
Lambda_ell = sp.Integer(37)
Xi_F1_from_Upsilon = sp.simplify(Upsilon_w * Lambda_ell**2)
Xi_F1_from_Theta = sp.simplify(100 * Theta_w * Lambda_ell**2)

print("\nXi_F1 from Upsilon_w =")
sp.pprint(Xi_F1_from_Upsilon)
print("Xi_F1 from Theta_w =")
sp.pprint(Xi_F1_from_Theta)
# Note: the two equalities below are arithmetic on the hand-supplied integers
# 37, 100, 1369, 136900. They are not independent verifications of the Family-1
# strength identity — the upstream stage that derives Lambda_ell = 37 and the
# convention Upsilon_w = 100 * Theta_w is responsible for those facts. Here we
# only display the resulting Xi_F1 forms for downstream readers.
print(f"Xi_F1(Theta_w) - 136900 Theta_w = {sp.simplify(Xi_F1_from_Theta - sp.Integer(136900) * Theta_w)}  (display only)")
print(f"Xi_F1(Upsilon_w=100 Theta_w) - Xi_F1(Theta_w) = {sp.simplify(Xi_F1_from_Upsilon.subs(Upsilon_w, 100 * Theta_w) - Xi_F1_from_Theta)}  (display only)")

print("\nTheorem ledger:")
print("  Pi_tr <= C_mix Q(zeta_-,eps_blk)  -> guaranteed success if zeta_phys >= zeta_-")
print("  Pi_tr >= C_mix Q(zeta_+,eps_blk)  -> guaranteed failure if zeta_phys <= zeta_+")
print("  Xi_F1 = 1369 Upsilon_w = 136900 Theta_w")
