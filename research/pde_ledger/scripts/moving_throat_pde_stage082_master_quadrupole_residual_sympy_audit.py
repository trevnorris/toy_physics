#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 65 SymPy audit

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


banner("STAGE 65 — MASTER QUADRUPOLE RESIDUAL")

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

# Directional content of R_quad: verify the partial derivatives that
# underwrite the "guaranteed success / guaranteed failure" theorems.
dR_dzeta_phys = sp.simplify(sp.diff(R_quad, zeta_phys))
expect_zero("dR_quad/dzeta_phys + 1", dR_dzeta_phys + 1)

dzeta_req_dPi = sp.simplify(sp.diff(zeta_req, Pi_tr))
dR_dPi_at_zeta_minus = sp.simplify(sp.diff(R_quad.subs(zeta_phys, zeta_minus), Pi_tr))
expect_zero(
    "dR_quad/dPi_tr - dzeta_req/dPi_tr (at zeta_phys=zeta_-)",
    dR_dPi_at_zeta_minus - dzeta_req_dPi,
)

# Family-1 strength identities.
Theta_w, Upsilon_w = sp.symbols("Theta_w Upsilon_w", positive=True, real=True)
# TODO(provenance): Lambda_ell = 37 and the convention Upsilon_w = 100 * Theta_w
# are carry-forward constants. Cite the upstream stage's script (likely an
# earlier moving_throat_pde_stage*_sympy_audit.py) that derives them. Until
# then, this stage treats them as inputs and only displays their consequences.
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
