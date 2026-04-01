#!/usr/bin/env python3
"""
SymPy audit for Stage 64.

Translates the explicit Family-1 zeta-demand thresholds into the selected-branch
product variable Pi_tr / C_mix.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS")

zeta, Pi, Cmix, eps_blk = sp.symbols('zeta Pi Cmix eps_blk', real=True)

# Exact Stage-35 support-demand law and its inversion.
zeta_expr = (Pi - Cmix) / (Cmix - eps_blk * (2 * Cmix - Pi))
Pi_of_zeta = sp.solve(sp.Eq(zeta, zeta_expr), Pi)[0]
Pi_of_zeta = sp.simplify(Pi_of_zeta)
Q = sp.simplify(Pi_of_zeta / Cmix)

print("Pi_of_zeta =", Pi_of_zeta)
print("Q(zeta;eps_blk) =", Q)
expect_zero("Q(0)-1", sp.simplify(Q.subs(zeta, 0) - 1))
expect_zero("Q(1)-2", sp.simplify(Q.subs(zeta, 1) - 2))
print("dQ/dzeta =", sp.simplify(sp.diff(Q, zeta)))

# Numerical zeta thresholds from Stage 63 at lambda_mu = 1.
zeta_suff_chi_1 = sp.Float('2.46622291347846')
zeta_fail_chi_1 = sp.Float('2.46752913273870')
zeta_suff_J_1 = sp.Float('2.44257571477179')
zeta_fail_J_1 = sp.Float('2.46752736855058')
zeta_max_F1 = sp.Float('2.46752922945601')

Pi_suff_chi_over_C = sp.simplify(Q.subs({zeta: zeta_suff_chi_1}))
Pi_fail_chi_over_C = sp.simplify(Q.subs({zeta: zeta_fail_chi_1}))
Pi_suff_J_over_C = sp.simplify(Q.subs({zeta: zeta_suff_J_1}))
Pi_fail_J_over_C = sp.simplify(Q.subs({zeta: zeta_fail_J_1}))
Pi_max_over_C = sp.simplify(Q.subs({zeta: zeta_max_F1}))

print("Pi_suff^(chi)/C_mix =", Pi_suff_chi_over_C)
print("Pi_fail^(chi)/C_mix =", Pi_fail_chi_over_C)
print("Pi_suff^(J)/C_mix   =", Pi_suff_J_over_C)
print("Pi_fail^(J)/C_mix   =", Pi_fail_J_over_C)
print("Pi_max^(F1)/C_mix   =", Pi_max_over_C)

# Unblocked numerical illustration.
print("At eps_blk = 0:")
print("  Pi_suff^(chi)/C_mix =", sp.N(Pi_suff_chi_over_C.subs(eps_blk, 0), 20))
print("  Pi_fail^(chi)/C_mix =", sp.N(Pi_fail_chi_over_C.subs(eps_blk, 0), 20))
print("  Pi_suff^(J)/C_mix   =", sp.N(Pi_suff_J_over_C.subs(eps_blk, 0), 20))
print("  Pi_fail^(J)/C_mix   =", sp.N(Pi_fail_J_over_C.subs(eps_blk, 0), 20))
print("  Pi_max^(F1)/C_mix   =", sp.N(Pi_max_over_C.subs(eps_blk, 0), 20))

# Hard blocking ceiling for denominator positivity at the Family-1 ceiling.
eps_ceiling = sp.simplify(1 / zeta_max_F1)
print("Blocking ceiling eps_blk <", sp.N(eps_ceiling, 20))

banner("FINAL LEDGER")
print("Pi_tr = C_mix Q(zeta_req;eps_blk)")
print("Q(zeta;eps_blk) = [1 + (1-2 eps_blk) zeta] / [1 - eps_blk zeta]")
print("The explicit Family-1 hard ceiling is Pi_max^(F1) = C_mix Q(zeta_max^(F1);eps_blk).")
print("At eps_blk = 0 and lambda_mu = 1 the natural success threshold is Pi_tr/C_mix <= 3.46622291347846.")
