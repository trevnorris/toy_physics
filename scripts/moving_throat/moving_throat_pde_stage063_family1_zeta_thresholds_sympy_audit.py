#!/usr/bin/env python3
"""
SymPy audit for Stage 63.

Converts the Stage-61 Family-1 Pe_req thresholds into explicit quadrupole-demand
thresholds zeta_req using the Stage-62 Family-1 demand map.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


banner("STAGE 63 — FAMILY-1 ZETA THRESHOLDS")

lam = sp.symbols('lambda_mu', positive=True, real=True)
Pe = sp.symbols('Pe', positive=True, real=True)
y = sp.symbols('y', real=True)

# Family-1 constants from Stage 62.
y_F1 = sp.nsolve(y * sp.tan(y) - 37, 1.53, tol=1e-34, maxsteps=100)
kappa_F1 = sp.Rational(12321, 5)
A_F1 = (kappa_F1 + sp.pi**2 / 4) / (kappa_F1 + y_F1**2)
Omega = sp.simplify(sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi) / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1)))
zeta = sp.simplify(A_F1 * Omega**2)
zeta_max = sp.simplify(sp.limit(zeta, Pe, sp.oo))

# Stage-61 explicit Pe thresholds.
Pe_suff_chi = sp.Float('96.5285247264386') * lam**2
Pe_fail_chi = sp.Float('11220.5441626259') * lam**2
Pe_suff_J = sp.Float('22.0062226330754') * lam**2
Pe_fail_J = sp.Float('2558.01892349205') * lam**2

zeta_suff_chi = sp.simplify(zeta.subs(Pe, Pe_suff_chi))
zeta_fail_chi = sp.simplify(zeta.subs(Pe, Pe_fail_chi))
zeta_suff_J = sp.simplify(zeta.subs(Pe, Pe_suff_J))
zeta_fail_J = sp.simplify(zeta.subs(Pe, Pe_fail_J))

print("zeta_max^(F1) =", sp.N(zeta_max, 20))
print("zeta_suff^(chi)(lambda_mu) = zeta_F1(96.5285247264386 lambda_mu^2)")
print("zeta_fail^(chi)(lambda_mu) = zeta_F1(11220.5441626259 lambda_mu^2)")
print("zeta_suff^(J)(lambda_mu)   = zeta_F1(22.0062226330754 lambda_mu^2)")
print("zeta_fail^(J)(lambda_mu)   = zeta_F1(2558.01892349205 lambda_mu^2)")

# Numerical values at lambda_mu = 1.
vals = {
    "zeta_suff^(chi)(1)": sp.N(zeta_suff_chi.subs(lam, 1), 20),
    "zeta_fail^(chi)(1)": sp.N(zeta_fail_chi.subs(lam, 1), 20),
    "zeta_suff^(J)(1)": sp.N(zeta_suff_J.subs(lam, 1), 20),
    "zeta_fail^(J)(1)": sp.N(zeta_fail_J.subs(lam, 1), 20),
}
for k, v in vals.items():
    print(k, '=', v)

# Large-lambda limits saturate at the hard ceiling.
for name, expr in [
    ("limit zeta_suff^(chi)", zeta_suff_chi),
    ("limit zeta_fail^(chi)", zeta_fail_chi),
    ("limit zeta_suff^(J)", zeta_suff_J),
    ("limit zeta_fail^(J)", zeta_fail_J),
]:
    lim = sp.simplify(sp.limit(expr, lam, sp.oo))
    print(name, '=', sp.N(lim, 20))
    if abs(complex(sp.N(lim - zeta_max, 30))) > 1e-10:
        raise AssertionError(f"{name} does not saturate to zeta_max^(F1)")

banner("FINAL LEDGER")
print("Guaranteed success: zeta_req <= zeta_suff^(chi)(lambda_mu)")
print("Guaranteed failure: zeta_req >= zeta_fail^(chi)(lambda_mu)")
print("Conservative floor success/failure obtained by replacing (chi) with (J)")
print("For lambda_mu = 1 the natural Family-1 branch already reaches zeta_req ≈ 2.46622 before the indeterminate band.")
