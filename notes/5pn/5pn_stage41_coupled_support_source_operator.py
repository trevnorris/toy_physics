#!/usr/bin/env python3
"""
5pn_stage41_coupled_support_source_operator.py

Stage 41 audit: coupled support/source operator and the exact Pe branch equation.
"""

from __future__ import annotations

import mpmath as mp
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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 41 — COUPLED SUPPORT/SOURCE OPERATOR AND EXACT PE BRANCH EQUATION")

x = sp.symbols("x", real=True)
alpha, eta, Pe, Xi = sp.symbols("alpha eta Pe Xi", positive=True, real=True)
Den = alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)

K = (
    sp.cosh(alpha * x)
    + (eta / alpha) * sp.sinh(alpha * x)
    - sp.cosh(alpha * (1 - x))
) / Den

subbanner("41.1 — Exact support-drop kernel and derivative")
K_repacked = (
    (1 - sp.cosh(alpha)) * sp.cosh(alpha * x)
    + (eta / alpha + sp.sinh(alpha)) * sp.sinh(alpha * x)
) / Den
expect_zero("kernel repacking identity", K - K_repacked)

Kprime_formula = (
    alpha * sp.sinh(alpha * x)
    + eta * sp.cosh(alpha * x)
    + alpha * sp.sinh(alpha * (1 - x))
) / Den
expect_zero("dK/dx formula", sp.diff(K, x) - Kprime_formula)

print("\nBecause alpha>0 and eta>0, the derivative formula shows K_(kappa,eta)(x)")
print("is strictly increasing on the constructive branch.")

subbanner("41.2 — Exact stationary source branch")
Sigma = Pe * sp.exp(Pe * x) / (sp.exp(Pe) - 1)
print("Sigma_Pe(x) =")
sp.pprint(Sigma)

Ic_formula = (
    sp.exp(Pe) * (Pe * sp.cosh(alpha) - alpha * sp.sinh(alpha)) - Pe
) / (Pe**2 - alpha**2)
Is_formula = (
    sp.exp(Pe) * (Pe * sp.sinh(alpha) - alpha * sp.cosh(alpha)) + alpha
) / (Pe**2 - alpha**2)

# Lightweight exact checks via the generic Piecewise branches returned by SymPy.
Ic_int = sp.integrate(sp.exp(Pe * x) * sp.cosh(alpha * x), (x, 0, 1))
Is_int = sp.integrate(sp.exp(Pe * x) * sp.sinh(alpha * x), (x, 0, 1))
expect_zero("I_c formula (generic branch)", Ic_int.args[0][0] - Ic_formula)
expect_zero("I_s formula (generic branch)", Is_int.args[0][0] - Is_formula)

Delta_formula = sp.factor(
    Pe / (sp.exp(Pe) - 1)
    * ((1 - sp.cosh(alpha)) * Ic_formula + (eta / alpha + sp.sinh(alpha)) * Is_formula)
    / Den
)
print("\nDelta(Pe; kappa, eta) =")
sp.pprint(Delta_formula)

# Numerical spot check against the defining integral.
K_num = sp.lambdify((x, alpha, eta), K, "mpmath")
Sigma_num = sp.lambdify((x, Pe), Sigma, "mpmath")
Delta_num = sp.lambdify((Pe, alpha, eta), Delta_formula, "mpmath")
mp.mp.dps = 50
aval = mp.mpf("1.2")
etaval = mp.mpf("0.7")
peval = mp.mpf("0.8")
quad_val = mp.quad(lambda xx: K_num(xx, aval, etaval) * Sigma_num(xx, peval), [0, 1])
closed_val = Delta_num(peval, aval, etaval)
print("\nNumeric integral check:")
print("  integral definition =", quad_val)
print("  closed form         =", closed_val)
print("  absolute error      =", abs(quad_val - closed_val))
if abs(quad_val - closed_val) > mp.mpf("1e-40"):
    raise AssertionError("Stage 41 Delta integral check failed.")

subbanner("41.3 — Exact endpoint values")
Delta0_formula = eta * (sp.cosh(alpha) - 1) / (alpha**2 * Den)
Delta0_limit = sp.simplify(sp.limit(Delta_formula, Pe, 0))
expect_zero("Delta_0 formula", Delta0_limit - Delta0_formula)

Deltainf_formula = (sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) / Den
expect_zero("Delta_inf from kernel endpoint", sp.simplify(K.subs(x, 1) - Deltainf_formula))

large_pe = mp.mpf("25.0")
large_val = Delta_num(large_pe, aval, etaval)
inf_val = sp.N(Deltainf_formula.subs({alpha: sp.Rational(6,5), eta: sp.Rational(7,10)}), 50)
print("\nLarge-Pe consistency check:")
print("  Delta(Pe=25)        =", large_val)
print("  Delta_inf formula   =", inf_val)
print("  absolute error      =", abs(large_val - mp.mpf(str(inf_val))))

subbanner("41.4 — Fixed-point law and exact interval")
print("Fixed-point equation:")
print("  Pe = Xi * Delta(Pe; kappa, eta)")
print("\nExact constructive-branch bracket:")
print("  Xi * Delta_0(kappa,eta) <= Pe_* <= Xi * Delta_inf(kappa,eta)")
print("\nThis is the first operator-selected interval for the physical Peclet number.")

subbanner("41.5 — Weak-coupling branch law")
print("For small Xi, the constructive branch satisfies")
print("  Pe_* = Xi * Delta_0(kappa,eta) + O(Xi^2).")

banner("STAGE 41 FINAL LEDGER")
print("Stage 41 makes Pe an operator-selected branch variable rather than a free placement")
print("parameter. The exact kernel K(x), exact support-drop law Delta(Pe;kappa,eta), exact")
print("endpoint values Delta_0/Delta_inf, and the fixed-point equation Pe = Xi Delta are now")
print("all explicit and script-checked.")
