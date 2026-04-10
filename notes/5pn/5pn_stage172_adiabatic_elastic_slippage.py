#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = '-' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        simp = expr.applyfunc(lambda e: sp.simplify(sp.expand(e)))
        print(f"{name} =")
        sp.pprint(simp)
        if any(e != 0 for e in simp):
            raise AssertionError(f"{name} is not zero")
    else:
        simp = sp.simplify(sp.expand(expr))
        print(f"{name} = {simp}")
        if simp != 0:
            raise AssertionError(f"{name} is not zero")


def lower_branch_g(r):
    return sp.simplify(r - sp.sqrt(1 + r**2) / 2)


def root1p(r):
    return sp.sqrt(1 + r**2)


"""
5pn_stage172_adiabatic_elastic_slippage.py

Stage 172 — rewrite the Stage-151 off-bundle slippages on the adiabatic wall
and impose the elastic/no-fraying branch law.

What this script does
---------------------
1. Specializes the Stage-151 slippage definitions to the adiabatic wall
      delta ln Theta_w = 0.
2. Shows the exact reduced slippage formulas
      epsilon_L, epsilon_v, epsilon_T
   relative to the adiabatic lower-branch transport law.
3. Builds the weighted scalar
      epsilon_perp
   that carries the off-family normal motion.
4. Proves that an elastic/no-fraying boundary law
      epsilon_L = epsilon_v = epsilon_T = 0
   forces epsilon_perp = 0 exactly.

Interpretation
--------------
The adiabatic wall removes isotropic thermal-fraying drift, while the elastic
boundary law removes the first scalar off-bundle source entirely. Any remaining
first-order obstruction must therefore come from the direct grouped outlet side
or from quotient motion in the Stage-169/170 orbit variables.
"""

banner("STAGE 172 — ADIABATIC ELASTIC OFF-BUNDLE SLIPPAGE")

# Adiabatic isotropic bundle observables

dlnKs, dlnKq, dlnP0 = sp.symbols('dlnKs dlnKq dlnP0', real=True)

# Actual first-order drifts of the transported variables

dlnLW, dlnvw0, dlnTm = sp.symbols('dlnLW dlnvw0 dlnTm', real=True)

epsL, epsv, epsT = sp.symbols('epsilon_L epsilon_v epsilon_T', real=True)
eps_perp = sp.symbols('epsilon_perp', real=True)

r = sp.symbols('r', positive=True, real=True)
gstar = lower_branch_g(r)

subbanner("1. Adiabatic lower-branch transport laws")

LW_expected = sp.Rational(1, 2) * dlnKs
vw0_expected = -sp.Rational(3, 4) * dlnKs + sp.Rational(1, 2) * dlnKq
Tm_expected = -sp.Rational(5, 4) * dlnKs + sp.Rational(1, 2) * dlnKq - sp.Rational(2, 5) * dlnP0

print("Expected adiabatic lower-branch laws:")
print("delta ln L_W =", LW_expected)
print("delta ln v_{w0} =", vw0_expected)
print("delta ln T_m =", Tm_expected)

subbanner("2. Exact adiabatic slippage variables")

# Stage-151 slippages specialized to dlnTheta = 0
slip_defs = {
    epsL: dlnLW - LW_expected,
    epsv: dlnvw0 - vw0_expected,
    epsT: dlnTm - Tm_expected,
}

print("epsilon_L =", slip_defs[epsL])
print("epsilon_v =", slip_defs[epsv])
print("epsilon_T =", slip_defs[epsT])

subbanner("3. Weighted scalar epsilon_perp on the Family-1 compensated branch")

coef_T = sp.simplify(gstar)
coef_v = sp.simplify(gstar + 1 / (2 * root1p(r)))
coef_L = sp.simplify(2 * gstar + 3 / (4 * root1p(r)))

eps_perp_expr = sp.simplify(coef_T * epsT + coef_v * epsv + coef_L * epsL)
print("epsilon_perp =")
sp.pprint(eps_perp_expr)

subbanner("4. Elastic/no-fraying boundary law kills the scalar off-bundle source")

elastic = {epsL: 0, epsv: 0, epsT: 0}
expect_zero("epsilon_perp under elastic/no-fraying boundary", eps_perp_expr.subs(elastic))

# Equivalently, imposing the exact adiabatic lower-branch transport law on the actual drifts
elastic_transport = {
    dlnLW: LW_expected,
    dlnvw0: vw0_expected,
    dlnTm: Tm_expected,
}
expect_zero(
    "epsilon_perp after substituting exact adiabatic transport",
    eps_perp_expr.subs({
        epsL: slip_defs[epsL].subs(elastic_transport),
        epsv: slip_defs[epsv].subs(elastic_transport),
        epsT: slip_defs[epsT].subs(elastic_transport),
    })
)

subbanner("5. Numerical Family-1 weights")

rF1 = sp.Float('1.77799353547498')
gF1 = sp.Float('0.758035078944663')
coef_T_num = sp.N(coef_T.subs(r, rF1), 16)
coef_v_num = sp.N(coef_v.subs(r, rF1), 16)
coef_L_num = sp.N(coef_L.subs(r, rF1), 16)
print("g_* =", gF1)
print("coef_T ≈", coef_T_num)
print("coef_v ≈", coef_v_num)
print("coef_L ≈", coef_L_num)

banner("FINAL STAGE-172 LEDGER")
print("1. On the adiabatic wall, the Stage-151 scalar slippages reduce to")
print("      epsilon_L = delta ln L_W - 1/2 delta ln K_s,")
print("      epsilon_v = delta ln v_{w0} + 3/4 delta ln K_s - 1/2 delta ln K_q,")
print("      epsilon_T = delta ln T_m + 5/4 delta ln K_s - 1/2 delta ln K_q + 2/5 delta ln P_0.")
print("2. The full scalar off-family source is still one weighted combination")
print("      epsilon_perp = g_* epsilon_T + (g_* + 1/(2 sqrt(1+r_*^2))) epsilon_v + (2 g_* + 3/(4 sqrt(1+r_*^2))) epsilon_L.")
print("3. If the topological boundary deforms elastically with no thermal fraying, so that")
print("      epsilon_L = epsilon_v = epsilon_T = 0,")
print("   then epsilon_perp = 0 exactly.")
print("4. Therefore the adiabatic-elastic boundary law removes the first scalar off-bundle source entirely.")
print("5. Any remaining first-order obstruction must come from the direct grouped outlet side or from Stage-169/170 quotient motion.")
