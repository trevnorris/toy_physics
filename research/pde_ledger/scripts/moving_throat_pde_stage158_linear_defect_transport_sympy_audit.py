#!/usr/bin/env python3
"""
moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py

SymPy-backed audit for Stage 158.

Checks:
1. Linear transport delta R(delta g) about the lower compensated Family-1 branch.
2. Linear gain transport for M_s = Sigma0 and M_q = -Sigma0 R.
3. Linear slope transport for Pi = Sigma0 (1 - R S).
4. Linearized outgoing-normalization defect Delta_Q = 5 b + a0/3 + 9 a5.
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


banner("STAGE 158 — LINEAR DEFECT TRANSPORT AROUND THE RENORMALIZED CANONICAL POINT")

# ---------------------------------------------------------------------------
# 1. delta R from delta g
# ---------------------------------------------------------------------------
g, r, dg = sp.symbols("g r dg", real=True)
g_star = r - sp.sqrt(1 + r**2) / 2
R = (g - r) ** 2 / (1 + r**2)
R_shift = sp.expand(R.subs(g, g_star + dg))
R_lin = sp.expand(sp.series(R_shift, dg, 0, 2).removeO())
R_expected = sp.expand(sp.Rational(1, 4) - dg / sp.sqrt(1 + r**2))
expect_zero("linear delta R law", R_lin - R_expected)

# ---------------------------------------------------------------------------
# 2. Mouth gains
# ---------------------------------------------------------------------------
Sigma0, dSigma0, Rstar, dR = sp.symbols("Sigma0 dSigma0 Rstar dR", real=True)
Mq = -(Sigma0 + dSigma0) * (Rstar + dR)
Mq_lin = sp.expand(Mq).subs({dSigma0 * dR: 0})
Mq0 = -Sigma0 * Rstar
expect_zero("delta Mq law", (Mq_lin - Mq0) - (-Rstar * dSigma0 - Sigma0 * dR))

# ---------------------------------------------------------------------------
# 3. Mouth slope / bias transport
# ---------------------------------------------------------------------------
Sstar, dS = sp.symbols("Sstar dS", real=True)
Pi = (Sigma0 + dSigma0) * (1 - (Rstar + dR) * (Sstar + dS))
Pi_lin = sp.expand(Pi).subs({dSigma0 * dR: 0, dSigma0 * dS: 0, dR * dS: 0})
Pi0 = Sigma0 * (1 - Rstar * Sstar)
dPi_expected = (1 - Rstar * Sstar) * dSigma0 - Sigma0 * (Rstar * dS + Sstar * dR)
expect_zero("delta Pi law", (Pi_lin - Pi0) - dPi_expected)

# ---------------------------------------------------------------------------
# 3b. Composed boxed identities (notes §3-§4)
# ---------------------------------------------------------------------------
dg_sym, r_sym = sp.symbols("dg_sym r_sym", real=True, positive=True)
dR_from_dg = -dg_sym / sp.sqrt(1 + r_sym**2)

dMq_composed = -sp.Rational(1, 4) * dSigma0 - Sigma0 * dR_from_dg
dMq_boxed = -sp.Rational(1, 4) * dSigma0 + Sigma0 / sp.sqrt(1 + r_sym**2) * dg_sym
expect_zero("composed delta Mq law", sp.expand(dMq_composed - dMq_boxed))

dPi_composed = (1 - sp.Rational(1, 4) * Sstar) * dSigma0 \
    - Sigma0 * (sp.Rational(1, 4) * dS + Sstar * dR_from_dg)
dPi_boxed = (1 - Sstar / 4) * dSigma0 \
    - (Sigma0 / 4) * dS \
    + (Sigma0 * Sstar) / sp.sqrt(1 + r_sym**2) * dg_sym
expect_zero("composed delta Pi law", sp.expand(dPi_composed - dPi_boxed))

# ---------------------------------------------------------------------------
# 4. Linear outgoing-normalization defect
# ---------------------------------------------------------------------------
eps, s, b, a0, a5 = sp.symbols("eps s b a0 a5", real=True)
S = 1 + eps * s
beta = 1 + eps * b
Sigma_0 = eps * a0
Sigma_5 = eps * a5
chi = sp.simplify(3 * (S * beta**5 + 9 * Sigma_5) / (3 * S - Sigma_0))
chi_lin = sp.expand(sp.series(chi, eps, 0, 2).removeO())
chi_expected = 1 + eps * (5 * b + a0 / 3 + 9 * a5)
expect_zero("linear Delta_Q law", chi_lin - chi_expected)

# ---------------------------------------------------------------------------
# 5. Numerical coefficients at the Stage 156 point
# ---------------------------------------------------------------------------
banner("Numerical coefficients at the renormalized canonical point")

rF1 = sp.Float("1.77799353547498")
Sigma0_can = sp.Float("4.651033550168876")
S_can = sp.Float("0.6703621156734617")
T_can = sp.Float("1.4467083664567624")
sqrt1 = sp.sqrt(1 + rF1**2)

coef_dR_dg = sp.N(-1 / sqrt1, 20)
coef_dMq_dSigma = sp.Rational(-1, 4)
coef_dMq_dg = sp.N(Sigma0_can / sqrt1, 20)
coef_dPi_dSigma = sp.N(1 - S_can / 4, 20)
coef_dPi_dS = sp.N(-Sigma0_can / 4, 20)
coef_dPi_dg = sp.N(Sigma0_can * S_can / sqrt1, 20)
coef_dSigma_dT = sp.N(sp.Rational(40, 9) * T_can, 20)
coef_dPi_dT = sp.N(coef_dPi_dSigma * coef_dSigma_dT, 20)

print("dR/dg        =", coef_dR_dg)
print("dMq/dSigma0  =", coef_dMq_dSigma)
print("dMq/dg       =", coef_dMq_dg)
print("dPi/dSigma0  =", coef_dPi_dSigma)
print("dPi/dS       =", coef_dPi_dS)
print("dPi/dg       =", coef_dPi_dg)
print("dSigma0/dThat=", coef_dSigma_dT)
print("dPi/dThat    =", coef_dPi_dT)

print("\nCarry-forward summary:")
print("  delta R  = -(delta g)/sqrt(1+r_F1^2) + O(delta g^2)")
print("  delta Mq = -(1/4) delta Sigma0 + Sigma0_can/sqrt(1+r_F1^2) delta g + O(2)")
print("  delta Pi = (1-S_can/4) delta Sigma0 - (Sigma0_can/4) delta S + Sigma0_can*S_can/sqrt(1+r_F1^2) delta g + O(2)")
print("  Delta_Q  = 5 b + a0/3 + 9 a5 + O(2)")
