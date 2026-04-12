
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 146 — OFF-FAMILY NORMAL COORDINATE")

r = sp.symbols('r', positive=True, real=True)
dg, dr = sp.symbols('deltag deltar', real=True)
g = sp.symbols('g', real=True)
root = root1p(r)
gstar = lower_branch_g(r)
gprime = lower_branch_gprime(r)
delta_perp = dg - gprime * dr

subbanner("1. Parent compensation-defect scalar")
F = 1 + r**2 - 4 * (g - r)**2
dF = sp.simplify(sp.diff(F, g).subs(g, gstar) * dg + sp.diff(F, r).subs(g, gstar) * dr)
expect_zero("deltaF - 4 sqrt(1+r^2) delta_perp", dF - 4 * root * delta_perp)
print("deltaF =", dF)

subbanner("2. Transport into the compensation ratio")
Rq = (g - r)**2 / (1 + r**2)
dRq = sp.simplify(sp.diff(Rq, g).subs(g, gstar) * dg + sp.diff(Rq, r).subs(g, gstar) * dr)
expect_zero("deltaRq + delta_perp/sqrt(1+r^2)", dRq + delta_perp / root)
print("deltaRq =", dRq)

subbanner("3. Exact microscopic formula for the off-family scalar")
dln_lambda, dlnKs, dlnKq = sp.symbols('dlnlambda dlnKs dlnKq', real=True)
dln_gs, dln_gq = sp.symbols('dln_gs dln_gq', real=True)
delta_r = sp.simplify(r * (dln_lambda - sp.Rational(1, 2) * dlnKs - sp.Rational(1, 2) * dlnKq))
delta_g = sp.simplify(gstar * (dln_gq - dln_gs + sp.Rational(1, 2) * dlnKs - sp.Rational(1, 2) * dlnKq))
delta_perp_expr = sp.simplify(delta_g - gprime * delta_r)
target = sp.simplify(
    gstar * (dln_gq + dlnKs - dln_gs - dln_lambda)
    + (dlnKs + dlnKq - 2 * dln_lambda) / (4 * root)
)
expect_zero("delta_perp microscopic identity", delta_perp_expr - target)
print("delta_perp =", target)

subbanner("4. Tangent/normal decomposition of the mouth-bias transport")
Sigma0, S0, dSigma0, dS = sp.symbols('Sigma0 S0 deltaSigma0 deltaS', real=True)
dPi = (1 - S0 / 4) * dSigma0 - Sigma0 * dS / 4 + Sigma0 * S0 * delta_perp / root
print("deltaPi =", dPi)
print("Family-1 coefficient of delta_perp =", sp.N((Sigma0_can * S_can / root).subs(r, rF1), 16))
