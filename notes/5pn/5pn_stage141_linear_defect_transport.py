
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 141 — LINEAR DEFECT TRANSPORT FROM THE RENORMALIZED FAMILY-1 POINT")

r = sp.symbols('r', positive=True, real=True)
Sigma0, S0, Xi = sp.symbols('Sigma0 S0 Xi', positive=True, real=True)
dSigma0, dg, dS = sp.symbols('deltaSigma0 deltag deltaS', real=True)
b, a0, a5 = sp.symbols('b a0 a5', real=True)

root = root1p(r)
gstar = lower_branch_g(r)
g = sp.symbols('g', real=True)

subbanner("1. Exact transport delta g -> delta R")
R = ((g - r) ** 2) / (1 + r**2)
dR = sp.simplify(sp.diff(R, g).subs(g, gstar) * dg)
expect_zero("deltaR + dg/sqrt(1+r^2)", dR + dg / root)
print("deltaR =", dR)

subbanner("2. Mouth gains Ms, Mq and traction transport")
Ms = Sigma0
Mq = -Sigma0 * R
dMq = sp.simplify(
    sp.diff(Mq, Sigma0).subs(g, gstar) * dSigma0
    + sp.diff(Mq, g).subs(g, gstar) * dg
)
target_dMq = -sp.Rational(1, 4) * dSigma0 + Sigma0 * dg / root
expect_zero("deltaMq - target", dMq - target_dMq)
print("deltaMs =", dSigma0)
print("deltaMq =", dMq)

dSigma0_from_dThat = sp.Symbol('alpha_T')  # placeholder for linear traction coefficient
print("Stage-141 traction law template: deltaSigma0 = alpha_T * delta T_hat")

subbanner("3. Slope / mouth-bias transport")
Pi = Sigma0 * (1 - R * S0)
dPi = sp.simplify(
    sp.diff(Pi, Sigma0).subs(g, gstar) * dSigma0
    + sp.diff(Pi, g).subs(g, gstar) * dg
    + sp.diff(Pi, S0).subs(g, gstar) * dS
)
target_dPi = (1 - sp.Rational(1, 4) * S0) * dSigma0 - Sigma0 * dS / 4 + Sigma0 * S0 * dg / root
expect_zero("deltaPi - target", dPi - target_dPi)
print("deltaPi =", dPi)

subbanner("4. Reduced 2.5PN bridge")
DeltaQ = 5 * b + a0 / 3 + 9 * a5
print("Delta_Q =", DeltaQ)
print("N_Q - 1 = -Delta_Q + O(Delta_Q^2) on the natural source-map branch")

subbanner("5. Numeric Family-1 coefficients")
numsubs = {
    r: rF1,
    Sigma0: Sigma0_can,
    S0: S_can,
}
print("deltaR coefficient =", sp.N((-1 / root).subs(numsubs), 16))
print("deltaMq coefficient of delta g =", sp.N((Sigma0 / root).subs(numsubs), 16))
print("deltaPi coefficient of deltaSigma0 =", sp.N((1 - S0 / 4).subs(numsubs), 16))
print("deltaPi coefficient of deltaS =", sp.N((-Sigma0 / 4).subs(numsubs), 16))
print("deltaPi coefficient of delta g =", sp.N((Sigma0 * S0 / root).subs(numsubs), 16))
