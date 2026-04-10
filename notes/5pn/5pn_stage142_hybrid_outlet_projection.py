
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 142 — HYBRID OUTLET PROJECTION")

r = sp.symbols('r', positive=True, real=True)
Sigma0, sigma = sp.symbols('Sigma0 sigma', positive=True, real=True)
Xi = sp.symbols('Xi', positive=True, real=True)
dSigma0, dg, dkappaW, dgammaW = sp.symbols('deltaSigma0 deltag deltaKappaW deltaGammaW', real=True)
deltaC = sp.symbols('deltaC', real=True)

root = root1p(r)

subbanner("1. Direct loading defect on the compensated hybrid surface")
# Use Stage-141 linearized M_s and M_q
dR = -dg / root
dMs = dSigma0
dMq = -sp.Rational(1, 4) * dSigma0 + Sigma0 * dg / root

rhoR = Xi * dMs
sigmaW = -Xi * dMq
dC = sp.simplify(rhoR - 4 * sigmaW)
sigma_star = Xi * Sigma0 / 4
target_dC = sp.simplify(-16 * sigma_star * dR)
expect_zero("deltaC - (-16 sigma_* deltaR)", dC - target_dC)
print("deltaC =", dC)

subbanner("2. Exact linear outlet algebra")
deltaE2 = (deltaC - 9 * sigma * dkappaW) / (27 * (1 - sigma))
deltaE4 = (5 * deltaC - 72 * sigma * dkappaW) / (243 * (1 - sigma))
DeltaQ = (deltaC - 27 * sigma * dgammaW) / (3 * (1 - sigma))
print("deltaE2 =", deltaE2)
print("deltaE4 =", deltaE4)
print("DeltaQ =", DeltaQ)

subbanner("3. Canonical-even gate")
sol = sp.solve([sp.Eq(deltaE2, 0), sp.Eq(deltaE4, 0)], [deltaC, dkappaW], dict=True)
print("solution to deltaE2 = deltaE4 = 0:", sol)
expect_zero("canonical-even deltaC", sol[0][deltaC])
expect_zero("canonical-even deltaKappaW", sol[0][dkappaW])

# If sigma != 0 and deltaC = 16 sigma_* dg / root, then delta g must vanish.
print("On a nontrivial loaded branch sigma_* != 0, canonical-even preservation implies delta g = 0.")

subbanner("4. Collapse of the final linear 2.5PN defect")
DeltaQ_even = sp.simplify(DeltaQ.subs({deltaC: 0, dkappaW: 0}))
target_even = -9 * sigma * dgammaW / (1 - sigma)
expect_zero("DeltaQ_even - target", DeltaQ_even - target_even)
print("DeltaQ under the canonical-even gate =", DeltaQ_even)
print("N_Q - 1 =", -DeltaQ_even)

subbanner("5. Numeric coefficient on the Family-1 branch")
numsubs = {r: rF1}
print("deltaC coefficient / sigma_* =", sp.N((16 / root).subs(numsubs), 16))
