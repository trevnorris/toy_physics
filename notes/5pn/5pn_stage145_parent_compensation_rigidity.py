
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage141_150_common import *

banner("STAGE 145 — PARENT COMPENSATION RIGIDITY")

r, dr = sp.symbols('r dr', positive=True, real=True)
root = root1p(r)
g = lower_branch_g(r)
gprime = lower_branch_gprime(r)

subbanner("1. Exact parent compensation family")
F = 1 + r**2 - 4 * (sp.Symbol('g') - r)**2
print("F(g,r) =", F)
print("g_-(r) =", g)

subbanner("2. Automatic similarity preservation on the parent family")
LW_over_a = sp.pi / 2 * sp.sqrt((1 + r**2) / 3)
gamma0 = (1 + r**2) / 9
dln_gamma0 = sp.simplify(sp.diff(sp.log(gamma0), r) * dr)
dln_Lratio = sp.simplify(sp.diff(sp.log(LW_over_a), r) * dr)
expect_zero("delta ln gamma0 - 2 delta ln(L_W/a)", dln_gamma0 - 2 * dln_Lratio)
print("delta ln gamma0 =", dln_gamma0)
print("delta ln(L_W/a) =", dln_Lratio)

subbanner("3. Lower-branch rigidity under delta g = 0")
positive_form = sp.simplify((4 + 3 * r**2) / (2 * root * (2 * root + r)))
expect_zero("gprime - positive form", gprime - positive_form)
print("g_-'(r) =", gprime)
print("Positive factor form =", positive_form)
print("Since g_-'(r) > 0 for all real r, delta g = 0 implies delta r = 0 on the lower branch.")

subbanner("4. Collapse of all first-order D/N similarity defects")
print("If delta r = 0, then delta r_c = 0, delta ln(L_W/a) = 0, delta ln gamma0 = 0,")
print("delta kappa_0 = 0, delta gamma_0 = 0, delta gamma_W = 0, and therefore Delta_Q = 0.")
print("Family-1 value of g_-'(r_F1) =", sp.N(gprime.subs(r, rF1), 16))
print("Inverse rigidity coefficient dr/dg =", sp.N((1 / gprime).subs(r, rF1), 16))
