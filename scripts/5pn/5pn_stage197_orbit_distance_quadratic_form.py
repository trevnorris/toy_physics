
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage193_198_common import *

"""
Stage 197 — canonical orbit-distance quadratic form.

What this script does
---------------------
1. Defines the exact logarithmic residual coordinates
      t = ln m_T, k = ln m_K, mu = ln m_mu.
2. Builds the exact quotient-coordinate linear map q = A r, where
      q = (q_tr,q_nt,q_eta), r = (t,k,mu).
3. Forms the canonical positive semidefinite scalar
      D^2 = q_tr^2 + q_nt^2 + q_eta^2 = r^T Q r.
4. Proves that Q is positive definite on the constructive branch by exact principal minors.
5. Gives the pure-channel orbit distances and the local Hessian weights at orbit lock.
"""

banner("STAGE 197 — CANONICAL ORBIT-DISTANCE QUADRATIC FORM")

star = base_symbols("")
chi0s = star["chi0s"]
Fstar = star["Fstar"]
res = residual_log_symbols()
A = orbit_distance_quadratic_form(star)["A"]
Q = orbit_distance_quadratic_form(star)["Q"]

subbanner("I. Exact linear map from residual logs to quotient coordinates")
print("A =")
sp.pprint(A)

qtr = sp.simplify((1 + chi0s) * res["t"])
qnt = sp.simplify(res["m"] - res["k"] - Fstar * res["t"])
qeta = sp.simplify(-res["k"])

expect_zero("A r - q", A * sp.Matrix([res["t"], res["k"], res["m"]]) - sp.Matrix([qtr, qnt, qeta]))

subbanner("II. Canonical scalar distance D^2 = r^T Q r")
print("Q = A^T A =")
sp.pprint(Q)

D2 = sp.simplify(
    qtr**2 + qnt**2 + qeta**2
    - (sp.Matrix([res["t"], res["k"], res["m"]]).T * Q * sp.Matrix([res["t"], res["k"], res["m"]]))[0]
)
expect_zero("D^2 identity", D2)

print("Expanded D^2 =")
sp.pprint(sp.expand(qtr**2 + qnt**2 + qeta**2))

subbanner("III. Exact positivity")
minor1 = sp.simplify(Q[:1,:1].det())
minor2 = sp.simplify(Q[:2,:2].det())
minor3 = sp.simplify(Q.det())

print("Principal minor 1 =")
sp.pprint(minor1)
print("Principal minor 2 =")
sp.pprint(minor2)
print("Principal minor 3 =")
sp.pprint(minor3)

# Symbolic positivity statements on the constructive branch.
print("\nOn chi0_* > 0,")
print("  minor1 = (1+chi0_*)^2 + F_*^2 > 0,")
print("  minor2 = 2(1+chi0_*)^2 + F_*^2 > 0,")
print("  minor3 = (1+chi0_*)^2 > 0.")
print("So Q is positive definite and D^2 = 0 iff t = k = mu = 0.")

subbanner("IV. Pure-channel orbit distances")
D2_T = sp.simplify((q_from_mismatch(star, sp.exp(res["t"]), 1, 1)["qtr"])**2 + (q_from_mismatch(star, sp.exp(res["t"]), 1, 1)["qnt"])**2)
D2_K = sp.simplify((q_from_mismatch(star, 1, sp.exp(res["k"]), 1)["qnt"])**2 + (q_from_mismatch(star, 1, sp.exp(res["k"]), 1)["qeta"])**2)
D2_M = sp.simplify((q_from_mismatch(star, 1, 1, sp.exp(res["m"]))["qnt"])**2)

print("Pure T_U channel D_T^2 =")
sp.pprint(D2_T)
print("Pure K_eta channel D_K^2 =")
sp.pprint(D2_K)
print("Pure mu_W channel D_mu^2 =")
sp.pprint(D2_M)

banner("STAGE 197 LEDGER")
print("1. The exact logarithmic residual triple carries a canonical positive-definite orbit")
print("   distance D^2 = q_tr^2 + q_nt^2 + q_eta^2.")
print("2. This distance vanishes iff the pairwise orbit-lock condition holds.")
print("3. So the orbit-lock failure of an actual branch can be quantified by one scalar that")
print("   is exactly reference-independent once written in pairwise residual coordinates.")
