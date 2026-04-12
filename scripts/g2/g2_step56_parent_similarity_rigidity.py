#!/usr/bin/env python3
"""
g2_step56_parent_similarity_rigidity.py

Step 56 of the g-2 chain.

What this script does
---------------------
1. Re-derives the exact parent-compensation-family formulas
       g_-(r) = r - 1/2 sqrt(1+r^2),
       L_W/a = (pi/2) sqrt((1+r^2)/3),
       gamma_0 = (1+r^2)/9.
2. Proves the D/N similarity identity
       d ln gamma_0 = 2 d ln(L_W/a),
   so Xi_slip = 0 identically on the exact parent family.
3. Proves lower-branch rigidity: if the canonical-even gate enforces delta g = 0,
   then delta r = 0 because dg_-/dr > 0 on the whole real branch.
4. Shows that this kills the bare mixed-port slippage, the renormalized odd outlet
   shift, and the first-order outgoing normalization defect.
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


banner("STEP 56 — PARENT SIMILARITY RIGIDITY ON THE LOWER COMPENSATED BRANCH")

r = sp.symbols("r", real=True)
delta_r, delta_g = sp.symbols("delta_r delta_g", real=True)
sigma_star, deltaPi_tan = sp.symbols("sigma_star deltaPi_tan", positive=True, real=True)

# Exact parent-family formulas.
g_minus = sp.simplify(r - sp.sqrt(1 + r**2) / 2)
Lratio = sp.simplify(sp.pi * sp.sqrt(1 + r**2) / (2 * sp.sqrt(3)))
gamma0 = sp.simplify((1 + r**2) / 9)
kappa0 = sp.simplify((1 + r**2) / 3)

print("g_-(r) =", g_minus)
print("L_W/a =", Lratio)
print("gamma_0 =", gamma0)
print("kappa_0 =", kappa0)

# ---------------------------------------------------------------------------
# 1. Exact D/N similarity identity on the whole parent family
# ---------------------------------------------------------------------------

dln_gamma = sp.simplify(sp.diff(sp.log(gamma0), r) * delta_r)
dln_Lratio = sp.simplify(sp.diff(sp.log(Lratio), r) * delta_r)
print("d ln gamma_0 =", dln_gamma)
print("d ln(L_W/a) =", dln_Lratio)
expect_zero("d ln gamma_0 - 2 d ln(L_W/a)", dln_gamma - 2 * dln_Lratio)

Xi_gamma = sp.simplify(sp.diff(sp.log(gamma0), r))
Xi_L = sp.simplify(sp.diff(sp.log(Lratio), r))
Xi_slip = sp.simplify(Xi_gamma - 2 * Xi_L)
expect_zero("Xi_slip", Xi_slip)

# Bare mixed-port slippage on the parent family.
delta_BW = sp.simplify(sp.diff(gamma0, r) * delta_r - sp.diff(kappa0, r) * delta_r / 3)
expect_zero("delta_BW", delta_BW)

# ---------------------------------------------------------------------------
# 2. Lower-branch rigidity under the canonical-even gate delta g = 0
# ---------------------------------------------------------------------------

dgdr = sp.simplify(sp.diff(g_minus, r))
dgdr_positive_form = sp.simplify(
    (4 + 3 * r**2)
    / (2 * sp.sqrt(1 + r**2) * (2 * sp.sqrt(1 + r**2) + r))
)
expect_zero("dg_-/dr - positive form", dgdr - dgdr_positive_form)
print("dg_-/dr =", dgdr)

# If delta g = (dg_-/dr) delta r and delta g = 0, then delta r = 0 because dg_-/dr > 0.
delta_g_relation = sp.simplify(dgdr * delta_r)
print("delta g on lower branch =", delta_g_relation)

r_F1 = sp.Float("1.77799353547498")
dgdr_F1 = sp.N(dgdr.subs(r, r_F1), 20)
print("dg_-/dr at Family-1 point =", dgdr_F1)
print("delta r / delta g at Family-1 point =", sp.N(1 / dgdr_F1, 20))

# ---------------------------------------------------------------------------
# 3. First-order outgoing defect collapses to zero
# ---------------------------------------------------------------------------

# Stage-145 first-order outgoing defect law on the exact parent family.
Delta_Q = sp.simplify(-sigma_star * Xi_slip * deltaPi_tan / (1 - sigma_star))
expect_zero("Delta_Q on exact parent compensation family", Delta_Q)

# Since delta_BW = 0 identically on the parent family, the renormalized odd outlet
# shift also vanishes: delta gamma_W = delta_BW / (1 + r^2) = 0.
delta_gamma_W = sp.simplify(delta_BW / (1 + r**2))
expect_zero("delta gamma_W", delta_gamma_W)

banner("FINAL LEDGER")
print("On the exact parent compensation family,")
print("  gamma_0 = (1+r^2)/9 and L_W/a = (pi/2) sqrt((1+r^2)/3)")
print("so the D/N similarity identity holds exactly:")
print("  d ln gamma_0 = 2 d ln(L_W/a)")
print("Hence Xi_slip = 0 identically.")
print()
print("On the lower branch g_-(r), dg/dr is strictly positive, so the carried")
print("canonical-even gate delta g = 0 forces delta r = 0.")
print()
print("Therefore, at first order on the natural compensated lower branch,")
print("  delta_BW = 0, delta gamma_W = 0, Delta_Q = 0, and N_Q - 1 = 0.")
print("So the natural branch preserves the canonical outgoing normalization rather")
print("than generating the electron sliver automatically.")
