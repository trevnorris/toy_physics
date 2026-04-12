#!/usr/bin/env python3
"""
g2_step55_compensated_branch_rigidity_conclusion.py

Step 55 of the g-2 chain.

What this script does
---------------------
1. Starts from the compensated hybrid canonical-even outlet law for the residual outgoing defect.
2. Solves the electron target for the required odd mixed-channel renormalization delta gamma_W.
3. Applies the bare mixed-port slippage theorem delta gamma_W = deltaB_W / (1 + r_c*).
4. Applies the similarity-slippage law Delta_Q = -(sigma_*/(1-sigma_*)) Xi_slip deltaPi_tan.
5. Shows that on the exact parent compensation family Xi_slip = 0, hence Delta_Q = 0 at first order.
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


banner("STEP 55 — COMPENSATED BRANCH RIGIDITY AND THE LAST LIVE SLIPPAGE SCALAR")

sigma_star, rc_star = sp.symbols("sigma_star rc_star", positive=True, real=True)
delta_gamma_W, deltaB_W = sp.symbols("delta_gamma_W deltaB_W", real=True)
Xi_slip, deltaPi_tan = sp.symbols("Xi_slip deltaPi_tan", real=True)
delta = sp.symbols("delta", positive=True, real=True)

# Electron-point outgoing defect from the g-2 chain:
Delta_e = sp.simplify(-delta / (1 + delta))
print("Delta_e =", Delta_e)

# ---------------------------------------------------------------------------
# 1. Compensated hybrid canonical-even outlet law
# ---------------------------------------------------------------------------

Delta_Q_hyb = sp.simplify(-9 * sigma_star * delta_gamma_W / (1 - sigma_star))
print("Delta_Q(hybrid) =", Delta_Q_hyb)

delta_gamma_req = sp.solve(sp.Eq(Delta_Q_hyb, Delta_e), delta_gamma_W)[0]
print("delta_gamma_W required for electron target =", sp.simplify(delta_gamma_req))

# ---------------------------------------------------------------------------
# 2. Bare mixed-port slippage theorem
# ---------------------------------------------------------------------------

delta_gamma_from_slip = sp.simplify(deltaB_W / (1 + rc_star))
print("delta_gamma_W(slippage) =", delta_gamma_from_slip)

deltaB_req = sp.solve(sp.Eq(delta_gamma_from_slip, delta_gamma_req), deltaB_W)[0]
print("deltaB_W required for electron target =", sp.simplify(deltaB_req))

# ---------------------------------------------------------------------------
# 3. Similarity-slippage law
# ---------------------------------------------------------------------------

Delta_Q_similarity = sp.simplify(-sigma_star * Xi_slip * deltaPi_tan / (1 - sigma_star))
print("Delta_Q(similarity) =", Delta_Q_similarity)

XiPi_req = sp.solve(sp.Eq(Delta_Q_similarity, Delta_e), Xi_slip * deltaPi_tan)[0]
print("Xi_slip * deltaPi_tan required for electron target =", sp.simplify(XiPi_req))

# ---------------------------------------------------------------------------
# 4. Exact parent compensation-family rigidity
# ---------------------------------------------------------------------------

# On the exact parent compensation family the later moving-throat notes give
# Xi_slip = Xi_gamma - 2 Xi_L = 0 identically.
Xi_slip_parent = 0
Delta_parent = sp.simplify(Delta_Q_similarity.subs(Xi_slip, Xi_slip_parent))
print("Delta_Q on exact parent compensation family =", Delta_parent)
if Delta_parent != 0:
    raise AssertionError("Parent compensation-family rigidity did not give Delta_Q = 0.")

banner("FINAL LEDGER")
print("On the compensated hybrid canonical-even branch,")
print("  Delta_Q = -(9 sigma_*/(1-sigma_*)) delta_gamma_W.")
print("So the electron-point sliver fixes one odd renormalization:")
print("  delta_gamma_W = {}".format(sp.simplify(delta_gamma_req)))
print()
print("Using the bare mixed-port slippage theorem, this is equivalent to")
print("  deltaB_W = {}".format(sp.simplify(deltaB_req)))
print()
print("Using the similarity-slippage form, it is also equivalent to")
print("  Xi_slip * deltaPi_tan = {}".format(sp.simplify(XiPi_req)))
print()
print("But on the exact parent compensation family Xi_slip = 0 identically, so")
print("  Delta_Q = 0 at first order.")
print("Therefore the natural compensated branch does not generate the electron")
print("sliver automatically; a nonzero sliver would require a genuine departure")
print("from the parent compensation family / pure-scale lock.")
