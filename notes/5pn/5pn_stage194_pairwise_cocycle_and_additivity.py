
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage193_198_common import *

"""
Stage 194 — pairwise cocycle laws and additive quotient coordinates.

What this script does
---------------------
1. Introduces three arbitrary states x, y, z on the same positive microscopic state space.
2. Proves that the exact pairwise residual mismatch ratios compose multiplicatively:
      m_(x->z) = m_(x->y) m_(y->z).
3. Proves that the logarithmic quotient coordinates compose additively:
      q_(x->z) = q_(x->y) + q_(y->z).
4. Shows the inverse laws under swapping the order of comparison.
5. Concludes that the residual triple forms an exact multiplicative cocycle over the
   microscopic branch trajectory, while the quotient coordinates are its additive log chart.
"""

banner("STAGE 194 — PAIRWISE COCYCLE LAWS AND ADDITIVE QUOTIENT COORDINATES")

star = base_symbols("")
x = base_symbols("x_")
y = base_symbols("y_")
z = base_symbols("z_")
for state in (x, y, z):
    state["chi0s"] = star["chi0s"]
    state["deltaUs"] = star["deltaUs"]
    state["Estar"] = star["Estar"]
    state["Fstar"] = star["Fstar"]
    state["sigma"] = star["sigma"]
    state["L"] = star["L"]

xy = pairwise_mismatch_ratios(x, y, star)
yz = pairwise_mismatch_ratios(y, z, star)
xz = pairwise_mismatch_ratios(x, z, star)

subbanner("I. Exact multiplicative cocycle law for mismatch ratios")
expect_zero("ln m_T(x->z) - ln m_T(x->y) - ln m_T(y->z)", sp.log(xz["mT"]) - sp.log(xy["mT"]) - sp.log(yz["mT"]))
expect_zero("ln m_K(x->z) - ln m_K(x->y) - ln m_K(y->z)", sp.log(xz["mK"]) - sp.log(xy["mK"]) - sp.log(yz["mK"]))
expect_zero("ln m_mu(x->z) - ln m_mu(x->y) - ln m_mu(y->z)", sp.log(xz["mMu"]) - sp.log(xy["mMu"]) - sp.log(yz["mMu"]))

print("m_T(x->z) = m_T(x->y) m_T(y->z)")
print("m_K(x->z) = m_K(x->y) m_K(y->z)")
print("m_mu(x->z) = m_mu(x->y) m_mu(y->z)")

subbanner("II. Additive law for exact quotient coordinates")
q_xy = q_from_mismatch(star, xy["mT"], xy["mK"], xy["mMu"])
q_yz = q_from_mismatch(star, yz["mT"], yz["mK"], yz["mMu"])
q_xz = q_from_mismatch(star, xz["mT"], xz["mK"], xz["mMu"])

expect_zero("q_tr(x->z) - q_tr(x->y) - q_tr(y->z)", q_xz["qtr"] - q_xy["qtr"] - q_yz["qtr"])
expect_zero("q_nt(x->z) - q_nt(x->y) - q_nt(y->z)", q_xz["qnt"] - q_xy["qnt"] - q_yz["qnt"])
expect_zero("q_eta(x->z) - q_eta(x->y) - q_eta(y->z)", q_xz["qeta"] - q_xy["qeta"] - q_yz["qeta"])

print("q_(x->z) = q_(x->y) + q_(y->z)")

subbanner("III. Inverse laws")
yx = pairwise_mismatch_ratios(y, x, star)
expect_zero("ln m_T(x->y) + ln m_T(y->x)", sp.log(xy["mT"]) + sp.log(yx["mT"]))
expect_zero("ln m_K(x->y) + ln m_K(y->x)", sp.log(xy["mK"]) + sp.log(yx["mK"]))
expect_zero("ln m_mu(x->y) + ln m_mu(y->x)", sp.log(xy["mMu"]) + sp.log(yx["mMu"]))

q_yx = q_from_mismatch(star, yx["mT"], yx["mK"], yx["mMu"])
expect_zero("q_tr(y->x) + q_tr(x->y)", q_yx["qtr"] + q_xy["qtr"])
expect_zero("q_nt(y->x) + q_nt(x->y)", q_yx["qnt"] + q_xy["qnt"])
expect_zero("q_eta(y->x) + q_eta(x->y)", q_yx["qeta"] + q_xy["qeta"])

banner("STAGE 194 LEDGER")
print("1. The residual mismatch triple is an exact multiplicative cocycle over branch comparisons.")
print("2. The quotient coordinates are the exact additive logarithmic chart of that cocycle.")
print("3. So for a sequence of PDE snapshots, orbit-lock transport can be accumulated exactly")
print("   either multiplicatively in (m_T,m_K,m_mu) or additively in (q_tr,q_nt,q_eta).")
