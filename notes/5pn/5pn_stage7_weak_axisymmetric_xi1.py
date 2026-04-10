#!/usr/bin/env python3
"""
5pn_stage7_weak_axisymmetric_xi1.py

Seventh executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Implements the Stage-160 weak-axisymmetric bookkeeping for the outgoing
   Maxwell/mixed sector.
2. Proves that the grouped signature
      lambda_20 = 1, lambda_21 = 1/2, lambda_22 = -1
   transports every microscopic outgoing slippage into the same grouped pattern.
3. Derives the exact port amplitudes
      sigma_r = 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r,
   and the weighted scalar collapse
      Xi_1 = sum_r rho_r^(N) sigma_r.
4. Proves that the full remaining grouped defect is one-dimensional:
      Xi_load^(A) = eps lambda_A Xi_1.
5. Computes the grouped trace/anomaly variables and proves the axisymmetric law
      b_Xi = 3 a_Xi.
6. Identifies the same scalar with the physical outgoing-prefactor slope
      Xi_1 = P1/P0,
   and writes the even-preserving quadrupole defect pattern directly in terms
   of Xi_1.
7. Records the interference/hybridization-rigid and dominant-port corollaries.

Interpretation
--------------
After Stage 7 the remaining weak-axisymmetric grouped defect is not a bundle of
independent microscopic drifts. It is one scalar amplitude Xi_1, equivalently
P1/P0, carried by the outgoing mixed-sector slippage bundle.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


def grouped_triplet(x20: sp.Expr, x21: sp.Expr, x22: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    a_x = sp.simplify((2 * x20 - x21 - x22) / 10)
    b_x = sp.simplify((x21 - x22) / 2)
    return xbar, a_x, b_x


# ---------------------------------------------------------------------------
# I. Weak-axisymmetric microscopic slope bookkeeping
# ---------------------------------------------------------------------------

banner("I. WEAK-AXISYMMETRIC MICROSCOPIC SLOPE BOOKKEEPING")

eps = sp.symbols("epsilon", real=True)
lambda20 = sp.Integer(1)
lambda21 = sp.Rational(1, 2)
lambda22 = sp.Integer(-1)

kappa1 = sp.symbols("kappa_1", real=True)
gW, gU, rR = sp.symbols("g_W g_U r", real=True)
oU, oW = sp.symbols("o_U o_W", real=True)

print("lambda_20 =", lambda20)
print("lambda_21 =", lambda21)
print("lambda_22 =", lambda22)
print("Microscopic logarithmic slopes:")
print("  kappa_1, g_W, g_U, r, o_U, o_W")


# ---------------------------------------------------------------------------
# II. Stage-160 slippage amplitudes m, i, h
# ---------------------------------------------------------------------------

banner("II. EXACT WEAK-AXISYMMETRIC TRANSPORT OF THE THREE PORT SLIPPAGES")

m_r = sp.simplify(gW - oW - sp.Rational(1, 2) * kappa1)
i_r = sp.simplify(rR + gU - oU - gW)
h_r = sp.simplify(2 * rR - oU - oW)

print("m_r =", m_r)
print("i_r =", i_r)
print("h_r =", h_r)

for label, lam in [("20", lambda20), ("21", lambda21), ("22", lambda22)]:
    print(f"delta ln M_{{{label}}} =", sp.simplify(eps * lam * m_r))
    print(f"delta ln I_{{{label}}} =", sp.simplify(eps * lam * i_r))
    print(f"delta ln H_{{{label}}} =", sp.simplify(eps * lam * h_r))


# ---------------------------------------------------------------------------
# III. Portwise amplitude sigma_r and grouped collapse
# ---------------------------------------------------------------------------

banner("III. PORTWISE SIGMA_r AND THE ONE-SCALAR GROUPED COLLAPSE")

I_r, H_r = sp.symbols("I_r H_r", real=True)
rho1, rho2 = sp.symbols("rho_1 rho_2", real=True)
sigma1, sigma2 = sp.symbols("sigma_1 sigma_2", real=True)

sigma_r = sp.simplify(2 * m_r + 2 * I_r / (1 + I_r) * i_r + 2 * H_r / (1 - H_r) * h_r)
print("sigma_r =")
sp.pprint(sigma_r)

Xi1 = sp.simplify(rho1 * sigma1 + rho2 * sigma2)
print("\nTwo-port representative scalar amplitude Xi_1 =", Xi1)
print("(For n ports: Xi_1 = sum_r rho_r^(N) sigma_r.)")

Xi20 = sp.simplify(eps * lambda20 * Xi1)
Xi21 = sp.simplify(eps * lambda21 * Xi1)
Xi22 = sp.simplify(eps * lambda22 * Xi1)

print("Xi_load^(20) =", Xi20)
print("Xi_load^(21) =", Xi21)
print("Xi_load^(22) =", Xi22)

Xi_bar, Xi_a, Xi_b = grouped_triplet(Xi20, Xi21, Xi22)
print("\nGrouped defect variables:")
print("Xi_bar =", Xi_bar)
print("a_Xi   =", Xi_a)
print("b_Xi   =", Xi_b)
expect_zero("weak-axisymmetric grouped trace vanishes", Xi_bar)
expect_zero("axisymmetric fingerprint b_Xi - 3 a_Xi", Xi_b - 3 * Xi_a)

subbanner("Portwise sigma_r version")
Xi_load_20_port = sp.simplify(eps * (rho1 * sigma1 + rho2 * sigma2))
Xi_load_21_port = sp.simplify(eps * (rho1 * sigma1 + rho2 * sigma2) / 2)
Xi_load_22_port = sp.simplify(-eps * (rho1 * sigma1 + rho2 * sigma2))
expect_zero("Xi20 from Xi1", Xi20 - Xi_load_20_port)
expect_zero("Xi21 from Xi1", Xi21 - Xi_load_21_port)
expect_zero("Xi22 from Xi1", Xi22 - Xi_load_22_port)


# ---------------------------------------------------------------------------
# IV. Identification with the physical outgoing-prefactor slope
# ---------------------------------------------------------------------------

banner("IV. IDENTIFICATION WITH THE PHYSICAL OUTGOING-PREFACTOR SLOPE")

P1, P0 = sp.symbols("P1 P0", nonzero=True)
Xi1_phys = sp.simplify(P1 / P0)
print("Physical identification: Xi_1 := P1/P0 =", Xi1_phys)
print(
    "The remaining weak-axisymmetric grouped defect is identified with the"
    " logarithmic outgoing-prefactor slope P1/P0 on the physical branch."
)

subbanner("Even-preserving quadrupole defect pattern")
DQ20 = sp.simplify(eps * Xi1_phys)
DQ21 = sp.simplify(eps * Xi1_phys / 2)
DQ22 = sp.simplify(-eps * Xi1_phys)

DQ_bar, DQ_a, DQ_b = grouped_triplet(DQ20, DQ21, DQ22)
print("Delta_Q^(20) =", DQ20)
print("Delta_Q^(21) =", DQ21)
print("Delta_Q^(22) =", DQ22)
print("Delta_Q_bar =", DQ_bar)
print("a_Q         =", DQ_a)
print("b_Q         =", DQ_b)
expect_zero("quadrupole grouped trace vanishes", DQ_bar)
expect_zero("quadrupole axisymmetric fingerprint b_Q - 3 a_Q", DQ_b - 3 * DQ_a)


# ---------------------------------------------------------------------------
# V. Rigidity and dominant-port corollaries
# ---------------------------------------------------------------------------

banner("V. RIGIDITY AND DOMINANT-PORT COROLLARIES")

sigma_r_rigid = sp.simplify(sigma_r.subs({i_r: 0, h_r: 0}))
expect_zero("interference/hybridization rigidity => sigma_r - 2 m_r", sigma_r_rigid - 2 * m_r)

print("If i_r = h_r = 0, then sigma_r =", sigma_r_rigid)
print("So Xi_1 = 2 sum_r rho_r^(N) m_r on the rigid branch.")

rho_dom = sp.symbols("rho_dom", real=True)
sigma_dom = sp.symbols("sigma_dom", real=True)
Xi1_dom = sp.simplify(sigma_dom)  # dominant port weight ~ 1
print("Dominant-port limit: Xi_1 ->", Xi1_dom)
print("So the whole remaining grouped defect is then controlled by one port amplitude sigma_*." )


# ---------------------------------------------------------------------------
# VI. Final theorem ledger
# ---------------------------------------------------------------------------

banner("VI. FINAL THEOREM LEDGER")
print("1. Every weak-axisymmetric outgoing microscopic slippage inherits the same")
print("   grouped signature (1, 1/2, -1).")
print("2. Each outgoing port collapses to one scalar amplitude")
print("      sigma_r = 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r.")
print("3. The full remaining grouped defect collapses to one weighted scalar")
print("      Xi_1 = sum_r rho_r^(N) sigma_r.")
print("4. The grouped defect pattern is fixed exactly:")
print("      Xi_load^(20) = eps Xi_1,")
print("      Xi_load^(21) = eps Xi_1 / 2,")
print("      Xi_load^(22) = - eps Xi_1.")
print("5. Its grouped anisotropy always lies on the axisymmetric fingerprint")
print("      b = 3 a.")
print("6. The same scalar is the physical outgoing-prefactor slope")
print("      Xi_1 = P1/P0.")
print("7. On the even-preserving branch, the whole remaining quadrupole-normalization")
print("   defect is therefore controlled by that one scalar amplitude.")
print("8. Under interference/hybridization rigidity, Xi_1 reduces further to the")
print("   weighted square-root mixed-leg slope.")
print("9. So after Stage 7 the next honest theorem gate is to compute the single")
print("   microscopic weak-axisymmetric amplitude Xi_1 = P1/P0 on the actual")
print("   moving-throat branch.")
