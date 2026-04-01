#!/usr/bin/env python3
"""
4pn_quartic_legendre_audit.py

Exact quartic-order perturbative Legendre compiler for the 4PN lift.

What this script does
---------------------
1. Considers a generic perturbative ordinary Lagrangian
       L = L0 + eps L1 + eps^2 L2 + eps^3 L3 + eps^4 L4
   with quadratic L0.
2. Solves p = dL/dv perturbatively through O(eps^3).
3. Verifies the exact quartic Hamiltonian coefficient formula

   H1 = -L1(v0)
   H2 = -L2(v0) + 1/2 A0 M^{-1} A0
   H3 = -L3(v0) + A0 M^{-1} B0 - 1/2 A0 M^{-1} C0 M^{-1} A0
   H4 = -L4(v0)
        + A0 M^{-1} D0
        + 1/2 B0 M^{-1} B0
        - B0 M^{-1} C0 M^{-1} A0
        - 1/2 A0 M^{-1} E0 M^{-1} A0
        + 1/2 A0 M^{-1} C0 M^{-1} C0 M^{-1} A0
        + 1/6 T0[M^{-1}A0, M^{-1}A0, M^{-1}A0]

   where
       v0 = M^{-1} p,
       A0 = (dL1/dv)|_{v0},
       B0 = (dL2/dv)|_{v0},
       D0 = (dL3/dv)|_{v0},
       C0 = (d^2 L1/dv^2)|_{v0},
       E0 = (d^2 L2/dv^2)|_{v0},
       T0 = (d^3 L1/dv^3)|_{v0}.

This is the exact compiler identity needed before any trustworthy 4PN comparable-mass
solve can begin.
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


# ---------------------------------------------------------------------------
# Generic 1D test model (enough to verify the exact quartic identity)
# ---------------------------------------------------------------------------

banner("QUARTIC-ORDER PERTURBATIVE LEGENDRE COMPILER")

eps = sp.symbols("eps", real=True)
m, p, v = sp.symbols("m p v", positive=True, real=True)

# Generic enough polynomial data to excite every quartic compiler term.
a1, a2, a3 = sp.symbols("a1 a2 a3", real=True)
b1, b2 = sp.symbols("b1 b2", real=True)
c1, c2 = sp.symbols("c1 c2", real=True)
d1, d2 = sp.symbols("d1 d2", real=True)

L0 = sp.Rational(1, 2) * m * v**2
L1 = a1 * v + a2 * v**2 + a3 * v**3
L2 = b1 * v + b2 * v**2
L3 = c1 * v + c2 * v**2
L4 = d1 * v + d2 * v**2
L = L0 + eps * L1 + eps**2 * L2 + eps**3 * L3 + eps**4 * L4

v0 = p / m
w1, w2, w3 = sp.symbols("w1 w2 w3", real=True)
v_series = v0 + eps * w1 + eps**2 * w2 + eps**3 * w3

peq = sp.expand(sp.diff(L, v).subs(v, v_series) - p)
eq1 = sp.expand(peq.coeff(eps, 1))
eq2 = sp.expand(peq.coeff(eps, 2))
eq3 = sp.expand(peq.coeff(eps, 3))

sol_w1 = sp.solve(sp.Eq(eq1, 0), w1)[0]
sol_w2 = sp.solve(sp.Eq(eq2.subs(w1, sol_w1), 0), w2)[0]
sol_w3 = sp.solve(sp.Eq(eq3.subs({w1: sol_w1, w2: sol_w2}), 0), w3)[0]

print("v1 =", sp.simplify(sol_w1))
print("v2 =", sp.simplify(sol_w2))
print("v3 =", sp.simplify(sol_w3))

# Exact Hamiltonian coefficients from the direct perturbative solve.
H_exact = p * (v0 + eps * sol_w1 + eps**2 * sol_w2 + eps**3 * sol_w3) - L.subs(
    v, v0 + eps * sol_w1 + eps**2 * sol_w2 + eps**3 * sol_w3
)
H_series = sp.expand(sp.series(H_exact, eps, 0, 5).removeO())

H1_exact = sp.expand(H_series.coeff(eps, 1))
H2_exact = sp.expand(H_series.coeff(eps, 2))
H3_exact = sp.expand(H_series.coeff(eps, 3))
H4_exact = sp.expand(H_series.coeff(eps, 4))

A0 = sp.diff(L1, v).subs(v, v0)
B0 = sp.diff(L2, v).subs(v, v0)
D0 = sp.diff(L3, v).subs(v, v0)
C0 = sp.diff(L1, v, 2).subs(v, v0)
E0 = sp.diff(L2, v, 2).subs(v, v0)
T0 = sp.diff(L1, v, 3).subs(v, v0)

H1_formula = sp.expand(-L1.subs(v, v0))
H2_formula = sp.expand(-L2.subs(v, v0) + sp.Rational(1, 2) * A0**2 / m)
H3_formula = sp.expand(-L3.subs(v, v0) + A0 * B0 / m - sp.Rational(1, 2) * A0**2 * C0 / m**2)
H4_formula = sp.expand(
    -L4.subs(v, v0)
    + A0 * D0 / m
    + sp.Rational(1, 2) * B0**2 / m
    - B0 * C0 * A0 / m**2
    - sp.Rational(1, 2) * A0**2 * E0 / m**2
    + sp.Rational(1, 2) * A0**2 * C0**2 / m**3
    + sp.Rational(1, 6) * A0**3 * T0 / m**3
)

print("\nExact coefficients from direct perturbative solve:")
print("H1 =", H1_exact)
print("H2 =", H2_exact)
print("H3 =", H3_exact)
print("H4 =", H4_exact)

print("\nClosed quartic compiler formulas:")
print("H1 = -L1(v0)")
print("H2 = -L2(v0) + 1/2 A0 M^{-1} A0")
print("H3 = -L3(v0) + A0 M^{-1} B0 - 1/2 A0 M^{-1} C0 M^{-1} A0")
print("H4 = -L4(v0) + A0 M^{-1} D0 + 1/2 B0 M^{-1} B0")
print("     - B0 M^{-1} C0 M^{-1} A0 - 1/2 A0 M^{-1} E0 M^{-1} A0")
print("     + 1/2 A0 M^{-1} C0 M^{-1} C0 M^{-1} A0")
print("     + 1/6 T0[M^{-1}A0, M^{-1}A0, M^{-1}A0]")

expect_zero("H1 exact - formula", H1_exact - H1_formula)
expect_zero("H2 exact - formula", H2_exact - H2_formula)
expect_zero("H3 exact - formula", H3_exact - H3_formula)
expect_zero("H4 exact - formula", H4_exact - H4_formula)

banner("FINAL QUARTIC COMPILER LEDGER")
print("The quartic perturbative Legendre compiler is exact.")
print("For a quadratic L0 with constant mass matrix M and v0 = M^{-1} p,")
print("the 4PN comparable-mass lift is controlled entirely by:")
print("  A0 = (∂L1/∂v)|_{v0}")
print("  B0 = (∂L2/∂v)|_{v0}")
print("  D0 = (∂L3/∂v)|_{v0}")
print("  C0 = (∂²L1/∂v²)|_{v0}")
print("  E0 = (∂²L2/∂v²)|_{v0}")
print("  T0 = (∂³L1/∂v³)|_{v0}")
print()
print("So before any 4PN residual solve is trusted, this compiler must be frozen and used")
print("as the unique ordinary-Lagrangian -> Hamiltonian bridge at quartic order.")
