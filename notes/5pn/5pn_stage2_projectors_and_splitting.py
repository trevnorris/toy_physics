#!/usr/bin/env python3
"""
5pn_stage2_projectors_and_splitting.py

Second executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Encodes the exact weighted grouped projectors for the {20,21,22} bundle.
2. Verifies projector idempotence, orthogonality, and reconstruction.
3. Derives the weak-axisymmetric splitting law b = 3 a.
4. Derives the first-order transport of conservative anisotropy into the
   normalized grouped response u2^(A).
5. Derives the first-order transport of conservative/transfer anisotropy into the
   outgoing static prefactor P0^(A).

This script is useful once actual PDE overlap data begin producing grouped-lane
triples, because it turns any non-isotropic output into precise diagnostics.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda x: sp.simplify(sp.expand(x)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("I. EXACT GROUPED PROJECTOR CALCULUS")

Ggrp = sp.diag(1, 2, 2)
Pbar = sp.Matrix([[1, 2, 2], [1, 2, 2], [1, 2, 2]]) / 5
Pa = sp.Matrix([[16, -8, -8], [-4, 2, 2], [-4, 2, 2]]) / 20
Pb = sp.Matrix([[0, 0, 0], [0, 2, -2], [0, -2, 2]]) / 4
I3 = sp.eye(3)

print("Ggrp =")
sp.pprint(Ggrp)
print("Pbar =")
sp.pprint(Pbar)
print("Pa =")
sp.pprint(Pa)
print("Pb =")
sp.pprint(Pb)

expect_zero("Pbar + Pa + Pb - I", Pbar + Pa + Pb - I3)
expect_zero("Pbar^2 - Pbar", Pbar * Pbar - Pbar)
expect_zero("Pa^2 - Pa", Pa * Pa - Pa)
expect_zero("Pb^2 - Pb", Pb * Pb - Pb)
expect_zero("Pbar*Pa", Pbar * Pa)
expect_zero("Pbar*Pb", Pbar * Pb)
expect_zero("Pa*Pb", Pa * Pb)

subx20, subx21, subx22 = sp.symbols("x20 x21 x22", real=True)
x = sp.Matrix([subx20, subx21, subx22])
xbar = sp.simplify((subx20 + 2 * subx21 + 2 * subx22) / 5)
a = sp.simplify((2 * subx20 - subx21 - subx22) / 10)
b = sp.simplify((subx21 - subx22) / 2)

xbar_part = sp.simplify(Pbar * x)
a_part = sp.simplify(Pa * x)
b_part = sp.simplify(Pb * x)
print("\nPbar x =")
sp.pprint(xbar_part)
print("Pa x =")
sp.pprint(a_part)
print("Pb x =")
sp.pprint(b_part)
expect_zero("reconstruction", xbar_part + a_part + b_part - x)
expect_zero("Pbar x - xbar*(1,1,1)", xbar_part - sp.Matrix([xbar, xbar, xbar]))
expect_zero("Pa x - a*(4,-1,-1)", a_part - sp.Matrix([4 * a, -a, -a]))
expect_zero("Pb x - b*(0,1,-1)", b_part - sp.Matrix([0, b, -b]))

banner("II. WEAK-AXISYMMETRIC SPLITTING LAW")

x0, x1, eps = sp.symbols("x0 x1 eps", real=True)
xa20 = x0 + eps * x1
xa21 = x0 + eps * x1 / 2
xa22 = x0 - eps * x1

xbar_ax = sp.simplify((xa20 + 2 * xa21 + 2 * xa22) / 5)
a_ax = sp.simplify((2 * xa20 - xa21 - xa22) / 10)
b_ax = sp.simplify((xa21 - xa22) / 2)
A2_ax = sp.simplify(4 * a_ax**2 + sp.Rational(4, 5) * b_ax**2)

print("xbar(axisymmetric) =", xbar_ax)
print("a(axisymmetric)    =", a_ax)
print("b(axisymmetric)    =", b_ax)
print("A^2(axisymmetric)  =", A2_ax)
expect_zero("axisymmetric xbar - x0", xbar_ax - x0)
expect_zero("axisymmetric b - 3 a", b_ax - 3 * a_ax)
expect_zero("axisymmetric A^2 - 7 eps^2 x1^2 / 10", A2_ax - sp.Rational(7, 10) * eps**2 * x1**2)

banner("III. FIRST-ORDER TRANSPORT INTO u2^(A)")

D0, u2_iso = sp.symbols("D0 u2")
dD0, dD2 = sp.symbols("dD0 dD2")

delta_u2 = sp.simplify(-(dD2 + u2_iso * dD0) / D0)
print("delta u2^(A) =", delta_u2)

ad0, bd0, ad2, bd2 = sp.symbols("aD0 bD0 aD2 bD2")
au2 = sp.simplify(-(ad2 + u2_iso * ad0) / D0)
bu2 = sp.simplify(-(bd2 + u2_iso * bd0) / D0)
print("a_(u,2) =", au2)
print("b_(u,2) =", bu2)

banner("IV. FIRST-ORDER TRANSPORT INTO P0^(A)")

N0, P0_iso = sp.symbols("N0 P0")
dN0 = sp.symbols("dN0")

delta_P0 = sp.simplify((dN0 - P0_iso * dD0) / D0)
print("delta P0^(A) =", delta_P0)

an0, bn0 = sp.symbols("aN0 bN0")
aP0 = sp.simplify((an0 - P0_iso * ad0) / D0)
bP0 = sp.simplify((bn0 - P0_iso * bd0) / D0)
print("a_(P,0) =", aP0)
print("b_(P,0) =", bP0)

banner("V. DIAGNOSTIC LEDGER")
print("1. Any exact O(3)-invariant grouped bundle must satisfy Pa x = Pb x = 0.")
print("2. A weak axisymmetric Y20 contamination forces b = 3 a.")
print("3. Conservative grouped anisotropy reaches u2^(A) only through")
print("      delta D_(A,2) + u2 * delta D_(A,0).")
print("4. Outgoing-prefactor anisotropy reaches P0^(A) only through")
print("      delta N_(A,0) - P0 * delta D_(A,0).")
print("5. So once PDE overlap data arrive, these formulas separate")
print("   support-side anisotropy from transfer-side anisotropy immediately.")
