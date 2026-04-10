#!/usr/bin/env python3
"""
3pn_sigma_collapse_and_unique_geometry_completion_audit.py

Audit the next exact 3PN step after the apparent P0/geometry one-parameter family.

Main result
-----------
The sigma-family introduced at COM level is algebraically redundant in the
full generic-frame mass polynomial.  The key identity

    nu * (p**3 + q**3) = (1 - 3*nu) * (p**2*q + p*q**2)

collapses the whole family to a unique pure-geometry static remainder.  The
result exactly matches the difference between the imported generic-frame static
3PN target coefficient and the static companion already forced by the richer
clustered-P2 middle-block fit.
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


banner("PART I — EXACT MASS-POLYNOMIAL IDENTITY")
p, q = sp.symbols("p q", positive=True)
M = p + q
mu = sp.simplify(p * q / M)
nu = sp.simplify(p * q / M**2)
U0 = p**3 + q**3
Ug = p**2 * q + p * q**2

identity = sp.simplify(nu * U0 - (1 - 3 * nu) * Ug)
expect_zero("nu*U0 - (1-3 nu)*Ug", identity)
print("nu =", nu)
print("U0 =", U0)
print("Ug =", Ug)


banner("PART II — SIGMA-FAMILY COLLAPSE")
sigma = sp.symbols("sigma", real=True)
f = sp.simplify((408 * nu**2 + 1232 * nu - 2080 + 63 * sp.pi**2) / 96)
expr_sigma = sp.simplify(sigma * nu * U0 + (f - sigma * (1 - 3 * nu)) * Ug)
expr_geom = sp.simplify(f * Ug)
expect_zero("sigma-family total - pure geometry form", expr_sigma - expr_geom)
print("collapsed coefficient f(nu) =", f)
print("unique generic-frame static polynomial =", sp.factor(expr_geom))

# Show that the sigma family is genuinely sigma-independent.
expr_dsigma = sp.diff(expr_sigma, sigma)
expect_zero("d/dsigma(total static polynomial)", expr_dsigma)


banner("PART III — DIRECT IMPORTED TARGET DECOMPOSITION")
# Direct imported generic-frame static residual coordinate from the earlier note.
c_target = sp.simplify(-sp.Rational(227, 24) + sp.Rational(21, 32) * sp.pi**2)
# Static companion already forced by the richer grouped-P2 middle block.
c_pred = sp.simplify((293 - 308 * nu - 102 * nu**2) / 24)
# Remaining unique geometry coefficient.
c_geom = sp.simplify(c_target - c_pred)

print("c_target =", c_target)
print("c_pred(P2 middle companion) =", c_pred)
print("c_geom(target - pred) =", c_geom)
expect_zero("c_geom - f(nu)", c_geom - f)
expect_zero("c_target - (c_pred + c_geom)", c_target - (c_pred + c_geom))

# Static generic-frame target split.
L_target = sp.simplify(c_target * Ug)
L_pred = sp.simplify(c_pred * Ug)
L_gap = sp.simplify(c_geom * Ug)
expect_zero("L_target - (L_pred + L_gap)", L_target - (L_pred + L_gap))
print("target static polynomial =", sp.factor(L_target))
print("grouped-P2 predicted static companion =", sp.factor(L_pred))
print("unique geometry remainder =", sp.factor(L_gap))


banner("PART IV — CONSISTENCY WITH THE COM REMAINDER")
G, r = sp.symbols("G r", positive=True)
U = sp.symbols("U", positive=True)
# COM reduction of the generic-frame basis element Ug.
com_map = sp.simplify((1 / mu) * (G**4 * p * q / r**4) * Ug)
com_target = sp.simplify(nu * U**4)
print("(1/mu) * (G^4 p q / r^4) * Ug =", com_map)
print("expected COM image = nu * U^4")
# After using U = G M / r, the equality is exact.
expect_zero(
    "COM map with U=GM/r",
    com_map.subs(U, G * M / r) - com_target.subs(U, G * M / r),
)

l15_gap = sp.simplify(nu * c_geom)
expected_l15_gap = sp.simplify(nu * (408 * nu**2 + 1232 * nu - 2080 + 63 * sp.pi**2) / 96)
expect_zero("l15_gap consistency", l15_gap - expected_l15_gap)
print("Delta l15^(g) =", l15_gap)


banner("PART V — HAMILTONIAN LIFT")
H_gap = -L_gap
print("ordinary static remainder =", sp.factor(L_gap))
print("Hamiltonian static remainder =", sp.factor(H_gap))
expect_zero("Hamiltonian sign-flip identity", H_gap + L_gap)


banner("PART VI — FINAL LEDGER")
print("1. The apparent COM one-parameter P0/g family is algebraically redundant.")
print("2. The exact identity nu*(p^3+q^3) = (1-3 nu)*(p^2 q + p q^2) collapses the")
print("   whole family to a unique pure-geometry static remainder.")
print("3. The direct imported generic-frame static target coefficient splits exactly into")
print("   the grouped-P2-predicted static companion plus this unique geometry remainder.")
print("4. Therefore sigma is not a physical ambiguity in the fixed ADM chart; it is only")
print("   a COM repartition of the same generic mass polynomial.")
print("5. The remaining conservative 3PN bottleneck is no longer the scalar split but the")
print("   separate pure-kinetic slot Delta l1.")
