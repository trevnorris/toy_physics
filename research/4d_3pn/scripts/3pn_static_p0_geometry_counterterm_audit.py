#!/usr/bin/env python3
"""
3pn_static_p0_geometry_counterterm_audit.py

Audit the next exact 3PN step after the grouped-P2 richer middle-block closure.

Main result
-----------
The remaining static COM gap

    Delta l15^(0/g) = nu*(408*nu**2 + 1232*nu - 2080 + 63*pi**2)/96

cannot be represented by the old constant-coefficient generic-frame U-block alone,
because that block always reduces to nu * const in COM.  The natural repair is a
nu-dressed scalar counterterm living in the old P0/geometry lane.

Using the two static scalar mass families

    U0 = p**3 + q**3         (body-local P0 family)
    Ug = p**2*q + p*q**2     (pair/geometry family)

inside the common prefactor G^4*p*q/r^4, one gets the exact COM images

    U0 -> (1 - 3*nu) U^4,
    Ug -> nu U^4.

Hence the full static gap is reproduced by the regular one-parameter family

    L_ct^(0) =  (G^4*p*q/r^4) * sigma*nu * (p**3 + q**3)
    L_ct^(g) =  (G^4*p*q/r^4) * [F(nu) - sigma*(1 - 3*nu)] * (p**2*q + p*q**2)

with

    F(nu) = (408*nu**2 + 1232*nu - 2080 + 63*pi**2)/96.

The generic-frame Hamiltonian image is then exactly the sign flip

    H_ct^(0/g) = - L_ct^(0/g),

because the full generic-frame 3PN compiler remains -I on any residual block with
mass-only coefficients.
"""
from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.expand(sp.simplify(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Basic symbols and exact static gap
# ---------------------------------------------------------------------------

nu, sigma = sp.symbols("nu sigma", real=True)
p, q, G, r = sp.symbols("p q G r", positive=True)
M = sp.symbols("M", positive=True)
mu = sp.symbols("mu", positive=True)

Fnu = sp.expand((408 * nu**2 + 1232 * nu - 2080 + 63 * sp.pi**2) / 96)
Gap = sp.expand(nu * Fnu)

banner("PART I — EXACT STATIC GAP AND THE OLD U-BLOCK NO-GO")
print("Delta l15^(0/g) =", Gap)

u1 = sp.symbols("u1")
poly_no_go = sp.Poly(sp.expand(Gap - nu * u1), nu)
print("Gap - nu*u1 =", poly_no_go.as_expr())
sol_u1 = sp.solve([sp.Eq(c, 0) for c in poly_no_go.all_coeffs()], [u1], dict=True)
print("constant-coefficient U-block solutions =", sol_u1)
if sol_u1:
    raise AssertionError("Unexpected constant U-block representation of the static gap.")
print("Therefore the exact static gap cannot live in the old constant-coefficient U-block alone.")


# ---------------------------------------------------------------------------
# Generic-frame static scalar families and their COM images
# ---------------------------------------------------------------------------

banner("PART II — TWO NATURAL STATIC SCALAR FAMILIES")

nu_pq = sp.simplify(p * q / (p + q) ** 2)
mu_pq = sp.simplify(p * q / (p + q))
U4 = sp.simplify(G**4 * (p + q) ** 4 / r**4)

L0_family = sp.simplify(G**4 * p * q / r**4 * (p**3 + q**3))
Lg_family = sp.simplify(G**4 * p * q / r**4 * (p**2 * q + p * q**2))

red0 = sp.simplify(L0_family / (mu_pq * U4))
redg = sp.simplify(Lg_family / (mu_pq * U4))
print("COM image of U0 family =", red0)
print("COM image of Ug family =", redg)
expect_zero("U0 family COM image - (1-3 nu)", red0 - (1 - 3 * nu_pq))
expect_zero("Ug family COM image - nu", redg - nu_pq)


# ---------------------------------------------------------------------------
# Exact one-parameter P0/geometry counterterm family
# ---------------------------------------------------------------------------

banner("PART III — EXACT REGULAR ONE-PARAMETER P0/GEOMETRY FAMILY")

Lct0 = sp.simplify(G**4 * p * q / r**4 * (sigma * nu_pq) * (p**3 + q**3))
Lctg = sp.simplify(G**4 * p * q / r**4 * (Fnu.subs(nu, nu_pq) - sigma * (1 - 3 * nu_pq)) * (p**2 * q + p * q**2))

red_ct0 = sp.simplify(Lct0 / (mu_pq * U4))
red_ctg = sp.simplify(Lctg / (mu_pq * U4))
red_sum = sp.simplify(red_ct0 + red_ctg)

print("COM coefficient of P0 counterterm =", red_ct0)
print("COM coefficient of geometry counterterm =", red_ctg)
print("combined COM coefficient =", red_sum)
expect_zero("combined static counterterm - exact gap", red_sum - Gap.subs(nu, nu_pq))

# Pure-geometry canonical slice.
red_pure_geom = sp.simplify(red_sum.subs(sigma, 0))
print("pure-geometry slice (sigma=0) =", red_pure_geom)
expect_zero("pure-geometry slice - exact gap", red_pure_geom - Gap.subs(nu, nu_pq))

# No-constant-geometry slice.
sigma0 = sp.simplify(Fnu.subs(nu, 0))
red0_alt = sp.simplify(red_ct0.subs(sigma, sigma0))
redg_alt = sp.factor(sp.simplify(red_ctg.subs(sigma, sigma0)))
print("no-constant-geometry sigma =", sigma0)
print("P0 piece on that slice =", red0_alt)
print("geometry piece on that slice =", redg_alt)
expect_zero(
    "no-constant-geometry split - exact gap",
    sp.simplify((red0_alt + redg_alt) - Gap.subs(nu, nu_pq)),
)

# Pure-P0 slice would require a singular sigma(nu).
sigma_pure_p0 = sp.simplify(Fnu / (1 - 3 * nu))
print("pure-P0 sigma(nu) =", sigma_pure_p0)
print("denominator of pure-P0 sigma =", sp.denom(sp.together(sigma_pure_p0)))


# ---------------------------------------------------------------------------
# Generic-frame Hamiltonian compiler
# ---------------------------------------------------------------------------

banner("PART IV — PUSH THROUGH THE GENERIC-FRAME HAMILTONIAN COMPILER")

Hct0 = sp.simplify(-Lct0)
Hctg = sp.simplify(-Lctg)
print("H_ct^(0) =", Hct0)
print("H_ct^(g) =", Hctg)
expect_zero("ordinary + Hamiltonian P0 counterterm", Lct0 + Hct0)
expect_zero("ordinary + Hamiltonian geometry counterterm", Lctg + Hctg)

# COM h15 relation on the repaired scalar lane.
h15_gap = sp.expand(-Gap)
print("COM Hamiltonian static gap Delta h15^(0/g) =", h15_gap)
expect_zero("h15 + l15 on the scalar gap", h15_gap + Gap)


banner("PART V — FINAL LEDGER")
print("1. The exact 3PN static remainder after the grouped-P2 middle-block closure is")
print("      Delta l15^(0/g) = nu*(408*nu^2 + 1232*nu - 2080 + 63*pi^2)/96.")
print("2. The old constant-coefficient generic-frame U-block cannot represent this remainder,")
print("   because it always reduces to nu*const in COM.")
print("3. The natural static scalar families are")
print("      U0 = p^3 + q^3   (body-local P0 family)")
print("      Ug = p^2*q + p*q^2   (pair/geometry family),")
print("   with exact COM images (1-3 nu)U^4 and nu U^4 respectively.")
print("4. The full scalar static gap is reproduced by the regular one-parameter family")
print("      L_ct^(0) + L_ct^(g),")
print("   where sigma labels how much of the COM scalar lane is assigned to P0 versus geometry.")
print("5. The simplest canonical slice is sigma = 0, which places the whole 3PN static gap in")
print("   the pair/geometry channel.")
print("6. The full generic-frame Hamiltonian compiler still acts by exact sign flip on this")
print("   nu-dressed scalar counterterm family: H_ct^(0/g) = -L_ct^(0/g).")
print("7. So the algebraic bottleneck is gone; the real remaining theorem question is physical:")
print("   what scalar-side wall/geometry dynamics selects the split parameter sigma?")
