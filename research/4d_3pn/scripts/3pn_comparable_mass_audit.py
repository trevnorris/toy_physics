#!/usr/bin/env python3
"""
3pn_comparable_mass_audit.py

First working scaffold for the 3PN comparable-mass lift.

What this script does
---------------------
1. Derives and verifies the cubic-order perturbative Legendre transform formula
   needed for a 3PN generic-frame residual solve.
2. Records the natural 3PN self/static two-body seed obtained by symmetric
   promotion of the exact one-body Schwarzschild target.
3. Generates an intentionally overcomplete exchange-symmetric residual basis for
   the new comparable-mass 3PN blocks:
      - G/r   times sextic velocity invariants,
      - G^2/r^2 times quartic velocity invariants with one extra mass power,
      - G^3/r^3 times quadratic velocity invariants with two extra mass powers,
      - G^4/r^4 static cross mass polynomial.
4. Verifies that every residual-basis element vanishes in the strict test-mass
   limit (body B fixed, body A test mass), so the one-body gate is cleanly
   separated from the comparable-mass lift.

Interpretation
--------------
This is not yet the full 3PN solve. It is the exact scaffold needed before that
solve is worth attempting:
  - a cubic Legendre-transform identity,
  - a frozen self/static seed,
  - and a concrete residual basis.
"""

from __future__ import annotations

import math
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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Part I. Cubic-order perturbative Legendre transform
# ---------------------------------------------------------------------------

def legendre_formula_check() -> None:
    banner("PART I — CUBIC-ORDER PERTURBATIVE LEGENDRE TRANSFORM")

    eps = sp.symbols("eps", real=True)
    m, p = sp.symbols("m p", positive=True, real=True)
    a1, a2, b1, b2, c1, c2 = sp.symbols("a1 a2 b1 b2 c1 c2", real=True)
    v = sp.symbols("v", real=True)

    # A deliberately generic 1D test model. 1D is enough to verify the exact
    # cubic-order identity because the formula is tensorial and the symbolic
    # cancellations already appear at this level.
    L0 = sp.Rational(1, 2) * m * v**2
    L1 = a1 * v**3 + a2 * v**4
    L2 = b1 * v**2 + b2 * v**3
    L3 = c1 * v + c2 * v**2
    L = L0 + eps * L1 + eps**2 * L2 + eps**3 * L3

    v0 = p / m
    w1, w2 = sp.symbols("w1 w2", real=True)
    v_series = v0 + eps * w1 + eps**2 * w2

    # Solve p = dL/dv perturbatively through O(eps^2), which is enough for H3.
    p_eq = sp.expand(sp.diff(L, v).subs(v, v_series) - p)
    eq1 = sp.expand(p_eq.coeff(eps, 1))
    eq2 = sp.expand(p_eq.coeff(eps, 2))

    sol_w1 = sp.solve(sp.Eq(eq1, 0), w1)[0]
    sol_w2 = sp.solve(sp.Eq(eq2.subs(w1, sol_w1), 0), w2)[0]

    print("v1 =", sp.simplify(sol_w1))
    print("v2 =", sp.simplify(sol_w2))

    H_exact = sp.expand(
        p * (v0 + eps * sol_w1 + eps**2 * sol_w2)
        - L.subs(v, v0 + eps * sol_w1 + eps**2 * sol_w2)
    )
    H_series = sp.expand(sp.series(H_exact, eps, 0, 4).removeO())

    H0_exact = sp.expand(H_series.coeff(eps, 0))
    H1_exact = sp.expand(H_series.coeff(eps, 1))
    H2_exact = sp.expand(H_series.coeff(eps, 2))
    H3_exact = sp.expand(H_series.coeff(eps, 3))

    A0 = sp.diff(L1, v).subs(v, v0)
    B0 = sp.diff(L2, v).subs(v, v0)
    C0 = sp.diff(L1, v, 2).subs(v, v0)

    H1_formula = sp.expand(-L1.subs(v, v0))
    H2_formula = sp.expand(-L2.subs(v, v0) + sp.Rational(1, 2) * A0**2 / m)
    H3_formula = sp.expand(-L3.subs(v, v0) + A0 * B0 / m - sp.Rational(1, 2) * A0**2 * C0 / m**2)

    print("\nExact cubic Legendre coefficients from direct solve:")
    print("H0 =", H0_exact)
    print("H1 =", H1_exact)
    print("H2 =", H2_exact)
    print("H3 =", H3_exact)

    print("\nClosed formulas:")
    print("H1 = -L1(v0)")
    print("H2 = -L2(v0) + 1/2 A0 M^{-1} A0")
    print("H3 = -L3(v0) + A0 M^{-1} B0 - 1/2 A0 M^{-1} C0 M^{-1} A0")

    expect_zero("H1 exact - formula", H1_exact - H1_formula)
    expect_zero("H2 exact - formula", H2_exact - H2_formula)
    expect_zero("H3 exact - formula", H3_exact - H3_formula)

    print("\nVector/tensor form carried forward to the two-body 3PN lift:")
    print("  v0   = M^{-1} p")
    print("  A0   = (∂L1/∂v)|_{v0}")
    print("  B0   = (∂L2/∂v)|_{v0}")
    print("  C0   = (∂²L1/∂v²)|_{v0}")
    print("  H1   = -L1(v0)")
    print("  H2   = -L2(v0) + 1/2 A0^T M^{-1} A0")
    print("  H3   = -L3(v0) + A0^T M^{-1} B0 - 1/2 A0^T M^{-1} C0 M^{-1} A0")


# ---------------------------------------------------------------------------
# Part II. Natural 3PN self/static seed
# ---------------------------------------------------------------------------

def self_static_seed() -> None:
    banner("PART II — NATURAL 3PN SELF/STATIC SEED")

    G, c, r = sp.symbols("G c r", positive=True, real=True)
    mA, mB = sp.symbols("mA mB", positive=True, real=True)
    vA2, vB2 = sp.symbols("vA2 vB2", nonnegative=True, real=True)

    # The exact one-body target slots promoted symmetrically to the two-body seed.
    L3_seed = (
        sp.Rational(5, 128) * (mA * vA2**4 + mB * vB2**4)
        + sp.Rational(11, 16) * G * mA * mB / r * (vA2**3 + vB2**3)
        + sp.Rational(47, 16) * G**2 * mA * mB / r**2 * (mB * vA2**2 + mA * vB2**2)
        + sp.Rational(13, 8) * G**3 * mA * mB / r**3 * (mB**2 * vA2 + mA**2 * vB2)
        - sp.Rational(1, 8) * G**4 * mA * mB / r**4 * (mB**3 + mA**3)
    )

    print("L3_seed (this multiplies c^{-6} in the full conservative ledger) =")
    sp.pprint(L3_seed)

    print("\nInterpretation:")
    print("  - v^8 slot     -> free 3PN kinematics")
    print("  - G/r block    -> exact one-body U v^6 seed")
    print("  - G^2/r^2 block-> exact one-body U^2 v^4 seed")
    print("  - G^3/r^3 block-> exact one-body U^3 v^2 seed")
    print("  - G^4/r^4 block-> exact one-body U^4 static seed")
    print("\nThis is the natural 3PN analogue of the frozen 2PN self/static seed.")


# ---------------------------------------------------------------------------
# Part III. Overcomplete comparable-mass residual basis
# ---------------------------------------------------------------------------

def swap_full(expr: sp.Expr, a: sp.Symbol, b: sp.Symbol, c: sp.Symbol, d: sp.Symbol, e: sp.Symbol, p: sp.Symbol, q: sp.Symbol) -> sp.Expr:
    return sp.expand(expr.xreplace({a: b, b: a, c: c, d: e, e: d, p: q, q: p}))


def canonical_sym(expr: sp.Expr, vars_order: tuple[sp.Symbol, ...], a: sp.Symbol, b: sp.Symbol, c: sp.Symbol, d: sp.Symbol, e: sp.Symbol, p: sp.Symbol, q: sp.Symbol) -> sp.Expr:
    s = sp.expand(expr + swap_full(expr, a, b, c, d, e, p, q))
    poly = sp.Poly(s, *vars_order, domain="QQ")
    coeffs = poly.coeffs()

    den_lcm = 1
    for coef in coeffs:
        frac = sp.Rational(coef)
        den_lcm = sp.ilcm(den_lcm, frac.q)
    ints = [int(sp.Rational(coef) * den_lcm) for coef in coeffs]

    g = 0
    for n in ints:
        g = math.gcd(g, abs(n))
    if g == 0:
        g = 1

    s_norm = sp.expand(s * den_lcm / g)
    poly2 = sp.Poly(s_norm, *vars_order)
    terms = poly2.terms()
    if terms and terms[0][1] < 0:
        s_norm = -s_norm
    return sp.expand(s_norm)


def generate_basis(mass_deg: int, vel_deg: int) -> list[sp.Expr]:
    a, b, c, d, e, p, q = sp.symbols("a b c d e p q")
    vars_order = (p, q, a, b, c, d, e)
    basis: set[sp.Expr] = set()

    for mp in range(mass_deg + 1):
        mq = mass_deg - mp
        for pa in range(5):
            for pb in range(5):
                for pc in range(5):
                    for pd in range(vel_deg + 1):
                        for pe in range(vel_deg + 1):
                            this_deg = 2 * pa + 2 * pb + 2 * pc + pd + pe
                            if this_deg != vel_deg:
                                continue
                            expr = p**mp * q**mq * a**pa * b**pb * c**pc * d**pd * e**pe
                            sym = canonical_sym(expr, vars_order, a, b, c, d, e, p, q)

                            # Strict one-body limit used here:
                            #   body A is test mass (p -> 0),
                            #   body B is the fixed central source (q -> 1, b = c = e = 0).
                            # Any genuine comparable-mass residual must vanish on this branch.
                            tm = sp.expand(sym.subs({b: 0, c: 0, e: 0, p: 0, q: 1}))
                            if tm != 0:
                                continue
                            basis.add(sym)

    return sorted(basis, key=str)


def comparable_mass_basis() -> None:
    banner("PART III — OVERCOMPLETE COMPARABLE-MASS RESIDUAL BASIS")

    G1_basis = generate_basis(mass_deg=0, vel_deg=6)
    G2_basis = generate_basis(mass_deg=1, vel_deg=4)
    G3_basis = generate_basis(mass_deg=2, vel_deg=2)
    G4_basis = generate_basis(mass_deg=3, vel_deg=0)

    print("Basis counts before any contact/gauge reduction:")
    print(f"  G/r sextic block      : {len(G1_basis)}")
    print(f"  G^2/r^2 quartic block : {len(G2_basis)}")
    print(f"  G^3/r^3 quadratic block: {len(G3_basis)}")
    print(f"  G^4/r^4 static block  : {len(G4_basis)}")

    subbanner("III.1 — G/r sextic residual invariants")
    for i, term in enumerate(G1_basis, start=1):
        print(f"Q{i:02d} = {term}")

    subbanner("III.2 — G^2/r^2 quartic mass-weighted residual invariants")
    for i, term in enumerate(G2_basis, start=1):
        print(f"T{i:02d} = {term}")

    subbanner("III.3 — G^3/r^3 quadratic mass-weighted residual invariants")
    for i, term in enumerate(G3_basis, start=1):
        print(f"S{i:02d} = {term}")

    subbanner("III.4 — G^4/r^4 static cross polynomial")
    for i, term in enumerate(G4_basis, start=1):
        print(f"U{i:02d} = {term}")

    print("\nInterpretation:")
    print("  - This basis is intentionally overcomplete.")
    print("  - It is the clean starting point before any contact-transformation or")
    print("    gauge reduction is imposed at 3PN.")
    print("  - Every element vanishes in the strict test-mass limit, so the one-body")
    print("    exact gate is cleanly separated from the comparable-mass residual.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    legendre_formula_check()
    self_static_seed()
    comparable_mass_basis()
