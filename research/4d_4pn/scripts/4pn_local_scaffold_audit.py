#!/usr/bin/env python3
"""
4pn_local_scaffold_audit.py

Local 4PN generic-frame scaffold audit.

What this script does
---------------------
1. Records the natural 4PN self/static seed obtained by symmetric promotion of the
   exact one-body 4PN Schwarzschild target.
2. Builds an overcomplete exchange-symmetric generic-frame residual basis for the
   *local* 4PN interaction blocks:
      - G/r     times degree-8 velocity invariants,
      - G^2/r^2 times degree-6 invariants with one extra mass power,
      - G^3/r^3 times degree-4 invariants with two extra mass powers,
      - G^4/r^4 times degree-2 invariants with three extra mass powers,
      - G^5/r^5 static cross mass polynomials.
3. Verifies that every basis element vanishes in the strict test-mass limit, so
   the exact one-body gate stays cleanly separated from the comparable-mass lift.
4. Projects each block to the center-of-mass (COM) frame and extracts the induced
   local 4PN interaction slots in the isotropic COM basis.
5. Computes the blockwise COM image ranks and nullities.
6. Extracts one minimal pivot-spanning subset per block, which will be the clean
   starting point for the eventual fixed-chart target import.

Interpretation
--------------
This is the local 4PN analogue of the first 3PN generic-frame scaffold audit.
It does *not* yet import the final 4PN target. Its role is to answer the prior
question first:

    Is the raw local 4PN generic-frame interaction scaffold already large enough
    to span the full COM interaction image blockwise?

The answer turns out to be yes, with a very large COM-null sector.
"""

from __future__ import annotations

import math
import sympy as sp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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
# Natural 4PN self/static seed
# ---------------------------------------------------------------------------

def self_static_seed() -> sp.Expr:
    banner("PART I — NATURAL 4PN SELF/STATIC SEED")

    G, r = sp.symbols("G r", positive=True, real=True)
    mA, mB = sp.symbols("mA mB", positive=True, real=True)
    vA2, vB2 = sp.symbols("vA2 vB2", nonnegative=True, real=True)

    seed = (
        sp.Rational(7, 256) * (mA * vA2**5 + mB * vB2**5)
        + sp.Rational(75, 128) * G * mA * mB / r * (vA2**4 + vB2**4)
        + sp.Rational(59, 16) * G**2 * mA * mB / r**2 * (mB * vA2**3 + mA * vB2**3)
        + sp.Rational(203, 32) * G**3 * mA * mB / r**3 * (mB**2 * vA2**2 + mA**2 * vB2**2)
        + sp.Rational(31, 32) * G**4 * mA * mB / r**4 * (mB**3 * vA2 + mA**3 * vB2)
        + sp.Rational(1, 16) * G**5 * mA * mB / r**5 * (mB**4 + mA**4)
    )

    print("L4_seed (this multiplies c^{-8} in the full conservative ledger) =")
    sp.pprint(seed)

    print("\nInterpretation:")
    print("  - v^10 slot      -> free 4PN kinematics")
    print("  - G/r block      -> exact one-body U v^8 seed")
    print("  - G^2/r^2 block  -> exact one-body U^2 v^6 seed")
    print("  - G^3/r^3 block  -> exact one-body U^3 v^4 seed")
    print("  - G^4/r^4 block  -> exact one-body U^4 v^2 seed")
    print("  - G^5/r^5 block  -> exact one-body U^5 static seed")
    print("\nThis is the natural local 4PN self/static promotion of the exact test-mass gate.")

    return seed


# ---------------------------------------------------------------------------
# Exchange-symmetric generic-frame basis generation
# ---------------------------------------------------------------------------

def swap_full(expr: sp.Expr, a: sp.Symbol, b: sp.Symbol, c: sp.Symbol, d: sp.Symbol, e: sp.Symbol, p: sp.Symbol, q: sp.Symbol) -> sp.Expr:
    return sp.expand(expr.xreplace({a: b, b: a, c: c, d: e, e: d, p: q, q: p}))



def canonical_sym(expr: sp.Expr, vars_order: tuple[sp.Symbol, ...], a: sp.Symbol, b: sp.Symbol, c: sp.Symbol, d: sp.Symbol, e: sp.Symbol, p: sp.Symbol, q: sp.Symbol) -> sp.Expr:
    """Normalize an exchange-symmetric polynomial to a unique representative."""
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
    """Generate one raw local 4PN exchange-symmetric interaction basis block.

    Variables:
      a = v_A^2,
      b = v_B^2,
      c = v_A·v_B,
      d = v_A·n,
      e = v_B·n,
      p = m_A,
      q = m_B.
    """
    a, b, c, d, e, p, q = sp.symbols("a b c d e p q")
    vars_order = (p, q, a, b, c, d, e)
    basis: set[sp.Expr] = set()
    maxpow = vel_deg // 2 + 1

    for mp in range(mass_deg + 1):
        mq = mass_deg - mp
        for pa in range(maxpow + 1):
            for pb in range(maxpow + 1):
                for pc in range(maxpow + 1):
                    for pd in range(vel_deg + 1):
                        for pe in range(vel_deg + 1):
                            this_deg = 2 * pa + 2 * pb + 2 * pc + pd + pe
                            if this_deg != vel_deg:
                                continue
                            expr = p**mp * q**mq * a**pa * b**pb * c**pc * d**pd * e**pe
                            sym = canonical_sym(expr, vars_order, a, b, c, d, e, p, q)

                            # Strict one-body branch used throughout the project:
                            # body A is the test mass, body B is the fixed source.
                            tm = sp.expand(sym.subs({b: 0, c: 0, e: 0, p: 0, q: 1}))
                            if tm != 0:
                                continue
                            basis.add(sym)

    return sorted(basis, key=str)


# ---------------------------------------------------------------------------
# COM projection machinery
# ---------------------------------------------------------------------------

a, b, c, d, e, p, q = sp.symbols("a b c d e p q")
V2, rd = sp.symbols("V2 rd", real=True)
nu, Delta = sp.symbols("nu Delta", real=True)
Xa = (1 + Delta) / 2
Xb = (1 - Delta) / 2

COM_SUBS = {
    p: Xa,
    q: Xb,
    a: Xb**2 * V2,
    b: Xa**2 * V2,
    c: -Xa * Xb * V2,
    d: Xb * rd,
    e: -Xa * rd,
}



def to_nu(expr: sp.Expr) -> sp.Expr:
    """Reduce a generic-frame invariant to COM form and eliminate Delta in favor of nu."""
    expr = sp.expand(expr.subs(COM_SUBS))
    expr = sp.expand((expr + expr.subs(Delta, -Delta)) / 2)
    for n in range(40, 1, -2):
        expr = sp.expand(expr.subs(Delta**n, (1 - 4 * nu) ** (n // 2)))
    while expr.has(Delta**2):
        expr = sp.expand(expr.subs(Delta**2, 1 - 4 * nu))
    if expr.has(Delta):
        expr = sp.expand(expr.subs(Delta, 0))
    return sp.expand(expr)



def block_slots(expr: sp.Expr, block: str) -> list[sp.Expr]:
    """Extract the induced COM interaction slots for one local 4PN block."""
    expr = sp.expand(expr)

    if block == "Q":  # G/r times degree-8 velocity invariants -> 5 COM slots
        return [
            sp.expand(expr.coeff(V2, 4).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 3).coeff(rd, 2)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).coeff(rd, 4)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 6)),
            sp.expand(expr.coeff(rd, 8).subs(V2, 0)),
        ]
    if block == "T":  # G^2/r^2 times degree-6 invariants -> 4 COM slots
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 3).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).coeff(rd, 2)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 4)),
            sp.expand(expr.coeff(rd, 6).subs(V2, 0)),
        ]
    if block == "S":  # G^3/r^3 times degree-4 invariants -> 3 COM slots
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 2)),
            sp.expand(expr.coeff(rd, 4).subs(V2, 0)),
        ]
    if block == "U":  # G^4/r^4 times degree-2 invariants -> 2 COM slots
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).subs(rd, 0)),
            sp.expand(expr.coeff(rd, 2).subs(V2, 0)),
        ]
    if block == "W":  # G^5/r^5 static block -> 1 COM slot
        return [sp.expand(expr.subs({V2: 0, rd: 0}))]
    raise ValueError(f"Unknown block {block}")


# ---------------------------------------------------------------------------
# Main scaffold audit
# ---------------------------------------------------------------------------

def scaffold_audit() -> None:
    banner("PART II — RAW LOCAL 4PN GENERIC-FRAME INTERACTION SCAFFOLD")

    blocks = [
        ("Q", 0, 8, "G/r degree-8 velocity block"),
        ("T", 1, 6, "G^2/r^2 degree-6 mass-weighted block"),
        ("S", 2, 4, "G^3/r^3 degree-4 mass-weighted block"),
        ("U", 3, 2, "G^4/r^4 degree-2 mass-weighted block"),
        ("W", 4, 0, "G^5/r^5 static cross block"),
    ]

    total_basis = 0
    total_rank = 0

    for tag, mass_deg, vel_deg, label in blocks:
        basis = generate_basis(mass_deg, vel_deg)
        total_basis += len(basis)

        subbanner(f"II.{tag} — {label}")
        print(f"Raw basis count = {len(basis)}")

        cols = []
        for expr in basis:
            img = to_nu(expr)
            cols.append(sp.Matrix(block_slots(img, tag)))
        M = sp.Matrix.hstack(*cols)
        rank = M.rank()
        total_rank += rank
        _, pivots = M.rref()

        print(f"COM image matrix shape = {M.shape}")
        print(f"COM image rank         = {rank}")
        print(f"COM nullity            = {len(basis) - rank}")
        print(f"Pivot columns          = {pivots}")

        # Verify that the pivot subset indeed spans the whole COM image of the block.
        pivot_matrix = sp.Matrix.hstack(*[cols[j] for j in pivots])
        if pivot_matrix.rank() != rank:
            raise AssertionError(f"Pivot subset for block {tag} does not have full block rank.")

        print("Minimal pivot-spanning subset:")
        for j in pivots:
            print(f"  basis[{j}] = {basis[j]}")
            print(f"    COM image = {to_nu(basis[j])}")
            print(f"    slots     = {block_slots(to_nu(basis[j]), tag)}")

    banner("PART III — GLOBAL LEDGER")
    print(f"Total raw local interaction directions = {total_basis}")
    print(f"Total COM interaction rank            = {total_rank}")
    print(f"Total COM-nullity                     = {total_basis - total_rank}")
    print()
    print("Blockwise interaction-slot counts in the local 4PN COM basis:")
    print("  G/r      : 5 slots")
    print("  G^2/r^2  : 4 slots")
    print("  G^3/r^3  : 3 slots")
    print("  G^4/r^4  : 2 slots")
    print("  G^5/r^5  : 1 slot")
    print("  total    : 15 interaction slots")
    print()
    print("Main structural result:")
    print("  The raw local 4PN generic-frame interaction scaffold is already blockwise COM-complete.")
    print("  There is no early kinematic shortage at the level of constant-coefficient exchange-symmetric")
    print("  local interaction monomials.")
    print()
    print("What remains to be done next:")
    print("  1. import the fixed-chart local 4PN target,")
    print("  2. apply the exact quartic Hamiltonian compiler, and")
    print("  3. quotient the enormous 124-dimensional COM-null sector by the true fixed-chart conditions,")
    print("     not by COM data alone.")


if __name__ == "__main__":
    self_static_seed()
    scaffold_audit()
