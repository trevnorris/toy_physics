#!/usr/bin/env python3
"""
3pn_seed_repair_and_com_null_ideal_audit.py

Standalone SymPy audit for the next exact step of the 3PN generic-frame lift.

What this script does
---------------------
1. Rebuild the overcomplete 24/17/8/1 generic-frame residual basis.
2. Rebuild the blockwise COM projection maps.
3. Construct the minimal nu-dressed seed-repair blocks that absorb the middle-
   block nu^3 obstructions found in the earlier COM projection audit.
4. Verify that those repair blocks:
   - reproduce exactly the obstruction COM tails,
   - vanish in the strict one-body limit.
5. Compute the exact sparse nullspace generators of the Q/T/S COM projections.
6. Verify that
   - T and S null generators lie in the linear COM-momentum ideal,
   - Q null generators lie in the full COM-null ideal.
7. Build one explicit canonical-gauge generic-frame representative and verify
   that it reproduces the exact solved COM 3PN target blockwise.

The point is to turn the remaining generic-frame ambiguity into something sharp:
not an unexplained mismatch, but a tiny nu-dressed seed repair plus a clean
27-dimensional COM-null gauge/contact ideal.
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
    expr = sp.expand(sp.simplify(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Basis generation
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
                            tm = sp.expand(sym.subs({b: 0, c: 0, e: 0, p: 0, q: 1}))
                            if tm != 0:
                                continue
                            basis.add(sym)

    return sorted(basis, key=str)


# ---------------------------------------------------------------------------
# COM reduction machinery
# ---------------------------------------------------------------------------

a, b, c, d, e, p, q = sp.symbols("a b c d e p q")
V2, rd = sp.symbols("V2 rd")
nu, Delta = sp.symbols("nu Delta")
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
    expr = sp.expand(expr.subs(COM_SUBS))
    expr = sp.expand((expr + expr.subs(Delta, -Delta)) / 2)
    for n in range(20, 1, -2):
        expr = sp.expand(expr.subs(Delta**n, (1 - 4 * nu) ** (n // 2)))
    while expr.has(Delta**2):
        expr = sp.expand(expr.subs(Delta**2, 1 - 4 * nu))
    if expr.has(Delta):
        expr = sp.expand(expr.subs(Delta, 0))
    return sp.expand(expr)


def block_slots(expr: sp.Expr, block: str) -> list[sp.Expr]:
    expr = sp.expand(expr)
    if block == "Q":
        return [
            sp.expand(expr.coeff(V2, 3).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).coeff(rd, 2)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 4)),
            sp.expand(expr.coeff(rd, 6).subs(V2, 0)),
        ]
    if block == "T":
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 2)),
            sp.expand(expr.coeff(rd, 4).subs(V2, 0)),
        ]
    if block == "S":
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).subs(rd, 0)),
            sp.expand(expr.coeff(rd, 2).subs(V2, 0)),
        ]
    if block == "U":
        return [sp.expand(expr)]
    raise ValueError(block)


def image_matrix_polynomial(block: str, basis: list[sp.Expr]) -> sp.Matrix:
    rows: list[list[sp.Expr]] = []
    for expr in basis:
        slots = block_slots(to_nu(expr), block)
        row: list[sp.Expr] = []
        if block == "Q":
            for s in slots:
                P = sp.Poly(s, nu)
                row.extend([P.coeff_monomial(nu), P.coeff_monomial(nu**2), P.coeff_monomial(nu**3)])
        elif block == "T":
            for s in slots:
                P = sp.Poly(s, nu)
                row.extend([P.coeff_monomial(nu), P.coeff_monomial(nu**2)])
        elif block == "S":
            for s in slots:
                P = sp.Poly(s, nu)
                row.extend([P.coeff_monomial(nu), P.coeff_monomial(nu**2)])
        else:
            P = sp.Poly(slots[0], nu)
            row.append(P.coeff_monomial(nu))
        rows.append(row)
    return sp.Matrix(rows).T


def pivot_cols(M: sp.Matrix) -> tuple[int, ...]:
    return M.rref()[1]


# ---------------------------------------------------------------------------
# Earlier exact COM targets and obstruction pieces
# ---------------------------------------------------------------------------

DL6 = sp.expand(sp.Rational(38, 16) * nu - sp.Rational(116, 16) * nu**2 - sp.Rational(57, 16) * nu**3)
DL7 = sp.expand(sp.Rational(20, 16) * nu**2 - sp.Rational(69, 16) * nu**3)
DL8 = sp.expand(sp.Rational(9, 16) * nu**2 - sp.Rational(33, 16) * nu**3)
DL9 = sp.expand(sp.Rational(5, 16) * nu**3)

DL10 = sp.expand(sp.Rational(129, 16) * nu - sp.Rational(98, 16) * nu**2 + sp.Rational(52, 16) * nu**3)
DL11 = sp.expand(-sp.Rational(3, 16) * nu + sp.Rational(52, 16) * nu**2 + sp.Rational(124, 16) * nu**3)
DL12 = sp.expand(-sp.Rational(5, 12) * nu + sp.Rational(11, 12) * nu**2 + 4 * nu**3)

DL13 = sp.expand(-sp.Rational(244, 192) * nu - sp.Rational(3, 192) * sp.pi**2 * nu - sp.Rational(1272, 192) * nu**2 - sp.Rational(96, 192) * nu**3)
DL14 = sp.expand(sp.Rational(452, 64) * nu + sp.Rational(3, 64) * sp.pi**2 * nu - 6 * nu**2 - sp.Rational(7, 2) * nu**3)
DL15 = sp.expand((-sp.Rational(908, 96) + sp.Rational(63, 96) * sp.pi**2) * nu)

TARGETS = {
    "Q": [DL6, DL7, DL8, DL9],
    "T": [DL10, DL11, DL12],
    "S": [DL13, DL14],
    "U": [DL15],
}

T_OBS = [sp.Rational(13, 4) * nu**3, sp.Rational(31, 4) * nu**3, 4 * nu**3]
S_OBS = [-sp.Rational(1, 2) * nu**3, -sp.Rational(7, 2) * nu**3]


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def vec_terms(vec: sp.Matrix, basis: list[sp.Expr]) -> list[tuple[sp.Expr, sp.Expr]]:
    out: list[tuple[sp.Expr, sp.Expr]] = []
    for coeff, expr in zip(vec, basis):
        coeff = sp.simplify(coeff)
        if coeff != 0:
            out.append((coeff, expr))
    return out


def main() -> None:
    Qbasis = generate_basis(0, 6)
    Tbasis = generate_basis(1, 4)
    Sbasis = generate_basis(2, 2)
    Ubasis = generate_basis(3, 0)

    MQ = image_matrix_polynomial("Q", Qbasis)
    MT = image_matrix_polynomial("T", Tbasis)
    MS = image_matrix_polynomial("S", Sbasis)
    MU = image_matrix_polynomial("U", Ubasis)

    banner("PART I — BLOCKWISE BASIS SIZES, RANKS, AND QUOTIENT DIMENSIONS")
    print(f"Q block size = {len(Qbasis)}, rank = {MQ.rank()}, nullity = {len(Qbasis) - MQ.rank()}")
    print(f"T block size = {len(Tbasis)}, rank = {MT.rank()}, nullity = {len(Tbasis) - MT.rank()}")
    print(f"S block size = {len(Sbasis)}, rank = {MS.rank()}, nullity = {len(Sbasis) - MS.rank()}")
    print(f"U block size = {len(Ubasis)}, rank = {MU.rank()}, nullity = {len(Ubasis) - MU.rank()}")
    print("Q pivot columns:", pivot_cols(MQ))
    print("T pivot columns:", pivot_cols(MT))
    print("S pivot columns:", pivot_cols(MS))

    banner("PART II — EXACT MINIMAL NU-DRESSED SEED REPAIR")
    nu_mass = p * q / (p + q) ** 2

    # Compact obstruction repair blocks.
    DeltaT_nu = sp.expand(nu_mass * (
        sp.Rational(13, 4) * Tbasis[1] +
        sp.Rational(31, 4) * Tbasis[12] +
        4 * Tbasis[13]
    ))
    DeltaS_nu = sp.expand(nu_mass * (
        sp.Rational(1, 2) * Sbasis[3] +
        sp.Rational(7, 2) * Sbasis[7]
    ))

    print("DeltaT_nu =", sp.factor(DeltaT_nu))
    print("DeltaS_nu =", sp.factor(DeltaS_nu))

    # One-body branch must vanish.
    tm_subs = {b: 0, c: 0, e: 0, p: 0, q: 1}
    expect_zero("strict test-mass limit of DeltaT_nu", DeltaT_nu.subs(tm_subs))
    expect_zero("strict test-mass limit of DeltaS_nu", DeltaS_nu.subs(tm_subs))

    # COM slots must reproduce the exact obstruction pieces.
    Tcorr_slots = block_slots(to_nu(DeltaT_nu), "T")
    Scorr_slots = block_slots(to_nu(DeltaS_nu), "S")
    for i, (lhs, rhs) in enumerate(zip(Tcorr_slots, T_OBS), start=10):
        expect_zero(f"DeltaT_nu COM slot dL{i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(Scorr_slots, S_OBS), start=13):
        expect_zero(f"DeltaS_nu COM slot dL{i}", lhs - rhs)

    banner("PART III — THE COM-NULL IDEAL")
    C1 = p * a + q * c
    C2 = p * c + q * b
    C3 = p * d + q * e
    C4 = a * b - c**2
    C5 = a * e - c * d
    C6 = b * d - c * e

    print("C1 =", C1)
    print("C2 =", C2)
    print("C3 =", C3)
    print("C4 =", C4)
    print("C5 =", C5)
    print("C6 =", C6)

    G_lin = sp.groebner([C1, C2, C3], a, b, c, d, e, p, q, order="lex", domain=sp.QQ)
    G_full = sp.groebner([C1, C2, C3, C4, C5, C6], a, b, c, d, e, p, q, order="lex", domain=sp.QQ)

    T_null = MT.nullspace()
    S_null = MS.nullspace()
    Q_null = MQ.nullspace()

    subbanner("Sparse null generators")
    print("T-block sparse null generators:")
    for k, vec in enumerate(T_null, start=1):
        expr = sp.expand(sum(coeff * term for coeff, term in vec_terms(vec, Tbasis)))
        print(f"  N_T{k} = {expr}")
        rem = G_lin.reduce(expr)[1]
        expect_zero(f"N_T{k} in <C1,C2,C3>", rem)

    print("\nS-block sparse null generators:")
    for k, vec in enumerate(S_null, start=1):
        expr = sp.expand(sum(coeff * term for coeff, term in vec_terms(vec, Sbasis)))
        print(f"  N_S{k} = {expr}")
        rem = G_lin.reduce(expr)[1]
        expect_zero(f"N_S{k} in <C1,C2,C3>", rem)

    print("\nQ-block sparse null generators:")
    for k, vec in enumerate(Q_null, start=1):
        expr = sp.expand(sum(coeff * term for coeff, term in vec_terms(vec, Qbasis)))
        print(f"  N_Q{k} = {expr}")
        rem = G_full.reduce(expr)[1]
        expect_zero(f"N_Q{k} in full COM-null ideal", rem)

    banner("PART IV — CANONICAL-GAUGE GENERIC-FRAME REPRESENTATIVE")

    # Carried exact representatives from the earlier COM-projection note.
    Q_can = sp.expand(
        sp.Rational(9, 4) * Qbasis[0]
        - sp.Rational(19, 8) * Qbasis[1]
        + sp.Rational(5, 4) * Qbasis[3]
        + sp.Rational(61, 16) * Qbasis[4]
        + sp.Rational(29, 16) * Qbasis[6]
        + sp.Rational(9, 16) * Qbasis[11]
        + sp.Rational(15, 32) * Qbasis[13]
        - sp.Rational(5, 16) * Qbasis[21]
    )

    T_can = sp.expand(
        sp.Rational(129, 16) * Tbasis[0]
        + sp.Rational(289, 16) * Tbasis[1]
        - sp.Rational(3, 16) * Tbasis[4]
        - sp.Rational(43, 16) * Tbasis[5]
        - sp.Rational(1, 3) * Tbasis[13]
        + sp.Rational(5, 12) * Tbasis[15]
        + DeltaT_nu
    )

    S_can = sp.expand(
        -(3 * sp.pi**2 + 880) * Sbasis[0] / 192
        -(3 * sp.pi**2 + 244) * Sbasis[1] / 192
        + (3 * sp.pi**2 + 260) * Sbasis[4] / 64
        + (3 * sp.pi**2 + 452) * Sbasis[5] / 64
        + DeltaS_nu
    )

    U_can = sp.expand((-908 + 63 * sp.pi**2) * Ubasis[0] / 96)

    print("Q_can =", sp.factor(Q_can))
    print("T_can =", sp.factor(T_can))
    print("S_can =", sp.factor(S_can))
    print("U_can =", sp.factor(U_can))

    # Verify exact COM target recovery.
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(Q_can), "Q"), TARGETS["Q"]), start=6):
        expect_zero(f"Q_can COM slot dL{i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(T_can), "T"), TARGETS["T"]), start=10):
        expect_zero(f"T_can COM slot dL{i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(S_can), "S"), TARGETS["S"]), start=13):
        expect_zero(f"S_can COM slot dL{i}", lhs - rhs)
    expect_zero("U_can COM slot dL15", block_slots(to_nu(U_can), "U")[0] - TARGETS["U"][0])

    banner("PART V — THEOREM LEDGER")
    print("1. The middle-block nu^3 obstructions are absorbed exactly by a tiny nu-dressed")
    print("   seed-repair sector DeltaT_nu + DeltaS_nu.")
    print("2. That repair vanishes in the strict one-body branch.")
    print("3. The remaining 27 unfixed generic-frame directions are precisely COM-null")
    print("   directions: T/S lie in the linear COM-momentum ideal, while Q lies in the")
    print("   full COM-null ideal including collinearity relations.")
    print("4. Setting the COM-null directions to zero defines a canonical gauge slice.")
    print("5. On that slice, the explicit representative (Q_can,T_can,S_can,U_can)")
    print("   reproduces the exact solved COM 3PN target blockwise.")


if __name__ == "__main__":
    main()
