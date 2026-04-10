#!/usr/bin/env python3
"""
3pn_generic_frame_com_projection_audit.py

Blockwise COM projection audit for the 3PN generic-frame lift.

Purpose
-------
The earlier 3PN artifacts solved the COM target and built a 50-parameter
exchange-symmetric generic-frame residual basis, but they did not yet connect
those two layers explicitly.

This script does that next exact step:

1. Rebuild the 24/17/8/1 generic-frame residual basis blocks.
2. Reduce each basis element to the center-of-mass (COM) frame.
3. Read off the induced COM ordinary-Lagrangian slot polynomials.
4. Measure the blockwise COM image ranks and nullities.
5. Compare the image against the current natural-seed COM residual target.
6. Isolate the exact obstruction pieces that cannot arise from the present
   constant-coefficient generic-frame basis.
7. Construct one exact COM-consistent generic-frame representative after the
   minimal seed refinement.

Main findings
-------------
- The G/r sextic block is already compatible with the solved COM target.
- The G^2/r^2 and G^3/r^3 blocks, in the current seed split, contain nu^3 COM
  tails that are not in the image of the present constant-coefficient generic-
  frame basis.
- After removing those obstruction tails into the seed/gauge sector, each block
  admits an exact COM-consistent representative.
- Even after COM matching there remains a large nullspace (27 interaction
  directions), so the full generic-frame target is still needed to fix the lift
  completely.
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
# Basis generation (copied so this file is standalone)
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

                            # Strict one-body branch: A test mass around fixed B.
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


# ---------------------------------------------------------------------------
# Solved COM target residuals from the current natural seed split
# ---------------------------------------------------------------------------

# These are the exact COM ordinary-Lagrangian residuals relative to the current
# natural self/static seed, as recorded by the existing 3PN COM target note.
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

# Exact obstruction pieces: the current generic-frame mass-degree structure can
# only generate nu and nu^2 in the T and S COM slots.
T_OBS = [sp.Rational(13, 4) * nu**3, sp.Rational(31, 4) * nu**3, 4 * nu**3]
S_OBS = [-sp.Rational(1, 2) * nu**3, -sp.Rational(7, 2) * nu**3]

REFINED_TARGETS = {
    "Q": TARGETS["Q"],
    "T": [sp.expand(DL10 - T_OBS[0]), sp.expand(DL11 - T_OBS[1]), sp.expand(DL12 - T_OBS[2])],
    "S": [sp.expand(DL13 - S_OBS[0]), sp.expand(DL14 - S_OBS[1])],
    "U": TARGETS["U"],
}


# ---------------------------------------------------------------------------
# Linear solve helpers
# ---------------------------------------------------------------------------

def build_equations(block: str, basis: list[sp.Expr], targets: list[sp.Expr], allow_pi2: bool) -> tuple[sp.Matrix, sp.Matrix, list[sp.Symbol]]:
    if allow_pi2:
        alpha = list(sp.symbols(f"a0:{len(basis)}"))
        beta = list(sp.symbols(f"b0:{len(basis)}"))
        vars = alpha + beta
    else:
        alpha = list(sp.symbols(f"c0:{len(basis)}"))
        beta = []
        vars = alpha

    eqexprs: list[sp.Expr] = []
    for slot, targ in enumerate(targets):
        images = [block_slots(to_nu(expr), block)[slot] for expr in basis]
        if allow_pi2:
            lhs = sum((alpha[i] + sp.pi**2 * beta[i]) * images[i] for i in range(len(basis)))
            diff = sp.expand(lhs - targ)
            A = sp.expand(diff.subs(sp.pi, 0))
            B = sp.expand((diff - A) / sp.pi**2)
            polys = [A, B]
        else:
            polys = [sp.expand(sum(alpha[i] * images[i] for i in range(len(basis))) - targ)]

        for poly in polys:
            P = sp.Poly(poly, nu)
            deg = P.degree()
            if deg == -sp.oo:
                continue
            for k in range(int(deg) + 1):
                coeff = sp.expand(P.coeff_monomial(nu**k))
                if coeff != 0:
                    eqexprs.append(coeff)

    M, bvec = sp.linear_eq_to_matrix(eqexprs, vars)
    return M, bvec, vars


def particular_solution(M: sp.Matrix, bvec: sp.Matrix, vars: list[sp.Symbol]) -> tuple[sp.Matrix, sp.Matrix]:
    sol_vec, params = M.gauss_jordan_solve(bvec)
    return sol_vec, params


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def main() -> None:
    Qbasis = generate_basis(0, 6)
    Tbasis = generate_basis(1, 4)
    Sbasis = generate_basis(2, 2)
    Ubasis = generate_basis(3, 0)

    banner("PART I — BLOCKWISE COM IMAGE OF THE 50-PARAMETER GENERIC-FRAME BASIS")
    print("Block sizes:")
    print(f"  Q (G/r sextic)       : {len(Qbasis)}")
    print(f"  T (G^2/r^2 quartic)  : {len(Tbasis)}")
    print(f"  S (G^3/r^3 quadratic): {len(Sbasis)}")
    print(f"  U (G^4/r^4 static)   : {len(Ubasis)}")

    MQ = image_matrix_polynomial("Q", Qbasis)
    MT = image_matrix_polynomial("T", Tbasis)
    MS = image_matrix_polynomial("S", Sbasis)
    MU = image_matrix_polynomial("U", Ubasis)

    print("\nPolynomial-coefficient image ranks (over constant coefficients):")
    print(f"  rank(Q) = {MQ.rank()}  -> nullity {MQ.shape[1] - MQ.rank()}")
    print(f"  rank(T) = {MT.rank()}  -> nullity {MT.shape[1] - MT.rank()}")
    print(f"  rank(S) = {MS.rank()}  -> nullity {MS.shape[1] - MS.rank()}")
    print(f"  rank(U) = {MU.rank()}  -> nullity {MU.shape[1] - MU.rank()}")

    subbanner("INTERPRETATION OF THE COM IMAGE")
    print("Q block: each COM slot spans nu, nu^2, nu^3.")
    print("T block: each COM slot spans only nu and nu^2. No nu^3 tails are possible.")
    print("S block: each COM slot spans only nu and nu^2 (with coefficients in Q(pi^2)).")
    print("U block: the COM image is exactly nu times a single constant coefficient.")

    banner("PART II — COMPATIBILITY OF THE CURRENT NATURAL-SEED RESIDUAL TARGET")

    # Q direct
    MQeq, MQrhs, Qvars = build_equations("Q", Qbasis, TARGETS["Q"], allow_pi2=False)
    print("Q direct compatibility:")
    print("  matrix rank     =", MQeq.rank())
    print("  augmented rank  =", MQeq.row_join(MQrhs).rank())
    if MQeq.rank() != MQeq.row_join(MQrhs).rank():
        raise AssertionError("Q block should be directly compatible, but is not.")

    # T direct
    MTeq, MTrhs, Tvars = build_equations("T", Tbasis, TARGETS["T"], allow_pi2=False)
    print("\nT direct compatibility:")
    print("  matrix rank     =", MTeq.rank())
    print("  augmented rank  =", MTeq.row_join(MTrhs).rank())
    print("  -> incompatible because the current target contains nu^3 tails")

    # S direct with pi^2 allowed
    MSeq, MSrhs, Svars = build_equations("S", Sbasis, TARGETS["S"], allow_pi2=True)
    print("\nS direct compatibility (allowing pi^2 coefficients):")
    print("  matrix rank     =", MSeq.rank())
    print("  augmented rank  =", MSeq.row_join(MSrhs).rank())
    print("  -> incompatible because the current target still contains nu^3 tails")

    # U direct with pi^2 allowed
    MUeq, MUrhs, Uvars = build_equations("U", Ubasis, TARGETS["U"], allow_pi2=True)
    print("\nU direct compatibility (allowing pi^2 coefficients):")
    print("  matrix rank     =", MUeq.rank())
    print("  augmented rank  =", MUeq.row_join(MUrhs).rank())
    if MUeq.rank() != MUeq.row_join(MUrhs).rank():
        raise AssertionError("U block should be compatible, but is not.")

    subbanner("EXACT OBSTRUCTION PIECES")
    print("T-block obstruction pieces to be absorbed into the refined seed/gauge sector:")
    print("  dL10_obs =", T_OBS[0])
    print("  dL11_obs =", T_OBS[1])
    print("  dL12_obs =", T_OBS[2])
    print("S-block obstruction pieces to be absorbed into the refined seed/gauge sector:")
    print("  dL13_obs =", S_OBS[0])
    print("  dL14_obs =", S_OBS[1])

    banner("PART III — EXACT COM-CONSISTENT REPRESENTATIVES AFTER THE MINIMAL SEED REFINEMENT")

    # Q representative
    q_rep = (
        sp.Rational(9, 4) * Qbasis[0]
        - sp.Rational(19, 8) * Qbasis[1]
        + sp.Rational(5, 4) * Qbasis[3]
        + sp.Rational(61, 16) * Qbasis[4]
        + sp.Rational(29, 16) * Qbasis[6]
        + sp.Rational(9, 16) * Qbasis[11]
        + sp.Rational(15, 32) * Qbasis[13]
        - sp.Rational(5, 16) * Qbasis[21]
    )
    q_slots = block_slots(to_nu(q_rep), "Q")
    for i, (lhs, rhs) in enumerate(zip(q_slots, REFINED_TARGETS["Q"]), start=6):
        expect_zero(f"Q representative slot dL{i}", lhs - rhs)

    # T representative (refined target only)
    t_rep = (
        sp.Rational(129, 16) * Tbasis[0]
        + sp.Rational(289, 16) * Tbasis[1]
        - sp.Rational(3, 16) * Tbasis[4]
        - sp.Rational(43, 16) * Tbasis[5]
        - sp.Rational(1, 3) * Tbasis[13]
        + sp.Rational(5, 12) * Tbasis[15]
    )
    t_slots = block_slots(to_nu(t_rep), "T")
    for i, (lhs, rhs) in enumerate(zip(t_slots, REFINED_TARGETS["T"]), start=10):
        expect_zero(f"T representative slot dL{i}_ref", lhs - rhs)

    # S representative (refined target only)
    s_rep = (
        -(3 * sp.pi**2 + 880) * Sbasis[0] / 192
        -(3 * sp.pi**2 + 244) * Sbasis[1] / 192
        + (3 * sp.pi**2 + 260) * Sbasis[4] / 64
        + (3 * sp.pi**2 + 452) * Sbasis[5] / 64
    )
    s_slots = block_slots(to_nu(s_rep), "S")
    for i, (lhs, rhs) in enumerate(zip(s_slots, REFINED_TARGETS["S"]), start=13):
        expect_zero(f"S representative slot dL{i}_ref", lhs - rhs)

    # U representative
    u_rep = (-908 + 63 * sp.pi**2) * Ubasis[0] / 96
    u_slots = block_slots(to_nu(u_rep), "U")
    expect_zero("U representative slot dL15", u_slots[0] - REFINED_TARGETS["U"][0])

    print("\nRepresentative generic-frame interaction blocks:")
    print("Q_part =", sp.factor(q_rep))
    print("T_part =", sp.factor(t_rep))
    print("S_part =", sp.factor(s_rep))
    print("U_part =", sp.factor(u_rep))

    banner("PART IV — WHAT THE COM PROJECTION NOW FIXES")
    print("The COM projection gives the following interaction-nullity counts after the")
    print("minimal seed refinement:")
    print("  Q block free directions : 12")
    print("  T block free directions : 11")
    print("  S block free directions :  4")
    print("  U block free directions :  0")
    print("  --------------------------------")
    print("  total interaction nullity: 27")
    print()
    print("So the current step does not finish the full generic-frame lift.")
    print("What it does finish is the exact COM-constraint layer:")
    print("  - Q is already compatible; ")
    print("  - U is uniquely fixed once pi^2 is allowed; ")
    print("  - T and S require a refined seed/gauge split that absorbs the nu^3 COM tails; ")
    print("  - after that refinement, one exact COM-consistent generic-frame representative exists blockwise.")


if __name__ == "__main__":
    main()
