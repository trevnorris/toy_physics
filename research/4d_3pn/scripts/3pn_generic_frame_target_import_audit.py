#!/usr/bin/env python3
"""
3pn_generic_frame_target_import_audit.py

Import the full generic-frame ordinary ADM-type 3PN Lagrangian target from
de Andrade, Blanchet, and Faye Eq. (4.11), substitute the GR value
lambda = -1987/3080, decompose the result on the 24/17/8/1 generic-frame
interaction scaffold, and compare its naive COM reduction with the previously
solved COM ordinary target.

Main point
----------
The full generic-frame target is imported exactly, but its naive COM reduction
does not reproduce the earlier COM ordinary target.  This is a genuine
reduction-ordering subtlety, consistent with the literature warning that one
cannot obtain the COM relative ordinary Lagrangian by naive substitution into
the general-frame one.
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
    expr = sp.expand(sp.simplify(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


a, b, c, d, e, p, q, r, G = sp.symbols("a b c d e p q r G")
nu, Delta, V2, rd = sp.symbols("nu Delta V2 rd")
P = sp.Symbol("P")


def swap_even(expr: sp.Expr) -> sp.Expr:
    return sp.expand(expr.xreplace({a: b, b: a, c: c, d: e, e: d, p: q, q: p}))


def canonical_sym(expr: sp.Expr, vars_order: tuple[sp.Symbol, ...]) -> sp.Expr:
    s = sp.expand(expr + swap_even(expr))
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
    vars_order = (p, q, a, b, c, d, e)
    basis: set[sp.Expr] = set()

    for mp in range(mass_deg + 1):
        mq = mass_deg - mp
        for pa in range(10):
            for pb in range(10):
                for pc in range(10):
                    for pd in range(vel_deg + 1):
                        for pe in range(vel_deg + 1):
                            this_deg = 2 * pa + 2 * pb + 2 * pc + pd + pe
                            if this_deg != vel_deg:
                                continue
                            expr = p**mp * q**mq * a**pa * b**pb * c**pc * d**pd * e**pe
                            sym = canonical_sym(expr, vars_order)
                            tm = sp.expand(sym.subs({b: 0, c: 0, e: 0, p: 0, q: 1}))
                            if tm != 0:
                                continue
                            basis.add(sym)

    return sorted(basis, key=str)


Qbasis = generate_basis(0, 6)
Tbasis = generate_basis(1, 4)
Sbasis = generate_basis(2, 2)
Ubasis = generate_basis(3, 0)


def coeff_dict(expr: sp.Expr) -> dict[tuple[int, ...], sp.Expr]:
    poly = sp.Poly(sp.expand(expr).subs(sp.pi**2, P), p, q, a, b, c, d, e, domain="EX")
    return {mon: coef for mon, coef in poly.terms()}


def coordinate_matrix(basis: list[sp.Expr]) -> tuple[sp.Matrix, dict[tuple[int, ...], int]]:
    monmap: dict[tuple[int, ...], int] = {}
    for expr in basis:
        for mon in coeff_dict(expr):
            monmap.setdefault(mon, len(monmap))
    M = sp.zeros(len(monmap), len(basis))
    for j, expr in enumerate(basis):
        for mon, coef in coeff_dict(expr).items():
            M[monmap[mon], j] = sp.expand(coef)
    return M, monmap


def coords_in_basis(expr: sp.Expr, basisM: sp.Matrix, monmap: dict[tuple[int, ...], int]) -> sp.Matrix:
    vec = sp.zeros(len(monmap), 1)
    for mon, coef in coeff_dict(expr).items():
        if mon not in monmap:
            raise ValueError(f"monomial {mon} not in basis set")
        vec[monmap[mon], 0] = sp.expand(coef)
    sol, params = basisM.gauss_jordan_solve(vec)
    if params.rows:
        sol = sol.subs({params[i, 0]: 0 for i in range(params.rows)})
    return sp.Matrix(sol)


QB, Qmap = coordinate_matrix(Qbasis)
TB, Tmap = coordinate_matrix(Tbasis)
SB, Smap = coordinate_matrix(Sbasis)
UB, Umap = coordinate_matrix(Ubasis)


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


def nonzero_terms(coords: sp.Matrix, basis: list[sp.Expr]) -> list[tuple[int, sp.Expr, sp.Expr]]:
    out: list[tuple[int, sp.Expr, sp.Expr]] = []
    for i, (coef, expr) in enumerate(zip(coords, basis)):
        coef = sp.expand(coef).subs(P, sp.pi**2)
        if coef != 0:
            out.append((i, sp.simplify(coef), expr))
    return out


def main() -> None:
    banner("PART I — IMPORT THE EXACT GENERIC-FRAME 3PN TARGET")

    lam = sp.Rational(-1987, 3080)
    print("lambda =", lam)

    Q_disp = (
        -sp.Rational(5, 32) * d**3 * e**3
        + sp.Rational(3, 16) * d**2 * e**2 * a
        + sp.Rational(9, 16) * d * e**3 * a
        - sp.Rational(3, 16) * d * e * a**2
        - sp.Rational(5, 16) * e**2 * a**2
        + sp.Rational(11, 16) * a**3
        - sp.Rational(15, 32) * d**2 * e**2 * c
        + sp.Rational(3, 4) * d * e * a * c
        - sp.Rational(1, 16) * e**2 * a * c
        - sp.Rational(21, 16) * a**2 * c
        + sp.Rational(5, 16) * d * e * c**2
        + sp.Rational(1, 8) * a * c**2
        + sp.Rational(1, 16) * c**3
        - sp.Rational(5, 16) * d**2 * a * b
        - sp.Rational(9, 32) * d * e * a * b
        + sp.Rational(7, 8) * a**2 * b
        - sp.Rational(15, 32) * a * b * c
    )

    T_disp = (
        -sp.Rational(5, 12) * d**4
        - sp.Rational(13, 8) * d**3 * e
        - sp.Rational(23, 24) * d**2 * e**2
        + sp.Rational(13, 16) * d**2 * a
        + sp.Rational(1, 4) * d * e * a
        + sp.Rational(5, 6) * e**2 * a
        + sp.Rational(21, 16) * a**2
        - sp.Rational(1, 2) * d**2 * c
        + sp.Rational(1, 3) * d * e * c
        - sp.Rational(97, 16) * a * c
        + sp.Rational(341, 48) * c**2
        + sp.Rational(29, 24) * d**2 * b
        - d * e * b
        + sp.Rational(43, 12) * a * b
        - sp.Rational(71, 8) * b * c
        + sp.Rational(47, 16) * b**2
    )

    S22_disp = (
        (sp.Rational(73, 16) + sp.Rational(3, 64) * sp.pi**2) * d**2
        + (-11 - sp.Rational(3, 64) * sp.pi**2) * d * e
        + (-sp.Rational(265, 48) - sp.Rational(1, 64) * sp.pi**2) * a
        + (sp.Rational(59, 8) + sp.Rational(1, 64) * sp.pi**2) * c
    )
    S31_disp = (
        -5 * d**2 - sp.Rational(1, 8) * d * e + sp.Rational(173, 48) * a
        - sp.Rational(27, 8) * c + sp.Rational(13, 8) * b
    )
    U41_disp = -sp.Rational(1, 8)
    U32_disp = sp.simplify(-sp.Rational(993, 140) + sp.Rational(11, 3) * lam + sp.Rational(21, 32) * sp.pi**2)

    displayed = (
        sp.Rational(5, 128) * p * a**4
        + G * p * q / r * Q_disp
        + G**2 * p**2 * q / r**2 * T_disp
        + G**3 * p**2 * q**2 / r**3 * S22_disp
        + G**3 * p**3 * q / r**3 * S31_disp
        + G**4 * p**4 * q / r**4 * U41_disp
        + G**4 * p**3 * q**2 / r**4 * U32_disp
    )
    L3_full = sp.expand(displayed + swap_even(displayed))

    subbanner("Split into 3PN blocks and subtract the frozen natural seed")
    def split_blocks(expr: sp.Expr) -> dict[tuple[int, int], sp.Expr]:
        expr = sp.expand(expr)
        out: dict[tuple[int, int], sp.Expr] = {}
        for term in sp.Add.make_args(expr):
            pd = term.as_powers_dict()
            gpow = int(pd.get(G, 0))
            rpow = int(pd.get(r, 0))
            coeff = sp.simplify(term / (G**gpow * r**rpow))
            out[(gpow, rpow)] = sp.expand(out.get((gpow, rpow), 0) + coeff)
        return out

    blocks = split_blocks(L3_full)
    Q_full = sp.expand(blocks[(1, -1)] / (p * q))
    T_full = sp.expand(blocks[(2, -2)] / (p * q))
    S_full = sp.expand(blocks[(3, -3)] / (p * q))
    U_full = sp.expand(blocks[(4, -4)] / (p * q))

    Q_seed = sp.Rational(11, 16) * (a**3 + b**3)
    T_seed = sp.Rational(47, 16) * (q * a**2 + p * b**2)
    S_seed = sp.Rational(13, 8) * (q**2 * a + p**2 * b)
    U_seed = -sp.Rational(1, 8) * (q**3 + p**3)

    Q_res = sp.expand(Q_full - Q_seed)
    T_res = sp.expand(T_full - T_seed)
    S_res = sp.expand(S_full - S_seed)
    U_res = sp.expand(U_full - U_seed)

    Qcoords = coords_in_basis(Q_res, QB, Qmap)
    Tcoords = coords_in_basis(T_res, TB, Tmap)
    Scoords = coords_in_basis(S_res, SB, Smap)
    Ucoords = coords_in_basis(U_res, UB, Umap)

    print("Nonzero Q residual coefficients:")
    for i, coef, expr in nonzero_terms(Qcoords, Qbasis):
        print(f"  Q[{i}] = {sp.simplify(coef)}   on   {expr}")

    print("\nNonzero T residual coefficients:")
    for i, coef, expr in nonzero_terms(Tcoords, Tbasis):
        print(f"  T[{i}] = {sp.simplify(coef)}   on   {expr}")

    print("\nNonzero S residual coefficients:")
    for i, coef, expr in nonzero_terms(Scoords, Sbasis):
        print(f"  S[{i}] = {sp.simplify(coef)}   on   {expr}")

    print("\nNonzero U residual coefficient:")
    for i, coef, expr in nonzero_terms(Ucoords, Ubasis):
        print(f"  U[{i}] = {sp.simplify(coef)}   on   {expr}")

    banner("PART II — STRICT ONE-BODY GATE")
    tm_subs = {b: 0, c: 0, e: 0, p: 0, q: 1}
    expect_zero("Q_res test-mass limit", Q_res.subs(tm_subs))
    expect_zero("T_res test-mass limit", T_res.subs(tm_subs))
    expect_zero("S_res test-mass limit", S_res.subs(tm_subs))
    expect_zero("U_res test-mass limit", U_res.subs(tm_subs))

    banner("PART III — NAIVE COM REDUCTION MISMATCH")
    targets = {
        "Q": [
            sp.Rational(38, 16) * nu - sp.Rational(116, 16) * nu**2 - sp.Rational(57, 16) * nu**3,
            sp.Rational(20, 16) * nu**2 - sp.Rational(69, 16) * nu**3,
            sp.Rational(9, 16) * nu**2 - sp.Rational(33, 16) * nu**3,
            sp.Rational(5, 16) * nu**3,
        ],
        "T": [
            sp.Rational(129, 16) * nu - sp.Rational(98, 16) * nu**2 + sp.Rational(52, 16) * nu**3,
            -sp.Rational(3, 16) * nu + sp.Rational(52, 16) * nu**2 + sp.Rational(124, 16) * nu**3,
            -sp.Rational(5, 12) * nu + sp.Rational(11, 12) * nu**2 + 4 * nu**3,
        ],
        "S": [
            -sp.Rational(244, 192) * nu - sp.Rational(3, 192) * sp.pi**2 * nu - sp.Rational(1272, 192) * nu**2 - sp.Rational(96, 192) * nu**3,
            sp.Rational(452, 64) * nu + sp.Rational(3, 64) * sp.pi**2 * nu - 6 * nu**2 - sp.Rational(7, 2) * nu**3,
        ],
        "U": [
            (-sp.Rational(908, 96) + sp.Rational(63, 96) * sp.pi**2) * nu,
        ],
    }

    q_slots = block_slots(to_nu(Q_res), "Q")
    t_slots = block_slots(to_nu(T_res), "T")
    s_slots = block_slots(to_nu(S_res), "S")
    u_slots = block_slots(to_nu(U_res), "U")

    for i, (lhs, rhs) in enumerate(zip(q_slots, targets["Q"]), start=6):
        print(f"dL{i} mismatch =", sp.simplify(lhs - rhs))
    for i, (lhs, rhs) in enumerate(zip(t_slots, targets["T"]), start=10):
        print(f"dL{i} mismatch =", sp.simplify(lhs - rhs))
    for i, (lhs, rhs) in enumerate(zip(s_slots, targets["S"]), start=13):
        print(f"dL{i} mismatch =", sp.simplify(lhs - rhs))
    print("dL15 mismatch =", sp.simplify(u_slots[0] - targets["U"][0]))

    banner("PART IV — THEOREM LEDGER")
    print("1. Eq. (4.11) imports the exact generic-frame ordinary ADM-type 3PN target.")
    print("2. After subtracting the frozen natural one-body/self-static seed, the target lies exactly")
    print("   in the compact 24/17/8/1 generic-frame interaction scaffold.")
    print("3. The strict test-mass gate remains clean.")
    print("4. However, the naive COM reduction of this generic-frame ordinary target does not reproduce")
    print("   the previously solved COM ordinary target in the Q,T,S blocks.")
    print("5. Therefore the remaining 3PN quotient cannot be fixed by naive COM substitution alone:")
    print("   one needs either the full Hamiltonian-level lift or the true generic-to-COM ordinary")
    print("   reduction map.")
    print("6. This is consistent with the literature warning that the COM relative ordinary Lagrangian")
    print("   does not straightforwardly follow from the general-frame one by naive reduction.")


if __name__ == "__main__":
    main()
