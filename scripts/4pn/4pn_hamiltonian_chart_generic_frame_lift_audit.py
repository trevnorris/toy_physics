#!/usr/bin/env python3
"""
4pn_hamiltonian_chart_generic_frame_lift_audit.py

Hamiltonian-chart generic-frame local 4PN lift audit.

What this script does
---------------------
1. Rebuilds the exchange-symmetric local 4PN generic-frame interaction scaffold:
      - G/r     times degree-8 momentum invariants,
      - G^2/r^2 times degree-6 invariants with one extra mass power,
      - G^3/r^3 times degree-4 invariants with two extra mass powers,
      - G^4/r^4 times degree-2 invariants with three extra mass powers,
      - G^5/r^5 static cross polynomials.
   The formal basis is the same as in the ordinary scaffold, but interpreted on the
   Hamiltonian side in Newtonian-order momentum variables.

2. Imports the fixed-chart *local* 4PN COM Hamiltonian target and subtracts the natural
   Hamiltonian self/static seed built from the exact one-body gate.

3. Decomposes each blockwise COM image into explicit nu-polynomial coefficient space:
      - G/r     : 5 slots x (nu,nu^2,nu^3,nu^4) = 20 coefficients,
      - G^2/r^2 : 4 slots x (nu,nu^2,nu^3)      = 12 coefficients,
      - G^3/r^3 : 3 slots x (nu,nu^2,nu^3)      = 9 coefficients,
      - G^4/r^4 : 2 slots x (nu,nu^2)           = 4 coefficients,
      - G^5/r^5 : 1 slot  x (nu,nu^2)           = 2 coefficients.

4. Computes the exact blockwise image ranks/nullities.

5. Solves for one exact generic-frame Hamiltonian representative of the local 4PN
   comparable-mass interaction residual and verifies that its COM reduction reproduces
   the imported target residual exactly.

Main result
-----------
In the Hamiltonian chart, the constant-coefficient local 4PN generic-frame scaffold is
not merely slot-complete after COM projection; it is polynomially complete against the
full fixed-chart local COM interaction target. Blockwise image ranks are maximal:

    G/r     : rank 20 / 20   (nullity 32)
    G^2/r^2 : rank 12 / 12   (nullity 34)
    G^3/r^3 : rank  9 /  9   (nullity 20)
    G^4/r^4 : rank  4 /  4   (nullity  6)
    G^5/r^5 : rank  2 /  2   (nullity  0)

So, unlike the ordinary-chart translation, the Hamiltonian chart shows no nu-tail
obstruction and no need for seed refinement before the local 4PN generic-frame lift.
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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.expand(sp.simplify(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.expand(sp.simplify(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Basis generation (standalone)
# ---------------------------------------------------------------------------

a, b, c, d, e, p, q = sp.symbols("a b c d e p q")
nu, Delta = sp.symbols("nu Delta", real=True)
V2, rd = sp.symbols("V2 rd", real=True)
P = sp.Symbol("P")  # placeholder for pi^2 in coordinate dictionaries

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


def swap_full(expr: sp.Expr) -> sp.Expr:
    return sp.expand(expr.xreplace({a: b, b: a, c: c, d: e, e: d, p: q, q: p}))


def canonical_sym(expr: sp.Expr, vars_order: tuple[sp.Symbol, ...]) -> sp.Expr:
    s = sp.expand(expr + swap_full(expr))
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
                            sym = canonical_sym(expr, vars_order)

                            # Strict one-body branch: A test mass orbiting fixed B.
                            tm = sp.expand(sym.subs({b: 0, c: 0, e: 0, p: 0, q: 1}))
                            if tm != 0:
                                continue
                            basis.add(sym)

    return sorted(basis, key=str)


def to_nu(expr: sp.Expr) -> sp.Expr:
    expr = sp.expand(expr.subs(COM_SUBS))
    expr = sp.expand((expr + expr.subs(Delta, -Delta)) / 2)
    for n in range(60, 1, -2):
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
            sp.expand(expr.coeff(V2, 4).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 3).coeff(rd, 2)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).coeff(rd, 4)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 6)),
            sp.expand(expr.coeff(rd, 8).subs(V2, 0)),
        ]
    if block == "T":
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 3).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).coeff(rd, 2)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 4)),
            sp.expand(expr.coeff(rd, 6).subs(V2, 0)),
        ]
    if block == "S":
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 2).subs(rd, 0)),
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).coeff(rd, 2)),
            sp.expand(expr.coeff(rd, 4).subs(V2, 0)),
        ]
    if block == "U":
        return [
            sp.expand(sp.collect(expr, V2).coeff(V2, 1).subs(rd, 0)),
            sp.expand(expr.coeff(rd, 2).subs(V2, 0)),
        ]
    if block == "W":
        return [sp.expand(expr.subs({V2: 0, rd: 0}))]
    raise ValueError(block)


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


def image_matrix_polynomial(block: str, basis: list[sp.Expr], maxdeg: int) -> sp.Matrix:
    rows: list[list[sp.Expr]] = []
    for expr in basis:
        slots = block_slots(to_nu(expr), block)
        row: list[sp.Expr] = []
        for s in slots:
            poly = sp.Poly(sp.expand(s), nu)
            for k in range(1, maxdeg + 1):
                row.append(sp.expand(poly.coeff_monomial(nu**k)))
        rows.append(row)
    return sp.Matrix(rows).T


def target_vector(residual_slots: list[sp.Expr], maxdeg: int) -> sp.Matrix:
    vec: list[sp.Expr] = []
    for s in residual_slots:
        poly = sp.Poly(sp.expand(s), nu)
        for k in range(1, maxdeg + 1):
            vec.append(sp.expand(poly.coeff_monomial(nu**k)))
    return sp.Matrix(vec)


def particular_solution(M: sp.Matrix, vec: sp.Matrix) -> sp.Matrix:
    sol, params = M.gauss_jordan_solve(vec)
    if params.rows:
        sol = sol.subs({params[i, 0]: 0 for i in range(params.rows)})
    return sp.Matrix(sol)


def nonzero_terms(coords: sp.Matrix, basis: list[sp.Expr]) -> list[tuple[int, sp.Expr, sp.Expr]]:
    out: list[tuple[int, sp.Expr, sp.Expr]] = []
    for i, (coef, expr) in enumerate(zip(coords, basis)):
        coef = sp.expand(coef).subs(P, sp.pi**2)
        if coef != 0:
            out.append((i, sp.simplify(coef), expr))
    return out


# ---------------------------------------------------------------------------
# Imported fixed local COM Hamiltonian target
# ---------------------------------------------------------------------------

def local_adm_target() -> dict[str, sp.Expr]:
    pi = sp.pi
    return {
        "K": sp.Rational(7, 256) - sp.Rational(63, 256) * nu + sp.Rational(189, 256) * nu**2 - sp.Rational(105, 128) * nu**3 + sp.Rational(63, 256) * nu**4,

        "Q1": sp.Rational(45, 128) - sp.Rational(45, 16) * nu + sp.Rational(423, 64) * nu**2 - sp.Rational(1013, 256) * nu**3 - sp.Rational(35, 128) * nu**4,
        "Q2": -sp.Rational(3, 32) * nu**2 + sp.Rational(23, 64) * nu**3 - sp.Rational(5, 32) * nu**4,
        "Q3": -sp.Rational(9, 64) * nu**2 + sp.Rational(69, 128) * nu**3 - sp.Rational(9, 64) * nu**4,
        "Q4": -sp.Rational(5, 64) * nu**3 - sp.Rational(5, 32) * nu**4,
        "Q5": sp.Rational(35, 256) * nu**3 - sp.Rational(35, 128) * nu**4,

        "T1": sp.Rational(13, 8) - sp.Rational(791, 64) * nu + sp.Rational(4857, 256) * nu**2 + sp.Rational(2335, 256) * nu**3,
        "T2": sp.Rational(49, 16) * nu - sp.Rational(545, 64) * nu**2 + sp.Rational(1135, 256) * nu**3,
        "T3": -sp.Rational(889, 192) * nu + sp.Rational(9475, 768) * nu**2 - sp.Rational(1649, 768) * nu**3,
        "T4": sp.Rational(369, 160) * nu - sp.Rational(1151, 128) * nu**2 + sp.Rational(10353, 1280) * nu**3,

        "S1": sp.Rational(105, 32) + (sp.Rational(2749, 8192) * pi**2 - sp.Rational(589189, 19200)) * nu + (sp.Rational(18491, 16384) * pi**2 - sp.Rational(1189789, 28800)) * nu**2 - sp.Rational(553, 128) * nu**3,
        "S2": (sp.Rational(63347, 1600) - sp.Rational(1059, 1024) * pi**2) * nu + (-sp.Rational(127, 3) - sp.Rational(4035, 2048) * pi**2) * nu**2 - sp.Rational(225, 64) * nu**3,
        "S3": (sp.Rational(375, 8192) * pi**2 - sp.Rational(23533, 1280)) * nu + (sp.Rational(57563, 1920) - sp.Rational(38655, 16384) * pi**2) * nu**2 - sp.Rational(381, 128) * nu**3,

        "U1": sp.Rational(105, 32) + (sp.Rational(185761, 19200) - sp.Rational(21837, 8192) * pi**2) * nu + (sp.Rational(672811, 19200) - sp.Rational(158177, 49152) * pi**2) * nu**2,
        "U2": (sp.Rational(3401779, 57600) - sp.Rational(28691, 24576) * pi**2) * nu + (sp.Rational(110099, 49152) * pi**2 - sp.Rational(21827, 3840)) * nu**2,

        "W1": -sp.Rational(1, 16) + (sp.Rational(6237, 1024) * pi**2 - sp.Rational(169199, 2400)) * nu + (sp.Rational(7403, 3072) * pi**2 - sp.Rational(1256, 45)) * nu**2,
    }


# ---------------------------------------------------------------------------
# Natural Hamiltonian seed and residuals
# ---------------------------------------------------------------------------

def to_nu_even(expr: sp.Expr) -> sp.Expr:
    expr = sp.expand((expr + expr.subs(Delta, -Delta)) / 2)
    for n in range(60, 1, -2):
        expr = sp.expand(expr.subs(Delta**n, (1 - 4 * nu) ** (n // 2)))
    while expr.has(Delta**2):
        expr = sp.expand(expr.subs(Delta**2, 1 - 4 * nu))
    if expr.has(Delta):
        expr = sp.expand(expr.subs(Delta, 0))
    return sp.expand(expr)


def hamiltonian_seed() -> dict[str, sp.Expr]:
    return {
        "K": to_nu_even(sp.Rational(7, 256) * (Xa**9 + Xb**9)),
        "Q1": to_nu_even(sp.Rational(45, 128) * (Xa**8 + Xb**8)),
        "Q2": sp.Integer(0),
        "Q3": sp.Integer(0),
        "Q4": sp.Integer(0),
        "Q5": sp.Integer(0),
        "T1": to_nu_even(sp.Rational(13, 8) * (Xa**7 + Xb**7)),
        "T2": sp.Integer(0),
        "T3": sp.Integer(0),
        "T4": sp.Integer(0),
        "S1": to_nu_even(sp.Rational(105, 32) * (Xa**6 + Xb**6)),
        "S2": sp.Integer(0),
        "S3": sp.Integer(0),
        "U1": to_nu_even(sp.Rational(105, 32) * (Xa**5 + Xb**5)),
        "U2": sp.Integer(0),
        "W1": to_nu_even(-sp.Rational(1, 16) * (Xa**4 + Xb**4)),
    }


def residual_slots(target: dict[str, sp.Expr], seed: dict[str, sp.Expr]) -> dict[str, list[sp.Expr]]:
    res = {k: sp.expand(target[k] - seed[k]) for k in target}
    return {
        "free": [res["K"]],
        "Q": [res["Q1"], res["Q2"], res["Q3"], res["Q4"], res["Q5"]],
        "T": [res["T1"], res["T2"], res["T3"], res["T4"]],
        "S": [res["S1"], res["S2"], res["S3"]],
        "U": [res["U1"], res["U2"]],
        "W": [res["W1"]],
    }


# ---------------------------------------------------------------------------
# Main theorem audit
# ---------------------------------------------------------------------------

def main() -> None:
    banner("PART I — GENERIC-FRAME LOCAL 4PN HAMILTONIAN SCAFFOLD")

    Qbasis = generate_basis(0, 8)
    Tbasis = generate_basis(1, 6)
    Sbasis = generate_basis(2, 4)
    Ubasis = generate_basis(3, 2)
    Wbasis = generate_basis(4, 0)

    print("Basis sizes:")
    print(f"  Q (G/r, degree 8)      : {len(Qbasis)}")
    print(f"  T (G^2/r^2, degree 6)  : {len(Tbasis)}")
    print(f"  S (G^3/r^3, degree 4)  : {len(Sbasis)}")
    print(f"  U (G^4/r^4, degree 2)  : {len(Ubasis)}")
    print(f"  W (G^5/r^5, static)    : {len(Wbasis)}")
    print(f"  total interaction dirs : {len(Qbasis)+len(Tbasis)+len(Sbasis)+len(Ubasis)+len(Wbasis)}")

    subbanner("I.1 — Hamiltonian-side formal invariant interpretation")
    print("The same formal invariant basis is used as in the ordinary scaffold, but with")
    print("  a = P_A^2 / m_A^2,  b = P_B^2 / m_B^2,  c = P_A.P_B / (m_A m_B),")
    print("  d = n.P_A / m_A,    e = n.P_B / m_B.")
    print("Because the quartic COM map showed that the local Hamiltonian chart matches the")
    print("constant-coefficient scaffold nu-degree ceilings blockwise, this is the correct")
    print("chart in which to attempt the local generic-frame lift.")

    banner("PART II — FIXED LOCAL COM HAMILTONIAN TARGET, NATURAL SEED, AND RESIDUAL")
    target = local_adm_target()
    seed = hamiltonian_seed()
    residual = residual_slots(target, seed)

    print("Seed slot coefficients:")
    for key in ["K", "Q1", "T1", "S1", "U1", "W1"]:
        print(f"  {key}^(seed) = {sp.factor(seed[key])}")

    subbanner("II.1 — Free slot check")
    expect_zero("free Hamiltonian residual", residual["free"][0])

    subbanner("II.2 — Interaction residual slots")
    for block, keys in [
        ("Q", ["Q1", "Q2", "Q3", "Q4", "Q5"]),
        ("T", ["T1", "T2", "T3", "T4"]),
        ("S", ["S1", "S2", "S3"]),
        ("U", ["U1", "U2"]),
        ("W", ["W1"]),
    ]:
        print(f"\n{block}-block residuals:")
        for key, expr in zip(keys, residual[block]):
            print(f"  {key}^(res) = {sp.factor(expr)}")

    banner("PART III — EXACT COEFFICIENT-SPACE IMAGE RANKS")
    blocks = {
        "Q": (Qbasis, 4),
        "T": (Tbasis, 3),
        "S": (Sbasis, 3),
        "U": (Ubasis, 2),
        "W": (Wbasis, 2),
    }

    image_data = {}
    total_rank = 0
    total_cols = 0
    total_rows = 0
    for block, (basis, maxdeg) in blocks.items():
        M = image_matrix_polynomial(block, basis, maxdeg)
        rank = M.rank()
        rows, cols = M.shape
        nullity = cols - rank
        image_data[block] = (M, rank, rows, cols, nullity)
        total_rank += rank
        total_cols += cols
        total_rows += rows
        print(f"{block}: shape {rows} x {cols}, rank = {rank}, nullity = {nullity}")

    print(f"\nTotal coefficient-space dimensions : {total_rows}")
    print(f"Total interaction directions       : {total_cols}")
    print(f"Total rank                         : {total_rank}")
    print(f"Total nullity                      : {total_cols - total_rank}")

    expected = {"Q": 20, "T": 12, "S": 9, "U": 4, "W": 2}
    for block, (_, rank, rows, _, _) in image_data.items():
        if rank != expected[block] or rows != expected[block]:
            raise AssertionError(f"{block}-block did not achieve full coefficient-space rank.")

    banner("PART IV — EXACT BLOCKWISE GENERIC-FRAME HAMILTONIAN LIFT")
    QB, Qmap = coordinate_matrix(Qbasis)
    TB, Tmap = coordinate_matrix(Tbasis)
    SB, Smap = coordinate_matrix(Sbasis)
    UB, Umap = coordinate_matrix(Ubasis)
    WB, Wmap = coordinate_matrix(Wbasis)

    mats = {"Q": QM if False else None}
    # Build polynomial coefficient matrices explicitly for the solves.
    QM = image_matrix_polynomial("Q", Qbasis, 4)
    TM = image_matrix_polynomial("T", Tbasis, 3)
    SM = image_matrix_polynomial("S", Sbasis, 3)
    UM = image_matrix_polynomial("U", Ubasis, 2)
    WM = image_matrix_polynomial("W", Wbasis, 2)

    Qvec = target_vector(residual["Q"], 4)
    Tvec = target_vector(residual["T"], 3)
    Svec = target_vector(residual["S"], 3)
    Uvec = target_vector(residual["U"], 2)
    Wvec = target_vector(residual["W"], 2)

    Qcoords = particular_solution(QM, Qvec)
    Tcoords = particular_solution(TM, Tvec)
    Scoords = particular_solution(SM, Svec)
    Ucoords = particular_solution(UM, Uvec)
    Wcoords = particular_solution(WM, Wvec)

    expect_zero("Q image solve", QM * Qcoords - Qvec)
    expect_zero("T image solve", TM * Tcoords - Tvec)
    expect_zero("S image solve", SM * Scoords - Svec)
    expect_zero("U image solve", UM * Ucoords - Uvec)
    expect_zero("W image solve", WM * Wcoords - Wvec)

    subbanner("IV.1 — One exact representative by block")
    for label, coords, basis in [
        ("Q", Qcoords, Qbasis),
        ("T", Tcoords, Tbasis),
        ("S", Scoords, Sbasis),
        ("U", Ucoords, Ubasis),
        ("W", Wcoords, Wbasis),
    ]:
        nz = nonzero_terms(coords, basis)
        print(f"\n{label}-block: {len(nz)} nonzero scaffold directions")
        for i, coef, expr in nz:
            print(f"  [{i:02d}] {coef}   *   {expr}")

    banner("PART V — DIRECT COM VERIFICATION OF THE REPRESENTATIVE")
    Qexpr = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Qcoords, Qbasis)))
    Texpr = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Tcoords, Tbasis)))
    Sexpr = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Scoords, Sbasis)))
    Uexpr = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Ucoords, Ubasis)))
    Wexpr = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Wcoords, Wbasis)))

    Qslots = [sp.expand(s) for s in block_slots(to_nu(Qexpr), "Q")]
    Tslots = [sp.expand(s) for s in block_slots(to_nu(Texpr), "T")]
    Sslots = [sp.expand(s) for s in block_slots(to_nu(Sexpr), "S")]
    Uslots = [sp.expand(s) for s in block_slots(to_nu(Uexpr), "U")]
    Wslots = [sp.expand(s) for s in block_slots(to_nu(Wexpr), "W")]

    for i, (lhs, rhs) in enumerate(zip(Qslots, residual["Q"]), start=1):
        expect_zero(f"Q slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(Tslots, residual["T"]), start=1):
        expect_zero(f"T slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(Sslots, residual["S"]), start=1):
        expect_zero(f"S slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(Uslots, residual["U"]), start=1):
        expect_zero(f"U slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(Wslots, residual["W"]), start=1):
        expect_zero(f"W slot {i}", lhs - rhs)

    banner("PART VI — FINAL HAMILTONIAN-CHART LIFT LEDGER")
    print("1. The natural Hamiltonian self/static seed already captures the entire free 4PN slot:")
    print("      Delta K = 0.")
    print("2. Every local interaction block of the fixed-chart COM Hamiltonian residual lies in the")
    print("   image of the constant-coefficient exchange-symmetric generic-frame scaffold.")
    print("3. Blockwise coefficient-space ranks are maximal:")
    print("      Q : 20/20,   T : 12/12,   S : 9/9,   U : 4/4,   W : 2/2.")
    print("4. Therefore there is no 4PN local Hamiltonian-chart analogue of the ordinary-chart")
    print("   nu^4-tail obstruction.")
    print("5. The remaining freedom is only the large generic-frame polynomial nullspace")
    print("      32 + 34 + 20 + 6 + 0 = 92 directions,")
    print("   not a failure of span.")
    print("6. So the fixed local 4PN program can proceed Hamiltonian-first:")
    print("      import/fix Hamiltonian chart  ->  solve generic-frame interaction lift  ->  translate")
    print("   back to the ordinary chart only afterward if desired.")


if __name__ == "__main__":
    main()
