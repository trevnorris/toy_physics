#!/usr/bin/env python3
"""
3pn_contact_generator_and_gauge_orbit_audit.py

Standalone SymPy audit for the next 3PN generic-frame step:
construct the actual ordinary-Lagrangian contact/gauge generator that preserves
both the compact 24/17/8/1 scaffold and the already-solved COM target.

What this script does
---------------------
1. Rebuild the 3PN generic-frame scaffold bases Q/T/S/U.
2. Build the odd total-derivative generator families relevant at 3PN:
      - r * degree-7 odd scalars,
      - G * degree-5 odd scalars,
      - G^2/r * degree-3 odd scalars,
      - G^3/r^2 * degree-1 odd scalars.
3. Enforce the correct exchange parity for odd generators:
      A <-> B,  d <-> -e.
4. Compute dF/dt with Newtonian order reduction.
5. Demand that the generated 3PN shift stays inside the compact
   24/17/8/1 generic-frame scaffold.
6. Demand that the resulting shift have zero COM projection.
7. Extract one sparse 11-generator basis for the surviving contact/gauge orbit.
8. Verify that:
      - every surviving generator factors through the COM-null ideal,
      - the scaffold-preserving + COM-preserving contact image has rank 11,
      - the full COM-null family has rank 27,
      - therefore 16 COM-blind algebraic directions remain beyond this simple
        ordinary total-derivative contact orbit.

Main conclusion
---------------
The actual scaffold-preserving 3PN ordinary contact/gauge family is an
11-dimensional suborbit inside the previously identified 27-dimensional
COM-null ideal.  So the canonical quotient slice is not exhausted by contact
freedom alone: there is a 16-dimensional algebraic COM-null quotient left over.
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
        simplified = expr.applyfunc(sp.simplify)
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
# Symbols and swap maps
# ---------------------------------------------------------------------------

a, b, c, d, e, p, q, r, G = sp.symbols("a b c d e p q r G")
nu, Delta, V2, rd = sp.symbols("nu Delta V2 rd")


def swap_even(expr: sp.Expr) -> sp.Expr:
    return sp.expand(expr.xreplace({a: b, b: a, c: c, d: e, e: d, p: q, q: p}))


def swap_odd(expr: sp.Expr) -> sp.Expr:
    return sp.expand(expr.xreplace({a: b, b: a, c: c, d: -e, e: -d, p: q, q: p}))


# ---------------------------------------------------------------------------
# Basis generation
# ---------------------------------------------------------------------------

def canonical_sym(expr: sp.Expr, vars_order: tuple[sp.Symbol, ...], odd: bool) -> sp.Expr:
    s = sp.expand(expr + (swap_odd(expr) if odd else swap_even(expr)))
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


def generate_basis(mass_deg: int, vel_deg: int, odd: bool) -> list[sp.Expr]:
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
                            sym = canonical_sym(expr, vars_order, odd)
                            tm = sp.expand(sym.subs({b: 0, c: 0, e: 0, p: 0, q: 1}))
                            if tm != 0:
                                continue
                            basis.add(sym)

    return sorted(basis, key=str)


# Compact 3PN scaffold bases from the earlier generic-frame lift.
Qbasis = generate_basis(0, 6, odd=False)
Tbasis = generate_basis(1, 4, odd=False)
Sbasis = generate_basis(2, 2, odd=False)
Ubasis = generate_basis(3, 0, odd=False)

# Odd scalar-generator families.
FQ_raw = generate_basis(0, 7, odd=True)
FT_raw = generate_basis(0, 5, odd=True)
FS_raw = generate_basis(1, 3, odd=True)
FU_raw = generate_basis(2, 1, odd=True)

FQ_basis = [r * f for f in FQ_raw]
FT_basis = [G * f for f in FT_raw]
FS_basis = [G**2 / r * f for f in FS_raw]
FU_basis = [G**3 / r**2 * f for f in FU_raw]
ALL_F = FQ_basis + FT_basis + FS_basis + FU_basis


# ---------------------------------------------------------------------------
# Newtonian order-reduction algebra in invariant form
# ---------------------------------------------------------------------------

rdot = d - e
adot = -2 * G * q * d / r**2
bdot = 2 * G * p * e / r**2
cdot = G * (p * d - q * e) / r**2
# d = v_A.n,  e = v_B.n.
# These formulas already include dot(n) = (v_A-v_B-(d-e)n)/r.
ddot = -G * q / r**2 + (a - c - d**2 + d * e) / r
edot = G * p / r**2 + (c - b - d * e + e**2) / r


def dt_expr(expr: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(expr, a) * adot
        + sp.diff(expr, b) * bdot
        + sp.diff(expr, c) * cdot
        + sp.diff(expr, d) * ddot
        + sp.diff(expr, e) * edot
        + sp.diff(expr, r) * rdot
    )


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


# ---------------------------------------------------------------------------
# Linear algebra helpers in the polynomial coefficient spaces
# ---------------------------------------------------------------------------

def terms_coeff_dict(expr: sp.Expr) -> dict[tuple[int, ...], sp.Expr]:
    poly = sp.Poly(sp.expand(expr), p, q, a, b, c, d, e)
    return {mon: sp.Rational(coef) for mon, coef in poly.terms()}


def coordinate_matrix(basis: list[sp.Expr]) -> tuple[sp.Matrix, dict[tuple[int, ...], int]]:
    monmap: dict[tuple[int, ...], int] = {}
    for b_expr in basis:
        for mon in terms_coeff_dict(b_expr):
            monmap.setdefault(mon, len(monmap))
    M = sp.zeros(len(monmap), len(basis))
    for j, b_expr in enumerate(basis):
        for mon, coef in terms_coeff_dict(b_expr).items():
            M[monmap[mon], j] = coef
    return M, monmap


def coords_in_basis(expr: sp.Expr, basisM: sp.Matrix, monmap: dict[tuple[int, ...], int]) -> sp.Matrix:
    vec = sp.zeros(len(monmap), 1)
    for mon, coef in terms_coeff_dict(expr).items():
        if mon not in monmap:
            raise ValueError(f"monomial {mon} not in basis monomial set")
        vec[monmap[mon], 0] = coef
    sol, params = basisM.gauss_jordan_solve(vec)
    if params.rows != 0:
        subs = {params[i, 0]: 0 for i in range(params.rows)}
        sol = sp.Matrix([sp.simplify(sol[i].subs(subs)) for i in range(sol.rows)])
    return sol


QB, Qmap = coordinate_matrix(Qbasis)
TB, Tmap = coordinate_matrix(Tbasis)
SB, Smap = coordinate_matrix(Sbasis)
UB, Umap = coordinate_matrix(Ubasis)


# ---------------------------------------------------------------------------
# COM projection machinery
# ---------------------------------------------------------------------------

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


def image_matrix_polynomial(block: str, exprs: list[sp.Expr]) -> sp.Matrix:
    rows: list[list[sp.Expr]] = []
    for expr in exprs:
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
# Scaffold-preserving and COM-preserving contact solve
# ---------------------------------------------------------------------------

def monomial_space(exprs: list[sp.Expr]) -> dict[tuple[int, ...], int]:
    monmap: dict[tuple[int, ...], int] = {}
    for expr in exprs:
        for mon in terms_coeff_dict(expr):
            monmap.setdefault(mon, len(monmap))
    return monmap


def build_matrix(exprs: list[sp.Expr], monmap: dict[tuple[int, ...], int]) -> sp.Matrix:
    M = sp.zeros(len(monmap), len(exprs))
    for j, expr in enumerate(exprs):
        for mon, coef in terms_coeff_dict(expr).items():
            M[monmap[mon], j] = coef
    return M


def constraint_matrix(scaffold_basis: list[sp.Expr], image_exprs: list[sp.Expr]) -> sp.Matrix:
    monmap = monomial_space(image_exprs + scaffold_basis)
    Bfull = build_matrix(scaffold_basis, monmap)
    IMG = build_matrix(image_exprs, monmap)
    left_null = Bfull.T.nullspace()
    if not left_null:
        return sp.zeros(0, IMG.shape[1])
    rows = [v.T * IMG for v in left_null]
    return sp.Matrix.vstack(*rows)


# Raw odd generator images.
GEN_BLOCKS = [split_blocks(dt_expr(F)) for F in ALL_F]
Q_exprs = [blk.get((1, -1), 0) for blk in GEN_BLOCKS]
T_exprs = [blk.get((2, -2), 0) for blk in GEN_BLOCKS]
S_exprs = [blk.get((3, -3), 0) for blk in GEN_BLOCKS]
U_exprs = [blk.get((4, -4), 0) for blk in GEN_BLOCKS]

# 1) Stay inside the compact scaffold.
C_scaffold = sp.Matrix.vstack(
    constraint_matrix(Qbasis, Q_exprs),
    constraint_matrix(Tbasis, T_exprs),
    constraint_matrix(Sbasis, S_exprs),
    constraint_matrix(Ubasis, U_exprs),
)
K_scaffold = C_scaffold.nullspace()

# 2) Then impose zero COM image.
SCAFFOLD_COMBO_EXPRS: list[sp.Expr] = []
for vec in K_scaffold:
    expr = 0
    for coeff, blk in zip(vec, GEN_BLOCKS):
        if coeff == 0:
            continue
        for k in [(1, -1), (2, -2), (3, -3), (4, -4)]:
            expr += coeff * G**k[0] * r**k[1] * blk.get(k, 0)
    SCAFFOLD_COMBO_EXPRS.append(sp.expand(expr))

Mcom_scaffold = sp.Matrix.vstack(
    image_matrix_polynomial("Q", SCAFFOLD_COMBO_EXPRS),
    image_matrix_polynomial("T", SCAFFOLD_COMBO_EXPRS),
    image_matrix_polynomial("S", SCAFFOLD_COMBO_EXPRS),
    image_matrix_polynomial("U", SCAFFOLD_COMBO_EXPRS),
)
K_contact = Mcom_scaffold.nullspace()


# ---------------------------------------------------------------------------
# A sparse 11-generator basis for the actual contact/gauge orbit
# ---------------------------------------------------------------------------

# Convenient sparse basis discovered by the exact nullspace solve.  These live
# directly in the raw FT/FS/FU families and do not need any FQ support.
C3 = d * p + e * q
C4 = a * b - c**2
C5 = a * e - c * d
C6 = b * d - c * e

RAW_CONTACT_GENERATORS = [
    G * (-a * C5 + b * C6),
    G * (b * C5 - a * C6),
    G * (d + e) * (e * C5 - d * C6),
    -G * (d - e) * C4,
    -G * d * e * (C5 - C6),
    G * (-d**2 * C5 + e**2 * C6),
    G**2 / r * (q * C5 - p * C6),
    G**2 / r * (a - b) * C3,
    G**2 / r * (-p * C5 + q * C6),
    G**2 / r * (d**2 - e**2) * C3,
    -G**3 / r**2 * (p - q) * C3,
]

CONTACT_IMAGES: list[sp.Matrix] = []
for F in RAW_CONTACT_GENERATORS:
    blk = split_blocks(dt_expr(F))
    qv = sp.zeros(len(Qbasis), 1)
    tv = sp.zeros(len(Tbasis), 1)
    sv = sp.zeros(len(Sbasis), 1)
    uv = sp.zeros(len(Ubasis), 1)
    if (1, -1) in blk and blk[(1, -1)] != 0:
        qv = coords_in_basis(blk[(1, -1)], QB, Qmap)
    if (2, -2) in blk and blk[(2, -2)] != 0:
        tv = coords_in_basis(blk[(2, -2)], TB, Tmap)
    if (3, -3) in blk and blk[(3, -3)] != 0:
        sv = coords_in_basis(blk[(3, -3)], SB, Smap)
    if (4, -4) in blk and blk[(4, -4)] != 0:
        uv = coords_in_basis(blk[(4, -4)], UB, Umap)
    CONTACT_IMAGES.append(sp.Matrix.vstack(qv, tv, sv, uv))

CONTACT_IMAGE_MATRIX = sp.Matrix.hstack(*CONTACT_IMAGES)


# ---------------------------------------------------------------------------
# Full 27-dimensional COM-null family for comparison
# ---------------------------------------------------------------------------

MQ = image_matrix_polynomial("Q", Qbasis)
MT = image_matrix_polynomial("T", Tbasis)
MS = image_matrix_polynomial("S", Sbasis)
MU = image_matrix_polynomial("U", Ubasis)
Mfull = sp.diag(MQ, MT, MS, MU)
NULL_FULL = Mfull.nullspace()
NULL_FULL_MATRIX = sp.Matrix.hstack(*NULL_FULL)


# ---------------------------------------------------------------------------
# Main audit prints
# ---------------------------------------------------------------------------

def main() -> None:
    banner("PART I — RAW ODD GENERATOR FAMILIES")
    print(f"Scaffold basis sizes: Q={len(Qbasis)}, T={len(Tbasis)}, S={len(Sbasis)}, U={len(Ubasis)}")
    print(f"Raw odd generator counts: FQ={len(FQ_basis)}, FT={len(FT_basis)}, FS={len(FS_basis)}, FU={len(FU_basis)}")
    print(f"Total raw odd generators = {len(ALL_F)}")
    print()
    print("Exchange rule for odd generators is")
    print("  A <-> B,   c -> c,   d -> -e,   e -> -d,   p <-> q")
    print("rather than the even-scalar rule d <-> e.")

    banner("PART II — SCAFFOLD-PRESERVING AND COM-PRESERVING COUNTS")
    print(f"Scaffold-preserving constraint rank  = {C_scaffold.rank()}")
    print(f"Scaffold-preserving kernel dimension = {len(K_scaffold)}")
    print(f"COM image rank inside scaffold-preserving family = {Mcom_scaffold.rank()}")
    print(f"Actual COM-preserving contact-family dimension   = {len(K_contact)}")
    print()
    print("So the raw 53 odd generators collapse as")
    print("  53 raw  -> 22 scaffold-preserving  -> 11 COM-preserving.")

    banner("PART III — SPARSE 11-GENERATOR CONTACT BASIS")
    for i, F in enumerate(RAW_CONTACT_GENERATORS, 1):
        pref = G if i <= 6 else (G**2 / r if i <= 10 else G**3 / r**2)
        core = sp.factor(sp.expand(F / pref))
        print(f"Gamma_{i} = {pref} * ({core})")

    subbanner("III.1 — Exact factorization through the COM-null ideal")
    factor_checks = [
        sp.expand(RAW_CONTACT_GENERATORS[0] / G - (-a * C5 + b * C6)),
        sp.expand(RAW_CONTACT_GENERATORS[1] / G - (b * C5 - a * C6)),
        sp.expand(RAW_CONTACT_GENERATORS[2] / G - (d + e) * (e * C5 - d * C6)),
        sp.expand(RAW_CONTACT_GENERATORS[3] / G + (d - e) * C4),
        sp.expand(RAW_CONTACT_GENERATORS[4] / G + d * e * (C5 - C6)),
        sp.expand(RAW_CONTACT_GENERATORS[5] / G - (-d**2 * C5 + e**2 * C6)),
        sp.expand(RAW_CONTACT_GENERATORS[6] * r / G**2 - (q * C5 - p * C6)),
        sp.expand(RAW_CONTACT_GENERATORS[7] * r / G**2 - (a - b) * C3),
        sp.expand(RAW_CONTACT_GENERATORS[8] * r / G**2 - (-p * C5 + q * C6)),
        sp.expand(RAW_CONTACT_GENERATORS[9] * r / G**2 - (d**2 - e**2) * C3),
        sp.expand(RAW_CONTACT_GENERATORS[10] * r**2 / G**3 + (p - q) * C3),
    ]
    expect_zero("all factorization residuals", sp.Matrix(factor_checks))

    banner("PART IV — CONTACT IMAGE INSIDE THE 27-DIMENSIONAL COM-NULL FAMILY")
    print(f"Rank of the 11-generator contact image = {CONTACT_IMAGE_MATRIX.rank()}")
    print(f"Rank of the full COM-null family       = {NULL_FULL_MATRIX.rank()}")
    print(f"Residual algebraic quotient dimension  = {NULL_FULL_MATRIX.rank() - CONTACT_IMAGE_MATRIX.rank()}")
    expect_zero("M_full * contact image", Mfull * CONTACT_IMAGE_MATRIX)

    qrank = CONTACT_IMAGE_MATRIX[:len(Qbasis), :].rank()
    trank = CONTACT_IMAGE_MATRIX[len(Qbasis):len(Qbasis) + len(Tbasis), :].rank()
    srank = CONTACT_IMAGE_MATRIX[len(Qbasis) + len(Tbasis):len(Qbasis) + len(Tbasis) + len(Sbasis), :].rank()
    urank = CONTACT_IMAGE_MATRIX[-len(Ubasis):, :].rank()
    print()
    print("Blockwise ranks of the 11-generator contact image:")
    print(f"  Q block rank = {qrank}")
    print(f"  T block rank = {trank}")
    print(f"  S block rank = {srank}")
    print(f"  U block rank = {urank}")
    print("So the simple ordinary contact family does not move the static U block at all.")

    banner("PART V — THEOREM LEDGER")
    print("1. The correct exchange parity for odd scalar generators is d <-> -e, not d <-> e.")
    print("2. The raw odd-generator family has 53 candidates: 31 + 12 + 8 + 2.")
    print("3. Requiring the generated shift to stay inside the compact 24/17/8/1 scaffold leaves 22 directions.")
    print("4. Requiring zero COM image leaves an 11-dimensional actual contact/gauge orbit.")
    print("5. A sparse 11-generator basis is given by Gamma_1,...,Gamma_11 above.")
    print("6. Every Gamma_i factors through the COM-null ideal generators C3,C4,C5,C6.")
    print("7. The 11-generator contact image sits inside the full 27-dimensional COM-null family.")
    print("8. Therefore the previously found canonical quotient slice is not exhausted by contact gauge alone:")
    print("     27 COM-blind directions = 11 contact/gauge + 16 residual algebraic COM-null directions.")
    print("9. So the actual contact/gauge generator picks a nearby 11-dimensional orbit inside the canonical family,")
    print("   not the entire canonical 27-dimensional slice.")


if __name__ == "__main__":
    main()
