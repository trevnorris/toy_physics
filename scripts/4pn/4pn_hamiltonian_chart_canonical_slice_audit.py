#!/usr/bin/env python3
"""
4pn_hamiltonian_chart_canonical_slice_audit.py

Canonical generic-frame Hamiltonian slice for the local 4PN lift.

What this script does
---------------------
1. Rebuilds the exchange-symmetric local 4PN Hamiltonian-chart generic-frame scaffold:
      - G/r     times degree-8 momentum invariants,
      - G^2/r^2 times degree-6 invariants with one extra mass power,
      - G^3/r^3 times degree-4 invariants with two extra mass powers,
      - G^4/r^4 times degree-2 invariants with three extra mass powers,
      - G^5/r^5 static cross polynomials.
2. Re-imports the fixed local COM Hamiltonian target and the natural Hamiltonian
   self/static seed, then reconstructs the exact comparable-mass interaction residual.
3. Computes the exact blockwise coefficient-space nullspaces:
      Q : 32,  T : 34,  S : 20,  U : 6,  W : 0.
4. Proves that those null directions are precisely COM-blind algebraic directions:
      - all Q-null basis vectors lie in the full COM ideal
            <pa+qc, pc+qb, pd+qe, ab-c^2, ae-cd, bd-ce>,
      - all T/S/U-null basis vectors already lie in the smaller linear-momentum ideal
            <pa+qc, pc+qb, pd+qe>.
5. Fixes a canonical quotient slice by setting all null coordinates to zero in the
   exact Gauss-Jordan lift.
6. Prints one explicit canonical generic-frame Hamiltonian representative
   (Q_can, T_can, S_can, U_can, W_can) and verifies that its COM reduction
   reproduces the imported fixed local COM Hamiltonian residual exactly.

Interpretation
--------------
This is the 4PN local analogue of the 3PN COM-null/canonical-slice step.

After the earlier existence theorem, the only local Hamiltonian ambiguity left was the
92-direction generic-frame null family. This audit shows that the family is purely
algebraic/COM-blind, identifies the corresponding ideals, and freezes a sparse
canonical representative in the chosen Hamiltonian chart.
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


def nullspace_dimension(M: sp.Matrix) -> int:
    return len(M.nullspace())


# ---------------------------------------------------------------------------
# Basis generation / COM projection machinery
# ---------------------------------------------------------------------------

a, b, c, d, e, p, q = sp.symbols("a b c d e p q")
nu, Delta = sp.symbols("nu Delta", real=True)
V2, rd = sp.symbols("V2 rd", real=True)
P = sp.Symbol("P")

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

                            # Strict one-body branch: A test mass around fixed B.
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


def null_polynomials(M: sp.Matrix, basis: list[sp.Expr]) -> list[sp.Expr]:
    polys: list[sp.Expr] = []
    for vec in M.nullspace():
        expr = sp.expand(sum(sp.expand(vec[i, 0]) * basis[i] for i in range(vec.rows)))
        polys.append(expr)
    return polys


def ideal_membership_counts(polys: list[sp.Expr], G: sp.GroebnerBasis) -> int:
    count = 0
    for expr in polys:
        rem = G.reduce(expr)[1]
        if sp.expand(rem) == 0:
            count += 1
    return count


# ---------------------------------------------------------------------------
# Imported fixed local COM Hamiltonian target and natural seed
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
    banner("PART I — BLOCKWISE NULLITIES OF THE LOCAL 4PN HAMILTONIAN-CHART LIFT")
    Qbasis = generate_basis(0, 8)
    Tbasis = generate_basis(1, 6)
    Sbasis = generate_basis(2, 4)
    Ubasis = generate_basis(3, 2)
    Wbasis = generate_basis(4, 0)

    QM = image_matrix_polynomial("Q", Qbasis, 4)
    TM = image_matrix_polynomial("T", Tbasis, 3)
    SM = image_matrix_polynomial("S", Sbasis, 3)
    UM = image_matrix_polynomial("U", Ubasis, 2)
    WM = image_matrix_polynomial("W", Wbasis, 2)

    print(f"Q nullity = {nullspace_dimension(QM)}")
    print(f"T nullity = {nullspace_dimension(TM)}")
    print(f"S nullity = {nullspace_dimension(SM)}")
    print(f"U nullity = {nullspace_dimension(UM)}")
    print(f"W nullity = {nullspace_dimension(WM)}")
    print("Total generic-frame interaction nullity =", sum(map(nullspace_dimension, [QM, TM, SM, UM, WM])))

    banner("PART II — THE COM-NULL IDEALS AT 4PN")
    C1 = p * a + q * c
    C2 = p * c + q * b
    C3 = p * d + q * e
    C4 = a * b - c**2
    C5 = a * e - c * d
    C6 = b * d - c * e

    print("Linear-momentum ideal generators:")
    print("  C1 = p a + q c")
    print("  C2 = p c + q b")
    print("  C3 = p d + q e")
    print("Full COM ideal adds:")
    print("  C4 = a b - c^2")
    print("  C5 = a e - c d")
    print("  C6 = b d - c e")

    Glin = sp.groebner([C1, C2, C3], a, b, c, d, e, p, q, order="grevlex", domain=sp.QQ)
    Gfull = sp.groebner([C1, C2, C3, C4, C5, C6], a, b, c, d, e, p, q, order="grevlex", domain=sp.QQ)

    Qnull = null_polynomials(QM, Qbasis)
    Tnull = null_polynomials(TM, Tbasis)
    Snull = null_polynomials(SM, Sbasis)
    Unull = null_polynomials(UM, Ubasis)

    q_full = ideal_membership_counts(Qnull, Gfull)
    q_lin = ideal_membership_counts(Qnull, Glin)
    t_full = ideal_membership_counts(Tnull, Gfull)
    t_lin = ideal_membership_counts(Tnull, Glin)
    s_full = ideal_membership_counts(Snull, Gfull)
    s_lin = ideal_membership_counts(Snull, Glin)
    u_full = ideal_membership_counts(Unull, Gfull)
    u_lin = ideal_membership_counts(Unull, Glin)

    print(f"Q-null basis vectors in full COM ideal     = {q_full} / {len(Qnull)}")
    print(f"Q-null basis vectors in linear ideal       = {q_lin} / {len(Qnull)}")
    print(f"T-null basis vectors in full COM ideal     = {t_full} / {len(Tnull)}")
    print(f"T-null basis vectors in linear ideal       = {t_lin} / {len(Tnull)}")
    print(f"S-null basis vectors in full COM ideal     = {s_full} / {len(Snull)}")
    print(f"S-null basis vectors in linear ideal       = {s_lin} / {len(Snull)}")
    print(f"U-null basis vectors in full COM ideal     = {u_full} / {len(Unull)}")
    print(f"U-null basis vectors in linear ideal       = {u_lin} / {len(Unull)}")

    if q_full != len(Qnull) or t_full != len(Tnull) or s_full != len(Snull) or u_full != len(Unull):
        raise AssertionError("Some null vectors failed the full COM-ideal test.")
    if q_lin != 0 or t_lin != len(Tnull) or s_lin != len(Snull) or u_lin != len(Unull):
        raise AssertionError("Unexpected linear-ideal classification of null blocks.")

    banner("PART III — CANONICAL QUOTIENT SLICE")
    target = local_adm_target()
    seed = hamiltonian_seed()
    residual = residual_slots(target, seed)

    expect_zero("free Hamiltonian residual", residual["free"][0])

    Qcoords = particular_solution(QM, target_vector(residual["Q"], 4))
    Tcoords = particular_solution(TM, target_vector(residual["T"], 3))
    Scoords = particular_solution(SM, target_vector(residual["S"], 3))
    Ucoords = particular_solution(UM, target_vector(residual["U"], 2))
    Wcoords = particular_solution(WM, target_vector(residual["W"], 2))

    print("Canonical slice rule: set every null coordinate to zero in the exact Gauss-Jordan lift.")
    print("Nonzero scaffold directions in the resulting representative:")
    print(f"  Q : {len(nonzero_terms(Qcoords, Qbasis))}")
    print(f"  T : {len(nonzero_terms(Tcoords, Tbasis))}")
    print(f"  S : {len(nonzero_terms(Scoords, Sbasis))}")
    print(f"  U : {len(nonzero_terms(Ucoords, Ubasis))}")
    print(f"  W : {len(nonzero_terms(Wcoords, Wbasis))}")

    Qcan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Qcoords, Qbasis)))
    Tcan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Tcoords, Tbasis)))
    Scan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Scoords, Sbasis)))
    Ucan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Ucoords, Ubasis)))
    Wcan = sp.expand(sum(sp.expand(coef) * expr for coef, expr in zip(Wcoords, Wbasis)))

    subbanner("III.1 — Canonical explicit blocks")
    print("Q_can =")
    sp.pprint(Qcan.subs(P, sp.pi**2))
    print("\nT_can =")
    sp.pprint(Tcan.subs(P, sp.pi**2))
    print("\nS_can =")
    sp.pprint(Scan.subs(P, sp.pi**2))
    print("\nU_can =")
    sp.pprint(Ucan.subs(P, sp.pi**2))
    print("\nW_can =")
    sp.pprint(Wcan.subs(P, sp.pi**2))

    banner("PART IV — DIRECT COM VERIFICATION OF THE CANONICAL REPRESENTATIVE")
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(Qcan), "Q"), residual["Q"]), start=1):
        expect_zero(f"Q slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(Tcan), "T"), residual["T"]), start=1):
        expect_zero(f"T slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(Scan), "S"), residual["S"]), start=1):
        expect_zero(f"S slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(Ucan), "U"), residual["U"]), start=1):
        expect_zero(f"U slot {i}", lhs - rhs)
    for i, (lhs, rhs) in enumerate(zip(block_slots(to_nu(Wcan), "W"), residual["W"]), start=1):
        expect_zero(f"W slot {i}", lhs - rhs)

    banner("PART V — FINAL 4PN LOCAL CANONICAL-SLICE LEDGER")
    print("1. The 92-direction generic-frame ambiguity is purely COM-blind algebraic freedom.")
    print("2. Q-null directions require the full COM ideal <C1,...,C6>.")
    print("3. T/S/U-null directions already lie in the smaller linear-momentum ideal <C1,C2,C3>.")
    print("4. The W block has zero nullity, so the top local static block is already fixed once")
    print("   the COM target is fixed.")
    print("5. Setting all null coordinates to zero defines a canonical Hamiltonian-chart quotient slice.")
    print("6. The explicit block tuple (Q_can,T_can,S_can,U_can,W_can) is an exact generic-frame")
    print("   representative of the fixed local 4PN comparable-mass Hamiltonian residual.")
    print("7. This does not yet prove fixed-chart uniqueness against the true generic-frame ADM")
    print("   target, but it freezes a sparse exact canonical representative for the local-first program.")


if __name__ == "__main__":
    main()
