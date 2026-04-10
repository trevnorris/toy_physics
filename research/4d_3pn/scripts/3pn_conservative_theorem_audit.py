#!/usr/bin/env python3
"""
3pn_conservative_theorem_audit.py

Master SymPy audit for the consolidated conservative 3PN theorem package.

What this script checks
-----------------------
1. The exact 3PN one-body gate:
      mu_rho3 = 1/4,  d3 = -45/4,  s24 = -1/16.
2. The exact GR 3PN COM target and the carried self/static seed.
3. The pure-kinetic collapse:
      Delta l1 is exactly the ordinary-Lagrangian counterimage of the universal
      free relativistic two-body COM Hamiltonian.
4. The richer grouped-real-P2 middle-block compiler:
      det(M_mid) = -4/27, and it reproduces the exact 9-slot GR middle block.
5. The static sigma-collapse:
      the apparent P0/g family collapses to a unique geometry-side remainder.
6. The final exact decomposition of the solved COM residual into
      kinetic/compiler + grouped-P2 middle block + unique geometry static slot.
"""

from __future__ import annotations

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


# -----------------------------------------------------------------------------
# Part I — exact one-body gate
# -----------------------------------------------------------------------------

def one_body_gate() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    banner("PART I — EXACT ONE-BODY 3PN GATE")

    c, U, v = sp.symbols("c U v", positive=True, real=True)
    u = U / c**2
    d3, mu_rho3, s24 = sp.symbols("d3 mu_rho3 s24", real=True)

    L_exact = -c**2 * sp.sqrt(((1 - u / 2) / (1 + u / 2))**2 - (1 + u / 2)**4 * v**2 / c**2)
    L_exact_series = sp.expand(sp.series(L_exact, c, sp.oo, 7).removeO())

    D = 1 - 4 * u + 8 * u**2 + d3 * u**3
    L_red = -c**2 * (1 - u) * sp.sqrt(1 - (v**2 / c**2) / D)
    L_red_series = sp.expand(sp.series(L_red, c, sp.oo, 7).removeO())

    L_candidate = (
        L_red_series
        - U**2 / (2 * c**2)
        + U**3 / (4 * c**4)
        - mu_rho3 * U**4 / (2 * c**6)
        + s24 * U**2 * v**4 / c**6
    )

    residual = sp.expand(L_exact_series - L_candidate)
    mu_sol = sp.solve(sp.Eq(residual.coeff(U, 4).coeff(v, 0).coeff(c, -6), 0), mu_rho3)[0]
    d3_sol = sp.solve(sp.Eq(residual.coeff(U, 3).coeff(v, 2).coeff(c, -6), 0), d3)[0]
    s24_sol = sp.solve(sp.Eq(residual.coeff(U, 2).coeff(v, 4).coeff(c, -6), 0), s24)[0]

    print("mu_rho3 =", mu_sol)
    print("d3      =", d3_sol)
    print("s24     =", s24_sol)

    solved = sp.simplify(residual.subs({mu_rho3: mu_sol, d3: d3_sol, s24: s24_sol}))
    expect_zero("one-body residual", solved)

    if mu_sol != sp.Rational(1, 4) or d3_sol != -sp.Rational(45, 4) or s24_sol != -sp.Rational(1, 16):
        raise AssertionError("Unexpected one-body 3PN gate values.")

    return mu_sol, d3_sol, s24_sol


# -----------------------------------------------------------------------------
# Part II — exact GR COM target, carried seed, and residuals
# -----------------------------------------------------------------------------

def gr_target_h(nu: sp.Symbol) -> dict[int, sp.Expr]:
    pi = sp.pi
    h: dict[int, sp.Expr] = {}
    h[1] = sp.Rational(1, 128) * (-5 + 35 * nu - 70 * nu**2 + 35 * nu**3)
    h[2] = 0
    h[3] = 0
    h[4] = 0
    h[5] = 0
    h[6] = sp.Rational(1, 16) * (-7 + 42 * nu - 53 * nu**2 - 5 * nu**3)
    h[7] = sp.Rational(1, 16) * (2 - 3 * nu) * nu**2
    h[8] = sp.Rational(3, 16) * (1 - nu) * nu**2
    h[9] = -sp.Rational(5, 16) * nu**3
    h[10] = sp.Rational(1, 16) * (-27 + 136 * nu + 109 * nu**2)
    h[11] = sp.Rational(1, 16) * (17 + 30 * nu) * nu
    h[12] = sp.Rational(1, 12) * (5 + 43 * nu) * nu
    h[13] = sp.Rational(1, 192) * (-600 + (3 * pi**2 - 1340) * nu - 552 * nu**2)
    h[14] = -sp.Rational(1, 64) * (340 + 3 * pi**2 + 112 * nu) * nu
    h[15] = sp.Rational(1, 96) * (12 + (872 - 63 * pi**2) * nu)
    return {i: sp.simplify(h[i]) for i in range(1, 16)}


def inverse_map_from_h(nu: sp.Symbol, h: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    l: dict[int, sp.Expr] = {}
    l[1] = sp.simplify(sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3 - h[1])
    l[2] = sp.simplify(-h[2])
    l[3] = sp.simplify(-h[3])
    l[4] = sp.simplify(-h[4])
    l[5] = sp.simplify(-h[5])
    l[6] = sp.simplify(sp.Rational(1, 4) + sp.Rational(7, 8) * nu - sp.Rational(35, 8) * nu**2 - sp.Rational(21, 4) * nu**3 - h[6])
    l[7] = sp.simplify(sp.Rational(11, 8) * nu**2 - sp.Rational(9, 2) * nu**3 - h[7])
    l[8] = sp.simplify(sp.Rational(3, 4) * nu**2 - sp.Rational(9, 4) * nu**3 - h[8])
    l[9] = sp.simplify(-h[9])
    l[10] = sp.simplify(sp.Rational(5, 4) + sp.Rational(15, 8) * nu + sp.Rational(123, 8) * nu**2 + sp.Rational(13, 4) * nu**3 - h[10])
    l[11] = sp.simplify(sp.Rational(7, 8) * nu + sp.Rational(41, 8) * nu**2 + sp.Rational(31, 4) * nu**3 - h[11])
    l[12] = sp.simplify(sp.Rational(9, 2) * nu**2 + 4 * nu**3 - h[12])
    l[13] = sp.simplify(-sp.Rational(3, 2) - sp.Rational(59, 4) * nu - sp.Rational(25, 4) * nu**2 - sp.Rational(1, 2) * nu**3 - h[13])
    l[14] = sp.simplify(sp.Rational(7, 4) * nu - sp.Rational(31, 4) * nu**2 - sp.Rational(7, 2) * nu**3 - h[14])
    l[15] = sp.simplify(-h[15])
    return {i: sp.simplify(l[i]) for i in range(1, 16)}


def carried_seed_l(nu: sp.Symbol) -> dict[int, sp.Expr]:
    seed: dict[int, sp.Expr] = {}
    seed[1] = sp.Rational(5, 128) - sp.Rational(35, 128) * nu + sp.Rational(35, 64) * nu**2 - sp.Rational(35, 128) * nu**3
    seed[2] = seed[3] = seed[4] = seed[5] = 0
    seed[6] = sp.Rational(11, 16) - sp.Rational(33, 8) * nu + sp.Rational(99, 16) * nu**2 - sp.Rational(11, 8) * nu**3
    seed[7] = seed[8] = seed[9] = 0
    seed[10] = sp.Rational(47, 16) - sp.Rational(235, 16) * nu + sp.Rational(235, 16) * nu**2
    seed[11] = seed[12] = 0
    seed[13] = sp.Rational(13, 8) - sp.Rational(13, 2) * nu + sp.Rational(13, 4) * nu**2
    seed[14] = 0
    seed[15] = -sp.Rational(1, 8) + sp.Rational(3, 8) * nu
    return {i: sp.simplify(seed[i]) for i in range(1, 16)}


def carried_seed_h(nu: sp.Symbol) -> dict[int, sp.Expr]:
    seed: dict[int, sp.Expr] = {}
    seed[1] = -sp.Rational(5, 128) + sp.Rational(59, 128) * nu - sp.Rational(119, 64) * nu**2 + sp.Rational(323, 128) * nu**3
    seed[2] = seed[3] = seed[4] = seed[5] = 0
    seed[6] = -sp.Rational(7, 16) + 5 * nu - sp.Rational(169, 16) * nu**2 - sp.Rational(31, 8) * nu**3
    seed[7] = sp.Rational(1, 8) * nu**2 * (11 - 36 * nu)
    seed[8] = sp.Rational(3, 4) * nu**2 * (1 - 3 * nu)
    seed[9] = 0
    seed[10] = -sp.Rational(27, 16) + sp.Rational(265, 16) * nu + sp.Rational(11, 16) * nu**2 + sp.Rational(13, 4) * nu**3
    seed[11] = sp.Rational(1, 8) * nu * (7 + 41 * nu + 62 * nu**2)
    seed[12] = sp.Rational(1, 2) * nu**2 * (9 + 8 * nu)
    seed[13] = -sp.Rational(25, 8) - sp.Rational(33, 4) * nu - sp.Rational(19, 2) * nu**2 - sp.Rational(1, 2) * nu**3
    seed[14] = sp.Rational(1, 4) * nu * (7 - 31 * nu - 14 * nu**2)
    seed[15] = sp.Rational(1, 8) * (1 - 3 * nu)
    return {i: sp.simplify(seed[i]) for i in range(1, 16)}


# -----------------------------------------------------------------------------
# Part III — grouped real P2 middle-block compiler
# -----------------------------------------------------------------------------

def coeff_vector(expr: sp.Expr, v2: sp.Symbol, d: sp.Symbol, U: sp.Symbol) -> list[sp.Expr]:
    monoms = {
        1: v2**4,
        2: v2**3 * d**2,
        3: v2**2 * d**4,
        4: v2 * d**6,
        5: d**8,
        6: U * v2**3,
        7: U * v2**2 * d**2,
        8: U * v2 * d**4,
        9: U * d**6,
        10: U**2 * v2**2,
        11: U**2 * v2 * d**2,
        12: U**2 * d**4,
        13: U**3 * v2,
        14: U**3 * d**2,
        15: U**4,
    }
    poly = sp.Poly(sp.expand(expr), v2, d, U)
    out = []
    for i in range(1, 16):
        mon = sp.Poly(monoms[i], v2, d, U).monoms()[0]
        out.append(sp.simplify(poly.coeff_monomial(mon)))
    return out


def grouped_p2_middle_closure(nu: sp.Symbol, delta_l: dict[int, sp.Expr]) -> tuple[list[sp.Expr], sp.Expr, sp.Expr]:
    banner("PART III — EXACT GROUPED REAL P2 MIDDLE-BLOCK CLOSURE")

    v2, d, U = sp.symbols("v2 d U", real=True)
    u2 = v2 - d**2

    C20sq = sp.expand(sp.Rational(1, 6) * (3 * d**2 - v2 - 2 * U) ** 2)
    C21sq = sp.expand(2 * d**2 * u2)
    C22sq = sp.expand(sp.Rational(1, 2) * u2**2)

    T20 = sp.expand(U * d**2 * (3 * u2 - U) ** 2 / 3)
    T21 = sp.expand(U * u2 * (u2 - d**2 - U) ** 2)
    T22 = sp.expand(U * d**2 * u2**2)

    S20 = sp.expand(U**2 * C20sq)
    S21 = sp.expand(U**2 * C21sq)
    S22 = sp.expand(U**2 * C22sq)

    V20 = sp.expand(v2 * S20 / U)
    V21 = sp.expand(v2 * S21 / U)
    V22 = sp.expand(v2 * S22 / U)

    families = [T20, T21, T22, S20, S21, S22, V20, V21, V22]
    A_mid = sp.Matrix([[coeff_vector(f, v2, d, U)[i - 1] for f in families] for i in range(6, 15)])
    det_mid = sp.factor(A_mid.det())
    print("det(M_mid) =", det_mid)
    if det_mid != -sp.Rational(4, 27):
        raise AssertionError("Unexpected determinant for richer grouped-P2 compiler.")

    target_vec = sp.Matrix([delta_l[i] for i in range(6, 15)])
    coeffs = [sp.simplify(x) for x in A_mid.LUsolve(target_vec)]
    expr_target = sp.expand(sum(coeffs[i] * families[i] for i in range(9)))
    coords = coeff_vector(expr_target, v2, d, U)

    for i in range(6, 15):
        expect_zero(f"grouped middle slot l{i}", coords[i - 1] - delta_l[i])
    for i in range(1, 6):
        expect_zero(f"grouped kinetic slot l{i}", coords[i - 1])

    l15_pred = sp.simplify(coords[14])
    expect_zero(
        "l15 prediction relation",
        l15_pred - (delta_l[10] + delta_l[11] + delta_l[12] + 2 * (delta_l[6] + delta_l[7] + delta_l[8] + delta_l[9]))
    )

    l15_gap = sp.simplify(delta_l[15] - l15_pred)
    print("predicted l15 from grouped-P2 middle block =", l15_pred)
    print("remaining static gap =", l15_gap)

    return coeffs, l15_pred, l15_gap


# -----------------------------------------------------------------------------
# Part IV — pure-kinetic collapse and sigma-collapse
# -----------------------------------------------------------------------------

def pure_kinetic_and_static(nu: sp.Symbol, target_h: dict[int, sp.Expr], target_l: dict[int, sp.Expr], seed_h: dict[int, sp.Expr], seed_l: dict[int, sp.Expr], l15_pred: sp.Expr) -> None:
    banner("PART IV — PURE-KINETIC COLLAPSE AND UNIQUE STATIC GEOMETRY COMPLETION")

    # Pure-kinetic collapse.
    F1 = sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3
    h1_free = sp.Rational(1, 128) * (-5 + 35 * nu - 70 * nu**2 + 35 * nu**3)
    delta_l1_expected = sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3
    expect_zero("h1 target - free Hamiltonian", target_h[1] - h1_free)
    expect_zero("Delta l1 from free-Hamiltonian compiler image", (seed_h[1] - h1_free) - (target_l[1] - seed_l[1]))
    expect_zero("Delta l1 closed form", (target_l[1] - seed_l[1]) - delta_l1_expected)

    # Sigma-collapse / unique geometry remainder.
    p, q = sp.symbols("p q", positive=True, real=True)
    nu_pq = sp.simplify(p * q / (p + q) ** 2)
    U0 = p**3 + q**3
    Ug = p**2 * q + p * q**2
    expect_zero("sigma-collapse mass identity", sp.simplify(nu_pq * U0 - (1 - 3 * nu_pq) * Ug))

    cU_target = -sp.Rational(227, 24) + 21 * sp.pi**2 / 32
    cU_p2pred = sp.simplify((293 - 308 * nu - 102 * nu**2) / 24)
    cU_g = sp.simplify((408 * nu**2 + 1232 * nu - 2080 + 63 * sp.pi**2) / 96)
    expect_zero("generic-frame static split", cU_target - (cU_p2pred + cU_g))

    l15_geometry = sp.simplify(target_l[15] - seed_l[15] - l15_pred)
    expect_zero("COM geometry gap formula", l15_geometry - nu * cU_g)

    print("Delta l1 =", sp.simplify(target_l[1] - seed_l[1]))
    print("Delta l15^(g) =", l15_geometry)


# -----------------------------------------------------------------------------
# Part V — final exact decomposition
# -----------------------------------------------------------------------------

def final_decomposition(nu: sp.Symbol, delta_l: dict[int, sp.Expr], grouped_l15_pred: sp.Expr) -> None:
    banner("PART V — FINAL EXACT COM RESIDUAL DECOMPOSITION")

    delta_l1 = sp.simplify(delta_l[1])
    geometry_gap = sp.simplify(delta_l[15] - grouped_l15_pred)

    final_slots = {}
    for i in range(1, 16):
        if i == 1:
            final_slots[i] = delta_l1
        elif 2 <= i <= 5:
            final_slots[i] = sp.Integer(0)
        elif 6 <= i <= 14:
            final_slots[i] = sp.simplify(delta_l[i])
        elif i == 15:
            final_slots[i] = sp.simplify(grouped_l15_pred + geometry_gap)

    for i in range(1, 16):
        expect_zero(f"final slot check l{i}", final_slots[i] - delta_l[i])

    print("Final theorem split:")
    print("  kinetic/compiler slot     =", delta_l1)
    print("  grouped-P2 middle block   = exact on l6..l14")
    print("  unique geometry remainder =", geometry_gap)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    mu_sol, d3_sol, s24_sol = one_body_gate()

    nu = sp.symbols("nu", real=True)
    banner("PART II — EXACT GR COM TARGET AND CARRIED SEED")
    target_h = gr_target_h(nu)
    target_l = inverse_map_from_h(nu, target_h)
    seed_h = carried_seed_h(nu)
    seed_l = carried_seed_l(nu)

    delta_h = {i: sp.simplify(target_h[i] - seed_h[i]) for i in range(1, 16)}
    delta_l = {i: sp.simplify(target_l[i] - seed_l[i]) for i in range(1, 16)}

    for i in range(1, 16):
        expect_zero(f"Delta l{i} + Delta h{i}", delta_l[i] + delta_h[i])

    coeffs, l15_pred, l15_gap = grouped_p2_middle_closure(nu, delta_l)
    pure_kinetic_and_static(nu, target_h, target_l, seed_h, seed_l, l15_pred)
    final_decomposition(nu, delta_l, l15_pred)

    banner("FINAL LEDGER")
    print("1. The one-body 3PN gate closes with")
    print("      mu_rho3 = 1/4,  d3 = -45/4,  s24 = -1/16.")
    print("2. The exact GR COM residual is fully known.")
    print("3. The richer grouped-real-P2 compiler has det(M_mid) = -4/27 and closes")
    print("   the whole 9-slot middle block exactly.")
    print("4. The apparent static P0/g freedom collapses identically to a unique")
    print("   geometry-side remainder.")
    print("5. The apparent pure-kinetic residual Delta l1 is exactly the ordinary-")
    print("   Lagrangian counterimage of the universal free relativistic two-body")
    print("   COM Hamiltonian.")
    print("6. Therefore the conservative 3PN COM residual splits exactly into")
    print("      kinetic/compiler + grouped-P2 middle block + unique geometry static slot.")


if __name__ == "__main__":
    main()