
#!/usr/bin/env python3
"""
4pn_local_hamiltonian_to_ordinary_audit.py

Hamiltonian-to-ordinary local 4PN translation audit.

What this script does
---------------------
1. Carries the exact reduced COM ordinary Lagrangian blocks through 3PN.
2. Builds the full 21-slot reduced COM 4PN ordinary-Lagrangian basis:
      - free kinetic : 6 slots,
      - G/r          : 5 slots,
      - G^2/r^2      : 4 slots,
      - G^3/r^3      : 3 slots,
      - G^4/r^4      : 2 slots,
      - G^5/r^5      : 1 slot.
3. Applies the exact quartic perturbative Legendre compiler and extracts the induced
   21-slot COM Hamiltonian map.
4. Verifies that the map is diagonal with Jacobian -I.
5. Imports the fixed-chart *local* 4PN ADM Hamiltonian target and solves the exact
   ordinary-Lagrangian target slot-by-slot.
6. Verifies the strict one-body limit of the ordinary target.
7. Compares the nu-degree structure of:
      - the Hamiltonian target,
      - the translated ordinary target,
      - and the constant-coefficient generic-frame local scaffold.

Main structural result
----------------------
The fixed local Hamiltonian target has interaction-block nu-degree ceilings

    G/r : 4,   G^2/r^2 : 3,   G^3/r^3 : 3,   G^4/r^4 : 2,   G^5/r^5 : 2,

which match the constant-coefficient local generic-frame scaffold.

But after quartic Hamiltonian->ordinary translation, the local ordinary target acquires

    G^2/r^2 : 4,   G^3/r^3 : 4,   G^4/r^4 : 4,

while G^5/r^5 stays at degree 2.

So the clean local 4PN generic-frame lift should be performed in the Hamiltonian
chart first. Translating to the ordinary chart *before* the lift introduces extra
nu^4 tails in the middle/upper local blocks that are not present in the fixed-chart
Hamiltonian target and are not reachable from the present constant-coefficient local
ordinary scaffold.
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


def evenize(expr: sp.Expr, pr: sp.Symbol, pt: sp.Symbol, p2: sp.Symbol, pr2: sp.Symbol) -> sp.Expr:
    """Rewrite an even polynomial in (pr,pt) as a polynomial in (p2,pr2)."""
    expr = sp.expand(expr)

    def repl(node: sp.Basic) -> sp.Basic:
        if isinstance(node, sp.Pow) and node.exp.is_integer and int(node.exp) % 2 == 0:
            n = int(node.exp)
            if node.base == pr:
                return pr2 ** (n // 2)
            if node.base == pt:
                return (p2 - pr2) ** (n // 2)
        return node

    expr = expr.replace(
        lambda z: isinstance(z, sp.Pow) and z.exp.is_integer and int(z.exp) % 2 == 0 and (z.base == pr or z.base == pt),
        repl,
    )
    return sp.expand(expr)


def max_nu_degree(expr: sp.Expr, nu: sp.Symbol) -> int:
    poly = sp.Poly(sp.expand(expr), nu)
    return int(poly.degree())


# ---------------------------------------------------------------------------
# Carried reduced COM lower-order ordinary Lagrangian blocks through 3PN
# ---------------------------------------------------------------------------

def carried_lower_order_blocks() -> tuple[sp.Expr, sp.Expr, sp.Expr, dict[str, sp.Symbol]]:
    nu, u = sp.symbols("nu u", real=True)
    rd, vt = sp.symbols("rd vt", real=True)
    v2 = rd**2 + vt**2
    pi = sp.pi

    l1 = (
        (1 - 3 * nu) / 8 * v2**2
        + u * ((3 + nu) / 2 * v2 + nu / 2 * rd**2)
        - u**2 / 2
    )

    l2 = (
        (1 - 5 * nu + 5 * nu**2) / 16 * v2**3
        + u * (
            (sp.Rational(7, 8) - sp.Rational(7, 4) * nu - nu**2 / 8) * v2**2
            + (nu / 4 - nu**2 / 4) * rd**2 * v2
            + sp.Rational(3, 8) * nu**2 * rd**4
        )
        + u**2 * ((2 - sp.Rational(7, 8) * nu) * v2 + sp.Rational(15, 8) * nu * rd**2)
        + u**3 * (sp.Rational(1, 4) + sp.Rational(3, 4) * nu)
    )

    # Exact carried 3PN COM Hamiltonian target.
    h3 = {
        1: (-5 + 35 * nu - 70 * nu**2 + 35 * nu**3) / 128,
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: (-7 + 42 * nu - 53 * nu**2 - 5 * nu**3) / 16,
        7: (2 - 3 * nu) * nu**2 / 16,
        8: 3 * (1 - nu) * nu**2 / 16,
        9: -5 * nu**3 / 16,
        10: (-27 + 136 * nu + 109 * nu**2) / 16,
        11: (17 + 30 * nu) * nu / 16,
        12: (5 + 43 * nu) * nu / 12,
        13: (-600 + (3 * pi**2 - 1340) * nu - 552 * nu**2) / 192,
        14: -(340 + 3 * pi**2 + 112 * nu) * nu / 64,
        15: (12 + (872 - 63 * pi**2) * nu) / 96,
    }

    # Exact inverse map from the carried 3PN audit.
    l3c = {
        1: sp.simplify(sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3 - h3[1]),
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: sp.simplify(sp.Rational(1, 4) + sp.Rational(7, 8) * nu - sp.Rational(35, 8) * nu**2 - sp.Rational(21, 4) * nu**3 - h3[6]),
        7: sp.simplify(sp.Rational(11, 8) * nu**2 - sp.Rational(9, 2) * nu**3 - h3[7]),
        8: sp.simplify(sp.Rational(3, 4) * nu**2 - sp.Rational(9, 4) * nu**3 - h3[8]),
        9: -h3[9],
        10: sp.simplify(sp.Rational(5, 4) + sp.Rational(15, 8) * nu + sp.Rational(123, 8) * nu**2 + sp.Rational(13, 4) * nu**3 - h3[10]),
        11: sp.simplify(sp.Rational(7, 8) * nu + sp.Rational(41, 8) * nu**2 + sp.Rational(31, 4) * nu**3 - h3[11]),
        12: sp.simplify(sp.Rational(9, 2) * nu**2 + 4 * nu**3 - h3[12]),
        13: sp.simplify(-sp.Rational(3, 2) - sp.Rational(59, 4) * nu - sp.Rational(25, 4) * nu**2 - sp.Rational(1, 2) * nu**3 - h3[13]),
        14: sp.simplify(sp.Rational(7, 4) * nu - sp.Rational(31, 4) * nu**2 - sp.Rational(7, 2) * nu**3 - h3[14]),
        15: -h3[15],
    }

    l3 = (
        l3c[1] * v2**4
        + l3c[2] * v2**3 * rd**2
        + l3c[3] * v2**2 * rd**4
        + l3c[4] * v2 * rd**6
        + l3c[5] * rd**8
        + u * (l3c[6] * v2**3 + l3c[7] * v2**2 * rd**2 + l3c[8] * v2 * rd**4 + l3c[9] * rd**6)
        + u**2 * (l3c[10] * v2**2 + l3c[11] * v2 * rd**2 + l3c[12] * rd**4)
        + u**3 * (l3c[13] * v2 + l3c[14] * rd**2)
        + u**4 * l3c[15]
    )

    return l1, l2, l3, {"nu": nu, "u": u, "rd": rd, "vt": vt, "v2": v2}


# ---------------------------------------------------------------------------
# Exact quartic COM map
# ---------------------------------------------------------------------------

def quartic_com_map() -> tuple[dict[int, sp.Expr], tuple[sp.Symbol, ...], dict[str, sp.Symbol]]:
    banner("PART I — EXACT 21-SLOT COM QUARTIC MAP")

    l1, l2, l3, syms = carried_lower_order_blocks()
    nu, u, rd, vt, v2 = syms["nu"], syms["u"], syms["rd"], syms["vt"], syms["v2"]
    pr, pt = sp.symbols("pr pt", real=True)
    p2, pr2 = sp.symbols("p2 pr2", real=True)

    coeffs = sp.symbols("L1:22")
    # 6 free + 5 + 4 + 3 + 2 + 1 = 21 slots.
    l4 = (
        coeffs[0] * v2**5
        + coeffs[1] * v2**4 * rd**2
        + coeffs[2] * v2**3 * rd**4
        + coeffs[3] * v2**2 * rd**6
        + coeffs[4] * v2 * rd**8
        + coeffs[5] * rd**10
        + u * (
            coeffs[6] * v2**4
            + coeffs[7] * v2**3 * rd**2
            + coeffs[8] * v2**2 * rd**4
            + coeffs[9] * v2 * rd**6
            + coeffs[10] * rd**8
        )
        + u**2 * (
            coeffs[11] * v2**3
            + coeffs[12] * v2**2 * rd**2
            + coeffs[13] * v2 * rd**4
            + coeffs[14] * rd**6
        )
        + u**3 * (
            coeffs[15] * v2**2
            + coeffs[16] * v2 * rd**2
            + coeffs[17] * rd**4
        )
        + u**4 * (
            coeffs[18] * v2
            + coeffs[19] * rd**2
        )
        + u**5 * coeffs[20]
    )

    A0 = sp.Matrix([sp.diff(l1, rd), sp.diff(l1, vt)]).subs({rd: pr, vt: pt})
    B0 = sp.Matrix([sp.diff(l2, rd), sp.diff(l2, vt)]).subs({rd: pr, vt: pt})
    D0 = sp.Matrix([sp.diff(l3, rd), sp.diff(l3, vt)]).subs({rd: pr, vt: pt})
    C0 = sp.hessian(l1, (rd, vt)).subs({rd: pr, vt: pt})
    E0 = sp.hessian(l2, (rd, vt)).subs({rd: pr, vt: pt})

    # T0[A0,A0,A0]
    vars2 = (rd, vt)
    Acomp = [A0[0], A0[1]]
    Tcontr = 0
    for i, xi in enumerate(vars2):
        for j, xj in enumerate(vars2):
            for k, xk in enumerate(vars2):
                Tcontr += sp.diff(l1, xi, xj, xk).subs({rd: pr, vt: pt}) * Acomp[i] * Acomp[j] * Acomp[k]

    h4_expr = sp.expand(
        -l4.subs({rd: pr, vt: pt})
        + A0.dot(D0)
        + sp.Rational(1, 2) * B0.dot(B0)
        - (B0.T * C0 * A0)[0]
        - sp.Rational(1, 2) * (A0.T * E0 * A0)[0]
        + sp.Rational(1, 2) * (A0.T * C0 * C0 * A0)[0]
        + sp.Rational(1, 6) * Tcontr
    )

    poly = sp.Poly(evenize(h4_expr, pr, pt, p2, pr2), p2, pr2, u)
    terms = dict(poly.terms())

    index_to_monom: dict[int, tuple[int, int, int]] = {}
    idx = 1
    for a_pow, b_pow in [(5, 0), (4, 1), (3, 2), (2, 3), (1, 4), (0, 5)]:
        index_to_monom[idx] = (a_pow, b_pow, 0)
        idx += 1
    for a_pow, b_pow in [(4, 0), (3, 1), (2, 2), (1, 3), (0, 4)]:
        index_to_monom[idx] = (a_pow, b_pow, 1)
        idx += 1
    for a_pow, b_pow in [(3, 0), (2, 1), (1, 2), (0, 3)]:
        index_to_monom[idx] = (a_pow, b_pow, 2)
        idx += 1
    for a_pow, b_pow in [(2, 0), (1, 1), (0, 2)]:
        index_to_monom[idx] = (a_pow, b_pow, 3)
        idx += 1
    for a_pow, b_pow in [(1, 0), (0, 1)]:
        index_to_monom[idx] = (a_pow, b_pow, 4)
        idx += 1
    index_to_monom[idx] = (0, 0, 5)

    hmap = {i: sp.simplify(terms.get(index_to_monom[i], 0)) for i in range(1, 22)}

    subbanner("I.1 — Extracted slot map h_i(L_j)")
    for i in range(1, 22):
        print(f"h{i} = {sp.factor(hmap[i])}")

    subbanner("I.2 — Exact diagonal Jacobian")
    J = sp.Matrix([[sp.diff(hmap[i], coeffs[j]) for j in range(21)] for i in range(1, 22)])
    expect_zero("Jacobian + I", J + sp.eye(21))

    return hmap, coeffs, {**syms, "pr": pr, "pt": pt, "p2": p2, "pr2": pr2}


# ---------------------------------------------------------------------------
# Fixed-chart local 4PN Hamiltonian target (full 21-slot COM basis)
# ---------------------------------------------------------------------------

def local_hamiltonian_target() -> dict[int, sp.Expr]:
    banner("PART II — FIXED-CHART LOCAL 4PN HAMILTONIAN TARGET")

    nu = sp.symbols("nu", real=True)
    pi = sp.pi

    h = {
        # free block
        1: sp.Rational(7, 256) - sp.Rational(63, 256) * nu + sp.Rational(189, 256) * nu**2 - sp.Rational(105, 128) * nu**3 + sp.Rational(63, 256) * nu**4,
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: 0,

        # G/r
        7: sp.Rational(45, 128) - sp.Rational(45, 16) * nu + sp.Rational(423, 64) * nu**2 - sp.Rational(1013, 256) * nu**3 - sp.Rational(35, 128) * nu**4,
        8: -sp.Rational(3, 32) * nu**2 + sp.Rational(23, 64) * nu**3 - sp.Rational(5, 32) * nu**4,
        9: -sp.Rational(9, 64) * nu**2 + sp.Rational(69, 128) * nu**3 - sp.Rational(9, 64) * nu**4,
        10: -sp.Rational(5, 64) * nu**3 - sp.Rational(5, 32) * nu**4,
        11: sp.Rational(35, 256) * nu**3 - sp.Rational(35, 128) * nu**4,

        # G^2/r^2
        12: sp.Rational(13, 8) - sp.Rational(791, 64) * nu + sp.Rational(4857, 256) * nu**2 + sp.Rational(2335, 256) * nu**3,
        13: sp.Rational(49, 16) * nu - sp.Rational(545, 64) * nu**2 + sp.Rational(1135, 256) * nu**3,
        14: -sp.Rational(889, 192) * nu + sp.Rational(9475, 768) * nu**2 - sp.Rational(1649, 768) * nu**3,
        15: sp.Rational(369, 160) * nu - sp.Rational(1151, 128) * nu**2 + sp.Rational(10353, 1280) * nu**3,

        # G^3/r^3
        16: sp.Rational(105, 32) + (sp.Rational(2749, 8192) * pi**2 - sp.Rational(589189, 19200)) * nu + (sp.Rational(18491, 16384) * pi**2 - sp.Rational(1189789, 28800)) * nu**2 - sp.Rational(553, 128) * nu**3,
        17: (sp.Rational(63347, 1600) - sp.Rational(1059, 1024) * pi**2) * nu + (-sp.Rational(127, 3) - sp.Rational(4035, 2048) * pi**2) * nu**2 - sp.Rational(225, 64) * nu**3,
        18: (sp.Rational(375, 8192) * pi**2 - sp.Rational(23533, 1280)) * nu + (sp.Rational(57563, 1920) - sp.Rational(38655, 16384) * pi**2) * nu**2 - sp.Rational(381, 128) * nu**3,

        # G^4/r^4
        19: sp.Rational(105, 32) + (sp.Rational(185761, 19200) - sp.Rational(21837, 8192) * pi**2) * nu + (sp.Rational(672811, 19200) - sp.Rational(158177, 49152) * pi**2) * nu**2,
        20: (sp.Rational(3401779, 57600) - sp.Rational(28691, 24576) * pi**2) * nu + (sp.Rational(110099, 49152) * pi**2 - sp.Rational(21827, 3840)) * nu**2,

        # G^5/r^5
        21: -sp.Rational(1, 16) + (sp.Rational(6237, 1024) * pi**2 - sp.Rational(169199, 2400)) * nu + (sp.Rational(7403, 3072) * pi**2 - sp.Rational(1256, 45)) * nu**2,
    }

    for i in range(1, 22):
        print(f"h{i}^(loc4PN) = {sp.factor(h[i])}")

    return h


# ---------------------------------------------------------------------------
# Solve exact ordinary local target
# ---------------------------------------------------------------------------

def ordinary_target_from_hamiltonian(hmap: dict[int, sp.Expr], coeffs: tuple[sp.Symbol, ...], target_h: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    banner("PART III — EXACT ORDINARY LOCAL 4PN TARGET")

    feedback = {i: sp.simplify(hmap[i].subs({coeffs[j]: 0 for j in range(21)})) for i in range(1, 22)}
    target_l = {i: sp.simplify(feedback[i] - target_h[i]) for i in range(1, 22)}

    subbanner("III.1 — Ordinary slot coefficients L_i^(target)")
    for i in range(1, 22):
        print(f"L{i}^(target) = {sp.factor(target_l[i])}")

    subbanner("III.2 — Verify h_i = feedback_i - L_i")
    for i in range(1, 22):
        expect_zero(f"map check slot {i}", hmap[i].subs({coeffs[j]: target_l[j + 1] for j in range(21)}) - target_h[i])

    return target_l


# ---------------------------------------------------------------------------
# One-body checks and natural seed
# ---------------------------------------------------------------------------

def even_nu(expr: sp.Expr, nu: sp.Symbol, Delta: sp.Symbol) -> sp.Expr:
    expr = sp.expand((expr + expr.subs(Delta, -Delta)) / 2)
    for n in range(40, 1, -2):
        expr = sp.expand(expr.subs(Delta**n, (1 - 4 * nu) ** (n // 2)))
    while expr.has(Delta**2):
        expr = sp.expand(expr.subs(Delta**2, 1 - 4 * nu))
    if expr.has(Delta):
        expr = sp.expand(expr.subs(Delta, 0))
    return sp.expand(expr)


def natural_seed_local_ordinary() -> dict[int, sp.Expr]:
    nu, Delta = sp.symbols("nu Delta", real=True)
    Xa = (1 + Delta) / 2
    Xb = (1 - Delta) / 2

    seed = {i: sp.Integer(0) for i in range(1, 22)}
    seed[1] = sp.simplify(sp.Rational(7, 256) * even_nu(Xa**9 + Xb**9, nu, Delta))
    seed[7] = sp.simplify(sp.Rational(75, 128) * even_nu(Xa**8 + Xb**8, nu, Delta))
    seed[12] = sp.simplify(sp.Rational(59, 16) * even_nu(Xa**7 + Xb**7, nu, Delta))
    seed[16] = sp.simplify(sp.Rational(203, 32) * even_nu(Xa**6 + Xb**6, nu, Delta))
    seed[19] = sp.simplify(sp.Rational(31, 32) * even_nu(Xa**5 + Xb**5, nu, Delta))
    seed[21] = sp.simplify(sp.Rational(1, 16) * even_nu(Xa**4 + Xb**4, nu, Delta))
    return seed


def onebody_and_seed_checks(target_l: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    banner("PART IV — ONE-BODY GATE AND NATURAL LOCAL SEED")

    nu = sp.symbols("nu", real=True)

    # strict one-body gate
    expect_zero("L1 one-body gate", target_l[1].subs(nu, 0) - sp.Rational(7, 256))
    expect_zero("L7 one-body gate", target_l[7].subs(nu, 0) - sp.Rational(75, 128))
    expect_zero("L12 one-body gate", target_l[12].subs(nu, 0) - sp.Rational(59, 16))
    expect_zero("L16 one-body gate", target_l[16].subs(nu, 0) - sp.Rational(203, 32))
    expect_zero("L19 one-body gate", target_l[19].subs(nu, 0) - sp.Rational(31, 32))
    expect_zero("L21 one-body gate", target_l[21].subs(nu, 0) - sp.Rational(1, 16))

    for i in [2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 20]:
        expect_zero(f"subleading one-body slot {i}", sp.simplify(target_l[i].subs(nu, 0)))

    seed = natural_seed_local_ordinary()
    subbanner("IV.1 — Natural local self/static seed coefficients")
    for i in [1, 7, 12, 16, 19, 21]:
        print(f"L{i}^(seed) = {sp.factor(seed[i])}")

    subbanner("IV.2 — Residual beyond the natural local seed")
    delta = {i: sp.simplify(target_l[i] - seed[i]) for i in range(1, 22)}
    for i in range(1, 22):
        if delta[i] != 0:
            print(f"Delta L{i} = {sp.factor(delta[i])}")

    return seed


# ---------------------------------------------------------------------------
# Constant-coefficient local scaffold degree ceilings
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
                            tm = sp.expand(sym.subs({b: 0, c: 0, e: 0, p: 0, q: 1}))
                            if tm != 0:
                                continue
                            basis.add(sym)

    return sorted(basis, key=str)


def scaffold_degree_ceilings() -> dict[str, int]:
    nu, Delta = sp.symbols("nu Delta", real=True)
    a, b, c, d, e, p, q = sp.symbols("a b c d e p q")
    V2, rd = sp.symbols("V2 rd", real=True)
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
        for n in range(40, 1, -2):
            expr = sp.expand(expr.subs(Delta**n, (1 - 4 * nu) ** (n // 2)))
        while expr.has(Delta**2):
            expr = sp.expand(expr.subs(Delta**2, 1 - 4 * nu))
        if expr.has(Delta):
            expr = sp.expand(expr.subs(Delta, 0))
        return sp.expand(expr)

    banner("PART V — CONSTANT-COEFFICIENT SCAFFOLD nu-DEGREE CEILINGS")

    ceilings = {}
    for tag, mass_deg, vel_deg in [("Q", 0, 8), ("T", 1, 6), ("S", 2, 4), ("U", 3, 2), ("W", 4, 0)]:
        basis = generate_basis(mass_deg, vel_deg)
        degs = [max_nu_degree(to_nu(expr), nu) for expr in basis]
        ceilings[tag] = max(degs)
        print(f"{tag} block: basis count = {len(basis)}, max nu-degree = {ceilings[tag]}")

    return ceilings


# ---------------------------------------------------------------------------
# Degree comparison theorem
# ---------------------------------------------------------------------------

def degree_comparison(target_h: dict[int, sp.Expr], target_l: dict[int, sp.Expr], ceilings: dict[str, int]) -> None:
    banner("PART VI — DEGREE COMPARISON THEOREM")

    nu = sp.symbols("nu", real=True)

    h_blocks = {
        "Q": [7, 8, 9, 10, 11],
        "T": [12, 13, 14, 15],
        "S": [16, 17, 18],
        "U": [19, 20],
        "W": [21],
    }
    l_blocks = h_blocks

    print("Hamiltonian target degrees by block:")
    h_deg = {}
    for tag, inds in h_blocks.items():
        h_deg[tag] = max(max_nu_degree(target_h[i], nu) for i in inds)
        print(f"  {tag}: {h_deg[tag]}")

    print("\nOrdinary target degrees by block:")
    l_deg = {}
    for tag, inds in l_blocks.items():
        l_deg[tag] = max(max_nu_degree(target_l[i], nu) for i in inds)
        print(f"  {tag}: {l_deg[tag]}")

    print("\nConstant-coefficient scaffold ceilings by block:")
    for tag in ["Q", "T", "S", "U", "W"]:
        print(f"  {tag}: {ceilings[tag]}")

    # Exact match in the Hamiltonian chart.
    for tag in ["Q", "T", "S", "U", "W"]:
        if h_deg[tag] != ceilings[tag]:
            raise AssertionError(f"Hamiltonian block {tag} does not match the scaffold ceiling.")

    # Ordinary-chart obstructions.
    if not (l_deg["T"] > ceilings["T"] and l_deg["S"] > ceilings["S"] and l_deg["U"] > ceilings["U"]):
        raise AssertionError("Expected ordinary-chart nu^4 obstructions were not detected.")
    if l_deg["W"] != ceilings["W"]:
        raise AssertionError("G^5/r^5 block unexpectedly changed degree after translation.")

    print("\nTheorem:")
    print("  The fixed local Hamiltonian target matches the constant-coefficient generic-frame")
    print("  scaffold ceilings blockwise.")
    print("  After quartic Hamiltonian->ordinary translation, the ordinary target acquires")
    print("  extra nu^4 tails in the G^2/r^2, G^3/r^3, and G^4/r^4 blocks, while G^5/r^5")
    print("  remains at the same degree ceiling.")
    print("  Therefore the clean local 4PN generic-frame lift should be done Hamiltonian-first.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    hmap, coeffs, syms = quartic_com_map()
    target_h = local_hamiltonian_target()
    target_l = ordinary_target_from_hamiltonian(hmap, coeffs, target_h)
    onebody_and_seed_checks(target_l)
    ceilings = scaffold_degree_ceilings()
    degree_comparison(target_h, target_l, ceilings)

    banner("FINAL LEDGER")
    print("1. The full 21-slot COM quartic map is diagonal with Jacobian -I.")
    print("2. The fixed local 4PN Hamiltonian target translates uniquely into an exact ordinary target.")
    print("3. The ordinary target preserves the strict one-body 4PN gate exactly.")
    print("4. The Hamiltonian target matches the constant-coefficient generic-frame local scaffold")
    print("   blockwise in nu-degree ceilings.")
    print("5. The ordinary target introduces extra nu^4 tails in the G^2/r^2, G^3/r^3, and G^4/r^4")
    print("   blocks that are not present in the fixed Hamiltonian target and are not reachable from")
    print("   the present constant-coefficient local ordinary scaffold.")
    print("6. So the next clean local 4PN step is the Hamiltonian-chart generic-frame lift, not an")
    print("   ordinary-chart lift done before the fixed-chart conditions are imposed.")
