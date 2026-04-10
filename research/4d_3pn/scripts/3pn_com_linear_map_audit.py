#!/usr/bin/env python3
"""
3pn_com_linear_map_audit.py

Center-of-mass 3PN ordinary-Lagrangian -> Hamiltonian compiler.

What this script does
---------------------
1. Carries forward the reduced COM Newtonian + 1PN + 2PN Lagrangian blocks
   implied by the existing 1PN/2PN toy-model derivations.
2. Introduces the most general isotropic, time-translation invariant 3PN COM
   Lagrangian basis with 15 coefficients.
3. Uses the exact cubic Legendre-transform identity already derived for the 3PN
   scaffold,

       H3 = -L3(v0) + A0·B0 - 1/2 A0^T C0 A0,

   with v0 = p, A0 = ∂L1|_{v0}, B0 = ∂L2|_{v0}, C0 = ∂²L1|_{v0},
   to compute the induced 3PN COM Hamiltonian.
4. Extracts the 15 Hamiltonian coefficients h_i in the standard COM basis and
   proves that the map is diagonal in the chosen ordinary-Lagrangian basis:

       h_i = F_i(nu) - l_i   (up to the obvious zero-feedback slots).

5. Evaluates the carried self/static 3PN seed on this map and prints its image
   in the COM Hamiltonian basis.

Interpretation
--------------
This does NOT yet import the final GR 3PN target coefficients. What it gives is
exactly the piece that was still missing between the scaffold and a real solve:

  * once the target COM Hamiltonian coefficients h_i are supplied,
    the ordinary 3PN COM Lagrangian coefficients l_i follow immediately;
  * the remaining work is then the lift from COM back to a clean generic-frame
    comparable-mass derivation.
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


def evenize(expr: sp.Expr, pr: sp.Symbol, pt: sp.Symbol, p2: sp.Symbol, pr2: sp.Symbol) -> sp.Expr:
    """Rewrite an even polynomial in (pr, pt) as a polynomial in (p2, pr2),
    where p2 = pr^2 + pt^2 and pr2 = pr^2.
    """
    expr = sp.expand(expr)

    def repl_pow(node: sp.Basic) -> sp.Basic:
        if isinstance(node, sp.Pow) and node.exp.is_integer and int(node.exp) % 2 == 0:
            n = int(node.exp)
            if node.base == pr:
                return pr2 ** (n // 2)
            if node.base == pt:
                return (p2 - pr2) ** (n // 2)
        return node

    expr = expr.replace(
        lambda z: isinstance(z, sp.Pow) and z.exp.is_integer and int(z.exp) % 2 == 0 and (z.base == pr or z.base == pt),
        repl_pow,
    )
    expr = sp.expand(expr)
    return sp.simplify(expr)


# ---------------------------------------------------------------------------
# Carried COM lower-order Lagrangian blocks
# ---------------------------------------------------------------------------

def carried_lower_order_blocks() -> tuple[sp.Expr, sp.Expr, sp.Expr, dict[str, sp.Symbol]]:
    nu, u = sp.symbols("nu u", real=True)
    rd, vt = sp.symbols("rd vt", real=True)
    v2 = rd**2 + vt**2

    # Reduced COM Lagrangian blocks L/μ with u = GM/r in dimensionless units.
    l0 = sp.Rational(1, 2) * v2 + u

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

    return l0, l1, l2, {"nu": nu, "u": u, "rd": rd, "vt": vt, "v2": v2}


# ---------------------------------------------------------------------------
# Generic 3PN COM Lagrangian basis and exact H3 map
# ---------------------------------------------------------------------------

def h3_linear_map() -> tuple[dict[int, sp.Expr], tuple[sp.Symbol, ...], dict[str, sp.Symbol]]:
    banner("PART I — EXACT COM 3PN LINEAR MAP")

    l0, l1, l2, syms = carried_lower_order_blocks()
    nu, u, rd, vt, v2 = syms["nu"], syms["u"], syms["rd"], syms["vt"], syms["v2"]
    pr, pt = sp.symbols("pr pt", real=True)
    p2, pr2 = sp.symbols("p2 pr2", real=True)

    coeffs = sp.symbols("l1:16")
    l1c, l2c, l3c, l4c, l5c, l6c, l7c, l8c, l9c, l10c, l11c, l12c, l13c, l14c, l15c = coeffs

    l3 = (
        l1c * v2**4
        + l2c * v2**3 * rd**2
        + l3c * v2**2 * rd**4
        + l4c * v2 * rd**6
        + l5c * rd**8
        + u * (l6c * v2**3 + l7c * v2**2 * rd**2 + l8c * v2 * rd**4 + l9c * rd**6)
        + u**2 * (l10c * v2**2 + l11c * v2 * rd**2 + l12c * rd**4)
        + u**3 * (l13c * v2 + l14c * rd**2)
        + u**4 * l15c
    )

    # Exact cubic Legendre-transform identity at 3PN for unit mass matrix.
    A0 = sp.Matrix([sp.diff(l1, rd), sp.diff(l1, vt)]).subs({rd: pr, vt: pt})
    B0 = sp.Matrix([sp.diff(l2, rd), sp.diff(l2, vt)]).subs({rd: pr, vt: pt})
    C0 = sp.hessian(l1, (rd, vt)).subs({rd: pr, vt: pt})
    h3_expr = sp.expand(-l3.subs({rd: pr, vt: pt}) + A0.dot(B0) - sp.Rational(1, 2) * (A0.T * C0 * A0)[0])

    h3_poly = sp.Poly(evenize(h3_expr, pr, pt, p2, pr2), p2, pr2, u)
    terms = dict(h3_poly.terms())

    # Standard COM Hamiltonian basis:
    # p2^4, p2^3 pr2, p2^2 pr2^2, p2 pr2^3, pr2^4,
    # u p2^3, u p2^2 pr2, u p2 pr2^2, u pr2^3,
    # u^2 p2^2, u^2 p2 pr2, u^2 pr2^2,
    # u^3 p2, u^3 pr2,
    # u^4.
    index_to_monom = {
        1: (4, 0, 0),
        2: (3, 1, 0),
        3: (2, 2, 0),
        4: (1, 3, 0),
        5: (0, 4, 0),
        6: (3, 0, 1),
        7: (2, 1, 1),
        8: (1, 2, 1),
        9: (0, 3, 1),
        10: (2, 0, 2),
        11: (1, 1, 2),
        12: (0, 2, 2),
        13: (1, 0, 3),
        14: (0, 1, 3),
        15: (0, 0, 4),
    }
    h_map = {i: sp.simplify(terms[index_to_monom[i]]) for i in range(1, 16)}

    subbanner("I.1 — Extracted Hamiltonian coefficients h_i")
    for i in range(1, 16):
        print(f"h{i} = {h_map[i]}")

    subbanner("I.2 — Inverse map l_i(h_j)")
    # The chosen Lagrangian basis diagonalizes the map.
    hsyms = {i: sp.Symbol(f"h{i}") for i in range(1, 16)}
    inverse = {
        1: sp.simplify(sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3 - hsyms[1]),
        2: sp.simplify(-hsyms[2]),
        3: sp.simplify(-hsyms[3]),
        4: sp.simplify(-hsyms[4]),
        5: sp.simplify(-hsyms[5]),
        6: sp.simplify(sp.Rational(1, 4) + sp.Rational(7, 8) * nu - sp.Rational(35, 8) * nu**2 - sp.Rational(21, 4) * nu**3 - hsyms[6]),
        7: sp.simplify(sp.Rational(11, 8) * nu**2 - sp.Rational(9, 2) * nu**3 - hsyms[7]),
        8: sp.simplify(sp.Rational(3, 4) * nu**2 - sp.Rational(9, 4) * nu**3 - hsyms[8]),
        9: sp.simplify(-hsyms[9]),
        10: sp.simplify(sp.Rational(5, 4) + sp.Rational(15, 8) * nu + sp.Rational(123, 8) * nu**2 + sp.Rational(13, 4) * nu**3 - hsyms[10]),
        11: sp.simplify(sp.Rational(7, 8) * nu + sp.Rational(41, 8) * nu**2 + sp.Rational(31, 4) * nu**3 - hsyms[11]),
        12: sp.simplify(sp.Rational(9, 2) * nu**2 + 4 * nu**3 - hsyms[12]),
        13: sp.simplify(-sp.Rational(3, 2) - sp.Rational(59, 4) * nu - sp.Rational(25, 4) * nu**2 - sp.Rational(1, 2) * nu**3 - hsyms[13]),
        14: sp.simplify(sp.Rational(7, 4) * nu - sp.Rational(31, 4) * nu**2 - sp.Rational(7, 2) * nu**3 - hsyms[14]),
        15: sp.simplify(-hsyms[15]),
    }
    for i in range(1, 16):
        print(f"l{i} = {inverse[i]}")

    # Direct verification of the inverse map.
    inverse_subs = {coeffs[i - 1]: inverse[i] for i in range(1, 16)}
    for i in range(1, 16):
        expect_zero(f"inverse check h{i}", h_map[i].subs(inverse_subs) - hsyms[i])

    return h_map, coeffs, {**syms, "pr": pr, "pt": pt, "p2": p2, "pr2": pr2}


# ---------------------------------------------------------------------------
# Carried self/static seed image in the COM Hamiltonian basis
# ---------------------------------------------------------------------------

def self_static_seed_image(h_map: dict[int, sp.Expr], coeffs: tuple[sp.Symbol, ...], syms: dict[str, sp.Symbol]) -> None:
    banner("PART II — COM IMAGE OF THE CARRIED 3PN SELF/STATIC SEED")

    nu, u, rd, vt, v2 = syms["nu"], syms["u"], syms["rd"], syms["vt"], syms["v2"]

    seed = {
        coeffs[0]: sp.Rational(5, 128) - sp.Rational(35, 128) * nu + sp.Rational(35, 64) * nu**2 - sp.Rational(35, 128) * nu**3,
        coeffs[1]: 0,
        coeffs[2]: 0,
        coeffs[3]: 0,
        coeffs[4]: 0,
        coeffs[5]: sp.Rational(11, 16) - sp.Rational(33, 8) * nu + sp.Rational(99, 16) * nu**2 - sp.Rational(11, 8) * nu**3,
        coeffs[6]: 0,
        coeffs[7]: 0,
        coeffs[8]: 0,
        coeffs[9]: sp.Rational(47, 16) - sp.Rational(235, 16) * nu + sp.Rational(235, 16) * nu**2,
        coeffs[10]: 0,
        coeffs[11]: 0,
        coeffs[12]: sp.Rational(13, 8) - sp.Rational(13, 2) * nu + sp.Rational(13, 4) * nu**2,
        coeffs[13]: 0,
        coeffs[14]: -sp.Rational(1, 8) + sp.Rational(3, 8) * nu,
    }

    subbanner("II.1 — Seed coefficients l_i")
    for i in range(1, 16):
        print(f"l{i}^(seed) = {sp.simplify(seed[coeffs[i-1]])}")

    subbanner("II.2 — Hamiltonian image h_i^(seed)")
    h_seed = {i: sp.simplify(h_map[i].subs(seed)) for i in range(1, 16)}
    for i in range(1, 16):
        print(f"h{i}^(seed) = {h_seed[i]}")

    print("\nInterpretation:")
    print("  - The carried self/static seed populates only the v^8, u v^6, u^2 v^4, u^3 v^2, u^4")
    print("    ordinary-Lagrangian slots.")
    print("  - After the exact cubic Legendre map, lower-order feedback fills several additional")
    print("    Hamiltonian coefficients.")
    print("  - Therefore the genuine 3PN comparable-mass problem is the residual between the")
    print("    eventual target h_i and these seed images h_i^(seed).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    h_map, coeffs, syms = h3_linear_map()
    self_static_seed_image(h_map, coeffs, syms)
