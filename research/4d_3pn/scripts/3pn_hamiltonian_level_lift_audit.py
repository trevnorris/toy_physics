#!/usr/bin/env python3
"""
3pn_hamiltonian_level_lift_audit.py

Carry the imported generic-frame ordinary 3PN target through the exact generic-frame
cubic Legendre transform *before* reducing to the center-of-mass frame.

Main point
----------
The previous 3PN step showed that naive COM substitution into the generic-frame
ordinary Lagrangian does not reproduce the solved COM ordinary target.  This
script checks the real next move: transform first, reduce second.

Result:
-------
Using the carried generic-frame 1PN/2PN blocks together with the imported ordinary
3PN target, the exact Hamiltonian-level COM reduction reproduces the standard GR
15-slot 3PN COM Hamiltonian target *exactly*.

So the earlier ordinary mismatch is a reduction-ordering artifact, not a failure of
the generic-frame target itself.
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


# ---------------------------------------------------------------------------
# Ordinary generic-frame 1PN / 2PN / 3PN inputs
# ---------------------------------------------------------------------------

def ordinary_blocks() -> tuple[sp.Expr, sp.Expr, sp.Expr, dict[str, sp.Symbol]]:
    G, r = sp.symbols("G r", positive=True, real=True)
    mA, mB = sp.symbols("mA mB", positive=True, real=True)
    a, b, c, d, e = sp.symbols("a b c d e", real=True)

    L1 = (
        sp.Rational(1, 8) * (mA * a**2 + mB * b**2)
        + G * mA * mB / r * (sp.Rational(3, 2) * (a + b) - sp.Rational(7, 2) * c - sp.Rational(1, 2) * d * e)
        - G**2 * mA * mB * (mA + mB) / (2 * r**2)
    )

    L2 = (
        sp.Rational(1, 16) * (mA * a**3 + mB * b**3)
        + G * mA * mB / r
        * (
            sp.Rational(7, 8) * (a**2 + b**2)
            - sp.Rational(7, 4) * c * (a + b)
            - sp.Rational(1, 4) * d * e * (a + b)
            + sp.Rational(11, 8) * a * b
            + sp.Rational(1, 4) * c**2
            - sp.Rational(5, 8) * (a * e**2 + b * d**2)
            + sp.Rational(3, 2) * c * d * e
            + sp.Rational(3, 8) * d**2 * e**2
        )
        + G**2 * mA * mB / r**2
        * (
            (2 * mB + sp.Rational(11, 8) * mA) * a
            + (2 * mA + sp.Rational(11, 8) * mB) * b
            - sp.Rational(15, 4) * (mA + mB) * c
            + sp.Rational(15, 8) * (mA * d**2 + mB * e**2)
        )
        + G**3 * mA * mB * (mA**2 + 5 * mA * mB + mB**2) / (4 * r**3)
    )

    lam = sp.Rational(-1987, 3080)
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
        sp.Rational(5, 128) * mA * a**4
        + G * mA * mB / r * Q_disp
        + G**2 * mA**2 * mB / r**2 * T_disp
        + G**3 * mA**2 * mB**2 / r**3 * S22_disp
        + G**3 * mA**3 * mB / r**3 * S31_disp
        + G**4 * mA**4 * mB / r**4 * U41_disp
        + G**4 * mA**3 * mB**2 / r**4 * U32_disp
    )
    swap = {a: b, b: a, c: c, d: e, e: d, mA: mB, mB: mA}
    L3 = sp.expand(displayed + displayed.xreplace(swap))

    return L1, L2, L3, {"G": G, "r": r, "mA": mA, "mB": mB, "a": a, "b": b, "c": c, "d": d, "e": e}


# ---------------------------------------------------------------------------
# COM Hamiltonian lift using invariant directional calculus
# ---------------------------------------------------------------------------

def com_hamiltonian_from_generic_lift() -> tuple[sp.Expr, dict[int, sp.Expr]]:
    banner("PART I — HAMILTONIAN-LEVEL LIFT BEFORE COM REDUCTION")

    L1, L2, L3, syms = ordinary_blocks()
    G, r, mA, mB, a, b, c, d, e = (syms[k] for k in ["G", "r", "mA", "mB", "a", "b", "c", "d", "e"])

    Delta, nu, Mtot = sp.symbols("Delta nu Mtot", real=True)
    P2, pr, u = sp.symbols("P2 pr u", real=True)

    Xa = (1 + Delta) / 2
    Xb = (1 - Delta) / 2
    mu = Xa * Xb * Mtot

    # COM values of the ordinary invariants evaluated on v0 = p/m.
    com_subs = {
        mA: Xa * Mtot,
        mB: Xb * Mtot,
        a: Xb**2 * P2,
        b: Xa**2 * P2,
        c: -Xa * Xb * P2,
        d: Xb * pr,
        e: -Xa * pr,
    }

    # Vector algebra in the 2-component representation v = alpha * P + beta * n,
    # with P·P = P2, P·n = pr, n·n = 1.
    def vec_dot(v: tuple[sp.Expr, sp.Expr], w: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
        ap, an = v
        bp, bn = w
        return sp.expand(ap * bp * P2 + (ap * bn + an * bp) * pr + an * bn)

    vA = (Xb, sp.Integer(0))
    vB = (-Xa, sp.Integer(0))

    def grad_com(F: sp.Expr, body: str) -> tuple[sp.Expr, sp.Expr]:
        if body == "A":
            Fa = sp.diff(F, a).subs(com_subs)
            Fc = sp.diff(F, c).subs(com_subs)
            Fd = sp.diff(F, d).subs(com_subs)
            return sp.simplify(2 * Fa * Xb - Fc * Xa), sp.simplify(Fd)
        Fb = sp.diff(F, b).subs(com_subs)
        Fc = sp.diff(F, c).subs(com_subs)
        Fe = sp.diff(F, e).subs(com_subs)
        return sp.simplify(Fc * Xb - 2 * Fb * Xa), sp.simplify(Fe)

    def w1_body(F: sp.Expr, body: str) -> tuple[sp.Expr, sp.Expr]:
        gp, gn = grad_com(F, body)
        mass = Xa * Mtot if body == "A" else Xb * Mtot
        return sp.simplify(gp / mass), sp.simplify(gn / mass)

    w1A = w1_body(L1, "A")
    w1B = w1_body(L1, "B")
    g2A = grad_com(L2, "A")
    g2B = grad_com(L2, "B")

    subbanner("I.1 — Exact COM first correction v1 = M^{-1} A0")
    print("w1_A =", w1A)
    print("w1_B =", w1B)

    term2 = sp.simplify(vec_dot(w1A, g2A) + vec_dot(w1B, g2B))

    def directional_second(F: sp.Expr) -> sp.Expr:
        da = 2 * vec_dot(vA, w1A)
        db = 2 * vec_dot(vB, w1B)
        dc = vec_dot(vB, w1A) + vec_dot(vA, w1B)
        dd = w1A[0] * pr + w1A[1]
        de = w1B[0] * pr + w1B[1]
        d2a = 2 * vec_dot(w1A, w1A)
        d2b = 2 * vec_dot(w1B, w1B)
        d2c = 2 * vec_dot(w1A, w1B)

        vars_inv = [a, b, c, d, e]
        deltas = [da, db, dc, dd, de]
        total = sp.Integer(0)
        for i, xi in enumerate(vars_inv):
            Fi = sp.diff(F, xi).subs(com_subs)
            if i == 0:
                total += Fi * d2a
            elif i == 1:
                total += Fi * d2b
            elif i == 2:
                total += Fi * d2c
        for i, xi in enumerate(vars_inv):
            for j, xj in enumerate(vars_inv):
                Fij = sp.diff(F, xi, xj).subs(com_subs)
                total += Fij * deltas[i] * deltas[j]
        return sp.simplify(sp.expand(total))

    term3 = sp.simplify(-sp.Rational(1, 2) * directional_second(L1))
    L3_com = sp.expand(L3.subs(com_subs))

    H3 = sp.simplify(sp.expand((-L3_com + term2 + term3) / mu))
    H3 = sp.expand(H3.subs(G * Mtot / r, u).subs(r, G * Mtot / u))
    H3 = sp.expand((H3 + H3.subs(Delta, -Delta)) / 2)
    for n in range(20, 1, -2):
        H3 = sp.expand(H3.subs(Delta**n, (1 - 4 * nu) ** (n // 2)))
    while H3.has(Delta**2):
        H3 = sp.expand(H3.subs(Delta**2, 1 - 4 * nu))
    H3 = sp.expand(H3.subs(Delta, 0))

    subbanner("I.2 — Reduced COM Hamiltonian generated from the full generic-frame lift")
    print("H3_com / mu =")
    sp.pprint(H3)

    h = {
        1: sp.Rational(1, 128) * (-5 + 35 * nu - 70 * nu**2 + 35 * nu**3),
        2: sp.Integer(0),
        3: sp.Integer(0),
        4: sp.Integer(0),
        5: sp.Integer(0),
        6: sp.Rational(1, 16) * (-7 + 42 * nu - 53 * nu**2 - 5 * nu**3),
        7: sp.Rational(1, 16) * (2 - 3 * nu) * nu**2,
        8: sp.Rational(3, 16) * (1 - nu) * nu**2,
        9: -sp.Rational(5, 16) * nu**3,
        10: sp.Rational(1, 16) * (-27 + 136 * nu + 109 * nu**2),
        11: sp.Rational(1, 16) * (17 + 30 * nu) * nu,
        12: sp.Rational(1, 12) * (5 + 43 * nu) * nu,
        13: sp.Rational(1, 192) * (-600 + (3 * sp.pi**2 - 1340) * nu - 552 * nu**2),
        14: -sp.Rational(1, 64) * (340 + 3 * sp.pi**2 + 112 * nu) * nu,
        15: sp.Rational(1, 96) * (12 + (872 - 63 * sp.pi**2) * nu),
    }

    target = (
        h[1] * P2**4
        + u * (h[6] * P2**3 + h[7] * P2**2 * pr**2 + h[8] * P2 * pr**4 + h[9] * pr**6)
        + u**2 * (h[10] * P2**2 + h[11] * P2 * pr**2 + h[12] * pr**4)
        + u**3 * (h[13] * P2 + h[14] * pr**2)
        + u**4 * h[15]
    )

    banner("PART II — EXACT MATCH TO THE SOLVED GR COM HAMILTONIAN TARGET")
    expect_zero("Hamiltonian-level COM mismatch", H3 - target)

    subbanner("II.1 — Slot-by-slot coefficients")
    extracted = {
        1: sp.simplify(sp.expand(H3.coeff(P2, 4).subs({u: 0, pr: 0}))),
        6: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 1).coeff(P2, 3).subs(pr, 0))),
        7: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 1).coeff(P2, 2).coeff(pr, 2))),
        8: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 1).coeff(P2, 1).coeff(pr, 4))),
        9: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 1).coeff(pr, 6).subs(P2, 0))),
        10: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 2).coeff(P2, 2).subs(pr, 0))),
        11: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 2).coeff(P2, 1).coeff(pr, 2))),
        12: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 2).coeff(pr, 4).subs(P2, 0))),
        13: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 3).coeff(P2, 1).subs(pr, 0))),
        14: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 3).coeff(pr, 2).subs(P2, 0))),
        15: sp.simplify(sp.expand(sp.collect(H3, u).coeff(u, 4))),
    }
    for idx in [1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]:
        print(f"h{idx}^(lift) = {extracted[idx]}")
        expect_zero(f"slot h{idx}", extracted[idx] - h[idx])

    banner("PART III — THEOREM LEDGER")
    print("1. The generic-frame 3PN ordinary target imported earlier is consistent with the exact GR")
    print("   COM Hamiltonian target once one performs the correct generic-frame cubic Legendre transform")
    print("   before reducing to COM.")
    print("2. The earlier naive ordinary-Lagrangian COM mismatch is therefore a reduction-ordering")
    print("   artifact, not a failure of the generic-frame target itself.")
    print("3. The real remaining problem is now sharply identified: extract the ordinary generic-frame")
    print("   representative/canonical slice by matching to the Hamiltonian target, not by naive COM")
    print("   substitution.")

    return H3, h


if __name__ == "__main__":
    com_hamiltonian_from_generic_lift()
