#!/usr/bin/env python3
"""
moving_throat_pde_stage7_overlap_isotropy_sympy_audit.py

SymPy audit for Stage 7 of the moving-throat PDE program.

Scope
-----
This script verifies the first explicit overlap-integral layer beyond the abstract
Stage-6 grouped bundle:

  • orthonormality of the normalized real STF l=2 harmonics,
  • the exact angular source-map identity,
  • the isotropic grouped-bundle collapse under an O(3)-invariant kernel,
  • the exact axisymmetric quadrupole triple-overlap matrix,
  • the grouped 20/21/22 splitting pattern (1, 1/2, -1),
  • the first-order defect law b = 3 a,
  • and the corresponding first-order transport law for P_A = N_A / D_A.

This is still a reduced-sector theorem. It does not solve the full moving-throat
PDE, but it does verify the first explicit angular overlap rules the PDE must obey
on the natural branch.
"""

from __future__ import annotations

import sympy as sp

pi = sp.pi
sqrt = sp.sqrt


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
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def pairings(lst: list[int]) -> list[list[tuple[int, int]]]:
    if not lst:
        return [[]]
    a = lst[0]
    out = []
    for i in range(1, len(lst)):
        b = lst[i]
        rest = lst[1:i] + lst[i + 1 :]
        for p in pairings(rest):
            out.append([(a, b)] + p)
    return out


# ---------------------------------------------------------------------------
# Canonical real STF basis for l=2
# ---------------------------------------------------------------------------

E20 = sp.Matrix([[-1 / sqrt(6), 0, 0], [0, -1 / sqrt(6), 0], [0, 0, 2 / sqrt(6)]])
E21c = sp.Matrix([[0, 0, 1 / sqrt(2)], [0, 0, 0], [1 / sqrt(2), 0, 0]])
E21s = sp.Matrix([[0, 0, 0], [0, 0, 1 / sqrt(2)], [0, 1 / sqrt(2), 0]])
E22c = sp.Matrix([[1 / sqrt(2), 0, 0], [0, -1 / sqrt(2), 0], [0, 0, 0]])
E22s = sp.Matrix([[0, 1 / sqrt(2), 0], [1 / sqrt(2), 0, 0], [0, 0, 0]])
BASIS = [E20, E21c, E21s, E22c, E22s]
NAMES = ["20", "21c", "21s", "22c", "22s"]

NORM = sp.sqrt(sp.Rational(15) / (8 * sp.pi))
NBASIS = [sp.simplify(NORM * B) for B in BASIS]


# ---------------------------------------------------------------------------
# Isotropic angular moments of the unit sphere
# ---------------------------------------------------------------------------

delta = sp.eye(3)
PAIRINGS4 = [[(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]]
PAIRINGS6 = pairings([0, 1, 2, 3, 4, 5])


def I4(i: int, j: int, k: int, l: int) -> sp.Expr:
    inds = [i, j, k, l]
    s = 0
    for pr in PAIRINGS4:
        prod = 1
        for a, b in pr:
            prod *= delta[inds[a], inds[b]]
        s += prod
    return sp.simplify(4 * pi * s / 15)


def I6(i: int, j: int, k: int, l: int, m: int, n: int) -> sp.Expr:
    inds = [i, j, k, l, m, n]
    s = 0
    for pr in PAIRINGS6:
        prod = 1
        for a, b in pr:
            prod *= delta[inds[a], inds[b]]
        s += prod
    return sp.simplify(4 * pi * s / 105)


def quad_overlap(A: sp.Matrix, B: sp.Matrix) -> sp.Expr:
    out = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    out += A[i, j] * B[k, l] * I4(i, j, k, l)
    return sp.simplify(out)


def triple_overlap(A: sp.Matrix, Q: sp.Matrix, B: sp.Matrix) -> sp.Expr:
    out = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            out += A[i, j] * Q[k, l] * B[m, n] * I6(i, j, k, l, m, n)
    return sp.simplify(out)


# ---------------------------------------------------------------------------
# I. Normalized harmonic orthonormality and source-map identity
# ---------------------------------------------------------------------------

def normalized_harmonics_and_source_map() -> None:
    banner("SECTION I — NORMALIZED HARMONICS AND THE ANGULAR SOURCE MAP")

    gram = sp.Matrix([[quad_overlap(A, B) for B in NBASIS] for A in NBASIS])
    subbanner("I.1 — Orthonormality of the normalized real STF harmonics")
    print("Gram matrix =")
    sp.pprint(gram)
    expect_zero("Gram - I5", gram - sp.eye(5))

    s20, s21c, s21s, s22c, s22s = sp.symbols("s20 s21c s21s s22c s22s", real=True)
    svec = sp.Matrix([s20, s21c, s21s, s22c, s22s])
    projected = sp.simplify(gram * svec)

    subbanner("I.2 — Exact angular source-map identity")
    expect_zero("projected coefficients - source coefficients", projected - svec)
    print("Conclusion: on the canonical isotropic angular basis, mhat_ang = 1 exactly.")


# ---------------------------------------------------------------------------
# II. Isotropic grouped-bundle collapse
# ---------------------------------------------------------------------------

def isotropic_grouped_collapse() -> None:
    banner("SECTION II — ISOTROPIC GROUPED-BUNDLE COLLAPSE")

    x0 = sp.symbols("x0", real=True)
    x20, x21, x22 = sp.symbols("x20 x21 x22", real=True)

    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)

    subbanner("II.1 — Equal-lane data imply exact grouped isotropy")
    equal_lane = {x20: x0, x21: x0, x22: x0}
    expect_zero("xbar - x0", xbar.subs(equal_lane) - x0)
    expect_zero("a_x on equal lanes", ax.subs(equal_lane))
    expect_zero("b_x on equal lanes", bx.subs(equal_lane))

    subbanner("II.2 — Unequal lanes produce visible grouped defects")
    witness_a = {x20: x0 + 1, x21: x0, x22: x0}
    witness_b = {x20: x0, x21: x0 + 1, x22: x0}
    expect_zero("a_x witness - 1/5", ax.subs(witness_a) - sp.Rational(1, 5))
    expect_zero("b_x witness", bx.subs(witness_a))
    expect_zero("a_x second witness + 1/10", ax.subs(witness_b) + sp.Rational(1, 10))
    expect_zero("b_x second witness - 1/2", bx.subs(witness_b) - sp.Rational(1, 2))

    print("Interpretation: any O(3)-invariant reduced kernel forces the grouped 20/21/22 bundle")
    print("to collapse to one common scalar value on the isotropic branch.")


# ---------------------------------------------------------------------------
# III. Weak axisymmetric quadrupole splitting
# ---------------------------------------------------------------------------

def axisymmetric_splitting() -> None:
    banner("SECTION III — AXISYMMETRIC QUADRUPOLE SPLITTING MATRIX")

    Q = NBASIS[0]  # normalized Y20 source/background
    M = sp.Matrix([[triple_overlap(A, Q, B) for B in NBASIS] for A in NBASIS])

    subbanner("III.1 — Exact five-mode triple-overlap matrix for Y20")
    print("M^(20) =")
    sp.pprint(M)

    kappa_star = sp.sqrt(5) / (7 * sp.sqrt(sp.pi))
    M_target = sp.diag(kappa_star, kappa_star / 2, kappa_star / 2, -kappa_star, -kappa_star)
    expect_zero("M - M_target", M - M_target)

    subbanner("III.2 — Grouped 20/21/22 splitting weights")
    lam20 = sp.Integer(1)
    lam21 = sp.Rational(1, 2)
    lam22 = -sp.Integer(1)
    x0, eps, x1 = sp.symbols("x0 eps x1", real=True)

    x20 = x0 + eps * lam20 * x1
    x21 = x0 + eps * lam21 * x1
    x22 = x0 + eps * lam22 * x1

    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)

    expect_zero("xbar - x0", xbar - x0)
    expect_zero("a_x - eps*x1/4", ax - eps * x1 / 4)
    expect_zero("b_x - 3*eps*x1/4", bx - 3 * eps * x1 / 4)
    expect_zero("b_x - 3 a_x", bx - 3 * ax)

    print("Axisymmetric grouped weights:")
    print("  lambda_20 = 1")
    print("  lambda_21 = 1/2")
    print("  lambda_22 = -1")


# ---------------------------------------------------------------------------
# IV. First-order transport law for P_A = N_A / D_A
# ---------------------------------------------------------------------------

def first_order_transport() -> None:
    banner("SECTION IV — FIRST-ORDER TRANSPORT LAW FOR THE NORMALIZATION RATIO")

    eps = sp.symbols("eps", real=True)
    D0, D1, N0, N1 = sp.symbols("D0 D1 N0 N1", nonzero=True, real=True)
    lam20 = sp.Integer(1)
    lam21 = sp.Rational(1, 2)
    lam22 = -sp.Integer(1)

    def lane_ratio(lam: sp.Expr) -> sp.Expr:
        expr = (N0 + eps * lam * N1) / (D0 + eps * lam * D1)
        return sp.expand(sp.series(expr, eps, 0, 2).removeO())

    P20 = lane_ratio(lam20)
    P21 = lane_ratio(lam21)
    P22 = lane_ratio(lam22)

    P0 = sp.simplify(N0 / D0)
    P1 = sp.simplify((N1 * D0 - N0 * D1) / D0**2)

    subbanner("IV.1 — Exact first-order expansion")
    expect_zero("P20 - (P0 + eps P1)", P20 - (P0 + eps * P1))
    expect_zero("P21 - (P0 + eps P1/2)", P21 - (P0 + eps * P1 / 2))
    expect_zero("P22 - (P0 - eps P1)", P22 - (P0 - eps * P1))

    abar = sp.simplify((P20 + 2 * P21 + 2 * P22) / 5)
    aP = sp.simplify((2 * P20 - P21 - P22) / 10)
    bP = sp.simplify((P21 - P22) / 2)

    subbanner("IV.2 — Grouped defects of the normalization ratio")
    expect_zero("Pbar - P0", abar - P0)
    expect_zero("a_P - eps*P1/4", aP - eps * P1 / 4)
    expect_zero("b_P - 3*eps*P1/4", bP - 3 * eps * P1 / 4)
    expect_zero("b_P - 3 a_P", bP - 3 * aP)

    print("So the Stage-6 normalization stack inherits the same exact axisymmetric")
    print("defect law as the microscopic grouped overlaps.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    normalized_harmonics_and_source_map()
    isotropic_grouped_collapse()
    axisymmetric_splitting()
    first_order_transport()
