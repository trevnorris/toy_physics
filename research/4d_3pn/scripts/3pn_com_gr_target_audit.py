#!/usr/bin/env python3
"""
3pn_com_gr_target_audit.py

Import the exact GR 3PN center-of-mass ordinary ADM Hamiltonian target, solve the
corresponding ordinary-Lagrangian coefficients using the previously derived diagonal
COM linear map, and isolate the residual beyond the carried self/static seed.

Primary sources used for the target formulas
-------------------------------------------
1. Memmesheimer, Gopakumar, Schafer, gr-qc/0407049, Eq. (7d): fully determined
   reduced ordinary 3PN ADM Hamiltonian in the center-of-mass frame.
2. Jain et al., arXiv:2211.15580, Eq. (4.1): the standard 15-slot isotropic,
   time-translation invariant COM Hamiltonian basis.

What this script checks
-----------------------
1. The imported GR target coefficients h_i^GR reproduce the exact one-body Schwarzschild
   gate in the strict test-mass limit nu -> 0.
2. The previously derived diagonal inverse map yields exact COM ordinary coefficients l_i^GR.
3. The carried self/static seed is recovered exactly in the strict one-body limit.
4. The genuine comparable-mass residuals Delta h_i and Delta l_i vanish as nu -> 0.
5. Because the COM map is diagonal with the same lower-order feedback for target and seed,
   the ordinary residual is exactly the negative Hamiltonian residual slot by slot:
       Delta l_i = - Delta h_i.
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


nu = sp.symbols("nu", real=True)
pi = sp.pi


# ---------------------------------------------------------------------------
# Exact GR 3PN COM Hamiltonian target, compiled in the standard 15-slot basis.
# ---------------------------------------------------------------------------

def gr_target_h() -> dict[int, sp.Expr]:
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


# ---------------------------------------------------------------------------
# Exact diagonal COM inverse map l_i = F_i(nu) - h_i.
# ---------------------------------------------------------------------------

def inverse_map_from_h(h: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
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


# ---------------------------------------------------------------------------
# Carried self/static seed in the same COM ordinary and Hamiltonian bases.
# ---------------------------------------------------------------------------

def carried_seed_l() -> dict[int, sp.Expr]:
    seed: dict[int, sp.Expr] = {}

    seed[1] = sp.Rational(5, 128) - sp.Rational(35, 128) * nu + sp.Rational(35, 64) * nu**2 - sp.Rational(35, 128) * nu**3
    seed[2] = 0
    seed[3] = 0
    seed[4] = 0
    seed[5] = 0

    seed[6] = sp.Rational(11, 16) - sp.Rational(33, 8) * nu + sp.Rational(99, 16) * nu**2 - sp.Rational(11, 8) * nu**3
    seed[7] = 0
    seed[8] = 0
    seed[9] = 0

    seed[10] = sp.Rational(47, 16) - sp.Rational(235, 16) * nu + sp.Rational(235, 16) * nu**2
    seed[11] = 0
    seed[12] = 0

    seed[13] = sp.Rational(13, 8) - sp.Rational(13, 2) * nu + sp.Rational(13, 4) * nu**2
    seed[14] = 0
    seed[15] = -sp.Rational(1, 8) + sp.Rational(3, 8) * nu

    return {i: sp.simplify(seed[i]) for i in range(1, 16)}


def carried_seed_h() -> dict[int, sp.Expr]:
    seed: dict[int, sp.Expr] = {}

    seed[1] = -sp.Rational(5, 128) + sp.Rational(59, 128) * nu - sp.Rational(119, 64) * nu**2 + sp.Rational(323, 128) * nu**3
    seed[2] = 0
    seed[3] = 0
    seed[4] = 0
    seed[5] = 0

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


# ---------------------------------------------------------------------------
# Checks and ledgers
# ---------------------------------------------------------------------------

def verify_map(target_h: dict[int, sp.Expr], target_l: dict[int, sp.Expr]) -> None:
    """Re-apply the diagonal map h_i = F_i(nu) - l_i."""
    reconstructed = inverse_map_from_h({i: sp.Symbol(f"h{i}") for i in range(1, 16)})
    # Convert the symbolic inverse-map formulas back into h(F,l) checks by direct substitution.
    # Since the map is diagonal, it is enough to verify the stored explicit formulas slot by slot.
    back = {
        1: sp.simplify(sp.Rational(3, 16) * nu - sp.Rational(21, 16) * nu**2 + sp.Rational(9, 4) * nu**3 - target_l[1]),
        2: sp.simplify(-target_l[2]),
        3: sp.simplify(-target_l[3]),
        4: sp.simplify(-target_l[4]),
        5: sp.simplify(-target_l[5]),
        6: sp.simplify(sp.Rational(1, 4) + sp.Rational(7, 8) * nu - sp.Rational(35, 8) * nu**2 - sp.Rational(21, 4) * nu**3 - target_l[6]),
        7: sp.simplify(sp.Rational(11, 8) * nu**2 - sp.Rational(9, 2) * nu**3 - target_l[7]),
        8: sp.simplify(sp.Rational(3, 4) * nu**2 - sp.Rational(9, 4) * nu**3 - target_l[8]),
        9: sp.simplify(-target_l[9]),
        10: sp.simplify(sp.Rational(5, 4) + sp.Rational(15, 8) * nu + sp.Rational(123, 8) * nu**2 + sp.Rational(13, 4) * nu**3 - target_l[10]),
        11: sp.simplify(sp.Rational(7, 8) * nu + sp.Rational(41, 8) * nu**2 + sp.Rational(31, 4) * nu**3 - target_l[11]),
        12: sp.simplify(sp.Rational(9, 2) * nu**2 + 4 * nu**3 - target_l[12]),
        13: sp.simplify(-sp.Rational(3, 2) - sp.Rational(59, 4) * nu - sp.Rational(25, 4) * nu**2 - sp.Rational(1, 2) * nu**3 - target_l[13]),
        14: sp.simplify(sp.Rational(7, 4) * nu - sp.Rational(31, 4) * nu**2 - sp.Rational(7, 2) * nu**3 - target_l[14]),
        15: sp.simplify(-target_l[15]),
    }
    for i in range(1, 16):
        expect_zero(f"map check h{i}", back[i] - target_h[i])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    banner("PART I — IMPORT THE EXACT GR 3PN COM HAMILTONIAN TARGET")
    target_h = gr_target_h()
    for i in range(1, 16):
        print(f"h{i}^(GR) = {sp.factor(target_h[i])}")

    banner("PART II — SOLVE THE EXACT COM ORDINARY-LAGRANGIAN COEFFICIENTS")
    target_l = inverse_map_from_h(target_h)
    for i in range(1, 16):
        print(f"l{i}^(GR) = {sp.expand(target_l[i])}")

    subbanner("II.1 — Re-apply the diagonal map")
    verify_map(target_h, target_l)

    banner("PART III — ONE-BODY GATE AND SEED CHECKS")
    seed_h = carried_seed_h()
    seed_l = carried_seed_l()

    # Strict one-body values must agree with the exact one-body target and with the carried seed.
    one_body_slots = [1, 6, 10, 13, 15]
    for i in one_body_slots:
        expect_zero(f"one-body h{i} target - seed", sp.simplify(target_h[i].subs(nu, 0) - seed_h[i].subs(nu, 0)))
        expect_zero(f"one-body l{i} target - seed", sp.simplify(target_l[i].subs(nu, 0) - seed_l[i].subs(nu, 0)))
    for i in [2, 3, 4, 5, 7, 8, 9, 11, 12, 14]:
        expect_zero(f"one-body h{i}", sp.simplify(target_h[i].subs(nu, 0)))
        expect_zero(f"one-body l{i}", sp.simplify(target_l[i].subs(nu, 0) - seed_l[i].subs(nu, 0)))

    banner("PART IV — GENUINE COMPARABLE-MASS RESIDUALS")
    delta_h = {i: sp.simplify(target_h[i] - seed_h[i]) for i in range(1, 16)}
    delta_l = {i: sp.simplify(target_l[i] - seed_l[i]) for i in range(1, 16)}

    subbanner("IV.1 — Hamiltonian residual Delta h_i = h_i^(GR) - h_i^(seed)")
    for i in range(1, 16):
        print(f"Delta h{i} = {sp.factor(delta_h[i])}")

    subbanner("IV.2 — Ordinary-Lagrangian residual Delta l_i = l_i^(GR) - l_i^(seed)")
    for i in range(1, 16):
        print(f"Delta l{i} = {sp.factor(delta_l[i])}")

    subbanner("IV.3 — Residuals are pure comparable-mass data and satisfy Delta l_i = -Delta h_i")
    for i in range(1, 16):
        expect_zero(f"nu -> 0 residual h{i}", sp.simplify(delta_h[i].subs(nu, 0)))
        expect_zero(f"nu -> 0 residual l{i}", sp.simplify(delta_l[i].subs(nu, 0)))
        expect_zero(f"Delta l{i} + Delta h{i}", sp.simplify(delta_l[i] + delta_h[i]))

    banner("PART V — QUICK READOUTS")
    nu_eq = sp.Rational(1, 4)
    print("Equal-mass target coefficients (nu = 1/4):")
    for i in range(1, 16):
        print(f"  h{i}(1/4) = {sp.simplify(target_h[i].subs(nu, nu_eq))}")

    print("\nEqual-mass ordinary coefficients (nu = 1/4):")
    for i in range(1, 16):
        print(f"  l{i}(1/4) = {sp.simplify(target_l[i].subs(nu, nu_eq))}")

    print("\nInterpretation:")
    print("  - The exact GR target is now imported into the 15-slot COM basis.")
    print("  - The diagonal map solves the ordinary COM coefficients immediately.")
    print("  - The genuine comparable-mass content is isolated by Delta h_i or equivalently Delta l_i.")
    print("  - The next remaining task is to lift this solved COM answer back to the generic frame,")
    print("    or import a generic-frame target block directly into the Phase-C residual solver.")


if __name__ == "__main__":
    main()
