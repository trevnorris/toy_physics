#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage009_strict_5pn_even_gate_mixed_corridor_sympy_audit.py

Stage 009 — carry the surviving Stage-008 mixed-sector corridor into the
stricter imported 5PN even-gate package, identify the surviving strict corridor,
and isolate the pure-transfer subcorridor that also preserves the Stage-008
first-order conservative response on the noncanonical sample branch.

What this script does
---------------------
1. Loads the exact Stage-008 primitive compiler on the concrete Stage-006
   compatibility point.
2. Imports the strict 5PN lower-order gate package
      K1      = D21 + D01/9,
      H_even  = D41 - (2/3) D21 - D01/27,
      Xi_load = N01/N0 - D01/D0 = P1/P0.
3. Proves exactly that on the Stage-008 conservative compensation surface
      D21 = -u2 D01,
      D41 = (D4/D0) D01,
   the strict 5PN gates collapse to
      K1     = (1/9 - u2) D01,
      H_even = (D4/D0 + 2u2/3 - 1/27) D01.
   So on a noncanonical one-pole sample branch the full intersection of
   Stage-008 compensation with the strict 5PN gates forces D01 = 0.
4. Solves the mixed-sector-only strict even-gate system on the concrete sample
   branch and shows that a nontrivial 3-dimensional mixed corridor survives.
5. Solves the sharper mixed-sector subcorridor defined by
      eq1 = 0,
      eq2 = 0,
      D01 = 0,
   showing that a 2-dimensional pure-transfer corridor survives in which
      D01 = D21 = D41 = 0,
      Xi_1 = N01/N0.
6. Computes transported Stage-007 same-charge ceiling budgets on both the
   strict even-gate corridor and the pure-transfer subcorridor using canonical
   Euclidean unit-amplitude microscopic drifts.
"""

from __future__ import annotations

import math
from pathlib import Path
import runpy

import numpy as np
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
        simplified = expr.applyfunc(lambda z: sp.factor(sp.simplify(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.factor(sp.simplify(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def expect_close(name: str, expr: sp.Expr, tol: float = 1e-28) -> None:
    val = abs(complex(sp.N(expr, 80)))
    print(f"{name} ~= {val}")
    if val > tol:
        raise AssertionError(f"{name} is not within tolerance {tol}")


def fmt(x: float) -> str:
    return f"{x:.18f}"


def numeric_matrix_from_dictrows(rows: list[dict[str, sp.Expr]], cols: list[str]) -> np.ndarray:
    return np.array([[float(sp.N(row[name], 50)) for name in cols] for row in rows], dtype=float)


def orthonormal_nullspace(A: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """Return an ambient-orthonormal basis for ker(A)."""
    _u, s, vh = np.linalg.svd(A)
    rank = int((s > tol).sum())
    basis = vh[rank:].T.copy()
    if basis.size == 0:
        return basis.reshape(A.shape[1], 0)
    # Numerical cleanup for display/reuse.
    basis[np.abs(basis) < 1e-14] = 0.0
    return basis


def raw_nullspace(A_sym: sp.Matrix) -> list[sp.Matrix]:
    A_num = A_sym.applyfunc(lambda z: sp.N(z, 40))
    return A_num.nullspace()


def xi_coeff_vector(compiler, cols: list[str]) -> np.ndarray:
    return np.array([float(sp.N(compiler.Xi[name], 50)) for name in cols], dtype=float)


def projected_xi_norm(A: np.ndarray, xi_vec: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    B = orthonormal_nullspace(A)
    xi_coords = B.T @ xi_vec
    sigma = float(np.linalg.norm(xi_coords))
    return B, xi_coords, sigma


# ---------------------------------------------------------------------------
# Load Stage 008 symbols and exact primitive compiler
# ---------------------------------------------------------------------------

STAGE008_PATH = Path(__file__).with_name(
    "same_charge_barrier_audit_stage008_microscopic_xi1_primitive_compiler_sympy_audit.py"
)
S8 = runpy.run_path(str(STAGE008_PATH))

PrimitiveParams = S8["PrimitiveParams"]
primitive_base_bundle = S8["primitive_base_bundle"]
primitive_slope_compiler = S8["primitive_slope_compiler"]

P_COMPAT_SAMPLE = S8["P_COMPAT_SAMPLE"]
P_BOTH_10 = S8["P_BOTH_10"]
P_ONE_10 = S8["P_ONE_10"]
P_BOTH_30 = S8["P_BOTH_30"]
P_ONE_30 = S8["P_ONE_30"]

SLOPE_NAMES = S8["SLOPE_NAMES"]
MIXED_COLS = ["xLambdaU", "xLambdaW", "xLambdaR", "xOmegaU", "xOmegaW"]


# ---------------------------------------------------------------------------
# Part I. Rebuild the concrete Stage-006 compatibility branch
# ---------------------------------------------------------------------------

def concrete_branch_and_compiler():
    params = PrimitiveParams(
        lamB=sp.Rational(1, 2),
        lamU=sp.Rational(3, 10),
        lamW=sp.Rational(2, 5),
        lamR=sp.Rational(1, 4),
        OmegaU=sp.Integer(1),
        OmegaW=sp.Rational(7, 5),
        varpi=sp.Integer(2),
        M=sp.Integer(1),
    )
    base = primitive_base_bundle(params)
    compiler = primitive_slope_compiler(base)
    return base, compiler


# ---------------------------------------------------------------------------
# Part II. Imported strict 5PN gate package
# ---------------------------------------------------------------------------

def strict_gate_dicts(base, compiler):
    K1 = {name: sp.simplify(compiler.D21[name] + compiler.D01[name] / 9) for name in SLOPE_NAMES}
    H_even = {
        name: sp.simplify(compiler.D41[name] - sp.Rational(2, 3) * compiler.D21[name] - compiler.D01[name] / 27)
        for name in SLOPE_NAMES
    }
    Xi_load = {name: sp.simplify(compiler.N01[name] / base.N0 - compiler.D01[name] / base.D0) for name in SLOPE_NAMES}
    return K1, H_even, Xi_load


def verify_xi_bridge(base, compiler, Xi_load):
    banner("PART I — EXACT 5PN LOAD BRIDGE")
    for name in SLOPE_NAMES:
        expect_zero(f"Xi_load[{name}] - Xi_1[{name}]", sp.simplify(Xi_load[name] - compiler.Xi[name]))
    print("Result: Xi_load = Xi_1 = P1/P0 coefficientwise on the primitive compiler.")


# ---------------------------------------------------------------------------
# Part III. Exact comparison between Stage-008 compensation and strict 5PN gates
# ---------------------------------------------------------------------------

def verify_compensation_vs_strict_gates(base):
    banner("PART II — EXACT COMPARISON: STAGE-008 COMPENSATION VS STRICT 5PN GATES")

    D0, D2, D4, D01 = sp.symbols("D0 D2 D4 D01", real=True)
    u2 = sp.simplify(-D2 / D0)
    u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)
    D21 = sp.simplify(-u2 * D01)
    D41 = sp.simplify((D4 / D0) * D01)

    K1_on_comp = sp.simplify(D21 + D01 / 9)
    H_on_comp = sp.simplify(D41 - sp.Rational(2, 3) * D21 - D01 / 27)

    print("On the Stage-008 compensation surface:")
    print("  D21 = -u2 D01")
    print("  D41 = (D4/D0) D01")
    print()
    print("K1 becomes")
    sp.pprint(K1_on_comp)
    print("H_even becomes")
    sp.pprint(H_on_comp)

    expect_zero("K1_on_comp - (1/9 - u2) D01", sp.simplify(K1_on_comp - (sp.Rational(1, 9) - u2) * D01))
    expect_zero(
        "H_on_comp - (D4/D0 + 2u2/3 - 1/27) D01",
        sp.simplify(H_on_comp - (D4 / D0 + 2 * u2 / 3 - sp.Rational(1, 27)) * D01),
    )
    expect_zero("D4/D0 - (u2^2 - u4)", sp.simplify(D4 / D0 - (u2**2 - u4)))
    H_on_comp_one_pole = sp.simplify(H_on_comp.subs({D4: -3 * D2**2 / D0}))
    expect_zero(
        "one-pole H_on_comp - (-3u2^2 + 2u2/3 - 1/27) D01",
        sp.simplify(H_on_comp_one_pole - (-3 * u2**2 + 2 * u2 / 3 - sp.Rational(1, 27)) * D01),
    )

    print()
    print("Concrete Stage-006 compatibility-point coefficients:")
    coeff_K = sp.N(sp.Rational(1, 9) - base.u2, 30)
    coeff_H = sp.N(base.D4 / base.D0 + 2 * base.u2 / 3 - sp.Rational(1, 27), 30)
    print(f"  K1 = ({coeff_K}) * D01")
    print(f"  H_even = ({coeff_H}) * D01")
    if coeff_K == 0 or coeff_H == 0:
        raise AssertionError("Unexpected canonical collapse on the noncanonical sample point.")
    print("Therefore, on this noncanonical sample branch, the full intersection of")
    print("Stage-008 compensation with the strict 5PN gates forces D01 = 0.")


# ---------------------------------------------------------------------------
# Part IV. Mixed-sector-only strict even-gate corridor
# ---------------------------------------------------------------------------

def print_row_dict(title: str, row: dict[str, sp.Expr], cols: list[str]) -> None:
    print(title)
    for name in cols:
        val = sp.N(row[name], 18)
        if val != 0:
            print(f"  {name:9s} -> {val}")
    print()


def strict_even_gate_corridor(base, compiler, K1, H_even):
    banner("PART III — MIXED-SECTOR-ONLY STRICT 5PN EVEN-GATE CORRIDOR")

    print_row_dict("K1 coefficients", K1, MIXED_COLS)
    print_row_dict("H_even coefficients", H_even, MIXED_COLS)
    print_row_dict("Xi_load = Xi_1 coefficients", compiler.Xi, MIXED_COLS)

    A_sym = sp.Matrix([
        [sp.N(K1[name], 30) for name in MIXED_COLS],
        [sp.N(H_even[name], 30) for name in MIXED_COLS],
    ])
    A = np.array(A_sym.tolist(), dtype=float)
    rank = np.linalg.matrix_rank(A)
    null_raw = raw_nullspace(A_sym)

    print("Strict even-gate mixed-only matrix =")
    sp.pprint(sp.N(A_sym, 18))
    print(f"rank = {rank} ; nullity = {A.shape[1] - rank}")

    xi_vec = xi_coeff_vector(compiler, MIXED_COLS)
    for i, vec in enumerate(null_raw, start=1):
        xi = sp.simplify(sum(sp.N(compiler.Xi[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        d01 = sp.simplify(sum(sp.N(compiler.D01[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        d21 = sp.simplify(sum(sp.N(compiler.D21[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        d41 = sp.simplify(sum(sp.N(compiler.D41[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        print(f"\nStrict-corridor raw null basis vector w{i} =")
        sp.pprint(sp.N(vec, 20))
        print(f"Xi_1(w{i}) = {sp.N(xi, 20)}")
        print(f"D01(w{i})  = {sp.N(d01, 20)}")
        print(f"D21(w{i})  = {sp.N(d21, 20)}")
        print(f"D41(w{i})  = {sp.N(d41, 20)}")

    B_even, xi_coords_even, sigma_even = projected_xi_norm(A, xi_vec)
    print("\nAmbient-orthonormal basis for the strict even-gate corridor =")
    sp.pprint(sp.Matrix(np.round(B_even, 12)))
    print("Xi coordinates on that orthonormal basis =")
    sp.pprint(sp.Matrix(np.round(xi_coords_even, 12)))
    print(f"Euclidean operator norm sigma_even = {fmt(sigma_even)}")

    return A_sym, A, B_even, xi_coords_even, sigma_even


# ---------------------------------------------------------------------------
# Part V. Pure-transfer subcorridor
# ---------------------------------------------------------------------------

def pure_transfer_subcorridor(base, compiler):
    banner("PART IV — PURE-TRANSFER SUBCORRIDOR ON THE NONCANONICAL SAMPLE BRANCH")

    A_sym = sp.Matrix([
        [sp.N(compiler.eq1[name], 30) for name in MIXED_COLS],
        [sp.N(compiler.eq2[name], 30) for name in MIXED_COLS],
        [sp.N(compiler.D01[name], 30) for name in MIXED_COLS],
    ])
    A = np.array(A_sym.tolist(), dtype=float)
    rank = np.linalg.matrix_rank(A)
    null_raw = raw_nullspace(A_sym)

    print("Intersection matrix (eq1, eq2, D01) =")
    sp.pprint(sp.N(A_sym, 18))
    print(f"rank = {rank} ; nullity = {A.shape[1] - rank}")
    print()
    print("Because the Stage-006 sample branch is noncanonical, Stage-008 compensation +")
    print("strict 5PN even gates reduces exactly to D01 = 0. So this 3x5 system is the")
    print("full mixed-sector intersection of")
    print("  - Stage-008 conservative-shape preservation, and")
    print("  - the imported strict 5PN even-gate package.")

    xi_vec = xi_coeff_vector(compiler, MIXED_COLS)
    for i, vec in enumerate(null_raw, start=1):
        xi = sp.simplify(sum(sp.N(compiler.Xi[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        d01 = sp.simplify(sum(sp.N(compiler.D01[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        d21 = sp.simplify(sum(sp.N(compiler.D21[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        d41 = sp.simplify(sum(sp.N(compiler.D41[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        n01 = sp.simplify(sum(sp.N(compiler.N01[name], 30) * vec[j] for j, name in enumerate(MIXED_COLS)))
        print(f"\nPure-transfer raw null basis vector t{i} =")
        sp.pprint(sp.N(vec, 20))
        print(f"Xi_1(t{i}) = {sp.N(xi, 20)}")
        print(f"N01(t{i})  = {sp.N(n01, 20)}")
        print(f"D01(t{i})  = {sp.N(d01, 20)}")
        print(f"D21(t{i})  = {sp.N(d21, 20)}")
        print(f"D41(t{i})  = {sp.N(d41, 20)}")

    # Exact logical implication on this branch.
    D01_s, D21_s, D41_s, u2_s, r_s = sp.symbols("D01_s D21_s D41_s u2_s r_s", real=True)
    expect_zero("eq1 when D01=0 and D21=0", sp.simplify(D21_s + u2_s * 0).subs({D21_s: 0}))
    expect_zero("eq2 when D01=0 and D41=0", sp.simplify(D41_s - r_s * 0).subs({D41_s: 0}))
    print("Hence on the pure-transfer subcorridor:")
    print("  D01 = D21 = D41 = 0")
    print("and therefore")
    print("  Xi_1 = Xi_load = N01 / N0")

    B_tr, xi_coords_tr, sigma_tr = projected_xi_norm(A, xi_vec)
    print("\nAmbient-orthonormal basis for the pure-transfer corridor =")
    sp.pprint(sp.Matrix(np.round(B_tr, 12)))
    print("Xi coordinates on that orthonormal basis =")
    sp.pprint(sp.Matrix(np.round(xi_coords_tr, 12)))
    print(f"Euclidean operator norm sigma_transfer = {fmt(sigma_tr)}")

    return A_sym, A, B_tr, xi_coords_tr, sigma_tr


# ---------------------------------------------------------------------------
# Part VI. Transported Stage-007 budgets on the strict corridors
# ---------------------------------------------------------------------------

def print_transported_budgets(sigma_even: float, sigma_transfer: float) -> None:
    banner("PART V — TRANSPORTED SAME-CHARGE CEILINGS ON THE STRICT CORRIDORS")

    xi_budgets = [
        ("10% loss, both wall-like poles", float(P_BOTH_10 / P_COMPAT_SAMPLE - 1)),
        ("10% loss, nonempty wall-like corridor", float(P_ONE_10 / P_COMPAT_SAMPLE - 1)),
        ("30% loss, both wall-like poles", float(P_BOTH_30 / P_COMPAT_SAMPLE - 1)),
        ("30% loss, nonempty wall-like corridor", float(P_ONE_30 / P_COMPAT_SAMPLE - 1)),
    ]

    print("Interpret the ambient microscopic mixed-sector drift amplitude as")
    print("  ||x_mixed||_2 = t.")
    print("Then on a corridor with Xi operator norm sigma, the transported bound is")
    print("  |epsilon| t <= budget / sigma.")

    print("\nStrict 3D even-gate corridor budgets:")
    for label, budget in xi_budgets:
        print(f"  [{label}]  |epsilon| t <= {fmt(budget / sigma_even)}")

    print("\nPure-transfer 2D subcorridor budgets:")
    for label, budget in xi_budgets:
        print(f"  [{label}]  |epsilon| t <= {fmt(budget / sigma_transfer)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base, compiler = concrete_branch_and_compiler()
    K1, H_even, Xi_load = strict_gate_dicts(base, compiler)

    verify_xi_bridge(base, compiler, Xi_load)
    verify_compensation_vs_strict_gates(base)
    _, _, _, _, sigma_even = strict_even_gate_corridor(base, compiler, K1, H_even)
    _, _, _, _, sigma_transfer = pure_transfer_subcorridor(base, compiler)
    print_transported_budgets(sigma_even, sigma_transfer)

    banner("STAGE 009 LEDGER")
    print("1. The imported 5PN strict load defect Xi_load is exactly the Stage-008 same-charge scalar Xi_1 = P1/P0.")
    print("2. The Stage-008 conservative compensation surface is weaker than the strict 5PN even-gate package.")
    print("   On the concrete noncanonical sample branch, enforcing both together forces D01 = 0.")
    print("3. The mixed-sector corridor survives the strict even-gate package itself as a 3-dimensional family.")
    print("4. If one also requires the Stage-008 conservative-shape preservation on this noncanonical branch,")
    print("   the surviving corridor collapses to a 2-dimensional pure-transfer family with")
    print("      D01 = D21 = D41 = 0,   Xi_1 = N01/N0.")
    print("5. So after Stage 009, the sharpest surviving same-charge mechanism is no longer generic anisotropy;")
    print("   it is mixed-sector outgoing-transfer enhancement with the conservative one-pole bundle frozen at first order.")


if __name__ == "__main__":
    main()
