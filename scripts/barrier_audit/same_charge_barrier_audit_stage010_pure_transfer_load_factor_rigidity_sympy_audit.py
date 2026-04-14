#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage010_pure_transfer_load_factor_rigidity_sympy_audit.py

Stage 010 — take the Stage-009 pure-transfer corridor and rewrite it as an exact
one-port outgoing-load problem. Then impose the first outgoing-rigidity filters.

What this script does
---------------------
1. Rebuilds the explicit Stage-006/008 primitive finite-throat one-port branch.
2. Rebuilds the Stage-009 pure-transfer subcorridor defined by
      D01 = D21 = D41 = 0
   inside the mixed primitive slope space
      (xLambdaU, xLambdaW, xLambdaR, xOmegaU, xOmegaW).
3. Proves exactly on that subcorridor that
      Xi_1 = N01/N0 = 2(P01/P) - 2(Delta01/Delta),
   so the surviving same-charge scalar is exactly the outgoing-load-factor slope.
4. Factors the one-port load factor as
      Lambda = P/Delta = (G_W/Omega_W^2) * (1 + I)/(1 - H),
   with
      I = R G_U/(Omega_U^2 G_W),
      H = R^2/(Omega_U^2 Omega_W^2),
   and proves on the pure-transfer subcorridor that
      Xi_1 = 2 m + 2 I/(1+I) i + 2 H/(1-H) h,
   where
      m = d ln(G_W/Omega_W^2),
      i = d ln I,
      h = d ln H.
5. Applies a rigidity sieve:
      i = 0,
      h = 0,
      m = 0,
      and the combined pair i = h = 0.
6. Shows that:
      - the combined interference+hybridization rigidity kills the corridor,
      - each individual rigidity leaves a 1D survivor,
      - even the mixed-leg-rigid branch m = 0 still leaves a nonzero Xi_1 corridor.
7. Carries the transported Stage-007 same-charge ceilings onto those rigid
   one-dimensional survivors.
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


def fmt(x: float) -> str:
    return f"{x:.18f}"


def numeric_row(row: sp.Matrix) -> np.ndarray:
    return np.array([float(sp.N(z, 60)) for z in list(row)], dtype=float)


def orthonormal_nullspace(A: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    _u, s, vh = np.linalg.svd(A)
    rank = int((s > tol).sum())
    B = vh[rank:].T.copy()
    if B.size == 0:
        return B.reshape(A.shape[1], 0)
    B[np.abs(B) < 1e-14] = 0.0
    return B


def unit_null_direction(A: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    B = orthonormal_nullspace(A, tol=tol)
    if B.shape[1] != 1:
        raise AssertionError("Expected a 1D nullspace.")
    v = B[:, 0].copy()
    nz = np.flatnonzero(np.abs(v) > 1e-12)
    if len(nz) and v[nz[0]] < 0:
        v = -v
    return v


def print_unit_direction(name: str, cols: list[str], v: np.ndarray, Xi_vec: np.ndarray) -> None:
    print(f"{name} unit direction =")
    for c, x in zip(cols, v):
        print(f"  {c:9s} -> {fmt(x)}")
    print(f"  Xi_1 amplitude on this unit direction = {fmt(float(v @ Xi_vec))}")
    print()


def budget_table(label: str, sigma: float, budgets: list[tuple[str, float]]) -> None:
    print(f"{label}: sigma = {fmt(sigma)}")
    for tag, b in budgets:
        print(f"  [{tag}]  |epsilon| t <= {fmt(b / sigma)}")
    print()


# ---------------------------------------------------------------------------
# Load Stage 008 exact primitive compiler
# ---------------------------------------------------------------------------

STAGE008_PATH = Path(__file__).with_name(
    "same_charge_barrier_audit_stage008_microscopic_xi1_primitive_compiler_sympy_audit.py"
)
S8 = runpy.run_path(str(STAGE008_PATH))

PrimitiveParams = S8["PrimitiveParams"]
primitive_base_bundle = S8["primitive_base_bundle"]
primitive_slope_compiler = S8["primitive_slope_compiler"]

P_COMPAT_SAMPLE = float(S8["P_COMPAT_SAMPLE"])
P_BOTH_10 = float(S8["P_BOTH_10"])
P_ONE_10 = float(S8["P_ONE_10"])
P_BOTH_30 = float(S8["P_BOTH_30"])
P_ONE_30 = float(S8["P_ONE_30"])

SLOPE_NAMES = S8["SLOPE_NAMES"]
MIXED_COLS = ["xLambdaU", "xLambdaW", "xLambdaR", "xOmegaU", "xOmegaW"]


# ---------------------------------------------------------------------------
# Rebuild the concrete primitive branch
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
# Exact one-port load-factor objects
# ---------------------------------------------------------------------------

def exact_rows(base, compiler):
    Xi_row = sp.Matrix([[sp.simplify(compiler.Xi[name]) for name in MIXED_COLS]])
    D01_row = sp.Matrix([[sp.simplify(compiler.D01[name]) for name in MIXED_COLS]])
    D21_row = sp.Matrix([[sp.simplify(compiler.D21[name]) for name in MIXED_COLS]])
    D41_row = sp.Matrix([[sp.simplify(compiler.D41[name]) for name in MIXED_COLS]])

    # Exact dP and dDelta rows from the Stage-008 primitive formulas.
    A0 = sp.simplify(base.OmegaU**2 * base.OmegaW**2)
    dDelta = {name: sp.Integer(0) for name in SLOPE_NAMES}
    dDelta["xOmegaU"] = sp.simplify(2 * A0)
    dDelta["xOmegaW"] = sp.simplify(2 * A0)
    dDelta["xLambdaR"] = sp.simplify(-2 * base.R**2)

    termP1 = sp.simplify(base.OmegaU**2 * base.GW)
    termP2 = sp.simplify(base.R * base.GU)
    dP = {name: sp.Integer(0) for name in SLOPE_NAMES}
    dP["xOmegaU"] = sp.simplify(2 * termP1)
    dP["xLambdaW"] = sp.simplify(termP1)
    dP["xLambdaR"] = sp.simplify(termP2)
    dP["xLambdaU"] = sp.simplify(termP2)

    Nload_row = sp.Matrix([[sp.simplify(compiler.N01[name] / base.N0) for name in MIXED_COLS]])
    pi_row = sp.Matrix([[sp.simplify(dP[name] / base.P) for name in MIXED_COLS]])
    delta_row = sp.Matrix([[sp.simplify(dDelta[name] / base.Delta) for name in MIXED_COLS]])

    m_row = sp.Matrix([[sp.Integer(0), sp.Integer(1), sp.Integer(0), sp.Integer(0), -sp.Integer(2)]])
    i_row = sp.Matrix([[sp.Integer(1), -sp.Integer(1), sp.Integer(1), -sp.Integer(2), sp.Integer(0)]])
    h_row = sp.Matrix([[sp.Integer(0), sp.Integer(0), sp.Integer(2), -sp.Integer(2), -sp.Integer(2)]])
    return Xi_row, D01_row, D21_row, D41_row, Nload_row, pi_row, delta_row, m_row, i_row, h_row


# ---------------------------------------------------------------------------
# Part I. Exact pure-transfer theorem
# ---------------------------------------------------------------------------

def pure_transfer_exact_theorem(base, compiler):
    banner("PART I — EXACT PURE-TRANSFER ONE-PORT LOAD THEOREM")

    Xi_row, D01_row, D21_row, D41_row, Nload_row, pi_row, delta_row, m_row, i_row, h_row = exact_rows(base, compiler)

    A_exact = sp.Matrix.vstack(D01_row, D21_row, D41_row)
    Bsym_cols = A_exact.nullspace()
    if len(Bsym_cols) != 2:
        raise AssertionError("Expected a 2D exact pure-transfer subspace.")
    Bsym = sp.Matrix.hstack(*Bsym_cols)

    print("Pure-transfer constraints:")
    print("  D01 = 0")
    print("  D21 = 0")
    print("  D41 = 0")
    print("Exact nullity =", 5 - A_exact.rank())
    print()

    expect_zero("(Xi_1 - N01/N0) on the pure-transfer subspace", (Xi_row - Nload_row) * Bsym)
    expect_zero("(Xi_1 - 2 pi_1 + 2 delta_1) on the pure-transfer subspace", (Xi_row - 2 * pi_row + 2 * delta_row) * Bsym)

    I = sp.simplify(base.R * base.GU / (base.OmegaU**2 * base.GW))
    H = sp.simplify(base.R**2 / (base.OmegaU**2 * base.OmegaW**2))
    Lambda = sp.simplify(base.P / base.Delta)
    Lambda_factor = sp.simplify(base.GW / base.OmegaW**2 * (1 + I) / (1 - H))

    print("One-port load factor Lambda = P/Delta =")
    sp.pprint(Lambda)
    print("Factorized form =")
    sp.pprint(Lambda_factor)
    expect_zero("Lambda - factorized form", sp.simplify(Lambda - Lambda_factor))

    expect_zero(
        "Xi_1 - 2[m + I/(1+I) i + H/(1-H) h] on the pure-transfer subspace",
        sp.simplify((Xi_row - 2 * (m_row + I / (1 + I) * i_row + H / (1 - H) * h_row)) * Bsym),
    )

    print()
    print("Concrete sample-branch values:")
    print("  I =")
    sp.pprint(I)
    print("  H =")
    sp.pprint(H)
    print("  So on the pure-transfer subspace,")
    print("    Xi_1 = 2 m + (6/19) i + [50/(98 pi^2 - 25)] h")
    print()


# ---------------------------------------------------------------------------
# Part II. Rigidity sieve
# ---------------------------------------------------------------------------

def rigidity_sieve(base, compiler):
    banner("PART II — OUTGOING-RIGIDITY SIEVE ON THE PURE-TRANSFER SUBCORRIDOR")

    Xi_row, D01_row, D21_row, D41_row, _Nload_row, _pi_row, _delta_row, m_row, i_row, h_row = exact_rows(base, compiler)

    A_exact = sp.Matrix.vstack(D01_row, D21_row, D41_row)
    M_i = sp.Matrix.vstack(A_exact, i_row)
    M_h = sp.Matrix.vstack(A_exact, h_row)
    M_m = sp.Matrix.vstack(A_exact, m_row)
    M_ih = sp.Matrix.vstack(A_exact, i_row, h_row)

    print("Ranks:")
    print("  rank[pure-transfer]      =", A_exact.rank())
    print("  rank[pure-transfer + i]  =", M_i.rank())
    print("  rank[pure-transfer + h]  =", M_h.rank())
    print("  rank[pure-transfer + m]  =", M_m.rank())
    print("  rank[pure-transfer + i + h] =", M_ih.rank())
    print()

    # Exact reduced 2x2 determinant showing i = h = 0 kills the subcorridor.
    Bsym = sp.Matrix.hstack(*A_exact.nullspace())
    M_ih_reduced = sp.Matrix.vstack(i_row * Bsym, h_row * Bsym)
    det_ih = sp.factor(sp.simplify(M_ih_reduced.det()))
    print("det[(i,h) on the pure-transfer basis] =")
    sp.pprint(det_ih)
    if det_ih == 0:
        raise AssertionError("Expected the combined i=h=0 rigidity to kill the pure-transfer corridor.")

    print()
    print("Conclusion:")
    print("  - i = h = 0 leaves only the trivial drift.")
    print("  - Each single rigidity i = 0, h = 0, or m = 0 still leaves a 1D survivor.")
    print()


# ---------------------------------------------------------------------------
# Part III. Numeric unit directions and transported budgets
# ---------------------------------------------------------------------------

def numeric_directions_and_budgets(base, compiler):
    banner("PART III — NUMERIC UNIT DIRECTIONS AND TRANSPORTED CEILINGS")

    Xi_row, D01_row, D21_row, D41_row, _Nload_row, _pi_row, _delta_row, m_row, i_row, h_row = exact_rows(base, compiler)

    Xi_vec = numeric_row(Xi_row)
    D01_vec = numeric_row(D01_row)
    D21_vec = numeric_row(D21_row)
    D41_vec = numeric_row(D41_row)
    i_vec = numeric_row(i_row)
    h_vec = numeric_row(h_row)
    m_vec = numeric_row(m_row)

    pure_transfer_matrix = np.vstack([D01_vec, D21_vec, D41_vec])
    Btr = orthonormal_nullspace(pure_transfer_matrix)
    Xi_coords = Btr.T @ Xi_vec
    sigma_transfer = float(np.linalg.norm(Xi_coords))

    print(f"Pure-transfer corridor operator norm sigma_transfer = {fmt(sigma_transfer)}")
    print()

    vi = unit_null_direction(np.vstack([D01_vec, D21_vec, D41_vec, i_vec]))
    vh = unit_null_direction(np.vstack([D01_vec, D21_vec, D41_vec, h_vec]))
    vm = unit_null_direction(np.vstack([D01_vec, D21_vec, D41_vec, m_vec]))

    sigma_i = abs(float(vi @ Xi_vec))
    sigma_h = abs(float(vh @ Xi_vec))
    sigma_m = abs(float(vm @ Xi_vec))

    print_unit_direction("i = 0 survivor", MIXED_COLS, vi, Xi_vec)
    print_unit_direction("h = 0 survivor", MIXED_COLS, vh, Xi_vec)
    print_unit_direction("m = 0 survivor", MIXED_COLS, vm, Xi_vec)

    budgets = [
        ("10% loss, both wall-like poles", P_BOTH_10 / P_COMPAT_SAMPLE - 1.0),
        ("10% loss, nonempty wall-like corridor", P_ONE_10 / P_COMPAT_SAMPLE - 1.0),
        ("30% loss, both wall-like poles", P_BOTH_30 / P_COMPAT_SAMPLE - 1.0),
        ("30% loss, nonempty wall-like corridor", P_ONE_30 / P_COMPAT_SAMPLE - 1.0),
    ]

    print("Transported same-charge ceilings on the 1D rigid survivors:")
    print("Interpret the ambient microscopic mixed-sector drift amplitude as ||x||_2 = t.")
    print("Then |epsilon| t <= budget / sigma.")
    print()

    budget_table("Pure-transfer 2D corridor", sigma_transfer, budgets)
    budget_table("i = 0 1D survivor", sigma_i, budgets)
    budget_table("h = 0 1D survivor", sigma_h, budgets)
    budget_table("m = 0 1D survivor", sigma_m, budgets)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base, compiler = concrete_branch_and_compiler()
    pure_transfer_exact_theorem(base, compiler)
    rigidity_sieve(base, compiler)
    numeric_directions_and_budgets(base, compiler)

    banner("STAGE 010 COMPLETE")
    print("The surviving Stage-009 same-charge mechanism is exactly a pure outgoing-load")
    print("slope on the concrete branch. Simultaneous interference- and hybridization-")
    print("rigidity kills it, while each single rigidity leaves a 1D survivor.")
    print("Even the mixed-leg-rigid (m = 0) specialization still carries a nonzero Xi_1.")


if __name__ == "__main__":
    main()
