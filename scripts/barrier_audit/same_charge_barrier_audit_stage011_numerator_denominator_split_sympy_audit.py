#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage011_numerator_denominator_split_sympy_audit.py

Stage 011 — split the surviving Stage-010 pure-transfer load-factor corridor into
numerator-rigid and denominator-rigid subcorridors, then carry both into the
actual wall-like dynamic window.

What this script does
---------------------
1. Rebuilds the concrete Stage-006/008 primitive compatibility branch and the
   Stage-010 exact pure-transfer corridor
       D01 = D21 = D41 = 0.
2. Defines the exact static numerator and denominator load slopes
       pi_1    := P01 / P,
       delta_1 := Delta01 / Delta,
   so that on the pure-transfer corridor
       Xi_1 = 2(pi_1 - delta_1).
3. Proves that:
      - pure-transfer + pi_1 = 0 leaves a 1D numerator-rigid survivor,
      - pure-transfer + delta_1 = 0 leaves a 1D denominator-rigid survivor,
      - imposing both pi_1 = 0 and delta_1 = 0 kills the corridor.
4. Extracts positively-oriented unit directions on the two 1D survivors and
   reports their exact static Xi_1 amplitudes.
5. Carries those directions into the full pole census of the concrete branch
   and measures, by symmetric finite difference, the first-order log-slopes of
      - P0,
      - the lower wall-like residue/linewidth figure R_Q,-,
      - the upper wall-like residue/linewidth figure R_Q,+.
6. Compares the resulting dynamic ceilings against the already carried Stage-007
   static ceilings.

Main structural result
----------------------
The pure-transfer corridor survives both numerator-rigid and denominator-rigid
splits.  Simultaneous numerator and denominator rigidity kills it exactly.
On the concrete compatibility branch, however, neither surviving split is killed
by the actual wall-like dynamic window:

- numerator-rigid positive-Xi motion hurts the upper wall pole but *helps* the
  lower wall pole, so the nonempty dynamic corridor remains open at first order;
- denominator-rigid positive-Xi motion weakens both wall poles, but its dynamic
  ceiling still lies above the transported static ceiling.

So after Stage 011 the live question is no longer whether numerator- or
 denominator-rigidity kills the mechanism outright.  It is which split the real
PDE-selected mixed branch actually resembles.
"""

from __future__ import annotations

import contextlib
import io
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
    if math.isinf(x):
        return "inf"
    return f"{x:.15f}"


def numeric_row(row: sp.Matrix) -> np.ndarray:
    return np.array([float(sp.N(z, 80)) for z in list(row)], dtype=float)


def orthonormal_nullspace(A: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    _u, s, vh = np.linalg.svd(A)
    rank = int((s > tol).sum())
    B = vh[rank:].T.copy()
    if B.size == 0:
        return B.reshape(A.shape[1], 0)
    B[np.abs(B) < 1e-14] = 0.0
    return B


def unit_null_direction(A: np.ndarray, Xi_vec: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    B = orthonormal_nullspace(A, tol=tol)
    if B.shape[1] != 1:
        raise AssertionError("Expected a 1D nullspace.")
    v = B[:, 0].copy()
    if float(v @ Xi_vec) < 0:
        v = -v
    v[np.abs(v) < 1e-14] = 0.0
    return v


def print_unit_direction(name: str, cols: list[str], v: np.ndarray, Xi_vec: np.ndarray, pi_vec: np.ndarray, delta_vec: np.ndarray) -> None:
    print(f"{name} unit direction (oriented to Xi_1 > 0) =")
    for c, x in zip(cols, v):
        print(f"  {c:9s} -> {fmt(float(x))}")
    print(f"  Xi_1 amplitude     = {fmt(float(v @ Xi_vec))}")
    print(f"  pi_1 amplitude     = {fmt(float(v @ pi_vec))}")
    print(f"  delta_1 amplitude  = {fmt(float(v @ delta_vec))}")
    print()


# ---------------------------------------------------------------------------
# Load prior stages quietly
# ---------------------------------------------------------------------------

BASE = Path("/mnt/data")
with contextlib.redirect_stdout(io.StringIO()):
    S5 = runpy.run_path(str(BASE / "same_charge_barrier_audit_stage005_concrete_branch_residue_linewidth_test_sympy_audit.py"))
    S8 = runpy.run_path(str(BASE / "same_charge_barrier_audit_stage008_microscopic_xi1_primitive_compiler_sympy_audit.py"))
    S10 = runpy.run_path(str(BASE / "same_charge_barrier_audit_stage010_pure_transfer_load_factor_rigidity_sympy_audit.py"))

PrimitiveParams = S8["PrimitiveParams"]
primitive_base_bundle = S8["primitive_base_bundle"]
primitive_slope_compiler = S8["primitive_slope_compiler"]

BranchParams = S5["BranchParams"]
branch_couplings = S5["branch_couplings"]
pole_census = S5["pole_census"]

exact_rows = S10["exact_rows"]

P_COMPAT_SAMPLE = float(S8["P_COMPAT_SAMPLE"])
P_BOTH_10 = float(S8["P_BOTH_10"])
P_ONE_10 = float(S8["P_ONE_10"])
P_BOTH_30 = float(S8["P_BOTH_30"])
P_ONE_30 = float(S8["P_ONE_30"])

MIXED_COLS = ["xLambdaU", "xLambdaW", "xLambdaR", "xOmegaU", "xOmegaW"]


# ---------------------------------------------------------------------------
# Rebuild the concrete compatibility branch and exact rows
# ---------------------------------------------------------------------------

def concrete_branch_and_rows():
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
    return base, compiler, exact_rows(base, compiler)


# ---------------------------------------------------------------------------
# Part I. Exact numerator / denominator split theorem
# ---------------------------------------------------------------------------

def exact_split_theorem(base, compiler):
    banner("PART I — EXACT NUMERATOR / DENOMINATOR SPLIT THEOREM")

    Xi_row, D01_row, D21_row, D41_row, _Nload_row, pi_row, delta_row, _m_row, _i_row, _h_row = exact_rows(base, compiler)

    print("On the Stage-010 pure-transfer corridor:")
    print("  D01 = D21 = D41 = 0")
    print("  Xi_1 = 2(pi_1 - delta_1)")
    print()

    print("Exact static numerator slope pi_1 = P01/P =")
    sp.pprint(pi_row)
    print("Exact static denominator slope delta_1 = Delta01/Delta =")
    sp.pprint(delta_row)
    print("Exact Xi_1 row =")
    sp.pprint(Xi_row)

    A_pure = sp.Matrix.vstack(D01_row, D21_row, D41_row)
    Bsym = sp.Matrix.hstack(*A_pure.nullspace())
    expect_zero("(Xi_1 - 2(pi_1 - delta_1)) on the pure-transfer basis", sp.simplify((Xi_row - 2 * (pi_row - delta_row)) * Bsym))

    M_num = sp.Matrix.vstack(A_pure, pi_row)
    M_den = sp.Matrix.vstack(A_pure, delta_row)
    M_both = sp.Matrix.vstack(A_pure, pi_row, delta_row)

    print()
    print("Exact ranks / nullities in the 5D mixed primitive slope space:")
    print("  pure-transfer          : rank =", A_pure.rank(), ", nullity =", 5 - A_pure.rank())
    print("  + numerator rigidity   : rank =", M_num.rank(), ", nullity =", 5 - M_num.rank())
    print("  + denominator rigidity : rank =", M_den.rank(), ", nullity =", 5 - M_den.rank())
    print("  + both rigidities      : rank =", M_both.rank(), ", nullity =", 5 - M_both.rank())

    M_red = sp.Matrix.vstack(pi_row * Bsym, delta_row * Bsym)
    det_red = sp.factor(sp.simplify(M_red.det()))
    print()
    print("det[(pi_1, delta_1) on the pure-transfer basis] =")
    sp.pprint(det_red)
    if det_red == 0:
        raise AssertionError("Expected simultaneous numerator and denominator rigidity to kill the corridor.")

    print()
    print("Conclusion:")
    print("  - pure-transfer + pi_1 = 0 leaves a 1D numerator-rigid survivor;")
    print("  - pure-transfer + delta_1 = 0 leaves a 1D denominator-rigid survivor;")
    print("  - imposing both pi_1 = 0 and delta_1 = 0 kills the corridor exactly.")


# ---------------------------------------------------------------------------
# Part II. Positive-Xi unit survivors
# ---------------------------------------------------------------------------

def positive_xi_unit_survivors(base, compiler):
    banner("PART II — POSITIVE-XI UNIT SURVIVORS")

    Xi_row, D01_row, D21_row, D41_row, _Nload_row, pi_row, delta_row, _m_row, _i_row, _h_row = exact_rows(base, compiler)
    Xi_vec = numeric_row(Xi_row)
    D01_vec = numeric_row(D01_row)
    D21_vec = numeric_row(D21_row)
    D41_vec = numeric_row(D41_row)
    pi_vec = numeric_row(pi_row)
    delta_vec = numeric_row(delta_row)

    A_pure = np.vstack([D01_vec, D21_vec, D41_vec])
    v_num = unit_null_direction(np.vstack([A_pure, pi_vec]), Xi_vec)
    v_den = unit_null_direction(np.vstack([A_pure, delta_vec]), Xi_vec)

    print_unit_direction("Numerator-rigid (pi_1 = 0)", MIXED_COLS, v_num, Xi_vec, pi_vec, delta_vec)
    print_unit_direction("Denominator-rigid (delta_1 = 0)", MIXED_COLS, v_den, Xi_vec, pi_vec, delta_vec)

    sigma_num = abs(float(v_num @ Xi_vec))
    sigma_den = abs(float(v_den @ Xi_vec))

    print("Static Xi_1 amplitudes:")
    print(f"  numerator-rigid sigma_num = {fmt(sigma_num)}")
    print(f"  denominator-rigid sigma_den = {fmt(sigma_den)}")
    print()
    print("Interpretation:")
    print("  - The numerator-rigid branch generates a larger Xi_1 per unit mixed drift.")
    print("  - The denominator-rigid branch is the gentler static lever.")

    return v_num, v_den, sigma_num, sigma_den


# ---------------------------------------------------------------------------
# Part III. Dynamic-window audit on the concrete compatibility branch
# ---------------------------------------------------------------------------

def perturbed_params(base, v: np.ndarray, eps: float) -> BranchParams:
    return BranchParams(
        lamB=float(base.lamB),
        lamU=float(base.lamU) * math.exp(eps * float(v[0])),
        lamW=float(base.lamW) * math.exp(eps * float(v[1])),
        lamR=float(base.lamR) * math.exp(eps * float(v[2])),
        OmegaU=float(base.OmegaU) * math.exp(eps * float(v[3])),
        OmegaW=float(base.OmegaW) * math.exp(eps * float(v[4])),
        varpi=float(base.varpi),
        K=float(base.Kcompat),
        M=float(base.M),
        a=1.0,
        cs=1.0,
    )


def wall_like_summary(p: BranchParams):
    couplings = branch_couplings(p)
    poles = pole_census(p)
    wall = sorted([pd for pd in poles if pd.family == "wall-like"], key=lambda pd: pd.omega)
    if len(wall) != 2:
        raise AssertionError("Expected exactly two wall-like poles on the concrete branch.")
    lower, upper = wall[0], wall[1]
    return {
        "P0": couplings["P0"],
        "lower_RQ": lower.R_Q,
        "upper_RQ": upper.R_Q,
        "lower_omega": lower.omega,
        "upper_omega": upper.omega,
        "lower_N": lower.N_star,
        "upper_N": upper.N_star,
    }


def symmetric_log_slope(base, v: np.ndarray, eps: float = 1.0e-6):
    plus = wall_like_summary(perturbed_params(base, v, eps))
    minus = wall_like_summary(perturbed_params(base, v, -eps))

    def dlog(key: str) -> float:
        return (math.log(plus[key]) - math.log(minus[key])) / (2 * eps)

    return {
        "dlogP0": dlog("P0"),
        "lower_dlogR": dlog("lower_RQ"),
        "upper_dlogR": dlog("upper_RQ"),
        "lower_dlogOmega": dlog("lower_omega"),
        "upper_dlogOmega": dlog("upper_omega"),
        "lower_dlogN": dlog("lower_N"),
        "upper_dlogN": dlog("upper_N"),
    }


def ratio_threshold(eta: float) -> float:
    V_known_at_1 = 1.181909222592
    epsilon = 0.1
    DeltaV_req = V_known_at_1 - epsilon
    return 2.0 * DeltaV_req * (1.0 + eta**2) / eta


def dynamic_ceiling(R0: float, dlogR: float, Rreq: float) -> float:
    if dlogR >= 0:
        return math.inf
    return math.log(R0 / Rreq) / (-dlogR)


def dynamic_window_audit(base, v_num: np.ndarray, v_den: np.ndarray, sigma_num: float, sigma_den: float) -> None:
    banner("PART III — FIRST DYNAMIC-WINDOW AUDIT OF THE TWO RIGID SPLITS")

    base_summary = wall_like_summary(perturbed_params(base, np.zeros(5), 0.0))
    print("Concrete compatibility-point wall-like poles:")
    print(f"  lower wall: omega = {fmt(base_summary['lower_omega'])}, R_Q = {fmt(base_summary['lower_RQ'])}")
    print(f"  upper wall: omega = {fmt(base_summary['upper_omega'])}, R_Q = {fmt(base_summary['upper_RQ'])}")
    print(f"  static P0  = {fmt(base_summary['P0'])}")
    print()

    num_slopes = symmetric_log_slope(base, v_num)
    den_slopes = symmetric_log_slope(base, v_den)

    print("Positive-Xi log-slopes from symmetric finite difference:")
    print("Numerator-rigid branch:")
    for key in ("dlogP0", "upper_dlogR", "lower_dlogR", "upper_dlogOmega", "lower_dlogOmega"):
        print(f"  {key:16s} = {fmt(num_slopes[key])}")
    print("Denominator-rigid branch:")
    for key in ("dlogP0", "upper_dlogR", "lower_dlogR", "upper_dlogOmega", "lower_dlogOmega"):
        print(f"  {key:16s} = {fmt(den_slopes[key])}")
    print()

    # Static transported ceilings from Stage 007 / Stage 010.
    stat_both_10_num = (P_BOTH_10 / P_COMPAT_SAMPLE - 1.0) / sigma_num
    stat_one_10_num = (P_ONE_10 / P_COMPAT_SAMPLE - 1.0) / sigma_num
    stat_both_30_num = (P_BOTH_30 / P_COMPAT_SAMPLE - 1.0) / sigma_num
    stat_one_30_num = (P_ONE_30 / P_COMPAT_SAMPLE - 1.0) / sigma_num

    stat_both_10_den = (P_BOTH_10 / P_COMPAT_SAMPLE - 1.0) / sigma_den
    stat_one_10_den = (P_ONE_10 / P_COMPAT_SAMPLE - 1.0) / sigma_den
    stat_both_30_den = (P_BOTH_30 / P_COMPAT_SAMPLE - 1.0) / sigma_den
    stat_one_30_den = (P_ONE_30 / P_COMPAT_SAMPLE - 1.0) / sigma_den

    # Dynamic ceilings from the actual wall-like poles on the concrete branch.
    for eta, tag in ((0.1, "10% loss"), (0.3, "30% loss")):
        Rreq = ratio_threshold(eta)
        print(f"Dynamic threshold [{tag}]  R_Q,req = {fmt(Rreq)}")

        num_upper = dynamic_ceiling(base_summary["upper_RQ"], num_slopes["upper_dlogR"], Rreq)
        num_lower = dynamic_ceiling(base_summary["lower_RQ"], num_slopes["lower_dlogR"], Rreq)
        num_both = min(num_upper, num_lower)
        num_one = max(num_upper, num_lower)

        den_upper = dynamic_ceiling(base_summary["upper_RQ"], den_slopes["upper_dlogR"], Rreq)
        den_lower = dynamic_ceiling(base_summary["lower_RQ"], den_slopes["lower_dlogR"], Rreq)
        den_both = min(den_upper, den_lower)
        den_one = max(den_upper, den_lower)

        print("  Numerator-rigid dynamic ceilings:")
        print(f"    both wall-like poles : |epsilon| t <= {fmt(num_both)}")
        print(f"    nonempty wall corridor: |epsilon| t <= {fmt(num_one)}")
        print("  Denominator-rigid dynamic ceilings:")
        print(f"    both wall-like poles : |epsilon| t <= {fmt(den_both)}")
        print(f"    nonempty wall corridor: |epsilon| t <= {fmt(den_one)}")

        if eta == 0.1:
            print("  Transported static ceilings on the same branch:")
            print(f"    numerator-rigid  both/nonempty = {fmt(stat_both_10_num)} / {fmt(stat_one_10_num)}")
            print(f"    denominator-rigid both/nonempty = {fmt(stat_both_10_den)} / {fmt(stat_one_10_den)}")
            if not (num_both > stat_both_10_num and den_both > stat_both_10_den):
                raise AssertionError("Expected the 10% dynamic both-pole ceilings to remain above the 10% static ceilings.")
        if eta == 0.3:
            print("  Transported static ceilings on the same branch:")
            print(f"    numerator-rigid  both/nonempty = {fmt(stat_both_30_num)} / {fmt(stat_one_30_num)}")
            print(f"    denominator-rigid both/nonempty = {fmt(stat_both_30_den)} / {fmt(stat_one_30_den)}")
            if not (num_both > stat_both_30_num and den_both > stat_both_30_den):
                raise AssertionError("Expected the 30% dynamic both-pole ceilings to remain above the 30% static ceilings.")
        print()

    print("Interpretation:")
    print("  - Numerator-rigid positive-Xi motion decreases the upper wall-like R_Q but increases the lower wall-like R_Q.")
    print("    So at first order it has no finite nonempty dynamic kill on the concrete branch.")
    print("  - Denominator-rigid positive-Xi motion decreases both wall-like R_Q values,")
    print("    so it is the only split with a genuinely finite nonempty dynamic ceiling.")
    print("  - But on the actual sample compatibility point even that dynamic ceiling still sits above")
    print("    the transported static ceiling. So neither branch dies dynamically before the static budget is hit.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    base, compiler, _rows = concrete_branch_and_rows()
    exact_split_theorem(base, compiler)
    v_num, v_den, sigma_num, sigma_den = positive_xi_unit_survivors(base, compiler)
    dynamic_window_audit(base, v_num, v_den, sigma_num, sigma_den)

    banner("STAGE 011 COMPLETE")
    print("The surviving Stage-010 pure-transfer corridor splits cleanly into")
    print("numerator-rigid and denominator-rigid 1D subcorridors.  Imposing both")
    print("rigidities kills the corridor exactly.  On the concrete compatibility point,")
    print("however, neither surviving split is killed by the wall-like dynamic window;")
    print("the actual first ceiling remains the transported static window.")


if __name__ == "__main__":
    main()
