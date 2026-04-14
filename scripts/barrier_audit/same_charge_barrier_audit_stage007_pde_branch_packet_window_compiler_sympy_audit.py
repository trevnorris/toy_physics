#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage007_pde_branch_packet_window_compiler_sympy_audit.py

Stage 007 — exact compiler from the actual 5PN/moving-throat branch packet to
transported same-charge survival-window tests.

What this script does
---------------------
1. Rebuilds the exact Packet-A / Delta_branch compiler for the outgoing static
   prefactor lanes P0^(20), P0^(21), P0^(22).
2. Converts the Stage-174 normalization defect Delta_norm into the actual
   isotropic prefactor \bar P0 seen by the same-charge window test.
3. Compiles the Stage-006 primitive-family survival ceilings into exact
   inequalities on the real branch packet.
4. Specializes to the weak-axisymmetric grouped signature
      lambda_(20)=1, lambda_(21)=1/2, lambda_(22)=-1,
   verifies the exact anomaly relation b_{P0}=3 a_{P0}, and rewrites the
   survival test as a single inequality in the scalar amplitude Xi_1=P1/P0.
5. Evaluates the explicit Stage-006 compatibility point and records the
   remaining weak-axisymmetric headroom before the transported primitive-family
   window closes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

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


def grouped_trace_anomaly(x20: sp.Expr, x21: sp.Expr, x22: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    return xbar, ax, bx


def grouped_inverse(xbar: sp.Expr, ax: sp.Expr, bx: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    x20 = sp.simplify(xbar + 4 * ax)
    x21 = sp.simplify(xbar - ax + bx)
    x22 = sp.simplify(xbar - ax - bx)
    return x20, x21, x22


# ---------------------------------------------------------------------------
# Carried primitive-family ceilings from Stage 006
# ---------------------------------------------------------------------------

P_BOTH_10 = sp.Float("0.0028313316855593175")
P_ONE_10 = sp.Float("0.0035965105896846573")
P_BOTH_30 = sp.Float("0.00817339430971383")
P_ONE_30 = sp.Float("0.0116633929790174")

P_COMPAT_SAMPLE = sp.Float("0.002069792318062885")


@dataclass
class CeilingBudget:
    name: str
    Pcrit: float
    xi_budget: float
    aP_budget: float


# ---------------------------------------------------------------------------
# Part I. Exact Packet-A / Delta_branch compiler
# ---------------------------------------------------------------------------

def verify_exact_branch_packet_compiler() -> None:
    banner("PART I — EXACT BRANCH-PACKET COMPILER TO LANE PREFATORS")

    G, cs, a, c, mhat0 = sp.symbols("G c_s a c mhat_0", positive=True, real=True)
    Delta_norm = sp.symbols("Delta_norm", real=True)
    aP0, bP0 = sp.symbols("a_P0 b_P0", real=True)

    Tquad = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5))
    Pbar0 = sp.simplify((Delta_norm + Tquad) / mhat0**2)
    P20, P21, P22 = grouped_inverse(Pbar0, aP0, bP0)

    print("T_quad =")
    sp.pprint(Tquad)
    print("Pbar0 from Delta_norm =")
    sp.pprint(Pbar0)
    print("Lane prefactors P0^(20), P0^(21), P0^(22) =")
    sp.pprint(sp.Matrix([P20, P21, P22]))

    Pbar_chk, a_chk, b_chk = grouped_trace_anomaly(P20, P21, P22)
    expect_zero("grouped inverse/trace roundtrip: Pbar", Pbar_chk - Pbar0)
    expect_zero("grouped inverse/trace roundtrip: a_P0", a_chk - aP0)
    expect_zero("grouped inverse/trace roundtrip: b_P0", b_chk - bP0)


# ---------------------------------------------------------------------------
# Part II. Exact isotropic window compiler
# ---------------------------------------------------------------------------

def verify_isotropic_window_compiler() -> None:
    banner("PART II — EXACT ISOTROPIC WINDOW TESTS")

    G, cs, a, c, mhat0 = sp.symbols("G c_s a c mhat_0", positive=True, real=True)
    Delta_norm = sp.symbols("Delta_norm", real=True)
    Pcrit = sp.symbols("P_crit", positive=True, real=True)

    Tquad = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5))
    Pbar0 = sp.simplify((Delta_norm + Tquad) / mhat0**2)
    Delta_norm_max = sp.simplify(mhat0**2 * Pcrit - Tquad)
    mhat_bound = sp.simplify(Tquad / Pcrit)

    print("Isotropic prefactor =")
    sp.pprint(Pbar0)
    print("Maximum admissible Delta_norm at ceiling P_crit =")
    sp.pprint(Delta_norm_max)
    print("Calibrated-branch lower bound on mhat_0^2 when Delta_norm = 0 =")
    sp.pprint(mhat_bound)

    expect_zero(
        "ceiling condition equivalence",
        sp.simplify(Pbar0 - Pcrit).subs({Delta_norm: Delta_norm_max}),
    )
    expect_zero(
        "calibrated branch bound equivalence",
        sp.simplify(Pbar0.subs({Delta_norm: 0, mhat0**2: mhat_bound}) - Pcrit),
    )


# ---------------------------------------------------------------------------
# Part III. Weak-axisymmetric grouped signature and Xi_1 compiler
# ---------------------------------------------------------------------------

def verify_axisymmetric_signature_compiler() -> None:
    banner("PART III — WEAK-AXISYMMETRIC GROUPED SIGNATURE AND Xi_1 COMPILER")

    P0bar = sp.symbols("Pbar_0", positive=True, real=True)
    eps, Xi1 = sp.symbols("epsilon Xi_1", real=True)
    z = sp.symbols("z", real=True)

    P20 = sp.simplify(P0bar * (1 + eps * Xi1))
    P21 = sp.simplify(P0bar * (1 + sp.Rational(1, 2) * eps * Xi1))
    P22 = sp.simplify(P0bar * (1 - eps * Xi1))

    Pbar_chk, aP_chk, bP_chk = grouped_trace_anomaly(P20, P21, P22)

    print("Axisymmetric weak-anisotropy lane prefactors =")
    sp.pprint(sp.Matrix([P20, P21, P22]))
    print("Grouped trace/anomaly variables =")
    sp.pprint(sp.Matrix([Pbar_chk, aP_chk, bP_chk]))

    expect_zero("axisymmetric grouped trace", Pbar_chk - P0bar)
    expect_zero("axisymmetric a_P0", aP_chk - sp.simplify(eps * P0bar * Xi1 / 4))
    expect_zero("axisymmetric b_P0", bP_chk - sp.simplify(3 * eps * P0bar * Xi1 / 4))
    expect_zero("axisymmetric exact law b_P0 - 3 a_P0", bP_chk - 3 * aP_chk)

    # Rewrite with z = epsilon * Xi_1
    P20z = sp.simplify(P20.subs(eps * Xi1, z))
    P21z = sp.simplify(P21.subs(eps * Xi1, z))
    P22z = sp.simplify(P22.subs(eps * Xi1, z))

    print("With z := epsilon * Xi_1, lane prefactors become =")
    sp.pprint(sp.Matrix([P20z, P21z, P22z]))

    # Positive-z branch: the worst lane is 20.
    zpos = sp.symbols("zpos", nonnegative=True, real=True)
    expect_zero(
        "positive-z max lane minus Pbar0(1+z)",
        sp.simplify(P20z.subs(z, zpos) - P0bar * (1 + zpos)),
    )
    expect_zero(
        "positive-z min lane minus Pbar0(1-z)",
        sp.simplify(P22z.subs(z, zpos) - P0bar * (1 - zpos)),
    )

    # Negative-z branch: the worst lane is 22.
    upos = sp.symbols("upos", nonnegative=True, real=True)
    expect_zero(
        "negative-z max lane minus Pbar0(1+u)",
        sp.simplify(P22z.subs(z, -upos) - P0bar * (1 + upos)),
    )
    expect_zero(
        "negative-z min lane minus Pbar0(1-u)",
        sp.simplify(P20z.subs(z, -upos) - P0bar * (1 - upos)),
    )

    print("Robust all-lane ceiling condition therefore reduces to:")
    sp.pprint(sp.simplify(P0bar * (1 + sp.Abs(z))))
    print("Selective lowered-lane condition reduces to:")
    sp.pprint(sp.simplify(P0bar * (1 - sp.Abs(z))))


# ---------------------------------------------------------------------------
# Part IV. Concrete Stage-006 compatibility-point headroom
# ---------------------------------------------------------------------------

def compute_headroom_budgets() -> list[CeilingBudget]:
    out: list[CeilingBudget] = []
    for name, crit in [
        ("10% loss, both wall-like poles", float(P_BOTH_10)),
        ("10% loss, nonempty wall-like corridor", float(P_ONE_10)),
        ("30% loss, both wall-like poles", float(P_BOTH_30)),
        ("30% loss, nonempty wall-like corridor", float(P_ONE_30)),
    ]:
        xi_budget = crit / float(P_COMPAT_SAMPLE) - 1.0
        aP_budget = (crit - float(P_COMPAT_SAMPLE)) / 4.0
        out.append(CeilingBudget(name=name, Pcrit=crit, xi_budget=xi_budget, aP_budget=aP_budget))
    return out


def print_headroom_budgets() -> None:
    banner("PART IV — EXPLICIT HEADROOM AT THE STAGE-006 COMPATIBILITY POINT")
    print(f"P0_target,compat (Stage-006 sample) = {float(P_COMPAT_SAMPLE):.18f}")
    print()
    print("Transported robust weak-axisymmetric budgets:")
    print("  |epsilon * Xi_1| <= P_crit / Pbar0 - 1")
    print("  |a_P0| <= (P_crit - Pbar0) / 4")

    budgets = compute_headroom_budgets()
    for row in budgets:
        print()
        print(f"[{row.name}]")
        print(f"  P_crit     = {row.Pcrit:.18f}")
        print(f"  Xi budget  = {row.xi_budget:.18f}")
        print(f"  a_P0 budget= {row.aP_budget:.18f}")

    # Sanity check: the ordering of the budgets follows the ordering of the ceilings.
    xi_list = [row.xi_budget for row in budgets]
    if not (xi_list[0] < xi_list[1] < xi_list[2] < xi_list[3]):
        raise AssertionError("Headroom budgets are not ordered as expected.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    verify_exact_branch_packet_compiler()
    verify_isotropic_window_compiler()
    verify_axisymmetric_signature_compiler()
    print_headroom_budgets()

    banner("STAGE 007 LEDGER")
    print("1. The actual 5PN branch packet now compiles directly to the lane prefactors")
    print("      P0^(20), P0^(21), P0^(22).")
    print("2. The isotropic normalization defect Delta_norm already determines the mean")
    print("      prefactor Pbar0 seen by the same-charge window test.")
    print("3. On the weak-axisymmetric grouped signature, the transported ceiling test")
    print("      collapses to one scalar amplitude Xi_1 = P1/P0.")
    print("4. The Stage-006 compatibility point still has finite anisotropy headroom, but")
    print("      that headroom is explicit and not large on the stricter 10% benchmark.")


if __name__ == "__main__":
    main()
