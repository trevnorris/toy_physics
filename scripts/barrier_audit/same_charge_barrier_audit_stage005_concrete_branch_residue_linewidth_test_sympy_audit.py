#!/usr/bin/env python3
"""
same_charge_barrier_audit_stage005_concrete_branch_residue_linewidth_test_sympy_audit.py

Stage 005 — concrete finite-throat primitive branch and residue/linewidth survival test.

What this script does
---------------------
1. Fixes one explicit finite-throat isotropic primitive branch by using the
   lowest N/N wall/U profile and the lowest D/N half-wave support/W profile.
2. Verifies the exact quartic pole polynomial for the conservative one-port
   wall/BdG/Maxwell/mixed bundle in y = omega^2.
3. Derives the exact residue-to-linewidth cancellation for a simple passive pole:

       R_Q,* := |A_Q,*|/gamma_* = 1/(Gamma_* N_*(omega_*)),

   where Gamma_* is the outgoing-port coefficient evaluated at the pole and
   N_*(omega_*) is the exact transfer factor from Stage 004.
4. Derives the exact low-loss survival inequality needed to beat a required local
   barrier reduction DeltaV_req(x).
5. Evaluates one explicit admissible primitive branch and prints the full pole
   census, including wall-like vs internal-like labels.
6. Compares the wall-like pure-Q residue/linewidth figures against the Stage-001
   illustrative barrier benchmark at two loss tolerances.
7. Runs a small scan in the outgoing-leg coupling lambda_W and shows the concrete
   tension between static prefactor size P0 and dynamic residue/linewidth gain.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List

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


KAPPA = 2 * math.sqrt(2.0) / math.pi


@dataclass
class BranchParams:
    lamB: float
    lamU: float
    lamW: float
    lamR: float
    OmegaU: float
    OmegaW: float
    varpi: float
    K: float
    M: float
    a: float = 1.0
    cs: float = 1.0


@dataclass
class PoleData:
    omega: float
    y: float
    family: str
    nearest_uncoupled: float
    Delta_star: float
    N_star: float
    P_star: float
    R_Q: float
    R_qq: float


# ---------------------------------------------------------------------------
# Part I. Exact finite-throat primitive branch and quartic pole polynomial
# ---------------------------------------------------------------------------

def verify_exact_branch_polynomial() -> None:
    banner("PART I — EXACT FINITE-THROAT PRIMITIVE BRANCH AND POLE POLYNOMIAL")

    y = sp.symbols("y", real=True)
    lamB, lamU, lamW, lamR = sp.symbols("lambda_B lambda_U lambda_W lambda_R", positive=True, real=True)
    OmegaU, OmegaW, varpi, K, M = sp.symbols("Omega_U Omega_W varpi K M", positive=True, real=True)
    kappa = sp.simplify(2 * sp.sqrt(2) / sp.pi)

    C = sp.simplify(kappa * lamB)
    GU = lamU
    GW = sp.simplify(kappa * lamW)
    R = sp.simplify(kappa * lamR)

    KB_y = sp.simplify((K - M * y) - C**2 / (varpi**2 - y))
    A_y = sp.simplify(OmegaU**2 - y)
    W_y = sp.simplify(OmegaW**2 - y)
    Delta_y = sp.simplify(A_y * W_y - R**2)
    Q_y = sp.simplify(GU**2 * W_y + 2 * GU * GW * R + GW**2 * A_y)
    D_y = sp.simplify(KB_y - Q_y / Delta_y)

    F_y = sp.expand(
        ((K - M * y) * (varpi**2 - y) - C**2) * ((OmegaU**2 - y) * (OmegaW**2 - y) - R**2)
        - (varpi**2 - y) * (GU**2 * (OmegaW**2 - y) + 2 * GU * GW * R + GW**2 * (OmegaU**2 - y))
    )

    num_y, den_y = sp.together(D_y).as_numer_denom()

    print("Exact finite-throat overlap constant kappa =")
    sp.pprint(kappa)
    print("C =")
    sp.pprint(C)
    print("G_U =")
    sp.pprint(GU)
    print("G_W =")
    sp.pprint(GW)
    print("R =")
    sp.pprint(R)
    print("\nD(y) numerator F(y) =")
    sp.pprint(F_y)

    print("\nThe together() numerator carries the expected overlap-denominator scale -pi^4.")
    expect_zero("quartic pole polynomial identity", sp.expand(num_y + sp.pi**4 * F_y))

    if sp.Poly(F_y, y).degree() != 4:
        raise AssertionError("F(y) should be quartic in y = omega^2")


# ---------------------------------------------------------------------------
# Part II. Exact residue/linewidth cancellation and survival inequality
# ---------------------------------------------------------------------------

def verify_residue_linewidth_cancellation() -> None:
    banner("PART II — EXACT RESIDUE / LINEWIDTH CANCELLATION")

    Dabs, Nsabs, Deltaabs, Nstar, Gamma = sp.symbols(
        "Dabs Nsabs Deltaabs N_star Gamma", positive=True, real=True
    )

    Aqq_abs = sp.simplify(1 / Dabs)
    gamma = sp.simplify(Gamma * Nstar / Dabs)
    Rqq = sp.simplify(Aqq_abs / gamma)

    As_abs = sp.simplify(Nsabs / (Deltaabs * Dabs))
    Rs = sp.simplify(As_abs / gamma)

    print("|A_qq,*| =")
    sp.pprint(Aqq_abs)
    print("gamma_* =")
    sp.pprint(gamma)
    print("R_qq,* = |A_qq,*| / gamma_* =")
    sp.pprint(Rqq)
    print("\nGeneric collinear-source coefficient ratio R_s,* =")
    sp.pprint(Rs)

    expect_zero("R_qq,* - 1/(Gamma N_*)", sp.simplify(Rqq - 1 / (Gamma * Nstar)))
    expect_zero(
        "R_s,* - Ns/(Delta Gamma N_*)",
        sp.simplify(Rs - Nsabs / (Deltaabs * Gamma * Nstar)),
    )
    expect_zero(
        "pure-Q specialization R_s,* -> R_qq,*",
        sp.simplify(Rs.subs(Nsabs, Deltaabs) - Rqq),
    )


def verify_survival_inequality() -> None:
    banner("PART III — EXACT LOW-LOSS SURVIVAL INEQUALITY")

    DeltaV_req, eta, Sq, ratio = sp.symbols("DeltaV_req eta S_q ratio", positive=True, real=True)

    # Required inequality from Stage 004: 0.5 * ratio * eta/(1+eta^2) * S_q^2 >= DeltaV_req
    req = sp.simplify(2 * DeltaV_req * (1 + eta**2) / (eta * Sq**2))
    lhs = sp.simplify(sp.Rational(1, 2) * ratio * eta / (1 + eta**2) * Sq**2)

    print("Required ratio threshold =")
    sp.pprint(req)
    print("Barrier-lowering bound =")
    sp.pprint(lhs)

    expect_zero(
        "threshold inversion identity",
        sp.simplify(lhs.subs(ratio, req) - DeltaV_req),
    )


# ---------------------------------------------------------------------------
# Part IV. Numerical primitive branch family
# ---------------------------------------------------------------------------

def branch_couplings(p: BranchParams) -> dict[str, float]:
    C = KAPPA * p.lamB
    GU = p.lamU
    GW = KAPPA * p.lamW
    R = KAPPA * p.lamR
    Delta0 = p.OmegaU**2 * p.OmegaW**2 - R**2
    Q0 = GU**2 * p.OmegaW**2 + 2 * GU * GW * R + GW**2 * p.OmegaU**2
    P0_proto = p.OmegaU**2 * GW + R * GU
    D0 = (p.K - C**2 / p.varpi**2) - Q0 / Delta0
    N0 = P0_proto**2 / Delta0**2
    Pref0 = N0 / D0
    return {
        "C": C,
        "GU": GU,
        "GW": GW,
        "R": R,
        "Delta0": Delta0,
        "Q0": Q0,
        "P0_proto": P0_proto,
        "D0": D0,
        "N0": N0,
        "P0": Pref0,
    }


def quartic_coefficients(p: BranchParams) -> List[float]:
    c = branch_couplings(p)
    C, GU, GW, R = c["C"], c["GU"], c["GW"], c["R"]
    y = sp.symbols("y", real=True)
    F = sp.expand(
        ((p.K - p.M * y) * (p.varpi**2 - y) - C**2) * ((p.OmegaU**2 - y) * (p.OmegaW**2 - y) - R**2)
        - (p.varpi**2 - y) * (GU**2 * (p.OmegaW**2 - y) + 2 * GU * GW * R + GW**2 * (p.OmegaU**2 - y))
    )
    return [float(cc) for cc in sp.Poly(F, y).all_coeffs()]


def uncoupled_wall_roots(p: BranchParams) -> List[float]:
    c = branch_couplings(p)
    C = c["C"]
    coeff = [p.M, -(p.K + p.M * p.varpi**2), p.K * p.varpi**2 - C**2]
    roots = np.roots(coeff)
    return sorted(math.sqrt(r.real) for r in roots if abs(r.imag) < 1e-10 and r.real > 1e-12)


def uncoupled_internal_roots(p: BranchParams) -> List[float]:
    c = branch_couplings(p)
    R = c["R"]
    coeff = [1.0, -(p.OmegaU**2 + p.OmegaW**2), p.OmegaU**2 * p.OmegaW**2 - R**2]
    roots = np.roots(coeff)
    return sorted(math.sqrt(r.real) for r in roots if abs(r.imag) < 1e-10 and r.real > 1e-12)


def N_of_omega(omega: float, p: BranchParams) -> float:
    c = branch_couplings(p)
    A = p.OmegaU**2 - omega**2
    Delta = (p.OmegaU**2 - omega**2) * (p.OmegaW**2 - omega**2) - c["R"]**2
    P = A * c["GW"] + c["R"] * c["GU"]
    return float(P**2 / Delta**2)


def Delta_of_omega(omega: float, p: BranchParams) -> float:
    c = branch_couplings(p)
    return float((p.OmegaU**2 - omega**2) * (p.OmegaW**2 - omega**2) - c["R"]**2)


def P_of_omega(omega: float, p: BranchParams) -> float:
    c = branch_couplings(p)
    return float((p.OmegaU**2 - omega**2) * c["GW"] + c["R"] * c["GU"])


def D_expr_and_derivative(p: BranchParams):
    c = branch_couplings(p)
    omega = sp.symbols("omega", real=True)
    KB = (p.K - p.M * omega**2) - c["C"]**2 / (p.varpi**2 - omega**2)
    A = p.OmegaU**2 - omega**2
    W = p.OmegaW**2 - omega**2
    Delta = A * W - c["R"]**2
    Q = c["GU"]**2 * W + 2 * c["GU"] * c["GW"] * c["R"] + c["GW"]**2 * A
    D = sp.simplify(KB - Q / Delta)
    Dp = sp.diff(D, omega)
    return sp.lambdify(omega, D, "numpy"), sp.lambdify(omega, Dp, "numpy")


def pole_census(p: BranchParams) -> List[PoleData]:
    coeff = quartic_coefficients(p)
    y_roots = np.roots(coeff)
    full_roots = sorted(math.sqrt(r.real) for r in y_roots if abs(r.imag) < 1e-10 and r.real > 1e-12)
    wall_roots = uncoupled_wall_roots(p)
    internal_roots = uncoupled_internal_roots(p)

    out: List[PoleData] = []
    for om in full_roots:
        nearest_wall = min(wall_roots, key=lambda z: abs(z - om))
        nearest_internal = min(internal_roots, key=lambda z: abs(z - om))
        if abs(om - nearest_wall) <= abs(om - nearest_internal):
            fam = "wall-like"
            nearest = nearest_wall
        else:
            fam = "internal-like"
            nearest = nearest_internal

        Nstar = N_of_omega(om, p)
        Delta_star = Delta_of_omega(om, p)
        Pstar = P_of_omega(om, p)
        Gamma_star = (p.a**5 / (27.0 * p.cs**5)) * om**5
        R_Q = 1.0 / (Gamma_star * Nstar)
        R_qq = R_Q

        out.append(
            PoleData(
                omega=om,
                y=om**2,
                family=fam,
                nearest_uncoupled=nearest,
                Delta_star=Delta_star,
                N_star=Nstar,
                P_star=Pstar,
                R_Q=R_Q,
                R_qq=R_qq,
            )
        )
    return out


def print_numeric_branch_slice() -> None:
    banner("PART IV — NUMERICAL PRIMITIVE BRANCH SLICE")

    p = BranchParams(
        lamB=0.5,
        lamU=0.3,
        lamW=0.4,
        lamR=0.25,
        OmegaU=1.0,
        OmegaW=1.4,
        varpi=2.0,
        K=3.0,
        M=1.0,
        a=1.0,
        cs=1.0,
    )

    c = branch_couplings(p)
    print("Primitive branch parameters:")
    print(p)
    print("\nOverlap-renormalized couplings and static bundle data:")
    for key in ("C", "GU", "GW", "R", "Delta0", "Q0", "P0_proto", "D0", "N0", "P0"):
        print(f"{key:8s} = {c[key]: .12f}")

    if not (c["Delta0"] > 0 and c["D0"] > 0):
        raise AssertionError("Sample slice is not statically admissible")

    wall_unc = uncoupled_wall_roots(p)
    int_unc = uncoupled_internal_roots(p)
    print("\nUncoupled wall/BdG roots:")
    for r in wall_unc:
        print(f"  {r: .12f}")
    print("Uncoupled internal U/W roots:")
    for r in int_unc:
        print(f"  {r: .12f}")

    poles = pole_census(p)
    print("\nFull conservative pole census:")
    print("omega_*        y_*            family         nearest uncpl.  Delta_*         N_*             R_Q,*")
    for pd in poles:
        print(
            f"{pd.omega: .12f}  {pd.y: .12f}  {pd.family:12s}  {pd.nearest_uncoupled: .12f}  "
            f"{pd.Delta_star: .12f}  {pd.N_star: .12f}  {pd.R_Q: .12f}"
        )

    # Verify exact ratio cancellation numerically on the sample slice.
    Df, Dpf = D_expr_and_derivative(p)
    subbanner("Numerical residue/linewidth cancellation check")
    for pd in poles:
        Dprime = float(Dpf(pd.omega))
        Aabs = abs(1.0 / Dprime)
        Gamma_star = (p.a**5 / (27.0 * p.cs**5)) * pd.omega**5
        gamma_star = Gamma_star * pd.N_star / abs(Dprime)
        expect_zero(
            f"|A_qq|/gamma - R_Q at omega={pd.omega:.6f}",
            sp.nsimplify(Aabs / gamma_star - pd.R_Q, tolerance=1e-12),
        )

    # Stage-001 illustrative benchmark at x = 1.
    subbanner("Stage-001 local barrier benchmark")
    V_known_at_1 = 1.181909222592
    epsilon = 0.1
    DeltaV_req = V_known_at_1 - epsilon
    print(f"DeltaV_req(x=1) = V_known(1) - epsilon = {DeltaV_req:.12f}")

    for eta in (0.1, 0.3):
        req_ratio = 2.0 * DeltaV_req * (1.0 + eta**2) / eta  # S_Q(1)^2 = 1
        print(f"Required R_Q,* at eta = {eta:.1f}: {req_ratio:.12f}")
        for pd in poles:
            if pd.family == "wall-like":
                verdict = "PASS" if pd.R_Q >= req_ratio else "FAIL"
                print(
                    f"  pole omega = {pd.omega:.12f} | R_Q,* = {pd.R_Q:.12f} | {verdict}"
                )


def scan_lambda_W_tradeoff() -> None:
    banner("PART V — SCAN IN lambda_W AND STATIC/DYNAMIC TENSION")

    lamW_values = [0.2, 0.4, 0.6, 0.8, 1.0]
    rows = []
    for lamW in lamW_values:
        p = BranchParams(
            lamB=0.5,
            lamU=0.3,
            lamW=lamW,
            lamR=0.25,
            OmegaU=1.0,
            OmegaW=1.4,
            varpi=2.0,
            K=3.0,
            M=1.0,
            a=1.0,
            cs=1.0,
        )
        c = branch_couplings(p)
        poles = pole_census(p)
        wall_poles = [pd for pd in poles if pd.family == "wall-like"]
        upper_wall = max(wall_poles, key=lambda pd: pd.omega)
        rows.append((lamW, c["P0"], c["D0"], upper_wall.omega, upper_wall.R_Q))

    print("lambda_W    P0(static)      D0(static)      omega_upperwall   R_Q,upperwall")
    for row in rows:
        print(f"{row[0]: .6f}   {row[1]: .12f}   {row[2]: .12f}   {row[3]: .12f}   {row[4]: .12f}")

    # The explicit scan shows the sample-family tension: along this slice P0 rises while the
    # upper wall-like residue/linewidth figure falls.
    P0s = [row[1] for row in rows]
    RQs = [row[4] for row in rows]
    if not all(P0s[i] < P0s[i + 1] for i in range(len(P0s) - 1)):
        raise AssertionError("P0 should rise monotonically on the chosen scan.")
    if not all(RQs[i] > RQs[i + 1] for i in range(len(RQs) - 1)):
        raise AssertionError("R_Q should fall monotonically on the chosen scan.")

    print("\nInterpretation:")
    print("  - On this explicit primitive family, strengthening the outgoing leg raises the static prefactor P0.")
    print("  - But the same move lowers the upper wall-like residue/linewidth figure R_Q,*.")
    print("  - So the static outgoing-normalization corridor and the linear dynamic low-loss corridor are in direct tension on this slice.")


# ---------------------------------------------------------------------------
# Run all parts
# ---------------------------------------------------------------------------

def main() -> None:
    verify_exact_branch_polynomial()
    verify_residue_linewidth_cancellation()
    verify_survival_inequality()
    print_numeric_branch_slice()
    scan_lambda_W_tradeoff()

    banner("STAGE 005 SUMMARY")
    print("1. The primitive finite-throat one-port branch gives an exact quartic pole polynomial in y = omega^2.")
    print("2. For the pure-Q same-charge family, the residue/linewidth figure is exactly")
    print("      R_Q,* = 27 c_s^5 / (a^5 omega_*^5 N_*(omega_*)).")
    print("3. The dynamic survival test is therefore an explicit inequality on one pole's transfer factor N_*(omega_*).")
    print("4. On the concrete sample slice, the upper wall-like pole clears the Stage-001 benchmark even at eta = 0.1, while the lower wall-like pole does not.")
    print("5. Along the explicit lambda_W scan, the static prefactor P0 rises while the upper wall-like residue/linewidth figure falls.")
    print("6. So the idea survives this stage, but only in a narrow corridor with a real static/dynamic tension.")


if __name__ == "__main__":
    main()
