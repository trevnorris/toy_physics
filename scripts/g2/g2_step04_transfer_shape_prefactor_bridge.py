#!/usr/bin/env python3
"""
Step 4 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Uses the exact branch composites
       N_* = T^2 R_tr^{B_*},   D = epsilon_eta,
   to show that the raw coherent defect Xi_1 is the logarithmic slope of the
   transfer shape T^2, while the corrected nontracking coordinate is the
   logarithmic slope of N_*.
2. Uses the grouped weak-axisymmetric outgoing-prefactor law to prove
       Xi_1 = P_1 / P_0,
   and on a one-port branch
       Xi_1 = 2 tau_eff.
3. Inserts the exact one-port continuum transfer-shape formula
       T_A^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-epsilon_W)^2]
             = [27 pi^2 G c_s^5 / (20 a^5 c^5)] * (1-epsilon_eta)/R_target,
   and derives two exact microscopic slope ledgers for Xi_1.
4. Combines those identities with the Step-3 quartic target Lambda_1 to give
   the first direct microscopic matching equations for the missing O(f^4)
   anomaly layer.

Interpretation
--------------
This is the first step where the quartic anomaly scalar is tied directly to
actual moving-throat branch data rather than to an abstract quotient tangent.
The remaining g-2 problem is recast as a required logarithmic slope of the
outgoing prefactor / effective transfer shape on the actual branch.
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


def step4_transfer_shape_bridge() -> None:
    banner("STEP 4 — TRANSFER-SHAPE / OUTGOING-PREFACTOR BRIDGE")

    # ------------------------------------------------------------------
    # I. Exact direct branch variables
    # ------------------------------------------------------------------
    subbanner("IV.1 — Exact branch-composite bridge")

    chi, delta = sp.symbols("chi_0 delta_U", positive=True, real=True)
    qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

    Atr = sp.simplify(2 * chi / ((1 + chi) * (1 + delta)))
    Ctr = sp.simplify(
        chi * delta / ((1 + chi) * (1 + delta) * (1 + chi + delta))
    )
    Bstar = sp.simplify(2 * (1 + chi + delta) / delta)

    print("A_tr =")
    sp.pprint(Atr)
    print("C_tr =")
    sp.pprint(Ctr)
    print("B_* =")
    sp.pprint(Bstar)

    expect_zero("A_tr - B_* C_tr", Atr - Bstar * Ctr)

    # Coherent branch formulas.
    dln_Rtr = sp.simplify(-Ctr * qtr)
    Xi1 = sp.simplify(Atr * qtr + qnt)
    dln_N = sp.simplify(qnt)
    dln_T2 = sp.symbols("delta_ln_T2", real=True)

    # N_* = T^2 R_tr^{B_*}  =>  dln N_* = dln T^2 + B_* dln R_tr
    dln_T2_from_N = sp.simplify(dln_N - Bstar * dln_Rtr)

    print("delta ln R_tr =")
    sp.pprint(dln_Rtr)
    print("delta ln N_* =")
    sp.pprint(dln_N)
    print("delta ln T^2 inferred from N_* = T^2 R_tr^{B_*} =")
    sp.pprint(dln_T2_from_N)

    expect_zero("delta ln T^2 - Xi_1", dln_T2_from_N - Xi1)

    print("\nKey identity:")
    print("  raw transfer-shape slope      = Xi_1 = A_tr q_tr + q_nt")
    print("  corrected nontracking slope   = q_nt = delta ln N_*")

    # ------------------------------------------------------------------
    # II. Outgoing prefactor slope
    # ------------------------------------------------------------------
    subbanner("IV.2 — Outgoing-prefactor identity Xi_1 = P_1/P_0")

    eps, lamA = sp.symbols("epsilon lambda_A", real=True)
    P0, P1 = sp.symbols("P_0 P_1", positive=True, real=True)
    T0, T1 = sp.symbols("T_0 T_1", positive=True, real=True)
    tau_eff = sp.symbols("tau_eff", real=True)

    P_A = P0 + eps * lamA * P1
    T2_A = T0 + eps * lamA * T1
    T_A = sp.sqrt(T0) * sp.exp(eps * lamA * tau_eff)

    Xi_from_P = sp.simplify(sp.diff(sp.log(P_A), eps).subs(eps, 0) / lamA)
    Xi_from_T2 = sp.simplify(sp.diff(sp.log(T2_A), eps).subs(eps, 0) / lamA)
    Xi_from_T = sp.simplify(sp.diff(sp.log(T_A**2), eps).subs(eps, 0) / lamA)

    print("Xi_1 from prefactor P_A = P_0 + epsilon lambda_A P_1:")
    sp.pprint(Xi_from_P)
    print("Xi_1 from T_eff^2 = T_0 + epsilon lambda_A T_1:")
    sp.pprint(Xi_from_T2)
    print("Xi_1 on one-port branch T_A = T_0^(1/2) exp(epsilon lambda_A tau_eff):")
    sp.pprint(Xi_from_T)

    expect_zero("Xi_1 - P_1/P_0", Xi_from_P - P1 / P0)
    expect_zero("Xi_1 - T_1/T_0", Xi_from_T2 - T1 / T0)
    expect_zero("Xi_1 - 2 tau_eff", Xi_from_T - 2 * tau_eff)

    # ------------------------------------------------------------------
    # III. Exact one-port microscopic slope ledgers
    # ------------------------------------------------------------------
    subbanner("IV.3 — One-port microscopic slope formulas for Xi_1")

    # Positive variables are varied multiplicatively; rho and epsilons are
    # varied additively because they appear inside (1+rho) and (1-epsilon).
    Z0, Om0 = sp.symbols("Z_W0 Omega_W0", positive=True, real=True)
    rho0 = sp.symbols("rho0", real=True)
    eW0, eEta0 = sp.symbols("epsilonW0 epsilonEta0", real=True)
    cs0, a0, Rt0 = sp.symbols("c_s0 a0 R_target0", positive=True, real=True)
    G, c = sp.symbols("G c", positive=True, real=True)

    sZ, sOm, sCs, sa, sRt = sp.symbols(
        "sigma_Z sigma_Omega sigma_cs sigma_a sigma_R", real=True
    )
    drho, deW, deEta = sp.symbols("delta_rho delta_epsilonW delta_epsilonEta", real=True)
    sRho, sEW, sEEta = sp.symbols(
        "sigma_rho sigma_epsilonW sigma_epsilonEta", real=True
    )

    Z = Z0 * sp.exp(eps * lamA * sZ)
    Om = Om0 * sp.exp(eps * lamA * sOm)
    cs = cs0 * sp.exp(eps * lamA * sCs)
    a = a0 * sp.exp(eps * lamA * sa)
    Rt = Rt0 * sp.exp(eps * lamA * sRt)
    rho = rho0 + eps * lamA * drho
    eW = eW0 + eps * lamA * deW
    eEta = eEta0 + eps * lamA * deEta

    T2_port = sp.simplify(Z * (1 + rho) ** 2 / (Om**2 * (1 - eW) ** 2))
    T2_geom = sp.simplify(
        (27 * sp.pi**2 * G / (20 * c**5)) * cs**5 / a**5 * (1 - eEta) / Rt
    )

    Xi_port = sp.simplify(sp.diff(sp.log(T2_port), eps).subs(eps, 0) / lamA)
    Xi_geom = sp.simplify(sp.diff(sp.log(T2_geom), eps).subs(eps, 0) / lamA)

    Xi_port_sig = sp.simplify(
        Xi_port.subs(
            {
                drho: sRho * (1 + rho0),
                deW: sEW * (1 - eW0),
            }
        )
    )
    Xi_geom_sig = sp.simplify(
        Xi_geom.subs(
            {
                deEta: sEEta * (1 - eEta0),
            }
        )
    )

    print("Xi_1 from T_A^2 = Z_W (1+rho)^2 / [Omega_W^2 (1-epsilon_W)^2] =")
    sp.pprint(Xi_port_sig)
    print("Xi_1 from T_A^2 = [27 pi^2 G c_s^5/(20 a^5 c^5)] * (1-epsilon_eta)/R_target =")
    sp.pprint(Xi_geom_sig)

    expect_zero(
        "Xi_port - (sigma_Z + 2 sigma_rho - 2 sigma_Omega + 2 sigma_epsilonW)",
        Xi_port_sig - (sZ + 2 * sRho - 2 * sOm + 2 * sEW),
    )
    expect_zero(
        "Xi_geom - (5 sigma_cs - 5 sigma_a - sigma_epsilonEta - sigma_R)",
        Xi_geom_sig - (5 * sCs - 5 * sa - sEEta - sRt),
    )

    # ------------------------------------------------------------------
    # IV. Quartic anomaly target in microscopic slope language
    # ------------------------------------------------------------------
    subbanner("IV.4 — Quartic anomaly target in branch and prefactor variables")

    Lam1 = sp.symbols("Lambda_1", real=True)
    s_tr = sp.symbols("s_tr", real=True)

    Xi_direct_match = sp.Eq(P1 / P0, Lam1)
    Xi_corrected_match = sp.Eq(P1 / P0 - Atr * s_tr, Lam1)
    Xi_tracking_rigid = sp.Eq(P1 / P0, Lam1)

    print("Direct observable quartic match:")
    sp.pprint(Xi_direct_match)
    print("Corrected nontracking quartic match:")
    sp.pprint(Xi_corrected_match)
    print("Tracking-rigid branch (s_tr = 0):")
    sp.pprint(Xi_tracking_rigid)

    micro_port_eq = sp.Eq(sZ + 2 * sRho - 2 * sOm + 2 * sEW, Lam1 + Atr * s_tr)
    micro_geom_eq = sp.Eq(5 * sCs - 5 * sa - sEEta - sRt, Lam1 + Atr * s_tr)

    print("\nMicroscopic port-side matching ledger:")
    sp.pprint(micro_port_eq)
    print("Microscopic geometric/selected-side matching ledger:")
    sp.pprint(micro_geom_eq)

    Lam1_num = sp.Float("0.279605891931464", 30)
    print(f"\nStep-3 benchmark value: Lambda_1 = {Lam1_num}")
    print("So on a tracking-rigid branch the required outgoing-prefactor slope is")
    print(f"  P_1/P_0 = Xi_1 = {Lam1_num}")
    print("and the exact one-port microscopic slope laws become")
    print(f"  sigma_Z + 2 sigma_rho - 2 sigma_Omega + 2 sigma_epsilonW = {Lam1_num}")
    print(f"  5 sigma_cs - 5 sigma_a - sigma_epsilonEta - sigma_R = {Lam1_num}")

    subbanner("IV.5 — Interpretation")
    print("1. Step 3's coherent defect Xi_1 is not abstract anymore: it is the")
    print("   logarithmic slope of the effective transfer shape, and equally of")
    print("   the outgoing prefactor, on the actual moving-throat branch.")
    print("2. The corrected nontracking composite N_* removes the universal")
    print("   tracking feed-through exactly. Therefore the quartic anomaly target")
    print("   can be imposed either on Xi_1 directly or on q_nt = Xi_1 - A_tr q_tr.")
    print("3. The first place actual branch data can now enter is explicit: through")
    print("   the microscopic slopes of Z_W, rho, Omega_W, epsilon_W, c_s, a,")
    print("   epsilon_eta, and R_target.")


if __name__ == "__main__":
    step4_transfer_shape_bridge()
