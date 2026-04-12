#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 208 — isotropic Packet-A / Packet-B compiler from the extracted overlap state.

What this script does
---------------------
1. Takes the Stage-207 extracted isotropic branch state
      (K, M, C, varpi, Omega_U, Omega_W, chi_0, epsilon_eta,
       epsilon_W, Z_W, delta_U)
   together with mhat_0.
2. Compiles the exact isotropic Packet-A branch defects
      Delta_pole, Delta_norm
   directly from that state.
3. Compiles the exact Packet-B orbit data
      (C_tr, C_nt, C_eta), (R_tr, R_nt, R_eta), (q_tr, q_nt, q_eta)
   from the same extracted state.
4. Proves the clean packet-separation theorem on the isotropic overlap branch:
      - Packet A is blind to (epsilon_W, delta_U),
      - Packet B is blind to the explicit support pair (C, varpi)
        and to the wall inertia pair (M, Omega_U) once the extracted invariants
        are formed.
5. Freezes the corrected minimal data ledger:
      the combined isotropic endgame needs 11 dynamic branch scalars plus mhat_0.

Interpretation
--------------
This is the first direct compiler from actual overlap-model state variables to the
Stage-199/200/201 endgame packets. It sharpens the roadmap from “extract some
branch data” to “extract this exact 11-scalar isotropic state and one source factor.”
"""


def packet_a_operator_moments(
    *,
    K: sp.Expr,
    M: sp.Expr,
    C: sp.Expr,
    varpi: sp.Expr,
    Omega_U: sp.Expr,
    Omega_W: sp.Expr,
    chi_0: sp.Expr,
    epsilon_eta: sp.Expr,
    Z_W: sp.Expr,
) -> dict[str, sp.Expr]:
    epsilon_mix = sp.simplify(Z_W * chi_0**2 / epsilon_eta)
    Delta = sp.simplify(Omega_U**2 * Omega_W**2 * (1 - epsilon_mix))

    G_U = sp.sqrt(epsilon_eta * K * Omega_U**2 / M)
    G_W = sp.sqrt(Z_W * K * Omega_W**2 / M)
    R = sp.simplify(chi_0 * Omega_U**2 * G_W / G_U)
    Q = sp.simplify(G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2)
    H = sp.simplify(G_U**2 + G_W**2)
    S2 = sp.simplify(Omega_U**2 + Omega_W**2)
    P = sp.simplify(Omega_U**2 * G_W + R * G_U)

    B0 = sp.simplify(C**2 / varpi**2)
    B2 = sp.simplify(C**2 / varpi**4)
    B4 = sp.simplify(C**2 / varpi**6)

    Z0 = sp.simplify(Q / Delta)
    Z2 = sp.simplify((Q * S2 - H * Delta) / Delta**2)
    Z4 = sp.simplify((Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3)

    N0 = sp.simplify(P**2 / Delta**2)
    N2 = sp.simplify(2 * P * (P * S2 - Delta * G_W) / Delta**3)
    N4 = sp.simplify(
        (Delta**2 * G_W**2 - 2 * Delta * P**2 - 4 * Delta * P * S2 * G_W + 3 * P**2 * S2**2)
        / Delta**4
    )

    D0 = sp.simplify(K - B0 - Z0)
    D2 = sp.simplify(-(M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    return {
        "epsilon_mix": epsilon_mix,
        "Delta": Delta,
        "B0": B0,
        "B2": B2,
        "B4": B4,
        "Z0": Z0,
        "Z2": Z2,
        "Z4": Z4,
        "N0": N0,
        "N2": N2,
        "N4": N4,
        "D0": D0,
        "D2": D2,
        "D4": D4,
    }


def packet_b_invariants(
    *,
    chi_0: sp.Expr,
    delta_U: sp.Expr,
    epsilon_eta: sp.Expr,
    epsilon_W: sp.Expr,
    Z_W: sp.Expr,
    Omega_W: sp.Expr,
    chi0_star: sp.Expr,
    deltaU_star: sp.Expr,
    E_star: sp.Expr,
    F_star: sp.Expr,
    chi0_ref: sp.Expr,
    deltaU_ref: sp.Expr,
    epsilon_eta_ref: sp.Expr,
    epsilon_W_ref: sp.Expr,
    ZW_ref: sp.Expr,
    OmegaW_ref: sp.Expr,
) -> dict[str, sp.Expr]:
    C_tr = sp.simplify(chi_0 ** (1 + deltaU_star) * delta_U ** (1 + chi0_star))
    C_nt = sp.simplify((Z_W / Omega_W**2) * epsilon_W**E_star * delta_U ** (-F_star))
    C_eta = sp.simplify(epsilon_eta)

    C_tr_ref = sp.simplify(chi0_ref ** (1 + deltaU_star) * deltaU_ref ** (1 + chi0_star))
    C_nt_ref = sp.simplify((ZW_ref / OmegaW_ref**2) * epsilon_W_ref**E_star * deltaU_ref ** (-F_star))
    C_eta_ref = sp.simplify(epsilon_eta_ref)

    R_tr = sp.simplify(C_tr / C_tr_ref)
    R_nt = sp.simplify(C_nt / C_nt_ref)
    R_eta = sp.simplify(C_eta / C_eta_ref)

    q_tr = sp.simplify(sp.log(R_tr))
    q_nt = sp.simplify(sp.log(R_nt))
    q_eta = sp.simplify(sp.log(R_eta))

    return {
        "C_tr": C_tr,
        "C_nt": C_nt,
        "C_eta": C_eta,
        "R_tr": R_tr,
        "R_nt": R_nt,
        "R_eta": R_eta,
        "q_tr": q_tr,
        "q_nt": q_nt,
        "q_eta": q_eta,
    }


def free_symbol_names(exprs: list[sp.Expr]) -> list[str]:
    names: set[str] = set()
    for expr in exprs:
        names |= {s.name for s in expr.free_symbols}
    return sorted(names)


if __name__ == "__main__":
    banner("STAGE 208 — ISOTROPIC OVERLAP PACKET COMPILER")

    # Extracted isotropic state from Stage 207.
    G, c, c_s, a, mhat0 = sp.symbols("G c c_s a mhat_0", positive=True, real=True)
    K, M, C, varpi = sp.symbols("K M C varpi", positive=True, real=True)
    Omega_U, Omega_W = sp.symbols("Omega_U Omega_W", positive=True, real=True)
    chi_0 = sp.symbols("chi_0", positive=True, real=True)
    epsilon_eta = sp.symbols("epsilon_eta", positive=True, real=True)
    epsilon_W = sp.symbols("epsilon_W", positive=True, real=True)
    Z_W = sp.symbols("Z_W", positive=True, real=True)
    delta_U = sp.symbols("delta_U", positive=True, real=True)

    # Reference-branch / orbit data.
    chi0_star, deltaU_star, E_star, F_star = sp.symbols(
        "chi0_star deltaU_star E_star F_star", positive=True, real=True
    )
    chi0_ref, deltaU_ref = sp.symbols("chi0_ref deltaU_ref", positive=True, real=True)
    epsilon_eta_ref, epsilon_W_ref = sp.symbols(
        "epsilon_eta_ref epsilon_W_ref", positive=True, real=True
    )
    ZW_ref, OmegaW_ref = sp.symbols("ZW_ref OmegaW_ref", positive=True, real=True)

    A = packet_a_operator_moments(
        K=K,
        M=M,
        C=C,
        varpi=varpi,
        Omega_U=Omega_U,
        Omega_W=Omega_W,
        chi_0=chi_0,
        epsilon_eta=epsilon_eta,
        Z_W=Z_W,
    )

    Delta_pole = sp.simplify(-(3 * A["D2"]**2 + A["D0"] * A["D4"]) / A["D0"]**2)
    P0 = sp.simplify(A["N0"] / A["D0"])
    Delta_norm = sp.simplify(mhat0**2 * P0 - 54 * G * c_s**5 / (5 * a**5 * c**5))

    B = packet_b_invariants(
        chi_0=chi_0,
        delta_U=delta_U,
        epsilon_eta=epsilon_eta,
        epsilon_W=epsilon_W,
        Z_W=Z_W,
        Omega_W=Omega_W,
        chi0_star=chi0_star,
        deltaU_star=deltaU_star,
        E_star=E_star,
        F_star=F_star,
        chi0_ref=chi0_ref,
        deltaU_ref=deltaU_ref,
        epsilon_eta_ref=epsilon_eta_ref,
        epsilon_W_ref=epsilon_W_ref,
        ZW_ref=ZW_ref,
        OmegaW_ref=OmegaW_ref,
    )

    subbanner("I. Exact isotropic Packet-A branch defects")
    print("D0 =")
    sp.pprint(A["D0"])
    print("D2 =")
    sp.pprint(A["D2"])
    print("D4 =")
    sp.pprint(A["D4"])
    print("Delta_pole =")
    sp.pprint(Delta_pole)
    print("Delta_norm =")
    sp.pprint(Delta_norm)

    subbanner("II. Exact Packet-B orbit data from the same extracted state")
    for key in ["C_tr", "C_nt", "C_eta", "R_tr", "R_nt", "R_eta", "q_tr", "q_nt", "q_eta"]:
        print(f"{key} =")
        sp.pprint(B[key])

    subbanner("III. Packet-separation theorem on the isotropic overlap branch")
    packet_a_dynamic = [A["D0"], A["D2"], A["D4"], A["N0"], A["N2"], A["N4"], Delta_pole, P0]
    packet_b_dynamic = [B["C_tr"], B["C_nt"], B["C_eta"], B["R_tr"], B["R_nt"], B["R_eta"], B["q_tr"], B["q_nt"], B["q_eta"]]

    packet_a_names = free_symbol_names(packet_a_dynamic)
    packet_b_names = free_symbol_names(packet_b_dynamic)
    print("Packet-A free symbols =")
    sp.pprint(packet_a_names)
    print("Packet-B free symbols =")
    sp.pprint(packet_b_names)

    # Clean blindness checks.
    packet_a_blind = {"epsilon_W", "delta_U"}
    packet_b_blind = {"C", "varpi", "M", "Omega_U"}

    if any(name in packet_a_blind for name in packet_a_names):
        raise AssertionError("Packet A unexpectedly depends on epsilon_W or delta_U.")
    if any(name in packet_b_blind for name in packet_b_names):
        raise AssertionError("Packet B unexpectedly depends on C, varpi, M, or Omega_U.")

    print("Packet A is blind to (epsilon_W, delta_U).")
    print("Packet B is blind to (C, varpi, M, Omega_U).")

    subbanner("IV. Corrected minimal isotropic data ledger")
    print("Packet A sufficiency set:")
    print("  (K, M, C, varpi, Omega_U, Omega_W, chi_0, epsilon_eta, Z_W) plus mhat_0")
    print("Packet B sufficiency set:")
    print("  (chi_0, delta_U, epsilon_eta, epsilon_W, Z_W, Omega_W)")
    print("Combined dynamic isotropic state:")
    print("  (K, M, C, varpi, Omega_U, Omega_W, chi_0, epsilon_eta, epsilon_W, Z_W, delta_U)")
    print("  plus the source factor mhat_0.")

    banner("STAGE 208 LEDGER")
    print("1. The isotropic Stage-199 branch packet reduces exactly to the two scalar defects")
    print("      Delta_pole and Delta_norm.")
    print("2. The orbit packet is compiled from the same extracted overlap state by")
    print("      (C_tr, C_nt, C_eta) -> (R_tr, R_nt, R_eta) -> (q_tr, q_nt, q_eta).")
    print("3. Packet A and Packet B do not need the same microscopic information.")
    print("4. The corrected combined isotropic endgame state is 11 dynamic branch scalars")
    print("      (K, M, C, varpi, Omega_U, Omega_W, chi_0, epsilon_eta, epsilon_W, Z_W, delta_U)")
    print("   together with mhat_0.")
    print("5. So the next moving-throat theorem gate is no longer vague: extract exactly this state")
    print("   from the actual branch, then evaluate Delta_pole, Delta_norm, and the orbit packet.")
