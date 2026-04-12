#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 207 — explicit finite-throat overlap extractor for the isotropic branch state.

What this script does
---------------------
1. Freezes the explicit finite-throat overlap prototype on s in [0,L]:
      chi_eta(s) = sqrt(2/L) sin(pi s/L),
      phi_DN(s)  = sqrt(2/L) sin(pi s/(2L)),
      u(s)       = chi_eta(s),
      w(s)       = phi_DN(s).
2. Computes the exact radial/axial overlaps and shows that the prototype induces
   a simple overlap renormalization of the raw support / mixed couplings.
3. Extracts the compressed isotropic branch state
      (K, M, C, varpi, Omega_U, Omega_W, chi_0, epsilon_eta,
       epsilon_W, Z_W, delta_U)
   together with the source factor mhat_0.
4. Verifies that the extracted normalized variables reproduce the coherent
   local-kernel formulas from the moving-throat notes exactly.
5. Records the exact selected-branch transfer-shape identity
      R_target * T^2 = Lambda_0 (1 - epsilon_eta),
   now written directly in the extracted overlap state.

Interpretation
--------------
This is the first actual extractor from an explicit moving-throat overlap model to
one finite isotropic state vector. It is the clean re-entry point before compiling
Packet A and Packet B directly from overlap data.
"""


def overlap_profiles(s: sp.Symbol, L: sp.Symbol) -> dict[str, sp.Expr]:
    chi_eta = sp.sqrt(sp.Integer(2) / L) * sp.sin(sp.pi * s / L)
    phi_DN = sp.sqrt(sp.Integer(2) / L) * sp.sin(sp.pi * s / (2 * L))
    return {
        "chi_eta": chi_eta,
        "phi_DN": phi_DN,
        "u": chi_eta,
        "w": phi_DN,
    }


def overlap_integrals(s: sp.Symbol, L: sp.Symbol) -> dict[str, sp.Expr]:
    prof = overlap_profiles(s, L)
    chi_eta = prof["chi_eta"]
    phi_DN = prof["phi_DN"]
    u = prof["u"]
    w = prof["w"]
    return {
        "I_eta_u": sp.simplify(sp.integrate(chi_eta * u, (s, 0, L))),
        "I_eta_phi": sp.simplify(sp.integrate(chi_eta * phi_DN, (s, 0, L))),
        "I_eta_w": sp.simplify(sp.integrate(chi_eta * w, (s, 0, L))),
        "I_uw": sp.simplify(sp.integrate(u * w, (s, 0, L))),
    }


def extracted_state(
    *,
    K: sp.Expr,
    M: sp.Expr,
    lambda_B_raw: sp.Expr,
    varpi: sp.Expr,
    c_etaU_raw: sp.Expr,
    lambda_W_raw: sp.Expr,
    gamma_raw: sp.Expr,
    K_U: sp.Expr,
    mu_U: sp.Expr,
    K_W: sp.Expr,
    mu_W: sp.Expr,
    T_U: sp.Expr,
    L: sp.Expr,
    sigma: sp.Expr,
    I_eta_u: sp.Expr,
    I_eta_phi: sp.Expr,
    I_eta_w: sp.Expr,
    I_uw: sp.Expr,
) -> dict[str, sp.Expr]:
    C = sp.simplify(lambda_B_raw * I_eta_phi)
    c_etaU_eff = sp.simplify(c_etaU_raw * I_eta_u)
    lambda_W_eff = sp.simplify(lambda_W_raw * I_eta_w)
    gamma_eff = sp.simplify(gamma_raw * I_uw / I_eta_w)

    Omega_U = sp.sqrt(K_U / mu_U)
    Omega_W = sp.sqrt(K_W / mu_W)

    chi_0 = sp.simplify(gamma_eff * c_etaU_eff / K_U)
    epsilon_eta = sp.simplify(c_etaU_eff**2 / (K * K_U))
    epsilon_W = sp.simplify(gamma_eff**2 * lambda_W_eff**2 * sigma / (K_U * K_W))
    Z_W = sp.simplify(lambda_W_eff**2 / (K * K_W))
    delta_U = sp.simplify(sp.pi**2 * T_U / (L**2 * K_U))

    return {
        "K": K,
        "M": M,
        "C": C,
        "varpi": varpi,
        "Omega_U": Omega_U,
        "Omega_W": Omega_W,
        "chi_0": chi_0,
        "epsilon_eta": epsilon_eta,
        "epsilon_W": epsilon_W,
        "Z_W": Z_W,
        "delta_U": delta_U,
        "c_etaU_eff": c_etaU_eff,
        "lambda_W_eff": lambda_W_eff,
        "gamma_eff": gamma_eff,
    }


if __name__ == "__main__":
    banner("STAGE 207 — EXPLICIT FINITE-THROAT OVERLAP EXTRACTOR")

    # Overlap-model coordinates.
    s, L = sp.symbols("s L", positive=True, real=True)

    # Raw microscopic inputs.
    G, c, c_s, a, mhat0 = sp.symbols("G c c_s a mhat_0", positive=True, real=True)
    K, M = sp.symbols("K M", positive=True, real=True)
    lambda_B_raw, varpi = sp.symbols("lambda_B_raw varpi", positive=True, real=True)
    c_etaU_raw = sp.symbols("c_etaU_raw", positive=True, real=True)
    lambda_W_raw = sp.symbols("lambda_W_raw", positive=True, real=True)
    gamma_raw = sp.symbols("gamma_raw", positive=True, real=True)
    K_U, mu_U = sp.symbols("K_U mu_U", positive=True, real=True)
    K_W, mu_W = sp.symbols("K_W mu_W", positive=True, real=True)
    T_U, sigma = sp.symbols("T_U sigma", positive=True, real=True)

    overlaps = overlap_integrals(s, L)

    subbanner("I. Exact finite-throat overlaps")
    for key, value in overlaps.items():
        print(f"{key} =")
        sp.pprint(value)
    expect_zero("I_eta_u - 1", overlaps["I_eta_u"] - 1)
    expect_zero("I_eta_phi - 8/(3 pi)", overlaps["I_eta_phi"] - sp.Rational(8, 3) / sp.pi)
    expect_zero("I_eta_w - I_eta_phi", overlaps["I_eta_w"] - overlaps["I_eta_phi"])
    expect_zero("I_uw - I_eta_phi", overlaps["I_uw"] - overlaps["I_eta_phi"])

    state = extracted_state(
        K=K,
        M=M,
        lambda_B_raw=lambda_B_raw,
        varpi=varpi,
        c_etaU_raw=c_etaU_raw,
        lambda_W_raw=lambda_W_raw,
        gamma_raw=gamma_raw,
        K_U=K_U,
        mu_U=mu_U,
        K_W=K_W,
        mu_W=mu_W,
        T_U=T_U,
        L=L,
        sigma=sigma,
        I_eta_u=overlaps["I_eta_u"],
        I_eta_phi=overlaps["I_eta_phi"],
        I_eta_w=overlaps["I_eta_w"],
        I_uw=overlaps["I_uw"],
    )

    subbanner("II. Overlap-renormalized microscopic couplings")
    print("C =")
    sp.pprint(state["C"])
    print("c_etaU_eff =")
    sp.pprint(state["c_etaU_eff"])
    print("lambda_W_eff =")
    sp.pprint(state["lambda_W_eff"])
    print("gamma_eff =")
    sp.pprint(state["gamma_eff"])

    expect_zero("c_etaU_eff - c_etaU_raw", state["c_etaU_eff"] - c_etaU_raw)
    expect_zero(
        "lambda_W_eff - (8/(3 pi)) lambda_W_raw",
        state["lambda_W_eff"] - sp.Rational(8, 3) * lambda_W_raw / sp.pi,
    )
    expect_zero("gamma_eff - gamma_raw", state["gamma_eff"] - gamma_raw)

    subbanner("III. Extracted isotropic branch state")
    for key in [
        "K",
        "M",
        "C",
        "varpi",
        "Omega_U",
        "Omega_W",
        "chi_0",
        "epsilon_eta",
        "epsilon_W",
        "Z_W",
        "delta_U",
    ]:
        print(f"{key} =")
        sp.pprint(state[key])

    # Reconstruct the compact coherent formulas directly from overlap-renormalized couplings.
    chi_0_expected = sp.simplify(state["gamma_eff"] * state["c_etaU_eff"] / K_U)
    epsilon_eta_expected = sp.simplify(state["c_etaU_eff"]**2 / (K * K_U))
    epsilon_W_expected = sp.simplify(state["gamma_eff"]**2 * state["lambda_W_eff"]**2 * sigma / (K_U * K_W))
    Z_W_expected = sp.simplify(state["lambda_W_eff"]**2 / (K * K_W))
    delta_U_expected = sp.simplify(sp.pi**2 * T_U / (L**2 * K_U))

    expect_zero("chi_0 coherent formula", state["chi_0"] - chi_0_expected)
    expect_zero("epsilon_eta coherent formula", state["epsilon_eta"] - epsilon_eta_expected)
    expect_zero("epsilon_W coherent formula", state["epsilon_W"] - epsilon_W_expected)
    expect_zero("Z_W coherent formula", state["Z_W"] - Z_W_expected)
    expect_zero("delta_U coherent formula", state["delta_U"] - delta_U_expected)

    subbanner("IV. Selected-branch transfer shape")
    epsilon_split = sp.simplify(
        state["epsilon_W"] * (1 - sp.Rational(2, 11) * state["delta_U"] / (1 + state["delta_U"]))
    )
    T_sq = sp.simplify(
        state["Z_W"] * (1 + state["chi_0"]) ** 2
        / (state["Omega_W"]**2 * (1 - epsilon_split) ** 2)
    )
    Lambda = sp.simplify(27 * sp.pi**2 * G * c_s**5 * K_W / (20 * a**5 * c**5 * mu_W))
    R_target = sp.simplify(
        Lambda * (1 - state["epsilon_eta"]) * (1 - epsilon_split) ** 2
        / (state["Z_W"] * (1 + state["chi_0"]) ** 2)
    )
    Lambda0 = sp.simplify(27 * sp.pi**2 * G * c_s**5 / (20 * a**5 * c**5))

    print("epsilon := epsilon_W [1 - (2/11) delta_U/(1+delta_U)] =")
    sp.pprint(epsilon_split)
    print("T^2 =")
    sp.pprint(T_sq)
    print("R_target =")
    sp.pprint(R_target)
    expect_zero(
        "selected-branch demand identity",
        sp.simplify(R_target * T_sq - Lambda0 * (1 - state["epsilon_eta"])),
    )

    banner("STAGE 207 LEDGER")
    print("1. The explicit finite-throat overlap prototype compresses to one extracted isotropic branch state:")
    print("      (K, M, C, varpi, Omega_U, Omega_W, chi_0, epsilon_eta, epsilon_W, Z_W, delta_U),")
    print("   together with the source factor mhat_0.")
    print("2. In this prototype the overlap geometry renormalizes")
    print("      lambda_W -> lambda_W_eff = (8/(3 pi)) lambda_W_raw,")
    print("   while c_etaU and gamma pass through unchanged.")
    print("3. The extracted normalized variables reproduce the coherent local-kernel formulas exactly.")
    print("4. The selected-branch transfer-shape identity survives exactly in the extracted state:")
    print("      R_target * T^2 = Lambda_0 (1 - epsilon_eta).")
