#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import *

"""
Stage 205 — normalized coherent isotropic Packet-A bridge.

What this script does
---------------------
1. Rewrites the Stage-203 minimal isotropic Packet-A model in the normalized
   coherent variables used later in the 5PN notes:
      chi_0, epsilon_eta, Z_W,
   together with the physical scale data
      (K, M, C, varpi, Omega_U, Omega_W).
2. Derives exact closed forms for the isotropic conservative/outgoing packet
      (B0,B2,B4,Z0,Z2,Z4,N0,N2,N4,D0,D2,D4)
   in that normalized language.
3. Verifies an exact separation theorem:
      - the outgoing-transfer moments N_n are support-blind with respect to
        the explicit BdG support pair (C,varpi),
      - the conservative wall operator still depends on that support pair.
4. Produces the exact normalized compatibility surface obtained by combining
   the 5PN one-pole identity with the 2.5PN/4PN normalization gate.

Interpretation
--------------
This is the clean splice between the Packet-A prototype and the normalized
coherent variables from Stages 12–13 of the 5PN notes. It tells us exactly
which microscopic sectors feed the isotropic Packet-A test.
"""


def minimal_isotropic_bundle(
    K: sp.Expr,
    M: sp.Expr,
    C: sp.Expr,
    varpi: sp.Expr,
    GU: sp.Expr,
    GW: sp.Expr,
    R: sp.Expr,
    OmegaU: sp.Expr,
    OmegaW: sp.Expr,
) -> dict[str, sp.Expr]:
    Delta = sp.simplify(OmegaU**2 * OmegaW**2 - R**2)
    Q = sp.simplify(GU**2 * OmegaW**2 + 2 * GU * GW * R + GW**2 * OmegaU**2)
    P = sp.simplify(OmegaU**2 * GW + R * GU)
    H = sp.simplify(GU**2 + GW**2)
    S2 = sp.simplify(OmegaU**2 + OmegaW**2)

    B0 = sp.simplify(C**2 / varpi**2)
    B2 = sp.simplify(C**2 / varpi**4)
    B4 = sp.simplify(C**2 / varpi**6)

    Z0 = sp.simplify(Q / Delta)
    Z2 = sp.simplify((Q * S2 - H * Delta) / Delta**2)
    Z4 = sp.simplify((Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3)

    N0 = sp.simplify(P**2 / Delta**2)
    N2 = sp.simplify(2 * P * (P * S2 - Delta * GW) / Delta**3)
    N4 = sp.simplify(
        (Delta**2 * GW**2 - 2 * Delta * P**2 - 4 * Delta * P * S2 * GW + 3 * P**2 * S2**2)
        / Delta**4
    )

    D0 = sp.simplify(K - B0 - Z0)
    D2 = sp.simplify(-(M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    return {
        "Delta": Delta,
        "Q": Q,
        "P": P,
        "H": H,
        "S2": S2,
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


def normalized_substitution(
    *,
    K: sp.Expr,
    M: sp.Expr,
    OmegaU: sp.Expr,
    OmegaW: sp.Expr,
    chi0: sp.Expr,
    eps_eta: sp.Expr,
    ZW: sp.Expr,
) -> dict[sp.Expr, sp.Expr]:
    GU = sp.sqrt(eps_eta * K * OmegaU**2 / M)
    GW = sp.sqrt(ZW * K * OmegaW**2 / M)
    R = sp.simplify(chi0 * OmegaU**2 * GW / GU)
    return {sp.Symbol("G_U"): GU, sp.Symbol("G_W"): GW, sp.Symbol("R"): R}


if __name__ == "__main__":
    banner("STAGE 205 — NORMALIZED COHERENT ISOTROPIC PACKET-A BRIDGE")

    # Physical scales kept explicit.
    G, c, cs, a, mhat0 = sp.symbols("G c c_s a mhat_0", positive=True, real=True)
    K, M = sp.symbols("K M", positive=True, real=True)
    C, varpi = sp.symbols("C varpi", positive=True, real=True)
    OmegaU, OmegaW = sp.symbols("Omega_U Omega_W", positive=True, real=True)

    # Normalized coherent variables from the later 5PN notes.
    chi0 = sp.symbols("chi_0", real=True)
    eps_eta = sp.symbols("epsilon_eta", positive=True, real=True)
    ZW = sp.symbols("Z_W", positive=True, real=True)
    eps_mix = sp.symbols("epsilon_mix", real=True)

    GU = sp.sqrt(eps_eta * K * OmegaU**2 / M)
    GW = sp.sqrt(ZW * K * OmegaW**2 / M)
    R = sp.simplify(chi0 * OmegaU**2 * GW / GU)

    bundle = minimal_isotropic_bundle(K, M, C, varpi, GU, GW, R, OmegaU, OmegaW)
    moms = response_moments_from_D(bundle["D0"], bundle["D2"], bundle["D4"])
    pref = prefactor_moments(bundle["D0"], bundle["D2"], bundle["D4"], bundle["N0"], bundle["N2"], bundle["N4"])

    subbanner("I. Exact normalized coherent dictionary")
    print("G_U =")
    sp.pprint(GU)
    print("G_W =")
    sp.pprint(GW)
    print("R =")
    sp.pprint(R)

    expect_zero(
        "chi_0 recovered",
        sp.simplify(R * GU / (OmegaU**2 * GW) - chi0),
    )
    expect_zero(
        "epsilon_eta recovered",
        sp.simplify(M * GU**2 / (K * OmegaU**2) - eps_eta),
    )
    expect_zero(
        "Z_W recovered",
        sp.simplify(M * GW**2 / (K * OmegaW**2) - ZW),
    )

    subbanner("II. Exact isotropic Packet-A coefficients in normalized variables")
    eps_expr = sp.simplify(ZW * chi0**2 / eps_eta)
    Delta_expected = sp.simplify(OmegaU**2 * OmegaW**2 * (1 - eps_expr))
    Z0_expected = sp.simplify(K * (eps_eta + ZW * (1 + 2 * chi0)) / (M * (1 - eps_expr)))
    N0_expected = sp.simplify(K * ZW * (1 + chi0) ** 2 / (M * OmegaW**2 * (1 - eps_expr) ** 2))

    print("epsilon_mix := chi_0^2 Z_W / epsilon_eta =")
    sp.pprint(eps_expr)
    print("Delta =")
    sp.pprint(bundle["Delta"])
    print("Z0 =")
    sp.pprint(bundle["Z0"])
    print("N0 =")
    sp.pprint(bundle["N0"])

    expect_zero("Delta normalized form", sp.simplify(bundle["Delta"] - Delta_expected))
    expect_zero("Z0 normalized form", sp.simplify(bundle["Z0"] - Z0_expected))
    expect_zero("N0 normalized form", sp.simplify(bundle["N0"] - N0_expected))

    Sigma2 = sp.simplify(bundle["Z2"] / (K / M))
    Sigma4 = sp.simplify(bundle["Z4"] / (K / M))
    print("Sigma2 := (M/K) Z2 =")
    sp.pprint(Sigma2)
    print("Sigma4 := (M/K) Z4 =")
    sp.pprint(Sigma4)

    subbanner("III. Support-blind outgoing-transfer theorem on the isotropic branch")
    expect_zero("dN0/dC", sp.diff(bundle["N0"], C))
    expect_zero("dN0/dvarpi", sp.diff(bundle["N0"], varpi))
    expect_zero("dN2/dC", sp.diff(bundle["N2"], C))
    expect_zero("dN2/dvarpi", sp.diff(bundle["N2"], varpi))
    expect_zero("dN4/dC", sp.diff(bundle["N4"], C))
    expect_zero("dN4/dvarpi", sp.diff(bundle["N4"], varpi))
    print("dD0/dC =")
    sp.pprint(sp.simplify(sp.diff(bundle["D0"], C)))
    print("dD0/dvarpi =")
    sp.pprint(sp.simplify(sp.diff(bundle["D0"], varpi)))

    subbanner("IV. Packet-A moments after normalized substitution")
    print("u2 =")
    sp.pprint(moms["u2"])
    print("u4 =")
    sp.pprint(moms["u4"])
    print("P0 =")
    sp.pprint(pref["P0"])

    subbanner("V. Exact normalized compatibility surface")
    P0_target = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5 * mhat0**2))
    compatibility = sp.simplify(bundle["N0"] / P0_target - 3 * (M + bundle["B2"] + bundle["Z2"]) ** 2 / (bundle["B4"] + bundle["Z4"]))
    print("P0_target =")
    sp.pprint(P0_target)
    print("Compatibility surface =")
    sp.pprint(sp.factor(compatibility))

    # A slightly more structural form.
    structural_compat = sp.simplify(
        K * ZW * (1 + chi0) ** 2 / (M * OmegaW**2 * (1 - eps_expr) ** 2 * P0_target)
        - 3 * (M + C**2 / varpi**4 + (K / M) * Sigma2) ** 2 / (C**2 / varpi**6 + (K / M) * Sigma4)
    )
    expect_zero("normalized compatibility equivalence", sp.simplify(structural_compat - compatibility))

    banner("STAGE 205 LEDGER")
    print("1. The isotropic Packet-A prototype is now written in the normalized coherent variables")
    print("      chi_0, epsilon_eta, Z_W, together with the scale data (K,M,C,varpi,Omega_U,Omega_W).")
    print("2. The mixed conservative/outgoing sector enters through the exact combination")
    print("      epsilon_mix = chi_0^2 Z_W / epsilon_eta,")
    print("   so Delta = Omega_U^2 Omega_W^2 (1 - epsilon_mix).")
    print("3. The outgoing-transfer moments N0,N2,N4 are exactly support-blind in the explicit")
    print("   BdG support pair (C,varpi), while the conservative wall operator D0,D2,D4 is not.")
    print("4. The isotropic Packet-A theorem gate is therefore an explicit compatibility surface in")
    print("   the normalized coherent variables plus the support pair, not a diffuse PDE statement.")
