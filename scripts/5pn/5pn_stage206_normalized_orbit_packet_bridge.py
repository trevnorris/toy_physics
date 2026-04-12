#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import *

"""
Stage 206 — normalized coherent orbit-packet bridge.

What this script does
---------------------
1. Packages the Stage-12/13 monomials directly in the normalized coherent
   variables used by the primitive overlap model.
2. Builds the exact invariant-ratio packet
      (R_tr, R_nt, R_eta)
   and the equivalent quotient-coordinate packet
      (q_tr, q_nt, q_eta).
3. Re-derives the exact normalized monomial-drift matrix from first principles
   and verifies that it matches the Stage-13 formula.
4. Solves the zero-defect linear system triangularly, giving the exact required
   co-drifts of (delta_U, M, Omega_W) once a candidate branch supplies the drifts
   of (G_W, G_U, R, K, Omega_U).

Interpretation
--------------
This is the clean splice between the Packet-B orbit theorem (Stages 181–201) and
the normalized coherent variables from Stages 12–13. It gives a direct way to
extract Delta_orbit from the same microscopic variable set used by the Packet-A
prototype.
"""


def log_drift_matrix(monomials: list[sp.Expr], variables: list[sp.Symbol]) -> sp.Matrix:
    return sp.Matrix(
        [
            [sp.simplify(v * sp.diff(sp.log(m), v)) for v in variables]
            for m in monomials
        ]
    )


if __name__ == "__main__":
    banner("STAGE 206 — NORMALIZED COHERENT ORBIT-PACKET BRIDGE")

    # Normalized coherent variables and reference data.
    GW, GU, R, K, M, OmegaU, OmegaW, deltaU = sp.symbols(
        "G_W G_U R K M Omega_U Omega_W delta_U",
        positive=True,
        real=True,
    )
    sigma = sp.symbols("sigma", positive=True, real=True)

    chi0_star, deltaU_star, E_star, F_star = sp.symbols(
        "chi0_star deltaU_star E_star F_star",
        positive=True,
        real=True,
    )

    # Actual branch invariants.
    chi0 = sp.simplify(R * GU / (OmegaU**2 * GW))
    eps_eta = sp.simplify(M * GU**2 / (K * OmegaU**2))
    eps_W = sp.simplify(R**2 * sigma / (OmegaU**2 * OmegaW**2))
    C_tr = sp.simplify(chi0 ** (1 + deltaU_star) * deltaU ** (1 + chi0_star))
    C_nt = sp.simplify((M * GW**2 / (K * OmegaW**4)) * eps_W**E_star * deltaU ** (-F_star))
    C_eta = eps_eta

    # Reference branch data.
    chi0_ref, deltaU_ref, eps_eta_ref, epsW_ref, ZW_ref, OmegaW_ref = sp.symbols(
        "chi0_ref deltaU_ref epsilon_eta_ref epsilon_W_ref ZW_ref OmegaW_ref",
        positive=True,
        real=True,
    )

    C_tr_ref = sp.simplify(chi0_ref ** (1 + deltaU_star) * deltaU_ref ** (1 + chi0_star))
    C_nt_ref = sp.simplify((ZW_ref / OmegaW_ref**2) * epsW_ref**E_star * deltaU_ref ** (-F_star))
    C_eta_ref = eps_eta_ref

    R_tr = sp.simplify(C_tr / C_tr_ref)
    R_nt = sp.simplify(C_nt / C_nt_ref)
    R_eta = sp.simplify(C_eta / C_eta_ref)
    q_packet = q_from_invariants(R_tr, R_nt, R_eta)

    subbanner("I. Exact normalized monomials and invariant ratios")
    print("chi_0 =")
    sp.pprint(chi0)
    print("epsilon_eta =")
    sp.pprint(eps_eta)
    print("epsilon_W =")
    sp.pprint(eps_W)
    print("C_tr =")
    sp.pprint(C_tr)
    print("C_nt =")
    sp.pprint(C_nt)
    print("C_eta =")
    sp.pprint(C_eta)
    print("R_tr =")
    sp.pprint(R_tr)
    print("R_nt =")
    sp.pprint(R_nt)
    print("R_eta =")
    sp.pprint(R_eta)

    subbanner("II. Exact quotient-coordinate packet")
    print("q_tr =")
    sp.pprint(q_packet["qtr"])
    print("q_nt =")
    sp.pprint(q_packet["qnt"])
    print("q_eta =")
    sp.pprint(q_packet["qeta"])

    # Verify the roundtrip through the Stage-200 common compiler utilities.
    chi0s = chi0_star
    Fstar = F_star
    m_from_q = mismatch_from_q(chi0s, Fstar, q_packet["qtr"], q_packet["qnt"], q_packet["qeta"])
    q_from_m = q_from_mismatch(chi0s, Fstar, m_from_q["mT"], m_from_q["mK"], m_from_q["mMu"])

    expect_zero("roundtrip q_tr", sp.simplify(q_from_m["qtr"] - q_packet["qtr"]))
    expect_zero("roundtrip q_nt", sp.simplify(q_from_m["qnt"] - q_packet["qnt"]))
    expect_zero("roundtrip q_eta", sp.simplify(q_from_m["qeta"] - q_packet["qeta"]))

    subbanner("III. Exact normalized monomial-drift matrix")
    variables = [GW, GU, R, K, M, OmegaU, OmegaW, deltaU]
    monomials = [C_tr, C_nt, C_eta]
    M_norm = log_drift_matrix(monomials, variables)
    print("M_norm =")
    sp.pprint(M_norm)

    M_expected = sp.Matrix(
        [
            [-(1 + deltaU_star), 1 + deltaU_star, 1 + deltaU_star, 0, 0, -2 * (1 + deltaU_star), 0, 1 + chi0_star],
            [2, 0, 2 * E_star, -1, 1, -2 * E_star, -(4 + 2 * E_star), -F_star],
            [0, 2, 0, -1, 1, -2, 0, 0],
        ]
    )
    expect_zero("M_norm - expected Stage-13 matrix", M_norm - M_expected)

    subbanner("IV. Exact triangular zero-defect solve")
    dlnGW, dlnGU, dlnR, dlnK, dlnM, dlnOmegaU, dlnOmegaW, dlndeltaU = sp.symbols(
        "dlnGW dlnGU dlnR dlnK dlnM dlnOmegaU dlnOmegaW dlndeltaU",
        real=True,
    )
    drift_vec = sp.Matrix([dlnGW, dlnGU, dlnR, dlnK, dlnM, dlnOmegaU, dlnOmegaW, dlndeltaU])
    eqs = list(M_norm * drift_vec)

    dlndeltaU_sol = sp.simplify(
        sp.solve(sp.Eq(eqs[0], 0), dlndeltaU)[0]
    )
    dlnM_sol = sp.simplify(
        sp.solve(sp.Eq(eqs[2], 0), dlnM)[0]
    )
    dlnOmegaW_sol = sp.simplify(
        sp.solve(sp.Eq(eqs[1].subs({dlndeltaU: dlndeltaU_sol, dlnM: dlnM_sol}), 0), dlnOmegaW)[0]
    )

    print("dln delta_U =")
    sp.pprint(dlndeltaU_sol)
    print("dln M =")
    sp.pprint(dlnM_sol)
    print("dln Omega_W =")
    sp.pprint(dlnOmegaW_sol)

    expect_zero(
        "tracking equation",
        sp.simplify(eqs[0].subs({dlndeltaU: dlndeltaU_sol})),
    )
    expect_zero(
        "dressing equation",
        sp.simplify(eqs[2].subs({dlnM: dlnM_sol})),
    )
    expect_zero(
        "nontracking equation",
        sp.simplify(eqs[1].subs({dlndeltaU: dlndeltaU_sol, dlnM: dlnM_sol, dlnOmegaW: dlnOmegaW_sol})),
    )

    banner("STAGE 206 LEDGER")
    print("1. Packet B is now explicit in the normalized coherent variable set")
    print("      (G_W, G_U, R, K, M, Omega_U, Omega_W, delta_U).")
    print("2. The exact invariant packet is")
    print("      (R_tr, R_nt, R_eta) = (C_tr/C_tr,ref, C_nt/C_nt,ref, epsilon_eta/epsilon_eta,ref),")
    print("   with quotient coordinates (q_tr,q_nt,q_eta) = (ln R_tr, ln R_nt, ln R_eta).")
    print("3. The direct monomial-drift matrix exactly reproduces the Stage-13 normalized kernel.")
    print("4. Zero-defect tangency is triangular: once the branch gives drifts of")
    print("      (G_W, G_U, R, K, Omega_U),")
    print("   the required co-drifts of (delta_U, M, Omega_W) are fixed exactly.")
