#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 228 — exact reduction of the overlap-image test on the coherent branch-law manifold.

What this script does
---------------------
1. Imports the exact Stage-224 overlap-image parameterization.
2. Rebuilds the seven-component overlap-image residual vector on the full
   eleven-component overlap-state drift space.
3. Restricts that tester to the actual coherent moving-throat branch-law manifold:
      B_eta = 0,
      B_tr  = 0,
      B_nt  = 0.
4. Proves that, after that restriction, the seven residuals collapse to:
      - a support residual triple (R_K, R_M, R_OmegaU),
      - one orbit/shape selector
            S_shape = dln epsilon_W - (2 dln chi_0 + dln Z_W),
      - with the remaining residuals vanishing or becoming exact multiples of S_shape.
5. Records the exact residual vector on the branch-law manifold.

Interpretation
--------------
Stage 226 showed that the full overlap-image test is equivalent to the coherent
branch laws plus a selector quartet. Stage 228 sharpens this further:
inside the coherent branch-law manifold, overlap-image membership already reduces to
    (R_K, R_M, R_OmegaU, S_shape) = 0.
"""


def _load_stage224_module():
    stage224_path = Path(__file__).with_name("5pn_stage224_overlap_image_parameterization.py")
    spec = importlib.util.spec_from_file_location(
        "stage224_overlap_image_parameterization", stage224_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Stage 224 module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def image_residual_vector(
    dlnK: sp.Expr,
    dlnM: sp.Expr,
    dlnC: sp.Expr,
    dlnvarpi: sp.Expr,
    dlnOmegaU: sp.Expr,
    dlnOmegaW: sp.Expr,
    dlnchi0: sp.Expr,
    dlnepsilon_eta: sp.Expr,
    dlnZW: sp.Expr,
    dlnepsilon_W: sp.Expr,
    dlndeltaU: sp.Expr,
) -> sp.Matrix:
    stage224 = _load_stage224_module()
    param = stage224.overlap_image_parameterization()
    chi, z, c, v = param["coords"]
    pred_subs = {chi: dlnchi0, z: dlnZW, c: dlnC, v: dlnvarpi}

    R_K = sp.simplify(dlnK - param["dlnK"].subs(pred_subs))
    R_M = sp.simplify(dlnM - param["dlnM"].subs(pred_subs))
    R_OU = sp.simplify(dlnOmegaU - param["dlnOmegaU"].subs(pred_subs))
    R_OW = sp.simplify(dlnOmegaW - param["dlnOmegaW"].subs(pred_subs))
    R_eta = sp.simplify(dlnepsilon_eta - param["dlnepsilon_eta"].subs(pred_subs))
    R_epsW = sp.simplify(dlnepsilon_W - param["dlnepsilon_W"].subs(pred_subs))
    R_delta = sp.simplify(dlndeltaU - param["dlndeltaU"].subs(pred_subs))
    return sp.Matrix([R_K, R_M, R_OU, R_OW, R_eta, R_epsW, R_delta])


if __name__ == "__main__":
    banner("STAGE 228 — EXACT BRANCH-LAW REDUCTION OF THE OVERLAP-IMAGE TEST")

    stage224 = _load_stage224_module()
    data = stage224.constructive_packet_null_data()
    chi0 = sp.simplify(data["chi0"])
    deltaU = sp.simplify(data["deltaU"])
    E_star = sp.simplify(data["E_star"])
    F_star = sp.simplify(data["F_star"])

    dlnK_obs, dlnM_obs, dlnC_obs, dlnvarpi_obs = sp.symbols(
        "dlnK_obs dlnM_obs dlnC_obs dlnvarpi_obs", real=True
    )
    dlnOmegaU_obs, dlnOmegaW_obs = sp.symbols("dlnOmegaU_obs dlnOmegaW_obs", real=True)
    dlnchi0_obs, dlnepsilon_eta_obs = sp.symbols("dlnchi0_obs dlnepsilon_eta_obs", real=True)
    dlnZW_obs, dlnepsilon_W_obs, dlndeltaU_obs = sp.symbols(
        "dlnZW_obs dlnepsilon_W_obs dlndeltaU_obs", real=True
    )

    R = image_residual_vector(
        dlnK_obs,
        dlnM_obs,
        dlnC_obs,
        dlnvarpi_obs,
        dlnOmegaU_obs,
        dlnOmegaW_obs,
        dlnchi0_obs,
        dlnepsilon_eta_obs,
        dlnZW_obs,
        dlnepsilon_W_obs,
        dlndeltaU_obs,
    )
    R_K, R_M, R_OU, R_OW, R_eta, R_epsW, R_delta = list(R)

    subbanner("I. Actual coherent branch-law manifold")
    print("B_eta = dlnepsilon_eta")
    print("B_tr  = (1+chi_0) dlndeltaU + (1+delta_U) dlnchi0")
    print("B_nt  = (dlnZW - 2 dlnOmegaW) + E_* dlnepsilon_W - F_* dlndeltaU")
    print()
    print("Constructive branch constants:")
    print("chi_0   =", chi0)
    print("delta_U =", deltaU)
    print("E_*     =", E_star)
    print("F_*     =", F_star)

    dlndeltaU_branch = sp.simplify(-(1 + deltaU) / (1 + chi0) * dlnchi0_obs)
    dlnOmegaW_branch = sp.simplify(
        (dlnZW_obs + E_star * dlnepsilon_W_obs - F_star * dlndeltaU_branch) / 2
    )
    branch_subs = {
        dlnepsilon_eta_obs: sp.Integer(0),
        dlndeltaU_obs: dlndeltaU_branch,
        dlnOmegaW_obs: dlnOmegaW_branch,
    }

    subbanner("II. Reduced residual vector on the branch-law manifold")
    R_branch = sp.Matrix([sp.simplify(expr.subs(branch_subs)) for expr in R])
    sp.pprint(R_branch)

    S_shape = sp.simplify(dlnepsilon_W_obs - (2 * dlnchi0_obs + dlnZW_obs))

    expect_zero("R_eta on branch-law manifold", R_branch[4])
    expect_zero("R_deltaU on branch-law manifold", R_branch[6])
    expect_zero("R_OmegaW - S_shape/8", sp.simplify(R_branch[3] - S_shape / 8))
    expect_zero("R_epsilonW - S_shape", sp.simplify(R_branch[5] - S_shape))

    subbanner("III. Exact orbit/shape selector")
    print("S_shape =")
    sp.pprint(S_shape)
    print()
    print("So on the coherent branch-law manifold,")
    print("  R_OmegaW = S_shape / 8,")
    print("  R_epsilonW = S_shape,")
    print("while R_eta = R_deltaU = 0 identically.")

    subbanner("IV. Remaining support residual triple")
    print("(R_K, R_M, R_OmegaU) =")
    sp.pprint(sp.Matrix([R_branch[0], R_branch[1], R_branch[2]]))
    expect_zero("d R_K / dln epsilon_W", sp.diff(R_branch[0], dlnepsilon_W_obs))
    expect_zero("d R_M / dln epsilon_W", sp.diff(R_branch[1], dlnepsilon_W_obs))
    expect_zero("d R_OmegaU / dln epsilon_W", sp.diff(R_branch[2], dlnepsilon_W_obs))
    print("So the remaining support residuals are completely blind to the orbit/shape selector S_shape.")

    banner("STAGE 228 LEDGER")
    print("1. Restricting the overlap-image tester to the actual coherent branch-law manifold")
    print("   kills R_eta and R_deltaU identically.")
    print("2. The whole orbit/shape mismatch collapses to one exact selector")
    print("      S_shape = dlnepsilon_W - (2 dlnchi_0 + dln Z_W).")
    print("3. On that manifold, R_OmegaW is not independent: it is exactly S_shape / 8.")
    print("4. The remaining independent image-membership conditions are therefore")
    print("      (R_K, R_M, R_OmegaU, S_shape) = 0.")
    print("5. This is already sharper than the Stage-226 selector quartet because the coherent")
    print("   branch laws have absorbed the (R_eta, R_deltaU, R_OmegaW) side completely.")
