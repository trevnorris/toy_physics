#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

_stage224_path = Path(__file__).with_name("5pn_stage224_overlap_image_parameterization.py")
_spec = importlib.util.spec_from_file_location("stage224_overlap_image_parameterization", _stage224_path)
_mod = importlib.util.module_from_spec(_spec)
assert _spec is not None and _spec.loader is not None
_spec.loader.exec_module(_mod)
overlap_image_parameterization = _mod.overlap_image_parameterization
constructive_packet_null_data = _mod.constructive_packet_null_data

"""
Stage 226 — branch-selection theorem from actual moving-throat coherence laws.

What this script does
---------------------
1. Packages the Stage-224 overlap image as seven exact residual equations on the full
   eleven-component overlap-state drift vector.
2. Packages the actual moving-throat coherent branch laws from the compact PDE program as the
   three direct monomial-preservation residuals
      (Sigma_eta, Sigma_tr, Sigma_nt).
3. Proves that the overlap image satisfies those three branch laws identically.
4. Shows that the full overlap-image residual system is exactly equivalent to:
      - the 3 branch-law residuals,
      - plus a 4-equation selector quartet.
5. Computes the exact dimension counts.

Interpretation
--------------
The actual moving-throat coherent branch laws do not kill the packet-null compensation family,
so the compensation route is not forbidden. But those laws also do not determine the full
packet-null image. They underdetermine it by four exact scalar conditions. So the present
branch-selection verdict is:

  not forced,
  not forbidden,
  underdetermined by four selector equations.
"""


def image_residual_vector(
    dlnK,
    dlnM,
    dlnC,
    dlnvarpi,
    dlnOmegaU,
    dlnOmegaW,
    dlnchi0,
    dlnepsilon_eta,
    dlnZW,
    dlnepsilon_W,
    dlndeltaU,
):
    param = overlap_image_parameterization()
    chi = param["coords"][0]
    z = param["coords"][1]
    c = param["coords"][2]
    v = param["coords"][3]
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
    banner("STAGE 226 — BRANCH-SELECTION THEOREM FROM COHERENT MOVING-THROAT LAWS")

    data = constructive_packet_null_data()
    param = overlap_image_parameterization()
    chi0 = sp.simplify(data["chi0"])
    deltaU = sp.simplify(data["deltaU"])
    E_star = sp.simplify(data["E_star"])
    F_star = sp.simplify(data["F_star"])

    # Observed overlap-state tangent.
    dlnK_obs, dlnM_obs, dlnC_obs, dlnvarpi_obs = sp.symbols(
        "dlnK_obs dlnM_obs dlnC_obs dlnvarpi_obs", real=True
    )
    dlnOmegaU_obs, dlnOmegaW_obs = sp.symbols("dlnOmegaU_obs dlnOmegaW_obs", real=True)
    dlnchi0_obs, dlnepsilon_eta_obs = sp.symbols("dlnchi0_obs dlnepsilon_eta_obs", real=True)
    dlnZW_obs, dlnepsilon_W_obs, dlndeltaU_obs = sp.symbols(
        "dlnZW_obs dlnepsilon_W_obs dlndeltaU_obs", real=True
    )
    obs = [
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
    ]

    R = image_residual_vector(*obs)
    R_K, R_M, R_OU, R_OW, R_eta, R_epsW, R_delta = list(R)

    subbanner("I. Full overlap-image residual system")
    print("Residuals are ordered as")
    print("  (R_K, R_M, R_OmegaU, R_OmegaW, R_epsilon_eta, R_epsilon_W, R_deltaU)")
    sp.pprint(R)

    subbanner("II. Actual moving-throat coherent branch-law residuals")
    B_eta = sp.simplify(dlnepsilon_eta_obs)
    B_tr = sp.simplify((1 + chi0) * dlndeltaU_obs + (1 + deltaU) * dlnchi0_obs)
    B_nt = sp.simplify((dlnZW_obs - 2 * dlnOmegaW_obs) + E_star * dlnepsilon_W_obs - F_star * dlndeltaU_obs)
    B = sp.Matrix([B_eta, B_tr, B_nt])

    print("Branch-law residuals are ordered as")
    print("  (B_eta, B_tr, B_nt)")
    sp.pprint(B)

    subbanner("III. Exact compatibility of the overlap image with the branch laws")
    image_subs = {
        dlnK_obs: param["dlnK"],
        dlnM_obs: param["dlnM"],
        dlnC_obs: param["dlnC"],
        dlnvarpi_obs: param["dlnvarpi"],
        dlnOmegaU_obs: param["dlnOmegaU"],
        dlnOmegaW_obs: param["dlnOmegaW"],
        dlnchi0_obs: param["dlnchi0"],
        dlnepsilon_eta_obs: param["dlnepsilon_eta"],
        dlnZW_obs: param["dlnZW"],
        dlnepsilon_W_obs: param["dlnepsilon_W"],
        dlndeltaU_obs: param["dlndeltaU"],
    }
    expect_zero("branch laws on the overlap image", sp.simplify(B.subs(image_subs)))

    subbanner("IV. Exact decomposition into branch laws + selector quartet")
    S = sp.Matrix([R_K, R_M, R_OU, R_epsW])
    print("Selector quartet =")
    sp.pprint(S)

    expect_zero("R_eta - B_eta", sp.simplify(R_eta - B_eta))
    expect_zero("(1+chi_0) R_delta - B_tr", sp.simplify((1 + chi0) * R_delta - B_tr))
    expect_zero("B_nt - (-2 R_OmegaW + E_* R_epsilon_W - F_* R_delta)", sp.simplify(B_nt - (-2 * R_OW + E_star * R_epsW - F_star * R_delta)))

    print("So the full overlap-image residual system is exactly equivalent to")
    print("  (B_eta, B_tr, B_nt) = 0")
    print("together with")
    print("  (R_K, R_M, R_OmegaU, R_epsilon_W) = 0.")
    print("Once those hold, R_delta and R_OmegaW vanish automatically.")

    subbanner("V. Exact dimension counts")
    J_R = R.jacobian(obs)
    J_B = B.jacobian(obs)
    J_BS = sp.Matrix.vstack(B, S).jacobian(obs)
    print("rank(full image residual system) =", J_R.rank())
    print("rank(actual branch-law system) =", J_B.rank())
    print("rank(branch laws + selector quartet) =", J_BS.rank())
    print("dimension of full overlap-state space = 11")
    print("dimension of branch-law manifold = 11 - 3 =", 11 - J_B.rank())
    print("dimension of packet-null overlap image = 11 - 7 =", 11 - J_R.rank())
    print("codimension gap =", J_R.rank() - J_B.rank())

    banner("STAGE 226 LEDGER")
    print("1. The actual moving-throat coherent branch laws from the compact PDE program")
    print("   are exactly the three monomial-preservation residuals")
    print("      (B_eta, B_tr, B_nt).")
    print("2. The Stage-224 overlap image satisfies those three branch laws identically, so the")
    print("   compensation family is not forbidden by the current moving-throat dynamical relations.")
    print("3. But the branch laws are not selective enough to recover the full overlap image.")
    print("   They underdetermine it by four exact scalar conditions.")
    print("4. A convenient selector quartet is")
    print("      (R_K, R_M, R_OmegaU, R_epsilon_W).")
    print("5. Therefore the present branch-selection verdict is: the current dynamical relations")
    print("   do not force the compensation family, they do not kill it, and they underdetermine")
    print("   it by exactly four scalar conditions.")
