#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 233 — microscopic two-selector theorem and selected-spectral corollary.

What this script does
---------------------
1. Combines the Stage-230 reduced residual vector with the actual minimal coherent
   continuum branch substitution.
2. Proves the exact microscopic form of the reduced residual vector:
      (0, 0, Sigma_sup, -Sigma_eta/8, 0, -Sigma_eta, 0).
3. Shows that the exact overlap-image restoration map therefore has the direct
   microscopic meaning
      dlnOmegaU    -> dlnOmegaU - Sigma_sup,
      dlnepsilon_W -> dlnepsilon_W + Sigma_eta,
      dlnOmegaW    -> dlnOmegaW + Sigma_eta/8.
4. Applies the same selector formulas to the raw selected spectral branch and
   derives its exact two-scalar miss.

Interpretation
--------------
This is the first theorem-level bridge between:
  - the Stage-228/230 overlap-image selector language,
  - the Stage-165/166 microscopic coherent-branch slippages,
  - and the actual selected-spectral moving-throat branch.
"""


def _load_stage230_module():
    stage230_path = Path(__file__).with_name("5pn_stage230_reduced_branch_realizability_tester.py")
    spec = importlib.util.spec_from_file_location(
        "stage230_reduced_branch_realizability_tester", stage230_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Stage 230 module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if __name__ == "__main__":
    banner("STAGE 233 — MICROSCOPIC TWO-SELECTOR THEOREM")

    stage230 = _load_stage230_module()
    data = stage230.branch_law_tester_data()
    obs = data["observed"]

    # Microscopic coherent-branch variables.
    lam1, c1, gam1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
    kapU, kapEta, kapW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)

    Sigma_chi = sp.simplify(gam1 + c1 - kapU)
    zeta_Z = sp.simplify(2 * lam1 - kapEta - kapW)
    Sigma_eta = sp.simplify(2 * c1 - kapU - kapEta)

    coeff_K = sp.simplify(sp.diff(data["OU_pred_from_support"], obs["dlnK"]))
    coeff_chi = sp.simplify(-sp.diff(data["OU_pred_from_support"], obs["dlnchi0"]))
    coeff_Z = sp.simplify(-sp.diff(data["OU_pred_from_support"], obs["dlnZW"]))

    Sigma_sup = sp.simplify(kapU / 2 - coeff_K * kapEta + coeff_chi * Sigma_chi + coeff_Z * zeta_Z)

    subs_actual = {
        obs["dlnK"]: kapEta,
        obs["dlnM"]: sp.Integer(0),
        obs["dlnOmegaU"]: kapU / 2,
        obs["dlnchi0"]: Sigma_chi,
        obs["dlnZW"]: zeta_Z,
        obs["dlnepsilon_W"]: sp.simplify(2 * gam1 + 2 * lam1 - kapU - kapW),
    }

    residual_actual = sp.simplify(data["residual_vector"].subs(subs_actual))
    target_residual = sp.Matrix([
        sp.Integer(0),
        sp.Integer(0),
        Sigma_sup,
        -Sigma_eta / 8,
        sp.Integer(0),
        -Sigma_eta,
        sp.Integer(0),
    ])

    subbanner("I. Exact microscopic residual vector")
    print("Residual vector on the minimal coherent continuum branch =")
    sp.pprint(residual_actual)
    print()
    print("Target microscopic form =")
    sp.pprint(target_residual)
    expect_zero("residual_actual - target_residual", residual_actual - target_residual)

    subbanner("II. Exact microscopic restoration map")
    print("Because the reduced restoration map is")
    print("  dlnOmegaU    -> dlnOmegaU - S_support,")
    print("  dlnepsilon_W -> dlnepsilon_W - S_shape,")
    print("  dlnOmegaW    -> dlnOmegaW - S_shape/8,")
    print("the actual microscopic meaning on this branch is")
    print("  dlnOmegaU    -> dlnOmegaU - Sigma_sup,")
    print("  dlnepsilon_W -> dlnepsilon_W + Sigma_eta,")
    print("  dlnOmegaW    -> dlnOmegaW + Sigma_eta/8.")

    subbanner("III. Raw selected-spectral branch corollary")
    sigma_K_raw, sigma_kappa_raw = sp.symbols("sigma_K_raw sigma_kappa_raw", real=True)
    subs_selected = {
        obs["dlnK"]: sigma_K_raw,
        obs["dlnM"]: sp.Integer(0),
        obs["dlnOmegaU"]: sp.Integer(0),
        obs["dlnchi0"]: sp.Integer(0),
        obs["dlnZW"]: sp.simplify(2 * sigma_kappa_raw - sigma_K_raw),
        obs["dlnepsilon_W"]: sp.simplify(2 * sigma_kappa_raw),
    }
    S_shape_sel = sp.simplify(data["S_shape"].subs(subs_selected))
    S_support_sel = sp.simplify(data["S_support"].subs(subs_selected))
    S_support_sel_target = sp.simplify(2 * coeff_Z * sigma_kappa_raw - (coeff_K + coeff_Z) * sigma_K_raw)

    print("S_shape(selected spectral raw branch) =")
    sp.pprint(S_shape_sel)
    print("S_support(selected spectral raw branch) =")
    sp.pprint(S_support_sel)
    expect_zero("S_shape(selected) - sigma_K_raw", sp.simplify(S_shape_sel - sigma_K_raw))
    expect_zero("S_support(selected) - target", sp.simplify(S_support_sel - S_support_sel_target))

    print()
    print("Numerically,")
    print("  S_support(selected) ≈")
    sp.pprint(sp.N(S_support_sel_target, 16))
    print("  = 3.99575193547 * sigma_kappa_raw - 3.95028944679 * sigma_K_raw")

    subbanner("IV. Meaning of the selected-spectral miss")
    print("The raw selected spectral branch therefore misses the overlap image in exactly two ways:")
    print("  - a shape/dressing miss     S_shape   = sigma_K_raw,")
    print("  - a support-plane miss      S_support = 2 c_Z sigma_kappa_raw - (c_K + c_Z) sigma_K_raw.")
    print("Since the stable selected branch has sigma_K_raw < 0 and sigma_kappa_raw > 0,")
    print("both misses are real and require compensating co-evolution.")

    banner("STAGE 233 LEDGER")
    print("1. On the minimal coherent continuum branch, the full reduced overlap-image residual")
    print("   vector collapses exactly to")
    print("      (0, 0, Sigma_sup, -Sigma_eta/8, 0, -Sigma_eta, 0).")
    print("2. So overlap-image membership is equivalent to exactly two microscopic conditions:")
    print("      Sigma_eta = 0,   Sigma_sup = 0.")
    print("3. The raw selected spectral branch fails these as")
    print("      S_shape   = sigma_K_raw,"
          )
    print("      S_support = 2 c_Z sigma_kappa_raw - (c_K + c_Z) sigma_K_raw.")
    print("4. This is the sharpest current continuation point: the actual moving-throat branch")
    print("   must dynamically generate precisely those two compensations if it is to land in")
    print("   the 5PN packet-null image on the coherent branch-law manifold.")
