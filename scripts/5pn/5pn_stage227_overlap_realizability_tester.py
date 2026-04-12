#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

_stage226_path = Path(__file__).with_name("5pn_stage226_branch_selection_theorem.py")
_spec = importlib.util.spec_from_file_location("stage226_branch_selection_theorem", _stage226_path)
_mod = importlib.util.module_from_spec(_spec)
assert _spec is not None and _spec.loader is not None
_spec.loader.exec_module(_mod)
image_residual_vector = _mod.image_residual_vector

_stage224_path = Path(__file__).with_name("5pn_stage224_overlap_image_parameterization.py")
_spec2 = importlib.util.spec_from_file_location("stage224_overlap_image_parameterization", _stage224_path)
_mod2 = importlib.util.module_from_spec(_spec2)
assert _spec2 is not None and _spec2.loader is not None
_spec2.loader.exec_module(_mod2)
overlap_image_parameterization = _mod2.overlap_image_parameterization

"""
Stage 227 — exact overlap-branch realizability tester.

What this script does
---------------------
1. Packages the Stage-224 overlap image into a concrete seven-residual tester.
2. Uses the overlap-extracted first-order drift state
      (dlnK, dlnM, dlnC, dlnvarpi, dlnOmegaU, dlnOmegaW,
       dlnchi0, dlnepsilon_eta, dlnZW, dlnepsilon_W, dlndeltaU)
   as its input.
3. Returns the exact seven residuals that must vanish for the observed tangent to lie in the
   constructive support-restored packet-null family.
4. Verifies exact vanishing on the internal Stage-224 image.

Interpretation
--------------
This is the drop-in interface for any later expensive moving-throat overlap computation.
Once a PDE / overlap solve returns a candidate first-order tangent in the overlap variables,
this tester gives an immediate yes/no residual against the current 4D packet-null family.
"""


def realizability_distance_squared(residuals: sp.Matrix) -> sp.Expr:
    return sp.simplify(sum(sp.expand(r**2) for r in residuals))


if __name__ == "__main__":
    banner("STAGE 227 — EXACT OVERLAP-BRANCH REALIZABILITY TESTER")

    # Symbolic observed overlap-state tangent.
    dlnK_obs, dlnM_obs, dlnC_obs, dlnvarpi_obs = sp.symbols(
        "dlnK_obs dlnM_obs dlnC_obs dlnvarpi_obs", real=True
    )
    dlnOU_obs, dlnOW_obs = sp.symbols("dlnOmegaU_obs dlnOmegaW_obs", real=True)
    dlnchi_obs, dlneta_obs = sp.symbols("dlnchi0_obs dlnepsilon_eta_obs", real=True)
    dlnZW_obs, dlnepsW_obs, dlndelta_obs = sp.symbols(
        "dlnZW_obs dlnepsilon_W_obs dlndeltaU_obs", real=True
    )

    residuals = image_residual_vector(
        dlnK_obs,
        dlnM_obs,
        dlnC_obs,
        dlnvarpi_obs,
        dlnOU_obs,
        dlnOW_obs,
        dlnchi_obs,
        dlneta_obs,
        dlnZW_obs,
        dlnepsW_obs,
        dlndelta_obs,
    )

    subbanner("I. Exact realizability residual vector")
    print("Residuals are ordered as")
    print("  (R_K, R_M, R_OmegaU, R_OmegaW, R_epsilon_eta, R_epsilon_W, R_deltaU)")
    sp.pprint(residuals)

    dist2 = realizability_distance_squared(residuals)
    print("Exact quadratic realizability distance D_real^2 =")
    sp.pprint(dist2)

    subbanner("II. Internal self-check on the Stage-224 image")
    param = overlap_image_parameterization()
    image_subs = {
        dlnK_obs: param["dlnK"],
        dlnM_obs: param["dlnM"],
        dlnC_obs: param["dlnC"],
        dlnvarpi_obs: param["dlnvarpi"],
        dlnOU_obs: param["dlnOmegaU"],
        dlnOW_obs: param["dlnOmegaW"],
        dlnchi_obs: param["dlnchi0"],
        dlneta_obs: param["dlnepsilon_eta"],
        dlnZW_obs: param["dlnZW"],
        dlnepsW_obs: param["dlnepsilon_W"],
        dlndelta_obs: param["dlndeltaU"],
    }
    expect_zero("realizability residuals on the internal overlap image", sp.simplify(residuals.subs(image_subs)))
    expect_zero("realizability distance on the internal overlap image", sp.simplify(dist2.subs(image_subs)))

    subbanner("III. Practical usage")
    print("To test an actual moving-throat overlap tangent, provide the eleven extracted drifts")
    print("  (dlnK, dlnM, dlnC, dlnvarpi, dlnOmegaU, dlnOmegaW, dlnchi0, dlnepsilon_eta, dlnZW, dlnepsilon_W, dlndeltaU)")
    print("and evaluate the seven residuals above.")
    print("The tangent lies in the current constructive support-restored packet-null family iff")
    print("all seven residuals vanish.")

    banner("STAGE 227 LEDGER")
    print("1. The overlap-level realizability test is now explicit and finite-dimensional.")
    print("2. It uses eleven directly extractable overlap-state drifts as input and returns seven exact residuals.")
    print("3. Zero residual means the observed overlap tangent lands inside the current 4D support-restored packet-null family.")
    print("4. This is the right drop-in tester for any later expensive moving-throat branch computation.")
