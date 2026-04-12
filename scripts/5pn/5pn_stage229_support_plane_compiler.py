#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 229 — exact support-plane compiler on the coherent branch-law manifold.

What this script does
---------------------
1. Imports the exact Stage-224 overlap-image parameterization.
2. Works on the coherent branch-law manifold isolated in Stage 228.
3. Solves the support residual pair (R_K,R_M)=0 exactly for the latent support
   coordinates (dlnC, dlnvarpi).
4. Substitutes that solve into R_OmegaU and derives the exact codimension-1
   support compatibility plane
        S_support = 0.
5. Prints exact formulas and useful numerical approximations.

Interpretation
--------------
Once the actual moving-throat branch is already coherent-law clean, and once the
latent support pair is allowed to adjust, overlap-image realizability reduces to
one exact support-plane scalar.
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


if __name__ == "__main__":
    banner("STAGE 229 — EXACT SUPPORT-PLANE COMPILER")

    stage224 = _load_stage224_module()
    param = stage224.overlap_image_parameterization()

    chi, z, c, v = param["coords"]
    Kobs, Mobs, OUobs = sp.symbols("dlnK_obs dlnM_obs dlnOmegaU_obs", real=True)

    dlnK_pred = sp.simplify(param["dlnK"])
    dlnM_pred = sp.simplify(param["dlnM"])
    dlnOU_pred = sp.simplify(param["dlnOmegaU"])

    subbanner("I. Exact support residual system")
    R_K = sp.simplify(Kobs - dlnK_pred)
    R_M = sp.simplify(Mobs - dlnM_pred)
    R_OU = sp.simplify(OUobs - dlnOU_pred)
    print("R_K =")
    sp.pprint(R_K)
    print()
    print("R_M =")
    sp.pprint(R_M)
    print()
    print("R_OmegaU =")
    sp.pprint(R_OU)

    support_matrix = sp.Matrix([
        [sp.diff(dlnK_pred, c), sp.diff(dlnK_pred, v)],
        [sp.diff(dlnM_pred, c), sp.diff(dlnM_pred, v)],
        [sp.diff(dlnOU_pred, c), sp.diff(dlnOU_pred, v)],
    ])
    print()
    print("Support coefficient matrix in (dlnC,dlnvarpi) =")
    sp.pprint(support_matrix)
    print("rank =", support_matrix.rank())
    if support_matrix.rank() != 2:
        raise AssertionError("Expected the support coefficient matrix to have rank 2.")

    subbanner("II. Exact support solve from (dlnK,dlnM)")
    sol_list = sp.solve([sp.Eq(Kobs, dlnK_pred), sp.Eq(Mobs, dlnM_pred)], [c, v], dict=True)
    if len(sol_list) != 1:
        raise AssertionError("Expected a unique support solve from (dlnK,dlnM).")
    sol = {k: sp.simplify(val) for k, val in sol_list[0].items()}

    print("dlnC =")
    sp.pprint(sol[c])
    print()
    print("dlnvarpi =")
    sp.pprint(sol[v])

    expect_zero("R_K after support solve", sp.simplify(R_K.subs(sol)))
    expect_zero("R_M after support solve", sp.simplify(R_M.subs(sol)))

    subbanner("III. Exact support compatibility plane")
    OU_from_KM = sp.simplify(dlnOU_pred.subs(sol))
    S_support = sp.simplify(OUobs - OU_from_KM)

    print("Predicted dlnOmegaU =")
    sp.pprint(OU_from_KM)
    print()
    print("S_support =")
    sp.pprint(S_support)

    expect_zero("R_OmegaU - S_support after support solve", sp.simplify(R_OU.subs(sol) - S_support))

    subbanner("IV. Numerical form")
    print("dlnC ≈")
    sp.pprint(sp.N(sol[c], 16))
    print()
    print("dlnvarpi ≈")
    sp.pprint(sp.N(sol[v], 16))
    print()
    print("dlnOmegaU(pred) ≈")
    sp.pprint(sp.N(OU_from_KM, 16))
    print()
    print("S_support ≈")
    sp.pprint(sp.N(S_support, 16))

    banner("STAGE 229 LEDGER")
    print("1. On the coherent branch-law manifold, the latent support pair (dlnC,dlnvarpi) is fixed")
    print("   uniquely by the observed support-side pair (dlnK,dlnM) once (dlnchi_0,dlnZ_W) are given.")
    print("2. The remaining support-side realizability test is codimension-1:")
    print("      S_support = dlnOmegaU - dlnOmegaU^(pred)(dlnK,dlnM,dlnchi_0,dlnZ_W) = 0.")
    print("3. So after solving out the latent support pair, overlap-image membership no longer depends")
    print("   on three support residuals. It depends on one exact support-plane scalar.")
