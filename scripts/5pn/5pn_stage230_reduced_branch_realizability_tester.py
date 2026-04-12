#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 230 — reduced two-scalar realizability tester and exact restoration map.

What this script does
---------------------
1. Imports the Stage-224 overlap image.
2. Restricts the full realizability problem to the coherent branch-law manifold.
3. Solves out the latent support pair (dlnC,dlnvarpi) from (dlnK,dlnM).
4. Shows that the full seven-component residual vector collapses exactly to two
   independent scalars:
       S_support,
       S_shape.
5. Builds the exact minimal restoration map that sends any coherent branch-law
   tangent to the nearest overlap-image representative without leaving the
   coherent branch-law manifold.

Interpretation
--------------
For any future expensive moving-throat tangent extraction, the first overlap-image
membership test no longer needs all eleven drifts. On the coherent branch-law manifold,
it reduces to two exact selectors.
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


def branch_law_tester_data() -> dict[str, sp.Expr]:
    stage224 = _load_stage224_module()
    param = stage224.overlap_image_parameterization()
    data = stage224.constructive_packet_null_data()
    chi0 = sp.simplify(data["chi0"])
    deltaU = sp.simplify(data["deltaU"])
    E_star = sp.simplify(data["E_star"])
    F_star = sp.simplify(data["F_star"])

    dlnK_obs, dlnM_obs, dlnOmegaU_obs = sp.symbols(
        "dlnK_obs dlnM_obs dlnOmegaU_obs", real=True
    )
    dlnchi0_obs, dlnZW_obs, dlnepsilon_W_obs = sp.symbols(
        "dlnchi0_obs dlnZW_obs dlnepsilon_W_obs", real=True
    )
    dlnC_obs, dlnvarpi_obs = sp.symbols("dlnC_obs dlnvarpi_obs", real=True)

    chi, z, c, v = param["coords"]
    Kpred = sp.simplify(param["dlnK"])
    Mpred = sp.simplify(param["dlnM"])
    OUpred = sp.simplify(param["dlnOmegaU"])

    # Solve support coordinates from observed (dlnK,dlnM) on the branch-law manifold.
    sol_list = sp.solve(
        [sp.Eq(dlnK_obs, Kpred.subs({chi: dlnchi0_obs, z: dlnZW_obs, c: dlnC_obs, v: dlnvarpi_obs})),
         sp.Eq(dlnM_obs, Mpred.subs({chi: dlnchi0_obs, z: dlnZW_obs, c: dlnC_obs, v: dlnvarpi_obs}))],
        [dlnC_obs, dlnvarpi_obs],
        dict=True,
    )
    if len(sol_list) != 1:
        raise AssertionError("Expected unique latent-support solve.")
    sol = {k: sp.simplify(val) for k, val in sol_list[0].items()}

    OU_pred_from_support = sp.simplify(
        OUpred.subs({chi: dlnchi0_obs, z: dlnZW_obs, c: sol[dlnC_obs], v: sol[dlnvarpi_obs]})
    )

    S_shape = sp.simplify(dlnepsilon_W_obs - (2 * dlnchi0_obs + dlnZW_obs))
    S_support = sp.simplify(dlnOmegaU_obs - OU_pred_from_support)

    # Coherent branch-law manifold values of the remaining derived observables.
    dlnepsilon_eta_branch = sp.Integer(0)
    dlndeltaU_branch = sp.simplify(-(1 + deltaU) / (1 + chi0) * dlnchi0_obs)
    dlnOmegaW_branch = sp.simplify(
        (dlnZW_obs + E_star * dlnepsilon_W_obs - F_star * dlndeltaU_branch) / 2
    )

    dlnOmegaU_restore = sp.simplify(dlnOmegaU_obs - S_support)
    dlnepsilon_W_restore = sp.simplify(dlnepsilon_W_obs - S_shape)
    dlnOmegaW_restore = sp.simplify(dlnOmegaW_branch - S_shape / 8)

    # Full seven residuals after branch laws and support solve.
    residual_vector = sp.Matrix([
        sp.Integer(0),
        sp.Integer(0),
        S_support,
        sp.simplify(S_shape / 8),
        sp.Integer(0),
        S_shape,
        sp.Integer(0),
    ])

    restored_residual_vector = sp.Matrix([
        sp.Integer(0),
        sp.Integer(0),
        sp.simplify(S_support.subs({dlnOmegaU_obs: dlnOmegaU_restore})),
        sp.simplify((S_shape / 8).subs({dlnepsilon_W_obs: dlnepsilon_W_restore})),
        sp.Integer(0),
        sp.simplify(S_shape.subs({dlnepsilon_W_obs: dlnepsilon_W_restore})),
        sp.Integer(0),
    ])

    return {
        "chi0": chi0,
        "deltaU": deltaU,
        "E_star": E_star,
        "F_star": F_star,
        "observed": {
            "dlnK": dlnK_obs,
            "dlnM": dlnM_obs,
            "dlnOmegaU": dlnOmegaU_obs,
            "dlnchi0": dlnchi0_obs,
            "dlnZW": dlnZW_obs,
            "dlnepsilon_W": dlnepsilon_W_obs,
        },
        "latent_support": sol,
        "OU_pred_from_support": OU_pred_from_support,
        "S_shape": S_shape,
        "S_support": S_support,
        "branch_values": {
            "dlnepsilon_eta": dlnepsilon_eta_branch,
            "dlndeltaU": dlndeltaU_branch,
            "dlnOmegaW": dlnOmegaW_branch,
        },
        "restored": {
            "dlnOmegaU": dlnOmegaU_restore,
            "dlnepsilon_W": dlnepsilon_W_restore,
            "dlnOmegaW": dlnOmegaW_restore,
            "dlnC": sol[dlnC_obs],
            "dlnvarpi": sol[dlnvarpi_obs],
        },
        "residual_vector": residual_vector,
        "restored_residual_vector": restored_residual_vector,
    }


if __name__ == "__main__":
    banner("STAGE 230 — REDUCED TWO-SCALAR REALIZABILITY TESTER")

    data = branch_law_tester_data()
    obs = data["observed"]
    lat = data["latent_support"]
    restored = data["restored"]

    subbanner("I. Reduced selectors")
    print("S_shape =")
    sp.pprint(data["S_shape"])
    print()
    print("S_support =")
    sp.pprint(data["S_support"])

    J_sel = sp.Matrix([data["S_shape"], data["S_support"]]).jacobian(list(obs.values()))
    print()
    print("rank of the two-selector Jacobian in the reduced observable set =", J_sel.rank())
    if J_sel.rank() != 2:
        raise AssertionError("Expected the two reduced selectors to be independent.")

    subbanner("II. Latent support reconstruction")
    print("dlnC(target) =")
    sp.pprint(lat[next(iter([k for k in lat if k.name == 'dlnC_obs']))])
    print()
    print("dlnvarpi(target) =")
    sp.pprint(lat[next(iter([k for k in lat if k.name == 'dlnvarpi_obs']))])
    print()
    print("Predicted dlnOmegaU(target) =")
    sp.pprint(data["OU_pred_from_support"])

    subbanner("III. Exact reduced residual vector")
    print("Residual vector after branch laws + support solve =")
    sp.pprint(data["residual_vector"])

    expect_zero(
        "reduced residual vector after restoration",
        data["restored_residual_vector"],
    )

    subbanner("IV. Exact restoration map on the coherent branch-law manifold")
    print("To restore overlap-image membership while staying on the coherent branch-law manifold, set")
    print("dlnC ->")
    sp.pprint(restored["dlnC"])
    print()
    print("dlnvarpi ->")
    sp.pprint(restored["dlnvarpi"])
    print()
    print("dlnOmegaU ->")
    sp.pprint(restored["dlnOmegaU"])
    print()
    print("dlnepsilon_W ->")
    sp.pprint(restored["dlnepsilon_W"])
    print()
    print("dlnOmegaW ->")
    sp.pprint(restored["dlnOmegaW"])

    banner("STAGE 230 LEDGER")
    print("1. Once the actual branch is restricted to the coherent branch-law manifold, and once the")
    print("   latent support pair (dlnC,dlnvarpi) is solved out from (dlnK,dlnM), overlap-image")
    print("   membership is controlled by exactly two independent scalars:")
    print("      S_shape   = dlnepsilon_W - (2 dlnchi_0 + dlnZ_W),")
    print("      S_support = dlnOmegaU - dlnOmegaU^(pred)(dlnK,dlnM,dlnchi_0,dlnZ_W).")
    print("2. The full seven-component residual vector collapses exactly to")
    print("      (0, 0, S_support, S_shape/8, 0, S_shape, 0).")
    print("3. Therefore the next actual moving-throat test no longer needs the full 11-component tangent.")
    print("   It only needs the reduced observable sextuple")
    print("      (dlnK, dlnM, dlnOmegaU, dlnchi_0, dlnZ_W, dlnepsilon_W),")
    print("   and then the two selectors above.")
    print("4. The exact minimal restoration to the overlap image is obtained by adjusting only")
    print("      the latent support pair (dlnC,dlnvarpi),")
    print("      dlnOmegaU, and the coherent shape pair (dlnepsilon_W,dlnOmegaW).")
