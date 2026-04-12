#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from functools import lru_cache
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 224 — exact overlap-image parameterization of the constructive packet-null family.

What this script does
---------------------
1. Reuses the Stage-218 constructive support-restored packet master matrix.
2. Replays the Stage-221 overlap-state realization solve on the constructive branch
   E_* = 1/4, F_* = 5/6.
3. Converts the compensated primitive overlap drifts into the overlap-observable drift set
      (dlnK, dlnM, dlnC, dlnvarpi, dlnOmegaU, dlnOmegaW,
       dlnchi0, dlnepsilon_eta, dlnZW, dlnepsilon_W, dlndeltaU).
4. Proves that the image has exact dimension 4 even though it starts from five primitive
   overlap-side drifts.
5. Chooses the branch-adapted overlap-image coordinates
      (dlnchi0, dlnZW, dlnC, dlnvarpi)
   and rewrites the whole image exactly in those four coordinates.

Interpretation
--------------
This is the first reusable finite-dimensional overlap-level parameterization of the
current 4D packet-null family. It is the right bridge between the symbolic packet-null
solve and any future expensive moving-throat overlap extraction.
"""


def load_stage218_builder():
    stage218_path = Path(__file__).with_name("5pn_stage218_support_restored_master_matrix.py")
    spec = importlib.util.spec_from_file_location(
        "stage218_support_restored_master_matrix", stage218_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Stage 218 builder.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_master_system


@lru_cache(maxsize=1)
def constructive_packet_null_data() -> dict[str, sp.Expr | dict[sp.Symbol, sp.Expr] | list[sp.Symbol]]:
    """Return the constructive overlap-side packet-null solve and the induced overlap drifts."""
    build_master_system = load_stage218_builder()
    E_star = sp.Rational(1, 4)
    F_star = sp.Rational(5, 6)
    data = build_master_system(E_star, F_star)

    alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi = data["free"]
    u2_1 = data["u2_1"]
    u4_1 = data["u4_1"]
    Xi_load = data["Xi_load"]
    chi0 = sp.simplify(data["chi0"])
    deltaU = sp.Integer(1)

    # Primitive overlap-state drifts.
    sigma_K, sigma_kappa = sp.symbols("sigma_K sigma_kappa", real=True)
    ell_B, ell_U, ell_W, ell_R = sp.symbols("ell_B ell_U ell_W ell_R", real=True)
    ell_OU, ell_varpi = sp.symbols("ell_OmegaU ell_varpi", real=True)

    carrier_subs = {
        alpha_K: sigma_K,
        alpha_GW: ell_W + sigma_kappa,
        alpha_GU: ell_U,
        alpha_R: ell_R + sigma_kappa,
        alpha_OU: ell_OU,
        beta_B: ell_B + sigma_kappa,
        beta_varpi: ell_varpi,
    }

    u2_overlap = sp.simplify(u2_1.subs(carrier_subs))
    u4_overlap = sp.simplify(u4_1.subs(carrier_subs))
    Xi_overlap = sp.simplify(Xi_load.subs(carrier_subs))

    sol_list = sp.solve(
        [sp.Eq(u2_overlap, 0), sp.Eq(u4_overlap, 0), sp.Eq(Xi_overlap, 0)],
        [ell_U, ell_OU, ell_varpi],
        dict=True,
    )
    if len(sol_list) != 1:
        raise AssertionError("Expected a unique realization solve for (ell_U, ell_OmegaU, ell_varpi).")
    sol = {k: sp.simplify(v) for k, v in sol_list[0].items()}

    alpha_GW_eff = sp.simplify(ell_W + sigma_kappa)
    alpha_GU_eff = sp.simplify(sol[ell_U])
    alpha_R_eff = sp.simplify(ell_R + sigma_kappa)
    alpha_OU_eff = sp.simplify(sol[ell_OU])
    alpha_K_eff = sigma_K
    beta_B_eff = sp.simplify(ell_B + sigma_kappa)
    beta_varpi_eff = sp.simplify(sol[ell_varpi])

    dln_deltaU = sp.simplify(
        -(1 + deltaU) / (1 + chi0) * (alpha_R_eff + alpha_GU_eff - alpha_GW_eff - 2 * alpha_OU_eff)
    )
    dln_M = sp.simplify(alpha_K_eff - 2 * alpha_GU_eff + 2 * alpha_OU_eff)
    dln_OW = sp.simplify(
        (
            alpha_GW_eff
            - alpha_GU_eff
            + (1 - E_star) * alpha_OU_eff
            + E_star * alpha_R_eff
            - (F_star / 2) * dln_deltaU
        )
        / (E_star + 2)
    )

    overlap_drifts = {
        "dlnK": sp.simplify(alpha_K_eff),
        "dlnM": sp.simplify(dln_M),
        "dlnC": sp.simplify(beta_B_eff),
        "dlnvarpi": sp.simplify(beta_varpi_eff),
        "dlnOmegaU": sp.simplify(alpha_OU_eff),
        "dlnOmegaW": sp.simplify(dln_OW),
        "dlnchi0": sp.simplify(alpha_R_eff + alpha_GU_eff - alpha_GW_eff - 2 * alpha_OU_eff),
        "dlnepsilon_eta": sp.simplify(dln_M + 2 * alpha_GU_eff - alpha_K_eff - 2 * alpha_OU_eff),
        "dlnZW": sp.simplify(dln_M + 2 * alpha_GW_eff - alpha_K_eff - 2 * dln_OW),
        "dlnepsilon_W": sp.simplify(2 * (alpha_R_eff - alpha_OU_eff - dln_OW)),
        "dlndeltaU": sp.simplify(dln_deltaU),
    }

    return {
        "E_star": E_star,
        "F_star": F_star,
        "chi0": chi0,
        "deltaU": deltaU,
        "primitive": [sigma_K, sigma_kappa, ell_B, ell_W, ell_R],
        "compensated": sol,
        "overlap_drifts": overlap_drifts,
    }


@lru_cache(maxsize=1)
def overlap_image_parameterization() -> dict[str, sp.Expr | sp.Matrix | list[sp.Symbol]]:
    """Return the exact 4D overlap-image parameterization in branch-adapted coordinates.

    Free image coordinates are
      (dlnchi0, dlnZW, dlnC, dlnvarpi).
    """
    data = constructive_packet_null_data()
    primitive = data["primitive"]
    overlap_drifts = data["overlap_drifts"]

    dlnchi0, dlnZW, dlnC, dlnvarpi = sp.symbols("dlnchi0 dlnZW dlnC dlnvarpi", real=True)
    coords = [dlnchi0, dlnZW, dlnC, dlnvarpi]

    def row_coeffs(expr: sp.Expr) -> sp.Matrix:
        return sp.Matrix([[sp.simplify(sp.diff(expr, p)) for p in primitive]])

    coord_rows = [row_coeffs(overlap_drifts[name]) for name in ("dlnchi0", "dlnZW", "dlnC", "dlnvarpi")]
    coord_matrix = sp.Matrix(coord_rows)
    coord_rank = coord_matrix.rank()
    if coord_rank != 4:
        raise AssertionError(f"Expected overlap image coordinate rank 4, got {coord_rank}.")
    nulls = coord_matrix.nullspace()
    if len(nulls) != 1:
        raise AssertionError("Expected a one-dimensional primitive null direction.")
    primitive_null = sp.Matrix([sp.simplify(x) for x in nulls[0]])

    param_exprs: dict[str, sp.Expr] = {}
    coeff_matrix_rows: list[list[sp.Expr]] = []

    coeff_syms = sp.symbols("a0:4", real=True)
    for name, expr in overlap_drifts.items():
        target_row = row_coeffs(expr)
        combo = sum((coeff_syms[i] * coord_rows[i] for i in range(4)), sp.zeros(1, len(primitive)))
        sol_list = sp.solve(
            [sp.Eq(combo[0, j], target_row[0, j]) for j in range(len(primitive))],
            coeff_syms,
            dict=True,
        )
        if len(sol_list) != 1:
            raise AssertionError(f"Expected unique coefficient solve for {name}.")
        coeffs = [sp.simplify(sol_list[0][s]) for s in coeff_syms]
        coeff_matrix_rows.append(coeffs)
        param_exprs[name] = sp.simplify(sum(c * q for c, q in zip(coeffs, coords)))

        # Exact linear-image check: coefficient rows match and every expression vanishes at the origin.
        expect_zero(f"row-image check for {name}", target_row - sum((c * r for c, r in zip(coeffs, coord_rows)), sp.zeros(1, len(primitive))))
        if sp.simplify(expr.subs({p: 0 for p in primitive})) != 0:
            raise AssertionError(f"Expected {name} to vanish at the primitive origin.")

    coeff_matrix = sp.Matrix(coeff_matrix_rows)

    return {
        "E_star": data["E_star"],
        "F_star": data["F_star"],
        "chi0": data["chi0"],
        "deltaU": data["deltaU"],
        "primitive": primitive,
        "primitive_null": primitive_null,
        "coords": coords,
        "coord_rank": coord_rank,
        "coeff_matrix": coeff_matrix,
        **param_exprs,
    }


if __name__ == "__main__":
    banner("STAGE 224 — EXACT OVERLAP-IMAGE PARAMETERIZATION")

    data = constructive_packet_null_data()
    param = overlap_image_parameterization()

    sigma_K, sigma_kappa, ell_B, ell_W, ell_R = data["primitive"]
    coords = param["coords"]
    dlnchi0, dlnZW, dlnC, dlnvarpi = coords

    subbanner("I. Constructive compensation solve")
    print("Constructive branch constants:")
    print("  E_* =", data["E_star"])
    print("  F_* =", data["F_star"])
    print("  chi_0 =", data["chi0"])
    print("  delta_U =", data["deltaU"])
    print()
    print("Unique compensation law in primitive overlap drifts:")
    ell_U_key = next(k for k in data["compensated"] if k.name == "ell_U")
    ell_OmegaU_key = next(k for k in data["compensated"] if k.name == "ell_OmegaU")
    ell_varpi_key = next(k for k in data["compensated"] if k.name == "ell_varpi")
    print("ell_U =")
    sp.pprint(data["compensated"][ell_U_key])
    print("ell_OmegaU =")
    sp.pprint(data["compensated"][ell_OmegaU_key])
    print("ell_varpi =")
    sp.pprint(data["compensated"][ell_varpi_key])

    subbanner("II. Exact overlap-observable drift set")
    for key in [
        "dlnK",
        "dlnM",
        "dlnC",
        "dlnvarpi",
        "dlnOmegaU",
        "dlnOmegaW",
        "dlnchi0",
        "dlnepsilon_eta",
        "dlnZW",
        "dlnepsilon_W",
        "dlndeltaU",
    ]:
        print(f"{key} =")
        sp.pprint(data["overlap_drifts"][key])
        print()

    subbanner("III. Image dimension and primitive null direction")
    print("rank of (dlnchi0, dlnZW, dlnC, dlnvarpi) as functions of the five primitive drifts =", param["coord_rank"])
    print("primitive null direction =")
    sp.pprint(param["primitive_null"])
    print("Interpretation: the primitive shift")
    print("  (sigma_kappa, ell_B, ell_W, ell_R) -> (-1, +1, +1, +1)")
    print("is invisible to the overlap-image coordinates.")

    subbanner("IV. Branch-adapted 4D overlap-image coordinates")
    print("Free overlap-image coordinates are")
    print("  (dlnchi0, dlnZW, dlnC, dlnvarpi).")
    print()
    print("Simple exact image laws:")
    print("dlnOmegaW =")
    sp.pprint(param["dlnOmegaW"])
    print("dlnepsilon_eta =")
    sp.pprint(param["dlnepsilon_eta"])
    print("dlnepsilon_W =")
    sp.pprint(param["dlnepsilon_W"])
    print("dlndeltaU =")
    sp.pprint(param["dlndeltaU"])
    print()
    print("The remaining image formulas are linear in the same four coordinates. Their exact")
    print("coefficient matrix is:")
    sp.pprint(param["coeff_matrix"])
    print()
    print("Numerically, the nontrivial rows are approximately:")
    for key in ("dlnK", "dlnM", "dlnOmegaU"):
        expr = param[key]
        coeffs = [sp.N(sp.diff(expr, q), 12) for q in coords]
        print(f"  {key} ≈ {coeffs[0]}*dlnchi0 + {coeffs[1]}*dlnZW + {coeffs[2]}*dlnC + {coeffs[3]}*dlnvarpi")

    banner("STAGE 224 LEDGER")
    print("1. The constructive support-restored packet-null family induces an exact 11-component")
    print("   overlap-observable drift vector.")
    print("2. Although it begins from five primitive overlap-side drifts")
    print("      (sigma_K, sigma_kappa, ell_B, ell_W, ell_R),")
    print("   its observable image has exact dimension 4.")
    print("3. A convenient exact coordinate system on that image is")
    print("      (dlnchi0, dlnZW, dlnC, dlnvarpi).")
    print("4. In these coordinates, the image splits naturally into")
    print("      a 2D orbit/shape layer (dlnchi0, dlnZW) and")
    print("      a 2D support layer      (dlnC, dlnvarpi).")
