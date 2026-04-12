#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 221 — realizability map from moving-throat overlap drifts to the Stage-220
packet-null family.

What this script does
---------------------
1. Reuses the exact constructive Packet-A/B null family from Stage 220.
2. Introduces the first actual overlap-state primitive drifts for the coherent
   moving-throat branch:
      - sigma_K      = d ln K_geo,
      - sigma_kappa  = d ln kappa,
      - ell_B        = d ln lambda_B,
      - ell_W        = d ln lambda_W,
      - ell_R        = d ln lambda_R,
   together with the still-free support / upper-leg carriers
      - ell_U        = d ln lambda_U,
      - ell_OU       = d ln Omega_U,
      - ell_varpi    = d ln varpi.
3. Uses the exact overlap dictionary
      C = lambda_B kappa,
      G_U = lambda_U,
      G_W = lambda_W kappa,
      R = lambda_R kappa,
   to map these primitive drifts into the Stage-220 carrier space.
4. Solves the full packet-null equations exactly for
      (ell_U, ell_OU, ell_varpi)
   in terms of the genuine overlap-side tangent data
      (sigma_K, sigma_kappa, ell_B, ell_W, ell_R).
5. Computes the induced orbit-lock companion drifts
      dln(delta_U), dln(M), dln(Omega_W).

Interpretation
--------------
This is the first exact realizability compiler from moving-throat overlap-state
branch data into the Stage-220 4-parameter packet-null family. It turns the next
question into a concrete one:

  given the actual overlap-state tangent of the moving-throat branch,
  does it satisfy these exact linear compensation laws?
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


if __name__ == "__main__":
    banner("STAGE 221 — OVERLAP-STATE REALIZABILITY MAP")

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

    # Primitive overlap-state branch drifts.
    sigma_K, sigma_kappa = sp.symbols("sigma_K sigma_kappa", real=True)
    ell_B, ell_U, ell_W, ell_R = sp.symbols("ell_B ell_U ell_W ell_R", real=True)
    ell_OU, ell_varpi = sp.symbols("ell_OmegaU ell_varpi", real=True)

    subbanner("I. Exact overlap-state dictionary")
    print("C   = lambda_B * kappa")
    print("G_U = lambda_U")
    print("G_W = lambda_W * kappa")
    print("R   = lambda_R * kappa")
    print()
    print("Carrier map:")
    print("  alpha_K   = sigma_K")
    print("  alpha_GW  = ell_W + sigma_kappa")
    print("  alpha_GU  = ell_U")
    print("  alpha_R   = ell_R + sigma_kappa")
    print("  alpha_OU  = ell_OmegaU")
    print("  beta_B    = ell_B + sigma_kappa")
    print("  beta_varpi= ell_varpi")

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

    subbanner("II. Exact packet equations in overlap-state drifts")
    print("u2^(1) =")
    sp.pprint(u2_overlap)
    print("u4^(1) =")
    sp.pprint(u4_overlap)
    print("Xi_load =")
    sp.pprint(Xi_overlap)

    subbanner("III. Exact realization solve")
    sol_list = sp.solve(
        [sp.Eq(u2_overlap, 0), sp.Eq(u4_overlap, 0), sp.Eq(Xi_overlap, 0)],
        [ell_U, ell_OU, ell_varpi],
        dict=True,
    )
    if len(sol_list) != 1:
        raise AssertionError("Expected a unique realization solve for (ell_U, ell_OmegaU, ell_varpi).")
    sol = {k: sp.simplify(v) for k, v in sol_list[0].items()}

    print("ell_U =")
    sp.pprint(sol[ell_U])
    print("ell_OmegaU =")
    sp.pprint(sol[ell_OU])
    print("ell_varpi =")
    sp.pprint(sol[ell_varpi])
    print()
    print("Numerically:")
    print("ell_U ≈")
    sp.pprint(sp.N(sol[ell_U], 14))
    print("ell_OmegaU ≈")
    sp.pprint(sp.N(sol[ell_OU], 14))
    print("ell_varpi ≈")
    sp.pprint(sp.N(sol[ell_varpi], 14))

    expect_zero("u2^(1) after realization solve", sp.simplify(u2_overlap.subs(sol)))
    expect_zero("u4^(1) after realization solve", sp.simplify(u4_overlap.subs(sol)))
    expect_zero("Xi_load after realization solve", sp.simplify(Xi_overlap.subs(sol)))

    subbanner("IV. Induced orbit-lock companion drifts")
    alpha_GW_eff = sp.simplify(ell_W + sigma_kappa)
    alpha_GU_eff = sp.simplify(sol[ell_U])
    alpha_R_eff = sp.simplify(ell_R + sigma_kappa)
    alpha_OU_eff = sp.simplify(sol[ell_OU])
    alpha_K_eff = sigma_K

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

    print("dln(delta_U) =")
    sp.pprint(dln_deltaU)
    print("dln(M) =")
    sp.pprint(dln_M)
    print("dln(Omega_W) =")
    sp.pprint(dln_OW)
    print()
    print("Numerically:")
    print("dln(delta_U) ≈")
    sp.pprint(sp.N(dln_deltaU, 14))
    print("dln(M) ≈")
    sp.pprint(sp.N(dln_M, 14))
    print("dln(Omega_W) ≈")
    sp.pprint(sp.N(dln_OW, 14))

    banner("STAGE 221 LEDGER")
    print("1. The exact packet-null family can now be written directly in moving-throat overlap-state drifts.")
    print("2. The genuinely branch-level input data are")
    print("      (sigma_K, sigma_kappa, ell_B, ell_W, ell_R).")
    print("3. Once those are supplied, exact packet-null realizability fixes uniquely")
    print("      ell_U, ell_OmegaU, ell_varpi,")
    print("   and therefore also the orbit-lock companion drifts")
    print("      dln(delta_U), dln(M), dln(Omega_W).")
    print("4. So the realizability problem is no longer 'find some null direction somehow'.")
    print("   It is now the exact test of whether the actual overlap-state tangent of the")
    print("   moving-throat branch satisfies these linear compensation laws.")
