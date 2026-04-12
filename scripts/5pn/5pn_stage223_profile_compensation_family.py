#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 223 — compensated coherent-profile realization family.

What this script does
---------------------
1. Starts from the Stage-221 exact overlap-state realization map.
2. Specializes it to the first coherent moving-throat overlap family
      kappa(theta), K_geo(theta),
   with lambda_B, lambda_W, lambda_R held fixed.
3. Solves for the unique compensating drifts
      dln lambda_U, dln Omega_U, dln varpi
   required to restore packet-null tangency.
4. Computes the induced orbit-lock companion drifts
      dln delta_U, dln M, dln Omega_W.
5. Evaluates the compensated branch at two concrete profile points:
      - theta = 0,
      - theta = theta_max.

Interpretation
--------------
This is the first explicit compensated realizability family attached directly to
moving-throat overlap data. The pure profile branch is dead, but it does admit a
unique packet-null completion once support / upper-leg / frequency co-drifts are
added in the exact proportions derived here.
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
    banner("STAGE 223 — COMPENSATED COHERENT-PROFILE REALIZATION FAMILY")

    build_master_system = load_stage218_builder()
    E_star = sp.Rational(1, 4)
    F_star = sp.Rational(5, 6)
    data = build_master_system(E_star, F_star)

    alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi = data["free"]
    chi0 = sp.simplify(data["chi0"])
    deltaU = sp.Integer(1)
    u2_1 = data["u2_1"]
    u4_1 = data["u4_1"]
    Xi_load = data["Xi_load"]

    # Primitive overlap-state drifts.
    sigma_K, sigma_kappa = sp.symbols("sigma_K sigma_kappa", real=True)
    ell_U, ell_OU, ell_varpi = sp.symbols("ell_U ell_OmegaU ell_varpi", real=True)

    # Fixed-coupling coherent-profile specialization:
    #   ell_B = ell_W = ell_R = 0,
    #   beta_B = sigma_kappa,
    #   alpha_GW = sigma_kappa,
    #   alpha_R = sigma_kappa.
    carrier_subs = {
        alpha_K: sigma_K,
        alpha_GW: sigma_kappa,
        alpha_GU: ell_U,
        alpha_R: sigma_kappa,
        alpha_OU: ell_OU,
        beta_B: sigma_kappa,
        beta_varpi: ell_varpi,
    }

    u2_comp = sp.simplify(u2_1.subs(carrier_subs))
    u4_comp = sp.simplify(u4_1.subs(carrier_subs))
    Xi_comp = sp.simplify(Xi_load.subs(carrier_subs))

    subbanner("I. Unique compensation law for the fixed-coupling coherent profile")
    sol_list = sp.solve(
        [sp.Eq(u2_comp, 0), sp.Eq(u4_comp, 0), sp.Eq(Xi_comp, 0)],
        [ell_U, ell_OU, ell_varpi],
        dict=True,
    )
    if len(sol_list) != 1:
        raise AssertionError("Expected a unique compensation law for (ell_U, ell_OmegaU, ell_varpi).")
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

    expect_zero("u2^(1) after compensation", sp.simplify(u2_comp.subs(sol)))
    expect_zero("u4^(1) after compensation", sp.simplify(u4_comp.subs(sol)))
    expect_zero("Xi_load after compensation", sp.simplify(Xi_comp.subs(sol)))

    subbanner("II. Induced orbit-lock companion drifts on the compensated branch")
    alpha_deltaU = sp.simplify(
        -(1 + deltaU) / (1 + chi0) * ((sigma_kappa) + sol[ell_U] - (sigma_kappa) - 2 * sol[ell_OU])
    )
    alpha_M = sp.simplify(sigma_K - 2 * sol[ell_U] + 2 * sol[ell_OU])
    alpha_OW = sp.simplify(
        (
            sigma_kappa
            - sol[ell_U]
            + (1 - E_star) * sol[ell_OU]
            + E_star * sigma_kappa
            - (F_star / 2) * alpha_deltaU
        )
        / (E_star + 2)
    )

    print("dln(delta_U) =")
    sp.pprint(alpha_deltaU)
    print("dln(M) =")
    sp.pprint(alpha_M)
    print("dln(Omega_W) =")
    sp.pprint(alpha_OW)
    print()
    print("Numerically:")
    print("dln(delta_U) ≈")
    sp.pprint(sp.N(alpha_deltaU, 14))
    print("dln(M) ≈")
    sp.pprint(sp.N(alpha_M, 14))
    print("dln(Omega_W) ≈")
    sp.pprint(sp.N(alpha_OW, 14))

    subbanner("III. First coherent overlap family kappa(theta), K_geo(theta)")
    theta, dtheta = sp.symbols("theta dtheta", real=True)
    K0, T_grad = sp.symbols("K_0 T_grad", positive=True, real=True)

    kappa0 = sp.simplify(2 * sp.sqrt(2) / sp.pi)
    kappa1 = sp.simplify(-4 / (3 * sp.pi))
    kappa_theta = sp.simplify(kappa0 * sp.cos(theta) + kappa1 * sp.sin(theta))
    K_geo_theta = sp.simplify(K0 + T_grad * sp.sin(theta) ** 2)

    sigma_kappa_theta = sp.simplify(sp.diff(sp.log(kappa_theta), theta) * dtheta)
    sigma_K_theta = sp.simplify(sp.diff(sp.log(K_geo_theta), theta) * dtheta)

    print("sigma_kappa(theta) =")
    sp.pprint(sigma_kappa_theta)
    print("sigma_K(theta) =")
    sp.pprint(sigma_K_theta)

    comp_theta = {
        sigma_kappa: sigma_kappa_theta,
        sigma_K: sigma_K_theta,
    }

    print("Compensating drifts along the coherent profile family:")
    print("dln(lambda_U) =")
    sp.pprint(sp.simplify(sol[ell_U].subs(comp_theta)))
    print("dln(Omega_U) =")
    sp.pprint(sp.simplify(sol[ell_OU].subs(comp_theta)))
    print("dln(varpi) =")
    sp.pprint(sp.simplify(sol[ell_varpi].subs(comp_theta)))

    subbanner("IV. Concrete profile points")
    theta_max = sp.simplify(-sp.atan(sp.sqrt(2) / 3))

    eval_points = {
        "theta = 0": {theta: 0},
        "theta = theta_max": {theta: theta_max},
    }

    for label, point in eval_points.items():
        print(label)
        print("  sigma_kappa / dtheta =")
        sp.pprint(sp.simplify((sigma_kappa_theta / dtheta).subs(point)))
        print("  sigma_K / dtheta =")
        sp.pprint(sp.simplify((sigma_K_theta / dtheta).subs(point)))
        print("  dln(lambda_U) / dtheta ≈")
        sp.pprint(sp.N((sol[ell_U].subs(comp_theta) / dtheta).subs(point), 12))
        print("  dln(Omega_U) / dtheta ≈")
        sp.pprint(sp.N((sol[ell_OU].subs(comp_theta) / dtheta).subs(point), 12))
        print("  dln(varpi) / dtheta ≈")
        sp.pprint(sp.N((sol[ell_varpi].subs(comp_theta) / dtheta).subs(point), 12))
        print("  dln(delta_U) / dtheta ≈")
        sp.pprint(sp.N((alpha_deltaU.subs(comp_theta) / dtheta).subs(point), 12))
        print("  dln(M) / dtheta ≈")
        sp.pprint(sp.N((alpha_M.subs(comp_theta) / dtheta).subs(point), 12))
        print("  dln(Omega_W) / dtheta ≈")
        sp.pprint(sp.N((alpha_OW.subs(comp_theta) / dtheta).subs(point), 12))
        print()

    banner("STAGE 223 LEDGER")
    print("1. Once the coherent-profile branch is given by")
    print("      sigma_kappa(theta) = dln kappa(theta),   sigma_K(theta) = dln K_geo(theta),")
    print("   there is a unique packet-null completion with fixed couplings")
    print("      dln(lambda_U), dln(Omega_U), dln(varpi).")
    print("2. The induced orbit-lock companion drifts")
    print("      dln(delta_U), dln(M), dln(Omega_W)")
    print("   are then fixed automatically." )
    print("3. So the pure profile branch is dead, but it is not a dead end: it admits a unique")
    print("   compensated realization family inside the exact Stage-220 packet-null manifold.")
    print("4. The next actual PDE question is whether the moving-throat branch dynamically generates")
    print("   these compensation ratios, not whether such a compensated family exists algebraically.")
