#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 222 — pure coherent-profile branch no-go for first-order packet-null realization.

What this script does
---------------------
1. Reuses the exact constructive packet equations from Stage 218.
2. Inserts the first concrete moving-throat coherent nonconstant overlap family,
   but with all microscopic couplings/frequencies held fixed:
      - only kappa(theta) and K_geo(theta) are allowed to drift.
3. Proves that the induced branch tangent cannot lie in the Stage-220 packet-null
   family except at the trivial zero tangent.
4. Evaluates the two most important concrete profile points:
      - theta = 0  (flat/constant branch),
      - theta = theta_max  (max-coupling branch).

Interpretation
--------------
This is the first serious realizability result. A pure geometry/profile deformation
of the explicit moving-throat overlap family is not enough. Nontrivial packet-null
realization requires compensating support / upper-leg / frequency co-drifts.
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
    banner("STAGE 222 — PURE COHERENT-PROFILE BRANCH NO-GO")

    build_master_system = load_stage218_builder()
    data = build_master_system(sp.Rational(1, 4), sp.Rational(5, 6))
    alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi = data["free"]
    u2_1 = data["u2_1"]
    u4_1 = data["u4_1"]
    Xi_load = data["Xi_load"]

    sigma_K, sigma_kappa = sp.symbols("sigma_K sigma_kappa", real=True)

    # Pure profile drift: only K_geo and kappa move; all microscopic couplings and
    # frequencies are frozen.
    pure_profile_subs = {
        alpha_K: sigma_K,
        alpha_GW: sigma_kappa,
        alpha_GU: 0,
        alpha_R: sigma_kappa,
        alpha_OU: 0,
        beta_B: sigma_kappa,
        beta_varpi: 0,
    }

    u2_profile = sp.simplify(u2_1.subs(pure_profile_subs))
    u4_profile = sp.simplify(u4_1.subs(pure_profile_subs))
    Xi_profile = sp.simplify(Xi_load.subs(pure_profile_subs))

    subbanner("I. Exact packet defects on the pure profile branch")
    print("u2^(1) =")
    sp.pprint(u2_profile)
    print("u4^(1) =")
    sp.pprint(u4_profile)
    print("Xi_load =")
    sp.pprint(Xi_profile)

    subbanner("II. Exact no-go solve")
    sol = sp.solve(
        [sp.Eq(u2_profile, 0), sp.Eq(u4_profile, 0), sp.Eq(Xi_profile, 0)],
        [sigma_K, sigma_kappa],
        dict=True,
    )
    print("Solutions of the pure-profile packet system =")
    sp.pprint(sol)
    if sol != [{sigma_K: 0, sigma_kappa: 0}]:
        raise AssertionError("Expected only the trivial zero tangent on the pure profile branch.")

    subbanner("III. Concrete coherent family kappa(theta) and K_geo(theta)")
    theta, dtheta = sp.symbols("theta dtheta", real=True)
    K0, T_grad = sp.symbols("K_0 T_grad", positive=True, real=True)

    kappa0 = sp.simplify(2 * sp.sqrt(2) / sp.pi)
    kappa1 = sp.simplify(-4 / (3 * sp.pi))
    kappa_theta = sp.simplify(kappa0 * sp.cos(theta) + kappa1 * sp.sin(theta))
    K_geo_theta = sp.simplify(K0 + T_grad * sp.sin(theta) ** 2)

    sigma_kappa_theta = sp.simplify(sp.diff(sp.log(kappa_theta), theta) * dtheta)
    sigma_K_theta = sp.simplify(sp.diff(sp.log(K_geo_theta), theta) * dtheta)

    print("kappa(theta) =")
    sp.pprint(kappa_theta)
    print("sigma_kappa(theta) = dln kappa =")
    sp.pprint(sigma_kappa_theta)
    print("K_geo(theta) =")
    sp.pprint(K_geo_theta)
    print("sigma_K(theta) = dln K_geo =")
    sp.pprint(sigma_K_theta)

    # Two important branch points.
    theta_max = sp.simplify(-sp.atan(sp.sqrt(2) / 3))
    flat_point = {
        theta: 0,
    }
    max_point = {
        theta: theta_max,
    }

    subbanner("IV. Flat branch and max-coupling branch diagnostics")
    print("At theta = 0:")
    print("  sigma_kappa / dtheta =")
    sp.pprint(sp.simplify((sigma_kappa_theta / dtheta).subs(flat_point)))
    print("  sigma_K / dtheta =")
    sp.pprint(sp.simplify((sigma_K_theta / dtheta).subs(flat_point)))

    print("At theta = theta_max = -arctan(sqrt(2)/3):")
    print("  sigma_kappa / dtheta =")
    sp.pprint(sp.simplify((sigma_kappa_theta / dtheta).subs(max_point)))
    print("  sigma_K / dtheta =")
    sp.pprint(sp.simplify((sigma_K_theta / dtheta).subs(max_point)))

    # Show explicitly that the two nontrivial branches miss the zero set.
    flat_res = {
        sigma_kappa: sp.simplify((sigma_kappa_theta / dtheta).subs(flat_point)),
        sigma_K: sp.simplify((sigma_K_theta / dtheta).subs(flat_point)),
    }
    max_res = {
        sigma_kappa: sp.simplify((sigma_kappa_theta / dtheta).subs(max_point)),
        sigma_K: sp.simplify((sigma_K_theta / dtheta).subs(max_point)),
    }

    print("Packet defects at theta = 0 (normalized by dtheta) =")
    sp.pprint(
        sp.Matrix(
            [
                sp.simplify((u2_profile / dtheta).subs(flat_res)),
                sp.simplify((u4_profile / dtheta).subs(flat_res)),
                sp.simplify((Xi_profile / dtheta).subs(flat_res)),
            ]
        )
    )

    print("Packet defects at theta = theta_max (normalized by dtheta) =")
    sp.pprint(
        sp.Matrix(
            [
                sp.simplify((u2_profile / dtheta).subs(max_res)),
                sp.simplify((u4_profile / dtheta).subs(max_res)),
                sp.simplify((Xi_profile / dtheta).subs(max_res)),
            ]
        )
    )

    banner("STAGE 222 LEDGER")
    print("1. A pure coherent-profile deformation with frozen microscopic couplings/frequencies")
    print("   induces only two logarithmic drifts:")
    print("      sigma_kappa = dln kappa,   sigma_K = dln K_geo.")
    print("2. The exact Stage-220 packet equations then have only the trivial solution")
    print("      sigma_kappa = sigma_K = 0.")
    print("3. Therefore the moving-throat coherent profile family cannot realize a nontrivial")
    print("   packet-null tangent by geometry/profile motion alone.")
    print("4. On the flat branch theta = 0, one has sigma_K = 0 but sigma_kappa = -sqrt(2)/3 dtheta != 0.")
    print("5. On the max-coupling branch theta = theta_max, one has sigma_kappa = 0 but sigma_K != 0")
    print("   whenever the axial-gradient penalty T_grad is present.")
    print("6. So the next realizability gate is necessarily a compensated one: support / upper-leg /")
    print("   frequency co-drifts must accompany the overlap-profile motion.")
