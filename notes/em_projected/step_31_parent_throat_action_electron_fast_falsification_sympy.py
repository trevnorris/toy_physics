#!/usr/bin/env python3
"""Fast electron falsification screen for the moderate branch-B patch."""
from __future__ import annotations

import hashlib
import json
import math

from step_29_parent_throat_action_moderate_branchb_sector_sympy import export_step29_moderate_branchb_patch


def assert_close(label: str, actual: float, expected: float, tol: float) -> None:
    if abs(actual - expected) > tol:
        raise AssertionError(f"{label} failed: {actual} vs {expected} (tol={tol})")


def assert_nonzero(label: str, value: float, floor: float) -> None:
    if abs(value) <= floor:
        raise AssertionError(f"{label} failed: |{value}| <= {floor}")


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_electron_fast_falsification",
        "pre_target_freeze": True,
        "target_blind": False,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "conditional electron screen with an explicit SI-to-PDE port scale S_port",
    }

    # Moderate branch-B patch data from step29.
    step29_patch = export_step29_moderate_branchb_patch()
    patch_20 = step29_patch["same_scale_lambda20"]
    patch_50 = step29_patch["same_scale_lambda50"]
    P0_base = float(step29_patch["P0_base"])
    lambda_out_20 = float(patch_20["lambda_out"])
    lambda_out_50 = float(patch_50["lambda_out"])
    P0_20 = lambda_out_20 * P0_base
    P0_50 = lambda_out_50 * P0_base
    mhat0_20 = float(patch_20["mhat0_req"])
    mhat0_50 = float(patch_50["mhat0_req"])
    K_patch_20 = mhat0_20**2 * P0_20
    K_patch_50 = mhat0_50**2 * P0_50
    K_unscaled_50 = mhat0_50**2 * P0_base

    # CODATA-like electron constants.
    G = 6.67430e-11
    c = 299792458.0
    hbar = 1.054571817e-34
    m_e = 9.1093837015e-31
    r_e = 2.8179403262e-15

    omega_compton = m_e * c**2 / hbar
    S_port_direct = 1.0

    # Port equation:
    #     S_port * (mhat0^2 lambda_out P0_base)
    #       = 54 G c_s^5 / (5 a^5 c^5), with a = r_e.
    # The naive direct-SI closure is the special case S_port = 1.
    c_s_direct = ((5.0 * S_port_direct * K_patch_50 * r_e**5 * c**5) / (54.0 * G)) ** (1.0 / 5.0)
    Omega_Q_direct = 3.0 * c_s_direct / (2.0 * r_e)
    ratio_direct = Omega_Q_direct / omega_compton
    E_Q_direct = hbar * Omega_Q_direct
    E_ratio = E_Q_direct / (m_e * c**2)
    K_direct_from_constants = 54.0 * G * c_s_direct**5 / (5.0 * r_e**5 * c**5)

    # Reverse calibration: force the pole scale to the electron Compton scale.
    c_s_compton = 2.0 * r_e * omega_compton / 3.0
    K_required_compton = 54.0 * G * c_s_compton**5 / (5.0 * r_e**5 * c**5)
    S_port_required_compton = K_required_compton / K_patch_50
    S_port_required_compton_20 = K_required_compton / K_patch_20
    S_port_required_if_lambda50_omitted = K_required_compton / K_unscaled_50
    lambda_load_guard = S_port_required_if_lambda50_omitted / S_port_required_compton
    direct_S_port_mutation_residual = K_direct_from_constants - 2.0 * S_port_direct * K_patch_50
    compton_S_port_mutation_residual = K_required_compton - 1.01 * S_port_required_compton * K_patch_50

    branch_freeze_payload = {
        "metadata": branch_metadata,
        "patch_data": {
            "P0_base": P0_base,
            "lambda_out_20": lambda_out_20,
            "lambda_out_50": lambda_out_50,
            "mhat0_20": mhat0_20,
            "mhat0_50": mhat0_50,
        },
        "si_constants": {"G": G, "c": c, "hbar": hbar, "m_e": m_e, "r_e": r_e},
        "exports": {
            "K_patch_20": K_patch_20,
            "K_patch_50": K_patch_50,
            "c_s_direct": c_s_direct,
            "S_port_required_compton": S_port_required_compton,
        },
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_freeze_payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    assert_close("step29 lambda20 patch K", K_patch_20, 10.8, 1e-12)
    assert_close("step29 lambda50 patch K", K_patch_50, 10.8, 1e-12)
    assert_close("direct-SI c_s", c_s_direct, 9.159491211330665e-05, 1e-16)
    assert_close("direct-SI Omega_Q", Omega_Q_direct, 48756308603.324455, 1e-6)
    assert_close("Compton omega", omega_compton, 7.76344071105011e20, 1e8)
    assert_close("direct Omega/Compton ratio", ratio_direct, 6.280244857660478e-11, 1e-22)
    assert_close("direct energy ratio", E_ratio, 6.280244857660478e-11, 1e-22)
    assert_close("Compton-calibrated c_s", c_s_compton, 1458460.8433153937, 1e-6)
    assert_close("Compton-calibrated K", K_required_compton, 1.1054545016851398e52, 1e36)
    assert_close("Compton-calibrated S_port", S_port_required_compton, 1.023568983041796e51, 2e35)
    assert_close("lambda20/lambda50 S_port agreement", S_port_required_compton_20 / S_port_required_compton, 1.0, 1e-12)
    assert_close("omitting lambda50 changes required S_port by lambda50", lambda_load_guard, 50.0, 1e-12)
    assert_close("direct constants round-trip to patch K", K_direct_from_constants, K_patch_50, 1e-12)
    assert_nonzero("mutating direct S_port breaks the port equation", direct_S_port_mutation_residual, 1.0)
    assert_nonzero("mutating Compton S_port breaks the port equation", compton_S_port_mutation_residual, 1e49)

    print("STEP 31 ELECTRON FAST FALSIFICATION AUDIT")
    print("Tested the shortest conditional electron screen for the moderate branch-B patch using real constants.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Moderate branch-B patch data:")
    print("  P0_base =", P0_base)
    print("  lambda_out = 20 -> P0 =", P0_20, "mhat0 =", mhat0_20, "mhat0^2 P0 =", K_patch_20)
    print("  lambda_out = 50 -> P0 =", P0_50, "mhat0 =", mhat0_50, "mhat0^2 P0 =", K_patch_50)
    print("Electron constants:")
    print("  G =", G)
    print("  c =", c)
    print("  hbar =", hbar)
    print("  m_e =", m_e)
    print("  r_e =", r_e)
    print("  omega_Compton =", omega_compton)
    print("Dimensional port equation:")
    print("  S_port * (mhat0^2 lambda_out P0_base) = 54 G c_s^5 / (5 r_e^5 c^5)")
    print("  direct closure uses S_port =", S_port_direct)
    print("Direct-SI closure with a = r_e and step29 patch K:")
    print("  c_s =", c_s_direct)
    print("  Omega_Q =", Omega_Q_direct)
    print("  Omega_Q / omega_Compton =", ratio_direct)
    print("  hbar*Omega_Q / (m_e c^2) =", E_ratio)
    print("Reverse calibration forcing Omega_Q = omega_Compton:")
    print("  c_s =", c_s_compton)
    print("  required K =", K_required_compton)
    print("  required S_port =", S_port_required_compton)
    print("  required S_port if lambda_out=50 is omitted =", S_port_required_if_lambda50_omitted)
    print("  lambda load guard =", lambda_load_guard)
    print("  direct S_port mutation residual =", direct_S_port_mutation_residual)
    print("  Compton S_port mutation residual =", compton_S_port_mutation_residual)
    print("Interpretation:")
    print("  Under the naive direct-SI closure, the moderate branch-B patch predicts a pole scale about 1.59e10 times too small for an electron.")
    print("  If one instead forces the pole to the electron Compton scale, the required port scale is S_port ~= 1.02e51.")
    print("  So the naive direct electron identification is falsified quickly: it fails either on the pole scale or on the normalization scale.")
    print("  This is a conditional falsification of the direct-SI electron map, not of the whole analog program with an unresolved dimensionalization map.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
