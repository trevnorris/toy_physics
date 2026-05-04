#!/usr/bin/env python3
"""Self-test the solver snapshot adapters into the runtime-monitor schema."""
from __future__ import annotations

import hashlib
import json
import pathlib
import tempfile

import numpy as np

from cfd_runtime_monitor_postprocess import analyze_snapshot, divergence3, gradient3, laplacian3, solve_phi_from_velocity
from cfd_snapshot_adapters import adapt_monopole_3d, adapt_wavefunction_4d, load_npz, save_npz
from step_34_parent_throat_action_cfd_runtime_postprocessor_sympy import make_periodic_consistency_snapshot


def assert_small(label: str, value: float, tol: float) -> None:
    if not np.isfinite(value) or abs(value) > tol:
        raise AssertionError(f"{label} failed: {value} > {tol}")


def assert_rel_small(label: str, value: np.ndarray, target: np.ndarray, tol: float) -> float:
    scale = max(float(np.max(np.abs(target))), 1e-12)
    rel = float(np.max(np.abs(value - target)) / scale)
    if not np.isfinite(rel) or rel > tol:
        raise AssertionError(f"{label} failed: {rel} > {tol}")
    return rel


def make_wavefunction_input() -> dict[str, np.ndarray]:
    base = make_periodic_consistency_snapshot()
    phi = base["phi3"]
    rho0 = float(base["rho0"])
    nw = base["w"].size
    psi = np.sqrt(rho0) * np.exp(1j * phi[..., None]) * np.ones((1, 1, 1, nw), dtype=np.complex128)
    return {
        "psi": psi,
        "x": base["x"],
        "y": base["y"],
        "z": base["z"],
        "w": base["w"],
        "W": base["W"],
        "rho0": base["rho0"],
        "phi3": base["phi3"],
        "N_probe": base["N_probe"],
        "Phi_eff": base["Phi_eff"],
        "hbar": np.array(1.0),
        "m": np.array(1.0),
    }


def main() -> None:
    branch_metadata = {
        "branch_id": "v2_local_parent_background_snapshot_adapter",
        "pre_target_freeze": True,
        "target_blind": False,
        "no_post_residual_refit": True,
        "boundary_class": "open_impedance_demo",
        "interpretation": "self-test for the adapter layer from solver-side dumps into the runtime-monitor snapshot schemas",
    }

    periodic = make_periodic_consistency_snapshot()
    spacings = (
        float(periodic["x"][1] - periodic["x"][0]),
        float(periodic["y"][1] - periodic["y"][0]),
        float(periodic["z"][1] - periodic["z"][0]),
    )

    wave_payload = adapt_wavefunction_4d(make_wavefunction_input(), mass=1.0, hbar=1.0, periodic_xyz=True, periodic_w=False)
    projected_jx = np.trapz(wave_payload["W"].reshape(1, 1, 1, -1) * wave_payload["jx"], x=wave_payload["w"], axis=-1)
    projected_jy = np.trapz(wave_payload["W"].reshape(1, 1, 1, -1) * wave_payload["jy"], x=wave_payload["w"], axis=-1)
    projected_jz = np.trapz(wave_payload["W"].reshape(1, 1, 1, -1) * wave_payload["jz"], x=wave_payload["w"], axis=-1)
    rho_brane_wave = np.trapz(wave_payload["W"].reshape(1, 1, 1, -1) * wave_payload["rho"], x=wave_payload["w"], axis=-1)
    expected_grad = gradient3(periodic["phi3"], spacings, periodic=True)
    wave_jx_rel = assert_rel_small("wave adapter jx", projected_jx, periodic["rho0"] * expected_grad[0], 5e-4)
    wave_jy_rel = assert_rel_small("wave adapter jy", projected_jy, periodic["rho0"] * expected_grad[1], 5e-4)
    wave_jz_rel = assert_rel_small("wave adapter jz", projected_jz, periodic["rho0"] * expected_grad[2], 5e-4)
    vx = projected_jx / np.maximum(np.abs(rho_brane_wave), 1e-12)
    vy = projected_jy / np.maximum(np.abs(rho_brane_wave), 1e-12)
    vz = projected_jz / np.maximum(np.abs(rho_brane_wave), 1e-12)
    phi_recovered = solve_phi_from_velocity(vx, vy, vz, spacings)
    wave_payload["dt_rho"] = -divergence3(projected_jx, projected_jy, projected_jz, spacings, periodic=True)
    wave_payload["phi3"] = phi_recovered
    wave_payload["Phi_eff"] = phi_recovered
    wave_payload["N_probe"] = 1.0 - 2.0 * phi_recovered
    wave_summary = analyze_snapshot(wave_payload, c_probe=1.0, bins=18, tail_fraction=0.35, periodic_xyz=True)
    assert_small("wave adapter continuity residual", wave_summary["max_abs_R_cont"], 1e-11)

    exact_proj_jx = np.trapz(periodic["W"].reshape(1, 1, 1, -1) * periodic["jx"], x=periodic["w"], axis=-1)
    exact_proj_jy = np.trapz(periodic["W"].reshape(1, 1, 1, -1) * periodic["jy"], x=periodic["w"], axis=-1)
    exact_proj_jz = np.trapz(periodic["W"].reshape(1, 1, 1, -1) * periodic["jz"], x=periodic["w"], axis=-1)
    exact_source = float(periodic["rho0"]) * laplacian3(periodic["phi3"], spacings, periodic=True)
    mono_input = {
        "rho": np.trapz(periodic["W"].reshape(1, 1, 1, -1) * periodic["rho"], x=periodic["w"], axis=-1),
        "mx": exact_proj_jx,
        "my": exact_proj_jy,
        "mz": exact_proj_jz,
        "S_rho": exact_source,
        "x": periodic["x"],
        "y": periodic["y"],
        "z": periodic["z"],
        "dt_rho": np.zeros_like(exact_proj_jx),
        "phi3": periodic["phi3"],
        "rho0": periodic["rho0"],
        "N_probe": periodic["N_probe"],
        "Phi_eff": periodic["Phi_eff"],
    }
    mono_payload = adapt_monopole_3d(mono_input)
    mono_summary = analyze_snapshot(mono_payload, c_probe=1.0, bins=18, tail_fraction=0.35, periodic_xyz=True)
    assert_small("monopole adapter continuity residual", mono_summary["max_abs_R_cont"], 1e-11)
    assert_small("monopole adapter exact Poisson residual", mono_summary["max_abs_R_pois_exact"], 1e-11)

    wave_missing_W = make_wavefunction_input()
    wave_missing_W.pop("W")
    try:
        adapt_wavefunction_4d(wave_missing_W, mass=1.0, hbar=1.0, periodic_xyz=True, periodic_w=False)
        raise AssertionError("wavefunction adapter accepted missing W without explicit opt-in")
    except KeyError as exc:
        if "--allow-uniform-W" not in str(exc):
            raise
    wave_missing_W_nonunit = dict(wave_missing_W)
    wave_missing_W_nonunit["w"] = np.linspace(-1.0, 1.0, wave_missing_W["w"].size)
    wave_uniform_payload = adapt_wavefunction_4d(
        wave_missing_W_nonunit,
        mass=1.0,
        hbar=1.0,
        periodic_xyz=True,
        periodic_w=False,
        allow_uniform_weight=True,
    )
    w_axis = wave_uniform_payload["w"]
    uniform_W_integral_error = abs(float(np.trapz(wave_uniform_payload["W"], x=w_axis)) - 1.0)
    assert_small("uniform W opt-in normalization", uniform_W_integral_error, 1e-15)
    uniform_W_expected_value = 1.0 / float(w_axis[-1] - w_axis[0])
    uniform_W_value_error = abs(float(wave_uniform_payload["W"][0]) - uniform_W_expected_value)
    assert_small("uniform W opt-in nonunit-span value", uniform_W_value_error, 1e-15)

    dx, dy, dz = spacings
    V_domain = float(exact_source.size) * dx * dy * dz
    W3 = np.full_like(exact_source, 1.0 / V_domain)
    Mdot = 0.75
    lambda_bulk = 0.4
    reconstructed_expected = -Mdot * W3 + (lambda_bulk * Mdot / V_domain) * np.ones_like(W3)
    mono_reconstruct_input = {
        "rho": mono_input["rho"],
        "mx": mono_input["mx"],
        "my": mono_input["my"],
        "mz": mono_input["mz"],
        "W": W3,
        "Mdot": np.array(Mdot),
        "lambda_bulk": np.array(lambda_bulk),
        "V_domain": np.array(V_domain),
        "x": periodic["x"],
        "y": periodic["y"],
        "z": periodic["z"],
    }
    mono_reconstructed = adapt_monopole_3d(mono_reconstruct_input)
    reconstruct_rel = assert_rel_small(
        "monopole reconstructed source",
        mono_reconstructed["S_rho"],
        reconstructed_expected,
        1e-15,
    )
    missing_lambda = dict(mono_reconstruct_input)
    missing_lambda.pop("lambda_bulk")
    try:
        adapt_monopole_3d(missing_lambda)
        raise AssertionError("monopole adapter accepted reconstructed source without lambda")
    except KeyError as exc:
        if "lambda" not in str(exc):
            raise

    with tempfile.TemporaryDirectory(prefix="em_projected_adapter_") as tmp:
        tmpdir = pathlib.Path(tmp)
        raw_wave = tmpdir / "raw_wave.npz"
        adapted_wave = tmpdir / "adapted_wave.npz"
        save_npz(raw_wave, make_wavefunction_input())
        save_npz(adapted_wave, wave_payload)
        roundtrip = load_npz(adapted_wave)
        if "rho" not in roundtrip or "jx" not in roundtrip:
            raise AssertionError("adapter roundtrip did not preserve full 4D schema")

    branch_freeze_payload = {
        "metadata": branch_metadata,
        "wave_metrics": {
            "rel_jx": wave_jx_rel,
            "rel_jy": wave_jy_rel,
            "rel_jz": wave_jz_rel,
            "max_abs_R_cont": wave_summary["max_abs_R_cont"],
            "alpha_fit_tail_mean": wave_summary["alpha_fit_tail_mean"],
            "uniform_W_integral_error": uniform_W_integral_error,
            "uniform_W_value_error": uniform_W_value_error,
        },
        "monopole_metrics": {
            "max_abs_R_cont": mono_summary["max_abs_R_cont"],
            "max_abs_R_pois_exact": mono_summary["max_abs_R_pois_exact"],
            "alpha_fit_tail_mean": mono_summary["alpha_fit_tail_mean"],
            "reconstructed_source_rel_error": reconstruct_rel,
        },
    }
    branch_freeze_hash = hashlib.sha256(json.dumps(branch_freeze_payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]

    print("STEP 36 SNAPSHOT ADAPTER SELF-TEST")
    print("Validated the solver-side snapshot adapters into the runtime-monitor schemas.")
    print("V2 branch-freeze metadata:")
    print("  branch_id =", branch_metadata["branch_id"])
    print("  branch_freeze_hash =", branch_freeze_hash)
    print("  pre_target_freeze =", str(branch_metadata["pre_target_freeze"]).lower())
    print("  target_blind =", str(branch_metadata["target_blind"]).lower())
    print("  no_post_residual_refit =", str(branch_metadata["no_post_residual_refit"]).lower())
    print("  boundary_class =", branch_metadata["boundary_class"])
    print("  interpretation =", branch_metadata["interpretation"])
    print("Supported adapter families:")
    print("  4D wavefunction snapshot -> full_4d monitor schema")
    print("  3D monopole state dump   -> projected_3d monitor schema")
    print("Wavefunction adapter pipeline:")
    print("  source_schema =", wave_summary["schema"]["source_schema"])
    print("  rel jx error =", wave_jx_rel)
    print("  rel jy error =", wave_jy_rel)
    print("  rel jz error =", wave_jz_rel)
    print("  max |R_cont| =", wave_summary["max_abs_R_cont"])
    print("  max |R_Pois_exact| =", wave_summary["max_abs_R_pois_exact"])
    print("  alpha_fit tail =", wave_summary["alpha_fit_tail_mean"])
    print("Monopole adapter pipeline:")
    print("  source_schema =", mono_summary["schema"]["source_schema"])
    print("  max |R_cont| =", mono_summary["max_abs_R_cont"])
    print("  max |R_Pois_exact| =", mono_summary["max_abs_R_pois_exact"])
    print("  alpha_fit tail =", mono_summary["alpha_fit_tail_mean"])
    print("  reconstructed source rel error =", reconstruct_rel)
    print("Adapter guardrails:")
    print("  missing W without --allow-uniform-W = rejected")
    print("  uniform W opt-in integral error =", uniform_W_integral_error)
    print("  uniform W opt-in nonunit-span value error =", uniform_W_value_error)
    print("  reconstructed monopole source without lambda = rejected")
    print("CLI:")
    print("  wavefunction -> python cfd_snapshot_adapters.py wavefunction-4d raw_wave.npz runtime_snapshot.npz")
    print("  wavefunction fallback W -> add --allow-uniform-W explicitly")
    print("  monopole     -> python cfd_snapshot_adapters.py monopole-3d raw_state.npz runtime_snapshot.npz")
    print("Interpretation:")
    print("  We can now feed either a 4D complex-wave snapshot or a 3D monopole state dump directly into the runtime monitor and fail-fast classifier without hand-editing arrays.")
    print("STATUS: PASS")


if __name__ == "__main__":
    main()
