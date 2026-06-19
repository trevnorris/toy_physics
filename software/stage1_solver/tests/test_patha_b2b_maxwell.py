from __future__ import annotations

import copy

import numpy as np

from stage1_solver import patha_b2b_maxwell as b2b


def _toy_background() -> dict:
    r = np.linspace(0.1, 1.4, 5)
    w = np.linspace(0.1, 1.75, 5)
    rr, ww = np.meshgrid(r, w, indexing="ij")
    psi = 0.08 * np.exp(-0.4 * rr**2) * (1.0 + 0.02 * np.cos(np.pi * ww / 1.85))
    return {
        "content_hash": "toy-background",
        "geometry": {"a": 1.0, "L": 1.85},
        "constants": {
            "tau": 1.0,
            "gauge_charge": 0.35,
            "particle_mass": 1.0,
            "mu0": 1.0,
        },
        "grid": {
            "r_max": 1.6,
            "r_centers": r.tolist(),
            "w_centers": w.tolist(),
        },
        "fields": {
            "psi_R0": psi.tolist(),
            "psi_I0": np.zeros_like(psi).tolist(),
            "A_00": (1.0e-4 * (1.0 + rr + 0.5 * ww)).tolist(),
            "A_r0": np.zeros_like(psi).tolist(),
            "A_w0": np.zeros_like(psi).tolist(),
            "R0_w": (1.0 + 1.0e-3 * np.sin(np.pi * w / 1.85)).tolist(),
        },
        "derived": {"Z_w": (1.0 + 0.1 * np.cos(np.pi * w / 1.85)).tolist()},
        "residuals": {"closed_stationary_linf": 1.0e-10, "self_consistent": True},
    }


def _toy_bdg_packet() -> dict:
    points = np.linspace(0.0, 1.85, 9)
    return {
        "content_hash": "toy-bdg",
        "grid": {"profile_points": points.tolist()},
        "wall": {"K": 7.7, "M": 1.0},
        "bdg_moments": {"B0": 1.0e-5, "B2": 1.0e-6, "B4": 1.0e-7},
        "bdg_modes": [
            {
                "name": "m0",
                "varpi": 4.5,
                "coupling": 0.02,
                "profile": np.sin(np.pi * points / 1.85).tolist(),
            },
            {
                "name": "m1",
                "varpi": 5.5,
                "coupling": 0.01,
                "profile": np.sin(2.0 * np.pi * points / 1.85).tolist(),
            },
        ],
    }


def test_b2b_dual_engine_gate_requires_abs_and_rel() -> None:
    passing = {
        "coefficients": {"max_abs": 1.0e-12, "max_rel": 1.0e-10},
        "per_coefficient": {
            key: {"max_abs": 1.0e-12, "max_rel": 1.0e-10}
            for key in b2b.COEFF_KEYS
        },
    }
    assert b2b._dual_engine_gate([passing])["passed"]

    abs_only = {
        "coefficients": {"max_abs": 1.0e-12, "max_rel": 1.0e-1},
        "per_coefficient": {
            key: {"max_abs": 1.0e-12, "max_rel": 1.0e-10}
            for key in b2b.COEFF_KEYS
        },
    }
    abs_only["per_coefficient"]["N0"] = {"max_abs": 1.0e-12, "max_rel": 1.0e-1}
    assert not b2b._dual_engine_gate([abs_only])["passed"]


def test_b2b_dual_engine_detects_perturbed_coefficient() -> None:
    base = {"coefficients": {key: 1.0e-6 for key in b2b.COEFF_KEYS}}
    perturbed = copy.deepcopy(base)
    perturbed["coefficients"]["Z2"] *= 1.20
    diffs = b2b.compare_engine_packets(base, perturbed)
    gate = b2b._dual_engine_gate([{**diffs, "tau": 1.0, "grid": "toy"}])
    assert diffs["per_coefficient"]["Z2"]["max_rel"] > 1.0e-3
    assert not gate["passed"]


def test_b2b_basis_invariance_gate_can_fail_when_rotation_is_incomplete() -> None:
    background = _toy_background()
    bdg = _toy_bdg_packet()
    assert b2b.basis_invariance_check(background, bdg)["passed"]
    leaked = b2b.basis_invariance_check(background, bdg, inject_port_leak=True)
    assert not leaked["passed"]
    assert leaked["port_leak_negative_control"]["failed_as_expected"]
    assert leaked["branch_difference"]["different_from_invariant"]
    assert leaked["max_tested_extraction_relative_delta"] > 1.0e-3


def test_b2b_convergence_gate_rejects_unconverged_sequence() -> None:
    control = b2b.convergence_negative_control()
    assert control["failed_as_expected"]
    assert control["grow_one_shrink_another"]
    assert control["N4_truncation_increments"][-1] > control["N4_truncation_increments"][-2]


def test_b2b_sweep_convergence_is_per_coefficient_not_max_flip() -> None:
    base = {key: 1.0 for key in b2b.COEFF_KEYS}
    packets = [
        b2b._synthetic_coeff_packet(
            {**base, "Z0": 1.0, "N4": 1.0},
            grid=(31, 13),
            window=0.028,
            radial_scale=3.0,
        ),
        b2b._synthetic_coeff_packet(
            {**base, "Z0": 1.20, "N4": 1.05},
            grid=(39, 15),
            window=0.028,
            radial_scale=4.0,
        ),
        b2b._synthetic_coeff_packet(
            {**base, "Z0": 1.26, "N4": 1.16},
            grid=(47, 17),
            window=0.028,
            radial_scale=5.0,
        ),
    ]
    table = b2b._sweep_table(packets, label_key="radial_scale")
    summary = b2b._sweep_convergence_summary(table)
    assert not summary["passed"]
    assert summary["per_coefficient"]["Z0"]["shrinking"]
    assert not summary["per_coefficient"]["N4"]["shrinking"]


def test_b2b_self_consistency_negative_control_rejects_smoke_background() -> None:
    control = b2b.self_consistency_negative_control(_toy_background())
    assert control["failed_as_expected"]


def test_b2b_tau_sensitivity_negative_control_rejects_relabelled_bundle() -> None:
    packet = {
        "tau": 1.0,
        "constants": {"tau": 1.0},
        "coefficients": {key: 1.0e-6 for key in b2b.COEFF_KEYS},
    }
    control = b2b.tau_sensitivity_negative_control(packet)
    assert control["failed_as_expected"]


def test_b2b_physicality_gate_rejects_nonpositive_n0() -> None:
    packet = {
        "coefficients": {
            "Z0": 1.0e-6,
            "Z2": 1.0e-7,
            "Z4": 1.0e-8,
            "N0": -1.0e-6,
            "N2": 1.0e-7,
            "N4": 1.0e-8,
        },
        "rows": [
            {
                "sigmaCons": {"re": 1.0e-6, "im": 0.0},
                "minusImSigmaRet": 1.0e-9,
                "outgoingFlux": 1.0e-9,
            }
        ],
    }
    gate = b2b._outgoing_physicality_gate(packet)
    assert not gate["passed"]
    assert not gate["N0_positive"]


def test_b2b_consumer_compatibility_is_cross_engine_direct_shape() -> None:
    bdg = _toy_bdg_packet()
    primary = {"coefficients": {key: 2.0e-6 for key in b2b.COEFF_KEYS}}
    independent = {"coefficients": {key: 2.0e-6 for key in b2b.COEFF_KEYS}}
    assert b2b._consumer_compatibility_gate(
        primary_packet=primary,
        independent_packet=independent,
        bdg_packet=bdg,
    )["passed"]
    independent["coefficients"]["N4"] *= 1.25
    assert not b2b._consumer_compatibility_gate(
        primary_packet=primary,
        independent_packet=independent,
        bdg_packet=bdg,
    )["passed"]
