from __future__ import annotations

import math

from stage1_solver import patha_extraction as pe
from stage1_solver import patha_b2a_bdg as b2a


def test_b2a_bdg_modes_are_direct_b1_moment_input() -> None:
    modes = [
        {
            "name": "m0",
            "lambda_B": 2.0,
            "overlap_I_eta_phi": 0.3,
            "coupling": 0.6,
            "varpi": 4.0,
            "profile": [0.0, 1.0],
        },
        {
            "name": "m1",
            "lambda_B": 0.5,
            "overlap_I_eta_phi": 0.8,
            "coupling": 0.4,
            "varpi": 5.0,
            "profile": [1.0, 0.0],
        },
    ]
    moments = pe.bdg_moments(modes)
    assert moments["B0"] == (0.6**2 / 4.0**2) + (0.4**2 / 5.0**2)
    assert moments["B2"] == (0.6**2 / 4.0**4) + (0.4**2 / 5.0**4)
    assert moments["B4"] == (0.6**2 / 4.0**6) + (0.4**2 / 5.0**6)


def test_b2a_structural_sanity_oracle_can_fail() -> None:
    good = {"bdg_moments": {"B0": 2.0, "B2": 1.0, "B4": 0.6}}
    bad = {"bdg_moments": {"B0": 1.0, "B2": 2.0, "B4": 1.0}}
    assert good["bdg_moments"]["B0"] * good["bdg_moments"]["B4"] - good["bdg_moments"]["B2"] ** 2 > 0.0
    assert bad["bdg_moments"]["B0"] * bad["bdg_moments"]["B4"] - bad["bdg_moments"]["B2"] ** 2 < 0.0


def test_b2a_dual_diff_detects_changed_coupling() -> None:
    base = {
        "bdg_modes": [
            {"varpi": 4.0, "coupling": 0.6},
            {"varpi": 5.0, "coupling": 0.4},
        ],
        "bdg_moments": {"B0": 0.1, "B2": 0.01, "B4": 0.001},
    }
    changed = {
        "bdg_modes": [
            {"varpi": 4.0, "coupling": 0.6},
            {"varpi": 5.0, "coupling": 0.5},
        ],
        "bdg_moments": {"B0": 0.1, "B2": 0.01, "B4": 0.001},
    }
    diffs = b2a.compare_engine_packets(base, changed)
    assert math.isclose(diffs["varpi"]["max_abs"], 0.0)
    assert diffs["coupling"]["max_abs"] == 0.09999999999999998


def test_b2a_modal_truncation_gate_can_fail_on_fat_tail() -> None:
    modes = [
        {"name": f"m{i}", "coupling": 1.0, "varpi": 1.0 + i}
        for i in range(8)
    ]
    exported = modes[:3]
    packet = {
        "grid": {"nr": 2, "nw": 2},
        "bdg_modes": exported,
        "bdg_moments": b2a._moments_with_stieltjes(exported),
        "diagnostics": {
            "modal_truncation": b2a._build_modal_truncation_metadata(
                modes,
                modes_to_export=3,
                tolerance=1.0e-3,
            )
        },
    }

    gate = b2a._modal_truncation_gate(packet)
    assert not gate["passed"]
    assert gate["checks"][0]["max_truncation_error_at_export"] > 1.0e-3


def test_b2a_consumer_gate_is_cross_engine_not_self_comparison() -> None:
    mma_packet = {
        "bdg_modes": [
            {"coupling": 0.6, "varpi": 4.0},
            {"coupling": 0.4, "varpi": 5.0},
        ],
    }
    good_moments = pe.bdg_moments(mma_packet["bdg_modes"])
    python_packet = {"bdg_moments": {**good_moments}}
    assert b2a._b1_cross_engine_consumer_check(
        mma_packet=mma_packet,
        python_packet=python_packet,
    )["passed"]

    python_packet["bdg_moments"] = {**good_moments, "B0": good_moments["B0"] * 1.01}
    assert not b2a._b1_cross_engine_consumer_check(
        mma_packet=mma_packet,
        python_packet=python_packet,
    )["passed"]


def test_b2a_dual_engine_gate_requires_abs_and_rel() -> None:
    passing_row = {
        "varpi": {"max_abs": 1.0e-12, "max_rel": 1.0e-12},
        "coupling": {"max_abs": 1.0e-14, "max_rel": 1.0e-10},
        "moments": {"max_abs": 1.0e-18, "max_rel": 1.0e-11},
    }
    assert b2a._dual_engine_gate([passing_row])["passed"]

    abs_only_row = {
        "varpi": {"max_abs": 1.0e-12, "max_rel": 1.0e-12},
        "coupling": {"max_abs": 1.0e-14, "max_rel": 1.0e-2},
        "moments": {"max_abs": 1.0e-18, "max_rel": 1.0e-11},
    }
    assert not b2a._dual_engine_gate([abs_only_row])["passed"]
