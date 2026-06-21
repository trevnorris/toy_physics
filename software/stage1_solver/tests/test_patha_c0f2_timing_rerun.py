from __future__ import annotations

import json

import numpy as np

from stage1_solver import patha_c0f2_timing_rerun as c0f2


def test_c0f2_depth_sequence_contains_required_fine_targets() -> None:
    sequence = c0f2.DEFAULT_C0F2_DEPTH_SEQUENCE

    assert 0.029125 in sequence
    assert 0.0290625 in sequence
    assert 0.029 in sequence
    assert 0.0285 in sequence
    assert 0.028 in sequence
    assert min(sequence) < 0.028
    assert all(left > right for left, right in zip(sequence, sequence[1:]))


def test_c0f2_extrapolation_uses_deepest_recent_measured_steps() -> None:
    table = [
        {
            "tau": 0.03,
            "wall_clock_seconds": 100.0,
            "accepted_default_success": True,
            "timing_instrumented": True,
        },
        {
            "tau": 0.0295,
            "wall_clock_seconds": 120.0,
            "accepted_default_success": True,
            "timing_instrumented": True,
        },
        {
            "tau": 0.029,
            "wall_clock_seconds": 140.0,
            "accepted_default_success": True,
            "timing_instrumented": True,
        },
    ]

    estimate = c0f2._extrapolate(table, target_tau=0.028)

    assert estimate["status"] == "MEASURED"
    assert estimate["model"] == "recent_median_per_fine_tau_step"
    assert np.isclose(estimate["median_recent_step_seconds"], 120.0)
    assert np.isclose(estimate["median_tau_step"], 0.0005)
    assert estimate["remaining_steps"] == 3


def test_c0f2_seed_import_marks_prior_rows_uninstrumented(tmp_path) -> None:
    state_path = tmp_path / "state.npz"
    np.savez_compressed(state_path, state=np.zeros(3))
    checkpoint_path = tmp_path / "old.json"
    checkpoint_path.write_text(
        json.dumps(
            {
                "tau_attempts": [
                    {
                        "target_tau": 0.03,
                        "final_physical_converged": True,
                        "state_artifact": str(state_path),
                    },
                    {
                        "target_tau": 0.029,
                        "final_physical_converged": False,
                        "state_artifact": None,
                    },
                ]
            }
        )
        + "\n",
        encoding="utf-8",
    )
    config = c0f2.C0f2Config(
        run_root=tmp_path / "run",
        report_path=tmp_path / "report.md",
        json_path=tmp_path / "run" / "result.json",
        seed_from_c0f_checkpoint=True,
        c0f_checkpoint_path=checkpoint_path,
    )

    payload = c0f2._load_checkpoint(config)

    assert payload is not None
    assert payload["phase_status"] == "seeded_from_prior_genuine_c0f_checkpoint"
    assert len(payload["tau_attempts"]) == 1
    assert payload["tau_attempts"][0]["timing_instrumented"] is False
    assert payload["scope_guard"]["prefer_existing_b2c_background_predictor"] is False
