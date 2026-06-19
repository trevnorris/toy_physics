from __future__ import annotations

import copy
import json
import math
from pathlib import Path

from stage1_solver import patha_b2c_calibration as b2c
from stage1_solver import patha_extraction as pe


def _coeff_for_tau(tau: float, *, b0: float = 4.0e-5, n0: float = 2.0e-8) -> dict[str, float]:
    return {
        "K": pe.hooke_kappa_hat() * tau,
        "M": 1.0,
        "B0": b0,
        "B2": 1.0e-6,
        "B4": 1.0e-7,
        "Z0": 2.0e-8,
        "Z2": -1.0e-8,
        "Z4": 5.0e-9,
        "N0": n0,
        "N2": -5.0e-9,
        "N4": 3.0e-9,
    }


def _sample(tau: float, coeff: dict[str, float]) -> dict:
    return {
        "tau": float(tau),
        "direct_coefficients": coeff,
        "observables": b2c.assemble_primary(coeff),
    }


def _evaluation_record(tau: float, kind: str, coeff: dict[str, float]) -> dict:
    assembled = b2c.assemble_dual(coeff)
    return {
        "tau": float(tau),
        "kind": kind,
        "content_hash": f"test-{kind}-{tau}",
        "direct_coefficients": coeff,
        "observables": {
            key: float(assembled["primary"][key])
            for key in ("D0", "P0", "R_norm", "R_pole", "P2", "P4")
        },
        "assembly": assembled,
        "stable_side": assembled["stable_side"],
        "paths": {},
        "r1_trace": {
            "closed_background": {
                "content_hash": f"bg-{kind}-{tau}",
                "converged": True,
                "residual_linf": 1.0e-9,
                "solver_message": "test converged background",
            }
        },
    }


def test_b2c_secondary_assembly_matches_b1_without_calling_lane_extract() -> None:
    tau = 0.03
    coeff = _coeff_for_tau(tau)
    primary = b2c.assemble_primary(coeff)
    secondary = b2c.assemble_secondary(coeff)
    dual = b2c.compare_assembly(primary, secondary)
    assert dual["passed"]
    assert math.isclose(primary["P4"], secondary["P4"], rel_tol=1.0e-11, abs_tol=1.0e-15)
    assert math.isclose(primary["R_pole"], secondary["R_pole"], rel_tol=1.0e-11, abs_tol=1.0e-15)


def test_b2c_stable_side_filter_rejects_d0_nonpositive() -> None:
    control = b2c.stable_side_negative_control()
    assert control["failed_as_expected"]
    assert control["D0"] <= 0.0


def test_b2c_confirmation_gate_fails_mislocated_tau() -> None:
    control = b2c.confirmation_negative_control()
    assert control["failed_as_expected"]
    assert abs(control["R_norm"]) > 1.0e-6


def test_b2c_frozen_root_moves_under_injected_drift_negative_control() -> None:
    control = b2c.negative_control_injected_drift()
    assert control["failed_as_expected"]
    assert control["relative_shift"] > 1.0e-2


def test_b2c_root_find_requires_real_stable_sign_bracket() -> None:
    coeff_hi = _coeff_for_tau(1.0)
    coeff_lo = _coeff_for_tau(0.03)
    samples = [_sample(0.03, coeff_lo), _sample(1.0, coeff_hi)]
    assert samples[0]["observables"]["R_norm"] < 0.0
    assert samples[1]["observables"]["R_norm"] < 0.0
    try:
        b2c.root_find_primary(samples)
    except ValueError as exc:
        assert "sign bracket" in str(exc)
    else:
        raise AssertionError("root finder accepted an unbracketed same-sign sample set")


def test_b2c_calibration_status_does_not_come_from_frozen_prior(tmp_path: Path) -> None:
    run_root = tmp_path / "run"
    eval_dir = run_root / "evaluations"
    eval_dir.mkdir(parents=True)
    seed = _evaluation_record(1.0, "seed_validated", _coeff_for_tau(1.0))
    real_floor = _evaluation_record(0.03, "resolved_floor", _coeff_for_tau(0.03, b0=5.0e-4, n0=4.0e-7))
    for record in (seed, real_floor):
        path = eval_dir / f"patha_b2c_{record['kind']}_tau_{b2c._format_tau(record['tau'])}.json"
        path.write_text(json.dumps(record), encoding="utf-8")

    result = b2c.calibrate_from_records(
        run_root=run_root,
        report_path=tmp_path / "report.md",
    )
    bundle = json.loads(Path(result["bundle_path"]).read_text(encoding="utf-8"))
    assert bundle["frozen_negative_control"]["tau_frozen"] < bundle["tau_floor"]["tau_floor"]
    assert bundle["status"] != "root_below_tau_floor"
    assert "Frozen-background root prior" not in bundle["finding"]
    assert bundle["root_agreement"] is None


def test_b2c_independent_root_agrees_on_bracketed_synthetic_samples() -> None:
    base = _coeff_for_tau(1.0)
    frozen_root = b2c.frozen_background_root(base)
    low_tau = frozen_root * 0.99999
    high_tau = frozen_root * 1.00020
    coeff_low = _coeff_for_tau(low_tau)
    coeff_high = _coeff_for_tau(high_tau)
    samples = [_sample(low_tau, coeff_low), _sample(high_tau, coeff_high)]
    assert samples[0]["observables"]["D0"] > 0.0
    assert samples[0]["observables"]["R_norm"] > 0.0
    assert samples[1]["observables"]["R_norm"] < 0.0
    primary = b2c.root_find_primary(samples)
    secondary = b2c.root_find_secondary(samples)
    assert abs(primary["tau"] - secondary["tau"]) <= max(1.0e-12 * primary["tau"], 1.0e-15)
    assert abs(primary["observables"]["R_norm"]) < 1.0e-8


def test_b2c_coefficient_error_propagation_uses_recorded_shape() -> None:
    coeff = _coeff_for_tau(0.03)
    rel = {key: 0.0 for key in b2c.DIRECT_KEYS}
    rel["N0"] = 0.004
    rel["B0"] = 0.01
    propagated = b2c.propagate_observable_errors(coeff, rel)
    assert propagated["absolute_errors"]["R_norm"] > 0.0
    assert propagated["tau_root_absolute_error_local"] > 0.0
    assert "N0" in propagated["coefficient_contributions"]
