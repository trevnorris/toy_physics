"""Step 7 numerical error-budget composition for the coupled branch.

This module intentionally performs no new branch solves.  It composes the
target-blind Step 4 surrogate-observable table with the Step 5 boundary floor
and Step 6 conservation/Gauss floor.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .config import HarnessConfig
from .convergence import OBSERVABLE_LABELS


NULL_DIAGNOSTIC = "null diagnostic"
SOLVER_FLOOR_DIAGNOSTIC = "solver-floor diagnostic"
EXPECTED_ORDER_CONVERGENCE = "expected-order convergence"

STEP4_PROVENANCE = {
    "source_step": "Step 4",
    "source_function": "stage1_solver.convergence.run_step4",
    "source_report_path": "software/stage1_solver/reports/step4_convergence_study.md",
    "source_commit": "c9b8b2c",
    "source_table_path": (
        "software/stage1_solver/runs/step4_convergence_study/"
        "step4_coupled_branch_grid_convergence/convergence_table.json"
    ),
    "config_hash": "6ca33829da5aec8b",
}

STEP5_PROVENANCE = {
    "source_step": "Step 5",
    "source_function": "stage1_solver.boundary_characterization.run_step5",
    "source_report_path": "software/stage1_solver/reports/step5_boundary_characterization.md",
    "source_commit": "63cd885",
    "source_table_path": (
        "software/stage1_solver/runs/step5_boundary_characterization/"
        "step5_boundary_characterization/boundary_characterization_table.json"
    ),
    "config_hash": "fd9bc1b134624e96",
}

STEP6_PROVENANCE = {
    "source_step": "Step 6",
    "source_function": "stage1_solver.conservation_diagnostics.run_step6",
    "source_report_path": "software/stage1_solver/reports/step6_conservation_diagnostics.md",
    "source_commit": "4a03797",
    "source_table_path": (
        "software/stage1_solver/runs/step6_conservation_diagnostics/"
        "step6_conservation_diagnostics/conservation_diagnostics_table.json"
    ),
    "config_hash": "684c53f5562074e4",
    "diagnostics_digest": "9752e5eb227f21f1",
}

# Recorded from STEP4_PROVENANCE["source_function"] at STEP4_PROVENANCE["source_commit"].
# These are not target values; they are the target-blind Step 4 numerical diagnostics.
RECORDED_STEP4_NOISE_FLOOR = {
    "solver_residual_floor_linf": 1.068892085953621e-08,
    "mass_constraint_floor": 2.0846830928178406e-09,
    "field_self_difference_floor_reached": False,
    "last_raw_field_relative_l2_change": 3.477500567058835e-04,
    "preliminary_numerical_floor": 1.068892085953621e-08,
}

RECORDED_STEP4_OBSERVABLE_SUMMARY = [
    {
        "observable": "density_mass",
        "label": "density mass integral",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 1.0000000013968622,
        "last_observed_order": 2.0051090084989176,
        "richardson_estimate": 1.0000000013945571,
        "finest_error_estimate": 2.30504504372675e-12,
        "verdict": SOLVER_FLOOR_DIAGNOSTIC,
        "diagnosis": "Newton mass constraint; differences read the solver floor, not discretization.",
    },
    {
        "observable": "peak_density",
        "label": "peak density",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 0.28005452830904265,
        "last_observed_order": 2.00543272096225,
        "richardson_estimate": 0.280106587688765,
        "finest_error_estimate": 5.205937972235786e-05,
        "verdict": EXPECTED_ORDER_CONVERGENCE,
        "diagnosis": "Pointwise extrema are sensitive to cell-center placement and coarse throat resolution.",
    },
    {
        "observable": "min_density",
        "label": "minimum density",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 2.3093962312588454e-07,
        "last_observed_order": 2.020453223287057,
        "richardson_estimate": 2.4526154852439287e-09,
        "finest_error_estimate": 2.284870076406406e-07,
        "verdict": EXPECTED_ORDER_CONVERGENCE,
        "diagnosis": "Pointwise extrema are sensitive to cell-center placement and coarse throat resolution.",
    },
    {
        "observable": "raw_field_l2_norm",
        "label": "raw field L2 norm",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 1.0002731839015613,
        "last_observed_order": 2.008130403351871,
        "richardson_estimate": 1.0002730696184459,
        "finest_error_estimate": 1.1428311541550329e-07,
        "verdict": EXPECTED_ORDER_CONVERGENCE,
        "diagnosis": "Raw-field integral on the shared coupled branch.",
    },
    {
        "observable": "scalar_gauge_l2",
        "label": "A0 L2 norm",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 0.02337607827897189,
        "last_observed_order": 2.005846932072963,
        "richardson_estimate": 0.023371179333987156,
        "finest_error_estimate": 4.898944984732534e-06,
        "verdict": EXPECTED_ORDER_CONVERGENCE,
        "diagnosis": "Coupled gauge/current response can show reduced order from open Robin boundaries.",
    },
    {
        "observable": "spatial_gauge_l2",
        "label": "spatial gauge L2 norm",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 0.0,
        "last_observed_order": None,
        "richardson_estimate": 0.0,
        "finest_error_estimate": 0.0,
        "verdict": NULL_DIAGNOSTIC,
        "diagnosis": (
            "This raw channel is identically zero on the completed placeholder branch; "
            "no order is measured."
        ),
    },
    {
        "observable": "spatial_current_l2",
        "label": "spatial current L2 norm",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 0.0,
        "last_observed_order": None,
        "richardson_estimate": 0.0,
        "finest_error_estimate": 0.0,
        "verdict": NULL_DIAGNOSTIC,
        "diagnosis": (
            "This raw channel is identically zero on the completed placeholder branch; "
            "no order is measured."
        ),
    },
    {
        "observable": "field_energy_like_integral",
        "label": "field-energy-like integral",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 2.044752362966379,
        "last_observed_order": 2.1697277676106155,
        "richardson_estimate": 2.0451155868505992,
        "finest_error_estimate": 3.6322388422016516e-04,
        "verdict": EXPECTED_ORDER_CONVERGENCE,
        "diagnosis": "Gradient-weighted integral; boundary and coarse-grid throat terms can reduce order.",
    },
    {
        "observable": "chemical_potential",
        "label": "chemical potential",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 2.0935765184503525,
        "last_observed_order": 1.9985255488712903,
        "richardson_estimate": 2.0938252935259962,
        "finest_error_estimate": 2.487750756436924e-04,
        "verdict": EXPECTED_ORDER_CONVERGENCE,
        "diagnosis": "Raw-field integral on the shared coupled branch.",
    },
    {
        "observable": "final_residual_linf",
        "label": "final residual Linf",
        "finest_grid": "convergence_l4_nr_96_nw_64",
        "finest_dof": 30721,
        "finest_value": 9.284885083005179e-09,
        "last_observed_order": 1.757521551237741,
        "richardson_estimate": 9.276180167212441e-09,
        "finest_error_estimate": 8.704915792737904e-12,
        "verdict": SOLVER_FLOOR_DIAGNOSTIC,
        "diagnosis": "Newton/GMRES stopping floor; useful as the numerical floor diagnostic.",
    },
]

# Recorded from Step 5/6 outputs at the provenance commits.
RECORDED_STEP5_BOUNDARY_RELATIVE_FLOOR = 6.850289509215066e-05
RECORDED_STEP6_CONSERVATION_RELATIVE_FLOOR = 8.23773483448239e-03

BOUNDARY_APPLIES_TO = frozenset(
    OBSERVABLE_LABELS.keys()
) - {"spatial_gauge_l2", "spatial_current_l2", "final_residual_linf"}
CONSERVATION_APPLIES_TO = frozenset(
    {"density_mass", "scalar_gauge_l2", "field_energy_like_integral"}
)
NULL_OBSERVABLES = frozenset({"spatial_gauge_l2", "spatial_current_l2"})

OBSERVABLE_CLASSES = {
    "density_mass": "charge/mass integral",
    "peak_density": "mass/density pointwise",
    "min_density": "mass/density pointwise",
    "raw_field_l2_norm": "raw-field aggregate",
    "scalar_gauge_l2": "gauge/Maxwell coupled",
    "spatial_gauge_l2": "null spatial gauge",
    "spatial_current_l2": "null spatial current",
    "field_energy_like_integral": "energy coupled",
    "chemical_potential": "branch scalar",
    "final_residual_linf": "solver-floor diagnostic",
}


@dataclass(frozen=True)
class ErrorFloorScales:
    solver_abs: float
    discretization_relative_floor: float
    boundary_relative_floor: float
    conservation_relative_floor: float


@dataclass(frozen=True)
class CombinationSensitivity:
    material_ratio_threshold: float = 1.05


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def recorded_prior_results() -> dict[str, dict[str, Any]]:
    """Return provenance-tagged recorded Step 4/5/6 outputs used by default."""

    return {
        "step4": {
            "passed": True,
            "config_hash": STEP4_PROVENANCE["config_hash"],
            "provenance": {
                **STEP4_PROVENANCE,
                "recorded_input_kind": "code-pinned run_step4 output snapshot",
            },
            "noise_floor": dict(RECORDED_STEP4_NOISE_FLOOR),
            "observable_summary": [dict(row) for row in RECORDED_STEP4_OBSERVABLE_SUMMARY],
        },
        "step5": {
            "passed": True,
            "config_hash": STEP5_PROVENANCE["config_hash"],
            "provenance": {
                **STEP5_PROVENANCE,
                "recorded_input_kind": "code-pinned run_step5 output scalar",
            },
            "boundary_error_metric": {
                "max_relative_l2": RECORDED_STEP5_BOUNDARY_RELATIVE_FLOOR,
                "reference": "fixed-window interior raw-field L2 signal magnitude",
            },
        },
        "step6": {
            "passed": True,
            "config_hash": STEP6_PROVENANCE["config_hash"],
            "provenance": {
                **STEP6_PROVENANCE,
                "recorded_input_kind": "code-pinned run_step6 output scalar",
            },
            "conservation_error_metric": {
                "gauss_closure_relative_floor": RECORDED_STEP6_CONSERVATION_RELATIVE_FLOOR,
                "source": (
                    "max finest-level sponge_on independent Gauss-closure relative "
                    "residual across nested surfaces"
                ),
            },
        },
    }


def _provenance(results: dict[str, Any], fallback: dict[str, Any]) -> dict[str, Any]:
    provenance = results.get("provenance")
    if provenance is not None:
        return dict(provenance)
    return dict(fallback)


def solver_floor_from_step4(step4_results: dict[str, Any]) -> float:
    floor = step4_results["noise_floor"]["preliminary_numerical_floor"]
    if floor is None or float(floor) <= 0.0:
        raise ValueError("Step 4 preliminary numerical floor must be positive")
    return float(floor)


def discretization_relative_floor_from_step4(step4_results: dict[str, Any]) -> float:
    floor = step4_results["noise_floor"]["last_raw_field_relative_l2_change"]
    if floor is None or float(floor) < 0.0:
        raise ValueError("Step 4 raw-field self-convergence floor must be non-negative")
    return float(floor)


def boundary_relative_floor_from_step5(step5_results: dict[str, Any]) -> float:
    floor = step5_results["boundary_error_metric"]["max_relative_l2"]
    if floor is None or float(floor) < 0.0:
        raise ValueError("Step 5 boundary floor must be non-negative")
    return float(floor)


def conservation_relative_floor_from_step6(step6_results: dict[str, Any]) -> float:
    metric = step6_results.get("conservation_error_metric")
    if metric is not None:
        floor = metric["gauss_closure_relative_floor"]
        if floor is None or float(floor) < 0.0:
            raise ValueError("Step 6 conservation floor must be non-negative")
        return float(floor)

    rows = [
        row
        for row in step6_results["gauss_closure_rows"]
        if row["mode"] == "sponge_on" and row["reconstruction"] == "independent_center_gradient"
    ]
    if not rows:
        raise ValueError("Step 6 result has no sponge_on independent Gauss-closure rows")
    finest = max(row["level"] for row in rows)
    floor = max(abs(float(row["relative_residual"])) for row in rows if row["level"] == finest)
    return floor


def combine_uncertainty(
    *,
    solver_floor: float,
    discretization_abs: float,
    boundary_abs: float,
    conservation_abs: float,
) -> dict[str, float]:
    """RSS independent numerical axes, then enforce the solver floor."""

    axes = [discretization_abs, boundary_abs, conservation_abs]
    if any(value < 0.0 for value in [solver_floor, *axes]):
        raise ValueError("uncertainty components must be non-negative")
    rss_unfloored = math.sqrt(sum(value * value for value in axes))
    rss_total = max(rss_unfloored, solver_floor)
    max_alternative = max(solver_floor, *axes)
    sum_bound = solver_floor + sum(axes)
    return {
        "rss_unfloored": rss_unfloored,
        "rss_total": rss_total,
        "max_alternative": max_alternative,
        "sum_bound": sum_bound,
        "rss_over_max": rss_total / max(max_alternative, 1.0e-300),
    }


def compose_observable_budget(
    summary: dict[str, Any],
    floors: ErrorFloorScales,
    *,
    sensitivity: CombinationSensitivity = CombinationSensitivity(),
) -> dict[str, Any]:
    observable = summary["observable"]
    value = float(summary["finest_value"])
    abs_value = abs(value)
    verdict = summary["verdict"]
    is_null = verdict == NULL_DIAGNOSTIC or observable in NULL_OBSERVABLES

    if is_null:
        return {
            "observable": observable,
            "label": OBSERVABLE_LABELS.get(observable, summary.get("label", observable)),
            "observable_class": OBSERVABLE_CLASSES.get(observable, "surrogate observable"),
            "finest_grid": summary.get("finest_grid"),
            "finest_dof": summary.get("finest_dof"),
            "finest_value": value,
            "u_solver": 0.0,
            "u_disc": 0.0,
            "u_boundary": 0.0,
            "u_conservation": 0.0,
            "u_total": 0.0,
            "relative_uncertainty": None,
            "u_max_alternative": 0.0,
            "u_sum_bound": 0.0,
            "rss_over_max": 1.0,
            "dominant_component": "null",
            "verdict": NULL_DIAGNOSTIC,
            "floor_limited": False,
            "material_rss_vs_max": False,
            "applied_components": [],
            "budget_note": "",
        }

    disc_raw = float(summary["finest_error_estimate"] or 0.0)
    u_boundary = floors.boundary_relative_floor * abs_value if observable in BOUNDARY_APPLIES_TO else 0.0
    u_conservation = (
        floors.conservation_relative_floor * abs_value
        if observable in CONSERVATION_APPLIES_TO
        else 0.0
    )
    combination = combine_uncertainty(
        solver_floor=floors.solver_abs,
        discretization_abs=disc_raw,
        boundary_abs=u_boundary,
        conservation_abs=u_conservation,
    )
    component_values = {
        "solver": floors.solver_abs,
        "discretization": disc_raw,
        "boundary": u_boundary,
        "conservation": u_conservation,
    }
    dominant_component = max(component_values, key=component_values.__getitem__)
    applied_components = ["solver", "discretization"]
    if u_boundary > 0.0:
        applied_components.append("boundary")
    if u_conservation > 0.0:
        applied_components.append("conservation")

    floor_limited = disc_raw < floors.solver_abs
    relative_uncertainty = combination["rss_total"] / abs_value if abs_value > 0.0 else None
    return {
        "observable": observable,
        "label": OBSERVABLE_LABELS.get(observable, summary.get("label", observable)),
        "observable_class": OBSERVABLE_CLASSES.get(observable, "surrogate observable"),
        "finest_grid": summary.get("finest_grid"),
        "finest_dof": summary.get("finest_dof"),
        "finest_value": value,
        "u_solver": floors.solver_abs,
        "u_disc": disc_raw,
        "u_boundary": u_boundary,
        "u_conservation": u_conservation,
        "u_total": combination["rss_total"],
        "relative_uncertainty": relative_uncertainty,
        "u_max_alternative": combination["max_alternative"],
        "u_sum_bound": combination["sum_bound"],
        "rss_over_max": combination["rss_over_max"],
        "dominant_component": dominant_component,
        "verdict": verdict,
        "floor_limited": floor_limited,
        "material_rss_vs_max": (
            combination["rss_over_max"] >= sensitivity.material_ratio_threshold
        ),
        "applied_components": applied_components,
        "budget_note": (
            "Step 4 verdict/floor_limited describes discretization; composed "
            "budget is conservation-dominated."
            if observable == "density_mass"
            else ""
        ),
    }


def _ordered_observable_summary(step4_results: dict[str, Any]) -> list[dict[str, Any]]:
    rows_by_name = {row["observable"]: dict(row) for row in step4_results["observable_summary"]}
    expected = set(OBSERVABLE_LABELS)
    observed = set(rows_by_name)
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise ValueError(f"Step 4 observable set mismatch; missing={missing}, extra={extra}")
    return [rows_by_name[name] for name in OBSERVABLE_LABELS]


def _observable_set_matches_step4(step4_results: dict[str, Any]) -> bool:
    return (
        set(row["observable"] for row in step4_results["observable_summary"])
        == set(OBSERVABLE_LABELS)
    )


def component_floor_rows(
    *,
    floors: ErrorFloorScales,
    step4_provenance: dict[str, Any],
    step5_provenance: dict[str, Any],
    step6_provenance: dict[str, Any],
) -> list[dict[str, Any]]:
    return [
        {
            "component": "solver",
            "value": floors.solver_abs,
            "units": "absolute",
            "source_step": "Step 4",
            "source_field": "noise_floor.preliminary_numerical_floor",
            "source_report_path": step4_provenance["source_report_path"],
            "source_commit": step4_provenance["source_commit"],
        },
        {
            "component": "discretization",
            "value": floors.discretization_relative_floor,
            "units": "relative raw-field self-difference; per-observable u_disc uses Step 4 finest_error_estimate",
            "source_step": "Step 4",
            "source_field": (
                "noise_floor.last_raw_field_relative_l2_change and "
                "observable_summary[*].finest_error_estimate"
            ),
            "source_report_path": step4_provenance["source_report_path"],
            "source_commit": step4_provenance["source_commit"],
        },
        {
            "component": "boundary",
            "value": floors.boundary_relative_floor,
            "units": "relative interior raw-field L2",
            "source_step": "Step 5",
            "source_field": "boundary_error_metric.max_relative_l2",
            "source_report_path": step5_provenance["source_report_path"],
            "source_commit": step5_provenance["source_commit"],
        },
        {
            "component": "conservation",
            "value": floors.conservation_relative_floor,
            "units": "relative Gauss closure",
            "source_step": "Step 6",
            "source_field": (
                "max finest-level sponge_on independent Gauss-closure relative_residual"
            ),
            "source_report_path": step6_provenance["source_report_path"],
            "source_commit": step6_provenance["source_commit"],
        },
    ]


def _governing_floor_by_class(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["observable_class"], []).append(row)
    class_rows = []
    for observable_class in sorted(grouped):
        class_budget_rows = grouped[observable_class]
        dominant = max(class_budget_rows, key=lambda row: row["u_total"])
        class_rows.append(
            {
                "observable_class": observable_class,
                "governing_component": dominant["dominant_component"],
                "governing_observable": dominant["observable"],
                "governing_u_total": dominant["u_total"],
                "null_class": all(row["verdict"] == NULL_DIAGNOSTIC for row in class_budget_rows),
            }
        )
    return class_rows


def _diagnostics_digest(results: dict[str, Any]) -> str:
    payload = {
        "component_floors": results["component_floors"],
        "observable_budget": results["observable_budget"],
        "governing_floor_by_class": results["governing_floor_by_class"],
        "pass_checks": results["pass_checks"],
        "asserted_checks": results["asserted_checks"],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=_json_default)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def run_step7(
    config: HarnessConfig | None = None,
    *,
    step4_results: dict[str, Any] | None = None,
    step5_results: dict[str, Any] | None = None,
    step6_results: dict[str, Any] | None = None,
    sensitivity: CombinationSensitivity = CombinationSensitivity(),
) -> dict[str, Any]:
    if config is None:
        config = HarnessConfig(
            run_root="software/stage1_solver/runs/step7_error_budget",
            report_path="software/stage1_solver/reports/step7_error_budget.md",
        )

    recorded = recorded_prior_results()
    step4 = recorded["step4"] if step4_results is None else step4_results
    step5 = recorded["step5"] if step5_results is None else step5_results
    step6 = recorded["step6"] if step6_results is None else step6_results

    step4_provenance = _provenance(step4, STEP4_PROVENANCE)
    step5_provenance = _provenance(step5, STEP5_PROVENANCE)
    step6_provenance = _provenance(step6, STEP6_PROVENANCE)
    floors = ErrorFloorScales(
        solver_abs=solver_floor_from_step4(step4),
        discretization_relative_floor=discretization_relative_floor_from_step4(step4),
        boundary_relative_floor=boundary_relative_floor_from_step5(step5),
        conservation_relative_floor=conservation_relative_floor_from_step6(step6),
    )
    observable_set_matches_step4 = _observable_set_matches_step4(step4)
    observable_budget = [
        compose_observable_budget(row, floors, sensitivity=sensitivity)
        for row in _ordered_observable_summary(step4)
    ]
    component_floors = component_floor_rows(
        floors=floors,
        step4_provenance=step4_provenance,
        step5_provenance=step5_provenance,
        step6_provenance=step6_provenance,
    )
    non_null_rows = [row for row in observable_budget if row["verdict"] != NULL_DIAGNOSTIC]
    pass_checks = {
        "non_null_uncertainties_floor_at_solver": all(
            row["u_total"] >= floors.solver_abs for row in non_null_rows
        ),
        "solver_floor_limited_rows_flagged": all(
            row["floor_limited"]
            for row in non_null_rows
            if row["u_disc"] < floors.solver_abs
        ),
        "null_sectors_remain_null": all(
            row["u_total"] == 0.0 and row["dominant_component"] == "null"
            for row in observable_budget
            if row["verdict"] == NULL_DIAGNOSTIC
        ),
        "conservation_floor_scoped": all(
            (row["u_conservation"] > 0.0) == (row["observable"] in CONSERVATION_APPLIES_TO)
            for row in observable_budget
            if row["verdict"] != NULL_DIAGNOSTIC
        ),
        "boundary_floor_scoped": all(
            (row["u_boundary"] > 0.0) == (row["observable"] in BOUNDARY_APPLIES_TO)
            for row in observable_budget
            if row["verdict"] != NULL_DIAGNOSTIC
        ),
        "observable_set_matches_step4": observable_set_matches_step4,
    }
    asserted_checks = {
        "target_blind_surrogate_budget_only_not_a_physics_gate": True,
        "physical_export_permitted_is_false_not_a_physics_gate": True,
        "prior_step4_passed_not_a_physics_gate": bool(step4.get("passed", True)),
        "prior_step5_passed_not_a_physics_gate": bool(step5.get("passed", True)),
        "prior_step6_passed_not_a_physics_gate": bool(step6.get("passed", True)),
    }
    asserted_check_notes = {
        "target_blind_surrogate_budget_only_not_a_physics_gate": (
            "Step 7 composes only target-blind surrogate observables and imports no "
            "target/reference/export map; asserted by module construction, not an "
            "independent physics measurement."
        ),
        "physical_export_permitted_is_false_not_a_physics_gate": (
            "Step 7 emits no physical packet; the external export guard lives in the "
            "firewalled physical model and is untouched - asserted by construction, "
            "not read here."
        ),
        "prior_step4_passed_not_a_physics_gate": (
            f"Default pinned provenance records Step 4 PASS at commit "
            f"{step4_provenance['source_commit']}; this is not a live Step 7 rerun."
        ),
        "prior_step5_passed_not_a_physics_gate": (
            f"Default pinned provenance records Step 5 PASS at commit "
            f"{step5_provenance['source_commit']}; this is not a live Step 7 rerun."
        ),
        "prior_step6_passed_not_a_physics_gate": (
            f"Default pinned provenance records Step 6 PASS at commit "
            f"{step6_provenance['source_commit']}; this is not a live Step 7 rerun."
        ),
    }
    results: dict[str, Any] = {
        "config": config.to_dict(),
        "config_hash": config.config_hash(),
        "method": {
            "scope": (
                "Numerical-only error budget on Step 4 target-blind surrogate observables; "
                "no physical section-G extraction map, no section-H targets, no physical packet."
            ),
            "combination_rule": (
                "RSS of independent numerical axes (discretization, boundary, "
                "conservation where applicable), then floored at the Step 4 "
                "solver/Newton-GMRES floor."
            ),
            "independence_assumption": (
                "Grid self-convergence, boundary truncation, and Gauss/conservation "
                "closure are treated as independent numerical axes. The solver floor "
                "is a hard reporting floor, not an independent axis to RSS twice."
            ),
            "max_alternative": (
                "For sensitivity, each row also reports max(u_solver, u_disc, "
                "u_boundary, u_conservation)."
            ),
            "material_rss_vs_max_ratio": sensitivity.material_ratio_threshold,
            "boundary_applies_to": sorted(BOUNDARY_APPLIES_TO),
            "conservation_applies_to": sorted(CONSERVATION_APPLIES_TO),
            "null_observables": sorted(NULL_OBSERVABLES),
            "export_guard": {"physical_export_permitted": False},
        },
        "provenance": {
            "step4": step4_provenance,
            "step5": step5_provenance,
            "step6": step6_provenance,
        },
        "component_floors": component_floors,
        "observable_budget": observable_budget,
        "governing_floor_by_class": _governing_floor_by_class(observable_budget),
        "material_sensitivity_rows": [
            row for row in observable_budget if row["material_rss_vs_max"]
        ],
        "pass_checks": pass_checks,
        "asserted_checks": asserted_checks,
        "asserted_check_notes": asserted_check_notes,
        "passed": all(pass_checks.values()),
    }
    results["diagnostics_digest"] = _diagnostics_digest(results)
    table_path = Path(config.run_root) / "step7_error_budget" / "error_budget_table.json"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    results["machine_readable_table"] = str(table_path)
    table_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, default=_json_default),
        encoding="utf-8",
    )
    return results


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6e}"
    if isinstance(value, (tuple, list)):
        return "[" + ", ".join(str(item) for item in value) + "]"
    return str(value)


def _table(headers: list[str], rows: list[dict[str, Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(row.get(header)) for header in headers) + " |")
    return "\n".join(lines)


def _component_report_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "component": row["component"],
            "value": row["value"],
            "units": row["units"],
            "source_step": row["source_step"],
            "source_commit": row["source_commit"],
            "source_field": row["source_field"],
        }
        for row in rows
    ]


def _budget_report_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "label": row["label"],
            "observable_class": row["observable_class"],
            "finest_value": row["finest_value"],
            "u_solver": row["u_solver"],
            "u_disc": row["u_disc"],
            "u_boundary": row["u_boundary"],
            "u_conservation": row["u_conservation"],
            "u_total": row["u_total"],
            "relative_uncertainty": row["relative_uncertainty"],
            "dominant_component": row["dominant_component"],
            "verdict": row["verdict"],
            "floor_limited": row["floor_limited"],
            "u_max_alternative": row["u_max_alternative"],
            "rss_over_max": row["rss_over_max"],
        }
        for row in rows
    ]


def write_step7_report(results: dict[str, Any], path: str) -> Path:
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("# Step 7 Numerical Error Budget")
    lines.append("")
    lines.append(f"Overall engineering gate: {'PASS' if results['passed'] else 'FAIL'}")
    lines.append(f"Config hash: `{results['config_hash']}`")
    lines.append(f"Diagnostics digest: `{results['diagnostics_digest']}`")
    lines.append("")
    lines.append(
        "**Scope framing:** this is a numerical-validity budget on the Step 4 "
        "target-blind surrogate observables only. It is not a physical section-G "
        "observable extraction, uses no extraction map or section-H target, emits no "
        "physical packet, and keeps `physical_export_permitted = False`."
    )
    lines.append("")
    lines.append(
        "**Honest limitations:** the budget covers only discretization, "
        "Newton-GMRES solver floor, boundary truncation, and conservation/Gauss "
        "closure. It excludes free_choice / physical-parameter uncertainty "
        "(GATE A: none frozen, and this will dominate eventual real observables) "
        "and model-form / parent-action uncertainty. No precision is claimed "
        "below the solver floor; null spatial sectors remain null diagnostics."
    )
    lines.append("")
    lines.append("## Component Noise Floors")
    lines.append("")
    lines.append(
        _table(
            ["component", "value", "units", "source_step", "source_commit", "source_field"],
            _component_report_rows(results["component_floors"]),
        )
    )
    lines.append("")
    lines.append("Governing floor by observable class:")
    lines.append("")
    lines.append(
        _table(
            [
                "observable_class",
                "governing_component",
                "governing_observable",
                "governing_u_total",
                "null_class",
            ],
            results["governing_floor_by_class"],
        )
    )
    lines.append("")
    lines.append("## Combination Rule")
    lines.append("")
    lines.append(results["method"]["combination_rule"])
    lines.append(results["method"]["independence_assumption"])
    lines.append(results["method"]["max_alternative"])
    lines.append(
        "Boundary relative floors are converted as `u_boundary(O) = boundary_rel * |O|` "
        "for non-null solution observables except the residual diagnostic. Conservation "
        "relative floors are converted as `u_conservation(O) = conservation_rel * |O|` "
        "only for the charge/mass integral, scalar-gauge/Maxwell, and energy-coupled "
        "surrogates."
    )
    lines.append(
        "For `density_mass`, the conservation floor is scoped through the localized "
        "Gauss-law proportionality `surface_integral Z F^{i0} dA = "
        "mu0*q_star*volume_integral rho dV`: the charge density "
        "`J^0 = q_star*rho` uses the same `rho` field as the mass integral "
        "`volume_integral rho dV`, so the Gauss-closure relative residual genuinely "
        "bounds that integral rather than switching categories."
    )
    lines.append(
        "`floor_limited = true` means the Step 4 finite-grid estimate is below the "
        "solver floor; the composed total can still be larger when boundary or "
        "conservation components also apply."
    )
    density_mass_row = next(
        row for row in results["observable_budget"] if row["observable"] == "density_mass"
    )
    lines.append(
        "For `density_mass`, `verdict = solver-floor diagnostic` and "
        "`floor_limited = true` describe the Step 4 convergence/discretization "
        "behavior only; the composed Step 7 budget is conservation-dominated "
        f"(`u_total = {_fmt(density_mass_row['u_total'])}`), so the row is not a "
        "claim of ~1e-8 precision."
    )
    lines.append(
        "Near-zero pointwise/residual surrogates such as `min_density` and "
        "`final_residual_linf` can carry relative uncertainty near or above 1; "
        "this is the budget correctly declining to claim precision on a near-zero "
        "quantity, not a defect."
    )
    lines.append("")
    lines.append("## Per-Observable Budget")
    lines.append("")
    lines.append(
        _table(
            [
                "label",
                "observable_class",
                "finest_value",
                "u_solver",
                "u_disc",
                "u_boundary",
                "u_conservation",
                "u_total",
                "relative_uncertainty",
                "dominant_component",
                "verdict",
                "floor_limited",
                "u_max_alternative",
                "rss_over_max",
            ],
            _budget_report_rows(results["observable_budget"]),
        )
    )
    lines.append("")
    lines.append("## RSS-vs-Max Sensitivity")
    lines.append("")
    if results["material_sensitivity_rows"]:
        lines.append(
            "Rows where RSS exceeds the conservative max alternative by the configured "
            f"{_fmt(results['method']['material_rss_vs_max_ratio'])} ratio are surfaced:"
        )
        lines.append("")
        lines.append(
            _table(
                ["label", "u_total", "u_max_alternative", "rss_over_max", "dominant_component"],
                [
                    {
                        "label": row["label"],
                        "u_total": row["u_total"],
                        "u_max_alternative": row["u_max_alternative"],
                        "rss_over_max": row["rss_over_max"],
                        "dominant_component": row["dominant_component"],
                    }
                    for row in results["material_sensitivity_rows"]
                ],
            )
        )
    else:
        lines.append("No per-observable RSS-vs-max divergence exceeds the configured material ratio.")
    lines.append("")
    lines.append("## Pass Checks")
    lines.append("")
    lines.append("Counted gates only; `passed` is the conjunction of these checks.")
    lines.append("")
    for key, value in results["pass_checks"].items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")
    lines.append("")
    lines.append("Asserted-by-construction checks, reported but not counted as physics gates:")
    for key, value in results["asserted_checks"].items():
        note = results["asserted_check_notes"].get(key, "")
        suffix = f" - {note}" if note else ""
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}{suffix}")
    lines.append("")
    lines.append("## Provenance")
    lines.append("")
    for key in ("step4", "step5", "step6"):
        provenance = results["provenance"][key]
        lines.append(
            f"- {key}: `{provenance['source_function']}` at `{provenance['source_commit']}`, "
            f"report `{provenance['source_report_path']}`, recorded input "
            f"`{provenance.get('recorded_input_kind', 'run output')}`."
        )
    lines.append("")
    lines.append("## Reproduction")
    lines.append("")
    lines.append("```bash")
    lines.append(
        "timeout 600 env PYTHONPATH=software/stage1_solver/src "
        "python -m stage1_solver.error_budget_harness"
    )
    lines.append(
        "timeout 600 env PYTHONPATH=software/stage1_solver/src "
        "pytest software/stage1_solver/tests/test_stage1_solver.py"
    )
    lines.append("```")
    lines.append("")
    lines.append(f"Machine-readable table: `{results['machine_readable_table']}`.")
    lines.append(
        "Target-blindness: Step 7 imports no benchmark targets, no references, no "
        "extraction map, and no physical export path; relative scales are the "
        "observable's own magnitude or prior target-blind aggregate norms."
    )
    lines.append(
        "Export guard: `physical_export_permitted` remains false; no physical packet is emitted."
    )
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path
