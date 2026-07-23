#!/usr/bin/env python3
"""Ledger stage041 SymPy audit: the NG5 one-medium reducibility gate.

Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.  The
production result is a characterized departure at the bookkeeping/labelling
level.  It is neither a dynamical proof of one medium nor a discovery of a
second substance.

Tooth-local runtime ablation uses ``LEDGER_STAGE041_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
import os
from typing import Any, Iterable, Mapping

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE041_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

ROUTE_CONJUNCTS = (
    "required_fields_present",
    "named_solve_in_provenance",
    "result_status_allowed",
    "finite_listed_missing_objects",
    "target_identity",
    "joint_target_identity",
    "target_blind",
    "falsifiers_clear",
)

CONJUNCT_TEETH = (
    "CTRL_ROUTE_CONJUNCT_REQUIRED_FIELDS_PRESENT",
    "CTRL_ROUTE_CONJUNCT_NAMED_SOLVE_IN_PROVENANCE",
    "CTRL_ROUTE_CONJUNCT_RESULT_STATUS_ALLOWED",
    "CTRL_ROUTE_CONJUNCT_FINITE_LISTED_MISSING_OBJECTS",
    "CTRL_ROUTE_CONJUNCT_TARGET_IDENTITY",
    "CTRL_ROUTE_CONJUNCT_JOINT_TARGET_IDENTITY",
    "CTRL_ROUTE_CONJUNCT_TARGET_BLIND",
    "CTRL_ROUTE_CONJUNCT_FALSIFIERS_CLEAR",
)

TOOTH_ORDER = (
    "ORIGIN_CLASSIFICATION",
    "FULL_ORIGIN_LEDGER",
    "FULL_ROUTE_EVALUATION_LEDGER",
    "PRODUCTION_VERDICT",
    "DRIFT_COUNT",
    "FREEDOM_STATE_RHO_BR",
    "FREEDOM_STATE_C_HU",
    "LOCATION_CLOSURE",
    "ARENA_LABELS",
    "LINEAGE_FINDING",
    "REDUCTION_STATUS",
    "ANTI_ABSORPTION_DRIFT_SIDE",
    "CTRL_AB_DELETE_REGISTRY",
    "CTRL_ROUTE_BLANK_ROUTE_A",
    "CTRL_ROUTE_FIELD_BLANK_TARGET_BLIND",
    "CTRL_ROUTE_FIELD_BLANK_MISSING_OBJECTS",
    *CONJUNCT_TEETH,
    "CTRL_CALIBRATION_ABLATION_Q_E",
    "CTRL_IRREDUCIBLE_SYNTHETIC",
    "CTRL_REDUCIBLE_DERIVED_SYNTHETIC",
    "CTRL_CONTRADICTION",
    "CTRL_ROUTE_EVAL_RECORDED",
    "CTRL_RESIDUAL_MULTIPLIER_ABLATION",
    "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA",
    "DIMENSION_HOMOGENEITY",
    "DUAL_ENGINE_TERMS",
    "VERDICT_REDERIVATION",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "ORIGIN_CLASSIFICATION":
        "make Future-Embedding-Overlap a valid deferred route before classifying C_hu",
    "FULL_ORIGIN_LEDGER":
        "redirect B_eff's dependency to a non-irreducible row before classification",
    "FULL_ROUTE_EVALUATION_LEDGER":
        "give c_gamma an absent route hint so its computed reason changes while its origin does not",
    "PRODUCTION_VERDICT":
        "make C_hu reducible-deferred in the private verdict witness",
    "DRIFT_COUNT":
        "inject xi_active into the private classified ledger feeding the computed count",
    "FREEDOM_STATE_RHO_BR":
        "remove Route-A in the private rho_br freedom-state witness",
    "FREEDOM_STATE_C_HU":
        "make Future-Embedding-Overlap fully valid and REGISTERED_DEFERRED",
    "LOCATION_CLOSURE":
        "inject a production row labelled unassigned",
    "ARENA_LABELS":
        "relabel C_hu from the seam to the still-admissible 3D brane surface",
    "LINEAGE_FINDING":
        "replace the production lineage facts by one same active object with a free residual multiplier",
    "REDUCTION_STATUS":
        "drop c_E while constructing the reduction-status pending map",
    "ANTI_ABSORPTION_DRIFT_SIDE":
        "rename the surviving historical rho_br entry to rho_B0 before retirement and projection",
    "CTRL_AB_DELETE_REGISTRY":
        "neutralize the all-route deletion control",
    "CTRL_ROUTE_BLANK_ROUTE_A":
        "neutralize the Route-A deletion control",
    "CTRL_ROUTE_FIELD_BLANK_TARGET_BLIND":
        "neutralize the target_blind field-deletion control",
    "CTRL_ROUTE_FIELD_BLANK_MISSING_OBJECTS":
        "neutralize the missing_objects field-deletion control",
    "CTRL_ROUTE_CONJUNCT_REQUIRED_FIELDS_PRESENT":
        "neutralize deletion of the otherwise-unused required source field",
    "CTRL_ROUTE_CONJUNCT_NAMED_SOLVE_IN_PROVENANCE":
        "neutralize the named-solve=False isolating control",
    "CTRL_ROUTE_CONJUNCT_RESULT_STATUS_ALLOWED":
        "neutralize the rejected-status isolating control",
    "CTRL_ROUTE_CONJUNCT_FINITE_LISTED_MISSING_OBJECTS":
        "neutralize the non-list missing_objects isolating control",
    "CTRL_ROUTE_CONJUNCT_TARGET_IDENTITY":
        "neutralize dropping only c_E from Route-A targets",
    "CTRL_ROUTE_CONJUNCT_JOINT_TARGET_IDENTITY":
        "neutralize the sentinel joint target",
    "CTRL_ROUTE_CONJUNCT_TARGET_BLIND":
        "neutralize target_blind=False while retaining the field",
    "CTRL_ROUTE_CONJUNCT_FALSIFIERS_CLEAR":
        "neutralize insertion of one Route-A falsifier",
    "CTRL_CALIBRATION_ABLATION_Q_E":
        "neutralize removal of the Q_E calibration anchor",
    "CTRL_IRREDUCIBLE_SYNTHETIC":
        "suppress xi_active injection",
    "CTRL_REDUCIBLE_DERIVED_SYNTHETIC":
        "change p_syn local_status away from SYNTHETIC_EARNED_RELATION",
    "CTRL_CONTRADICTION":
        "omit the cone locks so the freedom_tie system has an exact positive SAT witness",
    "CTRL_ROUTE_EVAL_RECORDED":
        "drop c_gamma's RouteEvaluation record from the active-record projection",
    "CTRL_RESIDUAL_MULTIPLIER_ABLATION":
        "close the free branch residual multiplier to 1",
    "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA":
        "suppress loc_sentinel injection",
    "DIMENSION_HOMOGENEITY":
        "change the cited speed dimension from L T^-1 to L^2 T^-1",
    "DUAL_ENGINE_TERMS":
        "drop the locally computed lock_B witness from the canonical inventory",
    "VERDICT_REDERIVATION":
        "leave one future route's named-solve provenance false in the conditional witness",
    "SOURCE_TO_STAGE_MANIFEST":
        "drop one source predicate and mis-scope a second predicate",
}

if len(TOOTH_ORDER) != 35:
    raise RuntimeError("stage041 tooth declaration is not exactly 35")


class AuditFailure(AssertionError):
    """A named stage predicate failed."""

    def __init__(self, predicate: str, detail: str = "") -> None:
        super().__init__(predicate)
        self.predicate = predicate
        self.detail = detail


def expect_bool(name: str, condition: bool, evidence: Any = None) -> None:
    global PASS_COUNT, FAIL_COUNT
    if bool(condition):
        PASS_COUNT += 1
        print(f"PASS  {name}")
        return
    FAIL_COUNT += 1
    print(f"FIRST_FAILURE={name}")
    if ACTIVE_MUTATION == name:
        print(f"FIRED_AT_OWN_ASSERT={name}")
    print(f"FAIL  {name}: residual = 1")
    if evidence is not None:
        print(f"      evidence = {evidence}")
    raise AuditFailure(name)


def section(text: str) -> None:
    print("")
    print(text)
    print("-" * len(text))


ALLOWED_ROUTE_STATUSES = frozenset({"SOLVED_PASS", "REGISTERED_DEFERRED"})
REJECTED_ROUTE_STATUSES = frozenset(
    {"FAILED", "BY_TUNING", "ABSENT", "PROMISSORY_ONLY"}
)
ARENAS = frozenset(
    {"4D bulk", "3D brane surface", "throat/embedding seam"}
)
EXPECTED_TRIO = ("rho_B0", "chi_c", "C_hu")
EXPECTED_ROUTE_A_TARGETS = ("rho_br", "mu_R", "c_E")
EXPECTED_VERDICT = (
    "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu})"
)


@dataclass(frozen=True)
class Row:
    param: str
    local_status: str
    active: bool = False
    out_of_active: bool = False
    dependencies: tuple[str, ...] = ()
    route_hint: str | None = None
    calibration_anchor: str | None = None
    calibrated_geometry: bool = False
    location: str = "unassigned"


@dataclass(frozen=True)
class DriftEntry:
    key: str
    name: str


def base_routes() -> dict[str, dict[str, Any]]:
    return {
        "Route-A": {
            "name": "Route-A",
            "source": "pathA_40:ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
            "result_status": "REGISTERED_DEFERRED",
            "named_solve_in_provenance": True,
            "missing_objects": ["R1", "R2", "R3", "R4", "R5"],
            "targets": ["rho_br", "mu_R", "c_E"],
            "required_joint_targets": ["rho_br", "mu_R"],
            "target_blind": True,
            "falsifiers": [],
        },
        "Future-Compression-4D-to-3D": {
            "name": "Future-Compression-4D-to-3D",
            "source": "named_future_reduction_route:not_currently_registered",
            "result_status": "PROMISSORY_ONLY",
            "named_solve_in_provenance": False,
            "missing_objects": [
                "compression-sector 4D-to-3D nonlinear brane solve"
            ],
            "targets": ["rho_B0", "chi_c"],
            "required_joint_targets": ["rho_B0", "chi_c"],
            "target_blind": True,
            "falsifiers": [],
        },
        "Future-Embedding-Overlap": {
            "name": "Future-Embedding-Overlap",
            "source": "named_future_reduction_route:not_currently_registered",
            "result_status": "PROMISSORY_ONLY",
            "named_solve_in_provenance": False,
            "missing_objects": [
                "embedding-overlap nonlinear throat solve deriving C_hu"
            ],
            "targets": ["C_hu"],
            "required_joint_targets": ["C_hu"],
            "target_blind": True,
            "falsifiers": [],
        },
    }


def copied_routes() -> dict[str, dict[str, Any]]:
    return {
        name: {
            key: list(value) if isinstance(value, list) else value
            for key, value in fact.items()
        }
        for name, fact in base_routes().items()
    }


def mutate_route(
    routes: Mapping[str, Mapping[str, Any]],
    route_name: str,
    **changes: Any,
) -> dict[str, dict[str, Any]]:
    out = {
        name: {
            key: list(value) if isinstance(value, list) else value
            for key, value in fact.items()
        }
        for name, fact in routes.items()
    }
    out[route_name].update(changes)
    return out


def valid_future_routes(
    *, future_status: str, route_a_status: str = "REGISTERED_DEFERRED"
) -> dict[str, dict[str, Any]]:
    routes = copied_routes()
    for route_name in (
        "Future-Compression-4D-to-3D",
        "Future-Embedding-Overlap",
    ):
        routes[route_name]["named_solve_in_provenance"] = True
        routes[route_name]["result_status"] = future_status
    routes["Route-A"]["named_solve_in_provenance"] = True
    routes["Route-A"]["result_status"] = route_a_status
    return routes


def build_rows(
    *,
    b_eff_dependency_mutation: bool = False,
    c_gamma_route_mutation: bool = False,
    c_hu_arena_mutation: bool = False,
    add_xi: bool = False,
    add_p_syn: bool = False,
    p_syn_status_mutation: bool = False,
    add_location_sentinel: bool = False,
) -> tuple[Row, ...]:
    rows = [
        Row("rho", "BASE_SUBSTRATE", location="4D bulk"),
        Row("K", "BASE_SUBSTRATE", location="4D bulk"),
        Row("m", "BASE_SUBSTRATE", location="4D bulk"),
        Row("a", "BASE_SUBSTRATE", location="4D bulk"),
        Row(
            "ell_g",
            "ACCEPTED_GEOMETRY_SUBSTRATE",
            location="throat/embedding seam",
        ),
        Row(
            "g_ell(w)",
            "ACCEPTED_PROFILE_GIVEN_ell_g",
            location="throat/embedding seam",
        ),
        Row(
            "varrho_br[rho]",
            "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT",
            out_of_active=True,
            location="3D brane surface",
        ),
        Row(
            "Sigma_n[rho]",
            "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT",
            out_of_active=True,
            location="3D brane surface",
        ),
        Row(
            "delta_Sigma[rho]",
            "PATHA25_CLOSED_DENSITY_SMECTIC_OBJECT",
            out_of_active=True,
            location="3D brane surface",
        ),
        Row(
            "rho_br",
            "POSTULATED_SHEAR_SURFACE_REDUCTION_TARGET",
            active=True,
            route_hint="Route-A",
            location="3D brane surface",
        ),
        Row(
            "mu_R",
            "POSTULATED_SHEAR_SURFACE_REDUCTION_TARGET",
            active=True,
            route_hint="Route-A",
            location="3D brane surface",
        ),
        Row(
            "c_E",
            "ROUTE_A_CONE_LOCK_REDUCTION_TARGET",
            active=True,
            route_hint="Route-A",
            location="throat/embedding seam",
        ),
        Row(
            "c_gamma",
            "DEPENDENT_SPEED",
            active=True,
            dependencies=("rho_br", "mu_R"),
            route_hint="Absent-C-Gamma-Route"
            if c_gamma_route_mutation
            else None,
            location="3D brane surface",
        ),
        Row(
            "rho_B0",
            "NO_CURRENT_REGISTERED_REDUCTION",
            active=True,
            route_hint="Future-Compression-4D-to-3D",
            location="3D brane surface",
        ),
        Row(
            "chi_c",
            "NO_CURRENT_REGISTERED_REDUCTION",
            active=True,
            route_hint="Future-Compression-4D-to-3D",
            location="3D brane surface",
        ),
        Row(
            "B_eff",
            "DEPENDENT_ON_COMPRESSION_INPUTS",
            active=True,
            dependencies=("c_gamma",)
            if b_eff_dependency_mutation
            else ("rho_B0", "chi_c"),
            location="3D brane surface",
        ),
        Row(
            "C_hu",
            "NO_CURRENT_REGISTERED_REDUCTION",
            active=True,
            route_hint="Future-Embedding-Overlap",
            location="3D brane surface"
            if c_hu_arena_mutation
            else "throat/embedding seam",
        ),
        Row(
            "Q_E",
            "DECLARED_CALIBRATION_ANCHOR",
            active=True,
            calibration_anchor="Q_E",
            location="throat/embedding seam",
        ),
        Row(
            "ell",
            "CALIBRATED_GEOMETRY_SOURCE_INPUT",
            active=True,
            calibrated_geometry=True,
            location="throat/embedding seam",
        ),
        Row(
            "b",
            "CALIBRATED_GEOMETRY_SOURCE_INPUT",
            active=True,
            calibrated_geometry=True,
            location="throat/embedding seam",
        ),
        Row(
            "M_h",
            "CALIBRATED_GEOMETRY_SOURCE_INPUT",
            active=True,
            calibrated_geometry=True,
            location="throat/embedding seam",
        ),
        Row(
            "K_h",
            "DEPENDENT_STIFFNESS",
            active=True,
            dependencies=("M_h", "c_E"),
            location="throat/embedding seam",
        ),
        Row(
            "q_h",
            "DEPENDENT_PROJECTION",
            active=True,
            dependencies=("Q_E", "b", "ell"),
            location="throat/embedding seam",
        ),
    ]
    rows.extend(
        Row(
            param,
            "PATHA25_DRIVER_OUT_OF_ACTIVE_SHEAR_SURFACE_NG5",
            out_of_active=True,
            location="3D brane surface",
        )
        for param in (
            "c_L1",
            "c_L2",
            "A_R",
            "k_R",
            "lambda_Cdiv",
            "chi_Cpin",
            "J_Pu",
            "kappa_Pu",
        )
    )
    rows.extend(
        (
            Row(
                "lambda_Pu",
                "FROZEN_UNUSED_OUT_OF_ACTIVE_NG5",
                out_of_active=True,
                location="3D brane surface",
            ),
            Row(
                "Omega_w",
                "FROZEN_UNUSED_OUT_OF_ACTIVE_NG5",
                out_of_active=True,
                location="3D brane surface",
            ),
            Row(
                "lambda_N",
                "WALL_INTERNAL_EXCLUDED_FROM_ACTIVE_NG5",
                out_of_active=True,
                location="throat/embedding seam",
            ),
            Row(
                "lambda_tau",
                "WALL_INTERNAL_EXCLUDED_FROM_ACTIVE_NG5",
                out_of_active=True,
                location="throat/embedding seam",
            ),
            Row(
                "Nu",
                "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
                out_of_active=True,
                location="throat/embedding seam",
            ),
            Row(
                "a_T",
                "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
                out_of_active=True,
                location="throat/embedding seam",
            ),
            Row(
                "a_Tp",
                "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
                out_of_active=True,
                location="throat/embedding seam",
            ),
            Row(
                "a_L",
                "SOURCE_NORMALIZATION_OUTSIDE_ACTIVE_ORIGIN_COMPETITION",
                out_of_active=True,
                location="throat/embedding seam",
            ),
        )
    )
    if add_xi:
        rows.append(
            Row(
                "xi_active",
                "SYNTHETIC_ACTIVE_INPUT",
                active=True,
                location="3D brane surface",
            )
        )
    if add_p_syn:
        rows.append(
            Row(
                "p_syn",
                "SYNTHETIC_ACTIVE_INPUT"
                if p_syn_status_mutation
                else "SYNTHETIC_EARNED_RELATION",
                active=True,
                location="4D bulk",
            )
        )
    if add_location_sentinel:
        rows.append(Row("loc_sentinel", "BASE_SUBSTRATE", location="unassigned"))
    return tuple(rows)


def route_evaluate(
    row: Row, routes: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    route_name = row.route_hint
    if route_name is None:
        return {
            "route_name": None,
            "valid": False,
            "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
            "result_status": None,
            "checks": {},
        }
    fact = routes.get(route_name)
    if fact is None:
        return {
            "route_name": route_name,
            "valid": False,
            "reason": "ROUTE_BLANKED_OR_ABSENT",
            "result_status": None,
            "checks": {"mapped_route_exists": False},
        }
    required = (
        "name",
        "source",
        "result_status",
        "named_solve_in_provenance",
        "missing_objects",
        "targets",
        "target_blind",
        "falsifiers",
    )
    missing_fields = tuple(key for key in required if key not in fact)
    status = fact.get("result_status")
    missing_objects = fact.get("missing_objects")
    targets = fact.get("targets")
    required_joint = fact.get("required_joint_targets", [])
    falsifiers = fact.get("falsifiers")
    checks = {
        "required_fields_present": not missing_fields,
        "named_solve_in_provenance":
            fact.get("named_solve_in_provenance") is True,
        "result_status_allowed": status in ALLOWED_ROUTE_STATUSES,
        "finite_listed_missing_objects": isinstance(missing_objects, list),
        "target_identity":
            isinstance(targets, list) and row.param in targets,
        "joint_target_identity":
            isinstance(targets, list)
            and all(target in targets for target in required_joint),
        "target_blind": fact.get("target_blind") is True,
        "falsifiers_clear":
            isinstance(falsifiers, list) and not falsifiers,
    }
    valid = all(checks[name] for name in ROUTE_CONJUNCTS)
    if valid:
        reason = "VALID_REGISTERED_ROUTE"
    elif status in REJECTED_ROUTE_STATUSES:
        reason = f"REJECTED_RESULT_STATUS_{status}"
    elif missing_fields:
        reason = "MISSING_REQUIRED_ROUTE_FIELD"
    else:
        reason = "ROUTE_EVALUATION_FAILED"
    return {
        "route_name": route_name,
        "valid": valid,
        "reason": reason,
        "result_status": status,
        "checks": checks,
    }


def classify_rows(
    rows: Iterable[Row],
    routes: Mapping[str, Mapping[str, Any]],
    *,
    anchors: frozenset[str] = frozenset({"Q_E"}),
    geometry: frozenset[str] = frozenset({"ell", "b", "M_h"}),
) -> tuple[dict[str, Any], ...]:
    classified: list[dict[str, Any]] = []
    prior: dict[str, dict[str, Any]] = {}
    for row in rows:
        route_eval = route_evaluate(row, routes) if row.active else None
        if row.out_of_active:
            origin = "OUT_OF_ACTIVE_NG5"
        elif row.local_status in {
            "BASE_SUBSTRATE",
            "ACCEPTED_GEOMETRY_SUBSTRATE",
            "ACCEPTED_PROFILE_GIVEN_ell_g",
        }:
            origin = row.local_status
        elif row.dependencies:
            dependency_origins = {
                prior[dependency]["origin"]
                for dependency in row.dependencies
                if dependency in prior
            }
            origin = (
                "DEPENDENT_ON_IRREDUCIBLE"
                if row.param == "B_eff"
                and "IRREDUCIBLY_INDEPENDENT" in dependency_origins
                else "DEPENDENT"
            )
        elif row.calibration_anchor:
            origin = (
                "CALIBRATED_ANCHOR"
                if row.calibration_anchor in anchors
                else "IRREDUCIBLY_INDEPENDENT"
            )
        elif row.calibrated_geometry:
            origin = (
                "CALIBRATED_GEOMETRY_INPUT"
                if row.param in geometry
                else "IRREDUCIBLY_INDEPENDENT"
            )
        elif row.local_status == "SYNTHETIC_EARNED_RELATION":
            origin = "REDUCIBLE_DERIVED"
        elif (
            route_eval
            and route_eval["valid"]
            and route_eval["result_status"] == "SOLVED_PASS"
        ):
            origin = "REDUCIBLE_DERIVED"
        elif (
            route_eval
            and route_eval["valid"]
            and route_eval["result_status"] == "REGISTERED_DEFERRED"
        ):
            origin = "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED"
        else:
            origin = "IRREDUCIBLY_INDEPENDENT"
        record = {"row": row, "origin": origin, "route_evaluation": route_eval}
        classified.append(record)
        prior[row.param] = record
    return tuple(classified)


def run_case(
    *,
    rows: tuple[Row, ...] | None = None,
    routes: Mapping[str, Mapping[str, Any]] | None = None,
    anchors: frozenset[str] = frozenset({"Q_E"}),
    geometry: frozenset[str] = frozenset({"ell", "b", "M_h"}),
) -> tuple[dict[str, Any], ...]:
    return classify_rows(
        build_rows() if rows is None else rows,
        copied_routes() if routes is None else routes,
        anchors=anchors,
        geometry=geometry,
    )


def by_param(
    classified: Iterable[dict[str, Any]]
) -> dict[str, dict[str, Any]]:
    return {record["row"].param: record for record in classified}


def active_irreducible(
    classified: Iterable[dict[str, Any]]
) -> tuple[str, ...]:
    return tuple(
        record["row"].param
        for record in classified
        if record["row"].active
        and record["origin"] == "IRREDUCIBLY_INDEPENDENT"
    )


def sim_deferred_map(
    classified: Iterable[dict[str, Any]]
) -> dict[str, str]:
    result: dict[str, str] = {}
    for record in classified:
        evaluation = record["route_evaluation"]
        if (
            record["row"].active
            and record["origin"]
            == "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED"
            and evaluation
            and evaluation["route_name"]
        ):
            result[record["row"].param] = evaluation["route_name"]
    return result


def calibrated_map(
    classified: Iterable[dict[str, Any]]
) -> dict[str, str]:
    return {
        record["row"].param: record["origin"]
        for record in classified
        if record["origin"]
        in {"CALIBRATED_ANCHOR", "CALIBRATED_GEOMETRY_INPUT"}
    }


def decide_verdict(
    classified: tuple[dict[str, Any], ...],
    *,
    no_go: bool = False,
) -> dict[str, Any]:
    irreducible = active_irreducible(classified)
    sim_deferred = sim_deferred_map(classified)
    calibrated = calibrated_map(classified)
    if no_go:
        verdict = "NO_GO(cone-lock-feedback)"
        drift_count = 0
        reported_irreducible: tuple[str, ...] = ()
    elif irreducible:
        verdict = (
            "SECOND_MEDIUM_DRIFT(active_irreducible={"
            + ",".join(irreducible)
            + "})"
        )
        drift_count = len(irreducible)
        reported_irreducible = irreducible
    elif sim_deferred:
        sim_text = ",".join(
            f"{param}->{route}" for param, route in sim_deferred.items()
        )
        calibrated_text = ",".join(
            f"{param}->{origin}" for param, origin in calibrated.items()
        )
        verdict = (
            "ONE_MEDIUM_CONDITIONAL(sim-deferred={"
            + sim_text
            + "}; calibrated={"
            + calibrated_text
            + "})"
        )
        drift_count = 0
        reported_irreducible = ()
    else:
        verdict = "ONE_MEDIUM_CONSISTENT"
        drift_count = 0
        reported_irreducible = ()
    return {
        "verdict": verdict,
        "drift_count": drift_count,
        "active_irreducible": reported_irreducible,
        "sim_deferred": sim_deferred,
        "calibrated": calibrated,
    }


def freedom_state(record: dict[str, Any]) -> str:
    evaluation = record["route_evaluation"]
    param = record["row"].param
    if (
        evaluation
        and evaluation["valid"]
        and evaluation["result_status"] == "REGISTERED_DEFERRED"
    ):
        return f"FREEDOM_SIM_DEFERRED{{{evaluation['route_name']}}}"
    if (
        evaluation
        and evaluation["valid"]
        and evaluation["result_status"] == "SOLVED_PASS"
    ):
        return f"FREEDOM_REDUCED{{{evaluation['route_name']}}}"
    return f"FREEDOM_CERTIFIED_CURRENT_LEDGER{{{param}}}"


def location_closure(rows: Iterable[Row]) -> dict[str, Any]:
    offending = tuple(
        (row.param, row.location) for row in rows if row.location not in ARENAS
    )
    return {
        "no_fourth_arena": len(offending) == 0,
        "offending": offending,
    }


def arena_labels(rows: Iterable[Row]) -> dict[str, str]:
    trio = set(EXPECTED_TRIO)
    return {row.param: row.location for row in rows if row.param in trio}


@dataclass(frozen=True)
class LineageFacts:
    varrho_object: str
    rho_object: str
    varrho_active: bool
    rho_active: bool
    varrho_role: str
    rho_role: str
    varrho_dimension: tuple[int, int, int]
    rho_dimension: tuple[int, int, int]
    residual_multiplier: str | None
    route_a_pending: bool


def production_lineage_facts() -> LineageFacts:
    return LineageFacts(
        varrho_object="pathA_25_closed_density_smectic_varrho_br",
        rho_object="pathA_35_active_shear_surface_rho_br",
        varrho_active=False,
        rho_active=True,
        varrho_role="density_smectic_layer_kinetic",
        rho_role="shear_surface_kinetic",
        varrho_dimension=(1, -3, 0),
        rho_dimension=(1, -3, 0),
        residual_multiplier=None,
        route_a_pending=True,
    )


def residual_lineage_facts(multiplier: str) -> LineageFacts:
    return LineageFacts(
        varrho_object="same_active_shear_surface_object",
        rho_object="same_active_shear_surface_object",
        varrho_active=True,
        rho_active=True,
        varrho_role="shear_surface_kinetic",
        rho_role="shear_surface_kinetic",
        varrho_dimension=(1, -3, 0),
        rho_dimension=(1, -3, 0),
        residual_multiplier=multiplier,
        route_a_pending=True,
    )


def lineage_adjudication(facts: LineageFacts) -> tuple[str, str]:
    same_active_object = (
        facts.varrho_active
        and facts.rho_active
        and facts.varrho_object == facts.rho_object
    )
    same_role = facts.varrho_role == facts.rho_role
    same_dimension = facts.varrho_dimension == facts.rho_dimension
    residual_closed = facts.residual_multiplier in {None, "1"}
    same = (
        same_active_object
        and same_role
        and same_dimension
        and residual_closed
    )
    if same:
        return "SAME", "OVERCOUNT_OR_SMUGGLE_CONTROL"
    if same_active_object and not residual_closed:
        return (
            "DIFFERENT",
            "DIFFERENT_BRANE_INERTIA_OBJECTS_RESIDUAL_MULTIPLIER",
        )
    if (
        not facts.varrho_active
        and facts.rho_active
        and facts.varrho_object != facts.rho_object
        and facts.route_a_pending
    ):
        return "DIFFERENT", "NO_OVERCOUNT_ROUTE_A_PENDING"
    return "DIFFERENT", "DIFFERENT_BRANE_INERTIA_OBJECTS"


def lineage_from_ledger(
    classified: tuple[dict[str, Any], ...], facts: LineageFacts
) -> tuple[str, str]:
    route_a_pending = (
        by_param(classified)["rho_br"]["origin"]
        == "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED"
    )
    return lineage_adjudication(
        replace(facts, route_a_pending=route_a_pending)
    )


def reduction_status(
    classified: tuple[dict[str, Any], ...],
    *,
    excluded: frozenset[str] = frozenset(),
) -> dict[str, str]:
    return {
        param: f"REGISTERED_PENDING({route})"
        for param, route in sim_deferred_map(classified).items()
        if param not in excluded
    }


HISTORICAL_STRUCTURAL_POSTULATES = (
    (
        "postulate_1",
        "imposed w-hat axis + w=0 surface (conceded-wall)",
    ),
    (
        "postulate_2",
        "u^a same-medium surface collective, tangentially free-slip",
    ),
    (
        "postulate_3",
        "T0 P^i reused as the Cosserat micro-rotation reservoir",
    ),
    (
        "postulate_4",
        "baseline P^i spin-wave status = massless",
    ),
    (
        "postulate_5",
        "w-hat-dependent parity-even P-u operator structural cost",
    ),
    (
        "postulate_6",
        "no C5 phi analog / no longitudinal constraint",
    ),
)
RETIRED_DRIFT_KEYS = frozenset(
    {"lambda_Pu", "postulate_3", "postulate_4", "postulate_5"}
)
EXPECTED_OPERATIVE_KEYS = frozenset(
    {
        "rho_br",
        "mu_R",
        "Omega_w",
        "g_ell",
        "postulate_1",
        "postulate_2",
        "postulate_6",
    }
)


def freeze_history(*, rename_survivor: bool = False) -> tuple[DriftEntry, ...]:
    table = (
        DriftEntry("rho_br", "rho_br"),
        DriftEntry("mu_R", "mu_R"),
        DriftEntry("lambda_Pu", "lambda_Pu"),
        DriftEntry("Omega_w", "Omega_w"),
        DriftEntry("g_ell", "g_ell(w)"),
        *(
            DriftEntry(key, name)
            for key, name in HISTORICAL_STRUCTURAL_POSTULATES
        ),
    )
    if rename_survivor:
        table = tuple(
            replace(entry, name="rho_B0")
            if entry.key == "rho_br"
            else entry
            for entry in table
        )
    return table


def anti_absorption_state(*, rename_survivor: bool = False) -> dict[str, Any]:
    historical = freeze_history(rename_survivor=rename_survivor)
    historical_by_key = {entry.key: entry for entry in historical}
    operative_keys = frozenset(historical_by_key) - RETIRED_DRIFT_KEYS
    operative = tuple(
        entry for entry in historical if entry.key in operative_keys
    )
    operative_names = frozenset(entry.name for entry in operative)
    return {
        "operative_keys": operative_keys,
        "operative_count": len(operative),
        "operative_names": operative_names,
        "disjoint": set(EXPECTED_TRIO).isdisjoint(operative_names),
    }


def current_nonentailment_witness() -> tuple[bool, sp.Expr]:
    values = {
        "rho": sp.Integer(1),
        "rho_br": sp.Integer(1),
        "rho_B0": sp.Integer(2),
        "chi_c": sp.Integer(4),
        "mu_R": sp.Integer(2),
        "K": sp.Integer(1),
        "m": sp.Integer(5),
        "M_h": sp.Integer(1),
        "c_E": sp.Integer(3),
        "C_hu": sp.Rational(1, 2),
        "K_h": sp.Integer(9),
        "B_eff": sp.Integer(1),
        "sigma": sp.Rational(35, 4),
    }
    base_residuals = (
        values["B_eff"] * values["chi_c"] - values["rho_B0"] ** 2,
        values["K_h"] - values["M_h"] * values["c_E"] ** 2,
        values["B_eff"] * values["K_h"]
        - values["C_hu"] ** 2
        - values["sigma"],
    )
    positive = all(bool(value > 0) for value in values.values())
    lock_b = values["c_E"] ** 2 * values["rho_br"] - values["mu_R"]
    return all(sp.simplify(item) == 0 for item in base_residuals) and positive, lock_b


def freedom_tie_status(*, include_locks: bool) -> tuple[str, Any]:
    B, Kh, sigma = sp.symbols("B Kh sigma", positive=True)
    if include_locks:
        q_from_charge_and_locks = sp.expand(2 * B * Kh)
        c_squared_from_tie = q_from_charge_and_locks
        stability_residual = sp.factor(
            B * Kh - c_squared_from_tie - sigma
        )
        expected_negative = sp.factor(-(B * Kh + sigma))
        status = (
            "UNSAT"
            if stability_residual == expected_negative
            and bool(B.is_positive)
            and bool(Kh.is_positive)
            and bool(sigma.is_positive)
            else "UNRESOLVED"
        )
        return status, stability_residual
    witness = {
        "rho_br": sp.Integer(1),
        "mu_R": sp.Integer(1),
        "M_h": sp.Integer(1),
        "c_E": sp.Integer(2),
        "K_h": sp.Integer(4),
        "B_eff": sp.Integer(1),
        "C_hu": sp.sqrt(2),
        "sigma": sp.Integer(2),
        "q_h_sq": sp.Integer(2),
    }
    residuals = (
        witness["K_h"] - witness["M_h"] * witness["c_E"] ** 2,
        witness["B_eff"] * witness["K_h"]
        - witness["C_hu"] ** 2
        - witness["sigma"],
        witness["C_hu"] ** 2 - witness["q_h_sq"],
        witness["q_h_sq"] * witness["rho_br"]
        - 2
        * witness["B_eff"]
        * witness["M_h"]
        * witness["mu_R"],
    )
    positive = all(bool(sp.N(value) > 0) for value in witness.values())
    return (
        "SAT"
        if all(sp.simplify(item) == 0 for item in residuals) and positive
        else "UNRESOLVED",
        residuals,
    )


def dimension_state(*, corrupt: bool = False) -> tuple[tuple[int, int, int], ...]:
    dim_mass = sp.Matrix([1, 0, 0])
    dim_bulk_density = sp.Matrix([0, -4, 0])
    dim_layer = sp.Matrix([0, 1, 0])
    dim_rho_br_source = sp.Matrix([1, -3, 0])
    dim_speed = sp.Matrix([0, 2 if corrupt else 1, -1])
    varrho = dim_mass + dim_bulk_density + dim_layer
    cone_product = 2 * dim_speed + dim_rho_br_source
    return (
        tuple(int(value) for value in varrho),
        tuple(int(value) for value in dim_rho_br_source),
        tuple(int(value) for value in cone_product),
    )


EXPECTED_FULL_ORIGIN = {
    "rho": "BASE_SUBSTRATE",
    "K": "BASE_SUBSTRATE",
    "m": "BASE_SUBSTRATE",
    "a": "BASE_SUBSTRATE",
    "ell_g": "ACCEPTED_GEOMETRY_SUBSTRATE",
    "g_ell(w)": "ACCEPTED_PROFILE_GIVEN_ell_g",
    "varrho_br[rho]": "OUT_OF_ACTIVE_NG5",
    "Sigma_n[rho]": "OUT_OF_ACTIVE_NG5",
    "delta_Sigma[rho]": "OUT_OF_ACTIVE_NG5",
    "rho_br": "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED",
    "mu_R": "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED",
    "c_E": "REDUCIBLE_IN_PRINCIPLE_SIM_DEFERRED",
    "c_gamma": "DEPENDENT",
    "rho_B0": "IRREDUCIBLY_INDEPENDENT",
    "chi_c": "IRREDUCIBLY_INDEPENDENT",
    "B_eff": "DEPENDENT_ON_IRREDUCIBLE",
    "C_hu": "IRREDUCIBLY_INDEPENDENT",
    "Q_E": "CALIBRATED_ANCHOR",
    "ell": "CALIBRATED_GEOMETRY_INPUT",
    "b": "CALIBRATED_GEOMETRY_INPUT",
    "M_h": "CALIBRATED_GEOMETRY_INPUT",
    "K_h": "DEPENDENT",
    "q_h": "DEPENDENT",
    "c_L1": "OUT_OF_ACTIVE_NG5",
    "c_L2": "OUT_OF_ACTIVE_NG5",
    "A_R": "OUT_OF_ACTIVE_NG5",
    "k_R": "OUT_OF_ACTIVE_NG5",
    "lambda_Cdiv": "OUT_OF_ACTIVE_NG5",
    "chi_Cpin": "OUT_OF_ACTIVE_NG5",
    "J_Pu": "OUT_OF_ACTIVE_NG5",
    "kappa_Pu": "OUT_OF_ACTIVE_NG5",
    "lambda_Pu": "OUT_OF_ACTIVE_NG5",
    "Omega_w": "OUT_OF_ACTIVE_NG5",
    "lambda_N": "OUT_OF_ACTIVE_NG5",
    "lambda_tau": "OUT_OF_ACTIVE_NG5",
    "Nu": "OUT_OF_ACTIVE_NG5",
    "a_T": "OUT_OF_ACTIVE_NG5",
    "a_Tp": "OUT_OF_ACTIVE_NG5",
    "a_L": "OUT_OF_ACTIVE_NG5",
}

EXPECTED_ROUTE_LEDGER = {
    "rho_br": {
        "route_name": "Route-A",
        "valid": True,
        "reason": "VALID_REGISTERED_ROUTE",
    },
    "mu_R": {
        "route_name": "Route-A",
        "valid": True,
        "reason": "VALID_REGISTERED_ROUTE",
    },
    "c_E": {
        "route_name": "Route-A",
        "valid": True,
        "reason": "VALID_REGISTERED_ROUTE",
    },
    "c_gamma": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
    "rho_B0": {
        "route_name": "Future-Compression-4D-to-3D",
        "valid": False,
        "reason": "REJECTED_RESULT_STATUS_PROMISSORY_ONLY",
    },
    "chi_c": {
        "route_name": "Future-Compression-4D-to-3D",
        "valid": False,
        "reason": "REJECTED_RESULT_STATUS_PROMISSORY_ONLY",
    },
    "B_eff": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
    "C_hu": {
        "route_name": "Future-Embedding-Overlap",
        "valid": False,
        "reason": "REJECTED_RESULT_STATUS_PROMISSORY_ONLY",
    },
    "Q_E": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
    "ell": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
    "b": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
    "M_h": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
    "K_h": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
    "q_h": {
        "route_name": None,
        "valid": False,
        "reason": "NO_REGISTERED_ROUTE_FOR_PARAM",
    },
}


def origin_map(classified: Iterable[dict[str, Any]]) -> dict[str, str]:
    return {
        record["row"].param: record["origin"] for record in classified
    }


def route_ledger(
    classified: Iterable[dict[str, Any]]
) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for record in classified:
        if not record["row"].active:
            continue
        evaluation = record["route_evaluation"]
        result[record["row"].param] = {
            "route_name": evaluation["route_name"],
            "valid": evaluation["valid"],
            "reason": evaluation["reason"],
        }
    return result


def recorded_state(
    classified: Iterable[dict[str, Any]], *, drop_param: str | None = None
) -> tuple[frozenset[str], frozenset[str]]:
    active = frozenset(
        record["row"].param
        for record in classified
        if record["row"].active
    )
    recorded = frozenset(
        record["row"].param
        for record in classified
        if record["row"].active
        and record["row"].param != drop_param
        and record["route_evaluation"] is not None
    )
    return active, recorded


def route_conjunct_control(
    conjunct: str, *, neutralize: bool
) -> dict[str, Any]:
    routes = copied_routes()
    if not neutralize:
        fact = routes["Route-A"]
        if conjunct == "required_fields_present":
            del fact["source"]
        elif conjunct == "named_solve_in_provenance":
            fact["named_solve_in_provenance"] = False
        elif conjunct == "result_status_allowed":
            fact["result_status"] = "BY_TUNING"
        elif conjunct == "finite_listed_missing_objects":
            fact["missing_objects"] = "not-a-list"
        elif conjunct == "target_identity":
            fact["targets"] = ["rho_br", "mu_R"]
        elif conjunct == "joint_target_identity":
            fact["required_joint_targets"] = [
                "rho_br",
                "mu_R",
                "joint_sentinel",
            ]
        elif conjunct == "target_blind":
            fact["target_blind"] = False
        elif conjunct == "falsifiers_clear":
            fact["falsifiers"] = ["synthetic-falsifier"]
        else:
            raise KeyError(conjunct)
    classified = run_case(routes=routes)
    records = by_param(classified)
    covered = set(EXPECTED_ROUTE_A_TARGETS)
    expected_flips = {"c_E"} if conjunct == "target_identity" else covered
    actual_flips = {
        param
        for param in covered
        if records[param]["origin"] == "IRREDUCIBLY_INDEPENDENT"
    }
    isolated = True
    for param in covered:
        checks = records[param]["route_evaluation"]["checks"]
        false_named = {
            name for name in ROUTE_CONJUNCTS if not checks[name]
        }
        expected_false = (
            {conjunct} if param in expected_flips else set()
        )
        isolated = isolated and false_named == expected_false
    return {
        "fired": actual_flips == expected_flips and isolated,
        "actual_flips": actual_flips,
        "isolated": isolated,
    }


SOURCE_PREDICATE_TOTAL = 37
SOURCE_PREDICATE_UNIVERSE = (
    "pathA41.py::dim_derivations",
    "pathA41.py::registered_route_for",
    "pathA41.py::classify_origin",
    "pathA41.py::classify_rows",
    "pathA41.py::decide_verdict",
    "pathA41.py::location_closure",
    "pathA41.py::lineage_adjudication",
    "pathA41.py::freedom_states",
    "pathA41.py::interpretation.reduction_status",
    "pathA41.py::interpretation.physical_meaning",
    "pathA41.py::production_location_closure_raise",
    "pathA41.py::p40_freedom_tie_no_go",
    "pathA41.py::p40_current_nonentailment",
    "pathA41.py::control.AB_delete_registry",
    "pathA41.py::control.route_blank_Route_A",
    "pathA41.py::control.route_field_blank_Route_A_target_blind",
    "pathA41.py::control.route_field_blank_Route_A_missing_objects",
    "pathA41.py::control.calibration_ablation_Q_E",
    "pathA41.py::control.irreducible_synthetic",
    "pathA41.py::control.reducible_derived_synthetic",
    "pathA41.py::control.location_closure_out_of_arena",
    "pathA41.py::control.contradiction",
    "pathA41.py::control.residual_multiplier_ablation",
    "pathA41.py::control.route_eval_recorded_for_all_active_rows",
    "pathA41.wl::hardcoded_production",
    "pathA41.wl::literal_controls",
    "pathA41.py::argparse_compare_harness",
    "pathA41.py::file_writing_and_main",
    "pathA41.py::build_results_yaml_write_report",
    "pathA41.py::comparison_payload_digest_sha_count",
    "pathA41.py::compare_with_mathematica",
    "pathA41.py::filesystem_source_token_scans",
    "pathA41.wl::filesystem_assertContains_scans",
    "pathA41.py::importlib_pathA40",
    "pathA41.py::extraction_audit_narration",
    "pathA41.py::comparison_map_json_serialization",
    "pathA41.wl::Import_canon_cross_read",
)

PRESERVED_SOURCE_IDS = SOURCE_PREDICATE_UNIVERSE[:7] + (
    SOURCE_PREDICATE_UNIVERSE[10],
    SOURCE_PREDICATE_UNIVERSE[11],
    SOURCE_PREDICATE_UNIVERSE[12],
    *SOURCE_PREDICATE_UNIVERSE[13:24],
)
REPLACED_SOURCE_IDS = (
    SOURCE_PREDICATE_UNIVERSE[7],
    SOURCE_PREDICATE_UNIVERSE[8],
    SOURCE_PREDICATE_UNIVERSE[9],
    SOURCE_PREDICATE_UNIVERSE[24],
    SOURCE_PREDICATE_UNIVERSE[25],
)
SCOPED_SOURCE_IDS = SOURCE_PREDICATE_UNIVERSE[26:]
NEW_STAGE_IDS = (
    "stage041::ANTI_ABSORPTION_DRIFT_SIDE",
    "stage041::DUAL_ENGINE_TERMS",
    "stage041::SOURCE_TO_STAGE_MANIFEST",
    *(
        "stage041::" + tooth
        for tooth in CONJUNCT_TEETH
    ),
)


def manifest_rows() -> tuple[tuple[str, str], ...]:
    return (
        *((identifier, "preserved-folded") for identifier in PRESERVED_SOURCE_IDS),
        *((identifier, "replaced-by-stronger") for identifier in REPLACED_SOURCE_IDS),
        *((identifier, "scoped-out") for identifier in SCOPED_SOURCE_IDS),
        *((identifier, "newly-added") for identifier in NEW_STAGE_IDS),
    )


EXPECTED_MANIFEST_CATEGORY = dict(manifest_rows())


def manifest_guard(*, mutate: bool = False) -> dict[str, Any]:
    rows = list(manifest_rows())
    if mutate:
        dropped = SCOPED_SOURCE_IDS[0]
        rows = [row for row in rows if row[0] != dropped]
        moved = PRESERVED_SOURCE_IDS[0]
        rows = [
            (identifier, "replaced-by-stronger")
            if identifier == moved
            else (identifier, category)
            for identifier, category in rows
        ]
    identifiers = [identifier for identifier, _ in rows]
    source_rows = [
        row for row in rows if row[0] in SOURCE_PREDICATE_UNIVERSE
    ]
    source_identifiers = [row[0] for row in source_rows]
    categories = {
        category: {identifier for identifier, item_category in rows if item_category == category}
        for category in (
            "preserved-folded",
            "replaced-by-stronger",
            "newly-added",
            "scoped-out",
        )
    }
    category_names = tuple(categories)
    disjoint = all(
        categories[left].isdisjoint(categories[right])
        for index, left in enumerate(category_names)
        for right in category_names[index + 1:]
    )
    source_counts = Counter(source_identifiers)
    engine_prefixes = {
        identifier.split("::", 1)[0] for identifier in identifiers
    }
    return {
        "ok": (
            len(source_rows) == SOURCE_PREDICATE_TOTAL
            and set(source_identifiers) == set(SOURCE_PREDICATE_UNIVERSE)
            and not any(count > 1 for count in source_counts.values())
            and set(identifier for identifier in identifiers if identifier in NEW_STAGE_IDS)
            == set(NEW_STAGE_IDS)
            and disjoint
            and dict(rows) == EXPECTED_MANIFEST_CATEGORY
            and engine_prefixes
            == {"pathA41.py", "pathA41.wl", "stage041"}
        ),
        "source_count": len(source_rows),
        "partition": Counter(category for _, category in rows),
        "disjoint": disjoint,
    }


def run_assertions() -> str:
    production_rows = build_rows()
    production = run_case(rows=production_rows)
    production_by_param = by_param(production)
    production_irreducible = active_irreducible(production)
    production_result = decide_verdict(production)

    section("Computed production origin and RouteEvaluation ledgers")
    origin_case = production
    if ACTIVE_MUTATION == "ORIGIN_CLASSIFICATION":
        mutated_routes = copied_routes()
        mutated_routes["Future-Embedding-Overlap"].update(
            {
                "named_solve_in_provenance": True,
                "result_status": "REGISTERED_DEFERRED",
            }
        )
        origin_case = run_case(routes=mutated_routes)
    expect_bool(
        "ORIGIN_CLASSIFICATION",
        active_irreducible(origin_case) == EXPECTED_TRIO,
        {"computed": active_irreducible(origin_case), "ratified": EXPECTED_TRIO},
    )

    full_origin_case = (
        run_case(rows=build_rows(b_eff_dependency_mutation=True))
        if ACTIVE_MUTATION == "FULL_ORIGIN_LEDGER"
        else production
    )
    expect_bool(
        "FULL_ORIGIN_LEDGER",
        origin_map(full_origin_case) == EXPECTED_FULL_ORIGIN,
        origin_map(full_origin_case),
    )

    route_case = (
        run_case(rows=build_rows(c_gamma_route_mutation=True))
        if ACTIVE_MUTATION == "FULL_ROUTE_EVALUATION_LEDGER"
        else production
    )
    expect_bool(
        "FULL_ROUTE_EVALUATION_LEDGER",
        route_ledger(route_case) == EXPECTED_ROUTE_LEDGER,
        route_ledger(route_case),
    )

    section("Production verdict, count, freedom, locations, and lineage")
    verdict_case = production
    if ACTIVE_MUTATION == "PRODUCTION_VERDICT":
        routes = copied_routes()
        routes["Future-Embedding-Overlap"].update(
            {
                "named_solve_in_provenance": True,
                "result_status": "REGISTERED_DEFERRED",
            }
        )
        verdict_case = run_case(routes=routes)
    verdict_actual = decide_verdict(verdict_case)["verdict"]
    expect_bool(
        "PRODUCTION_VERDICT",
        verdict_actual == EXPECTED_VERDICT,
        verdict_actual,
    )

    drift_case = (
        run_case(rows=build_rows(add_xi=True))
        if ACTIVE_MUTATION == "DRIFT_COUNT"
        else production
    )
    drift_count = len(active_irreducible(drift_case))
    expect_bool("DRIFT_COUNT", drift_count == 3, drift_count)

    rho_freedom_case = (
        run_case(routes={})
        if ACTIVE_MUTATION == "FREEDOM_STATE_RHO_BR"
        else production
    )
    rho_freedom = freedom_state(by_param(rho_freedom_case)["rho_br"])
    expect_bool(
        "FREEDOM_STATE_RHO_BR",
        rho_freedom == "FREEDOM_SIM_DEFERRED{Route-A}",
        rho_freedom,
    )

    c_hu_freedom_case = production
    if ACTIVE_MUTATION == "FREEDOM_STATE_C_HU":
        routes = copied_routes()
        routes["Future-Embedding-Overlap"].update(
            {
                "named_solve_in_provenance": True,
                "result_status": "REGISTERED_DEFERRED",
            }
        )
        c_hu_freedom_case = run_case(routes=routes)
    c_hu_freedom = freedom_state(by_param(c_hu_freedom_case)["C_hu"])
    expect_bool(
        "FREEDOM_STATE_C_HU",
        c_hu_freedom
        == "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}",
        c_hu_freedom,
    )

    closure_rows = (
        build_rows(add_location_sentinel=True)
        if ACTIVE_MUTATION == "LOCATION_CLOSURE"
        else production_rows
    )
    closure = location_closure(closure_rows)
    expect_bool(
        "LOCATION_CLOSURE",
        closure["no_fourth_arena"] is True,
        closure,
    )

    label_rows = (
        build_rows(c_hu_arena_mutation=True)
        if ACTIVE_MUTATION == "ARENA_LABELS"
        else production_rows
    )
    labels = arena_labels(label_rows)
    expected_labels = {
        "rho_B0": "3D brane surface",
        "chi_c": "3D brane surface",
        "C_hu": "throat/embedding seam",
    }
    expect_bool("ARENA_LABELS", labels == expected_labels, labels)

    lineage_facts = (
        residual_lineage_facts("Xi_residual")
        if ACTIVE_MUTATION == "LINEAGE_FINDING"
        else production_lineage_facts()
    )
    lineage = lineage_from_ledger(production, lineage_facts)
    q_e_isolation_case = run_case(anchors=frozenset())
    q_e_isolation = lineage_from_ledger(
        q_e_isolation_case, production_lineage_facts()
    )
    expect_bool(
        "LINEAGE_FINDING",
        lineage == ("DIFFERENT", "NO_OVERCOUNT_ROUTE_A_PENDING")
        and q_e_isolation
        == ("DIFFERENT", "NO_OVERCOUNT_ROUTE_A_PENDING"),
        {"production": lineage, "Q_E_ablation_isolation": q_e_isolation},
    )

    pending = reduction_status(
        production,
        excluded=frozenset({"c_E"})
        if ACTIVE_MUTATION == "REDUCTION_STATUS"
        else frozenset(),
    )
    sim_keys = frozenset(sim_deferred_map(production))
    expected_pending = frozenset(EXPECTED_ROUTE_A_TARGETS)
    expect_bool(
        "REDUCTION_STATUS",
        frozenset(pending) == expected_pending
        and sim_keys == expected_pending,
        {"pending": pending, "sim_deferred_keys": sim_keys},
    )

    anti = anti_absorption_state(
        rename_survivor=ACTIVE_MUTATION == "ANTI_ABSORPTION_DRIFT_SIDE"
    )
    expect_bool(
        "ANTI_ABSORPTION_DRIFT_SIDE",
        anti["operative_keys"] == EXPECTED_OPERATIVE_KEYS
        and anti["operative_count"] == 7
        and anti["disjoint"] is True,
        anti,
    )

    section("Individually asserted controls")
    blank_all = (
        production
        if ACTIVE_MUTATION == "CTRL_AB_DELETE_REGISTRY"
        else run_case(routes={})
    )
    blank_all_irred = active_irreducible(blank_all)
    expect_bool(
        "CTRL_AB_DELETE_REGISTRY",
        set(production_irreducible) < set(blank_all_irred)
        and blank_all_irred
        == ("rho_br", "mu_R", "c_E", "rho_B0", "chi_c", "C_hu"),
        blank_all_irred,
    )

    no_route_a_routes = copied_routes()
    if ACTIVE_MUTATION != "CTRL_ROUTE_BLANK_ROUTE_A":
        del no_route_a_routes["Route-A"]
    no_route_a = run_case(routes=no_route_a_routes)
    no_route_a_records = by_param(no_route_a)
    expect_bool(
        "CTRL_ROUTE_BLANK_ROUTE_A",
        all(
            no_route_a_records[param]["origin"]
            == "IRREDUCIBLY_INDEPENDENT"
            for param in EXPECTED_ROUTE_A_TARGETS
        ),
        {param: no_route_a_records[param]["origin"] for param in EXPECTED_ROUTE_A_TARGETS},
    )

    for tooth, field in (
        ("CTRL_ROUTE_FIELD_BLANK_TARGET_BLIND", "target_blind"),
        ("CTRL_ROUTE_FIELD_BLANK_MISSING_OBJECTS", "missing_objects"),
    ):
        routes = copied_routes()
        if ACTIVE_MUTATION != tooth:
            del routes["Route-A"][field]
        records = by_param(run_case(routes=routes))
        fired = all(
            records[param]["origin"] == "IRREDUCIBLY_INDEPENDENT"
            and records[param]["route_evaluation"]["reason"]
            == "MISSING_REQUIRED_ROUTE_FIELD"
            for param in EXPECTED_ROUTE_A_TARGETS
        )
        expect_bool(
            tooth,
            fired,
            {
                param: records[param]["route_evaluation"]
                for param in EXPECTED_ROUTE_A_TARGETS
            },
        )

    conjunct_map = dict(zip(CONJUNCT_TEETH, ROUTE_CONJUNCTS))
    for tooth in CONJUNCT_TEETH:
        control = route_conjunct_control(
            conjunct_map[tooth],
            neutralize=ACTIVE_MUTATION == tooth,
        )
        expect_bool(tooth, control["fired"], control)

    q_e_case = run_case(
        anchors=frozenset({"Q_E"})
        if ACTIVE_MUTATION == "CTRL_CALIBRATION_ABLATION_Q_E"
        else frozenset()
    )
    q_e_irred = active_irreducible(q_e_case)
    q_e_record = by_param(q_e_case)["Q_E"]
    expect_bool(
        "CTRL_CALIBRATION_ABLATION_Q_E",
        q_e_record["origin"] == "IRREDUCIBLY_INDEPENDENT"
        and q_e_irred == ("rho_B0", "chi_c", "C_hu", "Q_E")
        and len(q_e_irred) == 4,
        {"origin": q_e_record["origin"], "active_irreducible": q_e_irred},
    )

    xi_case = run_case(
        rows=build_rows(
            add_xi=ACTIVE_MUTATION != "CTRL_IRREDUCIBLE_SYNTHETIC"
        )
    )
    xi_irred = active_irreducible(xi_case)
    expect_bool(
        "CTRL_IRREDUCIBLE_SYNTHETIC",
        "xi_active" in xi_irred
        and by_param(xi_case)["xi_active"]["route_evaluation"]["reason"]
        == "NO_REGISTERED_ROUTE_FOR_PARAM",
        xi_irred,
    )

    p_syn_case = run_case(
        rows=build_rows(
            add_p_syn=True,
            p_syn_status_mutation=(
                ACTIVE_MUTATION == "CTRL_REDUCIBLE_DERIVED_SYNTHETIC"
            ),
        )
    )
    p_syn_records = by_param(p_syn_case)
    expect_bool(
        "CTRL_REDUCIBLE_DERIVED_SYNTHETIC",
        p_syn_records["p_syn"]["origin"] == "REDUCIBLE_DERIVED"
        and active_irreducible(p_syn_case) == EXPECTED_TRIO,
        {
            "origin": p_syn_records["p_syn"]["origin"],
            "active_irreducible": active_irreducible(p_syn_case),
        },
    )

    witness_ok, lock_b = current_nonentailment_witness()
    tie_status, tie_evidence = freedom_tie_status(
        include_locks=ACTIVE_MUTATION != "CTRL_CONTRADICTION"
    )
    contradiction_verdict = decide_verdict(
        production, no_go=tie_status == "UNSAT"
    )["verdict"]
    expect_bool(
        "CTRL_CONTRADICTION",
        witness_ok
        and lock_b == 7
        and tie_status == "UNSAT"
        and contradiction_verdict == "NO_GO(cone-lock-feedback)",
        {
            "lock_B_witness": lock_b,
            "tie_status": tie_status,
            "certificate_or_residuals": tie_evidence,
            "verdict": contradiction_verdict,
        },
    )

    active_records, recorded_records = recorded_state(
        production,
        drop_param="c_gamma"
        if ACTIVE_MUTATION == "CTRL_ROUTE_EVAL_RECORDED"
        else None,
    )
    expect_bool(
        "CTRL_ROUTE_EVAL_RECORDED",
        active_records == recorded_records,
        {"active": active_records, "recorded": recorded_records},
    )

    closed_lineage = lineage_adjudication(residual_lineage_facts("1"))
    free_lineage = lineage_adjudication(
        residual_lineage_facts(
            "1"
            if ACTIVE_MUTATION == "CTRL_RESIDUAL_MULTIPLIER_ABLATION"
            else "Xi_residual"
        )
    )
    expect_bool(
        "CTRL_RESIDUAL_MULTIPLIER_ABLATION",
        closed_lineage == ("SAME", "OVERCOUNT_OR_SMUGGLE_CONTROL")
        and free_lineage
        == (
            "DIFFERENT",
            "DIFFERENT_BRANE_INERTIA_OBJECTS_RESIDUAL_MULTIPLIER",
        ),
        {"closed": closed_lineage, "free": free_lineage},
    )

    control_rows = build_rows(
        add_location_sentinel=(
            ACTIVE_MUTATION != "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA"
        )
    )
    production_closure = location_closure(production_rows)
    control_closure = location_closure(control_rows)
    location_control_fired = (
        production_closure["no_fourth_arena"] is True
        and control_closure["no_fourth_arena"] is False
        and bool(control_closure["offending"])
    )
    expect_bool(
        "CTRL_LOCATION_CLOSURE_OUT_OF_ARENA",
        location_control_fired,
        {"production": production_closure, "control": control_closure},
    )

    section("Build-global dimensions, local inventory, and verdict table")
    dimensions = dimension_state(
        corrupt=ACTIVE_MUTATION == "DIMENSION_HOMOGENEITY"
    )
    expected_dimensions = (
        (1, -3, 0),
        (1, -3, 0),
        (1, -1, -2),
    )
    expect_bool(
        "DIMENSION_HOMOGENEITY",
        dimensions == expected_dimensions,
        dimensions,
    )

    anti_production = anti_absorption_state()
    production_lineage = lineage_from_ledger(
        production, production_lineage_facts()
    )[1]
    local_inventory = {
        "active_irreducible": production_irreducible,
        "drift_count": production_result["drift_count"],
        "freedom_rho_br": freedom_state(production_by_param["rho_br"]),
        "freedom_C_hu": freedom_state(production_by_param["C_hu"]),
        "no_fourth_arena": production_closure["no_fourth_arena"],
        "lineage": production_lineage,
        "anti_absorption": anti_production["disjoint"],
        "dimension_state": dimension_state(),
        "lock_B_witness": lock_b,
    }
    if ACTIVE_MUTATION == "DUAL_ENGINE_TERMS":
        del local_inventory["lock_B_witness"]
    expected_inventory = {
        "active_irreducible": EXPECTED_TRIO,
        "drift_count": 3,
        "freedom_rho_br": "FREEDOM_SIM_DEFERRED{Route-A}",
        "freedom_C_hu": "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}",
        "no_fourth_arena": True,
        "lineage": "NO_OVERCOUNT_ROUTE_A_PENDING",
        "anti_absorption": True,
        "dimension_state": expected_dimensions,
        "lock_B_witness": sp.Integer(7),
    }
    expect_bool(
        "DUAL_ENGINE_TERMS",
        local_inventory == expected_inventory,
        local_inventory,
    )

    blank_verdict = decide_verdict(run_case(routes={}))["verdict"]
    q_e_verdict = decide_verdict(run_case(anchors=frozenset()))["verdict"]
    tie_table_status, _ = freedom_tie_status(include_locks=True)
    no_go_verdict = decide_verdict(
        production, no_go=tie_table_status == "UNSAT"
    )["verdict"]
    conditional_routes = valid_future_routes(future_status="SOLVED_PASS")
    if ACTIVE_MUTATION == "VERDICT_REDERIVATION":
        conditional_routes["Future-Embedding-Overlap"][
            "named_solve_in_provenance"
        ] = False
    conditional_verdict = decide_verdict(
        run_case(routes=conditional_routes)
    )["verdict"]
    consistent_verdict = decide_verdict(
        run_case(
            routes=valid_future_routes(
                future_status="SOLVED_PASS",
                route_a_status="SOLVED_PASS",
            )
        )
    )["verdict"]
    actual_verdicts = (
        production_result["verdict"],
        blank_verdict,
        q_e_verdict,
        no_go_verdict,
        conditional_verdict,
        consistent_verdict,
    )
    expected_conditional = (
        "ONE_MEDIUM_CONDITIONAL(sim-deferred={"
        "rho_br->Route-A,mu_R->Route-A,c_E->Route-A"
        "}; calibrated={"
        "Q_E->CALIBRATED_ANCHOR,"
        "ell->CALIBRATED_GEOMETRY_INPUT,"
        "b->CALIBRATED_GEOMETRY_INPUT,"
        "M_h->CALIBRATED_GEOMETRY_INPUT})"
    )
    expected_verdicts = (
        EXPECTED_VERDICT,
        "SECOND_MEDIUM_DRIFT(active_irreducible={rho_br,mu_R,c_E,rho_B0,chi_c,C_hu})",
        "SECOND_MEDIUM_DRIFT(active_irreducible={rho_B0,chi_c,C_hu,Q_E})",
        "NO_GO(cone-lock-feedback)",
        expected_conditional,
        "ONE_MEDIUM_CONSISTENT",
    )
    expect_bool(
        "VERDICT_REDERIVATION",
        actual_verdicts == expected_verdicts,
        {"actual": actual_verdicts, "expected": expected_verdicts},
    )

    section("Bounded source-to-stage predicate manifest")
    manifest = manifest_guard(
        mutate=ACTIVE_MUTATION == "SOURCE_TO_STAGE_MANIFEST"
    )
    expect_bool("SOURCE_TO_STAGE_MANIFEST", manifest["ok"], manifest)

    arena_groups: dict[str, list[str]] = {}
    for param, location in labels.items():
        arena_groups.setdefault(location, []).append(param)
    physical_meaning = "; ".join(
        f"{location}: {','.join(params)}"
        for location, params in arena_groups.items()
    )
    print("")
    print(f"VERDICT={production_result['verdict']}")
    print(
        "ACTIVE_IRREDUCIBLE={"
        + ",".join(production_irreducible)
        + "}"
    )
    print(f"DRIFT_COUNT={production_result['drift_count']}")
    print(
        "FREEDOM_STATES="
        + freedom_state(production_by_param["rho_br"])
        + ";"
        + freedom_state(production_by_param["C_hu"])
    )
    print(f"NO_FOURTH_ARENA={production_closure['no_fourth_arena']}")
    print(f"LINEAGE={production_lineage}")
    print("POST_D16_DRIFT_TOKEN=POST_D16_DRIFT(7)")
    print("ANTI_ABSORPTION_DRIFT_SIDE=DISJOINT")
    print("PHYSICAL_MEANING_FROM_ARENA_MAP=" + physical_meaning)
    print(
        "Q_E_NEAR_MISS=without calibration anchor, "
        "active_irreducible={rho_B0,chi_c,C_hu,Q_E}; drift_count=4"
    )
    print(
        "HONEST_CAVEAT=bookkeeping/labelling closure only; dynamical "
        "brane-genesis remains open and is handed to the deferred "
        "throat/embedding solve"
    )
    print("FRAMING=CHARACTERIZED-DEPARTURE; no_fourth_arena is labelling-level")
    print(f"SOURCE_PREDICATE_TOTAL={SOURCE_PREDICATE_TOTAL}")
    print("EXECUTABLE_TOOTH_TOTAL=35")
    return production_result["verdict"]


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage041_ng5_second_medium_reducibility SymPy audit")
    print(
        "ROUTE=Row dataclasses + source-order classifier cascade + "
        "runtime RouteEvaluation conjuncts"
    )
    print("FILE_IO=none; CROSS_ENGINE_COMPARE=none")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")

    verdict = run_assertions()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print("")
    print(f"TOOTH_COUNT={len(TOOTH_ORDER)}")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print(f"OVERALL PASS: SymPy independently reached {verdict}")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(
            "OVERALL FAIL: SymPy stage041 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
