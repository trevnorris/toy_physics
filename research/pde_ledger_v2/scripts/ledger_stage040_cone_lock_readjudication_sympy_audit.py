#!/usr/bin/env python3
"""Ledger stage040 SymPy audit: cone-lock re-adjudication.

Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.  The
production token remains ``CONE_LOCK_CALIBRATED``; its adjudicated meaning is
CALIBRATED / UNCOMMITTED.  Both cone locks are available calibration choices,
not earned equalities.

The algebraic route is intentionally SymPy-native: exact positive witnesses,
Groebner remainder ideal-membership tests, and grevlex initial-monomial
Krull dimension.  Tooth-local ablation uses ``LEDGER_STAGE040_MUTATION``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import itertools
import os
from typing import Any, Iterable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE040_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

TOOTH_ORDER = (
    "ROUTE_A_GRADE",
    "ROUTE_B_STATUS",
    "FREEDOM_STATUS",
    "SECTOR_SAT",
    "LOCK_SAT",
    "PROVENANCE_LOCK_A",
    "LOCK_A_WITNESS_VALUE",
    "PROVENANCE_LOCK_B",
    "LOCK_B_WITNESS_VALUE",
    "DELTA_R",
    "PRODUCTION_VERDICT",
    "FIELD_OVERLAY_DET_ON_CONE",
    "FIELD_OVERLAY_POLES",
    "OPEN_110_CARRY",
    "EARNED_VS_CALIBRATED_PARTITION",
    "CONDITIONALITY_FREEDOM",
    "CTRL_WELL_POSED",
    "CTRL_ABSENT",
    "CTRL_PARTIAL_INVENTORY",
    "CTRL_FORCED_LOCK",
    "CTRL_A_ONLY_PARTIAL",
    "CTRL_B_ONLY_PARTIAL",
    "CTRL_OVER_CONSTRAINED",
    "CTRL_FREEDOM_TIE",
    "VERDICT_REDERIVATION",
    "DIMENSION_HOMOGENEITY",
    "DUAL_ENGINE_TERMS",
    "SOURCE_TO_STAGE_MANIFEST",
)

ABLATION_DESCRIPTIONS = {
    "ROUTE_A_GRADE": "mark R1 present in the private production inventory",
    "ROUTE_B_STATUS": "inject the rejected thin-plate over-import relation",
    "FREEDOM_STATUS": "inject a freedom-tie relation into the private freedom inventory",
    "SECTOR_SAT": "inject the positive over-stability contradiction",
    "LOCK_SAT": "inject the freedom-tie contradiction into E plus both locks",
    "PROVENANCE_LOCK_A": "inject the target-blind force-A bridge",
    "LOCK_A_WITNESS_VALUE": "change mu_R in the exact A non-entailment witness",
    "PROVENANCE_LOCK_B": "inject the target-blind force-B bridge",
    "LOCK_B_WITNESS_VALUE": "change c_E and its dependent base witness values",
    "DELTA_R": "drop lock B from the private after ideal",
    "PRODUCTION_VERDICT": "feed decide the independently computed force-B case",
    "FIELD_OVERLAY_DET_ON_CONE": "use q=B_eff/rho_br; determinant residual survives but the shear guard fails",
    "FIELD_OVERLAY_POLES": "drop the C_hu coupling from the private pole polynomial",
    "OPEN_110_CARRY": "replace the on-cone residual by a C_hu-independent residual",
    "EARNED_VS_CALIBRATED_PARTITION": "restore the source bug by inserting L_A into earned equalities",
    "CONDITIONALITY_FREEDOM": "suppress the tied branch before recomputing its verdict",
    "CTRL_WELL_POSED": "remove R5 from the synthetic well-posed inventory",
    "CTRL_ABSENT": "complete R1..R5 in the absent inventory",
    "CTRL_PARTIAL_INVENTORY": "complete R3..R5 in the partial inventory",
    "CTRL_FORCED_LOCK": "drop the synthetic bridge",
    "CTRL_A_ONLY_PARTIAL": "drop the force-A bridge",
    "CTRL_B_ONLY_PARTIAL": "drop the force-B bridge",
    "CTRL_OVER_CONSTRAINED": "drop the over-stability contradiction",
    "CTRL_FREEDOM_TIE": "drop the freedom-tie contradiction",
    "VERDICT_REDERIVATION": "replace only the forced-lock table row by a force-A computed case",
    "DIMENSION_HOMOGENEITY": "change [c_E] from L T^-1 to L^2 T^-1",
    "DUAL_ENGINE_TERMS": "drop the locally computed lock-B witness lane",
    "SOURCE_TO_STAGE_MANIFEST": "drop one source predicate, mis-scope another, and drop one live registry tooth",
}


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


# ---------------------------------------------------------------------------
# Neutral source facts and exact polynomial system
# ---------------------------------------------------------------------------

rho, rho_br, rho_B0, chi_c, mu_R, K, m = sp.symbols(
    "rho rho_br rho_B0 chi_c mu_R K m"
)
M_h, c_E, C_hu, K_h, B_eff, sigma = sp.symbols(
    "M_h c_E C_hu K_h B_eff sigma"
)
tau_forced, tau_A, tau_B, eta_over, q_h_sq = sp.symbols(
    "tau_forced tau_A tau_B eta_over q_h_sq"
)
k, q = sp.symbols("k q")

BASE_VARS = (
    rho,
    rho_br,
    rho_B0,
    chi_c,
    mu_R,
    K,
    m,
    M_h,
    c_E,
    C_hu,
    K_h,
    B_eff,
    sigma,
)
BASE_POSITIVE = (
    rho,
    rho_br,
    chi_c,
    mu_R,
    K,
    m,
    M_h,
    c_E,
    K_h,
    B_eff,
    sigma,
)
BASE_EQS = (
    B_eff * chi_c - rho_B0**2,
    K_h - M_h * c_E**2,
    B_eff * K_h - C_hu**2 - sigma,
)
LOCK_A = m * mu_R - 5 * K * rho**4 * rho_br
LOCK_B = c_E**2 * rho_br - mu_R
LOCKS = {"A": LOCK_A, "B": LOCK_B}

R_KINDS = {
    "R1": "nonlinear_action_with_gnls_and_u",
    "R2": "in_plane_shear_profile_fu",
    "R3": "common_normalization_rho_mu",
    "R4": "reduction_mu_over_rho_equals_cs",
    "R5": "no_residual_geometric_factor",
}


@dataclass(frozen=True)
class SourceFact:
    kind: str
    param: str | None = None


def facts(*kinds: str) -> tuple[SourceFact, ...]:
    return tuple(SourceFact(kind) for kind in kinds)


BASELINE_FACTS = (
    *facts(
        "candidate_bridge",
        "h_goldstone_profile_imported",
        "postulated_mu_rho",
        "surface_shear_postulated",
        "missing_closed_fu",
        "route_b_distinct_dof",
        "over_import_guard",
    ),
    SourceFact("freedom_free_parameter", "C_hu"),
    SourceFact("freedom_free_parameter", "rho_br"),
)


def add_facts(
    base: Iterable[SourceFact], *new_facts: SourceFact
) -> tuple[SourceFact, ...]:
    return tuple(base) + tuple(new_facts)


def kinds_of(source_facts: Iterable[SourceFact]) -> frozenset[str]:
    return frozenset(item.kind for item in source_facts)


def route_a(source_facts: tuple[SourceFact, ...]) -> dict[str, Any]:
    kinds = kinds_of(source_facts)
    present = {
        key: R_KINDS[key] in kinds
        for key in ("R1", "R2", "R3", "R4")
    }
    present["R5"] = (
        R_KINDS["R5"] in kinds if present["R4"] else False
    )
    missing = tuple(key for key in R_KINDS if not present[key])
    grade = (
        "ROUTE_A_WELL_POSED"
        if not missing
        else (
            "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT"
            if "candidate_bridge" in kinds
            else "ROUTE_A_ABSENT"
        )
    )
    return {"grade": grade, "missing": missing}


def route_b(source_facts: tuple[SourceFact, ...]) -> str:
    kinds = kinds_of(source_facts)
    if (
        "route_b_distinct_dof" in kinds
        and "thin_plate_over_import_relation" not in kinds
    ):
        return "ROUTE_B_CLOSED_CHECKED_NEGATIVE"
    return "ROUTE_B_GUARD_FAIL"


def freedom(source_facts: tuple[SourceFact, ...]) -> dict[str, Any]:
    if "freedom_tie_relation" in kinds_of(source_facts):
        return {"status": "FREEDOM_TIED", "free_params": ()}
    free_params = tuple(
        sorted(
            item.param
            for item in source_facts
            if item.kind == "freedom_free_parameter"
            and item.param is not None
        )
    )
    status = (
        "FREEDOM_UNCONSTRAINED{" + ",".join(free_params) + "}"
        if free_params
        else "FREEDOM_TIED"
    )
    return {"status": status, "free_params": free_params}


def case_facts(case: str) -> tuple[SourceFact, ...]:
    route_support = facts("route_b_distinct_dof", "over_import_guard")
    free_support = (
        SourceFact("freedom_free_parameter", "C_hu"),
        SourceFact("freedom_free_parameter", "rho_br"),
    )
    if case == "production":
        return BASELINE_FACTS
    if case == "well_posed":
        return (
            *facts("candidate_bridge", *R_KINDS.values()),
            *route_support,
            *free_support,
        )
    if case == "absent":
        return (*facts("postulated_mu_rho"), *route_support, *free_support)
    if case == "partial_inventory":
        return (
            *facts(
                "candidate_bridge",
                R_KINDS["R1"],
                R_KINDS["R2"],
                "postulated_mu_rho",
            ),
            *route_support,
            *free_support,
        )
    relation = {
        "forced_lock": "forced_lock_fake_relation",
        "a_only_partial": "force_A_fake_relation",
        "b_only_partial": "force_B_fake_relation",
        "over_constrained": "over_stability_relation",
        "freedom_tie": "freedom_tie_relation",
        "a_inconclusive": "inconclusive_A_witness",
    }.get(case)
    if relation is None:
        raise ValueError(case)
    return add_facts(BASELINE_FACTS, SourceFact(relation))


def source_relations(
    source_facts: tuple[SourceFact, ...],
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    kinds = kinds_of(source_facts)
    eqs: list[sp.Expr] = []
    variables: list[sp.Symbol] = []
    positives: list[sp.Symbol] = []
    if "forced_lock_fake_relation" in kinds:
        eqs.extend(
            (
                mu_R - rho_br * tau_forced,
                c_E**2 - tau_forced,
                m * tau_forced - 5 * K * rho**4,
            )
        )
        variables.append(tau_forced)
        positives.append(tau_forced)
    if "force_A_fake_relation" in kinds:
        eqs.extend(
            (
                mu_R - rho_br * tau_A,
                m * tau_A - 5 * K * rho**4,
            )
        )
        variables.append(tau_A)
        positives.append(tau_A)
    if "force_B_fake_relation" in kinds:
        eqs.extend((mu_R - rho_br * tau_B, c_E**2 - tau_B))
        variables.append(tau_B)
        positives.append(tau_B)
    if "over_stability_relation" in kinds:
        eqs.append(C_hu**2 - B_eff * K_h - eta_over)
        variables.append(eta_over)
        positives.append(eta_over)
    if "freedom_tie_relation" in kinds:
        eqs.extend(
            (
                C_hu**2 - q_h_sq,
                q_h_sq * rho_br - 2 * B_eff * M_h * mu_R,
            )
        )
        variables.append(q_h_sq)
        positives.append(q_h_sq)
    return tuple(eqs), tuple(variables), tuple(positives)


def variables_for(source_facts: tuple[SourceFact, ...]) -> tuple[sp.Symbol, ...]:
    return BASE_VARS + source_relations(source_facts)[1]


def positive_for(source_facts: tuple[SourceFact, ...]) -> tuple[sp.Symbol, ...]:
    return BASE_POSITIVE + source_relations(source_facts)[2]


def equations_for(
    source_facts: tuple[SourceFact, ...],
    lock_names: tuple[str, ...] = (),
) -> tuple[sp.Expr, ...]:
    relations = source_relations(source_facts)[0]
    return tuple(
        sp.factor(sp.expand(expr))
        for expr in (
            *BASE_EQS,
            *relations,
            *(LOCKS[name] for name in lock_names),
        )
    )


def base_witness() -> dict[sp.Symbol, sp.Expr]:
    return {
        rho: sp.Integer(1),
        rho_br: sp.Integer(1),
        rho_B0: sp.Integer(2),
        chi_c: sp.Integer(4),
        mu_R: sp.Integer(1),
        K: sp.Integer(1),
        m: sp.Integer(5),
        M_h: sp.Integer(1),
        c_E: sp.Integer(1),
        C_hu: sp.Rational(1, 2),
        K_h: sp.Integer(1),
        B_eff: sp.Integer(1),
        sigma: sp.Rational(3, 4),
    }


def unlocked_witness() -> dict[sp.Symbol, sp.Expr]:
    witness = base_witness()
    witness.update(
        {
            mu_R: sp.Integer(2),
            c_E: sp.Integer(3),
            K_h: sp.Integer(9),
            sigma: sp.Rational(35, 4),
        }
    )
    return witness


def witness_for(
    source_facts: tuple[SourceFact, ...],
    lock_names: tuple[str, ...],
) -> dict[sp.Symbol, sp.Expr] | None:
    kinds = kinds_of(source_facts)
    if "over_stability_relation" in kinds:
        return None
    if "freedom_tie_relation" in kinds:
        if lock_names:
            return None
        witness = base_witness()
        witness.update(
            {
                mu_R: sp.Integer(1),
                c_E: sp.Integer(2),
                K_h: sp.Integer(4),
                C_hu: sp.sqrt(2),
                sigma: sp.Integer(2),
                q_h_sq: sp.Integer(2),
            }
        )
        return witness
    if lock_names == ("A", "B") or "forced_lock_fake_relation" in kinds:
        witness = base_witness()
    elif "force_A_fake_relation" in kinds:
        witness = base_witness()
        witness.update(
            {
                c_E: sp.Integer(2),
                K_h: sp.Integer(4),
                sigma: sp.Rational(15, 4),
            }
        )
    elif "force_B_fake_relation" in kinds:
        witness = base_witness()
        witness.update(
            {
                mu_R: sp.Integer(4),
                c_E: sp.Integer(2),
                K_h: sp.Integer(4),
                sigma: sp.Rational(15, 4),
            }
        )
    else:
        witness = unlocked_witness()
    if "forced_lock_fake_relation" in kinds:
        witness[tau_forced] = sp.Integer(1)
    if "force_A_fake_relation" in kinds:
        witness[tau_A] = sp.Integer(1)
    if "force_B_fake_relation" in kinds:
        witness[tau_B] = witness[c_E] ** 2
    return witness


def verify_witness(
    source_facts: tuple[SourceFact, ...],
    witness: dict[sp.Symbol, sp.Expr],
    lock_names: tuple[str, ...],
) -> tuple[bool, tuple[str, ...]]:
    failures = [
        sp.sstr(expr)
        for expr in equations_for(source_facts, lock_names)
        if sp.simplify(expr.subs(witness)) != 0
    ]
    for variable in positive_for(source_facts):
        if variable not in witness or not bool(sp.N(witness[variable]) > 0):
            failures.append(f"{variable} not positive")
    if sp.simplify(witness.get(rho_B0, 0)) == 0:
        failures.append("rho_B0 is zero")
    return not failures, tuple(failures)


def unsat_certificate(
    source_facts: tuple[SourceFact, ...],
    lock_names: tuple[str, ...],
) -> str | None:
    kinds = kinds_of(source_facts)
    if "over_stability_relation" in kinds:
        return (
            "UNSAT: stability plus over-stability gives "
            "-(sigma+eta_over)=0 with sigma,eta_over>0"
        )
    if "freedom_tie_relation" in kinds and "B" in lock_names:
        return (
            "UNSAT: tie plus L_B and K_h=M_h*c_E^2 gives "
            "C_hu^2=2*B_eff*K_h, contradicting positive stability slack"
        )
    return None


@lru_cache(maxsize=None)
def sat_status(
    source_facts: tuple[SourceFact, ...],
    lock_names: tuple[str, ...],
) -> dict[str, Any]:
    certificate = unsat_certificate(source_facts, lock_names)
    if certificate is not None:
        return {"status": "UNSAT", "certificate": certificate}
    witness = witness_for(source_facts, lock_names)
    if witness is None:
        return {"status": "SAT_UNCERTIFIED"}
    valid, failures = verify_witness(source_facts, witness, lock_names)
    return {
        "status": "SAT" if valid else "SAT_UNCERTIFIED",
        "witness": witness,
        "failures": failures,
    }


def groebner_dimension(
    polynomials: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> int:
    if not polynomials:
        return len(variables)
    basis = sp.groebner(polynomials, *variables, order="grevlex")
    leading_supports: list[frozenset[int]] = []
    for polynomial in basis.polys:
        monomial = tuple(polynomial.LM(order=basis.order))
        support = frozenset(i for i, exponent in enumerate(monomial) if exponent)
        if not support:
            return -1
        leading_supports.append(support)
    dimension = -1
    variable_indices = range(len(variables))
    for size in range(len(variables) + 1):
        for subset in itertools.combinations(variable_indices, size):
            selected = frozenset(subset)
            if all(
                not support.issubset(selected)
                for support in leading_supports
            ):
                dimension = size
    return dimension


@lru_cache(maxsize=None)
def dimension_record(
    source_facts: tuple[SourceFact, ...],
    after_locks: tuple[str, ...] = ("A", "B"),
) -> dict[str, Any]:
    before_sat = sat_status(source_facts, ())
    after_sat = sat_status(source_facts, after_locks)
    if before_sat["status"] != "SAT" or after_sat["status"] != "SAT":
        return {"status": "NOT_RUN"}
    variables = variables_for(source_facts)
    before = groebner_dimension(equations_for(source_facts), variables)
    after = groebner_dimension(
        equations_for(source_facts, after_locks), variables
    )
    return {
        "status": "CERTIFIED",
        "dim_before": before,
        "dim_after": after,
        "delta_r": before - after,
    }


def nonentailment_witness(
    source_facts: tuple[SourceFact, ...], lock_name: str
) -> dict[sp.Symbol, sp.Expr] | None:
    kinds = kinds_of(source_facts)
    if "forced_lock_fake_relation" in kinds:
        return None
    if lock_name == "A" and (
        "force_A_fake_relation" in kinds
        or "inconclusive_A_witness" in kinds
    ):
        return None
    if lock_name == "B" and "force_B_fake_relation" in kinds:
        return None
    if (
        "over_stability_relation" in kinds
        or "freedom_tie_relation" in kinds
    ):
        return None
    if "force_A_fake_relation" in kinds or "force_B_fake_relation" in kinds:
        return witness_for(source_facts, ())
    return unlocked_witness()


@lru_cache(maxsize=None)
def provenance(
    source_facts: tuple[SourceFact, ...], lock_name: str
) -> dict[str, Any]:
    variables = variables_for(source_facts)
    basis = sp.groebner(
        equations_for(source_facts), *variables, order="grevlex"
    )
    remainder = sp.factor(basis.reduce(LOCKS[lock_name])[1])
    if remainder == 0:
        return {"status": "ENTAILED", "remainder": sp.Integer(0)}
    witness = nonentailment_witness(source_facts, lock_name)
    if witness is None:
        return {"status": "INCONCLUSIVE", "remainder": remainder}
    valid, failures = verify_witness(source_facts, witness, ())
    value = sp.simplify(LOCKS[lock_name].subs(witness))
    if valid and value != 0:
        return {
            "status": "WITNESSED",
            "remainder": remainder,
            "witness": witness,
            "value": value,
        }
    return {
        "status": "INCONCLUSIVE",
        "remainder": remainder,
        "value": value,
        "failures": failures,
    }


def decide(
    route_a_grade: str,
    sector_status: str,
    lock_status: str,
    prov_a: str,
    prov_b: str,
    dimension_status: str,
    delta_r: int | None,
) -> tuple[str, tuple[str, ...]]:
    riders: list[str] = []
    if sector_status == "UNSAT":
        return "NO_GO(sector-ledger)", tuple(riders)
    if route_a_grade == "ROUTE_A_WELL_POSED":
        return "HALT_ROUTE_A_WELL_POSED", tuple(riders)
    if lock_status == "UNSAT":
        return "NO_GO(cone-lock)", tuple(riders)
    if lock_status == "SAT_UNCERTIFIED":
        riders.append("SAT_UNCERTIFIED")
    if prov_a == "INCONCLUSIVE":
        riders.append("ENTAILMENT_INCONCLUSIVE(L_A)")
    if prov_b == "INCONCLUSIVE":
        riders.append("ENTAILMENT_INCONCLUSIVE(L_B)")
    if riders:
        return "CONE_LOCK_PROVENANCE_INCONCLUSIVE", tuple(riders)
    if dimension_status != "CERTIFIED" or lock_status != "SAT":
        return "CONE_LOCK_PROVENANCE_INCONCLUSIVE", tuple(riders)
    if prov_a == "ENTAILED" and prov_b == "ENTAILED" and delta_r == 0:
        return "CONE_LOCK_DERIVED", tuple(riders)
    if prov_a == "ENTAILED" and prov_b == "WITNESSED" and delta_r == 1:
        return (
            "CONE_LOCK_PARTIAL(derived=A, calibrated=B)",
            tuple(riders),
        )
    if prov_b == "ENTAILED" and prov_a == "WITNESSED" and delta_r == 1:
        return (
            "CONE_LOCK_PARTIAL(derived=B, calibrated=A)",
            tuple(riders),
        )
    if prov_a == "WITNESSED" and prov_b == "WITNESSED" and delta_r == 2:
        return "CONE_LOCK_CALIBRATED", tuple(riders)
    if prov_a == "WITNESSED" and prov_b == "WITNESSED" and delta_r == 1:
        return (
            "CONE_LOCK_SHARED_CALIBRATION(delta_r=1, derived=none)",
            tuple(riders),
        )
    return "CONE_LOCK_PROVENANCE_INCONCLUSIVE", tuple(riders)


@lru_cache(maxsize=None)
def run_case(source_facts: tuple[SourceFact, ...]) -> dict[str, Any]:
    route_a_result = route_a(source_facts)
    route_b_result = route_b(source_facts)
    freedom_result = freedom(source_facts)
    if route_a_result["grade"] == "ROUTE_A_WELL_POSED":
        sector = "NOT_RUN"
        locked = "NOT_RUN"
        prov_a = "NOT_RUN"
        prov_b = "NOT_RUN"
        dimension = {"status": "NOT_RUN"}
    else:
        sector = sat_status(source_facts, ())["status"]
        locked = sat_status(source_facts, ("A", "B"))["status"]
        if sector == "SAT" and locked == "SAT":
            prov_a = provenance(source_facts, "A")["status"]
            prov_b = provenance(source_facts, "B")["status"]
            dimension = dimension_record(source_facts)
        else:
            prov_a = "NOT_RUN"
            prov_b = "NOT_RUN"
            dimension = {"status": "NOT_RUN"}
    verdict, riders = decide(
        route_a_result["grade"],
        sector,
        locked,
        prov_a,
        prov_b,
        dimension["status"],
        dimension.get("delta_r"),
    )
    return {
        "route_a": route_a_result,
        "route_b": route_b_result,
        "freedom": freedom_result,
        "sector_sat": sector,
        "lock_sat": locked,
        "prov_A": prov_a,
        "prov_B": prov_b,
        "dimension": dimension,
        "verdict": verdict,
        "riders": riders,
    }


# ---------------------------------------------------------------------------
# Field overlay, dimensions, partition, and bounded manifest
# ---------------------------------------------------------------------------


def field_overlay(
    *,
    competing_cone: bool = False,
    drop_coupling_for_poles: bool = False,
    open_condition_mutation: bool = False,
) -> dict[str, Any]:
    q_gamma = B_eff / rho_br if competing_cone else mu_R / rho_br
    determinant_q = (
        (rho_br * q - B_eff) * (M_h * q - K_h) - C_hu**2
    ) * k**4
    det_on_cone = sp.factor(
        determinant_q.subs(
            {
                q: q_gamma,
                K_h: M_h * mu_R / rho_br,
            }
        )
    )
    shear_residual = sp.factor(rho_br * q_gamma - mu_R)

    coupling = sp.Integer(0) if drop_coupling_for_poles else C_hu**2
    pole_polynomial = (
        (rho_br * q - B_eff) * (M_h * q - K_h) - coupling
    )
    coefficients = sp.Poly(pole_polynomial, q)
    a_coeff, b_coeff, c_coeff = coefficients.all_coeffs()
    discriminant = sp.factor(b_coeff**2 - 4 * a_coeff * c_coeff)
    roots = (
        sp.factor((-b_coeff - sp.sqrt(discriminant)) / (2 * a_coeff)),
        sp.factor((-b_coeff + sp.sqrt(discriminant)) / (2 * a_coeff)),
    )
    deltas_under_b = tuple(
        sp.factor(
            (root - mu_R / rho_br).subs(
                mu_R, K_h * rho_br / M_h
            )
        )
        for root in roots
    )
    expected_discriminant = (
        B_eff**2 * M_h**2
        - 2 * B_eff * K_h * M_h * rho_br
        + 4 * C_hu**2 * M_h * rho_br
        + K_h**2 * rho_br**2
    )
    expected_poles = (
        (
            B_eff * M_h
            - K_h * rho_br
            - sp.sqrt(expected_discriminant)
        )
        / (2 * M_h * rho_br),
        (
            B_eff * M_h
            - K_h * rho_br
            + sp.sqrt(expected_discriminant)
        )
        / (2 * M_h * rho_br),
    )
    coincidence_residual = (
        -k**4 if open_condition_mutation else det_on_cone
    )
    coincidence_solutions = tuple(
        sp.solve(
            sp.Eq(coincidence_residual.subs(k, 1), 0),
            C_hu,
        )
    )
    return {
        "det_on_cone": det_on_cone,
        "shear_residual": shear_residual,
        "deltas_under_b": deltas_under_b,
        "expected_poles": expected_poles,
        "coincidence_solutions": coincidence_solutions,
    }


Dim = tuple[int, int, int]


def dadd(left: Dim, right: Dim) -> Dim:
    return tuple(a + b for a, b in zip(left, right, strict=True))


def dscale(power: int, value: Dim) -> Dim:
    return tuple(power * item for item in value)


def dimension_state(*, mutate_c_e: bool = False) -> dict[str, Any]:
    dims: dict[str, Dim] = {
        "rho": (1, -3, 0),
        "rho_br": (1, -3, 0),
        "mu_R": (1, -1, -2),
        "c_E": (0, 2 if mutate_c_e else 1, -1),
        "c_gamma": (0, 1, -1),
        "m": (1, 0, 0),
        "B_eff": (1, -1, -2),
        "M_h": (1, -1, 0),
        "K_h": (1, 1, -2),
        "C_hu": (1, 0, -2),
        "sigma": (2, 0, -4),
    }
    speed_squared = dscale(2, dims["c_gamma"])
    dims["K"] = dadd(
        dadd(speed_squared, dims["m"]),
        dscale(-4, dims["rho"]),
    )
    lock_b_dims = (
        dadd(dscale(2, dims["c_E"]), dims["rho_br"]),
        dims["mu_R"],
    )
    lock_a_dims = (
        dadd(dims["m"], dims["mu_R"]),
        dadd(
            dadd(dims["K"], dscale(4, dims["rho"])),
            dims["rho_br"],
        ),
    )
    stability_dims = (
        dadd(dims["B_eff"], dims["K_h"]),
        dscale(2, dims["C_hu"]),
        dims["sigma"],
    )
    return {
        "lock_A": lock_a_dims,
        "lock_B": lock_b_dims,
        "stability": stability_dims,
        "homogeneous": (
            lock_a_dims[0] == lock_a_dims[1]
            and lock_b_dims[0] == lock_b_dims[1]
            and len(set(stability_dims)) == 1
        ),
    }


BASE_EARNED = frozenset(
    {
        "c_s^2=5*K*rho^4/m",
        "c_gamma^2=mu_R/rho_br",
        "B_eff=rho_B0^2/chi_c",
        "K_h=M_h*c_E^2",
        "B_eff*K_h-C_hu^2=sigma",
    }
)


def partition_from_provenance(
    statuses: dict[str, str],
) -> tuple[frozenset[str], frozenset[str]]:
    earned = set(BASE_EARNED)
    calibrated: set[str] = set()
    for name in ("A", "B"):
        if statuses[name] == "ENTAILED":
            earned.add(f"L_{name}")
        if statuses[name] == "WITNESSED":
            calibrated.add(f"L_{name}")
    return frozenset(earned), frozenset(calibrated)


def partition_is_biconditional(
    statuses: dict[str, str],
    earned: frozenset[str],
    calibrated: frozenset[str],
) -> bool:
    return all(
        ((f"L_{name}" in earned) == (statuses[name] == "ENTAILED"))
        and (
            (f"L_{name}" in calibrated)
            == (statuses[name] == "WITNESSED")
        )
        for name in ("A", "B")
    )


SOURCE_PREDICATE_TOTAL = 22
SOURCE_PREDICATE_UNIVERSE = (
    "pathA40.py::route_a_inventory",
    "pathA40.py::route_b_check",
    "pathA40.py::freedom_check",
    "pathA40.py::sat_decision(sector)",
    "pathA40.py::sat_decision(locks)",
    "pathA40.py::entailment_status(A)",
    "pathA40.py::entailment_status(B)",
    "pathA40.py::dimension_delta",
    "pathA40.py::decide",
    "pathA40.py+pathA40.wl::field_overlay",
    "pathA40.py::control.well_posed",
    "pathA40.py::control.absent",
    "pathA40.py::control.partial_inventory",
    "pathA40.py::control.forced_lock",
    "pathA40.py::control.a_only_partial",
    "pathA40.py::control.b_only_partial",
    "pathA40.py::control.over_constrained",
    "pathA40.py::control.freedom_tie",
    "pathA40.py::ledger.earned_equalities",
    "pathA40.py+pathA40.wl::filesystem_token_scans",
    "pathA40.py+pathA40.wl::harness_and_file_writing",
    "pathA40.py+pathA40.wl::cross_engine_payload_and_cross_read",
)
PRESERVED_SOURCE_IDS = SOURCE_PREDICATE_UNIVERSE[:18]
REPLACED_SOURCE_IDS = (SOURCE_PREDICATE_UNIVERSE[18],)
SCOPED_SOURCE_IDS = SOURCE_PREDICATE_UNIVERSE[19:]
STAGE_ONLY_IDS = (
    "stage040::CONDITIONALITY_FREEDOM",
    "stage040::DIMENSION_HOMOGENEITY",
    "stage040::DUAL_ENGINE_TERMS",
    "stage040::SOURCE_TO_STAGE_MANIFEST",
)


def manifest_rows() -> tuple[tuple[str, str], ...]:
    return (
        *((item, "preserved-folded") for item in PRESERVED_SOURCE_IDS),
        *((item, "replaced-by-stronger") for item in REPLACED_SOURCE_IDS),
        *((item, "scoped-out") for item in SCOPED_SOURCE_IDS),
        *((item, "newly-added") for item in STAGE_ONLY_IDS),
    )


EXPECTED_MANIFEST_CATEGORY = {
    "pathA40.py::route_a_inventory": "preserved-folded",
    "pathA40.py::route_b_check": "preserved-folded",
    "pathA40.py::freedom_check": "preserved-folded",
    "pathA40.py::sat_decision(sector)": "preserved-folded",
    "pathA40.py::sat_decision(locks)": "preserved-folded",
    "pathA40.py::entailment_status(A)": "preserved-folded",
    "pathA40.py::entailment_status(B)": "preserved-folded",
    "pathA40.py::dimension_delta": "preserved-folded",
    "pathA40.py::decide": "preserved-folded",
    "pathA40.py+pathA40.wl::field_overlay": "preserved-folded",
    "pathA40.py::control.well_posed": "preserved-folded",
    "pathA40.py::control.absent": "preserved-folded",
    "pathA40.py::control.partial_inventory": "preserved-folded",
    "pathA40.py::control.forced_lock": "preserved-folded",
    "pathA40.py::control.a_only_partial": "preserved-folded",
    "pathA40.py::control.b_only_partial": "preserved-folded",
    "pathA40.py::control.over_constrained": "preserved-folded",
    "pathA40.py::control.freedom_tie": "preserved-folded",
    "pathA40.py::ledger.earned_equalities": "replaced-by-stronger",
    "pathA40.py+pathA40.wl::filesystem_token_scans": "scoped-out",
    "pathA40.py+pathA40.wl::harness_and_file_writing": "scoped-out",
    "pathA40.py+pathA40.wl::cross_engine_payload_and_cross_read":
        "scoped-out",
    "stage040::CONDITIONALITY_FREEDOM": "newly-added",
    "stage040::DIMENSION_HOMOGENEITY": "newly-added",
    "stage040::DUAL_ENGINE_TERMS": "newly-added",
    "stage040::SOURCE_TO_STAGE_MANIFEST": "newly-added",
}


def manifest_state(
    rows: tuple[tuple[str, str], ...],
    live_registry: tuple[str, ...],
) -> dict[str, Any]:
    source_rows = tuple(
        row for row in rows if row[0] in SOURCE_PREDICATE_UNIVERSE
    )
    identifiers = tuple(row[0] for row in rows)
    source_identifiers = tuple(row[0] for row in source_rows)
    prefixes = {
        identifier.split("::", 1)[0] for identifier in identifiers
    }
    return {
        "valid": (
            len(source_rows) == SOURCE_PREDICATE_TOTAL
            and set(source_identifiers) == set(SOURCE_PREDICATE_UNIVERSE)
            and len(source_identifiers) == len(set(source_identifiers))
            and len(identifiers) == len(set(identifiers))
            and dict(rows) == EXPECTED_MANIFEST_CATEGORY
            and prefixes
            == {
                "pathA40.py",
                "pathA40.py+pathA40.wl",
                "stage040",
            }
            and len(live_registry) == 28
        ),
        "source_count": len(source_rows),
        "registry_size": len(live_registry),
        "partition": Counter(category for _, category in rows),
        "prefixes": prefixes,
    }


# ---------------------------------------------------------------------------
# The 28 executable teeth
# ---------------------------------------------------------------------------


def run_assertions() -> str:
    production_facts = case_facts("production")
    production = run_case(production_facts)
    production_prov_a = provenance(production_facts, "A")
    production_prov_b = provenance(production_facts, "B")
    production_dimension = dimension_record(production_facts)
    production_field = field_overlay()

    section("Computed production inventory and algebra")
    route_a_facts = (
        add_facts(
            production_facts,
            SourceFact("nonlinear_action_with_gnls_and_u"),
        )
        if ACTIVE_MUTATION == "ROUTE_A_GRADE"
        else production_facts
    )
    route_a_actual = route_a(route_a_facts)
    expect_bool(
        "ROUTE_A_GRADE",
        route_a_actual
        == {
            "grade": "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
            "missing": ("R1", "R2", "R3", "R4", "R5"),
        },
        route_a_actual,
    )

    route_b_facts = (
        add_facts(
            production_facts,
            SourceFact("thin_plate_over_import_relation"),
        )
        if ACTIVE_MUTATION == "ROUTE_B_STATUS"
        else production_facts
    )
    route_b_actual = route_b(route_b_facts)
    expect_bool(
        "ROUTE_B_STATUS",
        route_b_actual == "ROUTE_B_CLOSED_CHECKED_NEGATIVE",
        route_b_actual,
    )

    freedom_facts = (
        add_facts(production_facts, SourceFact("freedom_tie_relation"))
        if ACTIVE_MUTATION == "FREEDOM_STATUS"
        else production_facts
    )
    freedom_actual = freedom(freedom_facts)
    expect_bool(
        "FREEDOM_STATUS",
        freedom_actual
        == {
            "status": "FREEDOM_UNCONSTRAINED{C_hu,rho_br}",
            "free_params": ("C_hu", "rho_br"),
        },
        freedom_actual,
    )

    sector_facts = (
        add_facts(production_facts, SourceFact("over_stability_relation"))
        if ACTIVE_MUTATION == "SECTOR_SAT"
        else production_facts
    )
    sector_actual = sat_status(sector_facts, ())["status"]
    expect_bool("SECTOR_SAT", sector_actual == "SAT", sector_actual)

    lock_sat_facts = (
        add_facts(production_facts, SourceFact("freedom_tie_relation"))
        if ACTIVE_MUTATION == "LOCK_SAT"
        else production_facts
    )
    lock_sat_actual = sat_status(lock_sat_facts, ("A", "B"))["status"]
    expect_bool("LOCK_SAT", lock_sat_actual == "SAT", lock_sat_actual)

    prov_a_facts = (
        add_facts(production_facts, SourceFact("force_A_fake_relation"))
        if ACTIVE_MUTATION == "PROVENANCE_LOCK_A"
        else production_facts
    )
    prov_a_actual = provenance(prov_a_facts, "A")
    expect_bool(
        "PROVENANCE_LOCK_A",
        prov_a_actual["status"] == "WITNESSED"
        and prov_a_actual["remainder"] != 0,
        prov_a_actual,
    )

    witness_a = unlocked_witness()
    if ACTIVE_MUTATION == "LOCK_A_WITNESS_VALUE":
        witness_a[mu_R] = sp.Integer(1)
    witness_a_valid, witness_a_failures = verify_witness(
        production_facts, witness_a, ()
    )
    witness_a_value = sp.simplify(LOCK_A.subs(witness_a))
    expect_bool(
        "LOCK_A_WITNESS_VALUE",
        witness_a_valid and witness_a_value == 5,
        {"value": witness_a_value, "failures": witness_a_failures},
    )

    prov_b_facts = (
        add_facts(production_facts, SourceFact("force_B_fake_relation"))
        if ACTIVE_MUTATION == "PROVENANCE_LOCK_B"
        else production_facts
    )
    prov_b_actual = provenance(prov_b_facts, "B")
    expect_bool(
        "PROVENANCE_LOCK_B",
        prov_b_actual["status"] == "WITNESSED"
        and prov_b_actual["remainder"] != 0,
        prov_b_actual,
    )

    witness_b = unlocked_witness()
    if ACTIVE_MUTATION == "LOCK_B_WITNESS_VALUE":
        witness_b.update(
            {
                c_E: sp.Integer(2),
                K_h: sp.Integer(4),
                sigma: sp.Rational(15, 4),
            }
        )
    witness_b_valid, witness_b_failures = verify_witness(
        production_facts, witness_b, ()
    )
    witness_b_value = sp.simplify(LOCK_B.subs(witness_b))
    witness_r_cone = sp.simplify(
        witness_b[c_E] ** 2 * witness_b[rho_br] / witness_b[mu_R]
    )
    witness_r_cone_lock_value = sp.simplify(
        witness_b[mu_R] * (witness_r_cone - 1)
    )
    expect_bool(
        "LOCK_B_WITNESS_VALUE",
        witness_b_valid
        and witness_b_value == 7
        and witness_r_cone == sp.Rational(9, 2)
        and witness_r_cone != 1
        and witness_r_cone_lock_value == 7,
        {
            "value": witness_b_value,
            "r_cone": witness_r_cone,
            "mu_R*(r_cone-1)": witness_r_cone_lock_value,
            "failures": witness_b_failures,
        },
    )

    dimension_actual = dimension_record(
        production_facts,
        ("A",)
        if ACTIVE_MUTATION == "DELTA_R"
        else ("A", "B"),
    )
    expect_bool(
        "DELTA_R",
        dimension_actual["status"] == "CERTIFIED"
        and dimension_actual["dim_before"] == 10
        and dimension_actual["dim_after"] == 8
        and dimension_actual["delta_r"] == 2,
        dimension_actual,
    )

    verdict_case = (
        run_case(case_facts("b_only_partial"))
        if ACTIVE_MUTATION == "PRODUCTION_VERDICT"
        else production
    )
    verdict_pair = (verdict_case["verdict"], verdict_case["riders"])
    expect_bool(
        "PRODUCTION_VERDICT",
        verdict_pair == ("CONE_LOCK_CALIBRATED", ()),
        verdict_pair,
    )

    section("Computed field overlay and OPEN_110")
    det_field = field_overlay(
        competing_cone=(
            ACTIVE_MUTATION == "FIELD_OVERLAY_DET_ON_CONE"
        )
    )
    det_residual_ok = (
        sp.factor(det_field["det_on_cone"] + C_hu**2 * k**4) == 0
    )
    shear_guard_ok = sp.factor(det_field["shear_residual"]) == 0
    expect_bool(
        "FIELD_OVERLAY_DET_ON_CONE",
        det_residual_ok and shear_guard_ok,
        {
            "det_on_cone": det_field["det_on_cone"],
            "det_residual_passed": det_residual_ok,
            "direct_shear_residual": det_field["shear_residual"],
        },
    )

    pole_field = field_overlay(
        drop_coupling_for_poles=(
            ACTIVE_MUTATION == "FIELD_OVERLAY_POLES"
        )
    )
    poles_ok = all(
        sp.simplify(actual - expected) == 0
        for actual, expected in zip(
            pole_field["deltas_under_b"],
            pole_field["expected_poles"],
            strict=True,
        )
    )
    expect_bool(
        "FIELD_OVERLAY_POLES",
        poles_ok,
        pole_field["deltas_under_b"],
    )

    open_field = field_overlay(
        open_condition_mutation=(ACTIVE_MUTATION == "OPEN_110_CARRY")
    )
    expect_bool(
        "OPEN_110_CARRY",
        open_field["coincidence_solutions"] == (sp.Integer(0),),
        open_field["coincidence_solutions"],
    )

    section("Earned/calibrated partition and freedom conditionality")
    production_statuses = {
        "A": production_prov_a["status"],
        "B": production_prov_b["status"],
    }
    earned, calibrated = partition_from_provenance(production_statuses)
    if ACTIVE_MUTATION == "EARNED_VS_CALIBRATED_PARTITION":
        earned = frozenset((*earned, "L_A"))
    forced_statuses = {
        "A": provenance(case_facts("forced_lock"), "A")["status"],
        "B": provenance(case_facts("forced_lock"), "B")["status"],
    }
    forced_earned, forced_calibrated = partition_from_provenance(
        forced_statuses
    )
    partition_ok = (
        earned.isdisjoint(calibrated)
        and partition_is_biconditional(
            production_statuses, earned, calibrated
        )
        and partition_is_biconditional(
            forced_statuses, forced_earned, forced_calibrated
        )
        and calibrated == frozenset({"L_A", "L_B"})
        and earned == BASE_EARNED
        and forced_earned == frozenset((*BASE_EARNED, "L_A", "L_B"))
        and not forced_calibrated
    )
    expect_bool(
        "EARNED_VS_CALIBRATED_PARTITION",
        partition_ok,
        {
            "production": (production_statuses, earned, calibrated),
            "forced": (
                forced_statuses,
                forced_earned,
                forced_calibrated,
            ),
        },
    )

    tied_facts = (
        production_facts
        if ACTIVE_MUTATION == "CONDITIONALITY_FREEDOM"
        else case_facts("freedom_tie")
    )
    tied_case = run_case(tied_facts)
    conditionality = (
        production["freedom"]["status"]
        == "FREEDOM_UNCONSTRAINED{C_hu,rho_br}"
        and production["verdict"] == "CONE_LOCK_CALIBRATED"
        and tied_case["freedom"]["status"] == "FREEDOM_TIED"
        and tied_case["verdict"] == "NO_GO(cone-lock)"
    )
    expect_bool(
        "CONDITIONALITY_FREEDOM",
        conditionality,
        {
            "free": (
                production["freedom"]["status"],
                production["verdict"],
            ),
            "tied": (
                tied_case["freedom"]["status"],
                tied_case["verdict"],
            ),
        },
    )

    section("Individually computed falsifiability controls")
    well_posed_facts = case_facts("well_posed")
    if ACTIVE_MUTATION == "CTRL_WELL_POSED":
        well_posed_facts = tuple(
            item
            for item in well_posed_facts
            if item.kind != R_KINDS["R5"]
        )
    well_posed_case = run_case(well_posed_facts)
    expect_bool(
        "CTRL_WELL_POSED",
        well_posed_case["verdict"] == "HALT_ROUTE_A_WELL_POSED",
        well_posed_case["verdict"],
    )

    absent_facts = case_facts("absent")
    if ACTIVE_MUTATION == "CTRL_ABSENT":
        absent_facts = add_facts(
            absent_facts,
            *(SourceFact(kind) for kind in R_KINDS.values()),
        )
    absent_case = run_case(absent_facts)
    absent_tuple = (
        absent_case["route_a"]["grade"],
        absent_case["route_a"]["missing"],
        absent_case["verdict"],
        absent_case["dimension"].get("delta_r"),
    )
    expect_bool(
        "CTRL_ABSENT",
        absent_tuple
        == (
            "ROUTE_A_ABSENT",
            ("R1", "R2", "R3", "R4", "R5"),
            "CONE_LOCK_CALIBRATED",
            2,
        ),
        absent_tuple,
    )

    partial_facts = case_facts("partial_inventory")
    if ACTIVE_MUTATION == "CTRL_PARTIAL_INVENTORY":
        partial_facts = add_facts(
            partial_facts,
            SourceFact(R_KINDS["R3"]),
            SourceFact(R_KINDS["R4"]),
            SourceFact(R_KINDS["R5"]),
        )
    partial_case = run_case(partial_facts)
    partial_tuple = (
        partial_case["route_a"]["grade"],
        partial_case["route_a"]["missing"],
        partial_case["verdict"],
        partial_case["dimension"].get("delta_r"),
    )
    expect_bool(
        "CTRL_PARTIAL_INVENTORY",
        partial_tuple
        == (
            "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
            ("R3", "R4", "R5"),
            "CONE_LOCK_CALIBRATED",
            2,
        ),
        partial_tuple,
    )

    control_specs = (
        (
            "CTRL_FORCED_LOCK",
            "forced_lock",
            "CONE_LOCK_DERIVED",
            0,
        ),
        (
            "CTRL_A_ONLY_PARTIAL",
            "a_only_partial",
            "CONE_LOCK_PARTIAL(derived=A, calibrated=B)",
            1,
        ),
        (
            "CTRL_B_ONLY_PARTIAL",
            "b_only_partial",
            "CONE_LOCK_PARTIAL(derived=B, calibrated=A)",
            1,
        ),
        (
            "CTRL_OVER_CONSTRAINED",
            "over_constrained",
            "NO_GO(sector-ledger)",
            None,
        ),
        (
            "CTRL_FREEDOM_TIE",
            "freedom_tie",
            "NO_GO(cone-lock)",
            None,
        ),
    )
    for tooth, case_name, expected_verdict, expected_delta in control_specs:
        control_facts = (
            production_facts
            if ACTIVE_MUTATION == tooth
            else case_facts(case_name)
        )
        control = run_case(control_facts)
        delta = control["dimension"].get("delta_r")
        delta_ok = (
            True if expected_delta is None else delta == expected_delta
        )
        expect_bool(
            tooth,
            control["verdict"] == expected_verdict and delta_ok,
            {
                "verdict": control["verdict"],
                "delta_r": delta,
                "sector_sat": control["sector_sat"],
                "lock_sat": control["lock_sat"],
            },
        )

    section("Verdict table, dimensions, local inventory, and manifest")
    verdict_rows = (
        (
            "production",
            "production",
            ("CONE_LOCK_CALIBRATED", ()),
        ),
        (
            "well_posed",
            "well_posed",
            ("HALT_ROUTE_A_WELL_POSED", ()),
        ),
        ("absent", "absent", ("CONE_LOCK_CALIBRATED", ())),
        (
            "partial_inventory",
            "partial_inventory",
            ("CONE_LOCK_CALIBRATED", ()),
        ),
        ("forced_lock", "forced_lock", ("CONE_LOCK_DERIVED", ())),
        (
            "a_only_partial",
            "a_only_partial",
            ("CONE_LOCK_PARTIAL(derived=A, calibrated=B)", ()),
        ),
        (
            "b_only_partial",
            "b_only_partial",
            ("CONE_LOCK_PARTIAL(derived=B, calibrated=A)", ()),
        ),
        (
            "over_constrained",
            "over_constrained",
            ("NO_GO(sector-ledger)", ()),
        ),
        (
            "freedom_tie",
            "freedom_tie",
            ("NO_GO(cone-lock)", ()),
        ),
        (
            "a_inconclusive",
            "a_inconclusive",
            (
                "CONE_LOCK_PROVENANCE_INCONCLUSIVE",
                ("ENTAILMENT_INCONCLUSIVE(L_A)",),
            ),
        ),
    )
    computed_verdict_rows = []
    for row_name, case_name, expected in verdict_rows:
        actual_case_name = (
            "a_only_partial"
            if (
                ACTIVE_MUTATION == "VERDICT_REDERIVATION"
                and row_name == "forced_lock"
            )
            else case_name
        )
        computed = run_case(case_facts(actual_case_name))
        computed_verdict_rows.append(
            (
                row_name,
                (computed["verdict"], computed["riders"]),
                expected,
            )
        )
    expect_bool(
        "VERDICT_REDERIVATION",
        all(actual == expected for _, actual, expected in computed_verdict_rows),
        computed_verdict_rows,
    )

    dimension_actual = dimension_state(
        mutate_c_e=(ACTIVE_MUTATION == "DIMENSION_HOMOGENEITY")
    )
    expect_bool(
        "DIMENSION_HOMOGENEITY",
        dimension_actual["homogeneous"],
        dimension_actual,
    )

    local_inventory = {
        "verdict": production["verdict"],
        "riders": production["riders"],
        "delta_r": production_dimension["delta_r"],
        "dim_before": production_dimension["dim_before"],
        "dim_after": production_dimension["dim_after"],
        "prov_A": production_prov_a["status"],
        "prov_B": production_prov_b["status"],
        "value_A": production_prov_a["value"],
        "value_B": production_prov_b["value"],
        "freedom": production["freedom"]["status"],
        "route_A": (
            production["route_a"]["grade"],
            production["route_a"]["missing"],
        ),
        "route_B": production["route_b"],
        "det_on_cone": production_field["det_on_cone"],
        "poles": production_field["deltas_under_b"],
    }
    if ACTIVE_MUTATION == "DUAL_ENGINE_TERMS":
        del local_inventory["value_B"]
    expected_inventory = {
        "verdict": "CONE_LOCK_CALIBRATED",
        "riders": (),
        "delta_r": 2,
        "dim_before": 10,
        "dim_after": 8,
        "prov_A": "WITNESSED",
        "prov_B": "WITNESSED",
        "value_A": sp.Integer(5),
        "value_B": sp.Integer(7),
        "freedom": "FREEDOM_UNCONSTRAINED{C_hu,rho_br}",
        "route_A": (
            "ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT",
            ("R1", "R2", "R3", "R4", "R5"),
        ),
        "route_B": "ROUTE_B_CLOSED_CHECKED_NEGATIVE",
        "det_on_cone": -C_hu**2 * k**4,
        "poles": production_field["expected_poles"],
    }
    inventory_ok = (
        set(local_inventory) == set(expected_inventory)
        and all(
            (
                all(
                    sp.simplify(a - b) == 0
                    for a, b in zip(actual, expected, strict=True)
                )
                if key == "poles"
                else actual == expected
            )
            for key, expected in expected_inventory.items()
            for actual in (local_inventory.get(key),)
        )
    )
    expect_bool(
        "DUAL_ENGINE_TERMS",
        inventory_ok,
        local_inventory,
    )

    live_manifest = manifest_rows()
    live_registry = TOOTH_ORDER
    if ACTIVE_MUTATION == "SOURCE_TO_STAGE_MANIFEST":
        live_manifest = tuple(
            row
            for row in live_manifest
            if row[0] != SOURCE_PREDICATE_UNIVERSE[0]
        )
        live_manifest = tuple(
            (
                identifier,
                (
                    "scoped-out"
                    if identifier == SOURCE_PREDICATE_UNIVERSE[1]
                    else category
                ),
            )
            for identifier, category in live_manifest
        )
        live_registry = TOOTH_ORDER[:-1]
    manifest_actual = manifest_state(live_manifest, live_registry)
    expect_bool(
        "SOURCE_TO_STAGE_MANIFEST",
        manifest_actual["valid"],
        manifest_actual,
    )

    print("")
    print(f"VERDICT={production['verdict']}")
    print(
        "INTERPRETATION=CALIBRATED/UNCOMMITTED; NO committed cone lock; "
        "both locks are available calibration choices, not earned facts"
    )
    print(
        "DIMENSION="
        f"{production_dimension['dim_before']}->"
        f"{production_dimension['dim_after']}; "
        f"delta_r={production_dimension['delta_r']}"
    )
    print(
        "LOCK_PROVENANCE="
        f"L_A:{production_prov_a['status']}@{production_prov_a['value']};"
        f"L_B:{production_prov_b['status']}@{production_prov_b['value']}"
    )
    print("R_CONE=OPEN_R71; witness=9/2; c_E=c_gamma is NOT settled")
    print(
        "OPEN_110=OFF_CONE_under_AB_proportional_to_C_hu_squared_OPEN_110"
    )
    print(
        "CONDITIONAL_PROVENANCE="
        "FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu} [stage041 citation only];"
        "FREEDOM_SIM_DEFERRED{rho_br}; no stage041 cross-read"
    )
    print(
        "ROUTE_A_DEBT=REGISTERED_DEFERRED nonlinear throat solve "
        "missing {R1,R2,R3,R4,R5}"
    )
    print("SOURCE_PREDICATE_TOTAL=22")
    print("EXECUTABLE_TOOTH_TOTAL=28")
    return production["verdict"]


def main() -> None:
    if len(TOOTH_ORDER) != 28:
        raise AuditFailure("TOOTH_REGISTRY_DECLARATION")
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage040_cone_lock_readjudication SymPy audit")
    print(
        "ROUTE=exact witnesses + Groebner remainder membership + "
        "grevlex initial-monomial Krull dimension"
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
            "OVERALL FAIL: SymPy stage040 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
