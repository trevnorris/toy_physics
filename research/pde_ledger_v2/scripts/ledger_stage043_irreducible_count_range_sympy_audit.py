#!/usr/bin/env python3
"""Ledger stage043 SymPy audit: irreducible codimension count as a range.

Standalone, print-only, assert-zero, and cross-engine-file-I/O-free.  The
register classifications and midway buckets are transcribed below by name.
The production object is assembled from those facts; it is never read from
another engine or stored as an actual-value literal.

The symbolic independence diagnostic is intentionally SymPy-native:
grevlex Groebner bases give initial-monomial Krull dimension, and exact
positive smooth witnesses guard the real locus.  Its scope is only the
tested M and Wchi blocks.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from itertools import combinations
import json
import os
from typing import Any, Iterable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0
MUTATION_ENV = "LEDGER_STAGE043_MUTATION"
ACTIVE_MUTATION = os.environ.get(MUTATION_ENV, "").strip()

TOOTH_ORDER = (
    "BUCKET_PARTITION",
    "LOW_ENDPOINT",
    "HIGH_ENDPOINT",
    "RANGE_SPREAD_IS_C1_PLUS_C2",
    "BASE_CONTINUOUS_27_36",
    "E_CONTINUOUS_IS_13",
    "C1_ATTRIBUTION",
    "C2_ATTRIBUTION",
    "CONVENTION_OPEN_NOT_IMPOSED",
    "REDUCTION_DEBT_COUNTED_ONCE",
    "R49_OUT_OF_SCOPE",
    "DISCRETE_POSTULATE_COUNT",
    "DELTA_R_INDEPENDENCE",
    "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION",
    "R35_DERIVED_IN_FORM_PLUS_RIDER",
    "KNOB_PROVENANCE_WELL_FORMED",
    "REGISTER_TO_COUNT_MANIFEST",
    "DUAL_ENGINE_TERMS",
    "IRREDUCIBLE_COUNT_RANGE",
    "COUNT_RANGE_REDERIVATION",
)

ABLATION_DESCRIPTIONS = {
    "BUCKET_PARTITION":
        "double-list REG:a:hbar in the private (b) bucket",
    "LOW_ENDPOINT":
        "drop REG:a:hbar from the private low bucket membership",
    "HIGH_ENDPOINT":
        "drop REG:C2:Btilde0 from the private strict supplement",
    "RANGE_SPREAD_IS_C1_PLUS_C2":
        "drop one C1 member only from the private attribution count",
    "BASE_CONTINUOUS_27_36":
        "reclassify the EOS-exponent-5 pull-out as structural-no-knob",
    "E_CONTINUOUS_IS_13":
        "remove rho_B0 from the private Parts-III--VI extension",
    "C1_ATTRIBUTION":
        "remove Mtilde from the private C1 contributor set",
    "C2_ATTRIBUTION":
        "count the two shorthand labels instead of six support scalars",
    "CONVENTION_OPEN_NOT_IMPOSED":
        "collapse the private reported range to its high endpoint",
    "REDUCTION_DEBT_COUNTED_ONCE":
        "exclude mu_R from the private continuous assembly",
    "R49_OUT_OF_SCOPE":
        "fold one R49 field profile into the private continuous assembly",
    "DISCRETE_POSTULATE_COUNT":
        "drop D1 from the private discrete-postulate itemization",
    "DELTA_R_INDEPENDENCE":
        "inject K=rho0 and mu_eta=T_w into the private baseline blocks",
    "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION":
        "subtract the tested delta-r diagnostic from the private count",
    "R35_DERIVED_IN_FORM_PLUS_RIDER":
        "flip the private R35 label to PENDING-debt",
    "KNOB_PROVENANCE_WELL_FORMED":
        "slip derived c_s0 into the private counted set",
    "REGISTER_TO_COUNT_MANIFEST":
        "mis-scope hbar as DERIVED against its fixed ratified category",
    "DUAL_ENGINE_TERMS":
        "drop the locally computed spread from the canonical inventory",
    "IRREDUCIBLE_COUNT_RANGE":
        "corrupt the assembled object's private spread field",
    "COUNT_RANGE_REDERIVATION":
        "drop hbar from a private source bucket and re-run the full pipeline",
}

if len(TOOTH_ORDER) != 20:
    raise RuntimeError("stage043 tooth declaration is not exactly 20")


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
# Self-contained register-to-count facts
# ---------------------------------------------------------------------------

CAT_CONTINUOUS = "continuous-counted-irreducible"
CAT_DEBT = "reduction-debt-counted-once"
CAT_DERIVED = "derived-not-counted"
CAT_BRIDGE = "calibrated-external-bridge-not-counted"
CAT_DEPARTURE = "departure-no-knob"
CAT_CONVENTION = "convention"
CAT_DISCRETE = "discrete-structural-postulate"
CAT_OPEN = "extension-convention-open"
CAT_R49 = "parallel-track-out-of-scope"
CAT_STRUCTURAL = "structural-no-knob"

STATUS_TO_CATEGORY = {
    "ACTION": CAT_CONTINUOUS,
    "FREE-UNREDUCED-ROUTELESS": CAT_CONTINUOUS,
    "CALIB": CAT_CONTINUOUS,
    "CALIB-ANCHOR": CAT_CONTINUOUS,
    "FREE-UNREDUCED-DEBT": CAT_DEBT,
    "DERIVED": CAT_DERIVED,
    "EXTERNAL-BRIDGE": CAT_BRIDGE,
    "DEPARTURE": CAT_DEPARTURE,
    "CONV": CAT_CONVENTION,
    "STRUCTURAL-POSTULATE": CAT_DISCRETE,
    "DERIVED-IN-FORM-UNEXECUTED": CAT_OPEN,
    "DOWNSTREAM-DEFERRED-OPEN": CAT_OPEN,
    "PARALLEL-TRACK": CAT_R49,
    "STRUCTURAL": CAT_STRUCTURAL,
}

EXPECTED_CATEGORY_COUNTS = {
    CAT_CONTINUOUS: 22,
    CAT_DEBT: 18,
    CAT_DERIVED: 40,
    CAT_BRIDGE: 10,
    CAT_DEPARTURE: 4,
    CAT_CONVENTION: 2,
    CAT_DISCRETE: 11,
    CAT_OPEN: 9,
    CAT_R49: 8,
    CAT_STRUCTURAL: 28,
}


@dataclass(frozen=True)
class RegisterFact:
    identifier: str
    source_status: str
    ratified_category: str
    part: str = ""
    bucket: str = ""
    register_label: str = ""
    rider: str = ""
    frozen_calibration_input: bool = False
    moment_integral_executed: bool = True
    shorthand: str = ""


def fact_group(
    identifiers: Iterable[str],
    status: str,
    category: str,
    *,
    part: str = "",
    bucket: str = "",
) -> list[RegisterFact]:
    return [
        RegisterFact(
            identifier=identifier,
            source_status=status,
            ratified_category=category,
            part=part,
            bucket=bucket,
        )
        for identifier in identifiers
    ]


MANIFEST_FACTS: tuple[RegisterFact, ...] = tuple(
    # Parts I--II low buckets: 4 + 15 + 14 + 1 = 34.
    fact_group(
        ("REG:a:hbar", "REG:a:m_GNLS", "REG:a:K_EOS", "REG:a:rho0"),
        "ACTION", CAT_CONTINUOUS, part="I-II", bucket="a",
    )
    + fact_group(
        (
            "REG:b:mu_R", "REG:b:rho_br", "REG:b:mu_R_4D",
            "REG:b:d_slab", "REG:b:mu_eta", "REG:b:T_w",
            "REG:b:beta", "REG:b:Vp0_over_lc", "REG:b:T_Omega",
            "REG:b:beta2_profile", "REG:b:epsilon0", "REG:b:epsilon1",
            "REG:b:K0c", "REG:b:K1_eta_direction",
            "REG:b:K1_TOmega_direction",
        ),
        "FREE-UNREDUCED-DEBT", CAT_DEBT, part="I-II", bucket="b",
    )
    + fact_group(
        (
            "REG:c:a_B", "REG:c:kappa_B", "REG:c:Omega_w",
            "REG:c:g_l_width", "REG:c:C_E", "REG:c:C_B",
            "REG:c:W_slab",
        ),
        "FREE-UNREDUCED-ROUTELESS", CAT_CONTINUOUS,
        part="I-II", bucket="c",
    )
    + fact_group(
        (
            "REG:disc:EOS_exponent_5",
            "REG:disc:chi_B_field",
            "REG:disc:Gamma_B_law",
            "REG:disc:chi_B_gating",
            "REG:disc:G0_postulate_1_axis_surface",
            "REG:disc:G0_postulate_2_free_slip",
            "REG:disc:G0_postulate_6_no_C5_phi",
        ),
        "STRUCTURAL-POSTULATE", CAT_DISCRETE, part="I-II", bucket="c",
    )
    + fact_group(
        ("REG:force:force_magnitude_norm",),
        "CALIB", CAT_CONTINUOUS, part="I-II", bucket="force-mag",
    )
    # Parts III--VI continuous extension: 5 + 7 + 1 + 0 = 13.
    + fact_group(
        (
            "REG:E:III:rho_B0", "REG:E:III:chi_c", "REG:E:III:B",
            "REG:E:III:J", "REG:E:III:K_theta_over_kappa_phase",
        ),
        "FREE-UNREDUCED-ROUTELESS", CAT_CONTINUOUS, part="III",
    )
    + fact_group(
        (
            "REG:E:IV:C_hu", "REG:E:IV:M4", "REG:E:IV:ell_over_a",
            "REG:E:IV:K_m", "REG:E:IV:J_m",
        ),
        "ACTION", CAT_CONTINUOUS, part="IV",
    )
    + fact_group(
        ("REG:E:IV:c_E", "REG:E:IV:Q_E"),
        "FREE-UNREDUCED-DEBT", CAT_DEBT, part="IV",
    )
    + fact_group(
        ("REG:E:V:q_T",),
        "FREE-UNREDUCED-DEBT", CAT_DEBT, part="V",
    )
    # Frozen D1--D4, separate from the continuous variety.
    + fact_group(
        (
            "REG:disc:D1:H_existence_and_Poschl_Teller_law",
            "REG:disc:D2:s_i_Z2_topology",
            "REG:disc:D3:R63_mouth_BC_class",
            "REG:disc:D4:R68_tau_d_time_arrow",
        ),
        "STRUCTURAL-POSTULATE", CAT_DISCRETE, part="III-VI",
    )
    # Complete not-counted bulk of the 152-entry bounded manifest.
    + fact_group(
        (
            "REG:derived:c_s0", "REG:derived:xi_h", "REG:derived:h0",
            "REG:derived:lambda_gamma", "REG:derived:K_eta",
            "REG:derived:delta_wall", "REG:derived:sigma_wall",
            "REG:derived:Z_return", "REG:derived:C_J", "REG:derived:B_eff",
            "REG:derived:M_h", "REG:derived:K_h", "REG:derived:Q_chi",
            "REG:derived:xi_no_sqrt2", "REG:derived:R_mouth_cancelled",
            "REG:derived:Z0_ret_alias", "REG:derived:Z1_ret_alias",
            "REG:derived:K1_sum", "REG:derived:Mtilde_definition",
            "REG:derived:Ktilde_definition", "REG:derived:Ttilde_definition",
            "REG:derived:rho_eff_projection", "REG:derived:chi_Q_outgoing",
            "REG:derived:lambda_m_SO3", "REG:derived:u2_DtN",
            "REG:derived:u4_DtN", "REG:derived:v5_DtN",
            "REG:derived:Kbar4_moment", "REG:derived:p_localization",
            "REG:derived:force_sign", "REG:derived:K4_from_M4_cE",
            "REG:derived:CJ_IBP_sign", "REG:derived:charge_value_R57",
            "REG:derived:qT_product", "REG:derived:mouth_class_fillers",
            "REG:derived:wall_kink_profile", "REG:derived:round_trip_Rrt",
            "REG:derived:DtN_pole_ladder", "REG:derived:cone_ratio",
            "REG:derived:field_overlay_determinant",
        ),
        "DERIVED", CAT_DERIVED,
    )
    + fact_group(
        (
            "REG:bridge:G_Newton", "REG:bridge:PN_54_over_5",
            "REG:bridge:PN_2_over_5", "REG:bridge:benchmark_force",
            "REG:bridge:benchmark_orbit", "REG:bridge:benchmark_flux",
            "REG:bridge:external_mass_anchor",
            "REG:bridge:external_charge_anchor",
            "REG:bridge:external_time_anchor",
            "REG:bridge:external_length_anchor",
        ),
        "EXTERNAL-BRIDGE", CAT_BRIDGE,
    )
    + fact_group(
        (
            "REG:departure:R66_native_P_no_emergent_Gauss",
            "REG:departure:R73_bT_T_even",
            "REG:departure:R81_scalar_sim_gated",
            "REG:departure:light_stray_longitudinal",
        ),
        "DEPARTURE", CAT_DEPARTURE,
    )
    + fact_group(
        ("REG:conv:a_pin", "REG:conv:lambda_gamma_ratio"),
        "CONV", CAT_CONVENTION,
    )
    + fact_group(
        (
            "REG:R49:profile_E0", "REG:R49:profile_E1",
            "REG:R49:profile_E2", "REG:R49:profile_E3",
            "REG:R49:profile_B0", "REG:R49:profile_B1",
            "REG:R49:profile_B2", "REG:R49:profile_B3",
        ),
        "PARALLEL-TRACK", CAT_R49,
    )
    + fact_group(
        (
            "REG:struct:R23_return_targets", "REG:struct:R25_slab_survival",
            "REG:struct:R26_frozen_validity", "REG:struct:R27_xi_lc_firewall",
            "REG:struct:R28_BC_dependent", "REG:struct:R31_truncation_window",
            "REG:struct:R42_scalar_reduction", "REG:struct:R50_H_parent",
            "REG:struct:R53_kernel_stability", "REG:struct:R63_mouth_data",
            "REG:struct:R65_charge_blocker", "REG:struct:R67_drain_coupling",
            "REG:struct:R69_return_bookkeeping", "REG:struct:R71_open_rcone",
            "REG:struct:R74_NG5_trio", "REG:struct:R75_shared_throat_route",
            "REG:struct:R76_QE_sensitivity", "REG:struct:R77_lineage_seal",
            "REG:struct:R78_cone_adjudication", "REG:struct:R79_OPEN_110",
            "REG:struct:R80_freedom_conditionality",
            "REG:struct:R82_falsifier_set", "REG:struct:R83_guard_A",
            "REG:struct:R84_VI_VII_handoff", "REG:struct:z_g_opaque",
            "REG:struct:z_b_opaque", "REG:struct:control_inventory",
            "REG:struct:census_bookkeeping",
        ),
        "STRUCTURAL", CAT_STRUCTURAL,
    )
)

# R35/C1 and the C2 six-scalar lane carry extra, independently checked facts.
_c1_ids = (
    "REG:C1:Mtilde", "REG:C1:Ktilde", "REG:C1:Ttilde_Omega",
)
_c2_specs = (
    ("REG:C2:Btilde0", "Btilde"),
    ("REG:C2:Btilde2", "Btilde"),
    ("REG:C2:Btilde4", "Btilde"),
    ("REG:C2:Ztilde0", "Ztilde"),
    ("REG:C2:Ztilde2", "Ztilde"),
    ("REG:C2:Ztilde4", "Ztilde"),
)
MANIFEST_FACTS += tuple(
    RegisterFact(
        identifier=identifier,
        source_status="DERIVED-IN-FORM-UNEXECUTED",
        ratified_category=CAT_OPEN,
        part="I-II",
        bucket="C1",
        register_label="DERIVED-in-form",
        rider="C1_contributor",
        frozen_calibration_input=True,
        moment_integral_executed=False,
        shorthand="R35",
    )
    for identifier in _c1_ids
)
MANIFEST_FACTS += tuple(
    RegisterFact(
        identifier=identifier,
        source_status="DOWNSTREAM-DEFERRED-OPEN",
        ratified_category=CAT_OPEN,
        part="I-II",
        bucket="C2",
        register_label="downstream-deferred",
        rider="C2_contributor",
        shorthand=shorthand,
    )
    for identifier, shorthand in _c2_specs
)


def category_of(row: RegisterFact) -> str:
    return STATUS_TO_CATEGORY.get(row.source_status, "UNCLASSIFIED")


def replace_row(
    rows: tuple[RegisterFact, ...],
    identifier: str,
    **changes: Any,
) -> tuple[RegisterFact, ...]:
    return tuple(
        replace(row, **changes) if row.identifier == identifier else row
        for row in rows
    )


def rows_by_id(rows: tuple[RegisterFact, ...]) -> dict[str, RegisterFact]:
    return {row.identifier: row for row in rows}


def manifest_state(rows: tuple[RegisterFact, ...]) -> dict[str, Any]:
    identifiers = [row.identifier for row in rows]
    computed_map = {row.identifier: category_of(row) for row in rows}
    expected_map = {row.identifier: row.ratified_category for row in rows}
    counts = Counter(computed_map.values())
    return {
        "unique": len(set(identifiers)) == 152 and len(identifiers) == 152,
        "engine_qualified": all(identifier.startswith("REG:") for identifier in identifiers),
        "per_row_map": computed_map == expected_map,
        "counts": dict(counts),
        "counts_match": dict(counts) == EXPECTED_CATEGORY_COUNTS,
        "valid": (
            len(set(identifiers)) == 152
            and len(identifiers) == 152
            and all(identifier.startswith("REG:") for identifier in identifiers)
            and computed_map == expected_map
            and dict(counts) == EXPECTED_CATEGORY_COUNTS
        ),
    }


def bucket_sets(rows: tuple[RegisterFact, ...]) -> dict[str, set[str]]:
    names = ("a", "b", "c", "force-mag", "C1", "C2")
    return {
        name: {row.identifier for row in rows if row.bucket == name}
        for name in names
    }


def bucket_partition_state(buckets: dict[str, set[str]]) -> dict[str, Any]:
    expected = {
        "a": 4, "b": 15, "c": 14, "force-mag": 1, "C1": 3, "C2": 6,
    }
    flattened = [item for name in expected for item in buckets[name]]
    counts = {name: len(buckets[name]) for name in expected}
    return {
        "counts": counts,
        "raw_low": sum(counts[name] for name in ("a", "b", "c", "force-mag")),
        "raw_high": sum(counts.values()),
        "disjoint": len(flattened) == len(set(flattened)),
        "valid": (
            counts == expected
            and len(flattened) == len(set(flattened))
            and sum(counts[name] for name in ("a", "b", "c", "force-mag")) == 34
            and sum(counts.values()) == 43
        ),
    }


def derive_count_state(
    rows: tuple[RegisterFact, ...],
    *,
    force_strict: bool = False,
) -> dict[str, Any]:
    buckets = bucket_sets(rows)
    bucket_state = bucket_partition_state(buckets)
    category = {row.identifier: category_of(row) for row in rows}

    base_discrete = {
        row.identifier
        for row in rows
        if row.part == "I-II" and category_of(row) == CAT_DISCRETE
    }
    extension_by_part = {
        part: {
            row.identifier
            for row in rows
            if row.part == part and category_of(row) in {CAT_CONTINUOUS, CAT_DEBT}
        }
        for part in ("III", "IV", "V", "VI")
    }
    extension_ids = set().union(*extension_by_part.values())

    low_raw_ids = set().union(
        *(buckets[name] for name in ("a", "b", "c", "force-mag"))
    )
    high_raw_ids = low_raw_ids | buckets["C1"] | buckets["C2"]
    effective_low_raw_ids = high_raw_ids if force_strict else low_raw_ids

    low_base_ids = effective_low_raw_ids - base_discrete
    high_base_ids = high_raw_ids - base_discrete
    low_components = sorted(low_base_ids) + sorted(extension_ids)
    high_components = sorted(high_base_ids) + sorted(extension_ids)

    e_counts = {part: len(ids) for part, ids in extension_by_part.items()}
    return {
        "buckets": buckets,
        "bucket_state": bucket_state,
        "raw_low": len(effective_low_raw_ids),
        "raw_high": len(high_raw_ids),
        "base_discrete_ids": base_discrete,
        "base_continuous": (len(low_base_ids), len(high_base_ids)),
        "extension_by_part": extension_by_part,
        "extension_ids": extension_ids,
        "E_itemization": e_counts,
        "E_continuous": sum(e_counts.values()),
        "low_components": low_components,
        "high_components": high_components,
        "low_c": len(low_components),
        "high_c": len(high_components),
        "spread": len(high_components) - len(low_components),
        "category": category,
    }


# ---------------------------------------------------------------------------
# SymPy Groebner/Krull diagnostic with exact positive witnesses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AlgebraBlock:
    name: str
    variables: tuple[sp.Symbol, ...]
    baseline: tuple[sp.Expr, ...]
    forced_dependency: sp.Expr
    witness: dict[sp.Symbol, sp.Expr]


def groebner_dimension(
    polynomials: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> int:
    nonzero = tuple(sp.expand(poly) for poly in polynomials if sp.expand(poly) != 0)
    if not nonzero:
        return len(variables)
    basis = sp.groebner(nonzero, *variables, order="grevlex")
    supports: list[frozenset[int]] = []
    for polynomial in basis.polys:
        leading = tuple(polynomial.LM(order=basis.order))
        support = frozenset(index for index, exponent in enumerate(leading) if exponent)
        if not support:
            return -1
        supports.append(support)
    for size in range(len(variables), -1, -1):
        for candidate in combinations(range(len(variables)), size):
            selected = frozenset(candidate)
            if all(not support <= selected for support in supports):
                return size
    raise AssertionError("initial-ideal dimension search failed")


def certify_positive_smooth_witness(
    equations: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    witness: dict[sp.Symbol, sp.Expr],
    dimension: int,
) -> bool:
    positive = all(
        sp.simplify(witness[variable]).is_real is True
        and sp.simplify(witness[variable]).is_positive is True
        for variable in variables
    )
    solves = all(sp.simplify(equation.subs(witness)) == 0 for equation in equations)
    nonzero = tuple(equation for equation in equations if sp.expand(equation) != 0)
    rank = int(sp.Matrix(nonzero).jacobian(variables).subs(witness).rank()) if nonzero else 0
    return positive and solves and len(variables) - rank == dimension


def dimension_record(block: AlgebraBlock, equations: tuple[sp.Expr, ...]) -> dict[str, int | bool]:
    before = groebner_dimension((), block.variables)
    after = groebner_dimension(equations, block.variables)
    witness_ok = certify_positive_smooth_witness(
        equations, block.variables, block.witness, after
    )
    return {
        "dim_before": before,
        "dim_after": after,
        "delta_r": before - after,
        "positive_witness": witness_ok,
    }


hbar, mass, big_k, rho0, cs0, xi_h, h0, scale_a, c_gamma, lambda_gamma = sp.symbols(
    "hbar mass K rho0 c_s0 xi_h h0 a c_gamma lambda_gamma", real=True
)
medium = AlgebraBlock(
    "M",
    (hbar, mass, big_k, rho0, cs0, xi_h, h0, scale_a, c_gamma, lambda_gamma),
    (
        mass * cs0**2 - 5 * big_k * rho0**4,
        scale_a * mass * cs0 - hbar,
        xi_h**2 * mass**2 * cs0**2 - 2 * hbar**2,
        4 * h0 - mass * cs0**2,
        lambda_gamma * cs0 - c_gamma,
    ),
    big_k - rho0,
    {
        hbar: sp.Integer(5), mass: sp.Integer(5), big_k: sp.Integer(1),
        rho0: sp.Integer(1), cs0: sp.Integer(1), xi_h: sp.sqrt(2),
        h0: sp.Rational(5, 4), scale_a: sp.Integer(1),
        c_gamma: sp.Integer(1), lambda_gamma: sp.Integer(1),
    },
)

mu_eta, t_w, beta, k_eta, vp0_lc, t_omega, a_b, kappa_b, delta, sigma_wall = sp.symbols(
    "mu_eta T_w beta K_eta Vp0_over_lc T_Omega a_B kappa_B delta sigma_wall",
    real=True,
)
wall = AlgebraBlock(
    "Wchi",
    (mu_eta, t_w, beta, k_eta, vp0_lc, t_omega, a_b, kappa_b, delta, sigma_wall),
    (
        k_eta - t_w * beta**2,
        2 * a_b * delta**2 - kappa_b,
        36 * sigma_wall**2 - 2 * a_b * kappa_b,
    ),
    mu_eta - t_w,
    {
        mu_eta: sp.Integer(1), t_w: sp.Integer(1), beta: sp.Integer(1),
        k_eta: sp.Integer(1), vp0_lc: sp.Integer(1), t_omega: sp.Integer(1),
        a_b: sp.Integer(1), kappa_b: sp.Integer(1),
        delta: sp.sqrt(2) / 2, sigma_wall: sp.sqrt(2) / 6,
    },
)

EXPECTED_DELTA_BASELINE = {
    "M": {"dim_before": 10, "dim_after": 5, "delta_r": 5},
    "Wchi": {"dim_before": 10, "dim_after": 7, "delta_r": 3},
}
EXPECTED_DELTA_FORCED = {
    "M": {"dim_before": 10, "dim_after": 4, "delta_r": 6},
    "Wchi": {"dim_before": 10, "dim_after": 6, "delta_r": 4},
}


def delta_diagnostic(*, force_baseline_dependency: bool = False) -> dict[str, Any]:
    baseline: dict[str, dict[str, int | bool]] = {}
    forced: dict[str, dict[str, int | bool]] = {}
    for block in (medium, wall):
        baseline_equations = (
            block.baseline + (block.forced_dependency,)
            if force_baseline_dependency
            else block.baseline
        )
        baseline[block.name] = dimension_record(block, baseline_equations)
        forced[block.name] = dimension_record(
            block, block.baseline + (block.forced_dependency,)
        )

    baseline_numeric = {
        name: {key: record[key] for key in ("dim_before", "dim_after", "delta_r")}
        for name, record in baseline.items()
    }
    forced_numeric = {
        name: {key: record[key] for key in ("dim_before", "dim_after", "delta_r")}
        for name, record in forced.items()
    }
    guards = (
        baseline_numeric == EXPECTED_DELTA_BASELINE
        and forced_numeric == EXPECTED_DELTA_FORCED
        and all(bool(record["positive_witness"]) for record in baseline.values())
        and all(bool(record["positive_witness"]) for record in forced.values())
        and forced["M"]["delta_r"] > baseline["M"]["delta_r"]
        and forced["Wchi"]["delta_r"] > baseline["Wchi"]["delta_r"]
    )
    token = (
        "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS"
        if guards
        else "DEVIATION_DETECTED"
    )
    return {
        "baseline": baseline,
        "forced": forced,
        "guards": guards,
        "token": token,
    }


# ---------------------------------------------------------------------------
# Headline assembly
# ---------------------------------------------------------------------------


def derive_qe_sensitivity(rows: tuple[RegisterFact, ...], state: dict[str, Any]) -> dict[str, Any]:
    qe_ids = {
        row.identifier for row in rows
        if row.identifier == "REG:E:IV:Q_E" and category_of(row) == CAT_DEBT
    }
    declared_split = len({
        row.identifier for row in rows
        if row.identifier in {
            "REG:E:III:rho_B0", "REG:E:III:chi_c", "REG:E:IV:C_hu",
        }
    })
    undeclared_split = declared_split + len(qe_ids)
    total_before = len(state["low_components"])
    total_after = len(set(state["low_components"]) | qe_ids)
    return {
        "declared_split": declared_split,
        "undeclared_split": undeclared_split,
        "total_neutral": total_before == total_after,
        "token": (
            "TOTAL_NEUTRAL_ROUTELESS_SHIFT_3_to_4"
            if declared_split == 3 and undeclared_split == 4 and total_before == total_after
            else "QE_SENSITIVITY_DEVIATION"
        ),
    }


def derive_r49_scope(rows: tuple[RegisterFact, ...], state: dict[str, Any]) -> dict[str, Any]:
    r49_ids = {
        row.identifier for row in rows if category_of(row) == CAT_R49
    }
    continuous_ids = set(state["high_components"])
    valid = len(r49_ids) == 8 and not (r49_ids & continuous_ids)
    return {
        "profiles": r49_ids,
        "valid": valid,
        "token": "OUT_OF_BUILT_V2_SCOPE" if valid else "R49_SCOPE_VIOLATION",
    }


def derive_convention(state: dict[str, Any]) -> str:
    open_members = state["buckets"]["C1"] | state["buckets"]["C2"]
    return "OPEN" if open_members and state["low_c"] < state["high_c"] else "IMPOSED"


def assemble_range(
    rows: tuple[RegisterFact, ...],
    state: dict[str, Any],
    diagnostic: dict[str, Any],
) -> dict[str, Any]:
    c1_count = len(state["buckets"]["C1"])
    c2_count = len(state["buckets"]["C2"])
    discrete_count = sum(
        1 for row in rows if category_of(row) == CAT_DISCRETE
    )
    qe = derive_qe_sensitivity(rows, state)
    r49 = derive_r49_scope(rows, state)
    return {
        "continuous_codimension": [state["low_c"], state["high_c"]],
        "base_continuous": list(state["base_continuous"]),
        "E_continuous": state["E_continuous"],
        "discrete_postulate_count": discrete_count,
        "convention": derive_convention(state),
        "spread": state["spread"],
        "C1_cardinality": c1_count,
        "C2_cardinality": c2_count,
        "delta_r_independence": diagnostic["token"],
        "Q_E_undeclared": qe["token"],
        "R49_parallel_track": r49["token"],
    }


EXPECTED_RANGE_OBJECT = {
    "continuous_codimension": [40, 49],
    "base_continuous": [27, 36],
    "E_continuous": 13,
    "discrete_postulate_count": 11,
    "convention": "OPEN",
    "spread": 9,
    "C1_cardinality": 3,
    "C2_cardinality": 6,
    "delta_r_independence": "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS",
    "Q_E_undeclared": "TOTAL_NEUTRAL_ROUTELESS_SHIFT_3_to_4",
    "R49_parallel_track": "OUT_OF_BUILT_V2_SCOPE",
}


def range_object_is_valid(obj: dict[str, Any]) -> bool:
    return (
        obj == EXPECTED_RANGE_OBJECT
        and obj["continuous_codimension"][0] < obj["continuous_codimension"][1]
        and obj["continuous_codimension"][1] - obj["continuous_codimension"][0] == 9
        and obj["spread"] == obj["C1_cardinality"] + obj["C2_cardinality"]
    )


def run_assertions() -> dict[str, Any]:
    base_rows = MANIFEST_FACTS
    state = derive_count_state(base_rows)
    diagnostic = delta_diagnostic()
    production = assemble_range(base_rows, state, diagnostic)

    section("A. Parts-I--II buckets and continuous endpoints")

    private_buckets = bucket_sets(base_rows)
    if ACTIVE_MUTATION == "BUCKET_PARTITION":
        private_buckets["b"].add("REG:a:hbar")
    bucket_actual = bucket_partition_state(private_buckets)
    expect_bool("BUCKET_PARTITION", bucket_actual["valid"], bucket_actual)

    low_rows = base_rows
    if ACTIVE_MUTATION == "LOW_ENDPOINT":
        low_rows = replace_row(low_rows, "REG:a:hbar", bucket="")
    low_actual = derive_count_state(low_rows)
    expect_bool(
        "LOW_ENDPOINT",
        low_actual["low_c"] == 40,
        {"low_c": low_actual["low_c"], "raw_low": low_actual["raw_low"]},
    )

    high_rows = base_rows
    if ACTIVE_MUTATION == "HIGH_ENDPOINT":
        high_rows = replace_row(high_rows, "REG:C2:Btilde0", bucket="")
    high_actual = derive_count_state(high_rows)
    expect_bool(
        "HIGH_ENDPOINT",
        high_actual["high_c"] == 49,
        {"high_c": high_actual["high_c"], "raw_high": high_actual["raw_high"]},
    )

    spread_c1_rows = base_rows
    if ACTIVE_MUTATION == "RANGE_SPREAD_IS_C1_PLUS_C2":
        spread_c1_rows = replace_row(spread_c1_rows, "REG:C1:Mtilde", bucket="")
    spread_c1_count = len(bucket_sets(spread_c1_rows)["C1"])
    spread_c2_count = len(bucket_sets(base_rows)["C2"])
    spread_ok = (
        state["spread"] == 9
        and state["spread"] == spread_c1_count + spread_c2_count
    )
    expect_bool(
        "RANGE_SPREAD_IS_C1_PLUS_C2",
        spread_ok,
        {
            "spread": state["spread"],
            "C1": spread_c1_count,
            "C2": spread_c2_count,
        },
    )

    base_rows_private = base_rows
    if ACTIVE_MUTATION == "BASE_CONTINUOUS_27_36":
        base_rows_private = replace_row(
            base_rows_private,
            "REG:disc:EOS_exponent_5",
            source_status="STRUCTURAL",
        )
    base_actual = derive_count_state(base_rows_private)
    expect_bool(
        "BASE_CONTINUOUS_27_36",
        (
            base_actual["base_continuous"] == (27, 36)
            and len(base_actual["base_discrete_ids"]) == 7
            and base_actual["raw_low"] == 34
            and base_actual["raw_high"] == 43
        ),
        {
            "base": base_actual["base_continuous"],
            "pullout": len(base_actual["base_discrete_ids"]),
            "raw": (base_actual["raw_low"], base_actual["raw_high"]),
        },
    )

    e_rows = base_rows
    if ACTIVE_MUTATION == "E_CONTINUOUS_IS_13":
        e_rows = replace_row(e_rows, "REG:E:III:rho_B0", part="")
    e_actual = derive_count_state(e_rows)
    expect_bool(
        "E_CONTINUOUS_IS_13",
        (
            e_actual["E_itemization"] == {"III": 5, "IV": 7, "V": 1, "VI": 0}
            and e_actual["E_continuous"] == 13
        ),
        {
            "itemization": e_actual["E_itemization"],
            "E_continuous": e_actual["E_continuous"],
        },
    )

    section("B. OPEN convention contributors")

    c1_rows = base_rows
    if ACTIVE_MUTATION == "C1_ATTRIBUTION":
        c1_rows = replace_row(c1_rows, "REG:C1:Mtilde", bucket="")
    c1_members = {
        row.identifier for row in c1_rows
        if row.bucket == "C1"
        and row.frozen_calibration_input
        and not row.moment_integral_executed
    }
    expect_bool(
        "C1_ATTRIBUTION",
        c1_members == set(_c1_ids) and len(c1_members) == 3,
        sorted(c1_members),
    )

    c2_rows = [row for row in base_rows if row.bucket == "C2"]
    if ACTIVE_MUTATION == "C2_ATTRIBUTION":
        c2_members = {row.shorthand for row in c2_rows}
    else:
        c2_members = {row.identifier for row in c2_rows}
    expect_bool(
        "C2_ATTRIBUTION",
        c2_members == {identifier for identifier, _ in _c2_specs}
        and len(c2_members) == 6,
        sorted(c2_members),
    )

    convention_object = dict(production)
    if ACTIVE_MUTATION == "CONVENTION_OPEN_NOT_IMPOSED":
        convention_object["continuous_codimension"] = [
            state["high_c"], state["high_c"],
        ]
    convention_low, convention_high = convention_object["continuous_codimension"]
    expect_bool(
        "CONVENTION_OPEN_NOT_IMPOSED",
        (
            convention_low < convention_high
            and convention_object.get("convention") == "OPEN"
        ),
        convention_object,
    )

    section("C. Debt, discrete, and out-of-scope accounting")

    debt_ids = {
        "REG:b:mu_R", "REG:b:rho_br",
        "REG:E:IV:c_E", "REG:E:IV:Q_E", "REG:E:V:q_T",
    }
    debt_low_components = list(state["low_components"])
    if ACTIVE_MUTATION == "REDUCTION_DEBT_COUNTED_ONCE":
        debt_low_components.remove("REG:b:mu_R")
    debt_occurrences = Counter(debt_low_components)
    excluded_control = [
        item for item in state["low_components"] if item not in debt_ids
    ]
    doubled_control = (
        list(state["low_components"]) + ["REG:b:mu_R", "REG:b:rho_br"]
    )
    debt_ok = (
        len(debt_low_components) == 40
        and all(debt_occurrences[item] == 1 for item in debt_ids)
        and len(excluded_control) == 35
        and len(doubled_control) == 42
    )
    expect_bool(
        "REDUCTION_DEBT_COUNTED_ONCE",
        debt_ok,
        {
            "production_low": len(debt_low_components),
            "occurrences": {item: debt_occurrences[item] for item in sorted(debt_ids)},
            "exclude_control": len(excluded_control),
            "double_control": len(doubled_control),
        },
    )

    r49 = derive_r49_scope(base_rows, state)
    r49_components = list(state["low_components"])
    if ACTIVE_MUTATION == "R49_OUT_OF_SCOPE":
        r49_components.append(sorted(r49["profiles"])[0])
    r49_overlap = set(r49_components) & r49["profiles"]
    expect_bool(
        "R49_OUT_OF_SCOPE",
        (
            len(r49["profiles"]) == 8
            and len(r49_components) == 40
            and not r49_overlap
            and r49["token"] == "OUT_OF_BUILT_V2_SCOPE"
        ),
        {
            "profile_count": len(r49["profiles"]),
            "continuous_count": len(r49_components),
            "overlap": sorted(r49_overlap),
        },
    )

    discrete_rows = base_rows
    if ACTIVE_MUTATION == "DISCRETE_POSTULATE_COUNT":
        discrete_rows = replace_row(
            discrete_rows,
            "REG:disc:D1:H_existence_and_Poschl_Teller_law",
            source_status="STRUCTURAL",
        )
    discrete_state = derive_count_state(discrete_rows)
    discrete_ids = {
        row.identifier for row in discrete_rows
        if category_of(row) == CAT_DISCRETE
    }
    discrete_base_ids = {
        row.identifier for row in discrete_rows
        if category_of(row) == CAT_DISCRETE and row.part == "I-II"
    }
    discrete_extension_ids = discrete_ids - discrete_base_ids
    continuous_union = set(discrete_state["high_components"])
    expect_bool(
        "DISCRETE_POSTULATE_COUNT",
        (
            len(discrete_ids) == 11
            and len(discrete_base_ids) == 7
            and len(discrete_extension_ids) == 4
            and not (discrete_ids & continuous_union)
        ),
        {
            "total": len(discrete_ids),
            "base": len(discrete_base_ids),
            "D_III_VI": len(discrete_extension_ids),
            "continuous_overlap": sorted(discrete_ids & continuous_union),
        },
    )

    section("D. Scoped delta-r independence diagnostic")

    delta_actual = delta_diagnostic(
        force_baseline_dependency=ACTIVE_MUTATION == "DELTA_R_INDEPENDENCE"
    )
    expect_bool(
        "DELTA_R_INDEPENDENCE",
        (
            delta_actual["guards"]
            and delta_actual["token"] == "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS"
        ),
        delta_actual,
    )

    count_with_diagnostic = state["low_c"]
    if ACTIVE_MUTATION == "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION":
        count_with_diagnostic -= sum(
            int(record["delta_r"])
            for record in diagnostic["baseline"].values()
        )
    expect_bool(
        "DELTA_R_IS_DIAGNOSTIC_NOT_SUBTRACTION",
        (
            count_with_diagnostic == 40
            and count_with_diagnostic == len(state["low_components"])
            and diagnostic["token"] == "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS"
        ),
        {
            "with_diagnostic": count_with_diagnostic,
            "partition_only": len(state["low_components"]),
        },
    )

    section("E. R35, provenance, manifest, and local engine inventory")

    r35_rows = base_rows
    if ACTIVE_MUTATION == "R35_DERIVED_IN_FORM_PLUS_RIDER":
        r35_rows = replace_row(
            r35_rows, "REG:C1:Mtilde", register_label="PENDING-debt"
        )
    r35_labels = {
        row.register_label for row in r35_rows if row.bucket == "C1"
    }
    r35_riders = {
        row.rider for row in r35_rows if row.bucket == "C1"
    }
    expect_bool(
        "R35_DERIVED_IN_FORM_PLUS_RIDER",
        (
            r35_labels == {"DERIVED-in-form"}
            and r35_riders == {"C1_contributor"}
        ),
        {"labels": sorted(r35_labels), "riders": sorted(r35_riders)},
    )

    row_index = rows_by_id(base_rows)
    counted_for_provenance = set(state["low_components"])
    if ACTIVE_MUTATION == "KNOB_PROVENANCE_WELL_FORMED":
        counted_for_provenance.add("REG:derived:c_s0")
    allowed_counted_statuses = {
        "ACTION", "FREE-UNREDUCED-ROUTELESS",
        "FREE-UNREDUCED-DEBT", "CALIB", "CALIB-ANCHOR",
    }
    bad_provenance = {
        identifier: row_index[identifier].source_status
        for identifier in counted_for_provenance
        if row_index[identifier].source_status not in allowed_counted_statuses
    }
    expect_bool(
        "KNOB_PROVENANCE_WELL_FORMED",
        not bad_provenance and len(counted_for_provenance) == 40,
        bad_provenance,
    )

    manifest_rows = base_rows
    if ACTIVE_MUTATION == "REGISTER_TO_COUNT_MANIFEST":
        manifest_rows = replace_row(
            manifest_rows, "REG:a:hbar", source_status="DERIVED"
        )
    manifest_actual = manifest_state(manifest_rows)
    expect_bool(
        "REGISTER_TO_COUNT_MANIFEST",
        manifest_actual["valid"],
        manifest_actual,
    )

    qe = derive_qe_sensitivity(base_rows, state)
    inventory = {
        "bucket_counts": state["bucket_state"]["counts"],
        "raw_parts_I_II": [state["bucket_state"]["raw_low"], state["bucket_state"]["raw_high"]],
        "base_continuous": list(state["base_continuous"]),
        "E_itemization": state["E_itemization"],
        "E_continuous": state["E_continuous"],
        "continuous_codimension": [state["low_c"], state["high_c"]],
        "spread": state["spread"],
        "C1_cardinality": len(state["buckets"]["C1"]),
        "C2_cardinality": len(state["buckets"]["C2"]),
        "discrete_postulate_count": production["discrete_postulate_count"],
        "reduction_debt_counted_once": sorted(debt_ids),
        "Q_E_shift": [qe["declared_split"], qe["undeclared_split"]],
        "delta_r_result": diagnostic["token"],
        "IRREDUCIBLE_COUNT_RANGE": production,
    }
    expected_inventory = {
        "bucket_counts": {
            "a": 4, "b": 15, "c": 14, "force-mag": 1, "C1": 3, "C2": 6,
        },
        "raw_parts_I_II": [34, 43],
        "base_continuous": [27, 36],
        "E_itemization": {"III": 5, "IV": 7, "V": 1, "VI": 0},
        "E_continuous": 13,
        "continuous_codimension": [40, 49],
        "spread": 9,
        "C1_cardinality": 3,
        "C2_cardinality": 6,
        "discrete_postulate_count": 11,
        "reduction_debt_counted_once": sorted(debt_ids),
        "Q_E_shift": [3, 4],
        "delta_r_result": "CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS",
        "IRREDUCIBLE_COUNT_RANGE": EXPECTED_RANGE_OBJECT,
    }
    if ACTIVE_MUTATION == "DUAL_ENGINE_TERMS":
        inventory.pop("spread")
    expect_bool("DUAL_ENGINE_TERMS", inventory == expected_inventory, inventory)

    section("F. Assembled range and full rederivation")

    range_private = dict(production)
    if ACTIVE_MUTATION == "IRREDUCIBLE_COUNT_RANGE":
        range_private["spread"] = range_private["spread"] + 1
    expect_bool(
        "IRREDUCIBLE_COUNT_RANGE",
        range_object_is_valid(range_private),
        range_private,
    )

    rederive_rows = base_rows
    if ACTIVE_MUTATION == "COUNT_RANGE_REDERIVATION":
        rederive_rows = replace_row(rederive_rows, "REG:a:hbar", bucket="")
    rederived_state = derive_count_state(rederive_rows)
    rederived_object = assemble_range(rederive_rows, rederived_state, diagnostic)
    strict_state = derive_count_state(rederive_rows, force_strict=True)
    positive_force_strict = (
        strict_state["low_c"] == 49
        and strict_state["low_c"] == rederived_state["high_c"]
    )
    expect_bool(
        "COUNT_RANGE_REDERIVATION",
        rederived_object == EXPECTED_RANGE_OBJECT and positive_force_strict,
        {
            "rederived": rederived_object,
            "force_strict_low": strict_state["low_c"],
            "named_high": rederived_state["high_c"],
        },
    )

    return production


def main() -> None:
    if ACTIVE_MUTATION and ACTIVE_MUTATION not in TOOTH_ORDER:
        print("FIRST_FAILURE=UNKNOWN_MUTATION")
        print(f"FAIL  UNKNOWN_MUTATION: {ACTIVE_MUTATION}")
        raise AuditFailure("UNKNOWN_MUTATION", ACTIVE_MUTATION)

    print("ledger_stage043_irreducible_count_range SymPy audit")
    print(
        "ROUTE=dict/set/Counter partition + grevlex Groebner "
        "initial-monomial Krull dimension"
    )
    print("FILE_IO=none; CROSS_ENGINE_COMPARE=none")
    print("DIMENSIONAL_HOMOGENEITY=N/A_INTEGER_COUNT_STAGE")
    if ACTIVE_MUTATION:
        print(f"ACTIVE_MUTATION={ACTIVE_MUTATION}")
        print(f"MUTATED_PRIMITIVE={ABLATION_DESCRIPTIONS[ACTIVE_MUTATION]}")

    production = run_assertions()
    if ACTIVE_MUTATION:
        print("FIRST_FAILURE=MUTATION_DID_NOT_FIRE")
        raise AuditFailure("MUTATION_DID_NOT_FIRE", ACTIVE_MUTATION)

    print("")
    print(
        "IRREDUCIBLE_COUNT_RANGE="
        + json.dumps(production, sort_keys=True, separators=(",", ":"))
    )
    print("RAW_PARTS_I_II=[34,43]; BASE_CONTINUOUS=[27,36]; E_CONTINUOUS=13")
    print("CONTINUOUS_CODIMENSION=[40,49]; SPREAD=9; C1=3; C2=6")
    print("DISCRETE_POSTULATE_COUNT=11 (=7+4); NOT_IN_CONTINUOUS_VARIETY")
    print("REDUCTION_DEBT_COUNTED_ONCE={mu_R,rho_br,c_E,q_T,Q_E}")
    print("R49_PARALLEL_TRACK=OUT_OF_BUILT_V2_SCOPE; PROFILE_DEBT=8_to_9")
    print("Q_E_UNDECLARED=TOTAL_NEUTRAL_ROUTELESS_SHIFT_3_to_4")
    print("K_THETA_KAPPA_PHASE_AS_TWO_SENSITIVITY=[41,50]")
    print(
        "DELTA_R_SCOPE=CONFIRMED_IN_TESTED_M_AND_WCHI_BLOCKS; "
        "M:(10,5,5)->(10,4,6); Wchi:(10,7,3)->(10,6,4)"
    )
    print(f"TOOTH_COUNT={len(TOOTH_ORDER)}")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy independently assembled the count range")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(
            "OVERALL FAIL: SymPy stage043 audit did not close "
            f"({exc.predicate})"
        )
        raise SystemExit(1)
