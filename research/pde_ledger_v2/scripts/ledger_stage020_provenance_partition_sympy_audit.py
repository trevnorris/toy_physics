#!/usr/bin/env python3
"""Ledger stage020 SymPy audit: provenance partition + CALIBRATED label.

Standalone, print-only, no arguments, no file I/O.  This is the pathA_33
II-G4c algebra/provenance slice only: the strengthened a^-5 bridge scaling,
the Gamma5/chi_Q equivalence, the expression-bound 54/5 = 2*27/5 identity,
the four-way provenance classifier, and the provenance-driven local verdict.
The adjacent fingerprint, prefactor, and dimensional-closure stages are cited
as provenance only.
"""

from __future__ import annotations

import inspect
from typing import Any, Iterable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

PROVENANCE_PARTITION_CALIBRATED = "PROVENANCE_PARTITION_CALIBRATED"
QUAD_PASS = "QUAD_PASS"
QUAD_CALIBRATED = "QUAD_CALIBRATED"
FAIL_SCALING = "FAIL_SCALING"
FAIL_EQUIVALENCE = "FAIL_EQUIVALENCE"
FAIL_PROVENANCE_PARTITION = "FAIL_PROVENANCE_PARTITION"
NO_FAIL = "NO_FAIL"

DERIVED_IN_GATE = "derived_in_gate"
EXTERNAL_BRIDGE_INPUT = "external_bridge_input"
DEFERRED_BRANCH_DATA = "deferred_branch_data"
CONVENTION = "convention"
UNCLASSIFIED = "unclassified"


class AuditFailure(AssertionError):
    pass


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(expr: Any) -> Any:
    if isinstance(expr, sp.MatrixBase):
        return expr.applyfunc(compact)
    if not isinstance(expr, sp.Basic):
        return expr
    if isinstance(expr, (sp.Equality, sp.Order)):
        return expr
    return sp.factor(sp.cancel(sp.together(sp.simplify(expr))))


def fmt(expr: Any) -> str:
    if isinstance(expr, bool):
        return "True" if expr else "False"
    if isinstance(expr, str):
        return expr
    try:
        return sp.sstr(compact(expr))
    except Exception:
        return sp.sstr(expr)


def assert_no_float(name: str, expr: Any) -> None:
    if isinstance(expr, dict):
        for key, value in expr.items():
            assert_no_float(f"{name}.{key}", value)
        return
    if isinstance(expr, sp.MatrixBase):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
        return
    if isinstance(expr, (list, tuple, set, frozenset)):
        for index, value in enumerate(expr):
            assert_no_float(f"{name}[{index}]", value)
        return
    if isinstance(expr, (str, type(None))):
        return
    if isinstance(expr, bool):
        expr = sp.Integer(1) if expr else sp.Integer(0)
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
    if floats:
        message = f"{name}: Float atom(s) found in exact audit expression: {floats}"
        _record_fail(f"FAIL  {message}")
        raise AuditFailure(message)


def _record_pass(message: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(message)


def _record_fail(message: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(message)


def expect_zero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_pass(f"PASS  {name}")
        return
    _record_fail(f"FAIL  {name}: residual = {fmt(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if bool(condition) else sp.Integer(1))


def expect_fail(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = compact(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {fmt(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def bool_from_residual(residual: sp.Expr | int) -> bool:
    assert_no_float("bool_from_residual", residual)
    return compact(residual) == 0


def verdict_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


def class_residual(actual: str, expected: str) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


a, c_s, c, G = sp.symbols("a c_s c G", positive=True)
N0, D0, chi_Q, lambda_G = sp.symbols(
    "N0 D0 chi_Q lambda_G", real=True, nonzero=True
)


DERIVED_TAGS = frozenset(
    {
        "dtn_hankel_expansion",
        "dtn_radiative_slot",
        "prefactor_series_algebra",
        "target_rhs_algebra",
        "gamma_bridge_algebra",
        "emergent_outgoing_passivity",
    }
)
EXTERNAL_TAGS = frozenset(
    {"external_gr_constant", "external_pn_bridge", "einstein_bridge_identity"}
)
DEFERRED_TAGS = frozenset({"gate6_branch_solve", "deferred_nonlinear_pde"})
CONVENTION_TAGS = frozenset(
    {"normalization_convention", "unit_choice", "static_slot_convention"}
)


PROVENANCE_SOURCES: dict[str, list[str]] = {
    "fingerprint_u2": ["dtn_hankel_expansion"],
    "fingerprint_u4": ["dtn_hankel_expansion"],
    "fingerprint_v5": ["dtn_hankel_expansion", "dtn_radiative_slot"],
    "fingerprint_27": ["dtn_radiative_slot"],
    "prefactor_P0_P2_P4": ["prefactor_series_algebra"],
    "P0_target_scaling_minus5": ["target_rhs_algebra"],
    "chi_Q": ["dtn_radiative_slot"],
    "Gamma5_equivalence_chain": ["gamma_bridge_algebra"],
    "emergent_passivity": ["emergent_outgoing_passivity"],
    "G": ["external_gr_constant"],
    "PN_2_over_5": ["external_pn_bridge"],
    "Einstein_2G_over_5c5_identity": ["einstein_bridge_identity"],
    "assembled_54_over_5_magnitude": [
        "external_pn_bridge",
        "dtn_radiative_slot",
    ],
    "D_n_N_n_numeric_values": ["gate6_branch_solve"],
    "port_scalars": ["gate6_branch_solve"],
    "actual_branch_a_scaling": ["gate6_branch_solve"],
    "unit_choices": ["unit_choice"],
    "static_slot_minus3": ["static_slot_convention"],
}


def a_power(expr: sp.Expr) -> sp.Rational | None:
    powers = sp.factor(expr).as_powers_dict()
    power = powers.get(a, sp.Rational(0))
    return sp.Rational(power) if power.is_number else None


def classify_provenance(tags: Iterable[str]) -> str:
    """Four-way source classifier with source-faithful dominance."""

    tagset = set(tags)
    if tagset & DEFERRED_TAGS:
        return DEFERRED_BRANCH_DATA
    if tagset & EXTERNAL_TAGS:
        return EXTERNAL_BRIDGE_INPUT
    if tagset & DERIVED_TAGS:
        return DERIVED_IN_GATE
    if tagset & CONVENTION_TAGS:
        return CONVENTION
    return UNCLASSIFIED


def group_partition(
    items: dict[str, dict[str, Any]], class_field: str
) -> dict[str, list[str]]:
    groups = {
        DERIVED_IN_GATE: [],
        EXTERNAL_BRIDGE_INPUT: [],
        DEFERRED_BRANCH_DATA: [],
        CONVENTION: [],
        UNCLASSIFIED: [],
    }
    for name, item in items.items():
        groups.setdefault(item[class_field], []).append(name)
    return groups


def copy_sources(
    sources: dict[str, list[str]] = PROVENANCE_SOURCES,
) -> dict[str, list[str]]:
    return {name: list(tags) for name, tags in sources.items()}


def build_partition(
    overrides: dict[str, str] | None = None,
    sources: dict[str, list[str]] | None = None,
) -> dict[str, Any]:
    source_map = copy_sources(sources or PROVENANCE_SOURCES)
    emitted_overrides = overrides or {}
    items: dict[str, dict[str, Any]] = {}
    for name, tags in source_map.items():
        computed = classify_provenance(tags)
        emitted = emitted_overrides.get(name, computed)
        items[name] = {
            "provenance_tags": tags,
            "computed_class": computed,
            "emitted_class": emitted,
            "class_matches_computed": emitted == computed,
        }
    return {
        "items": items,
        "groups": group_partition(items, "computed_class"),
        "emitted_groups": group_partition(items, "emitted_class"),
        "ok": all(item["class_matches_computed"] for item in items.values()),
    }


def build_scaling() -> dict[str, Any]:
    # The frozen 018 radiative slot is the sole a^5 carrier in the derived path.
    v5_slot = a**5 / (sp.Integer(27) * c_s**5)
    chi = sp.Integer(1)
    gamma_target = sp.Integer(2) * G / (sp.Integer(5) * c**5)
    target_from_bridge = compact(gamma_target / (chi * v5_slot))

    # This independently assembled target is the object compared to the bridge.
    target_rhs = sp.Integer(54) * G * c_s**5 / (
        sp.Integer(5) * a**5 * c**5
    )
    mutated_target_rhs = sp.Integer(54) * G * c_s**5 / (
        sp.Integer(5) * a**4 * c**5
    )

    gamma_power = a_power(gamma_target)
    slot_power = a_power(v5_slot)
    derived_power = (
        gamma_power - slot_power
        if gamma_power is not None and slot_power is not None
        else None
    )
    bridge_power = a_power(target_from_bridge)
    assembled_bridge_residual = compact(target_from_bridge - target_rhs)
    mutated_bridge_residual = compact(target_from_bridge - mutated_target_rhs)
    ok = (
        gamma_power == 0
        and slot_power == 5
        and derived_power == -5
        and bridge_power == derived_power
        and bool_from_residual(assembled_bridge_residual)
    )
    return {
        "chi": chi,
        "v5_slot": v5_slot,
        "gamma_target": gamma_target,
        "target_from_bridge": target_from_bridge,
        "target_rhs": target_rhs,
        "mutated_target_rhs": mutated_target_rhs,
        "gamma_power": gamma_power,
        "slot_power": slot_power,
        "derived_power": derived_power,
        "bridge_power": bridge_power,
        "assembled_bridge_residual": assembled_bridge_residual,
        "mutated_bridge_residual": mutated_bridge_residual,
        "mutated_ok": bool_from_residual(mutated_bridge_residual),
        "ok": ok,
    }


def build_equivalence(scaling: dict[str, Any]) -> dict[str, Any]:
    chi = scaling["chi"]
    v5_slot = scaling["v5_slot"]
    target_rhs = scaling["target_rhs"]
    gamma_target = scaling["gamma_target"]

    p0_physical = (c_s / a) ** 2 * (N0 / D0)
    gamma5 = chi_Q * p0_physical * v5_slot
    gamma5_with_outgoing_chi = compact(gamma5.subs(chi_Q, chi))

    forward_general = compact(target_rhs * chi_Q * v5_slot - gamma_target)
    forward_expected = sp.Integer(2) * G * (chi_Q - 1) / (
        sp.Integer(5) * c**5
    )
    forward = compact(forward_general.subs(chi_Q, chi))
    reverse = compact(gamma_target / (chi * v5_slot) - target_rhs)
    wrong_gamma = sp.Integer(3) * G / (sp.Integer(5) * c**5)
    wrong_reverse = compact(wrong_gamma / (chi * v5_slot) - target_rhs)
    ok = bool_from_residual(forward) and bool_from_residual(reverse)
    return {
        "P0_physical": p0_physical,
        "Gamma5": gamma5,
        "Gamma5_with_outgoing_chi": gamma5_with_outgoing_chi,
        "forward_general": forward_general,
        "forward_expected": forward_expected,
        "forward": forward,
        "reverse": reverse,
        "wrong_gamma": wrong_gamma,
        "wrong_reverse": wrong_reverse,
        "wrong_gamma_probe_fires": not bool_from_residual(wrong_reverse),
        "ok": ok,
    }


def build_bound_identity(
    scaling: dict[str, Any], partition: dict[str, Any]
) -> dict[str, Any]:
    target_rhs = scaling["target_rhs"]
    v5_slot = scaling["v5_slot"]
    unit = G * c_s**5 / (a**5 * c**5)
    mag = compact(target_rhs / unit)
    twenty_seven_from_slot = compact(a**5 / (c_s**5 * v5_slot))
    right = compact(
        sp.Integer(2) * twenty_seven_from_slot / sp.Integer(5)
    )
    residual = compact(mag - right)
    mutated_right = compact(
        sp.Integer(2) * (twenty_seven_from_slot - 1) / sp.Integer(5)
    )
    mutated_residual = compact(mag - mutated_right)
    items = partition["items"]
    return {
        "unit": unit,
        "mag": mag,
        "twenty_seven_from_slot": twenty_seven_from_slot,
        "left": mag,
        "right": right,
        "residual": residual,
        "mutated_right": mutated_right,
        "mutated_residual": mutated_residual,
        "earned_factor_class": items["fingerprint_27"]["computed_class"],
        "calibrated_factor_class": items["PN_2_over_5"]["computed_class"],
        "assembled_magnitude_class": items[
            "assembled_54_over_5_magnitude"
        ]["computed_class"],
        "ok": bool_from_residual(residual),
    }


def build_classifier_proof(partition: dict[str, Any]) -> dict[str, Any]:
    truth_inputs = {
        "deferred_over_external": [
            "gate6_branch_solve",
            "external_gr_constant",
        ],
        "external_over_derived": [
            "external_pn_bridge",
            "dtn_hankel_expansion",
        ],
        "derived_over_convention": ["dtn_radiative_slot", "unit_choice"],
        "single_deferred": ["gate6_branch_solve"],
        "single_external": ["external_gr_constant"],
        "single_derived": ["dtn_radiative_slot"],
        "single_convention": ["unit_choice"],
    }
    truth_expected = {
        "deferred_over_external": DEFERRED_BRANCH_DATA,
        "external_over_derived": EXTERNAL_BRIDGE_INPUT,
        "derived_over_convention": DERIVED_IN_GATE,
        "single_deferred": DEFERRED_BRANCH_DATA,
        "single_external": EXTERNAL_BRIDGE_INPUT,
        "single_derived": DERIVED_IN_GATE,
        "single_convention": CONVENTION,
    }
    truth_results = {
        name: classify_provenance(tags) for name, tags in truth_inputs.items()
    }
    truth_ok = all(
        truth_results[name] == truth_expected[name] for name in truth_inputs
    )

    items = partition["items"]
    key_expected = {
        "G": EXTERNAL_BRIDGE_INPUT,
        "PN_2_over_5": EXTERNAL_BRIDGE_INPUT,
        "fingerprint_27": DERIVED_IN_GATE,
        "assembled_54_over_5_magnitude": EXTERNAL_BRIDGE_INPUT,
    }
    key_results = {
        name: items[name]["computed_class"] for name in key_expected
    }
    key_ok = all(key_results[name] == expected for name, expected in key_expected.items())

    stripped_tags = [
        tag
        for tag in PROVENANCE_SOURCES["assembled_54_over_5_magnitude"]
        if tag != "external_pn_bridge"
    ]
    baseline_computed_class = classify_provenance(
        PROVENANCE_SOURCES["assembled_54_over_5_magnitude"]
    )
    mutated_computed_class = classify_provenance(stripped_tags)
    mutation_fires = (
        baseline_computed_class == EXTERNAL_BRIDGE_INPUT
        and mutated_computed_class == DERIVED_IN_GATE
        and mutated_computed_class != baseline_computed_class
    )
    return {
        "truth_inputs": truth_inputs,
        "truth_expected": truth_expected,
        "truth_results": truth_results,
        "truth_ok": truth_ok,
        "key_expected": key_expected,
        "key_results": key_results,
        "key_ok": key_ok,
        "tag_mutation": {
            "baseline_tags": list(
                PROVENANCE_SOURCES["assembled_54_over_5_magnitude"]
            ),
            "stripped_tags": stripped_tags,
            "baseline_computed_class": baseline_computed_class,
            "mutated_computed_class": mutated_computed_class,
            "fires": mutation_fires,
        },
        "ok": truth_ok and key_ok and mutation_fires,
    }


def is_g_invariant(expr: sp.Expr) -> bool:
    return bool_from_residual(expr.subs(G, lambda_G * G) - expr)


def build_g_invariance_diagnostic(
    partition: dict[str, Any],
    scaling: dict[str, Any],
    identity: dict[str, Any],
) -> dict[str, Any]:
    items = partition["items"]
    raw = {
        "G": (G, items["G"]["computed_class"]),
        "target_2G_over_5c5": (
            scaling["gamma_target"],
            items["PN_2_over_5"]["computed_class"],
        ),
        "pure_54_over_5": (
            identity["mag"],
            items["assembled_54_over_5_magnitude"]["computed_class"],
        ),
        "fingerprint_27": (
            identity["twenty_seven_from_slot"],
            items["fingerprint_27"]["computed_class"],
        ),
    }
    diagnostics: dict[str, dict[str, Any]] = {}
    for name, (expr, provenance_class) in raw.items():
        transformed = compact(expr.subs(G, lambda_G * G))
        residual = compact(transformed - expr)
        diagnostics[name] = {
            "expression": expr,
            "transformed": transformed,
            "residual": residual,
            "g_invariant": bool_from_residual(residual),
            "provenance_class": provenance_class,
        }
    trap = (
        diagnostics["pure_54_over_5"]["g_invariant"]
        and diagnostics["pure_54_over_5"]["provenance_class"]
        == EXTERNAL_BRIDGE_INPUT
    )
    return {
        "diagnostics": diagnostics,
        "classifier": "provenance_not_g_invariance",
        "invariance_only_trap_catches_54_over_5": trap,
    }


def verdict_from_partition(
    scaling_ok: bool,
    equivalence_ok: bool,
    provenance_ok: bool,
    partition: dict[str, Any],
) -> str:
    """Pure 020-local verdict: three gates, then provenance classes."""

    if not scaling_ok:
        return FAIL_SCALING
    if not equivalence_ok:
        return FAIL_EQUIVALENCE
    if not provenance_ok:
        return FAIL_PROVENANCE_PARTITION
    items = partition["items"]
    g_class = items["G"]["computed_class"]
    mag_class = items["assembled_54_over_5_magnitude"]["computed_class"]
    if g_class == DERIVED_IN_GATE and mag_class == DERIVED_IN_GATE:
        return QUAD_PASS
    return QUAD_CALIBRATED


def local_verdict(gates: dict[str, bool], partition: dict[str, Any]) -> str:
    return verdict_from_partition(
        gates["scaling_ok"],
        gates["equivalence_ok"],
        gates["provenance_ok"],
        partition,
    )


def dynamic_ablation(
    baseline_gates: dict[str, bool],
    baseline_partition: dict[str, Any],
    mutated_gates: dict[str, bool],
    mutated_partition: dict[str, Any],
    expected_fail: str,
) -> dict[str, Any]:
    """Rerun the pure 020-local verdict on computed baseline/mutant inputs."""

    with_mutation = local_verdict(mutated_gates, mutated_partition)
    without_mutation = local_verdict(baseline_gates, baseline_partition)
    return {
        "rerun_gate_logic": True,
        "with_mutation": with_mutation,
        "without_mutation": without_mutation,
        "expected_fail": expected_fail,
        "fail_suppressed": (
            with_mutation == expected_fail
            and without_mutation != expected_fail
            and with_mutation != without_mutation
        ),
    }


def control_partition(g_class: str, mag_class: str) -> dict[str, Any]:
    sources = copy_sources()
    sources["G"] = [
        "dtn_radiative_slot"
        if g_class == DERIVED_IN_GATE
        else "external_gr_constant"
    ]
    sources["assembled_54_over_5_magnitude"] = [
        "dtn_radiative_slot"
        if mag_class == DERIVED_IN_GATE
        else "external_pn_bridge"
    ]
    return build_partition(sources=sources)


def inverted_verdict_control(
    scaling_ok: bool,
    equivalence_ok: bool,
    provenance_ok: bool,
    partition: dict[str, Any],
) -> str:
    """Deliberately wrong PASS-unless-both-external mutation."""

    if not scaling_ok:
        return FAIL_SCALING
    if not equivalence_ok:
        return FAIL_EQUIVALENCE
    if not provenance_ok:
        return FAIL_PROVENANCE_PARTITION
    items = partition["items"]
    both_external = (
        items["G"]["computed_class"] == EXTERNAL_BRIDGE_INPUT
        and items["assembled_54_over_5_magnitude"]["computed_class"]
        == EXTERNAL_BRIDGE_INPUT
    )
    return QUAD_CALIBRATED if both_external else QUAD_PASS


def build_baseline() -> dict[str, Any]:
    scaling = build_scaling()
    equivalence = build_equivalence(scaling)
    partition = build_partition()
    identity = build_bound_identity(scaling, partition)
    classifier_proof = build_classifier_proof(partition)
    g_invariance = build_g_invariance_diagnostic(partition, scaling, identity)

    factor_classes_ok = (
        identity["earned_factor_class"] == DERIVED_IN_GATE
        and identity["calibrated_factor_class"] == EXTERNAL_BRIDGE_INPUT
        and identity["assembled_magnitude_class"] == EXTERNAL_BRIDGE_INPUT
    )
    provenance_ok = (
        partition["ok"]
        and identity["ok"]
        and classifier_proof["ok"]
        and factor_classes_ok
    )
    gates = {
        "scaling_ok": scaling["ok"],
        "equivalence_ok": equivalence["ok"],
        "provenance_ok": provenance_ok,
    }

    all_derived = control_partition(DERIVED_IN_GATE, DERIVED_IN_GATE)
    mixed = control_partition(EXTERNAL_BRIDGE_INPUT, DERIVED_IN_GATE)
    controls = {
        "all_derived_partition": all_derived,
        "mixed_partition": mixed,
        "all_derived_verdict": verdict_from_partition(
            True, True, all_derived["ok"], all_derived
        ),
        "mixed_verdict": verdict_from_partition(True, True, mixed["ok"], mixed),
        "inverted_mixed_verdict": inverted_verdict_control(
            True, True, mixed["ok"], mixed
        ),
        "constant_calibrated_on_all_derived": QUAD_CALIBRATED,
    }
    data = {
        "scaling": scaling,
        "equivalence": equivalence,
        "partition": partition,
        "identity": identity,
        "classifier_proof": classifier_proof,
        "g_invariance": g_invariance,
        "gates": gates,
        "verdict": local_verdict(gates, partition),
        "verdict_controls": controls,
    }
    data["probes"] = build_probes(data)
    return data


def build_probes(data: dict[str, Any]) -> dict[str, Any]:
    baseline_gates = data["gates"]
    baseline_partition = data["partition"]
    scaling = data["scaling"]
    equivalence = data["equivalence"]

    scaling_mutated_gates = dict(baseline_gates)
    scaling_mutated_gates["scaling_ok"] = scaling["mutated_ok"]
    scaling_ablation = dynamic_ablation(
        baseline_gates,
        baseline_partition,
        scaling_mutated_gates,
        baseline_partition,
        FAIL_SCALING,
    )

    mutated_equivalence_ok = bool_from_residual(equivalence["wrong_reverse"])
    equivalence_mutated_gates = dict(baseline_gates)
    equivalence_mutated_gates["equivalence_ok"] = mutated_equivalence_ok
    equivalence_ablation = dynamic_ablation(
        baseline_gates,
        baseline_partition,
        equivalence_mutated_gates,
        baseline_partition,
        FAIL_EQUIVALENCE,
    )

    g_as_derived = build_partition({"G": DERIVED_IN_GATE})
    fingerprint_as_external = build_partition(
        {"fingerprint_27": EXTERNAL_BRIDGE_INPUT}
    )

    def partition_ablation(mutated_partition: dict[str, Any]) -> dict[str, Any]:
        mutated_gates = dict(baseline_gates)
        mutated_gates["provenance_ok"] = mutated_partition["ok"]
        return dynamic_ablation(
            baseline_gates,
            baseline_partition,
            mutated_gates,
            mutated_partition,
            FAIL_PROVENANCE_PARTITION,
        )

    g_ablation = partition_ablation(g_as_derived)
    fingerprint_ablation = partition_ablation(fingerprint_as_external)
    return {
        "3c_wrong_scaling": {
            "mutated_residual": scaling["mutated_bridge_residual"],
            "mutated_ok": scaling["mutated_ok"],
            "verdict": FAIL_SCALING if not scaling["mutated_ok"] else NO_FAIL,
            "correct_object_verdict": NO_FAIL if scaling["ok"] else FAIL_SCALING,
            "self_ablation": scaling_ablation,
        },
        "3e_equivalence_break": {
            "wrong_gamma": equivalence["wrong_gamma"],
            "wrong_reverse": equivalence["wrong_reverse"],
            "mutated_ok": mutated_equivalence_ok,
            "verdict": (
                FAIL_EQUIVALENCE if not mutated_equivalence_ok else NO_FAIL
            ),
            "correct_object_verdict": (
                NO_FAIL if equivalence["ok"] else FAIL_EQUIVALENCE
            ),
            "self_ablation": equivalence_ablation,
        },
        "3f_partition_mislabel": {
            "G_external_to_derived": {
                "partition": g_as_derived,
                "verdict": (
                    FAIL_PROVENANCE_PARTITION if not g_as_derived["ok"] else NO_FAIL
                ),
                "self_ablation": g_ablation,
            },
            "fingerprint_27_derived_to_external": {
                "partition": fingerprint_as_external,
                "verdict": (
                    FAIL_PROVENANCE_PARTITION
                    if not fingerprint_as_external["ok"]
                    else NO_FAIL
                ),
                "self_ablation": fingerprint_ablation,
            },
            "with_mutation_cases": {
                "G_external_to_derived": g_ablation["with_mutation"],
                "fingerprint_27_derived_to_external": fingerprint_ablation[
                    "with_mutation"
                ],
            },
            "g_invariance_only_would_miss_54_over_5": data["g_invariance"][
                "invariance_only_trap_catches_54_over_5"
            ],
        },
    }


def run_scaling_and_equivalence(data: dict[str, Any]) -> None:
    scaling = data["scaling"]
    equivalence = data["equivalence"]
    subbanner("Strengthened a⁻⁵ scaling from the frozen 018 v5 slot")
    print(f"  gamma_target = {fmt(scaling['gamma_target'])}")
    print(f"  a_power(gamma_target) = {fmt(scaling['gamma_power'])}")
    print(f"  v5_slot = {fmt(scaling['v5_slot'])}")
    print(f"  a_power(v5_slot) = {fmt(scaling['slot_power'])}")
    print(f"  target_from_bridge = {fmt(scaling['target_from_bridge'])}")
    print(
        "  derived_power = a_power(gamma_target) - a_power(v5_slot) = "
        + fmt(scaling["derived_power"])
    )
    print(f"  a_power(target_from_bridge) = {fmt(scaling['bridge_power'])} (a⁻⁵)")
    print(f"  independently assembled target_rhs = {fmt(scaling['target_rhs'])}")
    expect_zero("Burke-Thorne gamma_target is a-free", scaling["gamma_power"])
    expect_zero("frozen v5_slot supplies a-power +5", scaling["slot_power"] - 5)
    expect_zero("derived bridge power is -5", scaling["derived_power"] + 5)
    expect_zero(
        "target_from_bridge power equals the derived power",
        scaling["bridge_power"] - scaling["derived_power"],
    )
    expect_zero(
        "bridge reconstruction equals the independently assembled target_rhs",
        scaling["assembled_bridge_residual"],
    )
    expect_fail(
        "3c a^-4 mutation of only the assembled target breaks bridge cancellation",
        scaling["mutated_bridge_residual"],
    )
    expect_bool("strengthened scaling gate passes", scaling["ok"])

    subbanner("Gamma5/chi_Q equivalence bridge")
    print(f"  P0_physical (019 provenance definition only) = {fmt(equivalence['P0_physical'])}")
    print(f"  Gamma5 DEF = {fmt(equivalence['Gamma5'])}")
    print(f"  forward_general = {fmt(equivalence['forward_general'])}")
    print(f"  forward = {fmt(equivalence['forward'])}")
    print(f"  reverse = {fmt(equivalence['reverse'])}")
    print(f"  wrong_reverse (3/5 mutation) = {fmt(equivalence['wrong_reverse'])}")
    p0_symbols = {N0, D0}
    bridge_residual_symbols = (
        equivalence["forward"].free_symbols
        | equivalence["reverse"].free_symbols
        | equivalence["forward_general"].free_symbols
    )
    expect_bool(
        "P0=N0/D0 enters the Gamma5 DEF and is absent from bridge residuals",
        p0_symbols <= equivalence["Gamma5"].free_symbols
        and p0_symbols.isdisjoint(bridge_residual_symbols),
    )
    expect_zero(
        "forward general form is 2*G*(chi_Q-1)/(5*c^5)",
        equivalence["forward_general"] - equivalence["forward_expected"],
    )
    expect_zero("Gamma5 bridge forward=0 for cited chi_Q=+1", equivalence["forward"])
    expect_zero("Gamma5 bridge reverse=0", equivalence["reverse"])
    expect_fail(
        "3e wrong_gamma=3*G/(5*c^5) makes wrong_reverse nonzero",
        equivalence["wrong_reverse"],
    )
    expect_bool("equivalence gate passes", equivalence["ok"])


def run_identity_and_classifier(data: dict[str, Any]) -> None:
    identity = data["identity"]
    proof = data["classifier_proof"]
    partition = data["partition"]
    subbanner("Expression-bound 54/5 = 2·27/5 decomposition")
    print(f"  unit extracted from target_rhs = {fmt(identity['unit'])}")
    print(f"  mag = target_rhs/unit = {fmt(identity['mag'])}")
    print(
        "  27_from_slot = a^5/(c_s^5*v5_slot) = "
        + fmt(identity["twenty_seven_from_slot"])
    )
    print(f"  bound left = {fmt(identity['left'])}")
    print(f"  bound right = 2*27_from_slot/5 = {fmt(identity['right'])}")
    print(f"  54/5=2·27/5 TRUE bound residual = {fmt(identity['residual'])}")
    expect_zero(
        "BOUND 54/5=2*27/5 identity from target_rhs and v5_slot",
        identity["residual"],
    )
    expect_fail(
        "bound-identity factor mutation 27_from_slot->27_from_slot-1 fires",
        identity["mutated_residual"],
    )

    subbanner("Four-way provenance classifier proof")
    print("  dominance order = deferred > external > derived > convention")
    for name, tags in proof["truth_inputs"].items():
        result = proof["truth_results"][name]
        expected = proof["truth_expected"][name]
        print(f"  truth-table {name}: {tags} -> {result}")
        expect_zero(
            f"classifier truth-table {name} -> {expected}",
            class_residual(result, expected),
        )
    for name, expected in proof["key_expected"].items():
        result = proof["key_results"][name]
        print(f"  key class {name} -> {result}")
        expect_zero(
            f"key provenance class {name} -> {expected}",
            class_residual(result, expected),
        )
    mutation = proof["tag_mutation"]
    print(
        "  tag mutation assembled magnitude: "
        f"{mutation['baseline_tags']} -> {mutation['stripped_tags']} -> "
        f"{mutation['mutated_computed_class']}"
    )
    expect_zero(
        "stripping external_pn_bridge drops assembled magnitude to derived_in_gate",
        class_residual(mutation["mutated_computed_class"], DERIVED_IN_GATE),
    )
    expect_zero(
        "baseline assembled magnitude class is external_bridge_input before tag stripping",
        class_residual(
            mutation["baseline_computed_class"], EXTERNAL_BRIDGE_INPUT
        ),
    )
    expect_bool("classifier truth-table + key classes + tag mutation all pass", proof["ok"])
    expect_bool("real emitted partition matches computed classes", partition["ok"])
    print(f"  earned factor class (READ from partition) = {identity['earned_factor_class']}")
    print(f"  calibrated factor class (READ from partition) = {identity['calibrated_factor_class']}")
    print(f"  assembled class (READ from partition) = {identity['assembled_magnitude_class']}")
    expect_zero(
        "27 factor class is read as derived_in_gate",
        class_residual(identity["earned_factor_class"], DERIVED_IN_GATE),
    )
    expect_zero(
        "2/5 factor class is read as external_bridge_input",
        class_residual(
            identity["calibrated_factor_class"], EXTERNAL_BRIDGE_INPUT
        ),
    )
    expect_zero(
        "assembled 54/5 class is read as external_bridge_input",
        class_residual(
            identity["assembled_magnitude_class"], EXTERNAL_BRIDGE_INPUT
        ),
    )


def run_g_invariance_diagnostic(data: dict[str, Any]) -> None:
    diagnostic = data["g_invariance"]
    subbanner("SEPARATE non-verdict-driving G->lambda_G*G diagnostic")
    expected = {
        "G": False,
        "target_2G_over_5c5": False,
        "pure_54_over_5": True,
        "fingerprint_27": True,
    }
    for name, item in diagnostic["diagnostics"].items():
        print(
            f"  {name}: expr={fmt(item['expression'])}, "
            f"g_invariant={item['g_invariant']}, class={item['provenance_class']}"
        )
        expect_bool(
            f"G-invariance diagnostic for {name} is {expected[name]}",
            item["g_invariant"] is expected[name],
        )
    trap = diagnostic["invariance_only_trap_catches_54_over_5"]
    print(f"  invariance_only_trap_catches_54_over_5 = {trap}")
    expect_bool(
        "separate invariance-only trap catches G-invariant yet calibrated 54/5",
        trap,
    )
    expect_zero(
        "provenance verdict remains CALIBRATED independently of the diagnostic",
        verdict_residual(data["verdict"], QUAD_CALIBRATED),
    )


def run_verdict_controls(data: dict[str, Any]) -> None:
    controls = data["verdict_controls"]
    real_items = data["partition"]["items"]
    all_items = controls["all_derived_partition"]["items"]
    mixed_items = controls["mixed_partition"]["items"]
    subbanner("Source-faithful provenance verdict and positive controls")
    print(
        "  real classes: G="
        + real_items["G"]["computed_class"]
        + ", assembled54/5="
        + real_items["assembled_54_over_5_magnitude"]["computed_class"]
    )
    print(f"  real verdict = {data['verdict']}")
    print(
        "  QUAD_PASS control classes: G="
        + all_items["G"]["computed_class"]
        + ", assembled54/5="
        + all_items["assembled_54_over_5_magnitude"]["computed_class"]
        + f"; partition_ok={controls['all_derived_partition']['ok']}"
    )
    print(f"  QUAD_PASS control verdict = {controls['all_derived_verdict']}")
    print(
        "  MIXED control classes: G="
        + mixed_items["G"]["computed_class"]
        + ", assembled54/5="
        + mixed_items["assembled_54_over_5_magnitude"]["computed_class"]
        + f"; partition_ok={controls['mixed_partition']['ok']}"
    )
    print(f"  MIXED control verdict = {controls['mixed_verdict']}")
    expect_zero(
        "real both-external partition lands QUAD_CALIBRATED",
        verdict_residual(data["verdict"], QUAD_CALIBRATED),
    )
    expect_bool(
        "coherent all-derived control partition is internally consistent",
        controls["all_derived_partition"]["ok"],
    )
    expect_zero(
        "coherent all-derived positive control lands QUAD_PASS",
        verdict_residual(controls["all_derived_verdict"], QUAD_PASS),
    )
    expect_bool(
        "coherent mixed control partition is internally consistent",
        controls["mixed_partition"]["ok"],
    )
    expect_zero(
        "required MIXED control lands QUAD_CALIBRATED",
        verdict_residual(controls["mixed_verdict"], QUAD_CALIBRATED),
    )
    expect_fail(
        "constant-CALIBRATED verdict ablation fails at the QUAD_PASS control",
        verdict_residual(
            controls["constant_calibrated_on_all_derived"], QUAD_PASS
        ),
    )
    expect_fail(
        "inverted PASS-unless-both-external rule fails at the MIXED control",
        verdict_residual(controls["inverted_mixed_verdict"], QUAD_CALIBRATED),
    )


def assert_dynamic_ablation(
    label: str, ablation: dict[str, Any], expected_fail: str
) -> None:
    print(f"  {label} dynamic self-ablation = {ablation}")
    expect_bool(f"{label} reruns 020-local gate logic", ablation["rerun_gate_logic"])
    expect_zero(
        f"{label} dynamic rerun with mutation reaches {expected_fail}",
        verdict_residual(ablation["with_mutation"], expected_fail),
    )
    expect_zero(
        f"{label} dynamic rerun without mutation returns QUAD_CALIBRATED",
        verdict_residual(ablation["without_mutation"], QUAD_CALIBRATED),
    )
    expect_bool(
        f"{label} self-ablation suppresses the local failure",
        ablation["fail_suppressed"],
    )


def run_probes(data: dict[str, Any]) -> None:
    probes = data["probes"]
    subbanner("Probes 3c/3e/3f and DYNAMIC 020-local self-ablations")
    scaling = probes["3c_wrong_scaling"]
    print(f"  3c wrong-scaling residual = {fmt(scaling['mutated_residual'])}")
    print(f"  3c verdict = {scaling['verdict']}")
    expect_zero(
        "3c wrong assembled scaling reaches FAIL_SCALING",
        verdict_residual(scaling["verdict"], FAIL_SCALING),
    )
    expect_zero(
        "3c correct bridge scaling is NO_FAIL",
        verdict_residual(scaling["correct_object_verdict"], NO_FAIL),
    )
    assert_dynamic_ablation("3c scaling_ablation", scaling["self_ablation"], FAIL_SCALING)

    equivalence = probes["3e_equivalence_break"]
    print(f"  3e wrong_gamma = {fmt(equivalence['wrong_gamma'])}")
    print(f"  3e wrong_reverse = {fmt(equivalence['wrong_reverse'])}")
    print(f"  3e verdict = {equivalence['verdict']}")
    expect_zero(
        "3e wrong Gamma bridge reaches FAIL_EQUIVALENCE",
        verdict_residual(equivalence["verdict"], FAIL_EQUIVALENCE),
    )
    expect_zero(
        "3e correct Gamma bridge is NO_FAIL",
        verdict_residual(equivalence["correct_object_verdict"], NO_FAIL),
    )
    assert_dynamic_ablation(
        "3e equivalence_ablation",
        equivalence["self_ablation"],
        FAIL_EQUIVALENCE,
    )

    partition = probes["3f_partition_mislabel"]
    g_case = partition["G_external_to_derived"]
    fp_case = partition["fingerprint_27_derived_to_external"]
    print(
        "  3f G external->derived: partition_ok="
        f"{g_case['partition']['ok']}, verdict={g_case['verdict']}"
    )
    print(
        "  3f fingerprint_27 derived->external: partition_ok="
        f"{fp_case['partition']['ok']}, verdict={fp_case['verdict']}"
    )
    expect_zero(
        "3f G external->derived reaches FAIL_PROVENANCE_PARTITION",
        verdict_residual(g_case["verdict"], FAIL_PROVENANCE_PARTITION),
    )
    expect_zero(
        "3f fingerprint_27 derived->external reaches FAIL_PROVENANCE_PARTITION",
        verdict_residual(fp_case["verdict"], FAIL_PROVENANCE_PARTITION),
    )
    expect_bool(
        "3f aggregate requires BOTH mutation directions to fire",
        all(
            value == FAIL_PROVENANCE_PARTITION
            for value in partition["with_mutation_cases"].values()
        ),
    )
    assert_dynamic_ablation(
        "3f partition_ablation G external->derived",
        g_case["self_ablation"],
        FAIL_PROVENANCE_PARTITION,
    )
    assert_dynamic_ablation(
        "3f partition_ablation 27 derived->external",
        fp_case["self_ablation"],
        FAIL_PROVENANCE_PARTITION,
    )
    expect_bool(
        "3f records that G-invariance alone would miss calibrated 54/5",
        partition["g_invariance_only_would_miss_54_over_5"],
    )


def earned_expression_inventory(data: dict[str, Any]) -> dict[str, sp.Expr]:
    scaling = data["scaling"]
    equivalence = data["equivalence"]
    identity = data["identity"]
    diagnostics = data["g_invariance"]["diagnostics"]
    inventory = {
        "chi": scaling["chi"],
        "v5_slot": scaling["v5_slot"],
        "gamma_target": scaling["gamma_target"],
        "target_from_bridge": scaling["target_from_bridge"],
        "target_rhs": scaling["target_rhs"],
        "mutated_target_rhs": scaling["mutated_target_rhs"],
        "gamma_power": sp.sympify(scaling["gamma_power"]),
        "slot_power": sp.sympify(scaling["slot_power"]),
        "derived_power": sp.sympify(scaling["derived_power"]),
        "bridge_power": sp.sympify(scaling["bridge_power"]),
        "assembled_bridge_residual": scaling["assembled_bridge_residual"],
        "mutated_bridge_residual": scaling["mutated_bridge_residual"],
        "P0_physical": equivalence["P0_physical"],
        "Gamma5": equivalence["Gamma5"],
        "Gamma5_with_outgoing_chi": equivalence["Gamma5_with_outgoing_chi"],
        "forward_general": equivalence["forward_general"],
        "forward_expected": equivalence["forward_expected"],
        "forward": equivalence["forward"],
        "reverse": equivalence["reverse"],
        "wrong_gamma": equivalence["wrong_gamma"],
        "wrong_reverse": equivalence["wrong_reverse"],
        "unit": identity["unit"],
        "mag": identity["mag"],
        "twenty_seven_from_slot": identity["twenty_seven_from_slot"],
        "bound_right": identity["right"],
        "bound_residual": identity["residual"],
        "mutated_bound_right": identity["mutated_right"],
        "mutated_bound_residual": identity["mutated_residual"],
    }
    for name, item in diagnostics.items():
        inventory[f"gdiag_{name}_expression"] = item["expression"]
        inventory[f"gdiag_{name}_transformed"] = item["transformed"]
        inventory[f"gdiag_{name}_residual"] = item["residual"]
    return inventory


def run_scope_and_provenance(data: dict[str, Any]) -> None:
    subbanner("020 scope, provenance-only consumption, and PARTIAL landing")
    print(f"  020 gate booleans = {data['gates']}")
    print(f"  020 scoped verdict = {data['verdict']}")
    expect_zero(
        "020 scoped verdict lands the CALIBRATED partial component",
        verdict_residual(data["verdict"], QUAD_CALIBRATED),
    )

    inventory = earned_expression_inventory(data)
    assert_no_float("earned_expression_inventory", inventory)
    allowed_names = {"a", "c_s", "c", "G", "N0", "D0", "chi_Q", "lambda_G"}
    names_by_expression = {
        name: {symbol.name for symbol in expr.free_symbols}
        for name, expr in inventory.items()
    }
    unexpected_by_expression = {
        name: sorted(names - allowed_names)
        for name, names in names_by_expression.items()
        if not names <= allowed_names
    }
    live_names = set().union(*names_by_expression.values())
    print(f"  earned symbolic expression count under the NAME allowlist = {len(inventory)}")
    print(f"  live symbolic names across every earned expression = {sorted(live_names)}")
    expect_bool(
        "EVERY earned symbolic expression obeys the free-symbol NAME allowlist",
        not unexpected_by_expression,
    )
    expect_bool(
        "units-bearing algebra exposes exactly a,c_s,c,G,N0,D0,chi_Q,lambda_G",
        live_names == allowed_names,
    )

    verdict_signature = inspect.signature(verdict_from_partition)
    verdict_parameter_names = tuple(verdict_signature.parameters)
    verdict_code = verdict_from_partition.__code__
    verdict_dependency_names = set(verdict_parameter_names).union(
        verdict_code.co_names,
        verdict_code.co_varnames,
        verdict_code.co_freevars,
        verdict_code.co_cellvars,
    )
    verdict_dim_like_names = sorted(
        name for name in verdict_dependency_names if "dim" in name.casefold()
    )
    print(
        "  structural dimensional cut: verdict_from_partition parameters = "
        f"{verdict_parameter_names}; dim-like dependencies = {verdict_dim_like_names}"
    )
    expect_bool(
        "020-local verdict has exactly three gate+partition parameters and no dim-like dependency",
        verdict_parameter_names
        == ("scaling_ok", "equivalence_ok", "provenance_ok", "partition")
        and not verdict_dim_like_names,
    )

    print("  FOUR-WAY SPLIT / THIRD LEG: 020 carries the provenance-partition + CALIBRATED-verdict component 3/4 and lands it as PARTIAL; 018 fingerprint DONE, 019 prefactor DONE, 021 dim closure remaining.")
    print("  CONSUMED-PROVENANCE -- LABELED NON-CHECK: cites 018 chi_Q=+1 and 27 (=1/v5), 019 P0=N0/D0, and 017 l=2 port-kernel D-lanes; NO guard/dual-site.")
    print("  UNITS-BEARING-BUT-NO-DIM-GATE: c_s,a,c,G plus abstract N0,D0, chi_Q, and 27 occur in algebra/provenance only; the homogeneity closure belongs to 021.")
    print("  54/5=2*27/5-COMPUTED / PROVENANCE-PARTITION: bound symbolic residual, computed factor classes, provenance-driven verdict; the separate invariance trap is non-driving.")
    print("  STRENGTHENED-a^-5 / 3f-BOTH-DIRECTIONS: -5 flows from the frozen v5 slot; both G external->derived and 27 derived->external fire.")
    print("  EARNED-vs-CALIBRATED scope: scaling, bridge, identity, classifier and verdict rule are earned; 54/5 and G are calibrated; actual branch scaling and numerical port scalars are Gate-6 deferred.")
    print("  dropped-bookkeeping: scratch-YAML agreement handoff and report/feed writers are absent; tri-review plus the independent Wolfram route and transcript agreement replace them.")
    print("  register note: likely zero new counted knobs; registration decides the structural edge and confirms the G/c bridge accounting.")
    print("  QUAD_CALIBRATED (3/4) — the 54/5=2·27/5 provenance partition + the CALIBRATED verdict label EARNED (PARTIAL)")
    print("    = the a⁻⁵ target scaling (DERIVED via equivalence a-cancellation, not a typed target power)")
    print("    AND the Gamma5/chi_Q equivalence (forward=0, reverse=0; closes iff chi=1 and 54/27=2)")
    print("    AND the bound 54/5=2*27/5 identity with 27 derived_in_gate and 2/5 external_bridge_input")
    print("    AND the four-way provenance dominance that makes assembled 54/5 external_bridge_input -> QUAD_CALIBRATED.")
    print("  REMAINING: fingerprint=018 (DONE); prefactor algebra=019 (DONE); the mu_hat0-free [P0^phys]=1 dim closure=021 (COMPLETING leg).")
    print("  CAVEATS: actual branch a-scaling from N0/D0 and numerical (D_n,N_n) port scalars are Gate-6 deferred; 54/5 and G are calibrated (G=GENUINE_BLOCKED).")
    print("  CONSUMES (PROVENANCE only): 018 chi_Q=+1 and 27; 019 P0=N0/D0 in the Gamma5 DEF; 017 port-kernel D-lanes; no dual-site.")
    print("  EXPORTS: partition + CALIBRATED label + Gamma5/chi_Q equivalence + a^-5 scaling -> 021 + 022 + 027; no file artifact is written.")
    print("  reduction certificate: FROZEN-INPUT 018 chi_Q/27, 019 P0 definition, 017 port provenance, and GR G/c/2/5 bridge; COMPUTED scaling/equivalence/bound identity/partition/trap/verdict; DEFERRED actual branch scaling and numerical port data.")


def print_verdict_labels() -> None:
    subbanner("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): PROVENANCE_PARTITION_CALIBRATED  (the assembled quadrupole magnitude 54/5 decomposes as 54/5 = 2*27/5 [SymPy-verified rational identity]; the 27 is derived_in_gate [018's fingerprint 1/v5], the 2/5 + G are external_bridge_input [GR Burke-Thorne 2G/5c^5, G=GENUINE_BLOCKED]; a 4-way provenance partition [dominance deferred>external>derived>convention] classifies the assembled 54/5 as external_bridge_input, so the verdict lands QUAD_CALIBRATED not QUAD_PASS -- PROVENANCE-driven, NOT G->lambda*G invariance [the separate g-invariance diagnostic exposes the invariance-only trap: 54/5 is G-invariant yet calibrated])")
    print("  source top-line verdict: QUAD_CALIBRATED  (JOINT 4-stage; 020 carries the 54/5=2*27/5 provenance-partition + the CALIBRATED verdict label component 3/4 and lands the token as a PARTIAL)")
    print("  joint composition (020 = the THIRD leg; 018/019 DONE, 021 remaining): QUAD_CALIBRATED = (018: outgoing DtN Hankel fingerprint + chi_Q sign)[EARNED, PARTIAL, DONE] AND (019: squared-denominator prefactor algebra)[EARNED, PARTIAL, DONE] AND (020: 54/5=2*27/5 provenance partition + the CALIBRATED verdict label)[EARNED here, PARTIAL] AND (021: mu_hat0-free [P0^phys]=1 dim closure)[the COMPLETING leg]")
    print("  earned (the partition + bridge + scaling): the a^-5 target scaling (DERIVED via the equivalence a-cancellation, not a typed target power); the Gamma5/chi_Q equivalence 54*G*c_s^5/(5*a^5*c^5) <=> 2*G/(5*c^5) (closes iff chi=1 and 54/27=2); the 54/5=2*27/5 SymPy identity; the 4-way provenance partition (assembled 54/5 = external_bridge_input, external dominates); the CALIBRATED verdict (provenance-driven, NOT G-invariance)")
    print("  earned (able-to-fail): 3c_wrong_scaling (STRENGTHENED -- a^5->a^4 breaks the equivalence a-cancellation -> FAIL_SCALING); 3e_equivalence_break (2/5->3/5 -> FAIL_EQUIVALENCE); 3f_partition_mislabel (BOTH directions: G external->derived AND 27 derived->external -> FAIL_PROVENANCE_PARTITION), each with a DYNAMIC 020-local self-ablation (re-run, NOT a constant, NOT the joint base_verdict)")
    print("  calibrated / deferred (NOT 020): the fingerprint (018, DONE); the prefactor algebra (019, DONE); the mu_hat0-free dim closure (021); the 54/5 magnitude + G (external_bridge_input, G=GENUINE_BLOCKED); the ACTUAL branch a-scaling + the numerical (D_n,N_n) port scalars (Gate-6 sim-deferred, report :49)")
    print("  consumed (PROVENANCE only -- NO guard/dual-site): 018's chi_Q=+1 + the 27 (=1/v5) [enter 020's self-contained equivalence bridge] + 019's P0=N0/D0 [Gamma5 def] + 017's l=2 port kernel D-lanes; NO mu_hat0/dim gate")
    print("  exports: the 54/5=2*27/5 partition + the CALIBRATED label + the Gamma5/chi_Q equivalence + the a^-5 scaling => 021 (dim closure completes the joint) + 022 (non-regression) + 027 (pathA_43 closure)")
    print("  new symbols first-appearing (020): none new-counted (G=GENUINE_BLOCKED/external_bridge_input, already in the register; c=GR/PN light-speed bridge = c_gamma role, cited benchmark; the 2/5=GR bridge; the 27=018's derived fingerprint); units-bearing (c_s/a/c/G live) but NO dim gate; no counted knobs expected (an EARNED classification slice, like 018/019)")


def main() -> int:
    banner("ledger_stage020_provenance_partition_sympy_audit")
    print("Target stem confirmed: ledger_stage020_provenance_partition")
    print("Engine: SymPy exact symbolic/rational algebra and tag-dominance classification; no floats/tolerances; zero file I/O.")
    data = build_baseline()
    run_scaling_and_equivalence(data)
    run_identity_and_classifier(data)
    run_g_invariance_diagnostic(data)
    run_verdict_controls(data)
    run_probes(data)
    run_scope_and_provenance(data)
    print_verdict_labels()
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
    except Exception as exc:
        if not isinstance(exc, AuditFailure):
            _record_fail(f"UNCAUGHT exception: {exc!r}")
        banner("Tallies")
        total = PASS_COUNT + FAIL_COUNT
        print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
        print("OVERALL FAIL")
        raise SystemExit(1) from exc
    banner("Tallies")
    total = PASS_COUNT + FAIL_COUNT
    print(f"TALLY sympy: {PASS_COUNT} pass + {FAIL_COUNT} fail = {total} checks")
    if FAIL_COUNT == 0 and exit_code == 0:
        print("OVERALL PASS")
        raise SystemExit(0)
    print("OVERALL FAIL")
    raise SystemExit(1)
