#!/usr/bin/env python3
from __future__ import annotations

import itertools
import math

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


def expect_equal(name: str, lhs, rhs) -> None:
    expect_zero(name, lhs - rhs)


def expect_true(name: str, condition: bool) -> None:
    print(f"{name} = {condition}")
    if not condition:
        raise AssertionError(f"{name} failed")


def packet_interval(packet_windows: dict[str, tuple[int, ...]]) -> tuple[int, int]:
    return (
        min(min(values) for values in packet_windows.values()),
        min(max(values) for values in packet_windows.values()),
    )


def boundary_best(assignment: dict[str, int], boundary_packet_names: tuple[str, ...]) -> int:
    return min(assignment[name] for name in boundary_packet_names)


def full_best(assignment: dict[str, int], boundary_packet_names: tuple[str, ...]) -> int:
    return min(boundary_best(assignment, boundary_packet_names), assignment["support5_int"])


def classify_family(
    family_name: str,
    packet_windows: dict[str, tuple[int, ...]],
    boundary_packet_names: tuple[str, ...],
) -> dict[str, int]:
    keys = tuple(packet_windows)
    full_lo, full_hi = packet_interval(packet_windows)
    counts = {"support5": 0, "boundary": 0, "tie": 0}

    for values in itertools.product(*(packet_windows[key] for key in keys)):
        assignment = dict(zip(keys, values, strict=True))
        boundary_value = boundary_best(assignment, boundary_packet_names)
        support5_value = assignment["support5_int"]
        full_value = min(boundary_value, support5_value)
        if not (full_lo <= full_value <= full_hi):
            raise AssertionError(f"{family_name} interval splice failed")
        if support5_value < boundary_value:
            counts["support5"] += 1
        elif boundary_value < support5_value:
            counts["boundary"] += 1
        else:
            counts["tie"] += 1

    print(f"{family_name} exhaustive outcomes = {counts}")
    return counts


banner("STAGE 201 — FULL SUPPORT-<=5 COMPLETION AND LOCAL MIXED-RAY SEARCH CLOSURE")

subbanner("I. Exact boundary-identification via imported Stage 215 face packets")

axes = ("lambda", "c", "gamma", "U", "W")
quadruple_faces = {
    f"omit_{axis}": {
        "omitted_axis": axis,
        "support": tuple(candidate for candidate in axes if candidate != axis),
        "source_stage": 198,
    }
    for axis in axes
}
boundary_packets = {
    "support_le3": {
        "support": "<=3 global ledger",
        "source_stage": 198,
    },
    **quadruple_faces,
}
boundary_face_names = tuple(quadruple_faces)
proper_faces = []
for k in range(1, 5):
    proper_faces.extend(itertools.combinations(axes, k))

quadruple_supports = {
    tuple(sorted(packet["support"])) for packet in quadruple_faces.values()
}
expected_quadruples = {
    tuple(sorted(face)) for face in itertools.combinations(axes, 4)
}

print("primitive axes =", axes)
print("imported Stage 215 boundary packets =", tuple(boundary_packets))
print("quadruple face supports =", tuple(packet["support"] for packet in quadruple_faces.values()))
print("#proper nonempty support strata =", len(proper_faces))

expect_zero("primitive axes - 5", len(axes) - 5)
expect_zero("imported quadruple packet count - 5", len(quadruple_faces) - 5)
expect_true(
    "imported quadruple supports match the five simplex facets",
    quadruple_supports == expected_quadruples,
)
expect_zero("support-cardinality ceiling 5 - #axes", len(axes) - 5)

for subset in proper_faces:
    covering_faces = [
        name
        for name, packet in quadruple_faces.items()
        if set(subset).issubset(packet["support"])
    ]
    expected = len(axes) - len(subset)
    print(f"boundary coverage incidence {subset} -> {covering_faces}")
    if len(covering_faces) != expected:
        raise AssertionError("boundary-identification coverage count is incorrect")

subbanner("II. Exact imported support-<=4 and support-<=5 ledger splice")

tau_le3_best = sp.symbols("tau_le3_best", real=True)
tau_le3_lo, tau_le3_hi = sp.symbols("tau_le3_lo tau_le3_hi", real=True)
tau5_best_int = sp.symbols("tau_5_best_int", real=True)
tau5_lo_int, tau5_hi_int = sp.symbols("tau_5_lo_int tau_5_hi_int", real=True)

boundary_best_symbols = [tau_le3_best]
boundary_lo_symbols = [tau_le3_lo]
boundary_hi_symbols = [tau_le3_hi]

for face_name in boundary_face_names:
    symbol_stub = face_name.replace("omit_", "face_")
    best_sym = sp.symbols(f"tau_{symbol_stub}_best", real=True)
    lo_sym = sp.symbols(f"tau_{symbol_stub}_lo", real=True)
    hi_sym = sp.symbols(f"tau_{symbol_stub}_hi", real=True)
    quadruple_faces[face_name]["best_symbol"] = best_sym
    quadruple_faces[face_name]["lo_symbol"] = lo_sym
    quadruple_faces[face_name]["hi_symbol"] = hi_sym
    boundary_best_symbols.append(best_sym)
    boundary_lo_symbols.append(lo_sym)
    boundary_hi_symbols.append(hi_sym)

tau_le4_best = sp.Min(*boundary_best_symbols)
tau_le4_lo = sp.Min(*boundary_lo_symbols)
tau_le4_hi = sp.Min(*boundary_hi_symbols)

tau_le5_best = sp.Min(tau_le4_best, tau5_best_int)
tau_le5_lo = sp.Min(tau_le4_lo, tau5_lo_int)
tau_le5_hi = sp.Min(tau_le4_hi, tau5_hi_int)

tau_le5_best_flat = sp.Min(*boundary_best_symbols, tau5_best_int)
tau_le5_lo_flat = sp.Min(*boundary_lo_symbols, tau5_lo_int)
tau_le5_hi_flat = sp.Min(*boundary_hi_symbols, tau5_hi_int)

print("tau_{<=4,*}^{best} =")
sp.pprint(tau_le4_best)
print("tau_{<=5,*}^{best} =")
sp.pprint(tau_le5_best)
print("tau_{<=5,lo} =")
sp.pprint(tau_le5_lo)
print("tau_{<=5,hi} =")
sp.pprint(tau_le5_hi)

expect_equal("support<=5 best flattening over imported packets", tau_le5_best, tau_le5_best_flat)
expect_equal("support<=5 lower splice over imported packets", tau_le5_lo, tau_le5_lo_flat)
expect_equal("support<=5 upper splice over imported packets", tau_le5_hi, tau_le5_hi_flat)

subbanner("III. Exact improvement / no-improvement / overlap families on the actual finite ledger")

improvement_family = {
    "support_le3": (10, 11),
    "omit_lambda": (12, 13),
    "omit_c": (11, 12),
    "omit_gamma": (13, 14),
    "omit_U": (15, 16),
    "omit_W": (14, 15),
    "support5_int": (2, 3, 4),
}
no_improvement_family = {
    "support_le3": (2, 3),
    "omit_lambda": (4, 5),
    "omit_c": (3, 4),
    "omit_gamma": (5, 6),
    "omit_U": (6, 7),
    "omit_W": (4, 6),
    "support5_int": (9, 10, 11),
}
overlap_family = {
    "support_le3": (5, 6),
    "omit_lambda": (4, 8),
    "omit_c": (7, 8),
    "omit_gamma": (6, 9),
    "omit_U": (8, 9),
    "omit_W": (5, 7),
    "support5_int": (3, 7),
}

improvement_counts = classify_family(
    "genuine support-5 improvement family", improvement_family, boundary_face_names + ("support_le3",)
)
expect_zero("support-5 improvement family boundary wins", improvement_counts["boundary"])
expect_zero("support-5 improvement family ties", improvement_counts["tie"])
expect_true(
    "support-5 improvement family interior wins exist",
    improvement_counts["support5"] > 0,
)

no_improvement_counts = classify_family(
    "support-5 no-improvement family", no_improvement_family, boundary_face_names + ("support_le3",)
)
expect_zero("support-5 no-improvement family interior wins", no_improvement_counts["support5"])
expect_zero("support-5 no-improvement family ties", no_improvement_counts["tie"])
expect_true(
    "support-5 no-improvement family boundary wins exist",
    no_improvement_counts["boundary"] > 0,
)

overlap_counts = classify_family(
    "support-5 overlap family", overlap_family, boundary_face_names + ("support_le3",)
)
expect_true(
    "support-5 overlap family retains boundary winners",
    overlap_counts["boundary"] > 0,
)
expect_true(
    "support-5 overlap family retains interior winners",
    overlap_counts["support5"] > 0,
)

subbanner("IV. Exact support-five candidate filters and final budget ledger")

support_five_packet = {
    "source_stage": 200,
    "canonical_screens": ("gradient-optimal", "equal-mix"),
    "preferred_lifted_degree_pattern": (3, 3, 3, 3, 2),
    "fallback_projected_degree_pattern": (5, 5, 5, 6),
}
lifted_per_envelope = math.prod(support_five_packet["preferred_lifted_degree_pattern"])
projected_per_envelope = math.prod(support_five_packet["fallback_projected_degree_pattern"])

print("support-five canonical screens =", support_five_packet["canonical_screens"])
print("preferred lifted degree pattern =", support_five_packet["preferred_lifted_degree_pattern"])
print("fallback projected degree pattern =", support_five_packet["fallback_projected_degree_pattern"])

expect_zero("support-five canonical screen count - 2", len(support_five_packet["canonical_screens"]) - 2)
expect_zero("lifted compiler bound - 162", lifted_per_envelope - 162)
expect_zero("projected compiler bound - 750", projected_per_envelope - 750)

# Stage 215 carries the support-<=3 ledger budget 600 and the exact five
# primitive quadruple face packets, each with 54 interior evaluations per envelope.
support_le3_budget = 600
quadruple_eval_per_envelope = 54
support_le4_budget = support_le3_budget + len(quadruple_faces) * 2 * quadruple_eval_per_envelope

preferred_total = support_le4_budget + 2 * lifted_per_envelope
fallback_total = support_le4_budget + 2 * projected_per_envelope

print(f"support<=3 imported budget = {support_le3_budget}")
print(f"support<=4 rebuilt budget = {support_le4_budget}")
print(f"preferred full support<=5 budget = {preferred_total}")
print(f"fallback full support<=5 budget = {fallback_total}")

expect_zero("support<=4 rebuilt budget - 1140", support_le4_budget - 1140)
expect_zero("preferred full support<=5 budget - 1464", preferred_total - 1464)
expect_zero("fallback full support<=5 budget - 2640", fallback_total - 2640)

print("\nAll Stage 218 imported ledger, splice, classification, and budget checks verified.")
