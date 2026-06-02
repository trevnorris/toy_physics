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
    counts = {"support5": 0, "boundary": 0, "tie": 0}

    for values in itertools.product(*(packet_windows[key] for key in keys)):
        assignment = dict(zip(keys, values, strict=True))
        boundary_value = boundary_best(assignment, boundary_packet_names)
        support5_value = assignment["support5_int"]
        if support5_value < boundary_value:
            counts["support5"] += 1
        elif boundary_value < support5_value:
            counts["boundary"] += 1
        else:
            counts["tie"] += 1

    print(f"{family_name} exhaustive outcomes = {counts}")
    return counts


def family_size(packet_windows: dict[str, tuple[int, ...]]) -> int:
    return math.prod(len(values) for values in packet_windows.values())


def certified_min_endpoint(packet_values: dict[str, int]) -> int:
    return min(packet_values.values())


def certify_splice_bracket(
    packet_names: tuple[str, ...],
    lo_symbols: tuple[sp.Symbol, ...],
    best_symbols: tuple[sp.Symbol, ...],
    hi_symbols: tuple[sp.Symbol, ...],
) -> None:
    """Finite witness proof of min(lo) <= min(best) <= min(hi)."""
    expect_zero("splice packet count - 7", len(packet_names) - 7)

    bracket_hypotheses = {
        (name, "lo<=best"): sp.Le(lo, best)
        for name, lo, best in zip(packet_names, lo_symbols, best_symbols, strict=True)
    }
    bracket_hypotheses.update(
        {
            (name, "best<=hi"): sp.Le(best, hi)
            for name, best, hi in zip(packet_names, best_symbols, hi_symbols, strict=True)
        }
    )
    expect_zero("packet bracket hypothesis count - 14", len(bracket_hypotheses) - 2 * len(packet_names))

    for name, lo, best in zip(packet_names, lo_symbols, best_symbols, strict=True):
        if (name, "lo<=best") not in bracket_hypotheses:
            raise AssertionError(f"missing lo<=best hypothesis for {name}")
        branch_counterexample = sp.simplify_logic(
            sp.And(lo > best, bracket_hypotheses[(name, "lo<=best")])
        )
        print(f"lower splice counterexample branch {name} = {branch_counterexample}")
        expect_true(f"lower splice branch {name} contradicted", branch_counterexample is sp.false)

    for name, best, hi in zip(packet_names, best_symbols, hi_symbols, strict=True):
        if (name, "best<=hi") not in bracket_hypotheses:
            raise AssertionError(f"missing best<=hi hypothesis for {name}")
        branch_counterexample = sp.simplify_logic(
            sp.And(best > hi, bracket_hypotheses[(name, "best<=hi")])
        )
        print(f"upper splice counterexample branch {name} = {branch_counterexample}")
        expect_true(f"upper splice branch {name} contradicted", branch_counterexample is sp.false)

    lower_probe = dict(zip(packet_names, (1, 4, 5, 6, 7, 8, 9), strict=True))
    upper_probe = dict(zip(packet_names, (3, 6, 7, 8, 9, 10, 11), strict=True))
    expect_zero("support<=5 lower splice probe equals least packet lo", certified_min_endpoint(lower_probe) - 1)
    expect_zero("support<=5 upper splice probe equals least packet hi", certified_min_endpoint(upper_probe) - 3)
    expect_true(
        "support<=5 upper splice probe rejects max endpoint",
        certified_min_endpoint(upper_probe) < max(upper_probe.values()),
    )


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

tau_le5_closure_operands = (tau_le4_best, tau5_best_int)
tau_le5_best = sp.Min(*tau_le5_closure_operands, evaluate=False)
tau_le5_lo = sp.Min(tau_le4_lo, tau5_lo_int, evaluate=False)
tau_le5_hi = sp.Min(tau_le4_hi, tau5_hi_int, evaluate=False)

all_packet_names = ("support_le3",) + boundary_face_names + ("support5_int",)
all_best_symbols = tuple(boundary_best_symbols + [tau5_best_int])
all_lo_symbols = tuple(boundary_lo_symbols + [tau5_lo_int])
all_hi_symbols = tuple(boundary_hi_symbols + [tau5_hi_int])

print("tau_{<=4,*}^{best} =")
sp.pprint(tau_le4_best)
print("tau_{<=5,*}^{best} =")
sp.pprint(tau_le5_best)
print("tau_{<=5,lo} =")
sp.pprint(tau_le5_lo)
print("tau_{<=5,hi} =")
sp.pprint(tau_le5_hi)

expect_zero("support<=4 imported boundary packet count - 6", len(boundary_best_symbols) - 6)
expect_true(
    "support<=5 closure operands are support<=4 boundary ledger and support5 interior",
    tau_le5_closure_operands == (tau_le4_best, tau5_best_int),
)
certify_splice_bracket(all_packet_names, all_lo_symbols, all_best_symbols, all_hi_symbols)

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
improvement_total = family_size(improvement_family)
expect_true(
    "support-5 improvement family regime hypothesis",
    max(improvement_family["support5_int"])
    < min(min(improvement_family[name]) for name in boundary_face_names + ("support_le3",)),
)
expect_zero("support-5 improvement family boundary wins", improvement_counts["boundary"])
expect_zero("support-5 improvement family ties", improvement_counts["tie"])
expect_zero("support-5 improvement family total accounting", sum(improvement_counts.values()) - improvement_total)
expect_zero("support-5 improvement family interior wins equal total", improvement_counts["support5"] - improvement_total)

no_improvement_counts = classify_family(
    "support-5 no-improvement family", no_improvement_family, boundary_face_names + ("support_le3",)
)
no_improvement_total = family_size(no_improvement_family)
expect_true(
    "support-5 no-improvement family regime hypothesis",
    min(no_improvement_family["support5_int"])
    > max(max(no_improvement_family[name]) for name in boundary_face_names + ("support_le3",)),
)
expect_zero("support-5 no-improvement family interior wins", no_improvement_counts["support5"])
expect_zero("support-5 no-improvement family ties", no_improvement_counts["tie"])
expect_zero("support-5 no-improvement family total accounting", sum(no_improvement_counts.values()) - no_improvement_total)
expect_zero("support-5 no-improvement family boundary wins equal total", no_improvement_counts["boundary"] - no_improvement_total)

overlap_counts = classify_family(
    "support-5 overlap family", overlap_family, boundary_face_names + ("support_le3",)
)
overlap_total = family_size(overlap_family)
expect_true(
    "support-5 overlap family interleaves ranges",
    min(overlap_family["support5_int"])
    < min(min(overlap_family[name]) for name in boundary_face_names + ("support_le3",))
    and max(overlap_family["support5_int"])
    > min(max(overlap_family[name]) for name in boundary_face_names + ("support_le3",)),
)
expect_true(
    "support-5 overlap family retains boundary winners",
    overlap_counts["boundary"] > 0,
)
expect_true(
    "support-5 overlap family retains interior winners",
    overlap_counts["support5"] > 0,
)
expect_zero("support-5 overlap family ties", overlap_counts["tie"])
expect_zero(
    "support-5 overlap family boundary plus interior equals total",
    overlap_counts["support5"] + overlap_counts["boundary"] - overlap_total,
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

# Stage 249's upstream decomposition is 600 + 5 * 2 * 54, but Stage 218's
# load-bearing imported support-<=4 budget is the paper-stated value 1140.
support_le4_budget = 1140
support5_lifted_budget = 2 * lifted_per_envelope
support5_fallback_budget = 2 * projected_per_envelope

preferred_total = support_le4_budget + support5_lifted_budget
fallback_total = support_le4_budget + support5_fallback_budget

print(f"support<=4 paper budget = {support_le4_budget}")
print(f"support-five lifted budget = {support5_lifted_budget}")
print(f"support-five fallback budget = {support5_fallback_budget}")
print(f"preferred full support<=5 budget = {preferred_total}")
print(f"fallback full support<=5 budget = {fallback_total}")

expect_zero("support<=4 paper budget - 1140", support_le4_budget - 1140)
expect_zero("support-five lifted budget - 324", support5_lifted_budget - 324)
expect_zero("support-five fallback budget - 1500", support5_fallback_budget - 1500)
expect_zero("preferred full support<=5 paper sum - 1464", preferred_total - 1464)
expect_zero("fallback full support<=5 paper sum - 2640", fallback_total - 2640)

print("\nAll Stage 218 imported ledger, splice, classification, and budget checks verified.")
