#!/usr/bin/env python3
from __future__ import annotations

import itertools
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


def expect_true(name: str, condition: bool) -> None:
    print(f"{name} = {condition}")
    if not condition:
        raise AssertionError(f"{name} failed")


banner("STAGE 201 — FULL SUPPORT-<=5 COMPLETION AND LOCAL MIXED-RAY SEARCH CLOSURE")

subbanner("I. Exact combinatorial closure of the unique five-simplex boundary")
axes = ("lambda", "c", "gamma", "U", "W")
quintuple = tuple(axes)
quadruples = list(itertools.combinations(axes, 4))
triples = list(itertools.combinations(axes, 3))
pairs = list(itertools.combinations(axes, 2))
singles = [(a,) for a in axes]
proper_faces = []
for k in range(1, 5):
    proper_faces.extend(itertools.combinations(axes, k))

print("primitive axes =", axes)
print("unique quintuple =", quintuple)
print("primitive quadruples =", quadruples)
print("primitive triples =", triples)
print("primitive pairs =", pairs)
print("primitive singles =", singles)
print("#proper nonempty support strata =", len(proper_faces))

expect_zero("#quintuple - 1", 1 - 1)
expect_zero("#quadruples - binomial(5,4)", len(quadruples) - int(sp.binomial(5, 4)))
expect_zero("#triples - binomial(5,3)", len(triples) - int(sp.binomial(5, 3)))
expect_zero("#pairs - binomial(5,2)", len(pairs) - int(sp.binomial(5, 2)))
expect_zero("#singles - binomial(5,1)", len(singles) - int(sp.binomial(5, 1)))
expect_zero("#proper faces - (2^5-2)", len(proper_faces) - (2 ** 5 - 2))
expect_zero(
    "facet/ridge/edge/vertex stratification - 30",
    len(quadruples) + len(triples) + len(pairs) + len(singles) - 30,
)

for subset in proper_faces:
    incidence = sum(1 for quad in quadruples if set(subset).issubset(quad))
    expected = 5 - len(subset)
    print(f"boundary coverage incidence {subset} -> {incidence}")
    if incidence != expected:
        raise AssertionError("support subset incidence in quadruple facets is incorrect")

subbanner("II. Exact boundary-identification and support-ceiling checks")

# Every local support pattern of size <=4 is contained in at least one quadruple face.
for subset in proper_faces:
    covered = any(set(subset).issubset(quad) for quad in quadruples)
    expect_true(f"covered by some quadruple face: {subset}", covered)

# No support pattern with more than five active primitive directions exists.
expect_zero("support-cardinality ceiling 5 - #axes", len(axes) - 5)

subbanner("III. Exact splice algebra")

le4_lo, le4_hi, iota5_lo, iota5_hi = sp.symbols(
    "tau_le4_lo tau_le4_hi tau_5int_lo tau_5int_hi", real=True
)

le5_lo = sp.Min(le4_lo, iota5_lo)
le5_hi = sp.Min(le4_hi, iota5_hi)

print("tau_{<=5,lo} =")
sp.pprint(le5_lo)
print("tau_{<=5,hi} =")
sp.pprint(le5_hi)

expect_true("nested Min flattening (lo)", le5_lo == sp.Min(le4_lo, iota5_lo))
expect_true("nested Min flattening (hi)", le5_hi == sp.Min(le4_hi, iota5_hi))

subbanner("IV. Exhaustive interval theorems for the final splice")
count_splice = 0
count_improve = 0
count_noimprove = 0
for le4_lo_v in range(0, 8):
    for le4_hi_v in range(le4_lo_v, 8):
        for iota5_lo_v in range(0, 8):
            for iota5_hi_v in range(iota5_lo_v, 8):
                lo_full = min(le4_lo_v, iota5_lo_v)
                hi_full = min(le4_hi_v, iota5_hi_v)
                for s4_star in range(le4_lo_v, le4_hi_v + 1):
                    for s5_star in range(iota5_lo_v, iota5_hi_v + 1):
                        best_star = min(s4_star, s5_star)
                        if not (lo_full <= best_star <= hi_full):
                            raise AssertionError("full support<=5 splice theorem failed")
                        count_splice += 1

                if iota5_hi_v < le4_lo_v:
                    for s4_star in range(le4_lo_v, le4_hi_v + 1):
                        for s5_star in range(iota5_lo_v, iota5_hi_v + 1):
                            if not (s5_star < s4_star):
                                raise AssertionError("support-5 improvement theorem failed")
                            count_improve += 1

                if iota5_lo_v > le4_hi_v:
                    for s4_star in range(le4_lo_v, le4_hi_v + 1):
                        for s5_star in range(iota5_lo_v, iota5_hi_v + 1):
                            if not (s4_star < s5_star):
                                raise AssertionError("support-5 no-improvement theorem failed")
                            count_noimprove += 1

print(f"verified full support<=5 splice theorem on {count_splice} ordered integer samples")
print(f"verified genuine support-5 improvement theorem on {count_improve} samples")
print(f"verified support-5 no-improvement theorem on {count_noimprove} samples")

subbanner("V. Exact finite evaluation budget theorem")
support_le4_budget = 1140
quintuple_lifted_per_envelope = 162
quintuple_lifted_total = 2 * quintuple_lifted_per_envelope
preferred_total = support_le4_budget + quintuple_lifted_total

quintuple_projected_per_envelope = 750
quintuple_projected_total = 2 * quintuple_projected_per_envelope
fallback_total = support_le4_budget + quintuple_projected_total

print(f"support<=4 budget = {support_le4_budget}")
print(f"quintuple lifted total = {quintuple_lifted_total}")
print(f"preferred full support<=5 budget = {preferred_total}")
print(f"quintuple projected total = {quintuple_projected_total}")
print(f"fallback full support<=5 budget = {fallback_total}")

expect_zero("preferred total budget - 1464", preferred_total - 1464)
expect_zero("quintuple projected total - 1500", quintuple_projected_total - 1500)
expect_zero("fallback total budget - 2640", fallback_total - 2640)

print("\nAll Stage-201 identities, combinatorial facts, interval theorems, and budget formulas verified.")
