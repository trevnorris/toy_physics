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


banner("STAGE 215 — FULL PRIMITIVE-QUADRUPLE RANKING AND THE UP-TO-FOUR-COORDINATE SIEVE")

subbanner("I. Exact combinatorial ledger")
axes = ("lambda", "c", "gamma", "U", "W")
triples = list(itertools.combinations(axes, 3))
quadruples = list(itertools.combinations(axes, 4))

print("primitive axes =", axes)
print("primitive triples =", triples)
print("primitive quadruples =", quadruples)
print("#triples =", len(triples))
print("#quadruples =", len(quadruples))

expect_zero("#triples - binomial(5,3)", len(triples) - int(sp.binomial(5, 3)))
expect_zero("#quadruples - binomial(5,4)", len(quadruples) - int(sp.binomial(5, 4)))

for quad in quadruples:
    face_count = sum(1 for tri in triples if set(tri).issubset(quad))
    print(f"quadruple face incidence {quad} -> {face_count}")
    if face_count != 4:
        raise AssertionError("each primitive quadruple must have exactly four triple faces")

for tri in triples:
    incidence = sum(1 for quad in quadruples if set(tri).issubset(quad))
    print(f"triple incidence {tri} -> {incidence}")
    if incidence != 2:
        raise AssertionError("each primitive triple must belong to exactly two quadruples")

for axis in axes:
    incidence = sum(1 for quad in quadruples if axis in quad)
    print(f"axis incidence {axis} -> {incidence}")
    if incidence != 4:
        raise AssertionError("each primitive axis must belong to exactly four quadruples")

subbanner("II. Exact boundary-splice algebra")

tijk_lo, tijl_lo, tikl_lo, tjkl_lo, iota_lo = sp.symbols(
    "tau_ijk_lo tau_ijl_lo tau_ikl_lo tau_jkl_lo iota_lo", real=True
)
tijk_hi, tijl_hi, tikl_hi, tjkl_hi, iota_hi = sp.symbols(
    "tau_ijk_hi tau_ijl_hi tau_ikl_hi tau_jkl_hi iota_hi", real=True
)

beta_lo = sp.Min(tijk_lo, tijl_lo, tikl_lo, tjkl_lo)
beta_hi = sp.Min(tijk_hi, tijl_hi, tikl_hi, tjkl_hi)
full_lo = sp.Min(iota_lo, beta_lo)
full_hi = sp.Min(iota_hi, beta_hi)

print("beta_lo =")
sp.pprint(beta_lo)
print("beta_hi =")
sp.pprint(beta_hi)
print("full_lo =")
sp.pprint(full_lo)
print("full_hi =")
sp.pprint(full_hi)

expect_true("nested Min flattening (lo)", full_lo == sp.Min(iota_lo, tijk_lo, tijl_lo, tikl_lo, tjkl_lo))
expect_true("nested Min flattening (hi)", full_hi == sp.Min(iota_hi, tijk_hi, tijl_hi, tikl_hi, tjkl_hi))

subbanner("III. Exhaustive interval-theorem checks")
count_local = 0
for beta_lo_v in range(0, 6):
    for beta_hi_v in range(beta_lo_v, 6):
        for iota_lo_v in range(0, 6):
            for iota_hi_v in range(iota_lo_v, 6):
                lo_full = min(beta_lo_v, iota_lo_v)
                hi_full = min(beta_hi_v, iota_hi_v)
                for b_star in range(beta_lo_v, beta_hi_v + 1):
                    for i_star in range(iota_lo_v, iota_hi_v + 1):
                        simplex_best = min(b_star, i_star)
                        if not (lo_full <= simplex_best <= hi_full):
                            raise AssertionError("local quadruple full-simplex interval theorem failed")
                        count_local += 1
print(f"verified local full-simplex interval theorem on {count_local} ordered integer samples")

count_int = 0
count_bdry = 0
for beta_lo_v in range(1, 8):
    for beta_hi_v in range(beta_lo_v, 8):
        for iota_lo_v in range(0, 8):
            for iota_hi_v in range(iota_lo_v, 8):
                if iota_hi_v < beta_lo_v:
                    for b_star in range(beta_lo_v, beta_hi_v + 1):
                        for i_star in range(iota_lo_v, iota_hi_v + 1):
                            if not (i_star < b_star):
                                raise AssertionError("interior-certified quadruple ordering failed")
                            count_int += 1
                if iota_lo_v > beta_hi_v:
                    for b_star in range(beta_lo_v, beta_hi_v + 1):
                        for i_star in range(iota_lo_v, iota_hi_v + 1):
                            if not (b_star < i_star):
                                raise AssertionError("boundary-certified quadruple ordering failed")
                            count_bdry += 1
print(f"verified interior-certified local order on {count_int} samples")
print(f"verified boundary-certified local order on {count_bdry} samples")

count_rank = 0
for L1 in range(0, 6):
    for U1 in range(L1, 6):
        for L2 in range(0, 8):
            for U2 in range(L2, 8):
                if U1 < L2:
                    for x in range(L1, U1 + 1):
                        for y in range(L2, U2 + 1):
                            if not (x < y):
                                raise AssertionError("primitive-quadruple ranking theorem failed")
                            count_rank += 1
print(f"verified certified interval ordering on {count_rank} samples")

count_unique = 0
interval_samples = [(L, U) for L in range(0, 5) for U in range(L, 5)]
for intervals in itertools.product(interval_samples, repeat=5):
    for star in range(5):
        other_lows = [intervals[p][0] for p in range(5) if p != star]
        if intervals[star][1] < min(other_lows):
            for x_star in range(intervals[star][0], intervals[star][1] + 1):
                for p in range(5):
                    if p == star:
                        continue
                    for x_p in range(intervals[p][0], intervals[p][1] + 1):
                        if not (x_star < x_p):
                            raise AssertionError("unique primitive-quadruple certified winner theorem failed")
                        count_unique += 1
print(f"verified unique certified winner theorem on {count_unique} min-over-others samples")

count_global = 0
for le3_lo in range(0, 6):
    for le3_hi in range(le3_lo, 6):
        for quad_lo in range(0, 6):
            for quad_hi in range(quad_lo, 6):
                lo_global = min(le3_lo, quad_lo)
                hi_global = min(le3_hi, quad_hi)
                for s3_star in range(le3_lo, le3_hi + 1):
                    for q4_star in range(quad_lo, quad_hi + 1):
                        best_star = min(s3_star, q4_star)
                        if not (lo_global <= best_star <= hi_global):
                            raise AssertionError("up-to-four-coordinate splice theorem failed")
                        count_global += 1
print(f"verified global splice theorem on {count_global} samples")

count_improve = 0
count_noimprove = 0
for le3_lo in range(1, 8):
    for le3_hi in range(le3_lo, 8):
        for quad_lo in range(0, 8):
            for quad_hi in range(quad_lo, 8):
                if quad_hi < le3_lo:
                    for s3_star in range(le3_lo, le3_hi + 1):
                        for q4_star in range(quad_lo, quad_hi + 1):
                            if not (q4_star < s3_star):
                                raise AssertionError("global four-coordinate improvement theorem failed")
                            count_improve += 1
                if quad_lo > le3_hi:
                    for s3_star in range(le3_lo, le3_hi + 1):
                        for q4_star in range(quad_lo, quad_hi + 1):
                            if not (s3_star < q4_star):
                                raise AssertionError("global four-coordinate no-improvement theorem failed")
                            count_noimprove += 1
print(f"verified global four-coordinate improvement theorem on {count_improve} samples")
print(f"verified global four-coordinate no-improvement theorem on {count_noimprove} samples")

subbanner("IV. Exact finite evaluation budget theorem")
support_le3_budget = 10 * 12 + 10 * 48
quad_eval_per_envelope = (3 * 3 * 3) * 2
quad_eval_total = len(quadruples) * 2 * quad_eval_per_envelope
full_eval_total = support_le3_budget + quad_eval_total

print(f"imported support<=3 budget = {support_le3_budget}")
print(f"interior quadruple evaluation budget = {quad_eval_total}")
print(f"full support<=4 budget = {full_eval_total}")

expect_zero("quadruple interior budget - 540", quad_eval_total - 540)
expect_zero("full support<=4 budget - 1140", full_eval_total - 1140)

print("\nAll Stage 215 identities and interval theorems verified.")
