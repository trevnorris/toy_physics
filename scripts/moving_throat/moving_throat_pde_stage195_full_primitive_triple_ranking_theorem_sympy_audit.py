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


banner("STAGE 195 — FULL PRIMITIVE-TRIPLE RANKING AND THE UP-TO-THREE-COORDINATE SIEVE")

# ---------------------------------------------------------------------------
# I. Exact combinatorial ledger for the five primitive free directions
# ---------------------------------------------------------------------------
subbanner("I. Exact combinatorial ledger")

axes = ("lambda", "c", "gamma", "U", "W")
pairs = list(itertools.combinations(axes, 2))
triples = list(itertools.combinations(axes, 3))

print("primitive axes =", axes)
print("primitive pairs =", pairs)
print("primitive triples =", triples)
print("#pairs =", len(pairs))
print("#triples =", len(triples))

expect_zero("#pairs - binomial(5,2)", len(pairs) - int(sp.binomial(5, 2)))
expect_zero("#triples - binomial(5,3)", len(triples) - int(sp.binomial(5, 3)))

for pair in pairs:
    incidence = sum(1 for tri in triples if set(pair).issubset(tri))
    print(f"pair incidence {pair} -> {incidence}")
    if incidence != 3:
        raise AssertionError("each primitive pair must belong to exactly three triples")

for axis in axes:
    incidence = sum(1 for tri in triples if axis in tri)
    print(f"axis incidence {axis} -> {incidence}")
    if incidence != 6:
        raise AssertionError("each primitive axis must belong to exactly six triples")

# ---------------------------------------------------------------------------
# II. Exact boundary-splice algebra in Min form
# ---------------------------------------------------------------------------
subbanner("II. Exact boundary-splice algebra")

tij_lo, tik_lo, tjk_lo, iota_lo = sp.symbols("tau_ij_lo tau_ik_lo tau_jk_lo iota_lo", real=True)
tij_hi, tik_hi, tjk_hi, iota_hi = sp.symbols("tau_ij_hi tau_ik_hi tau_jk_hi iota_hi", real=True)

beta_lo = sp.Min(tij_lo, tik_lo, tjk_lo)
beta_hi = sp.Min(tij_hi, tik_hi, tjk_hi)
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

expect_true("nested Min flattening (lo)", full_lo == sp.Min(iota_lo, tij_lo, tik_lo, tjk_lo))
expect_true("nested Min flattening (hi)", full_hi == sp.Min(iota_hi, tij_hi, tik_hi, tjk_hi))

# ---------------------------------------------------------------------------
# III. Exact interval theorems verified by exhaustive finite order logic
# ---------------------------------------------------------------------------
subbanner("III. Exhaustive interval-theorem checks")

# 1) Local full-simplex interval theorem
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
                            raise AssertionError("local full-simplex interval theorem failed")
                        count_local += 1
print(f"verified local full-simplex interval theorem on {count_local} ordered integer samples")

# 2) Interior-certified and boundary-certified class tests
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
                                raise AssertionError("interior-certified ordering failed")
                            count_int += 1
                if iota_lo_v > beta_hi_v:
                    for b_star in range(beta_lo_v, beta_hi_v + 1):
                        for i_star in range(iota_lo_v, iota_hi_v + 1):
                            if not (b_star < i_star):
                                raise AssertionError("boundary-certified ordering failed")
                            count_bdry += 1
print(f"verified interior-certified local order on {count_int} samples")
print(f"verified boundary-certified local order on {count_bdry} samples")

# 3) Certified interval ordering between two triples
count_rank = 0
for L1 in range(0, 6):
    for U1 in range(L1, 6):
        for L2 in range(0, 8):
            for U2 in range(L2, 8):
                if U1 < L2:
                    for x in range(L1, U1 + 1):
                        for y in range(L2, U2 + 1):
                            if not (x < y):
                                raise AssertionError("primitive-triple ranking theorem failed")
                            count_rank += 1
print(f"verified certified interval ordering on {count_rank} samples")

# 4) Global up-to-three-coordinate splice theorem
count_global = 0
for pair_lo in range(0, 6):
    for pair_hi in range(pair_lo, 6):
        for tri_lo in range(0, 6):
            for tri_hi in range(tri_lo, 6):
                lo_global = min(pair_lo, tri_lo)
                hi_global = min(pair_hi, tri_hi)
                for p_star in range(pair_lo, pair_hi + 1):
                    for t_star in range(tri_lo, tri_hi + 1):
                        best_star = min(p_star, t_star)
                        if not (lo_global <= best_star <= hi_global):
                            raise AssertionError("up-to-three-coordinate splice theorem failed")
                        count_global += 1
print(f"verified global splice theorem on {count_global} samples")

# 5) Global three-coordinate improvement / no-improvement theorems
count_improve = 0
count_noimprove = 0
for pair_lo in range(1, 8):
    for pair_hi in range(pair_lo, 8):
        for tri_lo in range(0, 8):
            for tri_hi in range(tri_lo, 8):
                if tri_hi < pair_lo:
                    for p_star in range(pair_lo, pair_hi + 1):
                        for t_star in range(tri_lo, tri_hi + 1):
                            if not (t_star < p_star):
                                raise AssertionError("global three-coordinate improvement theorem failed")
                            count_improve += 1
                if tri_lo > pair_hi:
                    for p_star in range(pair_lo, pair_hi + 1):
                        for t_star in range(tri_lo, tri_hi + 1):
                            if not (p_star < t_star):
                                raise AssertionError("global three-coordinate no-improvement theorem failed")
                            count_noimprove += 1
print(f"verified global three-coordinate improvement theorem on {count_improve} samples")
print(f"verified global three-coordinate no-improvement theorem on {count_noimprove} samples")

# ---------------------------------------------------------------------------
# IV. Exact finite evaluation budget theorem
# ---------------------------------------------------------------------------
subbanner("IV. Exact finite evaluation budget theorem")

pair_eval_per_envelope = 6
pair_eval_total = len(pairs) * 2 * pair_eval_per_envelope
triple_eval_per_envelope = 24
triple_eval_total = len(triples) * 2 * triple_eval_per_envelope
full_eval_total = pair_eval_total + triple_eval_total

print(f"pairwise evaluation budget = {pair_eval_total}")
print(f"interior triple evaluation budget = {triple_eval_total}")
print(f"full up-to-three-coordinate budget = {full_eval_total}")

expect_zero("pairwise budget - 120", pair_eval_total - 120)
expect_zero("triple interior budget - 480", triple_eval_total - 480)
expect_zero("full budget - 600", full_eval_total - 600)

print("\nAll Stage-195 identities and interval theorems verified.")
