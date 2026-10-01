#!/usr/bin/env python3
"""Leg script 10: does restricting to R1 (spurion along direction 1 only) and class P (fields
independent of x3) BEFORE the §3a quotient-rank test change the number of independent
invariants?  Uses the S11c-b SymPy engine's own enumeration/signature/rank functions from a
read-only COPY (engine_copy/), with the S11c-a exports copy on the path.

Directive B §2.2 says: implement the group-fixed domain in each engine BEFORE the energy is
constructed; B §1 says keep the accepted 40-record object (10 + 15 + 15) as the coefficient
domain.  Print: uniform rank and per-source spurion rank, generic vs restricted.
Prints computed counts only.
"""
import sys
sys.path.insert(0, "/tmp/s11cd_clean_review_r4_claude/engine_copy")
sys.path.insert(0, "/tmp/s11cd_clean_review_r4_claude/import_probe")
sys.setrecursionlimit(100000)
import sympy as sp
import S11c_b_brane_operator_sympy_audit as E

def restriction_subs(restrict_spurion_to_1: bool, class_p: bool):
    subs = {}
    if restrict_spurion_to_1:
        subs[E.bg[1]] = 0
        subs[E.bg[2]] = 0
        for (i, j), sym in E.basis_background_second.items():
            if i != 0 or j != 0:
                subs[sym] = 0
    if class_p:
        for field, first, second in E.basis_fields:
            subs[first[2]] = 0
            for (i, j), sym in second.items():
                if i == 2 or j == 2:
                    subs[sym] = 0
    return subs

def rank_of(candidates, signatures, subs):
    restricted_c = tuple(sp.expand(c.subs(subs, simultaneous=True)) for c in candidates)
    restricted_s = tuple(tuple(sp.expand(e.subs(subs, simultaneous=True)) for e in sig) for sig in signatures)
    # drop candidates whose restricted signature is identically zero (they would be null on the domain)
    keep = [k for k, sig in enumerate(restricted_s) if any(e != 0 for e in sig)]
    if not keep:
        return 0, len(candidates), 0
    sel, _ = E.quotient_independent_indices(tuple(restricted_c[k] for k in keep), tuple(restricted_s[k] for k in keep))
    return len(sel), len(candidates), len(candidates) - len(keep)

uniform_c = tuple(expr for _, expr in E.UNIFORM_CANDIDATES)
uniform_s = E.UNIFORM_SIGNATURES
for label, (r1, p) in (("GENERIC", (False, False)), ("CLASS_P_ONLY", (False, True)), ("R1_AND_CLASS_P", (True, True))):
    rk, n, nnull = rank_of(uniform_c, uniform_s, restriction_subs(r1, p))
    print("UNIFORM", label, "CANDIDATES", n, "NULL_ON_DOMAIN", nnull, "QUOTIENT_RANK", rk)

cands = E.enumerate_new_candidates(tuple(E.bg))
exprs = tuple(expr for _, expr in cands)
sigs = E.basis_euler_signatures(exprs, E.basis_fields, background_first_jets=E.bg,
                                background_second_jets=E.basis_background_second,
                                background_depth=E.STRONG_ROW_JET_DEPTH)
for label, (r1, p) in (("GENERIC", (False, False)), ("R1_ONLY", (True, False)), ("CLASS_P_ONLY", (False, True)), ("R1_AND_CLASS_P", (True, True))):
    rk, n, nnull = rank_of(exprs, sigs, restriction_subs(r1, p))
    print("SPURION_PER_SOURCE", label, "CANDIDATES", n, "NULL_ON_DOMAIN", nnull, "QUOTIENT_RANK", rk)
