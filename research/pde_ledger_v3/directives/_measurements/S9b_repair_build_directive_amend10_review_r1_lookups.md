# Measurements — S9b directive amendment 10, review round 1 (generated 2026-10-13 02:21)

Generator: `_scratch/s9b_build/gen/s9b_amend10_r1_lookups.sh`. The legs' scripts and stdout are copied (files under
2 MB) to `_scratch/s9b_build/s9b_amend10_review_r1_codex_evidence/` (from `/tmp/s9b-amend10-r1-review-4qe91gs_/`) and `_scratch/s9b_build/s9b_amend10_review_r1_grok_evidence/` (from
`/tmp/s9b-amend10-r1-physics/`).

## The reviewed text is preserved byte-exact
```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r1.md
5da3c7e2da97d67ae8f5cf83124a319f6e98adda3efe7d89118b13d7121a0797  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r1.md
```

```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend10_review_baseline_r1.sha256
5da3c7e2da97d67ae8f5cf83124a319f6e98adda3efe7d89118b13d7121a0797  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

## Both verdicts
```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend10_review_r1_codex_final.txt
Verdict: CLEAR
```

```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend10_review_r1_grok.txt
Verdict: NEEDS REVISION
```

```
$ cat _scratch/s9b_build/s9b_amend10_review_r1_codex.DONE _scratch/s9b_build/s9b_amend10_review_r1_grok.DONE
EXIT=0
EXIT=0
```

## B1: the round-1 sentence
```
$ grep -n -A1 'A condition that binds variables, in a set or a quantifier, is' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r1.md
286:      such a point under some admissible profile. A condition that binds variables, in a set or a quantifier, is
287-      compared as the condition it states on its free variables, so its point gives values to those only.
```

## B1: Codex's set lines
```
$ grep -n 'R6_SET_DIFFERENCE\|BOUND_SET_MEMBERSHIP_AT_3_OVER_2_AFTER_A_1\|^BOUND_SET_OPERANDS\|^MEMBERSHIP_EVALUATIONS' _scratch/s9b_build/s9b_amend10_review_r1_codex_evidence/report.txt
2484:BOUND_SET_OPERANDS: (ConditionSet(r, r > a, Reals), ConditionSet(r, r > a + 1, Reals))
2487:BOUND_SET_MEMBERSHIP_AT_3_OVER_2_AFTER_A_1: (3/2 > a, 3/2 > a + 1)
2492:R6_SET_DIFFERENCE: ((VALUE_DIFFERENCE, (UNDECIDED, (BASELINE, ConditionSet(r, r > a, Reals)), (CORRUPTED, ConditionSet(r, r > a + 1, Reals)), (RESIDUAL, SymmetricDifference(ConditionSet(_0, _0 - a > 0, Reals), ConditionSet(_0, _0 - a - 1 > 0, Reals))))), (DOMAIN_DIFFERENCE, 0))
2524:BOUND_SET_OPERANDS: (ConditionSet(r, r > a, Reals), ConditionSet(r, r > a + 1, Reals))
2527:MEMBERSHIP_EVALUATIONS: (True, False)
2574:BOUND_SET_OPERANDS: {ImplicitRegion[r > a, {r}], ImplicitRegion[r > 1 + a, {r}]}
2576:BOUND_SET_MEMBERSHIP_AT_3_OVER_2_AFTER_A_1: {True, False}
```

## B1: Grok's set_decision stdout
```
$ cat _scratch/s9b_build/s9b_amend10_review_r1_grok_evidence/set_decision.py.stdout
equal_structural False
equal_is_subset_ab None
equal_is_subset_ba None
equal_complement Complement(ConditionSet(r, (r > 0) & (r > 1), Reals), ConditionSet(r, r > 1, Reals))
equal_rel_a (r > 0) & (r > 1)
equal_rel_b r > 1
differ_structural False
differ_is_subset_ab None
differ_is_subset_ba None
differ_complement Complement(ConditionSet(r, (r > 1) & (r < 3), Reals), ConditionSet(r, (r > 2) & (r < 4), Reals))
differ_rel_a (r > 1) & (r < 3)
differ_rel_b (r > 2) & (r < 4)
differ_contains_3_2 (True, False)
subs_repr ConditionSet(r, r > b, Interval(0, oo))
xreplace_repr ConditionSet(r, r > 1, Interval(0, oo))
xreplace_structural False
xreplace_is_subset_ab None
xreplace_is_subset_ba None
xreplace_complement Complement(ConditionSet(r, r > 1, Interval(0, oo)), ConditionSet(r, r > 2, Interval(0, oo)))
xreplace_rel_a (0 <= r) & (r > 1)
xreplace_rel_b (0 <= r) & (r > 2)
xreplace_contains (True, False)
bool_structural False
bool_P_subset_empty None
bool_Q_subset_empty None
bool_Q_contains_3_4 True
bool_P_contains_3_4 False
```

## B1: set-valued objects in the SymPy engine
```
$ grep -n 'roots=sp.ConditionSet\|traversed=sp.ConditionSet\|traversed=sp.Union' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
569:        roots=sp.ConditionSet(sp.Tuple(t,j),root_condition,
572:        traversed=sp.ConditionSet(r,sp.Ne(below,sp.EmptySet,evaluate=False),sp.Interval(0,upper,left_open=True))
695:    traversed=sp.Union(*(ray['radii'] for ray in rays.values()),evaluate=False)
```

## Grok's guarded run directories
```
$ ls -la --time-style=+%H:%M _scratch/s11c/s9b-amend10-review-r1/
total 64
drwxrwxr-x 16 trevnorris trevnorris 4096 01:45 .
drwxrwxr-x 56 trevnorris trevnorris 4096 23:37 ..
drwxrwxr-x  2 trevnorris trevnorris 4096 23:51 binding-py
drwxrwxr-x  2 trevnorris trevnorris 4096 23:52 binding-py-2
drwxrwxr-x  2 trevnorris trevnorris 4096 23:53 binding-py-3
drwxrwxr-x  2 trevnorris trevnorris 4096 23:56 binding-py-4
drwxrwxr-x  2 trevnorris trevnorris 4096 23:44 s9b-amend10-r1-review-4qe91gs_-rule-1
drwxrwxr-x  2 trevnorris trevnorris 4096 23:37 s9b-amend10-r1-review-4qe91gs_-sympy-1
drwxrwxr-x  2 trevnorris trevnorris 4096 23:40 s9b-amend10-r1-review-4qe91gs_-sympy-2
drwxrwxr-x  2 trevnorris trevnorris 4096 23:43 s9b-amend10-r1-review-4qe91gs_-sympy_bound_conditions-1
drwxrwxr-x  2 trevnorris trevnorris 4096 23:40 s9b-amend10-r1-review-4qe91gs_-wolfram-1
drwxrwxr-x  2 trevnorris trevnorris 4096 23:43 s9b-amend10-r1-review-4qe91gs_-wolfram_final-1
drwxrwxr-x  2 trevnorris trevnorris 4096 23:42 s9b-amend10-r1-review-4qe91gs_-wolfram_v2-1
drwxrwxr-x  2 trevnorris trevnorris 4096 00:00 sat-bound-py
drwxrwxr-x  2 trevnorris trevnorris 4096 01:45 sat-bound-py-2
drwxrwxr-x  2 trevnorris trevnorris 4096 23:58 set-decision-py
```

## The round-2 text
```
$ grep -n -A3 'Two sets differ at an element that one contains' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
286:      such a point under some admissible profile. Two sets differ at an element that one contains and the other does
287-      not, so a point for a set gives values to its free variables and to that element. A statement that binds
288-      variables, such as a quantifier or a test of whether a set is empty, is compared as the statement it makes on
289-      its free variables, so its point gives values to those only.
```

```
$ grep -c 'A condition that binds variables, in a set or a quantifier' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
0
```

