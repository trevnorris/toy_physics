# Measurements — S9b directive amendment 10, review round 2 (generated 2026-10-13 06:49)

Generator: `_scratch/s9b_build/gen/s9b_amend10_r2_lookups.sh`. The legs' scripts and stdout are copied (files under
2 MB) to `_scratch/s9b_build/s9b_amend10_review_r2_codex_evidence/` (from `/tmp/s9b_amend10_r2_review_20261013/`) and `_scratch/s9b_build/s9b_amend10_review_r2_grok_evidence/` (from
`/tmp/s9b-amend10-r2-physics/`).

## The reviewed text is the accepted text
```
$ sha256sum research/pde_ledger_v3/directives/S9b_repair_build_directive.md
faf60e3cdea94630b30b1f50aac1faace269260c43b959074ecedfaf7b801144  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend10_review_baseline_r2.sha256
faf60e3cdea94630b30b1f50aac1faace269260c43b959074ecedfaf7b801144  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

## Both verdicts
```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend10_review_r2_codex_final.txt
Verdict: CLEAR
```

```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend10_review_r2_grok.txt
Verdict: CLEAR
```

```
$ cat _scratch/s9b_build/s9b_amend10_review_r2_codex.DONE _scratch/s9b_build/s9b_amend10_review_r2_grok.DONE
EXIT=0
EXIT=0
```

## The round-0 and round-1 preservation commits, and the seconds the directive path held them
```
$ git log -3 --format='%h %ci %s' 192adc31
192adc31 2026-10-13 06:16:01 -0600 Preserve S9b directive amendment 10 as reviewed in round 1 (not accepted)
331e748d 2026-10-13 06:16:01 -0600 Preserve S9b directive amendment 10 as reviewed in round 0 (not accepted)
dfa07223 2026-10-13 06:15:24 -0600 Preserve the S9b repair build after repair round 6 (review r6 baseline; not accepted)
```

```
$ ls -la --time-style=+%H:%M:%S _scratch/s9b_build/s9b_amend10_review_r2_codex_evidence/final_file_evidence.txt _scratch/s9b_build/s9b_amend10_review_r2_codex_evidence/attempts.jsonl
-rw-rw-r-- 1 trevnorris trevnorris 13304 06:42:32 _scratch/s9b_build/s9b_amend10_review_r2_codex_evidence/attempts.jsonl
-rw-rw-r-- 1 trevnorris trevnorris 24847 06:47:24 _scratch/s9b_build/s9b_amend10_review_r2_codex_evidence/final_file_evidence.txt
```

```
$ grep -n 'S9b_repair_build_directive.md: ' _scratch/s9b_build/s9b_amend10_review_r2_codex_evidence/final_file_evidence.txt
2:research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

```
$ ls -la --time-style=+%H:%M:%S /tmp/s9b-amend10-r2-physics/evidence/
total 20
drwxrwxr-x 2 trevnorris trevnorris 4096 06:29:56 .
drwxrwxr-x 3 trevnorris trevnorris 4096 06:45:10 ..
-rw-rw-r-- 1 trevnorris trevnorris 6528 06:29:56 amendment10.diff
-rw-rw-r-- 1 trevnorris trevnorris 2124 06:29:56 diffs.txt
```

```
$ grep -c 'Two sets differ at an element that one contains' _scratch/s9b_build/s9b_amend10_review_r2_grok_evidence/evidence/amendment10.diff
1
```

## Grok's nested-set and binding stdout
```
$ cat _scratch/s9b_build/s9b_amend10_review_r2_grok_evidence/nested_sets.stdout
SYMPY 1.14.0
B1_FREE_SYMBOLS [a] bound_in_free False
B1_MEMBERSHIP True False print CHANGE
ROOTS_FREE []
ROOTS_CONTAINS_TRUE True
ROOTS_CONTAINS_FALSE Eq(ConditionSet(s9b_path, s9b_path <= 2, Interval(1, 10)), EmptySet)
GATE_TRUE_COUNTEREXAMPLES [] SAMPLES [True, True, True]
GATE_FALSE_COUNTEREXAMPLES [1] SAMPLES [False, True, True]
ROOTS_CHANGED_BY_CERTIFICATE True False print CHANGE
TWO_FALSE_GATES_COUNTEREXAMPLES [1, 5.5, 10] [1]
EMPTINESS_OF_ROOTS_BOTH_FALSE_GATES roots_nonempty False False emptiness_statements_differ False inner_counterexample_sets_differ True print_of_emptiness_part ZERO print_if_inner_set_bubbled CHANGE
TRAVERSED_MEMBERSHIP_BASE {1/2: False, 2: True, 6: True}
TRAVERSED_MEMBERSHIP_CORRUPT {1/2: False, 2: False, 6: True}
TRAVERSED_DIFFER_AT [2] print CHANGE
TRAVERSED_SAME_DESPITE_INNER_ROOTS True inner_roots_differ True print_outer ZERO print_if_inner_bubbled CHANGE
REAL_DOMAIN_UNCHANGED_WHEN_NESTED_TRAVERSED_CHANGES True True counters [] [] nested_traversed_differ True print ZERO
REAL_DOMAIN_CHANGED_BY_PREDICATE True False counters [] [3/2, 2] certificate_radius_not_a_free_coordinate [3/2, 2] print CHANGE
SAME_SYMBOL_FREE [] bound_symbol_still_free False
CAPTURED_SUBS_OF_BOUND_SYMBOL Eq(EmptySet, EmptySet)
CAPTURED_EQUALS_ORIGINAL True
INTERVAL_DOMAIN True False print CHANGE
INERT_STRUCTURAL_EQ False
INERT_DOIT_EQ True
INERT_PRINT ZERO
INEXACT_EQUAL_VALUE delta 3/2 exceeds_each_single_bound True exceeds_sum_bound False print_under_sum_bound UNDECIDED
AGREEING_NUMERIC -0.e-124 equals_zero_float False
AGREEING_DOES_NOT_PRINT_CHANGE UNDECIDED
SYMBOLIC_ZERO_SHOWS_UNCHANGED ZERO
ADMISSIBLE_PROFILE derivative 3 meets_negative_condition False values (7, 0) print CHANGE
EXCLUSIVITY [('equal_shown', 'ZERO'), ('changed_shown', 'CHANGE'), ('neither', 'UNDECIDED'), ('both_shown_zero_wins', 'ZERO')]
```

```
$ cat _scratch/s9b_build/s9b_amend10_review_r2_grok_evidence/binding.stdout
WL_START
INERT_SAMEQ False
INERT_ACTIVATE_SAMEQ True
INERT_ACTIVATED 0
FORALL_BOTH_FALSE {False, False}
MATRIX_AT_3_2 {True, False}
FORALL_TRUTH {True, False}
NESTED_EXISTS_TRUTH {False, True}
NESTED_EXISTS_SAME_TRUTH {False, False}
INNER_EXISTS_AT_R3 {True, True}
REDUCE_NO_REAL_VELOCITY yRoot < 0
REDUCE_STILL_HAS_EXISTS False
ELEMENT_VS_TRUE {False, True}
ELEMENT_DIFFERS_AT_I {True, False}
WL_DONE
```

