# Measurements — amendment 11 review, round 8 (generated 2026-10-16 08:28)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r8_lookups.sh` (sha256sum/grep/sed/cat/wc only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r8.md` are the reviewed round-8 text's (saved before launch). The legs' stdout is quoted from their evidence
directories.

## The reviewed round-8 text, its scope, and the baseline after both legs
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r8.sha256
930999fcf49d0c90286c9708d628a1aea897af94e288b332bb476b0a68dcbb84  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r8.md research/pde_ledger_v3/directives/S9b_repair_build_directive.md
930999fcf49d0c90286c9708d628a1aea897af94e288b332bb476b0a68dcbb84  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r8.md
930999fcf49d0c90286c9708d628a1aea897af94e288b332bb476b0a68dcbb84  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ git diff -U0 899458a1 -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep '^@@'
@@ -306,4 +306,5 @@ comparator, the production runs and the record.
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r8.sha256 | grep -c ': OK$'
22
```

```
$ git diff --stat 9af9f705 -- research/pde_ledger_v3/scripts research/pde_ledger_v3/mathematica | wc -l
0
```

## Verdicts
```
$ grep -m1 'Verdict' _scratch/s9b_build/s9b_amend11_review_r8_claude.md
**Verdict:** CLEAR
```

```
$ grep -m1 'Findings' _scratch/s9b_build/s9b_amend11_review_r8_claude.md
**Findings:** None.
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r8_grok.txt
Verdict: CLEAR
```

```
$ grep -c 'No findings' _scratch/s9b_build/s9b_amend11_review_r8_grok.txt
1
```

## The pairing paragraph and the legs' checks of it
```
$ sed -n '306,310p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r8.md
      For value witnesses and compared domain conditions, compute and print the intersection of the solutions each
      copy admits at the witness inputs, irrespective of whether the emitted fixing equations or constructions are
      the same expression. If the intersection is nonempty, compare at a common solution, or print and compare the
      complete sets of (solution, value) pairs, using condition evaluations for domains. Pair different solutions
      only when the intersection is decided empty; keep each copy's own solutions.
```

```
$ grep 'RULE\|AT_COMMON\|BARRED\|SAME_SET\|SAME_EXPRESSION' _scratch/s9b_build/s9b_amend11_review_r8_claude_evidence/py_pairing.stdout
SQUARED_SAME_EXPRESSION_AS_EMITTED False
RULE BRANCH 0 BASE EmptySet K4C EmptySet INTERSECTION EmptySet
RULE BRANCH 1 BASE {2, 4} K4C {2, 4} INTERSECTION {2, 4}
  AT_COMMON z 2 BASE -3*sqrt(2)/4 K4C -3*sqrt(2)/4 RESIDUAL 0
  AT_COMMON z 4 BASE -3/16 K4C -3/16 RESIDUAL 0
  COMPLETE_SETS BASE [(2, -3*sqrt(2)/4), (4, -3/16)] K4C [(2, -3*sqrt(2)/4), (4, -3/16)] SAME_SET True
  BARRED_DIFFERENT_SOLUTION_PAIR BASE z 2 -3*sqrt(2)/4 K4C z 4 -3/16 RESIDUAL 3/16 - 3*sqrt(2)/4
  BARRED_DIFFERENT_SOLUTION_PAIR BASE z 4 -3/16 K4C z 2 -3*sqrt(2)/4 RESIDUAL -3/16 + 3*sqrt(2)/4
RULE BRANCH 2 BASE {2, 4} K4C {2, 4} INTERSECTION {2, 4}
  AT_COMMON z 2 BASE 3*sqrt(2)/4 K4C 3*sqrt(2)/4 RESIDUAL 0
  AT_COMMON z 4 BASE 3/16 K4C 3/16 RESIDUAL 0
  COMPLETE_SETS BASE [(2, 3*sqrt(2)/4), (4, 3/16)] K4C [(2, 3*sqrt(2)/4), (4, 3/16)] SAME_SET True
  BARRED_DIFFERENT_SOLUTION_PAIR BASE z 2 3*sqrt(2)/4 K4C z 4 3/16 RESIDUAL -3/16 + 3*sqrt(2)/4
  BARRED_DIFFERENT_SOLUTION_PAIR BASE z 4 3/16 K4C z 2 3*sqrt(2)/4 RESIDUAL 3/16 - 3*sqrt(2)/4
RULE BRANCH 3 BASE EmptySet K4C EmptySet INTERSECTION EmptySet
RULE BRANCH 0 BASE EmptySet K4C EmptySet INTERSECTION EmptySet
RULE BRANCH 1 BASE {3, 5} K4C EmptySet INTERSECTION EmptySet
RULE BRANCH 2 BASE {3, 5} K4C EmptySet INTERSECTION EmptySet
RULE BRANCH 3 BASE EmptySet K4C EmptySet INTERSECTION EmptySet
```

```
$ grep 'RULE\|AT_COMMON\|BARRED\|SAME_EXPRESSION\|HAS_XI' _scratch/s9b_build/s9b_amend11_review_r8_claude_evidence/wl_pairing.stdout
SQUARED_VELOCITY_SAME_EXPRESSION_AS_EMITTED False
BASE_SQUARED_HAS_XI_SLOPE True
K4C_SQUARED_HAS_XI_SLOPE False
RULE INTERSECTION {2, 4}
  AT_COMMON case 1 rBoundary 2 BASE {-3/(2*Sqrt[2]), 3/(2*Sqrt[2])} K4C {-3/(2*Sqrt[2]), 3/(2*Sqrt[2])} RESIDUAL {0, 0}
  AT_COMMON case 1 rBoundary 4 BASE {-3/16, 3/16} K4C {-3/16, 3/16} RESIDUAL {0, 0}
  BARRED_DIFFERENT_SOLUTION_PAIR case 1 BASE 2 {-3/(2*Sqrt[2]), 3/(2*Sqrt[2])} K4C 4 {-3/16, 3/16} RESIDUAL {(3*(1 - 4*Sqrt[2]))/16, (3*(-1 + 4*Sqrt[2]))/16}
  BARRED_DIFFERENT_SOLUTION_PAIR case 1 BASE 4 {-3/16, 3/16} K4C 2 {-3/(2*Sqrt[2]), 3/(2*Sqrt[2])} RESIDUAL {(3*(-1 + 4*Sqrt[2]))/16, (3*(1 - 4*Sqrt[2]))/16}
  AT_COMMON case 2 rBoundary 2 BASE {3/(2*Sqrt[2]), -3/(2*Sqrt[2])} K4C {3/(2*Sqrt[2]), -3/(2*Sqrt[2])} RESIDUAL {0, 0}
  AT_COMMON case 2 rBoundary 4 BASE {3/16, -3/16} K4C {3/16, -3/16} RESIDUAL {0, 0}
  BARRED_DIFFERENT_SOLUTION_PAIR case 2 BASE 2 {3/(2*Sqrt[2]), -3/(2*Sqrt[2])} K4C 4 {3/16, -3/16} RESIDUAL {(3*(-1 + 4*Sqrt[2]))/16, (3*(1 - 4*Sqrt[2]))/16}
  BARRED_DIFFERENT_SOLUTION_PAIR case 2 BASE 4 {3/16, -3/16} K4C 2 {3/(2*Sqrt[2]), -3/(2*Sqrt[2])} RESIDUAL {(3*(1 - 4*Sqrt[2]))/16, (3*(-1 + 4*Sqrt[2]))/16}
```

```
$ grep 'RULE_ON\|ADMITTED\|COMPLETE_SETS\|INTERSECTION' _scratch/s9b_build/s9b_amend11_review_r8_grok_evidence/pairing_rule.stdout
ADMITTED_ZEROS_BASE: {2, 4}
ADMITTED_ZEROS_K4C: {2, 4}
RULE_ON_CROSSED: {'intersection': [2, 4], 'decided_empty': False, 'compare_at_this_common_element': False, 'pair_these_different_solutions': False, 'compare_complete_sets': True}
RULE_ON_COMMON: {'intersection': [2, 4], 'decided_empty': False, 'compare_at_this_common_element': True, 'pair_these_different_solutions': False, 'compare_complete_sets': True}
COMPLETE_SETS_EQUAL: True
TURN_baseline_ADMITTED: [3]
TURN_k3_ADMITTED: [3/2 + sqrt(21)/2]
TURN_ADMITTED_INTERSECTION_EMPTY: True
TURN_ALGEBRAIC_INTERSECTION: [1]
RULE_ON_DIFFERENT_ADMITTED_TURNING_RADII: {'intersection': [], 'decided_empty': True, 'compare_at_this_common_element': False, 'pair_these_different_solutions': True, 'compare_complete_sets': True}
RULE_IF_FAILED_ROOT_COUNTS_AS_ADMITTED: {'intersection': [1], 'decided_empty': False, 'compare_at_this_common_element': False, 'pair_these_different_solutions': False, 'compare_complete_sets': True}
```

## Word counts
```
$ sed -n '/^14\. \*\*Ablation harness/,/^15\. /p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r8.md | wc -w
1353
```

