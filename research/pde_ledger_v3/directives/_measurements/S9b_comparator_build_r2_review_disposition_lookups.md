# Measurements — S9b comparator build review, round 2 (generated 2026-10-14 15:28)

Generator: `_scratch/s9b_build/gen/s9b_comparator_r2_lookups.sh`. Leg evidence: `_scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/` (copied from
`/tmp/s9b_cmp_review_r2_claudeleg/`) and `_scratch/s9b_build/s9b_comparator_build_review_r2_grok_evidence/` (copied from `/tmp/s9b_r2_*`). Operand values the legs computed
stay in that evidence. Line numbers are those of the reviewed comparator `_scratch/s9b_build/comparator_r2/S9b_cross_engine_comparator.py`.

## The reviewed artifacts
```
$ sha256sum -c _scratch/s9b_build/s9b_comparator_build_review_baseline_r2.sha256
_scratch/s9b_build/comparator_r2/S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_r2/test_S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_r2/S9b_comparator_build_report.md: OK
_scratch/s9b_build/comparator_r2/function_name_table.md: OK
_scratch/s9b_build/comparator_r2/coverage_reason_report.md: OK
_scratch/s9b_build/comparator_r2/validate_run.py: OK
_scratch/s9b_build/comparator_output/repair-02/full-01/comparator_source.py: OK
_scratch/s9b_build/comparator_output/repair-02/full-01/accounting.jsonl: OK
_scratch/s9b_build/comparator_output/repair-02/full-01/output_summary.json: OK
_scratch/s9b_build/comparator_output/repair-02/full-01/coverage_reason_inventory.json: OK
_scratch/s9b_build/comparator_output/repair-02/full-01/guard/outcome.json: OK
_scratch/s9b_build/comparator_output/repair-02/tests-04/comparator_source.py: OK
_scratch/s9b_build/comparator_output/repair-02/tests-04/test_source.py: OK
_scratch/s9b_build/comparator_output/repair-02/tests-04/guard/outcome.json: OK
research/pde_ledger_v3/directives/S9b_comparator_build_directive.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out: OK
research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out: OK
_scratch/s9b_build/comparator_output/repair-02/full-01/comparison.jsonl: OK
_scratch/s9b_build/s9b_comparator_build_review_prompt_r2.md: OK
```

## Both verdicts
```
$ head -1 _scratch/s9b_build/s9b_comparator_build_review_r2_claude.md
**Verdict:** NEEDS REVISION
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_comparator_build_review_r2_grok.txt
Verdict: CLEAR
```

```
$ cat _scratch/s9b_build/s9b_comparator_build_review_r2_grok.DONE
EXIT=0
```

## F1 (Claude 1): the row-level formed flag counts marker and label comparisons
```
$ sed -n '1214p;1468,1469p' _scratch/s9b_build/comparator_r2/S9b_cross_engine_comparator.py
    result['comparison_formed']=formed
        row['comparison_formed']=any(used)
        if not row['comparison_formed']: summary['rows_without_formed_comparison']+=1
```

```
$ grep -c '"comparison_formed":true' _scratch/s9b_build/comparator_output/repair-02/full-01/accounting.jsonl
293
```

```
$ grep -c '"comparison_formed":false' _scratch/s9b_build/comparator_output/repair-02/full-01/accounting.jsonl
37
```

```
$ grep '^{"name":"A_ONE_WAY_ER_G021"' _scratch/s9b_build/comparator_output/repair-02/full-01/accounting.jsonl
{"name":"A_ONE_WAY_ER_G021","accounting":"joined","py":{"parsed_leaves":29489,"compared_leaves":1,"uncompared_leaves":29488,"attempted_leaves":10631,"attempted_unformed_leaves":10630,"unattempted_leaves":18858},"wl":{"parsed_leaves":3315,"compared_leaves":1,"uncompared_leaves":3314,"attempted_leaves
```

```
$ sed -n 6p _scratch/s9b_build/comparator_r2/S9b_comparator_build_report.md | grep -o 'Rows with no formed comparison: [0-9]*'
Rows with no formed comparison: 37
```

The round-1 comparator (preserved at 1af2209c) has no row-level formed flag (match count):
```
$ git show 1af2209c:research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py | grep -c 'comparison_formed\|rows_without_formed'
0
```

```
$ cat _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/out/row_value_coverage.stdout
216 accounting comparison_formed=True formed kinds= ('marker', 'value')
54 accounting comparison_formed=True formed kinds= ('marker',)
37 accounting comparison_formed=False formed kinds= ()
22 accounting comparison_formed=True formed kinds= ('domain-formed', 'label')
1 accounting comparison_formed=True formed kinds= ('marker', 'value-boolean')
joined rows: 330  rows with a VALUE comparison formed: 217  rows flagged comparison_formed=true: 293
rows flagged formed but with no value comparison formed: ['A_ONE_WAY_ER_G021', 'A_ONE_WAY_ER_G101', 'A_ONE_WAY_ER_G120', 'A_ONE_WAY_ER_G121', 'A_ONE_WAY_RE_G021', 'A_ONE_WAY_RE_G101', 'A_ONE_WAY_RE_G120', 'A_ONE_WAY_RE_G121', 'A_ROUND_TRIP_G021', 'A_ROUND_TRIP_G101', 'A_ROUND_TRIP_G120', 'A_ROUND_TR
```

```
$ grep 'C_FIXED_RATIO_CONDITION_DEFLECTION' _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/ablate/analyze_s1.stdout | head -2
  C_FIXED_RATIO_CONDITION_DEFLECTION   computed products changed: False | comparison_formed True -> True
```

## F2 (Claude 2): the exponential profile never forms on the improper-integral rows
```
$ sed -n '1015p;1064p;1073p' _scratch/s9b_build/comparator_r2/S9b_cross_engine_comparator.py
def integral_screen(integral):
            leading=transformed.as_leading_term(t,cdir=1)
                raise UnresolvedEvaluation('endpoint convergence not established at '+math_text(endpoint))
```

Lines (one per row record) and occurrences:
```
$ grep -c 'endpoint convergence not established' _scratch/s9b_build/comparator_output/repair-02/full-01/comparison.jsonl
34
```

```
$ grep -c 'exponential: formed=0' _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/out/exp_profile_rows.stdout
34
```

```
$ head -3 _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/out/exp_profile_rows.stdout
rows with exponential samples unresolved for endpoint convergence: 34
  A_DEFLECTION_G001                  exponential: formed=0 endpoint=3 other=0 | inverse formed=3 | shifted formed=3
  A_DEFLECTION_G020                  exponential: formed=0 endpoint=3 other=0 | inverse formed=3 | shifted formed=3
```

```
$ grep -o '^p[0-9] [a-z]* *screen=[A-Za-z]*[^|]*' _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/ablate/runs/screen_check2/stdout
p0 inverse     screen=passed 
p0 shifted     screen=passed 
p0 exponential screen=UnresolvedEvaluation: endpoint convergence not established at 7/3 
p1 inverse     screen=passed 
p1 shifted     screen=passed 
p1 exponential screen=UnresolvedEvaluation: endpoint convergence not established at 11/5 
p2 inverse     screen=passed 
p2 shifted     screen=passed 
p2 exponential screen=UnresolvedEvaluation: endpoint convergence not established at 13/7 
```

```
$ grep -n 'several test profiles of different' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
60:   - For an object that is a functional of a live profile, evaluate under several test profiles of different
```

## F3 (Claude 3): no formed sample on 36 mixed-grade one-way and round-trip rows
```
$ sed -n '954p;1026p' _scratch/s9b_build/comparator_r2/S9b_cross_engine_comparator.py
        return {'unresolved':'held derivative, substitution or limit; evaluation not established for the original object'}
        raise UnresolvedEvaluation('trigonometric pole set not established')
```

Lines (one per row record) and occurrences:
```
$ grep -c 'held derivative, substitution or limit' _scratch/s9b_build/comparator_output/repair-02/full-01/comparison.jsonl
36
```

```
$ grep -o 'held derivative, substitution or limit' _scratch/s9b_build/comparator_output/repair-02/full-01/comparison.jsonl | wc -l
1260
```

```
$ grep -c 'trigonometric pole set not established' _scratch/s9b_build/comparator_output/repair-02/full-01/comparison.jsonl
33
```

```
$ grep -o 'trigonometric pole set not established' _scratch/s9b_build/comparator_output/repair-02/full-01/comparison.jsonl | wc -l
813
```

```
$ grep -E '(ONE_WAY_E?R?_?[ER]*|ROUND_TRIP)_G(021|101|120|121) .*on/value/\(value\)' _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/out/residual_census.stdout | grep -c "{'unresolved coverage': [0-9]*}$"
36
```

```
$ grep -E '(ONE_WAY_E?R?_?[ER]*|ROUND_TRIP)_G(021|101|120|121) .*on/value/\(value\)' _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/out/residual_census.stdout | head -3
   A_ONE_WAY_ER_G021 /on/value/(value) {'unresolved coverage': 33}
   A_ONE_WAY_ER_G101 /on/value/(value) {'unresolved coverage': 33}
   A_ONE_WAY_ER_G120 /on/value/(value) {'unresolved coverage': 33}
```

```
$ grep -n 'A stratum the comparator cannot evaluate is printed as' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
58:     which a parameter, a profile or a density vanishes. A stratum the comparator cannot evaluate is printed as
```

```
$ grep -o '36 one-way and round-trip values at G021, G101, G120 and G121 are never evaluated on a live' research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r1_review_disposition.md
36 one-way and round-trip values at G021, G101, G120 and G121 are never evaluated on a live
```

## O1 (Claude, checked): the leaf cache key has no control
```
$ grep -n 'residual formed' _scratch/s9b_build/s9b_comparator_build_review_r2_claude_evidence/ablate/analyze_s2_cache_corrupt.stdout | head -2
5:      residual formed                                                           37 ->     0
6:      residual formed, |computed_30| > absolute_tolerance                       36 ->     0
```

## Grok's corresponding checks
```
$ grep -o 'rows_without_formed_comparison. 37' _scratch/s9b_build/s9b_comparator_build_review_r2_grok.txt
rows_without_formed_comparison` 37
```

```
$ grep -o 'Those are engine decomposition differences, printed as such.' _scratch/s9b_build/s9b_comparator_build_review_r2_grok.txt
Those are engine decomposition differences, printed as such.
```

