# Measurements — S9b comparator build review, round 4 (generated 2026-10-15 23:36)

Generator: `_scratch/s9b_build/gen/s9b_comparator_r4_lookups.sh`. Leg evidence: `_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/` (copied from
`/tmp/s9b_cmp_review_r4_claudeleg/`) and `_scratch/s9b_build/s9b_comparator_build_review_r4_grok_evidence/` (copied from `/tmp/s9b-r4-review/`). Operand values the legs computed
stay in that evidence. Line numbers are those of the reviewed comparator `_scratch/s9b_build/comparator_r4/S9b_cross_engine_comparator.py`. Earlier versions: round 0
`3ec82aa9`, round 1 `1af2209c`, round 2 `75676d01`, round 3 `52698c77`.

## The reviewed artifacts
```
$ sha256sum -c _scratch/s9b_build/s9b_comparator_build_review_baseline_r4.sha256 | grep -c ': OK$'
33
```

```
$ wc -l < _scratch/s9b_build/s9b_comparator_build_review_baseline_r4.sha256
33
```

```
$ sha256sum _scratch/s9b_build/comparator_r4/S9b_cross_engine_comparator.py _scratch/s9b_build/comparator_r4/test_S9b_cross_engine_comparator.py _scratch/s9b_build/s9b_comparator_build_review_prompt_r4.md
326d204bcf76f81d40c141803106e42aa8bbb71da91414aef07a3cac710bb1f7  _scratch/s9b_build/comparator_r4/S9b_cross_engine_comparator.py
efe7249bc18698fedb9d365aa1816a17014ab5967737af680744a5e31c989ad8  _scratch/s9b_build/comparator_r4/test_S9b_cross_engine_comparator.py
5479dab51b508381139484fcbad91234acc244e1e7fdb7ca2f3d212326669ea1  _scratch/s9b_build/s9b_comparator_build_review_prompt_r4.md
```

## Both verdicts
```
$ head -1 _scratch/s9b_build/s9b_comparator_build_review_r4_claude.md
## Verdict: CLEAR
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_comparator_build_review_r4_grok.txt
Verdict: CLEAR
```

```
$ grep -c '^## Findings' _scratch/s9b_build/s9b_comparator_build_review_r4_claude.md; grep -A1 '^## Findings' _scratch/s9b_build/s9b_comparator_build_review_r4_claude.md
1
## Findings
None.
```

```
$ grep -A2 '^## Findings' _scratch/s9b_build/s9b_comparator_build_review_r4_grok.txt
## Findings

None.
```

## Claude observation 1: quadrature formation where an integral cancels to about zero
```
$ sed -n '1205,1208p' _scratch/s9b_build/comparator_r4/S9b_cross_engine_comparator.py
        lo,hi=bounds()
        middle=(lo+hi)/2 if mp.isfinite(lo) and mp.isfinite(hi) else lo+1 if mp.isfinite(lo) else hi-1 if mp.isfinite(hi) else mp.mpf(0)
        second,error2=mp.quad(normalized,[lo,middle,hi],method='tanh-sinh',error=True)
        tolerance=mp.power(10,-digits+4)*abs(second)
```

The round-3 requirement (C2), as dispositioned:
```
$ grep -o 'Required: whether a quadrature sample forms does not depend on the magnitude of the integral[^.]*\.' research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r3_review_disposition.md
Required: whether a quadrature sample forms does not depend on the magnitude of the integral; a value established at the printed precision forms.
```

```
$ grep -h 'exact-zero by cancellation' _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/probe_quad/probe_quad.stdout
exact-zero by cancellation, int_0^2 (r-1) dr: outcome=symbolic samples={'unresolved coverage': 33}
```

Occurrences in the reviewed production output of the reason this probe prints, and of the two reasons the leg
lists for the unresolved value samples:
```
$ grep -o 'quadrature value not established' _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl | wc -l
0
```

```
$ grep -o 'interior denominator zero at' _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl | wc -l
183
```

```
$ grep -o 'nonpolynomial denominator zero set not established' _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl | wc -l
27
```

```
$ grep -h "^value-sample outcomes\|brane_zero_at_point {\|^formed but reliability" _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/out/value_samples.stdout
value-sample outcomes: {'residual formed': 5786, 'unresolved coverage': 105}
   brane_zero_at_point {'residual formed': 408, 'unresolved coverage': 105}
formed but reliability_established false: 0
```

## Claude observation 2: the tests and the WL reader's precedence and chained comparisons
Test totals: baseline, then the six ablations F1–F6 in order:
```
$ for f in _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/tests_run/tests_baseline.stderr _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F1/run/tests.stderr _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F2/run/tests.stderr _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F3/run/tests.stderr _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F4/run/tests.stderr _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F5/run/tests.stderr _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F6/run/tests.stderr; do echo "$f: $(grep -h '^OK\|^FAILED' $f)"; done
_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/tests_run/tests_baseline.stderr: OK
_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F1/run/tests.stderr: OK
_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F2/run/tests.stderr: OK
_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F3/run/tests.stderr: FAILED (failures=2)
_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F4/run/tests.stderr: FAILED (failures=3)
_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F5/run/tests.stderr: FAILED (failures=1)
_scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F6/run/tests.stderr: OK
```

F1 on the real streams, rows whose records differ from the unablated run (match count), and F2's:
```
$ grep -c 'records identical = False' _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F1/run/base1_vs_F1.stdout
4
```

```
$ grep -c 'records identical = True' _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F1/run/base1_vs_F1.stdout
0
```

```
$ grep -h 'records identical' _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F2/run/base1_vs_F2.stdout | sort | uniq -c
      1 == A_DEFLECTION_G100: records identical = False
      1 == A_RADAR_LOG_G020: records identical = False
      1 == A_ROUND_TRIP_G020: records identical = False
      1 == F_FLOW_A_ONE_WAY_ER_G020: records identical = False
```

Lines reporting a moved residual, F1 then F2 (match counts):
```
$ grep -c 'residual moved' _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F1/run/base1_vs_F1.stdout
3
```

```
$ grep -c 'residual moved' _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F2/run/base1_vs_F2.stdout
0
```

```
$ cat _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/ablate/F2/evidence/count_chains.stdout
chained comparison sites (same level, >=2 relations): 3276 in 333 tags
   66 WL_S9B_C_FIXED_RATIO_IMPLIED_JN_DEFLECTION
   66 WL_S9B_C_FIXED_RATIO_IMPLIED_JN_RADAR
   66 WL_S9B_C_POWER_IMPLIED_JN_DEFLECTION
   66 WL_S9B_C_POWER_IMPLIED_JN_RADAR
   36 WL_S9B_C_FIXED_RATIO_CONDITION_DEFLECTION
   36 WL_S9B_C_FIXED_RATIO_CONDITION_RADAR
   36 WL_S9B_C_POWER_CONDITION_DEFLECTION
   36 WL_S9B_C_POWER_CONDITION_RADAR
   32 WL_S9B_B_IMPLIED_JN_DEFLECTION
   32 WL_S9B_B_IMPLIED_JN_RADAR
   32 WL_S9B_C_CONSTANT_IMPLIED_JN_DEFLECTION
   32 WL_S9B_C_CONSTANT_IMPLIED_JN_RADAR
Inequality[ heads: 0 in 0 tags
```

The density-label rows' cited engine lines (main tree), on which the table's orientation rests (F6):
```
$ sed -n '91,92p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
                ('RHO_POSITIVE',sp.Gt(density,0)),
                ('RHO_NONPOSITIVE',sp.Le(density,0)),
```

```
$ sed -n '620p;622p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  "DensityParts" -> {<|"Label" -> "rho > 0", "LabelCondition" -> (rho0 (1+fProfile[r]) > 0),
    <|"Label" -> "rho <= 0", "LabelCondition" -> (rho0 (1+fProfile[r]) <= 0),
```

## Claude observation 3: order of the domain-stratum inventory
```
$ grep -h 'multiset equal' _scratch/s9b_build/s9b_comparator_build_review_r4_claude_evidence/runs/analysis/strata_sets_base1.stdout
A_DEFLECTION_G100 /children/on/children/domain: n=110/110 multiset equal=True order equal=False other differing keys=[]
A_DEFLECTION_G100 /children/on/children/value/unpaired_metadata/applicability/domain_comparison: n=1/1 multiset equal=True order equal=True other differing keys=[]
A_DEFLECTION_G100 /children/outside/children/domain: n=110/110 multiset equal=True order equal=False other differing keys=[]
A_RADAR_LOG_G020 /children/on/children/domain: n=110/110 multiset equal=True order equal=False other differing keys=[]
A_RADAR_LOG_G020 /children/on/children/value/unpaired_metadata/applicability/domain_comparison: n=1/1 multiset equal=True order equal=True other differing keys=[]
A_RADAR_LOG_G020 /children/outside/children/domain: n=110/110 multiset equal=True order equal=False other differing keys=[]
A_ROUND_TRIP_G020 /children/on/children/domain: n=110/110 multiset equal=True order equal=False other differing keys=[]
A_ROUND_TRIP_G020 /children/on/children/value/unpaired_metadata/applicability/domain_comparison: n=1/1 multiset equal=True order equal=True other differing keys=[]
A_ROUND_TRIP_G020 /children/outside/children/domain: n=110/110 multiset equal=True order equal=False other differing keys=[]
F_FLOW_A_ONE_WAY_ER_G020 /comparison/children/on/children/domain: n=97/97 multiset equal=True order equal=False other differing keys=[]
F_FLOW_A_ONE_WAY_ER_G020 /comparison/children/on/children/value/unpaired_metadata/applicability/domain_comparison: n=1/1 multiset equal=True order equal=True other differing keys=[]
F_FLOW_A_ONE_WAY_ER_G020 /comparison/children/outside/children/domain: n=97/97 multiset equal=True order equal=False other differing keys=[]
F_FLOW_A_ONE_WAY_ER_G020 /unpaired_metadata/premise/children/Flux/domain_comparison: n=1/1 multiset equal=True order equal=True other differing keys=[]
F_FLOW_A_ONE_WAY_ER_G020 /unpaired_metadata/premise/children/Measure/domain_comparison: n=2/2 multiset equal=True order equal=True other differing keys=[]
F_FLOW_A_ONE_WAY_ER_G020 /unpaired_metadata/premise/children/NormalExchange/domain_comparison: n=1/1 multiset equal=True order equal=True other differing keys=[]
```

## Grok: the right operand under another profile family
```
$ grep -h '^OK\|^FAILED\|^Ran ' _scratch/s9b_build/s9b_comparator_build_review_r4_grok_evidence/tests_form_ablation.stdout _scratch/s9b_build/s9b_comparator_build_review_r4_grok_evidence/tests_unpatched.stdout
Ran 65 tests in 124.028s
OK
Ran 65 tests in 124.556s
OK
```

The round-3 disposition of the same check:
```
$ grep -o 'evaluating the Wolfram operand under a different profile family[^|]*| \*\*NO REPAIR\.\*\*[^|]*' research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r3_review_disposition.md
evaluating the Wolfram operand under a different profile family passes the tests but moves the printed output on every profile-dependent row. | **NO REPAIR.** Not a required item-6 mutation, and the printed output catches it. 
```

## Record caveats, in the round-4 builder report and full-02
```
$ sed -n '19p' _scratch/s9b_build/comparator_r4/S9b_comparator_build_report.md
Rows with formed comparisons, by kind (overlapping): value 253, domain 22, marker 271, label 22. No formed value: 77 rows; no comparison of any kind: 37. Markers, labels and domains are not counted as compared values.
```

```
$ wc -l < _scratch/s9b_build/comparator_output/repair-04/full-02/rows_without_value.txt; grep -c 'B_GAMMA' _scratch/s9b_build/comparator_output/repair-04/full-02/rows_without_value.txt
77
18
```

```
$ grep -c 'B_GAMMA_[A-Z]*. | marker |' _scratch/s9b_build/comparator_r4/S9b_comparator_build_report_coverage_reasons.md
18
```

```
$ grep -o '| not_formed | container/scalar representation differs | [0-9]* |' _scratch/s9b_build/comparator_r4/S9b_comparator_build_report_coverage_reasons.md
| not_formed | container/scalar representation differs | 18 |
```

## The production output: the reviewed full-02 run
```
$ cat _scratch/s9b_build/comparator_output/repair-04/full-02/guard/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 11116.410525696992,
  "unit": "s11c-guard-c5425124acaf",
  "stderrBytes": 1732566,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 11116.238143314957,
    "childPid": 435791
  }
}
```

```
$ sha256sum _scratch/s9b_build/comparator_output/repair-04/full-02/source/S9b_cross_engine_comparator.py _scratch/s9b_build/comparator_output/repair-04/full-02/source/test_S9b_cross_engine_comparator.py
326d204bcf76f81d40c141803106e42aa8bbb71da91414aef07a3cac710bb1f7  _scratch/s9b_build/comparator_output/repair-04/full-02/source/S9b_cross_engine_comparator.py
efe7249bc18698fedb9d365aa1816a17014ab5967737af680744a5e31c989ad8  _scratch/s9b_build/comparator_output/repair-04/full-02/source/test_S9b_cross_engine_comparator.py
```

```
$ cmp _scratch/s9b_build/comparator_output/repair-04/full-02/source/S9b_cross_engine_comparator.py _scratch/s9b_build/comparator_r4/S9b_cross_engine_comparator.py && echo FULL02_SOURCE_IDENTICAL_TO_REVIEWED_COMPARATOR
FULL02_SOURCE_IDENTICAL_TO_REVIEWED_COMPARATOR
```

Input transcripts named in the configuration record, and the committed transcripts:
```
$ head -1 _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl | grep -o '"input_sha256":{[^}]*}'
"input_sha256":{"py":"b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142","wl":"6519d10c86c9b16741bd32814a57216116c33717a0aec4150f61903617e9ffdb"}
```

```
$ head -1 _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl | grep -o '"target":"[^"]*"'
"target":"no residual target was supplied"
```

```
$ sha256sum research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142  research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out
6519d10c86c9b16741bd32814a57216116c33717a0aec4150f61903617e9ffdb  research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
```

```
$ wc -c < _scratch/s9b_build/comparator_output/repair-04/full-02/comparison.jsonl; wc -l < _scratch/s9b_build/comparator_output/repair-04/full-02/accounting.jsonl
1348220210
423
```

