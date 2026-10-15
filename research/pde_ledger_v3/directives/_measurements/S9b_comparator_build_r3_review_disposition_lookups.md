# Measurements — S9b comparator build review, round 3 (generated 2026-10-15 15:00)

Generator: `_scratch/s9b_build/gen/s9b_comparator_r3_lookups.sh`. Leg evidence: `_scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/` (copied from
`/tmp/s9b_cmp_review_r3_claudeleg/`) and `_scratch/s9b_build/s9b_comparator_build_review_r3_grok_evidence/` (copied from `/tmp/s9b_r3_*`). Operand values the legs computed
stay in that evidence. Line numbers are those of the reviewed comparator `_scratch/s9b_build/comparator_r3/S9b_cross_engine_comparator.py`. Earlier versions: round 0
`3ec82aa9`, round 1 `1af2209c`, round 2 `75676d01`.

## The reviewed artifacts
```
$ sha256sum -c _scratch/s9b_build/s9b_comparator_build_review_baseline_r3.sha256
_scratch/s9b_build/comparator_r3/S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_r3/test_S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_r3/S9b_comparator_build_report.md: OK
_scratch/s9b_build/comparator_r3/S9b_comparator_build_report_function_name_table.md: OK
_scratch/s9b_build/comparator_r3/S9b_comparator_build_report_coverage_reasons.md: OK
_scratch/s9b_build/comparator_r3/probe.py: OK
_scratch/s9b_build/comparator_r3/summarize.py: OK
_scratch/s9b_build/comparator_r3/write_report.py: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/source/S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/source/test_S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/accounting.jsonl: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/output_summary.json: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/coverage_reason_inventory.json: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/rows_without_value.txt: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/validation.stdout: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/verification.json: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/guard/outcome.json: OK
_scratch/s9b_build/comparator_output/repair-03/tests-01/guard/outcome.json: OK
_scratch/s9b_build/comparator_output/repair-03/tests-01/guard/stderr: OK
research/pde_ledger_v3/directives/S9b_comparator_build_directive.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r2_review_disposition.md: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out: OK
research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out: OK
_scratch/s9b_build/comparator_output/repair-03/full-01/comparison.jsonl: OK
```

## Both verdicts
```
$ head -1 _scratch/s9b_build/s9b_comparator_build_review_r3_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_comparator_build_review_r3_grok.txt
Verdict: NEEDS REVISION
```

## C1 (Claude 1): the brane density has one functional form in every generic test profile
```
$ sed -n '758p;762,764p;771p' _scratch/s9b_build/comparator_r3/S9b_cross_engine_comparator.py
    if family=='exponential':
    else:
        d,v,xi,f=1/(101*x),c/(17*sp.sqrt(x)),2*sp.sqrt(x)/19,1/(23*x)
    rho=1+1/(13*x)
    if family=='brane_zero_at_point': rho=(x-anchor)**2/(1+(x-anchor)**2)
```

The same line in each earlier version:
```
$ for c in 3ec82aa9 1af2209c 75676d01; do echo $c; git show $c:research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py | grep -n 'rho=1+1/(13\*x)\|brane_zero_at_point.: rho='; done
3ec82aa9
661:    rho=1+1/(13*x)
668:    if family=='brane_zero_at_point': rho=(x-anchor)**2/(1+(x-anchor)**2)
1af2209c
716:    rho=1+1/(13*x)
723:    if family=='brane_zero_at_point': rho=(x-anchor)**2/(1+(x-anchor)**2)
75676d01
764:    rho=1+1/(13*x)
771:    if family=='brane_zero_at_point': rho=(x-anchor)**2/(1+(x-anchor)**2)
```

```
$ grep -n 'several test profiles of different' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
60:   - For an object that is a functional of a live profile, evaluate under several test profiles of different
```

```
$ git show 34a994d5:research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md | sed -n '56,57p;80p'
  S9's basis writes it as `μ_R` (`steps/S9_light_requires_shear.md:78–79`). Both `μ_⊥(x)` and `ρ_br(x)`
  remain live. The same `ρ_br(x)` appears in the mass balance below. There is one isotropic speed, the same
  `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is
```

```
$ grep -c 'formed_by_profile' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/out/rho_value_leaves.stdout
56
```

```
$ tail -1 _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/out/rho_value_leaves.stdout
TOTAL {'leaves': 56}
```

Value leaves among them, and value leaves with no formed brane_zero_at_point sample (match counts):
```
$ grep -c ' on/value symbolic ' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/out/rho_value_leaves.stdout
50
```

```
$ grep ' on/value symbolic ' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/out/rho_value_leaves.stdout | grep -vc 'brane_zero_at_point'
35
```

```
$ cat _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/ablate/rho_manifest.stdout
WL base is the verbatim production line: True
PY copy is the verbatim production line: True
base count of Derivative[1][rhoBr][integrationRadius] = 2 | corrupt count = 0 | corrupt count of ((1 - rhoBr[integrationRadius])/integrationRadius) = 2
corrupt == base with that literal replaced: True
```

```
$ grep -h '^LEAF' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/ablate/R4_rho_base.stdout _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/ablate/R5_rho_corrupt.stdout
LEAF rhobase F_FLOW_A_ONE_WAY_ER_G020 on/value value symbolic formed=True zero_class=unresolved samples=37 nonzero_or_true_diff=0 []
LEAF rhocorrupt F_FLOW_A_ONE_WAY_ER_G020 on/value value symbolic formed=True zero_class=unresolved samples=37 nonzero_or_true_diff=0 []
```

```
$ grep -h 'Ran \|^OK' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/T1_tests.stdout _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/ablate/T1_tests.stdout 2>/dev/null | head -2
Ran 60 tests in 103.718s
OK
```

## C2 (Claude 2): whether a quadrature sample forms depends on the integral's magnitude
```
$ sed -n '1196,1198p' _scratch/s9b_build/comparator_r3/S9b_cross_engine_comparator.py
        tolerance=mp.power(10,-digits+4)*abs(second)
        if not quadrature_consistent(first,second,error,error2,tolerance):
            raise UnresolvedEvaluation('quadrature error/precision/panel consistency check did not converge')
```

The same tolerance line in each earlier version:
```
$ for c in 3ec82aa9 1af2209c 75676d01; do echo $c; git show $c:research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py | grep -n 'tolerance=mp.power(10,-digits+4)\*abs(second)'; done
3ec82aa9
1af2209c
997:        tolerance=mp.power(10,-digits+4)*abs(second)
75676d01
1093:        tolerance=mp.power(10,-digits+4)*abs(second)
```

Lines (one per row record) carrying the reason, and occurrences:
```
$ grep -c 'quadrature error/precision/panel consistency check did not converge' _scratch/s9b_build/comparator_output/repair-03/full-01/comparison.jsonl
3
```

```
$ grep -o 'quadrature error/precision/panel consistency check did not converge' _scratch/s9b_build/comparator_output/repair-03/full-01/comparison.jsonl | wc -l
12
```

```
$ grep '^FILE' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/out/g111_across_rounds.stdout
FILE repair-01/full-03/comparison.jsonl
FILE repair-02/full-01/comparison.jsonl
FILE repair-03/full-01/comparison.jsonl
```

```
$ grep -c 'G111 children/on/children/value/comparison p[01]/exponential | cov: unresolved coverage | PY: quadrature UnresolvedEvaluation' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/out/g111_across_rounds.stdout
18
```

```
$ grep -o 'consistent=[A-Za-z]*' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/ablate/g111_diag.stdout | sort | uniq -c
      4 consistent=False
      8 consistent=True
```

## G1 (Grok 1): density-part labels are compared through an uncited label map
```
$ sed -n '733p;1440,1442p' _scratch/s9b_build/comparator_r3/S9b_cross_engine_comparator.py
LABEL_TABLE={'RHO_POSITIVE':'rho > 0','RHO_NONPOSITIVE':'rho <= 0'}
def text_comparison(left,right,role):
    lv=LABEL_TABLE.get(left,left) if role=='label' else left
    return lv==right
```

```
$ git show 3ec82aa9:research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py | grep -n 'LABEL_TABLE='
630:LABEL_TABLE={'RHO_POSITIVE':'rho > 0','RHO_NONPOSITIVE':'rho <= 0'}
```

Label entries in the builder's name-table appendix (match count):
```
$ grep -ci 'RHO_POSITIVE\|rho > 0\|label' _scratch/s9b_build/comparator_r3/S9b_comparator_build_report_function_name_table.md
0
```

```
$ sed -n '40,46p' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md

2. **Symbols and functions through an explicit table.** Any map between the two engines' symbol and function names
   (the coordinates, the impact parameter and endpoints, `GM`, `c₀`, the live profiles, the OPEN quantities, the
   Part C parameters) is an explicit, injective table in the comparator's source. Each entry cites the spec object
   both names denote and the line in each engine that introduces it. A mapped name is a substitution, never by
   itself a zero. A **repoint control** must show that pointing a mapped name at a different object moves the
   printed residual.
```

```
$ sed -n '48,52p' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
3. **Value and domain compared separately; labelled parts compared part by part.** Where a payload attaches a
   domain to its computed object (item 14 of the build directive defines both), the comparator prints a value
   residual and a domain comparison separately. Where a Part C row carries the labelled `ρ > 0` and `ρ ≤ 0` parts
   (item 9), each part is compared with its counterpart, and the labels themselves are compared. A
   `NOT_ESTABLISHED` marker is compared as a marker together with the domain it carries.
```

```
$ sed -n '72p;94p' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
   - A difference in representation is printed as a difference. It is never reconciled inside the comparator.
   - the two density labels swapped;
```

```
$ sed -n '126,127p' _scratch/s9b_build/comparator_r3/test_S9b_cross_engine_comparator.py
    def test_density_labels_swapped(self):
        self.mutation(py_parts('Integer(7)'),wl_parts('7'),wl_parts('7',labels=('rho <= 0','rho > 0')),channel='text_comparison')
```

```
$ grep -A3 '^label comparisons' _scratch/s9b_build/s9b_comparator_build_review_r3_grok_evidence/s9b_r3_fidelity_census.stdout
label comparisons 44
  equal strings 0 unequal strings 44 unequal but comparison true 44
LABEL_QUOTE {"name": "C_FIXED_RATIO_CONDITION_DEFLECTION", "comparison": true, "py_text": "RHO_NONPOSITIVE", "wl_text": "rho <= 0", "comparison_formed": true, "py_path": [2, 0], "wl_path": [1, 0, 1, 0, 0], "role": "label"}
LABEL_PAIRS Counter({('RHO_NONPOSITIVE', 'rho <= 0', True): 22, ('RHO_POSITIVE', 'rho > 0', True): 22})
```

```
$ grep '^LABEL_PAIRS' _scratch/s9b_build/s9b_comparator_build_review_r3_grok_evidence/s9b_r3_fidelity_census.stdout
LABEL_PAIRS Counter({('RHO_NONPOSITIVE', 'rho <= 0', True): 22, ('RHO_POSITIVE', 'rho > 0', True): 22})
```

## O1 (both legs): the leaf cache key has no control
```
$ sed -n '1514p' _scratch/s9b_build/comparator_r3/S9b_cross_engine_comparator.py
        key=(lv,rv,role)
```

```
$ grep -h 'ABLATED_TESTS_EXIT\|^Ran \|^OK' _scratch/s9b_build/s9b_comparator_build_review_r3_grok_evidence/s9b_r3_cache_ablation.stdout
ABLATED_TESTS_EXIT 0
Ran 60 tests in 102.333s
OK
```

```
$ grep -h 'Ran \|^OK' _scratch/s9b_build/s9b_comparator_build_review_r3_claude_evidence/ablate/T3_cache_tests.stdout
Ran 60 tests in 102.593s
OK
```

