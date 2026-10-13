# Measurements — S9b scoped build review, round 6 (generated 2026-10-12 22:16)

Generator: `_scratch/s9b_build/gen/s9b_build_r6_lookups.sh` (sha256sum/grep/sed/cat/wc/ls/awk only). `_scratch/s9b_build/s9b_repair_build_review_r6_claude_evidence/`
holds the Claude leg's scripts and stdout, copied from `/tmp/s9b_rev_r6_fc_1791856929/` (files under 2 MB; the
larger ones are listed in its `NOT_COPIED_over_2MB.txt`); `_scratch/s9b_build/s9b_repair_build_review_r6_grok_evidence/` holds Grok's, copied from `/tmp/s9b_r6_*`.

## The reviewed artifacts are the working-tree files, unchanged after both legs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r6.sha256
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/scripts/S9b_exports.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_not_evaluated_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_k11_evaluated_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_not_established_tags.txt: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_executable_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_payload_checks.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_validation.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_undecided_rows.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_changed_rows.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_difference_rows.json: OK
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round6_measurements.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: FAILED
sha256sum: WARNING: 1 computed checksum did NOT match
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_review_prompt_r6.md
7be1fe7ff347b3c25ca54d57a8655e31bf4ae0c354d078c129adad13ad411a82  _scratch/s9b_build/s9b_repair_build_review_prompt_r6.md
```

```
$ cat _scratch/s9b_build/r6_harness_output/py/final3/source-pins.json
{
  "research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py": "a230d052203aee81f08f942b1be3f8ed08e501776a862a179f8bcef91af94308",
  "research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py": "06818e19a8d29e36a7e9b8c86c6032d5d22d17f483ec503d4050718483f78dc6",
  "_scratch/s9b-repair-amend9/final2/method_probe.py": "2eafe0b8941b6a170260f4d7e3c5ea97222ca9c2c9deca1e393be30db986f144"
}
```

```
$ grep -o '"research[^"]*ablation.wl": "[0-9a-f]*"\|"research[^"]*audit.wl": "[0-9a-f]*"' _scratch/s9b_build/r6_harness_output/wl/builder_scratch/final_source_hashes.json
"research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl": "2028badd3772d68d43659e06d6ebc6644843dd0fb77563fc6fff37255f241561"
"research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl": "373219de8394ec0850630cb3df43d58f2c05334e3b6c3292111aead5ae5f003b"
```

```
$ cat _scratch/s9b_build/r6_harness_output/py/harness/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 10659.089266449,
  "unit": "s11c-guard-f4651997369b",
  "stderrBytes": 2082357,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 10658.951396039,
    "childPid": 6728
  }
}
```

## Both verdicts
```
$ grep -n -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r6_claude.md
1:Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r6_grok.txt
Verdict: NEEDS REVISION
```

```
$ cat _scratch/s9b_build/s9b_repair_build_review_r6_grok.DONE
EXIT=0
```

## The finding is round 5's C1, again
```
$ grep -n -o '| C1 | Both harnesses print UNDECIDED for differences the engine can decide.\{0,60\}' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r5_review_disposition.md
55:| C1 | Both harnesses print UNDECIDED for differences the engine can decide. SymPy prints a change only when `residual.is_zero is False
```

## Directive: item 14 asks for a point under some admissible profile; item 13 forbids choosing a profile
```
$ grep -n 'such a point under some admissible profile' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
287:      such a point under some admissible profile. An admissible profile is a general radial function (item 13)
```

```
$ grep -n 'Where the engine cannot decide whether a part is unchanged or' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ grep -n 'No choice may turn a live profile into a chosen one' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
256:    of these. No choice may turn a live profile into a chosen one or add an expected value.
```

```
$ grep -n 'No profile family is used for a shared-vocabulary object' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
267:    - No profile family is used for a shared-vocabulary object. An engine may also print family-restricted
```

## The harnesses' change tests
### Wolfram: the round-6 comment rules out samples and profile substitutions
```
$ sed -n '122,128p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
(* Amendment 9: equality is universal, change is existential. An unresolved
   pointwise equality is NOT evidence that no difference can be established.
   Exact real quantifier elimination proves existence on the common defined
   domain. Unsupported functional/analytic formulas retain their COMPLETE
   existential obligation and UNDECIDED status; no samples, term cutoffs,
   profile substitutions, timed fallbacks, or independent-jet assumptions.
   Real noninteger powers are tested on their nonnegative-base chart, which
```

```
$ grep -n 'realFormulaQ\[\|CompleteRealDecision' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | head -12
137:realFormulaQ[z_] := Which[
157:  supported = realFormulaQ[condition] &&
177:    "SupportedFormula" -> supported,"CompleteRealDecision" -> complete|>;
190:  If[TrueQ[change["CompleteRealDecision"]] && SameQ[change["ExistenceResult"],False],Return[0]];
```

```
$ grep -n '"Status" -> If' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
193:  <|"Status" -> If[TrueQ[change["ExistenceResult"]],"CHANGED","UNDECIDED"],
```

### SymPy: the change words it can print, and its existence test under an attached domain
```
$ grep -o "Str('[A-Z_]*')" research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | sort | uniq -c | grep -E 'RATIONAL_NONIDENTITY|EXISTENTIAL_DIFFERENCE|SET_DIFFERENCE_WITNESS|UNDECIDED'
      2 Str('EXISTENTIAL_DIFFERENCE')
      1 Str('RATIONAL_NONIDENTITY')
      1 Str('SET_DIFFERENCE_WITNESS')
      1 Str('UNDECIDED')
```

```
$ grep -n '_linear_existence' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
335:def _linear_existence(predicate):
399:        certificate=_linear_existence(sp.And(context,result))
427:            certificate=_linear_existence(sp.And(context,sp.Ne(result,0)))
```

```
$ grep -n 'if whole.has(sp.Integral,sp.Limit,sp.Subs,sp.Set)' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
289:    if whole.has(sp.Integral,sp.Limit,sp.Subs,sp.Set):
```

```
$ grep -n 'RADIAL_JET_PROFILES = ' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
283:RADIAL_JET_PROFILES = frozenset(('s9b_delta','s9b_V','s9b_xi',
```

## What the harnesses printed, counted per knife
### Wolfram (`_scratch/s9b_build/r6_harness_output/wl/harness/stdout`): Status words; and UNDECIDED on an inert integral of zero
```
$ for k in K1 K2 K3 K4a K4b K4c K5 K6 K7 K8 K9a K9b K10 K11; do printf '%s CHANGED=%s UNDECIDED=%s UNDECIDED_on_Inactive_Integrate_0=%s\n' $k $(grep "^WL_LOCAL_S9B_ABLATION_${k}_" _scratch/s9b_build/r6_harness_output/wl/harness/stdout | grep -o '"Status" -> "CHANGED"' | wc -l) $(grep "^WL_LOCAL_S9B_ABLATION_${k}_" _scratch/s9b_build/r6_harness_output/wl/harness/stdout | grep -o '"Status" -> "UNDECIDED"' | wc -l) $(grep "^WL_LOCAL_S9B_ABLATION_${k}_" _scratch/s9b_build/r6_harness_output/wl/harness/stdout | grep -o '"Status" -> "UNDECIDED", "Equality" -> Inactive\[Integrate\]\[0, {r, b, Infinity}\] == 0' | wc -l); done
K1 CHANGED=23 UNDECIDED=1924 UNDECIDED_on_Inactive_Integrate_0=16
K2 CHANGED=10 UNDECIDED=1222 UNDECIDED_on_Inactive_Integrate_0=18
K3 CHANGED=11 UNDECIDED=1507 UNDECIDED_on_Inactive_Integrate_0=18
K4a CHANGED=0 UNDECIDED=734 UNDECIDED_on_Inactive_Integrate_0=6
K4b CHANGED=126 UNDECIDED=1049 UNDECIDED_on_Inactive_Integrate_0=0
K4c CHANGED=0 UNDECIDED=1223 UNDECIDED_on_Inactive_Integrate_0=18
K5 CHANGED=0 UNDECIDED=44 UNDECIDED_on_Inactive_Integrate_0=0
K6 CHANGED=0 UNDECIDED=432 UNDECIDED_on_Inactive_Integrate_0=0
K7 CHANGED=36 UNDECIDED=1502 UNDECIDED_on_Inactive_Integrate_0=0
K8 CHANGED=0 UNDECIDED=770 UNDECIDED_on_Inactive_Integrate_0=0
K9a CHANGED=2 UNDECIDED=336 UNDECIDED_on_Inactive_Integrate_0=0
K9b CHANGED=0 UNDECIDED=325 UNDECIDED_on_Inactive_Integrate_0=0
K10 CHANGED=0 UNDECIDED=792 UNDECIDED_on_Inactive_Integrate_0=0
K11 CHANGED=11 UNDECIDED=6 UNDECIDED_on_Inactive_Integrate_0=0
```

### SymPy (`_scratch/s9b_build/r6_harness_output/py/harness/stdout`, 9.7 GB, one pass): change words and UNDECIDED
```
$ LC_ALL=C grep -o -E "^PY_LOCAL_S9B_ABLATION_K[0-9]+[a-c]?_|Str\('(RATIONAL_NONIDENTITY|EXISTENTIAL_DIFFERENCE|SET_DIFFERENCE_WITNESS|UNDECIDED)'\)" _scratch/s9b_build/r6_harness_output/py/harness/stdout | awk '/^PY_LOCAL/{split($0,a,"_"); k=a[5]; next} {c[k" "$0]++} END{for (x in c) print x, c[x]}' | sort -V
K1 Str('RATIONAL_NONIDENTITY') 1
K1 Str('UNDECIDED') 1736
K2 Str('RATIONAL_NONIDENTITY') 1
K2 Str('UNDECIDED') 972
K3 Str('RATIONAL_NONIDENTITY') 1
K3 Str('UNDECIDED') 1528
K4a Str('UNDECIDED') 454
K4b Str('RATIONAL_NONIDENTITY') 4
K4b Str('UNDECIDED') 1272
K4c Str('UNDECIDED') 972
K5 Str('UNDECIDED') 33
K6 Str('RATIONAL_NONIDENTITY') 5
K6 Str('UNDECIDED') 349
K7 Str('UNDECIDED') 506
K8 Str('UNDECIDED') 690
K9a Str('UNDECIDED') 167
K9b Str('UNDECIDED') 434
K10 Str('UNDECIDED') 1145
```

## The legs' evidence
```
$ cat _scratch/s9b_build/s9b_repair_build_review_r6_claude_evidence/count/wl_changed_by_knife.stdout
K1 tags_with_CHANGED=14 tags_with_UNDECIDED=314
K2 tags_with_CHANGED=4 tags_with_UNDECIDED=234
K3 tags_with_CHANGED=4 tags_with_UNDECIDED=334
K4a tags_with_CHANGED=0 tags_with_UNDECIDED=104
K4b tags_with_CHANGED=20 tags_with_UNDECIDED=104
K4c tags_with_CHANGED=0 tags_with_UNDECIDED=235
K5 tags_with_CHANGED=0 tags_with_UNDECIDED=38
K6 tags_with_CHANGED=0 tags_with_UNDECIDED=54
K7 tags_with_CHANGED=26 tags_with_UNDECIDED=47
K8 tags_with_CHANGED=0 tags_with_UNDECIDED=184
K9a tags_with_CHANGED=1 tags_with_UNDECIDED=13
K9b tags_with_CHANGED=0 tags_with_UNDECIDED=12
K10 tags_with_CHANGED=0 tags_with_UNDECIDED=325
K11 tags_with_CHANGED=4 tags_with_UNDECIDED=2
```

```
$ cat _scratch/s9b_build/s9b_repair_build_review_r6_claude_evidence/count/wl_inert_zero.stdout | head -20
UNDECIDED leaves with Equality Inactive[Integrate][0,...] == 0:
112
tags (lines) containing such a leaf, by knife:
     25 K1
     24 K2
     24 K3
      6 K4a
     24 K4c
      9 K5
SymPy harness: leaves UNDECIDED whose BASELINE or CORRUPTED is Integral(Integer(0), ...):
0
```

```
$ cat _scratch/s9b_build/s9b_repair_build_review_r6_claude_evidence/count/py_printed_changes.stdout | head -40
K1	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K1	PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G020	value=OTHER	domain=ZERO
K1	PY_LOCAL_S9B_VECTOR_DISPERSION	value=OTHER+CERT	domain=ZERO
K2	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K2	PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G001	value=OTHER	domain=ZERO
K2	PY_LOCAL_S9B_VECTOR_DISPERSION	value=OTHER+CERT	domain=ZERO
K3	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K3	PY_LOCAL_S9B_VECTOR_DISPERSION	value=OTHER+CERT	domain=ZERO
K4a	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K4b	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K4b	PY_LOCAL_S9B_OPTICAL_Q	value=OTHER+CERT	domain=ZERO
K4c	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K4c	PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G001	value=OTHER	domain=ZERO
K5	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K5	PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G010	value=OTHER	domain=ZERO
K6	PY_LOCAL_S9B_ABEL_REFERENCE_INVERSES	value=OTHER+CERT	domain=ZERO
K6	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K6	PY_LOCAL_S9B_RADAR_ANGULAR_OPERATOR_RATIO	value=OTHER	domain=ZERO
K6	PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G001	value=OTHER	domain=ZERO
K6	PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G010	value=OTHER	domain=ZERO
K6	PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G020	value=OTHER	domain=ZERO
K7	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K8	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K9a	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K9b	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
K10	PY_LOCAL_S9B_BUILD_INPUT_DIGESTS	value=OTHER	domain=ZERO
```

```
$ grep -n 'wl_form\|sympy_form' _scratch/s9b_build/s9b_repair_build_review_r6_grok.txt | head -5 | cut -c1-200
28:The same gate is what the form ablation moves. Each copy keeps the harness and replaces applied functions and their derivatives by an admissible chart. A residual rescale is not the mutation. Liter
```

## The two legs differ on the Wolfram inert-zero case
```
$ grep -n 'INERT_ZERO' _scratch/s9b_build/s9b_repair_build_review_r6_grok_evidence/s9b_r6_form_1791856918/wl_unmutated.stdout _scratch/s9b_build/s9b_repair_build_review_r6_grok_evidence/s9b_r6_form_1791856918/sympy_unmutated.stdout | head -4
_scratch/s9b_build/s9b_repair_build_review_r6_grok_evidence/s9b_r6_form_1791856918/wl_unmutated.stdout:4:INERT_ZERO: <|"ValueDifference" -> 0, "DomainDifference" -> 0|>
_scratch/s9b_build/s9b_repair_build_review_r6_grok_evidence/s9b_r6_form_1791856918/sympy_unmutated.stdout:4:INERT_ZERO: Tuple(Tuple(Str('VALUE_DIFFERENCE'), Integer(0)), Tuple(Str('DOMAIN_DIFFERENCE'), Integer(0)))
```

```
$ grep -n 'INERT_ZERO' _scratch/s9b_build/s9b_repair_build_review_r6_grok_evidence/s9b_r6_form_1791856918/drive_wl_unmutated.wl | head -4
9:row["INERT_ZERO", Inactive[Integrate][0, {r, 0, 1}], Inactive[Integrate][0, {r, 0, 1}]];
```

## Guard state after both legs
```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

```
$ systemctl --user list-units 's11c-guard-*' --no-legend | wc -l
1
```

