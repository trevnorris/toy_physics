# Measurements — S9b scoped build review, round 5 (generated 2026-10-12 07:47)

Generator: `_scratch/s9b_build/gen/s9b_build_r5_lookups.sh` (sha256sum/grep/sed/cat/wc/ls only). `_scratch/s9b_build/s9b_repair_build_review_r5_claude_evidence/` holds
the Claude leg's scripts and stdout, copied from `/tmp/s9b_rev5_leg/`; `_scratch/s9b_build/s9b_repair_build_review_r5_grok_evidence/` holds Grok's, copied from
`/tmp/s9b_r5_grok_review_20261012/` (files under 2 MB).

## The reviewed artifacts are the working-tree files, unchanged after both legs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r5.sha256
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py: OK
research/pde_ledger_v3/scripts/S9b_exports.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
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
research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round5_measurements.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_review_prompt_r5.md
51c70dbf2bdfbf07a98544440f70eb5b5759bdfaef21324a9837e2cbde733723  _scratch/s9b_build/s9b_repair_build_review_prompt_r5.md
```

```
$ cat _scratch/s9b_build/r5_harness_output/py/builder_state/source-pins.json
{
  "research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py": "a230d052203aee81f08f942b1be3f8ed08e501776a862a179f8bcef91af94308",
  "research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py": "6f94e6e34da4e3a1c3462d51dfbea9a5650765d1bcba0f9449d6c1fc2a2e7cc7"
}
```

```
$ grep -o '"research[^"]*ablation.wl": "[0-9a-f]*"\|"research[^"]*audit.wl": "[0-9a-f]*"' _scratch/s9b_build/r5_harness_output/wl/builder_scratch/final_source_hashes.json
"research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl": "2028badd3772d68d43659e06d6ebc6644843dd0fb77563fc6fff37255f241561"
"research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl": "2ba2b150e969ee673203a4a9eaf5e45d043862faf1d75a0815ac216295d759ed"
```

## Both verdicts
```
$ grep -n -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r5_claude.md
2:Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r5_grok.txt
Verdict: CLEAR
```

```
$ cat _scratch/s9b_build/build_r5_grok.DONE
EXIT=0
```

## Finding C1 (Claude 1): a difference the engine could decide prints as UNDECIDED
### Directive item 14
```
$ grep -n 'Where the engine cannot decide whether a part is' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
270:      including a condition carried inside a set or a quantifier. Where the engine cannot decide whether a part is
```

```
$ grep -n 'so an unchanged part prints as zero' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
267:      condition. Each part's difference is computed by the engine, so an unchanged part prints as zero. A part is
```

### The harnesses' change criteria
```
$ grep -n 'is_zero is False\|is_zero' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
127:        if residual.is_zero is True:return sp.S.Zero
128:        if residual.is_zero is False:return residual
264:                        if factor.is_zero is False:divisor*=factor
300:        if result.is_zero is True:return sp.S.Zero
301:        if result.is_zero is False:return result
```

```
$ grep -n 'NEVER proves a change' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
68:   A nonzero structural/algebraic residual alone NEVER proves a change. *)
```

```
$ grep -n '"UNDECIDED"\|"CHANGED"' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
131:  <|"Status" -> If[SameQ[proof,False],"CHANGED","UNDECIDED"],
```

### The round-4 harness printed the residual
```
$ git show 6b6e2ce4:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | grep -n 'return b-a'
249:        return b-a
251:        return b-a
```

### Quoted lines in the builders' harness output
```
$ LC_ALL=C grep -m1 -F 'PY_LOCAL_S9B_ABLATION_K6_PY_LOCAL_S9B_ABEL_REFERENCE_INVERSES' _scratch/s9b_build/r5_harness_output/py/k1-k10/stdout | grep -o 'UNDECIDED.\{0,300\}' | head -2
UNDECIDED'), Tuple(Str('BASELINE'), Mul(Integer(2), Symbol('GM', real=True), Pow(Symbol('c_0', positive=True), Integer(-2)), Pow(Symbol('s9b_r', positive=True), Integer(-1)))), Tuple(Str('CORRUPTED'), Mul(Integer(4), Symbol('GM', real=True), Pow(Symbol('c_0', positive=True), Integer(-2)), Pow(Symbol('s9b_r',
```

```
$ LC_ALL=C grep -m1 -F 'PY_LOCAL_S9B_ABLATION_K9b_PY_LOCAL_S9B_C_POWER_RESPONSE' _scratch/s9b_build/r5_harness_output/py/k1-k10/stdout | grep -o 'UNDECIDED.\{0,200\}' | head -2
UNDECIDED'), Tuple(Str('BASELINE'), Pow(Add(Symbol('s9b_fractional_bulk_change', real=True), Integer(1)), Symbol('s', real=True))), Tuple(Str('CORRUPTED'), Integer(1)), Tuple(Str('RESIDUAL'), Add(Integer(1), M
UNDECIDED'), Tuple(Str('BASELINE'), Mul(Symbol('s', real=True), Function('s9b_f', **{'complex': True, 'extended_real': True, 'finite': True, 'hermitian': True, 'imaginary': False, 'infinite': False, 'real': Tr
```

```
$ LC_ALL=C grep -m1 -F 'WL_LOCAL_S9B_ABLATION_K5_WL_S9B_F_FLOW_A_NONRECIPROCAL_G011' _scratch/s9b_build/r5_harness_output/wl/harness/stdout | grep -o '"Status" -> "[A-Z]*"' | sort | uniq -c
      1 "Status" -> "UNDECIDED"
```

```
$ LC_ALL=C grep -m1 -F 'WL_LOCAL_S9B_ABLATION_K5_WL_S9B_F_FLOW_A_NONRECIPROCAL_G011' _scratch/s9b_build/r5_harness_output/wl/harness/stdout | grep -o 'Inactive\[Integrate\]\[0, .\{0,60\}' | head -1
Inactive[Integrate][0, {integrationRadius, Sqrt[b^2 + ZE^2], Sqrt[b^2 + ZR^2]}], In
```

### The legs' counts
```
$ cat _scratch/s9b_build/s9b_repair_build_review_r5_claude_evidence/wl_harness_status.json | head -c 1500
[{"knife": "K1", "tag": "WL_LOCAL_S9B_ABEL_REDUCTION", "changed": 0, "undecided": 2, "allzero": false, "value0": false, "difflen": 2104}, {"knife": "K1", "tag": "WL_LOCAL_S9B_BRANCH_FREQUENCIES", "changed": 0, "undecided": 2, "allzero": false, "value0": false, "difflen": 1202}, {"knife": "K1", "tag": "WL_LOCAL_S9B_BRANCH_TYPES", "changed": 0, "undecided": 1, "allzero": false, "value0": false, "dif
```

```
$ cat _scratch/s9b_build/s9b_repair_build_review_r5_claude_evidence/py_decided_literal.stdout | head -30
K1 {'tags': 379, 'with_undecided': 336, 'with_decided_symbolic_change': 1} | decided-symbolic tags: ['PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G020']
K10 {'tags': 379, 'with_undecided': 319} | decided-symbolic tags: []
K2 {'tags': 379, 'with_undecided': 214, 'with_decided_symbolic_change': 1} | decided-symbolic tags: ['PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G001']
K3 {'tags': 379, 'with_undecided': 334} | decided-symbolic tags: []
K4a {'tags': 379, 'with_undecided': 96} | decided-symbolic tags: []
K4b {'tags': 379, 'with_undecided': 273} | decided-symbolic tags: []
K4c {'tags': 379, 'with_undecided': 213, 'with_decided_symbolic_change': 1} | decided-symbolic tags: ['PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G001']
K5 {'tags': 379, 'with_undecided': 31, 'with_decided_symbolic_change': 1} | decided-symbolic tags: ['PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G010']
K6 {'tags': 379, 'with_undecided': 66, 'with_decided_symbolic_change': 4} | decided-symbolic tags: ['PY_LOCAL_S9B_RADAR_ANGULAR_OPERATOR_RATIO', 'PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G001', 'PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G010', 'PY_LOCAL_S9B_RADAR_ENDPOINT_COUNTING_CERTIFICATE_G020']
K7 {'tags': 379, 'with_undecided': 50} | decided-symbolic tags: []
K8 {'tags': 379, 'with_undecided': 167} | decided-symbolic tags: []
K9a {'tags': 379, 'with_undecided': 14} | decided-symbolic tags: []
K9b {'tags': 379, 'with_undecided': 15} | decided-symbolic tags: []
```

```
$ grep -n 'UNDECIDED' _scratch/s9b_build/s9b_repair_build_review_r5_grok_evidence/harness_quotes.txt | head -5 | cut -c1-200
1:K9A_RESPONSE_LEN 18046 VALUE_AT 6316 DOMAIN_AT 8378 UNDECIDED 4 Integer0 48
2:VALUE_WINDOW VALUE_DIFFERENCE'), Tuple(Integer(0), Integer(0), Integer(0), Tuple(Integer(0), Tuple(Tuple(Integer(0), Tuple(Str('UNDECIDED'), Tuple(Str('BASELINE'), Mul(Pow(Symbol('rho_0', positive=T
3:DOMAIN_WINDOW DOMAIN_DIFFERENCE'), Tuple(Tuple(Tuple(Integer(3), Integer(1), Integer(0), Integer(2)), Tuple(Str('UNDECIDED'), Tuple(Str('BASELINE'), And(StrictGreaterThan(Mul(Symbol('K', real=True),
```

```
$ grep -n 'hash counts\|open comparison is' _scratch/s9b_build/s9b_repair_build_review_r5_grok.txt | cut -c1-200
35:**Harness.** Builder pins match the baseline hashes. From the builders’ stdout, an unchanged slot is `Integer(0)`, and an open comparison is `UNDECIDED`. K9a’s fixed-ratio response (`harness_qu
```

## Guard state after both legs
```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

```
$ systemctl --user list-units 's11c-guard-*' --no-legend | wc -l
0
```

```
$ cat _scratch/s11c/s9b-review-r5/legC_py_harness_decided/outcome.json
{
  "exitCode": 1,
  "wallSeconds": 1269.8706513009965,
  "unit": "s11c-guard-f26affde3cf5",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}
```

