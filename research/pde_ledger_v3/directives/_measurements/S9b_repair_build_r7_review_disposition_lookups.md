# Measurements — S9b scoped build review, round 7 (generated 2026-10-16 01:38)

Generator: `_scratch/s9b_build/gen/s9b_build_r7_lookups.sh` (sha256sum/grep/sed/cat/wc/git show, and literal
field retrieval from the named inventory JSON). `_scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/` holds the Claude leg's scripts and stdout, copied from
`/tmp/s9b_review_r7_claudeleg/`; `_scratch/s9b_build/s9b_repair_build_review_r7_grok_evidence/` holds Grok's, copied from `/tmp/s9b_r7_*`. Reviewer witness values stay
in that evidence. Line numbers are those of the reviewed working-tree harnesses.

## The reviewed artifacts are the working-tree files, unchanged after both legs
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r7.sha256 | grep -c ': OK$'
27
```

```
$ wc -l < _scratch/s9b_build/s9b_repair_build_review_baseline_r7.sha256
27
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_review_prompt_r7.md
010b53ba88b58c3aec1b511d3bde4d19238196085ff20e09598d697dfc97c6a1  _scratch/s9b_build/s9b_repair_build_review_prompt_r7.md
```

## Both verdicts
```
$ head -1 _scratch/s9b_build/s9b_repair_build_review_r7_claude.md
## Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r7_grok.txt
Verdict: NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_repair_build_review_r7_claude.md; grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_repair_build_review_r7_grok.txt
3
1
```

## R1 (Claude 1): a value prints zero when the two attached domains share no point
The shortcut in each harness, and at the round-6 build `dfa07223`:
```
$ grep -n 'context is sp.false' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
1020:    if a==b or context is sp.false:return sp.S.Zero
```

```
$ grep -n 'SameQ\[change\["ExistenceResult"\],False\],Return\[0\]' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
328:  If[TrueQ[change["CompleteRealDecision"]] && SameQ[change["ExistenceResult"],False],Return[0]];
```

```
$ git show dfa07223:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | grep -n 'context is sp.false'
116:    if a==b or context is sp.false:
385:    if a == b or context is sp.false:
```

```
$ git show dfa07223:research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | grep -n 'SameQ\[change\["ExistenceResult"\],False\],Return\[0\]'
190:  If[TrueQ[change["CompleteRealDecision"]] && SameQ[change["ExistenceResult"],False],Return[0]];
```

The builder's own Wolfram check: its input, and what the final run printed:
```
$ grep -n 'show\["DISJOINT_DOMAINS"' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation_checks.wl
33:show["DISJOINT_DOMAINS",<|"Value"->x,"Domain"->x>0|>,<|"Value"->2x,"Domain"->x<0|>];
```

```
$ grep -o '^DISJOINT_DOMAINS: <|"ValueDifference" -> [^,]*' _scratch/s9b_build/r7_harness_output/wl/comparison-final2/stdout
DISJOINT_DOMAINS: <|"ValueDifference" -> 0
```

The Claude leg's SymPy FORM ablation (mutA removes the shortcut): diff of its driver output:
```
$ grep -c '' _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/form_py/out/diff_driver_unmut_mutA.txt; grep -o '^[<>] FIXED [A-Z_0-9a-z]* | VALUE: [a-zA-Z]* change=[0-9]* undecided=[0-9]*' _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/form_py/out/diff_driver_unmut_mutA.txt
12
< FIXED DISJOINT_DOMAINS_VALUES_x_vs_2x | VALUE: ZERO change=0 undecided=0
< FIXED COMPLEMENT_DOMAINS_VALUES_x_vs_2x | VALUE: ZERO change=0 undecided=0
< FIXED FALSE_CORRUPTED_DOMAIN_VALUES_1_vs_2 | VALUE: ZERO change=0 undecided=0
> FIXED DISJOINT_DOMAINS_VALUES_x_vs_2x | VALUE: nonzero change=0 undecided=1
> FIXED COMPLEMENT_DOMAINS_VALUES_x_vs_2x | VALUE: nonzero change=0 undecided=1
> FIXED FALSE_CORRUPTED_DOMAIN_VALUES_1_vs_2 | VALUE: nonzero change=0 undecided=1
```

The builder's SymPy checks under the same mutation (lines of diff):
```
$ grep -c '' _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/form_py/out/diff_checks_unmut_mutA.txt
0
```

## R2 (Claude 2): a value witness is tried only where both attached domains evaluate to exactly true
```
$ sed -n '954,957p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
            allowed=_evaluate_at(context,profiles,points)
            allowed,allowed_error,allowed_certificate,reason=_certified_value(allowed)
            if reason:unresolved['COMMON_DOMAIN:'+reason]+=1
            if allowed is not sp.true:continue
```

```
$ sed -n '281,282p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
  cv=exactWitnessEvaluation[prepared[[3]]];
  If[TrueQ[cv],
```

```
$ git show dfa07223:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | grep -c 'allowed is not sp.true'; git show dfa07223:research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | grep -c 'TrueQ\[cv\]'
0
0
```

The Wolfram builder's statement of the rule:
```
$ grep -o 'Value comparisons require both attached domains\.' research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_measurements.md
Value comparisons require both attached domains.
```

The directive's text, item 14:
```
$ grep -n 'for a value, a point of its variables where both parts are defined' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
284:      baseline part at some point: for a value, a point of its variables where both parts are defined; for a domain,
```

```
$ grep -n 'the difference shows separately whether the value changed and whether' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
276:    - For a payload that carries a domain, the difference shows separately whether the value changed and whether
```

The Claude leg's Wolfram FORM ablation (mutB relaxes the gate): site lines that move, and the checks diff:
```
$ grep -o '^[<>] [A-Z]* [A-Za-z0-9_]* .*ValueDifference=.*' _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/form_wl/out/diff_site_unmut_mutB.txt | cut -c1-40 | sort | uniq -c
      1 < FIXED PROSE_IN_DOMAIN_VALUE_1_vs_2 Val
      1 > FIXED PROSE_IN_DOMAIN_VALUE_1_vs_2 Val
      1 < SITE K1_WL_LOCAL_S9B_ABEL_REDUCTION re
      1 > SITE K1_WL_LOCAL_S9B_ABEL_REDUCTION re
      1 < SITE K9a_WL_LOCAL_S9B_C_RESPONSES reco
      1 > SITE K9a_WL_LOCAL_S9B_C_RESPONSES reco
      1 < SITE K9b_WL_LOCAL_S9B_C_RESPONSES reco
      1 > SITE K9b_WL_LOCAL_S9B_C_RESPONSES reco
```

```
$ grep -c '' _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/form_wl/out/diff_checks_unmut_mutB.txt
0
```

SymPy K7, the knife's own relation part, unmutated (reason counts the leg's script printed):
```
$ grep -h 'reasons:\|statuses:' _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/form_py/out/k7_reasons_unmut.stdout
K7 B_CONDITION_DEFLECTION reasons: [('COMMON_DOMAIN:EVALUATION_UNRESOLVED_NATIVE_TREE_BUDGET', '36')]
statuses: ['UNDECIDED']
```

## R3 (Grok 1): Wolfram prints a part undecided after its equality test returned False
```
$ grep -n 'SameQ\[equality,False\],True,True' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
163:    SameQ[equality,False],True,True,aa != bb];
```

Leaves in the final Wolfram run, per knife, printed with this exact sequence:
```
$ P='"Status" -> "UNDECIDED", "Equality" -> False, "ChangeProof" -> <|"ExistenceCondition" -> Inactive\[Exists\]\[{}, True\], "ExistenceResult" -> Missing\["UNDECIDED"\]'; for k in K1 K2 K3 K4a K4b K4c K5 K6 K7 K8 K9a K9b K10 K11; do printf '%s %s\n' $k $(LC_ALL=C grep "^WL_LOCAL_S9B_ABLATION_${k}_" _scratch/s9b_build/r7_harness_output/wl/full-final/stdout | LC_ALL=C grep -o "$P" | wc -l); done; LC_ALL=C grep -o "$P" _scratch/s9b_build/r7_harness_output/wl/full-final/stdout | wc -l
K1 14
K2 12
K3 10
K4a 10
K4b 0
K4c 12
K5 0
K6 0
K7 14
K8 14
K9a 4
K9b 4
K10 0
K11 0
94
```

```
$ head -2 _scratch/s9b_build/s9b_repair_build_review_r7_grok_evidence/s9b_r7_equality_false_leaves.txt
leaves_by_knife {'K1': 14, 'K2': 12, 'K3': 10, 'K4a': 10, 'K4c': 12, 'K7': 14, 'K8': 14, 'K9a': 4, 'K9b': 4}
leaf_total 94
```

## R4 (Claude 3; Grok, not raised): SymPy K7 and K10 print a change only at the build-input digest
Inventory: harness tags of each knife whose markers include CHANGE, and the stored per-knife counts:
```
$ python3 -c "import json;d=json.load(open('research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json'));h=[r for r in d['runs'] if r['id']=='harness'][0];print(h['output'],len(h['tags']));[print(k,[t['name'] for t in h['tags'] if t['name'].startswith('PY_LOCAL_S9B_ABLATION_'+k+'_') and 'CHANGE' in t['markers']]) for k in ('K7','K10')]"
_scratch/s11c/s9b-amend10-full9/stdout 5321
K7 ['PY_LOCAL_S9B_ABLATION_K7_PY_LOCAL_S9B_BUILD_INPUT_DIGESTS']
K10 ['PY_LOCAL_S9B_ABLATION_K10_PY_LOCAL_S9B_BUILD_INPUT_DIGESTS']
```

```
$ python3 -c "import json;d=json.load(open('research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json'));[print(k,{x:d['difference_counts_by_knife'][k][x] for x in ('value_change_count','domain_change_count','value_undecided_count','domain_undecided_count')}) for k in ('K7','K10')]"
K7 {'value_change_count': 1, 'domain_change_count': 0, 'value_undecided_count': 130, 'domain_undecided_count': 376}
K10 {'value_change_count': 1, 'domain_change_count': 0, 'value_undecided_count': 0, 'domain_undecided_count': 1145}
```

The round-6 disposition listed SymPy K7 and K10 among the knives with no change at their site:
```
$ grep -o 'SymPy (no change word) K4a, K4c, K5, K7, K8, K9a, K9b and K10' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r6_review_disposition.md
SymPy (no change word) K4a, K4c, K5, K7, K8, K9a, K9b and K10
```

The legs' witness files for these sites (verdict lines only):
```
$ grep -c 'CHANGED (exact)' _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/relation_witness.stdout _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/k10_domain_witness.stdout
_scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/relation_witness.stdout:3
_scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/k10_domain_witness.stdout:3
```

## Observations
### Reuse: the Claude leg's guarded reruns of the reused baseline and K1 engine outputs
```
$ cat _scratch/s9b_build/s9b_repair_build_review_r7_claude_evidence/reuse_rerun_compare.txt
b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142  /tmp/s9b_review_r7_claudeleg/baseline_rerun/stdout
4cf4fcd8dcf7fd1fa8f62bdf5258007d7608c214e8bc8bfbf82fba569a2f689a  /tmp/s9b_review_r7_claudeleg/k1_rerun/stdout
baseline-identical
K1-identical
```

### K11: Wolfram tags reached and not evaluated; SymPy not-evaluated count
```
$ wc -l < research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_k11_evaluated_tags.txt; wc -l < research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round7_k11_not_evaluated_tags.txt
10
364
```

```
$ python3 -c "import json;d=json.load(open('research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_tag_inventory.json'));print(d['difference_counts_by_knife']['K11']['not_evaluated_count'])"
376
```

### Builder models for every harness round so far (Codex session headers)
```
$ grep -h '^model:\|^reasoning effort:' _scratch/s9b_build/repair_build/build_*.txt | sort | uniq -c
     38 model: gpt-6-astra
     38 reasoning effort: high
```

## Guard state after both legs
```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json
{}
```

```
$ systemctl --user list-units 's11c-guard-*' --no-legend | wc -l
0
```

