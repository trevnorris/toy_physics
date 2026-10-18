# Measurements — S9b build review, round 9 (generated 2026-10-17 18:34)

Generator: `_scratch/s9b_build/gen/s9b_build_r9_lookups.sh` (sha256sum/grep/sed/cat/head/awk line selection
only). Line numbers are the reviewed files' (baseline `_scratch/s9b_build/s9b_repair_build_review_baseline_r9.sha256`).

## The reviewed files and the prompt
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r9.sha256 | grep -c ': OK$'
33
```

```
$ wc -l < _scratch/s9b_build/s9b_repair_build_review_baseline_r9.sha256
33
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_review_prompt_r9.md | cut -c1-64
5d5d369becc652da4525ddc8c1ec192f4ce3e503a9321748da6a8dd897c1b39c
```

```
$ git diff --quiet 86097a44 -- research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py research/pde_ledger_v3/scripts/S9b_exports.py research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl && echo engines-and-export-unchanged-since-86097a44
engines-and-export-unchanged-since-86097a44
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_repair_build_review_r9_claude.md
## Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict:\*\* [A-Z][A-Z ]*' _scratch/s9b_build/s9b_repair_build_review_r9_grok.txt
Verdict:** NEEDS REVISION
```

## Claude 1: SymPy K10 attempts cannot return a witness at a mathematical part
```
$ awk 'NR>=782 && NR<=904 && /return base/ {print NR": "substr($0,1,160)}' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
796:         return base|{'status':'undecided','CAS_operations':('compact method does not reach this tag; no mathematical comparison',),
807:         return base|{'status':'no change' if result is sp.true else 'witness',
824:             return base|{'status':'no change','statement':statement,'CAS_result':decision,
883:     return base|{'status':'no change' if no_change else 'undecided','statement':statement,'CAS_result':decision,
```

```
$ sed -n '805,808p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    if knife!='K11' and all(p and p['role']=='text' for p in pair):
        result=sp.sympify(values[0]==values[1])
        return base|{'status':'no change' if result is sp.true else 'witness',
            'CAS_operations':('SymPy Str literal equality',),'texts':values,'CAS_result':result,
```

```
$ sed -n '893p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
            'This partial trial does not solve the joint fixing equations and their cross-copy intersection, decide every active quantified condition, or certify ordinary finite unequal values. No witness is inferred.')}
```

```
$ sed -n '904,905p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    if knife not in ('K7','K10'):
        trial=direct_search(knife,left,right,source)
```

```
$ git show 86097a44:research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py | grep -n 'cautious_condition_search' | head -3
439:def cautious_condition_search(knife,left_stream,right_stream,source):
818:                search=(cautious_condition_search if knife in ('K7','K10') else direct_search)(knife,baseline,corrupted,source)
```

```
$ sed -n '308,312p' research/pde_ledger_v3/directives/S9b_repair_build_directive.md

      For mathematical parts report **no change** only on an engine CAS decision of agreement between the copies in
      the part and its attached conditions at every admissible input. Report **cannot be tried** only on an
      engine CAS decision that no admissible input meets the witness requirements apart from difference. For every
      no-change or cannot-be-tried decision, print its statement, CAS result and full input coverage: all inputs and
```

```
$ grep -n 'K10' research/pde_ledger_v3/directives/_measurements/S9b_repair_sympy_build_report.md
7:**Searches:** witnesses: **K1, K2, K3, K4a, K4b, K4c, K5, K6, K7, K8, K9a, K9b, K11**; undecided: **K10**; none found: **none**. K7's witness compares a literal reduction status. K10 tried **4163** parts; **3022** lack a no-change decision and are listed with their unsettled requirements. Per-part decisions are printed; counts are inventoried. Marks and method residuals furnish no witness.
13:**Knife sites:** K1 advection L508; K2 metric L509; K3 speed L510; K4a–c derivative freezes L39; K5 return orientation L773; K6 slope source L820; K7 fixed impact L895; K8 mass density L903; K9a local sound speed L105; K9b bulk density L108; K10 gate L770; K11 azimuthal advection L506.
```

```
$ head -24 _scratch/s9b_build/s9b_repair_build_review_r9_claude_evidence/py_k10_gate_decide2.stdout
===== baseline
 conjunct 0 Unequality -> undecided by doit; root-set membership test:
   candidate (12, 12) base_set.contains -> True
   condition at candidate (inner ConditionSets resolved by solveset) -> True
      conj True -> True
 conjunct 1 Unequality -> undecided by doit; root-set membership test:
   candidate (12, 12) base_set.contains -> True
   condition at candidate (inner ConditionSets resolved by solveset) -> Eq(-asin(12/13) - asin(4/5) + pi, acos(33/65))
      Eq(trig) Eq(-asin(12/13) - asin(4/5) + pi, acos(33/65)) simplify(lhs-rhs)= -asin(12/13) - acos(33/65) - asin(4/5) + pi | exact cos(lhs)-cos(rhs) via expand_trig+simplify = 0 | lhs,rhs to 60 digits = 1.03829222849304590354290875910690322493668079869101658590098 1.03829222849304590354290875910690322493668079869101658590098 | both in (0,pi): True
 conjunct 2 Unequality -> undecided by doit; root-set membership test:
   candidate (12, 12) base_set.contains -> True
   condition at candidate (inner ConditionSets resolved by solveset) -> Eq(-asin(12/13) - asin(4/5) + pi, acos(33/65))
      Eq(trig) Eq(-asin(12/13) - asin(4/5) + pi, acos(33/65)) simplify(lhs-rhs)= -asin(12/13) - acos(33/65) - asin(4/5) + pi | exact cos(lhs)-cos(rhs) via expand_trig+simplify = 0 | lhs,rhs to 60 digits = 1.03829222849304590354290875910690322493668079869101658590098 1.03829222849304590354290875910690322493668079869101658590098 | both in (0,pi): True
 conjunct 3 Equality -> True
 conjunct 4 Equality -> True
===== corrupted
 conjunct 0 Unequality -> undecided by doit; root-set membership test:
   candidate (12, 12) base_set.contains -> True
   condition at candidate (inner ConditionSets resolved by solveset) -> True
      conj True -> True
 conjunct 1 Unequality -> undecided by doit; root-set membership test:
   candidate (12, 12) base_set.contains -> True
   condition at candidate (inner ConditionSets resolved by solveset) -> Eq(-asin(12/13) - asin(4/5) + pi, acos(33/65))
      Eq(trig) Eq(-asin(12/13) - asin(4/5) + pi, acos(33/65)) simplify(lhs-rhs)= -asin(12/13) - acos(33/65) - asin(4/5) + pi | exact cos(lhs)-cos(rhs) via expand_trig+simplify = 0 | lhs,rhs to 60 digits = 1.03829222849304590354290875910690322493668079869101658590098 1.03829222849304590354290875910690322493668079869101658590098 | both in (0,pi): True
```

```
$ sed -n '1,8p' _scratch/s9b_build/s9b_repair_build_review_r9_claude_evidence/py_conditionset_subs2.stdout
2-variable ConditionSet : Ne(ConditionSet((x, y), Eq(x, b/c_0), ProductSet(Interval.open(0, oo), Interval.open(0, oo))), EmptySet)
  subs(simultaneous=True): Ne(ConditionSet((x, y), Eq(x, b/c_0), ProductSet(Interval.open(0, oo), Interval.open(0, oo))), EmptySet)
  subs():                  Ne(ConditionSet((x, y), Eq(x, b/c_0), ProductSet(Interval.open(0, oo), Interval.open(0, oo))), EmptySet)
  xreplace():              Ne(ConditionSet((x, y), Eq(x, 12), ProductSet(Interval.open(0, oo), Interval.open(0, oo))), EmptySet)
1-variable ConditionSet : Ne(ConditionSet(x, Eq(x, b/c_0), Interval.open(0, oo)), EmptySet)
  subs(simultaneous=True): Ne(ConditionSet(x, Eq(x, b/c_0), Interval.open(0, oo)), EmptySet)
  subs():                  Ne(ConditionSet(x, Eq(x, b/c_0), Interval.open(0, oo)), EmptySet)
  xreplace():              Ne(ConditionSet(x, Eq(x, 12), Interval.open(0, oo)), EmptySet)
```

```
$ sed -n '1,4p' _scratch/s9b_build/s9b_repair_build_review_r9_claude_evidence/py_k10_attempt1819_inputs.stdout
chars 1094684
printed point: Tuple(Str('point'), Tuple(Tuple(Symbol('Z_E', positive=True), Integer(5)), Tuple(Symbol('Z_R', positive=True), Integer(9)), Tuple(Symbol('c_0', positive=True), Integer(1)), Tuple(Symbol('b', positive=True), Integer(12)))), Tuple(Str('ranged'), Tuple(Tuple(), T
evaluations segment chars 150546 | occurrences of Symbol('b', positive=True): 6 | Symbol('c_0', positive=True): 6
first occurrence context: Tuple(Str('evaluations'), Tuple(Tuple(Tuple(Str('part'), And(Unequality(ConditionSet(Tuple(Symbol('s9b_turn_FLYBY', positive=True), Symbol('s9b_momentum_FLYBY', positive=True)), And(Equality(Symbol('s9b_momentum_FLYBY', positive=True), Mul(Symbol('b', positive=True), Pow(Symbol('c_0', positive=True)
```

## Claude 2: the Wolfram K10 mutation is not the negated computed branch-existence condition
```
$ sed -n '393,394p' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
    - **K10, branch gate:** the computed branch-existence condition is replaced by its logical negation before the
      observables are gated on it.
```

```
$ sed -n '193p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
```

```
$ sed -n '165,166p;174p;203p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
propagatingLocal = Element[relativeNorm, Reals] && relativeNorm > 0 && Det[metric2] > 0;
branchGateLocal = propagatingLocal; (* K10 *)
physicalPropagation = branchGateLocal /. branchRules;
  "Domain" -> (r > 0 && physicalTraversal && physicalPropagation && speedChart),
```

```
$ sed -n '261,265p;275p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
traversedRadii = flybyRadii || radarRadii;
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
  physicalPropagation && speedChart && rhoBr[r] != 0 &&
  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];
traversalGate = rayConstruction && Inactive[ForAll][r, traversedRadii, physicalTraversal];
  emit[prefix <> "BRANCH_EXISTENCE", wrapper[propagationGate /. rules]];
```

```
$ sed -n '696p;720p;770,771p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    real_domain = everywhere(r,existence.subs(branch_map),traversed)
    emit('BRANCH_EXISTENCE',real_domain)
    observable_gate = real_domain  # K10: existence only
    branch_domain=sp.And(observable_gate,travel_domain)
```

```
$ head -6 _scratch/s9b_build/s9b_repair_build_review_r9_claude_evidence/wl_k10_branch_existence.stdout
baseline chars 4472 | ForAll occurrences 3 | Not[ occurrences 0 | "!(" occurrences 0
   ForAll body start: Inactive[ForAll][r, r >= rTurnFlyby, ((2*c0^2*r)/((c0 + c0*deltaProfile[r])^2 - velocityProfile[r]^2*(1 + Derivative[1][xiProfile][r]^2)) - (c0^2*r^2*(2*c0*(c0 + c0*deltaProfile[r])*Derivative[1][deltaProfile][r] - 2*velocityProfile[r]*Derivative[1][velocityProfile][r]*(1 + Derivative[1][xiProfile][r]^2) - 2*velocityProfile[r]^2*Derivative[1][xiProfile][r]*Derivative[2][xiPro
corrupted chars 4476 | ForAll occurrences 3 | Not[ occurrences 0 | "!(" occurrences 1
   ForAll body start: Inactive[ForAll][r, r >= rTurnFlyby, ((2*c0^2*r)/((c0 + c0*deltaProfile[r])^2 - velocityProfile[r]^2*(1 + Derivative[1][xiProfile][r]^2)) - (c0^2*r^2*(2*c0*(c0 + c0*deltaProfile[r])*Derivative[1][deltaProfile][r] - 2*velocityProfile[r]*Derivative[1][velocityProfile][r]*(1 + Derivative[1][xiProfile][r]^2) - 2*velocityProfile[r]^2*Derivative[1][xiProfile][r]*Derivative[2][xiPro
corrupted: negation inside ForAll starting at 4116 ; negation at 4216
   context: Inactive[ForAll][r, r >= rTurnFlyby || r ... TurnRadar <= r <= Max[Sqrt[b^2 + ZE^2], Sqrt[b^2 + ZR^2]],  !(Element[muPerp[r]/rhoBr[r], Reals] && muPerp[r]/rhoBr[r] > 0 && r^2*(1 + Derivative[1][xiProfile][r]^2) > 0) && c0*(1 + deltaProfile[r]) == Sqrt[muPerp[r]/rhoBr[r]] && rhoBr[r] != 0 && Element[velocityProfile[r] | xiProfile[r] | muPerp[r] | rhoBr[r], Reals]]
```

## Claude 3: a stale sentence in the rewritten round-8 measurements
```
$ sed -n '5p' research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round8_measurements.md
The [current round-9 amendment-12 report](S9b_repair_wolfram_build_report.md) and [measurements](S9b_repair_wolfram_round9_measurements.md) record the actual resource stop and incomplete coverage. No previous successful run supplies a current result.
```

```
$ sed -n '116p' research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round9_measurements.md
All four runs exited **0**; no guard refusal/kill, memory-exhaustion, task-limit or item-12 event occurred. Full-run peak cgroup memory **2863075328 bytes**, maximum sampled tasks **32**, sampled swap **0 bytes**. These are cgroup memory figures, not a native virtual-address-space peak. Effective cgroup and native address-space limits were both **13958643712 bytes**. Per-run resource summaries and
```

## Claude 4 and Grok 1: Profiles and IntegrationConstants
```
$ sed -n '303,305p' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
      Exclude exactly two additional kinds of part in both engines, printing each reason: lists of symbol, function
      or row names extracted from expressions, whether emitted as text or symbols; and text parts whose literal
      difference is the knife's replacement source. For text parts replace mathematical tests with engine CAS literal
```

```
$ sed -n '397p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
extractedNameKeys={"Profiles","FreeInputs","RemainingInputs"};
```

```
$ sed -n '401p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  "Profiles" -> {deltaProfile,velocityProfile,xiProfile,rhoBr,fProfile},
```

```
$ sed -n '537p;542p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
      "IntegrationConstants" -> If[homogeneous === 0,{}, {velocityIntegrationConstant}],
    "IntegrationConstants" -> If[homogeneous === 0,{}, {velocityIntegrationConstant}],
```

```
$ sed -n '499,500p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
    "FreeInputs" -> DeleteDuplicates[Join[{velocityProfile,rhoBr},
      Cases[residual,deltaProfile | xiProfile | fProfile,{0,Infinity},Heads -> True]]], "Constraint" -> residual == 0,
```

```
$ head -12 _scratch/s9b_build/s9b_repair_build_review_r9_claude_evidence/wl_profiles_entries_rerun.stdout
K1 50 ('"IntegrationConstants"', 'value', 'site part', '"emitted", "emitted"', '')
K1 42 ('"IntegrationConstants"', 'value', 'site part', '"emitted", "not emitted"', '')
K1 67 ('"Profiles"', 'domain', 'excluded', '"emitted", "emitted"', 'additional exclusion 1: lists of symbol/function/row names e')
K2 50 ('"IntegrationConstants"', 'value', 'site part', '"emitted", "emitted"', '')
K2 67 ('"Profiles"', 'domain', 'excluded', '"emitted", "emitted"', 'additional exclusion 1: lists of symbol/function/row names e')
K3 50 ('"IntegrationConstants"', 'value', 'site part', '"emitted", "emitted"', '')
K3 67 ('"Profiles"', 'domain', 'excluded', '"emitted", "emitted"', 'additional exclusion 1: lists of symbol/function/row names e')
K4a 50 ('"IntegrationConstants"', 'value', 'site part', '"emitted", "emitted"', '')
K4a 67 ('"Profiles"', 'domain', 'excluded', '"emitted", "emitted"', 'additional exclusion 1: lists of symbol/function/row names e')
K4b 50 ('"IntegrationConstants"', 'value', 'site part', '"emitted", "emitted"', '')
K4b 67 ('"Profiles"', 'domain', 'excluded', '"emitted", "emitted"', 'additional exclusion 1: lists of symbol/function/row names e')
K4c 50 ('"IntegrationConstants"', 'value', 'site part', '"emitted", "emitted"', '')
```

```
$ head -8 _scratch/s9b_build/s9b_repair_build_review_r9_grok_evidence/stdout/wl_integration_constants.txt
emitted_count 176
126 "IntegrationConstants" -> {velocityIntegrationConstant}
34 "IntegrationConstants" -> {}, "AlongRayDomain" -> Condi
16 "IntegrationConstants" -> {}, "AlongRayDomain" -> Inact
```

```
$ head -8 _scratch/s9b_build/s9b_repair_build_review_r9_grok_evidence/stdout/wl_site_fields.txt
{'head': 'WL_LOCAL_S9B_ABLATION_K1_SITE_INVENTORY', 'IntegrationConstants_site': 436, 'IntegrationConstants_excl': 0, 'Profiles_site': 0, 'Profiles_excl': 67, 'bytes': 8299008, 'BRANCH_EXISTENCE_site': 8, 'BRANCH_EXISTENCE_excl': 18, 'BRANCH_TYPE_site': 185, 'BRANCH_TYPE_excl': 34, 'PATH_TRAVERSAL_site': 8, 'PATH_TRAVERSAL_excl': 18, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
{'head': 'WL_LOCAL_S9B_ABLATION_K2_SITE_INVENTORY', 'IntegrationConstants_site': 352, 'IntegrationConstants_excl': 0, 'Profiles_site': 0, 'Profiles_excl': 67, 'bytes': 8062038, 'BRANCH_EXISTENCE_site': 8, 'BRANCH_EXISTENCE_excl': 18, 'BRANCH_TYPE_site': 205, 'BRANCH_TYPE_excl': 34, 'PATH_TRAVERSAL_site': 8, 'PATH_TRAVERSAL_excl': 18, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
{'head': 'WL_LOCAL_S9B_ABLATION_K3_SITE_INVENTORY', 'IntegrationConstants_site': 352, 'IntegrationConstants_excl': 0, 'Profiles_site': 0, 'Profiles_excl': 67, 'bytes': 8311450, 'BRANCH_EXISTENCE_site': 8, 'BRANCH_EXISTENCE_excl': 18, 'BRANCH_TYPE_site': 233, 'BRANCH_TYPE_excl': 34, 'PATH_TRAVERSAL_site': 8, 'PATH_TRAVERSAL_excl': 18, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
{'head': 'WL_LOCAL_S9B_ABLATION_K4a_SITE_INVENTORY', 'IntegrationConstants_site': 352, 'IntegrationConstants_excl': 0, 'Profiles_site': 0, 'Profiles_excl': 67, 'bytes': 8029231, 'BRANCH_EXISTENCE_site': 0, 'BRANCH_EXISTENCE_excl': 26, 'BRANCH_TYPE_site': 0, 'BRANCH_TYPE_excl': 176, 'PATH_TRAVERSAL_site': 0, 'PATH_TRAVERSAL_excl': 26, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
{'head': 'WL_LOCAL_S9B_ABLATION_K4b_SITE_INVENTORY', 'IntegrationConstants_site': 352, 'IntegrationConstants_excl': 0, 'Profiles_site': 0, 'Profiles_excl': 67, 'bytes': 7994391, 'BRANCH_EXISTENCE_site': 0, 'BRANCH_EXISTENCE_excl': 26, 'BRANCH_TYPE_site': 0, 'BRANCH_TYPE_excl': 176, 'PATH_TRAVERSAL_site': 0, 'PATH_TRAVERSAL_excl': 26, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
{'head': 'WL_LOCAL_S9B_ABLATION_K4c_SITE_INVENTORY', 'IntegrationConstants_site': 352, 'IntegrationConstants_excl': 0, 'Profiles_site': 0, 'Profiles_excl': 67, 'bytes': 8060069, 'BRANCH_EXISTENCE_site': 8, 'BRANCH_EXISTENCE_excl': 18, 'BRANCH_TYPE_site': 170, 'BRANCH_TYPE_excl': 54, 'PATH_TRAVERSAL_site': 8, 'PATH_TRAVERSAL_excl': 18, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
{'head': 'WL_LOCAL_S9B_ABLATION_K5_SITE_INVENTORY', 'IntegrationConstants_site': 176, 'IntegrationConstants_excl': 108, 'Profiles_site': 0, 'Profiles_excl': 25, 'bytes': 8330454, 'BRANCH_EXISTENCE_site': 0, 'BRANCH_EXISTENCE_excl': 2, 'BRANCH_TYPE_site': 0, 'BRANCH_TYPE_excl': 52, 'PATH_TRAVERSAL_site': 0, 'PATH_TRAVERSAL_excl': 0, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
{'head': 'WL_LOCAL_S9B_ABLATION_K6_SITE_INVENTORY', 'IntegrationConstants_site': 176, 'IntegrationConstants_excl': 104, 'Profiles_site': 0, 'Profiles_excl': 21, 'bytes': 8415000, 'BRANCH_EXISTENCE_site': 0, 'BRANCH_EXISTENCE_excl': 2, 'BRANCH_TYPE_site': 0, 'BRANCH_TYPE_excl': 52, 'PATH_TRAVERSAL_site': 0, 'PATH_TRAVERSAL_excl': 0, 'LOCAL_NAMES_site': 0, 'LOCAL_NAMES_excl': 2}
```

## Grok's K10 rows
```
$ grep -n '| K10 |' _scratch/s9b_build/s9b_repair_build_review_r9_grok.txt
27:| K10 | undecided: 4163 parts, 1141 no change, 3022 undecided | confirms | witness, local `Domain`, False vs True | confirms |
```

```
$ grep -n 'K10_CONJUNCTS\|RESULT K10' _scratch/s9b_build/s9b_repair_build_review_r9_grok.txt
37:K10_CONJUNCTS corrupted traversal,propagation,speed,r>0 {True, False, True, True}
38:K10_CONJUNCTS baseline traversal,propagation,speed,r>0 {True, True, True, True}
39:RESULT K10 Domain recomputed={False, True} printed={False, True} equal=True
```

