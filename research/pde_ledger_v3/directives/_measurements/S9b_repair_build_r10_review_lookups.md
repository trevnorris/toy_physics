# Measurements: S9b build review, round 10 (generated 2026-10-18 17:16)

Generator: `_scratch/s9b_build/gen/s9b_build_r10_lookups.sh` (sha256sum/grep/sed/head/awk line selection/wc
only). The legs' evidence is copied to `_scratch/s9b_build/s9b_repair_build_review_r10_{claude,grok}_evidence/`.

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_baseline_r10.sha256 | grep -c ': OK$')
34
```

```
$ wc -l < /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_baseline_r10.sha256
34
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_prompt_r10.md | cut -c1-64
5d43f2d85e4cb0708d906ede2ad8b247060d6f5c4dc1345789cf28698978af9d
```

```
$ git -C /var/projects/toy_physics diff --quiet 86097a44 -- /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_exports.py /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl && echo engines-and-export-unchanged-since-86097a44
engines-and-export-unchanged-since-86097a44
```

## Verdicts
```
$ head -1 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_grok.txt
Verdict: CLEAR
```

## The knife's definition
```
$ sed -n '393,394p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_repair_build_directive.md
    - **K10, branch gate:** the computed branch-existence condition is replaced by its logical negation before the
      observables are gated on it.
```

```
$ sed -n '193,196p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
**Branch existence.** First print both the condition under which the branch is real and propagating everywhere
along the ray and the condition under which a ray can traverse the full flyby path and both legs of the round
trip in the required directions. Outside either condition, report the branch type (growing, decaying, absent,
or unable to traverse in a required direction), print `NOT_ESTABLISHED` for the observables, and do not compute
```

```
$ sed -n '288,290p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_repair_build_directive.md
      uses the mutated construction's result, with its payload location and the source location where that result
      enters. Exclude the build-input digest and parts computed before the mutation enters. For K10, the result
      is the gate after negation. For K11, the witness compares compact corrupted/unmutated parts of a reached tag
```

```
$ grep -n 'the same object' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r9_review_disposition.md | cut -c1-220
56:- **Wolfram:** the K10 mutation is the negation of the engine's computed branch-existence condition, the same object
```

## Finding 1: each engine's K10 copy negates a different object
```
$ sed -n '29p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | cut -c1-260
 "K10" -> {"propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,\n  physicalPropagation && speedChart && rhoBr[r] != 0 &&\n  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];", "propagationGate = Not[rayConstruction && I
```

```
$ sed -n '166p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
branchGateLocal = propagatingLocal; (* K10 *)
```

```
$ sed -n '262,266p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
  physicalPropagation && speedChart && rhoBr[r] != 0 &&
  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];
traversalGate = rayConstruction && Inactive[ForAll][r, traversedRadii, physicalTraversal];
rayGate = propagationGate && traversalGate;
```

```
$ sed -n '275,276p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  emit[prefix <> "BRANCH_EXISTENCE", wrapper[propagationGate /. rules]];
  emit[prefix <> "PATH_TRAVERSAL", wrapper[traversalGate /. rules]];
```

```
$ sed -n '46p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    'K10': ('observable_gate = real_domain', 'observable_gate = sp.Not(real_domain)'),
```

```
$ sed -n '696,697p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py | cut -c1-200
    real_domain = everywhere(r,existence.subs(branch_map),traversed)
    travel_domain = sp.And(*(sp.Ne(ray['roots'],sp.EmptySet,evaluate=False) for ray in rays.values()),
```

```
$ sed -n '720p;770,771p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    emit('BRANCH_EXISTENCE',real_domain)
    observable_gate = real_domain  # K10: existence only
    branch_domain=sp.And(observable_gate,travel_domain)
```

```
$ grep -h 'corrupted = Not\|BRANCH_TYPE baseline SameQ\|PATH_TRAVERSAL baseline SameQ\|BRANCH_EXISTENCE {baseline' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_claude_evidence/wl_scripts/check_wl_k10*.stdout | head -6
corrupted = Not[baseline] literally: True; head of corrupted: Not
BRANCH_TYPE baseline SameQ corrupted: True
PATH_TRAVERSAL baseline SameQ corrupted: True
BRANCH_EXISTENCE {baseline, corrupted} whole emitted objects at the input: InputForm[{True, False}]; harness printed {corrupted, baseline} = InputForm[{False, True}]
```

```
$ grep -n 'K10 is not the same formula' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_grok.txt | cut -c1-300
61:K10 is not the same formula in the two engines, and each harness follows its own engine. SymPy sets `observable_gate = Not(real_domain)` and then `branch_domain = And(observable_gate, travel_domain)`. `travel_domain` stays in the conjunction. `BRANCH_EXISTENCE` and `BRANCH_TYPE` are emitted befor
```

## Finding 2: the SymPy K10 attempts cannot decide the gated domain parts
```
$ sed -n '343p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
        if isinstance(node,sp.ConditionSet) and isinstance(node.sym,sp.Symbol):
```

```
$ grep -n "check('scalar_solver_preserves_nested_ray_binders'" /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
1321:    check('scalar_solver_preserves_nested_ray_binders',{
```

```
$ sed -n '232,238p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
PROFILES = {'s9b_delta':1/(10*RADIUS**2), 's9b_V':1/(8*RADIUS**2),
            's9b_xi':1/(6*RADIUS), 'rho_br_profile':1+1/RADIUS,
            's9b_f':1/(10*RADIUS**2), 'mu_perp':(1+1/RADIUS)*(1+1/(10*RADIUS**2))**2}
POINT = {'s9b_X':2,'s9b_r':2,'c_0':1,'b':12,'Z_E':5,'Z_R':9,
         's9b_b_min':1,'s9b_fixed_b':12,'s9b_constant_density':1,
         'rho_0':1,'K':1,'m':1,'n':3,'s':1,'GM':1,
         'x_1':3,'x_2':4,'x_3':0}
```

```
$ head -12 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_claude_evidence/py_scripts/k10_harness_attempt_probe.stdout
(1) harness PROFILES: status: undecided
  copy 0 evaluated domain: type And; conjunct heads: ['Unequality', 'Unequality', 'Unequality', 'Equality'] ; tuple-binder ConditionSets remaining: 6
  copy 1 evaluated domain: type BooleanFalse; conjunct heads: ['BooleanFalse'] ; tuple-binder ConditionSets remaining: 0
  unmet: ('Ordinary finite exact evaluations not established in both copies; unevaluated binders/integrals/derivatives decide nothing.', 'Compared domain is not decided true/false in both copies.', 'No full-input agreement or impossibility decision.')
(2) flat admissible profile: status: undecided
  copy 0 evaluated domain: type And; conjunct heads: ['Unequality', 'Unequality', 'Unequality'] ; tuple-binder ConditionSets remaining: 3
  copy 1 evaluated domain: type BooleanFalse; conjunct heads: ['BooleanFalse'] ; tuple-binder ConditionSets remaining: 0
  unmet: ('Ordinary finite exact evaluations not established in both copies; unevaluated binders/integrals/derivatives decide nothing.', 'Compared domain is not decided true/false in both copies.', 'No full-input agreement or impossibility decision.')
```

```
$ grep -h '(a) flat\|conjunction at element\|decides d = 0' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_claude_evidence/py_scripts/k10_flat_step2.stdout | head -6
(a) flat: corrupted copy evaluates to False
     conjunction at element -> True
     conjunction at element -> Eq(-pi + asin(4/5) + asin(12/13), -acos(33/65))
     exact decision of the integral equality: d = lhs - rhs = -pi + asin(4/5) + acos(33/65) + asin(12/13) ; cos(d) by expand_trig+simplify = 1 ; d to 50 digits = -0.e-160 ; (cos(d)=1 exactly and |d|<2*pi) decides d = 0: True
     conjunction at element -> Eq(-pi + asin(4/5) + asin(12/13), -acos(33/65))
     exact decision of the integral equality: d = lhs - rhs = -pi + asin(4/5) + acos(33/65) + asin(12/13) ; cos(d) by expand_trig+simplify = 1 ; d to 50 digits = -0.e-160 ; (cos(d)=1 exactly and |d|<2*pi) decides d = 0: True
```

```
$ cat /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_claude_evidence/py_scripts/k10_tuple_binder_count.stdout | head -6
252 ('no change', 'evaluated objects field present', 'none remains')
889 ('no change', 'no evaluated objects field', 'none remains')
1763 ('undecided', 'evaluated objects field present', 'none remains')
1259 ('undecided', 'evaluated objects field present', 'tuple-binder ray ConditionSet remains')
```

```
$ cat /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_claude_evidence/py_scripts/k10_unmet_inputs_count.stdout | head -8
Unassigned independent inputs: {('s9b_velocity_zero_radius',): 672, ('Phi',): 710, ('s9b_gamma',): 8, ('Phi', 's9b_gamma'): 4} total attempts 1394
Unfixed profile functions: {('j_n',): 1046} total attempts 1046
```

```
$ grep -n 'undecided' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r10_grok.txt | head -3 | cut -c1-300
3:No finding changes what the harnesses or the SymPy tag inventory compute or show. Every printed witness that was recomputed matches. K10 on the SymPy side is undecided, which is what its attempt counts imply.
27:| SymPy | K10 | undecided; 4163 attempts; 3022 without a no-change decision | confirms the outcome label |
61:K10 is not the same formula in the two engines, and each harness follows its own engine. SymPy sets `observable_gate = Not(real_domain)` and then `branch_domain = And(observable_gate, travel_domain)`. `travel_domain` stays in the conjunction. `BRANCH_EXISTENCE` and `BRANCH_TYPE` are emitted befor
```

