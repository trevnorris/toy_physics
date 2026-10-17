# Measurements — amendment 12 review, round 3 (generated 2026-10-16 18:56)

Generator: `_scratch/s9b_build/gen/s9b_amend12_r3_lookups.sh` (sha256sum/grep/sed/cat/git show only). Line
numbers in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r3.md` are the reviewed round-3 text's (saved before launch); the engines' and harnesses' are the
committed files' at `86097a44`.

## The reviewed text and the baseline
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r3.sha256 | cut -c1-64
1a11a74f0b1e5a12d7f9b8ff9a9206a6fd5580c53ec3cacdcd1abed69b1b3052
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r3.md | cut -c1-64
1a11a74f0b1e5a12d7f9b8ff9a9206a6fd5580c53ec3cacdcd1abed69b1b3052
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r3.sha256 | grep -c ': OK$'
13
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend12_review_prompt_r3.md | cut -c1-64
04bbe3489abb16de2caf00eb861dc0572ab29245b13238cf16afdd9acf5d79af
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend12_review_r3_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend12_review_r3_grok.txt
Verdict: NEEDS REVISION
```

## The round-3 clauses
```
$ sed -n '288,289p;293p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r3.md
      enters. Exclude the build-input digest and parts computed before the mutation enters. For K10, the result
      is the gate after negation. For K11, the witness compares compact corrupted/unmutated parts of a reached tag
      Take the union across copies (K11's compact pair), without filtering by structural marks. Print the set and
```

```
$ sed -n '300,322p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r3.md
      Include reached text only when the engine writes it as a literal status or label with no mathematical reading.
      Exclude parts whose difference comes from expression names or structure, or the knife's replacement text; print
      the reason. For text parts replace mathematical witness requirements with literal comparison: print both texts
      and location, count the part tried whether equal or different, and count a literal difference as a witness in
      both engines. One-copy text parts follow the one-copy rule.

      For mathematical parts count a part tried only when all witness requirements except difference hold. Decide
      values and conditions by engine CAS; substitution alone and unresolved evaluations do not count. Count an attempt
      only after running its witness test's CAS operations at printed inputs: solve fixing equations and evaluate values
      and attached conditions. Before listing an untried part, print the actual unresolved evaluation or requirement the
      CAS decided cannot be met; without those evaluations, continue the search. Print each copy's CAS-evaluated output
      of the construction defined by the knife list at trial inputs; exclude inputs giving the same output.

      Try item 9's shared object on each sign even if stored once; label trials do not try it. Count a sign trial
      only at a printed concrete input or throughout a printed subrange of the part's stated interval, with printed
      density evaluations in that sign in each emitting copy. Respect the spec's profile definitions, including ρ₀
      as asymptotic bulk density. If no such input meets the sign on that interval, print the CAS decision and unmet
      requirement; retain the part as untried, with no witness.

      End each search after finitely many attempts with exactly one outcome: a witness; otherwise, after attempting
      every part, **none found** if all count as tried, or **undecided**, listing each untried part and unmet
      requirement. Print each attempt's tag, location, role, inputs, evaluations or **not emitted**, whether tried, and
      unmet requirements.
```

## The orchestrator's round-3 entry 2, which the clause follows
```
$ grep -n 'A part whose' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_directive_amend12_review_disposition.md
220:2. **Text entries (N).** A text witness is text the engine writes as a literal status or label. A part whose
```

## The exclusion's scope (Claude 1, Grok 1)
```
$ sed -n '29p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
 "K10" -> {"branchGateLocal = propagatingLocal;", "branchGateLocal = Not[propagatingLocal];"},
```

```
$ sed -n '42p;45p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    'K8': ('mass_density = rho', "mass_density = sp.Symbol('s9b_constant_density',positive=True)"),
    'K10': ('observable_gate = real_domain', 'observable_gate = sp.Not(real_domain)'),
```

```
$ sed -n '597p;770,771p;903,904p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    oneform = (oneform*(1-return_orientation)/2).applyfunc(sp.factor)
    observable_gate = real_domain  # K10: existence only
    branch_domain=sp.And(observable_gate,travel_domain)
    mass_density = rho  # K8
    exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)
```

```
$ sed -n '166p;498,500p;512,513p;555,556p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
branchGateLocal = propagatingLocal; (* K10 *)
  If[quantifierMode === "FIXED_B", Return[<|"Status" -> "velocity undetermined by one fixed-b condition",
    "FreeInputs" -> DeleteDuplicates[Join[{velocityProfile,rhoBr},
      Cases[residual,deltaProfile | xiProfile | fProfile,{0,Infinity},Heads -> True]]], "Constraint" -> residual == 0,
    "FreeInputs" -> DeleteDuplicates[Join[{velocityProfile,rhoBr},
      Cases[odeD,deltaProfile | xiProfile | fProfile,{0,Infinity},Heads -> True]]], "ConstraintOnOtherInputs" -> (odeD == 0),
    "RemainingInputs" -> DeleteDuplicates[Cases[{ySolution,massDensity},
      deltaProfile | xiProfile | fProfile | rhoBr, {0,Infinity}, Heads -> True]]|>];
```

```
$ grep -n 'K8_RHS\|NAMES_ONLY\|VALUE_BASELINE_AT_r2\|K10_ONLY_CONJUNCT\|IMPLIED_JN_DEFLECTION_Density' _scratch/s9b_build/s9b_amend12_review_r3_claude_evidence/out/exclusion_scope.stdout
3:K8_RHS s9b_constant_density*(-s9b_r*Derivative(s9b_V(s9b_r), s9b_r) - 2*s9b_V(s9b_r))/s9b_r
4:NAMES_ONLY_BASELINE ['rho_br_profile']
5:NAMES_ONLY_K8 ['s9b_constant_density']
6:K8_RHS_WITH_KNIFE_SYMBOL_RENAMED_TO_RHO_BR (-s9b_r*Derivative(s9b_V(s9b_r), s9b_r) - 2*s9b_V(s9b_r))*rho_br_profile(s9b_r)/s9b_r
8:VALUE_BASELINE_AT_r2 1/128 VALUE_K8_AT_r2 0 EXACT_DIFFERENCE 1/128
14:K10_ONLY_CONJUNCT_HEAD Unequality IS_NOT_OF_A_BASELINE_ONLY_CONJUNCT True
16:K10_ONLY_CONJUNCT_STR Ne(ConditionSet(s9b_r, mu_perp(s9b_r)/((Derivative(s9b_xi(s9b_r), s9b_r)**2 + 1)*rho_br_profile(s9b_r)) <= 0, Union(ConditionSet(s9b_r, Ne(Intersection(ConditionSet((s9b_turn_ER, s9b_momentum_ER), Eq((s9b_momentum_ER**2*(Derivative(s9b_xi(s9b_turn_ER), s9b_turn_ER)**2 + 1)*s9b_V(s9b_turn_ER)**2 - s9
18:K8_IMPLIED_JN_DEFLECTION_Density_BASELINE ['rhoBr[r]', 'rhoBr[r]', 'rhoBr[r]']
19:K8_IMPLIED_JN_DEFLECTION_Density_CORRUPTED ['rhoBrConstant', 'rhoBrConstant', 'rhoBrConstant']
```

```
$ grep -n 'HEADLIST_EXAMPLE' _scratch/s9b_build/s9b_amend12_review_r3_claude_evidence/out/headlists_inert.stdout | head -3
45:HEADLIST_EXAMPLE ('K1', 'FreeInputs') WL_S9B_B_IMPLIED_JN_DEFLECTION BASELINE ['{deltaProfile, xiProfile, rhoBr}'] CORRUPTED ['{velocityProfile, rhoBr, deltaProfile, xiProfile}']
46:HEADLIST_EXAMPLE ('K2', 'FreeInputs') WL_S9B_B_IMPLIED_JN_DEFLECTION BASELINE ['{deltaProfile, xiProfile, rhoBr}'] CORRUPTED ['{deltaProfile, rhoBr}']
47:HEADLIST_EXAMPLE ('K2', 'RemainingInputs') WL_S9B_B_CONDITION_DEFLECTION BASELINE ['{deltaProfile, xiProfile, rhoBr}'] CORRUPTED ['{deltaProfile, rhoBr}']
```

```
$ cat _scratch/s9b_build/s9b_amend12_review_r3_grok_evidence/exclusion_scope.stdout
PY K5 ('return_orientation = -1', 'return_orientation = 1')
WL K5 ('returnDirection = -1;', 'returnDirection = 1;')
PY K7 ('fixed_impact = None', "fixed_impact = sp.Symbol('s9b_fixed_b',positive=True)")
WL K7 ('quantifierMode = "EVERY_B";', 'quantifierMode = "FIXED_B";')
PY K8 ('mass_density = rho', "mass_density = sp.Symbol('s9b_constant_density',positive=True)")
WL K8 ('massDensity = rhoBr[r];', 'massDensity = rhoBrConstant;')
PY K10 ('observable_gate = real_domain', 'observable_gate = sp.Not(real_domain)')
WL K10 ('branchGateLocal = propagatingLocal;', 'branchGateLocal = Not[propagatingLocal];')
K5_EMITTED_FACTOR
  mathematical: True
  replacement: return_orientation = 1
  before: 1.00000000000000
  after: 0.0
  structure_differs True
  expression_names_differ False
  replacement_text_in_new_content False
  excluded_if_sentence_covers_every_part True
  excluded_if_sentence_covers_mathematical_parts True
  introduced: ["Float('0.0', precision=53)"]
  factor_values 1.00000000000000 0.0 differ True
K8_MASS_EXCHANGE
  mathematical: True
  replacement: mass_density = sp.Symbol('s9b_constant_density',positive=True)
  before: -s9b_V(s9b_r)*Derivative(s9b_rho(s9b_r), s9b_r) - s9b_rho(s9b_r)*Derivative(s9b_V(s9b_r), s9b_r)
  after: -s9b_constant_density*Derivative(s9b_V(s9b_r), s9b_r)
  structure_differs True
  expression_names_differ True
  replacement_text_in_new_content True
  excluded_if_sentence_covers_every_part True
  excluded_if_sentence_covers_mathematical_parts True
  introduced: ["Mul(Integer(-1), Symbol('s9b_constant_density', positive=True), Derivative(Function('s9b_V')(Symbol('s9b_r', positive=True)), Tuple(Symbol('s9b_r', positive=True), Integer(1))))", "Symbol('s9b_constant_density', positive=True)"]
  exchange_difference -s9b_constant_density*Derivative(s9b_V(s9b_r), s9b_r) + s9b_V(s9b_r)*Derivative(s9b_rho(s9b_r), s9b_r) + s9b_rho(s9b_r)*Derivative(s9b_V(s9b_r), s9b_r)
K8_NAME_LIST
  mathematical: True
  replacement: mass_density = sp.Symbol('s9b_constant_density',positive=True)
  before: (s9b_V, s9b_r, s9b_rho)
  after: (s9b_V, s9b_constant_density, s9b_r)
  structure_differs True
  expression_names_differ True
  replacement_text_in_new_content True
  excluded_if_sentence_covers_every_part True
  excluded_if_sentence_covers_mathematical_parts True
  introduced: ["Symbol('s9b_constant_density')", "Tuple(Symbol('s9b_V'), Symbol('s9b_constant_density'), Symbol('s9b_r'))"]
  only_in_copy_names ['s9b_constant_density']
K10_SYPY_GATE
  mathematical: True
  replacement: observable_gate = sp.Not(real_domain)
  before: real_domain & travel_domain
  after: travel_domain & ~real_domain
  structure_differs True
  expression_names_differ False
  replacement_text_in_new_content True
  excluded_if_sentence_covers_every_part True
  excluded_if_sentence_covers_mathematical_parts True
  introduced: ["And(Symbol('travel_domain'), Not(Symbol('real_domain')))", "Not(Symbol('real_domain'))"]
K10_WOLFRAM_GATE
  mathematical: True
  replacement: branchGateLocal = Not[propagatingLocal];
  before: propagatingLocal
  after: ~propagatingLocal
  structure_differs True
  expression_names_differ False
  replacement_text_in_new_content True
  excluded_if_sentence_covers_every_part True
  excluded_if_sentence_covers_mathematical_parts True
  introduced: ["Not(Symbol('propagatingLocal'))"]
K7_STATUS
  mathematical: False
  replacement: quantifierMode = "FIXED_B";
  before: velocity determined on signed components with smooth zero-set joining
  after: velocity undetermined by one fixed-b condition
  structure_differs True
  expression_names_differ False
  replacement_text_in_new_content False
  excluded_if_sentence_covers_every_part True
  excluded_if_sentence_covers_mathematical_parts False
  introduced: ["Str('velocity undetermined by one fixed-b condition')"]
  status_contains_FIXED_B False
```

## Inert inputs as attempts (Claude 2)
```
$ grep -n 'K8_CONSTRUCTION_OUTPUT' _scratch/s9b_build/s9b_amend12_review_r3_claude_evidence/out/headlists_inert.stdout
62:K8_CONSTRUCTION_OUTPUT round8_inputs rho_br=1+1/r BASELINE 1/128 COPY 0 DIFFERENCE 1/128
63:K8_CONSTRUCTION_OUTPUT_AS_FUNCTION_OF_r round8_inputs rho_br=1+1/r BASELINE 1/(8*s9b_r**4) COPY 0 DIFFERENCE 1/(8*s9b_r**4)
64:K8_CONSTRUCTION_OUTPUT rho_br constant 1 BASELINE 0 COPY 0 DIFFERENCE 0
65:K8_CONSTRUCTION_OUTPUT_AS_FUNCTION_OF_r rho_br constant 1 BASELINE 0 COPY 0 DIFFERENCE 0
```

