# Measurements — amendment 12 review, round 1 (generated 2026-10-16 17:42)

Generator: `_scratch/s9b_build/gen/s9b_amend12_r1_lookups.sh` (sha256sum/grep/sed/cat/git show only). Line
numbers in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r1.md` are the reviewed round-1 text's (saved before launch); the engines' are the committed files'
at `86097a44`.

## The reviewed text and the baseline
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r1.sha256 | cut -c1-64
7375a6d6ad693acfcd7adc96b355aa5d75b39384c723947295918acd51debd0f
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r1.md | cut -c1-64
7375a6d6ad693acfcd7adc96b355aa5d75b39384c723947295918acd51debd0f
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r1.sha256 | grep -c ': OK$'
11
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend12_review_prompt_r1.md | cut -c1-64
8872e3b9ef7ad7189efd7f25c98fbfcce5322eaaa4f74a812e5408a44f81d023
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend12_review_r1_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend12_review_r1_grok.txt
Verdict: NEEDS REVISION
```

## The round-1 clauses
```
$ sed -n '292,315p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r1.md
      Enumerate all such values and domains by mathematical object in both engines. Recurse through mappings'
      keyed values; matrices', vectors' and lists' components; solution sets' represented branches and domains;
      and cases' values and selectors. For equations, include the predicate and both sides. Keep other scalar
      expressions and predicates whole. Identify locations by tag and mathematical key, component, side, solution
      or case identity, never storage position in a sparse or reordered container. Include engine-local tags and
      both item 9 density-sign parts. Trace direct and transitive computational uses, and take the union of
      locations in the two copies (the compact pair for K11), without filtering by structural marks. Print the
      site-part set and every excluded value/domain part with its reason.

      For a part emitted by only one copy, evaluate that copy and its attached conditions; print **not emitted**
      for the other. Use and print an absent entry's reading only if the engine defines it, with its computation;
      supply no undefined value. Absence alone furnishes no witness; apply the witness tests below to any
      engine-defined readings.

      Count a part as tried only when its inputs, construction solutions, available values and attached conditions
      meet the witness rules below, apart from establishing a difference. Decide values and conditions by the engine
      CAS; substitution alone and unresolved evaluations do not count. Use no input choice that alone neutralizes
      the mutation at its construction; print the CAS computation checking this. For each density-sign part,
      compute and print the density at the trial inputs; count the trial only if it lies in that part in each
      emitting copy, throughout any ranged radius.

      Print one witness, or **none found** only after every site part counts as tried. For each attempt, print the
      tag, location, value/domain role, inputs, performed evaluations or **not emitted**, whether it counts as tried,
      and every unmet witness requirement.
```

## A part that cannot count as tried (Claude 1, Grok 1)
```
$ grep -n 'ConditionSet\|Integrals in its condition\|structurally equal' _scratch/s9b_build/s9b_amend12_review_r1_claude_evidence/out/py_labels_gates.stdout | head -8
28:  not a Str-keyed mapping in one copy; baseline head/len: Tuple 4 | K7 head/len: Tuple 4 | K7 value: (VELOCITY_FIXED_BY_RESTRICTION, Eq(s9b_V(s9b_r), 0), (ATTACHED_DOMAIN, Ne(ConditionSet((s9b_turn_FLYBY, s9b_momentum_FLYBY), Eq(s9b_momentum_FLYBY, b/c_0) & Eq((-s9b_momentum_FLYBY**2*mu_perp(s9b_turn_FLYBY)/rho_br_profile(s9b_turn_FLYBY) + s9b_turn_FLYBY**2)/s9b_turn_FLYBY**2, 0) & Eq(Condition
30:  not a Str-keyed mapping in one copy; baseline head/len: Tuple 4 | K7 head/len: Tuple 4 | K7 value: (VELOCITY_FIXED_BY_RESTRICTION, Eq(s9b_V(s9b_r), 0), (ATTACHED_DOMAIN, Ne(ConditionSet((s9b_turn_FLYBY, s9b_momentum_FLYBY), Eq(s9b_momentum_FLYBY, b/c_0) & Eq((-s9b_momentum_FLYBY**2*mu_perp(s9b_turn_FLYBY)/rho_br_profile(s9b_turn_FLYBY) + s9b_turn_FLYBY**2)/s9b_turn_FLYBY**2, 0) & Eq(Condition
32:  not a Str-keyed mapping in one copy; baseline head/len: Tuple 4 | K7 head/len: Tuple 4 | K7 value: (VELOCITY_FIXED_BY_RESTRICTION, Eq(s9b_V(s9b_r), 0), (ATTACHED_DOMAIN, Ne(ConditionSet((s9b_turn_FLYBY, s9b_momentum_FLYBY), Eq(s9b_momentum_FLYBY, b/c_0) & Eq((-c_0**2*s9b_momentum_FLYBY**2 + s9b_turn_FLYBY**2)/s9b_turn_FLYBY**2, 0) & (2/s9b_turn_FLYBY - 2*(-c_0**2*s9b_momentum_FLYBY**2 + s9b_t
34:  not a Str-keyed mapping in one copy; baseline head/len: Tuple 4 | K7 head/len: Tuple 4 | K7 value: (VELOCITY_FIXED_BY_RESTRICTION, Eq(s9b_V(s9b_r), 0), (ATTACHED_DOMAIN, Ne(ConditionSet((s9b_turn_FLYBY, s9b_momentum_FLYBY), Eq(s9b_momentum_FLYBY, b/c_0) & Eq((-c_0**2*s9b_momentum_FLYBY**2 + s9b_turn_FLYBY**2)/s9b_turn_FLYBY**2, 0) & (2/s9b_turn_FLYBY - 2*(-c_0**2*s9b_momentum_FLYBY**2 + s9b_t
76:  not a Str-keyed mapping in one copy; baseline head/len: Tuple 4 | K7 head/len: Tuple 4 | K7 value: (VELOCITY_FIXED_BY_RESTRICTION, Eq(s9b_V(s9b_r), 0), (ATTACHED_DOMAIN, Contains((s9b_f(s9b_r) + 1)**s, Reals) & ((s9b_f(s9b_r) + 1)**s > 0) & Ne(ConditionSet((s9b_turn_FLYBY, s9b_momentum_FLYBY), Eq(s9b_momentum_FLYBY, b/c_0) & Eq((-c_0**2*s9b_momentum_FLYBY**2*(s9b_f(s9b_turn_FLYBY) + 1)**(2*s)
78:  not a Str-keyed mapping in one copy; baseline head/len: Tuple 4 | K7 head/len: Tuple 4 | K7 value: (VELOCITY_FIXED_BY_RESTRICTION, Eq(s9b_V(s9b_r), 0), (ATTACHED_DOMAIN, Contains((s9b_f(s9b_r) + 1)**s, Reals) & ((s9b_f(s9b_r) + 1)**s > 0) & Ne(ConditionSet((s9b_turn_FLYBY, s9b_momentum_FLYBY), Eq(s9b_momentum_FLYBY, b/c_0) & Eq((-c_0**2*s9b_momentum_FLYBY**2*(s9b_f(s9b_turn_FLYBY) + 1)**(2*s)
82:entry 3 label: ATTACHED_DOMAIN | baseline and K7 entry 3 structurally equal: True
84:ConditionSet nodes: 12 | Integral nodes: 6
```

```
$ head -12 _scratch/s9b_build/s9b_amend12_review_r1_grok_evidence/wl_absent_reading.stdout
MISSING_INPUTFORM Missing["KeyAbsent", "absent"]
MISSING_HEAD Missing
MISSING_NUMERIC False
MISSING_FINITE_EXACT False
CE_HEAD ConditionalExpression
CE_ASSOCIATION False
CE_LIST False
CE_EQUAL False
CE_PIECEWISE False
CE_ARG_COUNT 2
```

```
$ grep -n 'JET_GETITEM\|EMITTED_CONTAINS_ZERO_GRADE\|TUPLE_KEY_READING' _scratch/s9b_build/s9b_amend12_review_r1_grok_evidence/engine_accessors.stdout
2:JET_GETITEM True
4:EMITTED_CONTAINS_ZERO_GRADE False
5:TUPLE_KEY_READING TypeError
```

```
$ sed -n '10,18p' _scratch/s9b_build/s9b_amend12_review_r1_grok_evidence/extract_nested.stdout
    IntegralKernels: ['100', '010', '020', '001']
LINE K7_WL_S9B_B_CONDITION_DEFLECTION:
  Baseline
    case ['Domain', 'GM', 'ProfileEquation', 'ProfileEquationRadius', 'InverseWeight', 'ReferenceInverse', 'VelocityReduction']
    VelocityReduction ['Status', 'SquaredVelocityEquation', 'LinearizationResidual', 'GeneralSquaredVelocity', 'GeneralConstantDomain', 'SquaredVelocity', 'SquaredVelocityResidual', 'HomogeneousCountingRatioLimit', 'ExcludedBranches', 'NoRealVelocityDomain', 'Cases', 'ZeroSetJoining', 'RemainingInputs']
  Corrupted
    case ['Domain', 'GM', 'VelocityReduction']
    VelocityReduction ['Status', 'FreeInputs', 'Constraint', 'Cases']
LINE K9b_WL_LOCAL_S9B_C_RESPONSES:
```

```
$ grep -n 'no time limit\|No time limit' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r1.md | head -3
253:    - There are no time limits (`AGENTS.md`).
```

## The density-sign parts and rho_0 (Claude 1, Grok 2)
```
$ head -8 _scratch/s9b_build/s9b_amend12_review_r1_claude_evidence/out/py_density_range.stdout
PY_S9B_C_POWER_CONDITION_DEFLECTION
entry 0: LOCAL_BULK_DENSITY_PARTS
entry 1 (label, predicate): (RHO_POSITIVE, rho_0*(s9b_f(s9b_r) + 1) > 0)
entry 2 (label, predicate): (RHO_NONPOSITIVE, rho_0*(s9b_f(s9b_r) + 1) <= 0)
entry 3 label: OBJECT_ON_BOTH_PARTS
object entry 0: RADIAL_DIFFERENTIAL_IDENTITY | entry 1 (radius): s9b_r | entry 2: (ATTACHED_DOMAIN, Interval.open(s9b_b_min, oo))
label predicate free symbols: ['rho_0', 's9b_r'] | contains the ranged radius symbol: True
label density: rho_0*(s9b_f(s9b_r) + 1)
```

```
$ grep -n 'asymptotic bulk' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
223:    `ρ₀` is the asymptotic bulk number density, `c_s0` is the corresponding asymptotic sound speed, and Part C
```

```
$ grep -n 'SIGN_PART_CONTAINS_RESPONSE' _scratch/s9b_build/s9b_amend12_review_r1_grok_evidence/engine_accessors.stdout
9:SIGN_PART_CONTAINS_RESPONSE [('RHO_POSITIVE', False), ('RHO_NONPOSITIVE', False)]
```

## "neutralizes the mutation at its construction" (Claude 2)
```
$ cat _scratch/s9b_build/s9b_amend12_review_r1_claude_evidence/out/py_neutralize.stdout | head -12
== K8: source line `mass_density = rho` -> `mass_density = s9b_constant_density`
K8 mutated source expression (mass density) at r=2, constant=1, V=0 inputs
  baseline  : 3/2
  corrupted : 1
  baseline - corrupted : 1/2
K8 MASS_BALANCE baseline payload: Eq(j_n(s9b_r), -rho_br_profile(s9b_r)*Derivative(s9b_V(s9b_r), s9b_r) - s9b_V(s9b_r)*Derivative(rho_br_profile(s9b_r), s9b_r) - 2*rho_br_profile(s9b_r)*s9b_V(s9b_r)/s9b_r)
K8 MASS_BALANCE corrupted payload: Eq(j_n(s9b_r), s9b_constant_density*(-s9b_r*Derivative(s9b_V(s9b_r), s9b_r) - 2*s9b_V(s9b_r))/s9b_r)
K8 MASS_BALANCE rhs, r8 harness profiles (V = 1/(8 r^2))
  baseline  : 1/128
  corrupted : 0
  baseline - corrupted : 1/128
K8 MASS_BALANCE rhs, V = 0 (other profiles as r8)
```

## Label and status entries (Claude 3)
```
$ grep -n -A4 'B_IMPLIED_JN_DEFLECTION' _scratch/s9b_build/s9b_amend12_review_r1_claude_evidence/out/py_labels_gates.stdout | head -8
3:B_IMPLIED_JN_DEFLECTION
4-  baseline keys: ['kind', 'definition', 'squared', 'equation', 'integrating_factor', 'general_solution', 'residual', 'constants', 'tail_domain', 'nonreal_branch']
5-  K7 keys      : ['kind', 'reason']
6-  kind baseline / K7: DIFFERENTIAL / UNDETERMINED
7-  K7 reason: Single-impact integral constraint leaves the radial V profile undetermined.
--
15:R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION
16-  baseline keys: ['kind', 'definition', 'squared', 'equation', 'integrating_factor', 'general_solution', 'residual', 'constants', 'tail_domain', 'nonreal_branch']
```

```
$ sed -n '498p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  If[quantifierMode === "FIXED_B", Return[<|"Status" -> "velocity undetermined by one fixed-b condition",
```

## Conditional values and case lists (Claude 4, Grok 3)
```
$ sed -n '611p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    return cas(((attached_domain(domain),value),(attached_domain(sp.Not(domain)),'NOT_ESTABLISHED')))
```

```
$ sed -n '271p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
gated[z_, rules_:{}] := <|"OnDomain" -> ConditionalExpression[z, rayGate /. rules],
```

```
$ head -6 _scratch/s9b_build/s9b_amend12_review_r1_grok_evidence/extract_k10_gate_carriers.stdout
TAG WL_LOCAL_S9B_BRANCH_TYPES:
  Baseline bytes=760 flags={'Not[': False, 'physicalPropagation': False, 'branchGate': False, 'propagatingLocal': False, 'ForAll': False, 'ConditionalExpression': False, 'Piecewise': True} head=Piecewise[{{{Piecewise[{{"grow
  Corrupted bytes=764 flags={'Not[': False, 'physicalPropagation': False, 'branchGate': False, 'propagatingLocal': False, 'ForAll': False, 'ConditionalExpression': False, 'Piecewise': True} head=Piecewise[{{{Piecewise[{{"grow
TAG WL_LOCAL_S9B_NONRECIPROCAL_CONSTRUCTION:
  Baseline bytes=827 flags={'Not[': False, 'physicalPropagation': False, 'branchGate': False, 'propagatingLocal': False, 'ForAll': False, 'ConditionalExpression': False, 'Piecewise': False} head=<|"OneForm" -> {-(velocityProf
  Corrupted bytes=831 flags={'Not[': False, 'physicalPropagation': False, 'branchGate': False, 'propagatingLocal': False, 'ForAll': False, 'ConditionalExpression': False, 'Piecewise': False} head=<|"OneForm" -> {-(velocityProf
```

```
$ head -12 _scratch/s9b_build/s9b_amend12_review_r1_grok_evidence/extract_case_counts.stdout
Baseline outer_case_assocs 6 vr_cases_assocs 3 vr_cases_head {<|"Domain" -> veloc
Corrupted outer_case_assocs 3 vr_cases_assocs 1 vr_cases_head {<|"Domain" -> bFixe
```

## Round 8's SymPy K1-K4c witness sites (Claude 5)
```
$ grep -n 'only baseline\|structurally equal' _scratch/s9b_build/s9b_amend12_review_r1_claude_evidence/out/py_optical_q_grades.stdout | head -6
3:   only baseline: [(0, 2, 0), (0, 2, 1), (1, 2, 0), (1, 2, 1)]  only copy: []
4:   common grade (0, 0, 0) | structurally equal: True | baseline - copy at r8 inputs: 0
5:   common grade (0, 0, 1) | structurally equal: True | baseline - copy at r8 inputs: 0
6:   common grade (1, 0, 0) | structurally equal: True | baseline - copy at r8 inputs: 0
7:   common grade (1, 0, 1) | structurally equal: True | baseline - copy at r8 inputs: 0
9:   only baseline: [(0, 0, 1), (0, 2, 1), (1, 0, 1), (1, 2, 1)]  only copy: []
```

```
$ sed -n '124p;127p;769p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
        self.c = {g: v for g, v in self.c.items() if v != 0}
        return self.c.get(g, sp.S.Zero)
    emit('OPTICAL_Q',optical_q.c)
```

