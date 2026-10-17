# Measurements — amendment 12 review, round 2 (generated 2026-10-16 18:31)

Generator: `_scratch/s9b_build/gen/s9b_amend12_r2_lookups.sh` (sha256sum/grep/sed/cat/git show only). Line
numbers in `_scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r2.md` are the reviewed round-2 text's (saved before launch); the engines' and harnesses' are the
committed files' at `86097a44`.

## The reviewed text and the baseline
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r2.sha256 | cut -c1-64
0cca29c88f6a8e6d0e2f691a9f263335e877288744e33c80b376cb65b45dd1fd
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r2.md | cut -c1-64
0cca29c88f6a8e6d0e2f691a9f263335e877288744e33c80b376cb65b45dd1fd
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend12_review_baseline_r2.sha256 | grep -c ': OK$'
12
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend12_review_prompt_r2.md | cut -c1-64
fafd9c2557b0f0f570b82f1223bc8882dda49ae3d209255e19cd9b795bec0d9b
```

## Verdicts
```
$ head -1 _scratch/s9b_build/s9b_amend12_review_r2_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend12_review_r2_grok.txt
Verdict: NEEDS REVISION
```

## The round-2 clauses
```
$ sed -n '292,317p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r2.md
      Enumerate all direct and transitive site parts, including engine-local tags and both item 9 signs.
      Take the union across copies (K11's compact pair), without filtering by structural marks. Print the set and
      all exclusions with reasons. Print every reading and its computation: parts versus containers; case lists as
      ordered first-match or branch sets; conditional values; location identities across copies; absent entries.
      Use only engine-defined readings, identical for the same mathematical object in both engines; never use
      storage positions. Evaluate one-copy parts and their conditions; print **not emitted** for the other.
      Cite source lines defining absent-entry readings; invent none. Absence alone is no witness.

      Include reached status, label and other text with no mathematical reading as site parts. Replace their
      mathematical witness requirements with literal comparison: print both texts and location. A literal difference
      is a witness in both engines.

      Otherwise count a part tried only when all witness requirements except difference hold. Decide values and
      conditions by engine CAS; substitution alone and unresolved evaluations do not count. Print each copy's
      CAS-evaluated output of the construction defined by the knife list at trial inputs; exclude inputs giving the
      same output.

      Try item 9's shared object on each sign even if stored once; label trials do not try it. Count sign trials
      only with printed density evaluations in that sign in each emitting copy throughout any ranged radius.
      Respect the spec's profile definitions, including ρ₀ as asymptotic bulk density. Print a CAS decision and unmet
      requirement for an unreachable sign part.

      End each search after finitely many attempts with exactly one outcome: a witness; otherwise, after attempting
      every part, **none found** if all count as tried, or **undecided**, listing each untried part and unmet
      requirement. Print each attempt's tag, location, role, inputs, evaluations or **not emitted**, whether tried, and
      unmet requirements.
```

```
$ sed -n '321,323p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend12_r2.md
      for the default. A branch label qualifies only its own branch. Item 9's density-sign labels partition the
      shared object; they add no attached condition. A witness may lie on either part, including `ρ ≤ 0`; print
      its point's part, or each part's subrange where the radius is ranged. Descriptive entries (`Profiles`,
```

## Claude 1: an untried part with no attempt
```
$ sed -n '217,225p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
  If[kind==="ray",
    (* Read the full candidate and fixing equations. The finite search stops
       at the unresolved full ray boundary problem; no guessed turning radius,
       impact, endpoint cutoff or inactive-domain decision is used. *)
    entry=<|"Result"->"none found","Site"->location,"Attempt"-><|
      "CandidateParts"->{aa,bb},"ProfilesTried"->Normal[profiles],
      "ConstructionEquationsAfterProfiles"->{a["State"]["RayConstruction"],b["State"]["RayConstruction"]}/.Normal[profiles],
      "UnassignedDependentQuantities"->{Global`rTurnFlyby,Global`rTurnRadar,Global`impactRadar},
      "Reason"->"no exact solution of the full endpoint/turning construction and all ranged conditions was obtained; no scalar point or bound-variable substitution was made"|>|>;
```

## Claude 2: name lists
```
$ sed -n '333,334p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    progress('EXPORT_CLOSURE')
    emit('EXPORT_CLOSURE', sorted(own), local=True)
```

```
$ sed -n '1085,1091p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    declarations={}
    for value in objects.values():
        for symbol in value.atoms(sp.Symbol):
            declarations[symbol.name]=(symbol,'COORDINATE' if symbol in (r,x,b,B,bmin) else 'KNOB')
        for function in value.atoms(AppliedUndef):
            classification='DERIVED' if function.func.__name__=='s9b_velocity_squared' else 'PREMISE'
            declarations[function.func.__name__]=(function.func,classification)
```

```
$ sed -n '41,42p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
    'K7': ('fixed_impact = None', "fixed_impact = sp.Symbol('s9b_fixed_b',positive=True)"),
    'K8': ('mass_density = rho', "mass_density = sp.Symbol('s9b_constant_density',positive=True)"),
```

```
$ grep -n 'EXPORT_CLOSURE' _scratch/s9b_build/s9b_amend12_review_r2_claude_evidence/out/py_text_tags.stdout | head -8
1:K1 EXPORT_CLOSURE mark Str('not compared') | baseline Str count 84 | copy Str count 76 | only baseline ['+', '-', 'Phi', 's9b_velocity_integration_constant', 's9b_velocity_integration_radius', 's9b_velocity_join_radius', 's9b_velocity_squared', 's9b_velocity_zero_radius'] | only copy []
7:K2 EXPORT_CLOSURE mark Str('not compared') | baseline Str count 84 | copy Str count 83 | only baseline ['s9b_xi'] | only copy []
13:K3 EXPORT_CLOSURE mark Str('identical') | baseline Str count 84 | copy Str count 84 | only baseline [] | only copy []
19:K4a EXPORT_CLOSURE mark Str('identical') | baseline Str count 84 | copy Str count 84 | only baseline [] | only copy []
25:K4b EXPORT_CLOSURE mark Str('not compared') | baseline Str count 84 | copy Str count 82 | only baseline ['s9b_velocity_integration_constant', 's9b_velocity_integration_radius'] | only copy []
31:K4c EXPORT_CLOSURE mark Str('not compared') | baseline Str count 84 | copy Str count 83 | only baseline ['s9b_xi'] | only copy []
37:K5 EXPORT_CLOSURE mark Str('identical') | baseline Str count 84 | copy Str count 84 | only baseline [] | only copy []
43:K6 EXPORT_CLOSURE mark Str('identical') | baseline Str count 84 | copy Str count 84 | only baseline [] | only copy []
```

```
$ cat _scratch/s9b_build/s9b_amend12_review_r2_claude_evidence/out/closure_names_vs_value.stdout | head -6
log pair
  operand A: Add(Mul(Integer(-1), log(Symbol('s9b_constant_density', positive=True))), log(Mul(Symbol('s9b_constant_density', positive=True), Symbol('s9b_r', positive=True))))
  operand B: log(Symbol('s9b_r', positive=True))
  referenced-atom names A: ['s9b_constant_density', 's9b_r']
  referenced-atom names B: ['s9b_r']
  name lists literally equal: False
```

## Grok 1: the density-sign trial over a ranged radius
```
$ sed -n '966p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
            condition=cas(('RADIAL_DIFFERENTIAL_IDENTITY',r,attached_domain(sp.Interval.open(bmin,sp.oo)),attached_domain(domain),eq(GM,local_mass)))
```

```
$ sed -n '400p;421p;620,622p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
regularityDomain = <|"Radius" -> r > bFar, "Anchor" -> rAnchor > bFar,
    "Cases" -> {<|"Domain" -> extra && r > bFar && rAnchor > bFar,
  "DensityParts" -> {<|"Label" -> "rho > 0", "LabelCondition" -> (rho0 (1+fProfile[r]) > 0),
      "Object" -> object|>,
    <|"Label" -> "rho <= 0", "LabelCondition" -> (rho0 (1+fProfile[r]) <= 0),
```

```
$ sed -n '223p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
    `ρ₀` is the asymptotic bulk number density, `c_s0` is the corresponding asymptotic sound speed, and Part C
```

```
$ cat _scratch/s9b_build/s9b_amend12_review_r2_grok_evidence/density_sign_trial.stdout
LIMIT rho_0 EQUALS_RHO0 True
AT_WELL -rho_0 SIMPLIFIES_TO_NEG_RHO0 True
WELL_OVER_RHO0 -1
EDGE_OVER_RHO0_AT_SQRT_LOG2 0
OUTSIDE_TWICE_HALF_OVER_RHO0 7/8 POSITIVE True
SOLVESET_RATIO_LE_0 ConditionSet(r, 1 - 2*exp(-(r - r0)**2) <= 0, Reals)
EPSILON_GAP rho_0/2 POSITIVE True
SUBRANGE_LEFT 4 - sqrt(log(2)) RIGHT sqrt(log(2)) + 4 LEFT_GT_BMIN True
FAR_RATIO_R0_PLUS_10 1 - 2*exp(-100) N 1.0000000000000000000
K9B_S1_BASELINE_RATIO -1
K9B_S1_CORRUPTED_RATIO 1
K9B_BASELINE_GT0 False SIMPLIFIES False
K9B_CORRUPTED_GT0 True SIMPLIFIES True
K9B_BASELINE_REAL True SIMPLIFIES True
K9B_CORRUPTED_REAL True SIMPLIFIES True
POINT_IN_RHO_LE_0 True
NOT_THROUGHOUT_RAY True
SOURCE_CONTAINS S9b_SHARED_PHYSICS.md True :: ρ₀` is the asymptotic bulk number density
SOURCE_CONTAINS S9b_light_bending_sympy_audit.py True :: sp.Interval.open(bmin,sp.oo)
SOURCE_CONTAINS S9b_light_bending_sympy_audit.py True :: ('RHO_NONPOSITIVE',sp.Le(density,0))
SOURCE_CONTAINS S9b_light_bending_sympy_audit.py True :: ('OBJECT_ON_BOTH_PARTS',value)
SOURCE_CONTAINS S9b_light_bending_sympy_audit.py True :: power_density = density  # K9b
SOURCE_CONTAINS S9b_light_bending_sympy_audit.py True :: power_response=(power_density/rho0)**s
SOURCE_CONTAINS S9b_light_bending_sympy_audit.py True :: sp.Gt(local,0)
SOURCE_CONTAINS S9b_light_bending_mathematica_audit.wl True :: Radius" -> r > bFar
SOURCE_CONTAINS S9b_light_bending_mathematica_audit.wl True :: "LabelCondition" -> (rho0 (1+fProfile[r]) <= 0)
SOURCE_CONTAINS S9b_light_bending_mathematica_audit.wl True :: extra && r > bFar && rAnchor > bFar
SOURCE_CONTAINS S9b_light_bending_mathematica_audit.wl True :: localBulkDensity = rho0 (1 + fLocal); (* K9b *)
SOURCE_CONTAINS S9b_light_bending_mathematica_audit.wl True :: c0 densityPower[localBulkDensity/rho0,s]
SOURCE_CONTAINS S9b_repair_build_directive.md True :: throughout any ranged radius
SOURCE_CONTAINS S9b_repair_build_directive.md True :: they add no attached condition
SOURCE_CONTAINS S9b_repair_build_directive.md True :: including `ρ ≤ 0`
```

## Grok 2: equal text entries
```
$ sed -n '444,446p;478p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
        return {'kind':'UNDETERMINED','reason':'Matching condition contains no V; V remains a free radial profile.'}
    if mass.has(sp.Integral):
        return {'kind':'UNDETERMINED','reason':'Single-impact integral constraint leaves the radial V profile undetermined.'}
    return {'kind':'DIFFERENTIAL','definition':eq(W,velocity**2),'squared':squared,'equation':eq(equation,0),
```

```
$ sed -n '498p;511p;544p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  If[quantifierMode === "FIXED_B", Return[<|"Status" -> "velocity undetermined by one fixed-b condition",
    "Status" -> "velocity undetermined by this optical condition",
  <|"Status" -> "velocity determined on signed components with smooth zero-set joining",
```

