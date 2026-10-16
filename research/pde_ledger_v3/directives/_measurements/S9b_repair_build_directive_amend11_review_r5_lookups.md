# Measurements — amendment 11 review, round 5 (generated 2026-10-16 06:51)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r5_lookups.sh` (sha256sum/grep/sed/cat/wc only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r5.md` are the reviewed round-5 text's (saved before launch); the engines' are the committed files' at
`9af9f705`. The legs' stdout is quoted from their evidence directories.

## The reviewed round-5 text, its scope, and the baseline after both legs
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r5.sha256
efa14a5109b830e91687c59b90af8c3f51fd919ccb2e664db084d265a824ae73  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r5.md
efa14a5109b830e91687c59b90af8c3f51fd919ccb2e664db084d265a824ae73  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r5.md
```

```
$ git diff -U0 46d32255 -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep '^@@'
@@ -292,29 +292,35 @@ comparator, the production runs and the record.
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r5.sha256 | grep -c ': OK$'
19
```

```
$ git diff --stat 9af9f705 -- research/pde_ledger_v3/scripts research/pde_ledger_v3/mathematica | wc -l
0
```

## Verdicts
```
$ grep -m1 'Verdict' _scratch/s9b_build/s9b_amend11_review_r5_claude.md
**Verdict:** NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r5_claude.md
2
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r5_grok.txt
Verdict: NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r5_grok.txt
1
```

## Claude 1: a construction with no value or several
```
$ sed -n '300,302p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r5.md
      Print a concrete real point assigning all independent inputs free in either copy or its conditions after
      accounting for binders and stated ranges; never substitute bound variables. Each quantity fixed by a copy's
      construction takes that copy's own computed value; print it and its computation. Print a concrete admissible
```

```
$ sed -n '252,260p;348p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
radarEndpointEquation = Total[fi[radarAngleIntegrand,{r,rTurnRadar,#}] & /@ endpointRadii] == angularSeparation;
turningCondition[j_, t_] := (rayRadius /. r -> t) == j && t > bFar &&
  (D[rayRadius, r] /. r -> t) > 0;
flybyRadii = r >= rTurnFlyby;
radarRadii = rTurnRadar <= r <= Max[endpointRadii];
rayConstruction = radarEndpointEquation && turningCondition[b, rTurnFlyby] &&
  turningCondition[impactRadar, rTurnRadar] && rTurnRadar <= Min[endpointRadii] &&
  Inactive[ForAll][r, r >= rTurnFlyby, D[rayRadius, r] > 0] &&
  Inactive[ForAll][r, radarRadii, D[rayRadius, r] > 0];
  "RadarImpact" -> (impactRadar == b + (shift /. flatJetRules /. jetRules /. Thread[markers -> 1]))|>];
```

```
$ sed -n '944p;953p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
        z=sp.Symbol('s9b_velocity_zero_radius',positive=True)
                joining=sp.And(common,sp.Gt(z,bmin),eq(squared.subs(r,z),0),eq(left,right),
```

```
$ grep '^F1a \|^F1b ' _scratch/s9b_build/s9b_amend11_review_r5_claude_evidence/fixed_and_ranged.stdout
F1a profile "delta(r)=r^2, V=0, xi=0; c0=1, bFar=1, b=2"
F1a rayRadius base | K3 {r/(1 + r^2), r}
F1a turningCondition[2,t] solutions base False
F1a turningCondition[2,t] solutions K3 rev`tt == 2
F1a Exists t turningCondition base | K3 {False, True}
F1b profile "delta(r)=1/(r^2/3-3r+8)-1, V=0, xi=0; c0=1, bFar=1, b=6"
F1b rayRadius base (r*(24 - 9*r + r^2))/3
F1b turningCondition[6,t] solutions base rev`tt == 3 - Sqrt[3] || rev`tt == 3 + Sqrt[3]
F1b roots numeric {1.26794919243112270647255365849412763306`12., 4.73205080756887729352744634150587236694`12.}
F1b D[rayRadius,r] at profile 8 - 6*r + r^2
F1b engine flyby monotonicity (simplified derivative) per root {False, True}
```

```
$ grep 'ZERO_SET\|^BRANCH' _scratch/s9b_build/s9b_amend11_review_r5_claude_evidence/zero_join_and_labels.stdout
BRANCH_COUNT 4
ZERO_SET_OF_OWN_CONSTRUCTION {2, 4}
BRANCH 0 signs (-1, -1) z 2 attached_domain False j_n value 4
BRANCH 0 signs (-1, -1) z 4 attached_domain False j_n value 8
BRANCH 1 signs (-1, 1) z 2 attached_domain True j_n value -4
BRANCH 1 signs (-1, 1) z 4 attached_domain True j_n value -8
BRANCH 2 signs (1, -1) z 2 attached_domain True j_n value 4
BRANCH 2 signs (1, -1) z 4 attached_domain True j_n value 8
BRANCH 3 signs (1, 1) z 2 attached_domain False j_n value -4
BRANCH 3 signs (1, 1) z 4 attached_domain False j_n value -8
```

## Claude 2: which variables are ranges
```
$ sed -n '309,312p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r5.md
      named in that engine's symbols; existence and convergence are ordinary. Relations are decided over their
      whole stated ranges, whether carried by attached sets or domain formulas, including the radius of a reduced
      every-`b` condition in both engines (item 6). The exterior radial half-line is `s9b_r > s9b_b_min` (SymPy) or
      `r > bFar` (Wolfram).
```

```
$ sed -n '411,413p;421p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
    Return[<|"Domain" -> extra && bFixed > bFar,
      "BranchDomain" -> (rayGate /. usedGateRules),
      "Cases" -> {<|"Domain" -> extra && bFixed > bFar,
    "Cases" -> {<|"Domain" -> extra && r > bFar && rAnchor > bFar,
```

```
$ grep '^F2 \|^(delta' _scratch/s9b_build/s9b_amend11_review_r5_claude_evidence/fixed_and_ranged.stdout
F2 base top Domain <|"Radius" -> r > bFar, "Anchor" -> rAnchor > bFar, "Profiles" -> {deltaProfile, velocityProfile, xiProfile, rhoBr, fProfile}, "FunctionalDomain" -> {"required derivatives exist", "displayed improper integrals converge", "Abel inverse exists on the exterior half-line"}, "OpticalCo
F2 base case Domain r > bFar && rAnchor > bFar
F2 base ProfileEquationRadius r
F2 K7 top Domain bFixed > bFar
F2 K7 case Domain bFixed > bFar
F2 K7 case keys {"Domain", "GM", "VelocityReduction"}
F2 K7 VelocityReduction Status "velocity undetermined by one fixed-b condition"
F2 K7 Constraint head Equal
F2 profile delta(r) r^(-2)
F2 K7 Constraint at profile 4*bFixed*GM + Pi == 0
F2 GM solving the K7 constraint at bFixed=2 -1/8*Pi
F2 point reading: K7 Constraint at bFixed=2, GM=GMstar True
F2 range reading: ForAll bFixed>bFar K7 Constraint at GM=GMstar False
F2 Theta(b) at profile (engine comparisonTheta activated) -(Pi/b^2)
F2 baseline-quantifier statement ForAll b>bFar Theta(b)==ref(b) at GM=GMstar False
F2 base ProfileEquation at profile False
F2 base ProfileEquation decided over r>bFar False
F2 K7 GM entry at profile (value), bFixed unassigned {{GM -> -1/4*Pi/bFixed}}
F2 profile delta(r) r^(-1)
F2 K7 Constraint at profile 1 + 2*GM == 0
F2 GM solving the K7 constraint at bFixed=2 -1/2
F2 point reading: K7 Constraint at bFixed=2, GM=GMstar True
F2 range reading: ForAll bFixed>bFar K7 Constraint at GM=GMstar True
F2 Theta(b) at profile (engine comparisonTheta activated) -2/b
F2 baseline-quantifier statement ForAll b>bFar Theta(b)==ref(b) at GM=GMstar True
F2 base ProfileEquation at profile True
F2 base ProfileEquation decided over r>bFar True
F2 K7 GM entry at profile (value), bFixed unassigned {{GM -> -1/2}}
```

```
$ grep -i 'counterexample\|label' _scratch/s9b_build/s9b_amend11_review_r5_claude_evidence/zero_join_and_labels.stdout
RANGE_READING counterexample set of Eq(squared(z),0) over z>bmin Union(Interval.open(1, 2), Interval.open(2, 4), Interval.open(4, oo)) | decided: False
LABEL f(r)= -2 + 3/s9b_r on (b_min,oo) with b_min=1
LABEL RHO_POSITIVE set Interval.open(1, 3)
LABEL RHO_NONPOSITIVE set Interval(3, oo)
LABEL RHO_POSITIVE for every r in range Eq(ConditionSet(s9b_r, rho_0*(-1 + 3/s9b_r) <= 0, Interval.open(1, oo)), EmptySet)
LABEL RHO_NONPOSITIVE for every r in range Eq(ConditionSet(s9b_r, rho_0*(-1 + 3/s9b_r) > 0, Interval.open(1, oo)), EmptySet)
LABEL RHO_POSITIVE counterexample set over range Interval(3, oo) | decided for every r: False
LABEL RHO_NONPOSITIVE counterexample set over range Interval.open(1, 3) | decided for every r: False
```

## Grok 1: "required derivatives exist"
```
$ sed -n '306,309p;322,324p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r5.md
      All conditions require exact decisions of their mathematical content for that copy's object under this
      directive and the spec, over their whole stated ranges, with quantifiers active. Inactive or unevaluated
      forms decide nothing. For each worded condition, print its mathematical statement with variables and ranges
      named in that engine's symbols; existence and convergence are ordinary. Relations are decided over their
      An admissible profile assigns each live profile in the part or its conditions one specific radial function
      with fixed coefficients, including `δ`, `V`, `ξ_w`, `ρ_br`, `f`, `μ_⊥` and `j_n` where present, defined with
      all required derivatives at the witness point. Witness inputs are chosen for this test only, enter no engine
```

```
$ sed -n '400,404p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
regularityDomain = <|"Radius" -> r > bFar, "Anchor" -> rAnchor > bFar,
  "Profiles" -> {deltaProfile,velocityProfile,xiProfile,rhoBr,fProfile},
  "FunctionalDomain" -> {"required derivatives exist", "displayed improper integrals converge",
    "Abel inverse exists on the exterior half-line"},
  "OpticalCounting" -> {"delta O(1/r)", "V/c0 O(1/Sqrt[r])", "xiPrimeSquared O(1/r)"}|>;
```

```
$ grep -A6 'WORDED required-derivatives' _scratch/s9b_build/s9b_amend11_review_r5_grok_evidence/py-witness.stdout
== WORDED required-derivatives: point vs whole half-line
profile Abs(r-2) at r=3 left/right/exists 1 1 True
same profile at r=2 left/right/exists -1 1 False
Radius formula r>bFar at the witness point True
kink lies in r>bFar True
exists at every r>bFar False
point-reading condition true True
```

