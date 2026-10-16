# Measurements — amendment 11 review, round 4 (generated 2026-10-16 06:21)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r4_lookups.sh` (sha256sum/grep/sed/cat/wc only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md` are the reviewed round-4 text's (saved before launch); the engines' are the committed files' at
`9af9f705`. The legs' stdout is quoted from their evidence directories.

## The reviewed round-4 text, its scope, and the baseline after both legs
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r4.sha256
166c108deaece96d98d959c82a623ffbcf2d35b41d700a0b36b74bc3466db81a  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
166c108deaece96d98d959c82a623ffbcf2d35b41d700a0b36b74bc3466db81a  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
```

```
$ git diff -U0 a65f752c -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep '^@@'
@@ -281,3 +281,4 @@ comparator, the production runs and the record.
@@ -288,2 +289,2 @@ comparator, the production runs and the record.
@@ -291,6 +292,15 @@ comparator, the production runs and the record.
@@ -301,2 +311,4 @@ comparator, the production runs and the record.
@@ -305,5 +317,4 @@ comparator, the production runs and the record.
@@ -343,5 +354,7 @@ comparator, the production runs and the record.
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r4.sha256 | grep -c ': OK$'
18
```

```
$ git diff --stat 9af9f705 -- research/pde_ledger_v3/scripts research/pde_ledger_v3/mathematica | wc -l
0
```

## Verdicts
```
$ grep -m1 'Verdict' _scratch/s9b_build/s9b_amend11_review_r4_claude.md
**Verdict:** NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r4_claude.md
3
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r4_grok.txt
Verdict: NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r4_grok.txt
3
```

## Claude 1: ordered cases
```
$ sed -n '292,295p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
      Attached conditions are entries stating where the part's value holds: domains, gates and premises, at every
      containing payload level and in its object labels. A branch label qualifies only its own branch. Descriptive
      entries such as `Profiles` and `Measure`, and computed values such as `Closedness`, are not conditions. Print
      which entries were evaluated as conditions.
```

```
$ sed -n '178,183p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
branchTypes = Piecewise[{
  {Map[Function[root, Piecewise[{{"growing", Im[root] > 0},
       {"decaying", Im[root] < 0}}, "absent"]], frequencyRoots], frequencyRadicand < 0},
  {"absent", relativeNorm == 0 || Det[metric2] <= 0 || rhoBr[r] == 0},
  {"unable to traverse in a required direction", branchGateLocal && Not[subcritical]}},
  "real propagating"];
```

```
$ grep '^A[1-4] ' _scratch/s9b_build/s9b_amend11_review_r4_claude_evidence/witness_masking.stdout
A1 branch 1 condition base | K10 {True, True}
A1 branch 2 condition base | K10 {False, False}
A1 branch 3 condition base | K10 {False, True}
A2 third-branch condition base Element[cSquared, Reals] && cSquared > 0 && aMetric*r^2 > 0 && cSquared <= aMetric*vLocal^2
A2 third-branch condition K10  !(Element[cSquared, Reals] && cSquared > 0 && aMetric*r^2 > 0) && cSquared <= aMetric*vLocal^2
A3 emitted Piecewise at point base {"decaying", "growing"}
A3 emitted Piecewise at point K10 {"decaying", "growing"}
A3 SameQ of the two evaluated objects True
A4 shared LocalClassification profile "muPerp(x)=-1, rhoBr(x)=1, V(x)=0, xi(x)=0; r=1, kr=1, kphi=0, c0=1"
A4 branch 1 condition base | K10 {True, True}
A4 branch 2 condition base | K10 {False, False}
A4 branch 3 condition base | K10 {False, True}
A4 LocalClassification at point base | K10 | SameQ {{"decaying", "growing"}, {"decaying", "growing"}, True}
```

## Claude 2: symbols each copy's construction fixes
```
$ sed -n '297,298p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
      Print a concrete real point assigning every symbol free in either copy or its attached conditions after
      accounting for binders; never substitute bound variables. Print a concrete admissible profile for every live
```

```
$ sed -n '157,158p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
   - Every along-ray gate quantifies over a constructed domain: the radii the flyby and both legs of the round trip
     traverse, each turning radius fixed by that ray's own construction. This covers both branch-existence
```

```
$ sed -n '252,260p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
radarEndpointEquation = Total[fi[radarAngleIntegrand,{r,rTurnRadar,#}] & /@ endpointRadii] == angularSeparation;
turningCondition[j_, t_] := (rayRadius /. r -> t) == j && t > bFar &&
  (D[rayRadius, r] /. r -> t) > 0;
flybyRadii = r >= rTurnFlyby;
radarRadii = rTurnRadar <= r <= Max[endpointRadii];
rayConstruction = radarEndpointEquation && turningCondition[b, rTurnFlyby] &&
  turningCondition[impactRadar, rTurnRadar] && rTurnRadar <= Min[endpointRadii] &&
  Inactive[ForAll][r, r >= rTurnFlyby, D[rayRadius, r] > 0] &&
  Inactive[ForAll][r, radarRadii, D[rayRadius, r] > 0];
```

```
$ grep '^B[0-4] ' _scratch/s9b_build/s9b_amend11_review_r4_claude_evidence/witness_masking.stdout
B0 ForAll binders in baseline rayGate {r}
B0 turning/impact symbols in baseline rayGate {impactRadar, rTurnFlyby, rTurnRadar}
B0 turningCondition[b,rTurnFlyby] is a conjunct of rayConstruction (base | K3) {True, True}
B1 profile and parameters {"delta(x)=1/x, V(x)=0, xi(x)=0", {c0 -> 1, b -> 3, bFar -> 1}}
B1 rayRadius base | K3 {r^2/(1 + r), r}
B2 turning radius from own construction base | K3 {r == (3 + Sqrt[21])/2, r == 3}
B3 turningCondition at shared rTurnFlyby = base root: base | K3 {True, False}
B3 turningCondition at K3 own root: K3 True
B4 K3 rayGate with rTurnFlyby = base root (other symbols unassigned) False
```

```
$ grep -m1 'B_IMPLIED_JN_DEFLECTION' research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out | grep -o 's9b_velocity_zero_radius' | wc -l
76
```

## Claude 3: the range of a reduced condition in Wolfram
```
$ sed -n '303,304p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
      worded condition, also print that mathematical statement. A set attached to a relation over a variable is
      that variable's range; decide the relation over the whole range (item 6). The exterior half-line is
```

```
$ sed -n '966p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
            condition=cas(('RADIAL_DIFFERENTIAL_IDENTITY',r,attached_domain(sp.Interval.open(bmin,sp.oo)),attached_domain(domain),eq(GM,local_mass)))
```

```
$ sed -n '400,404p;420,423p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
regularityDomain = <|"Radius" -> r > bFar, "Anchor" -> rAnchor > bFar,
  "Profiles" -> {deltaProfile,velocityProfile,xiProfile,rhoBr,fProfile},
  "FunctionalDomain" -> {"required derivatives exist", "displayed improper integrals converge",
    "Abel inverse exists on the exterior half-line"},
  "OpticalCounting" -> {"delta O(1/r)", "V/c0 O(1/Sqrt[r])", "xiPrimeSquared O(1/r)"}|>;
  <|"Domain" -> regularityDomain, "BranchDomain" -> (rayGate /. usedGateRules),
    "Cases" -> {<|"Domain" -> extra && r > bFar && rAnchor > bFar,
      "GM" -> (rule /. r -> rAnchor), "ProfileEquation" -> equation,
      "ProfileEquationRadius" -> r, "InverseWeight" -> weight,
```

```
$ grep '^P[12] ' _scratch/s9b_build/s9b_amend11_review_r4_claude_evidence/range_reading.stdout
P1 delta=(r-2)/r**2 | identities: Eq(1/r - (r - 2)/r**2, 0) | Eq((r - 2)/r**2, 0)
P1 delta=(r-2)/r**2 | pointwise at r=2: baseline False | corrupted True | domain True
P1 delta=(r-2)/r**2 | every-r counterexample sets: baseline Interval.open(1, oo) | corrupted Union(Interval.open(1, 2), Interval.open(2, oo)) | statements: False False
P2 delta=1/r | identities: True | False
P2 delta=1/r | pointwise at r=2: baseline True | corrupted False | domain True
P2 delta=1/r | every-r counterexample sets: baseline EmptySet | corrupted Interval.open(1, oo) | statements: True False
```

## Grok 1: the density-sign labels
```
$ sed -n '292,293p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
      Attached conditions are entries stating where the part's value holds: domains, gates and premises, at every
      containing payload level and in its object labels. A branch label qualifies only its own branch. Descriptive
```

```
$ sed -n '90,93p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    return cas(('LOCAL_BULK_DENSITY_PARTS',
                ('RHO_POSITIVE',sp.Gt(density,0)),
                ('RHO_NONPOSITIVE',sp.Le(density,0)),
                ('OBJECT_ON_BOTH_PARTS',value)))
```

```
$ sed -n '616,623p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
(* Labels partition the supplied physical density, even when a knife replaces
   a density at a response site. They are not added to any response domain. *)
densityParts[i_, object_] := If[i == 1, object, <|
  "SuppliedBulkDensityPremise" -> (rho0 > 0),
  "DensityParts" -> {<|"Label" -> "rho > 0", "LabelCondition" -> (rho0 (1+fProfile[r]) > 0),
      "Object" -> object|>,
    <|"Label" -> "rho <= 0", "LabelCondition" -> (rho0 (1+fProfile[r]) <= 0),
      "Object" -> object|>}|>];
```

```
$ grep '^DENSITY\|^RESPONSE\|^SYMPY\|^WOLFRAM\|^POINT' _scratch/s9b_build/s9b_amend11_review_r4_grok_evidence/witness_tests.stdout
DENSITY both labels: False
DENSITY solveset of both: EmptySet
DENSITY f=0 pos/nonpos: True False
DENSITY f=-2 pos/nonpos: False True
RESPONSE n=1,f=-2 (rho<0): density=-rho_0 square=K/m defined=True value=c_0 pos=False nonpos=True
RESPONSE n=1,f=0 (rho>0): density=rho_0 square=K/m defined=True value=c_0 pos=True nonpos=False
RESPONSE n=2,f=0 (rho>0): density=rho_0 square=2*K*rho_0/m defined=True value=c_0 pos=True nonpos=False
RESPONSE n=2,f=-2 (rho<0): density=-rho_0 square=-2*K*rho_0/m defined=False value=NOT_A_VALUE pos=False nonpos=True
SYMPY containing conditions: [('entry', 'RHO_POSITIVE'), ('entry', 'RHO_NONPOSITIVE')]
WOLFRAM positive-part containing conditions: [('LabelCondition', rho_0*(s9b_f + 1) > 0)]
POINT f=0 response value exists: True
POINT f=0 accepted if both sibling labels are conditions: False
POINT f=0 accepted if only the positive label is a condition: True
POINT f=-2 n=1 value exists while positive label is false: True
```

## Grok 2: the worded Wolfram domain
```
$ sed -n '300,306p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
      computation establishing each. All conditions follow one test: an exact decision of their mathematical
      content for that copy's object as defined by this directive and the spec, over their whole stated ranges,
      with quantifiers active and bound variables bound. Inactive or unevaluated forms decide nothing. For a
      worded condition, also print that mathematical statement. A set attached to a relation over a variable is
      that variable's range; decide the relation over the whole range (item 6). The exterior half-line is
      `(b_far, ∞)`, using the engine's cutoff (`s9b_b_min` or `bFar`). Neither `OpticalCounting` nor the spec's
      order counting is evaluated.
```

```
$ sed -n '390,395p;409p;269p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  "FlatConnectedBranch" -> Thread[markers -> 0], "FarZone" -> b > bFar|>;
abelInverse[ref_] := Module[{transformed, primitive},
  transformed = (ref /. b -> Sqrt[t])/(2 Sqrt[t]);
  primitive = Integrate[transformed/Sqrt[t-u], {t,u,Infinity},
    Assumptions -> u > 0 && c0 > 0 && Element[{GM,gamma},Reals], GenerateConditions -> False];
  FullSimplify[(-2 u D[primitive,u]/abelKernel) /. u -> r^2, r > 0 && c0 > 0]];
quantifierMode = "EVERY_B"; (* K7 *)
  If[quantifierMode === "FIXED_B",
```

```
$ grep '^ABEL' _scratch/s9b_build/s9b_amend11_review_r4_grok_evidence/worded_domain.stdout
ABEL GenerateConditions False: InputForm[0]
ABEL with u>0: InputForm[Integrate[1/Sqrt[t - u], {t, u, Infinity}, Assumptions -> u > 0]]
ABEL limit form: InputForm[Infinity]
```

```
$ grep '^HALFLINE\|^IMPROPER\|^MIXED' _scratch/s9b_build/s9b_amend11_review_r4_grok_evidence/witness_tests.stdout
HALFLINE r>bFar at point: False
HALFLINE b>bFar at point: True
IMPROPER partial antiderivative: 2*sqrt(R + 1) - 2
IMPROPER limit R->oo: oo
MIXED formula r>bFar at r=2: True
MIXED ordinary-convergence of displayed 1/sqrt(r): False
MIXED integrand-finite reading: True
MIXED optical-counting evaluated: not evaluated
```

## Grok 3: approximate numbers
```
$ sed -n '281,284p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r4.md
      **identical** exactly when they are the same expression as emitted, under a named literal structural
      comparison, with no evaluation, simplification or canonicalization first. Approximate numbers must be exactly
      equal; their printing preserves every distinction. Otherwise mark it **not compared**. Marks claim only
      structural identity; any claim of change requires a witness. Claim no mathematical equality. Print no verdict.
```

```
$ grep '^SAMEQ\|^SUBTRACT\|^INPUTFORM' _scratch/s9b_build/s9b_amend11_review_r4_grok_evidence/worded_domain.stdout
SAMEQ machine step: True
SUBTRACT: InputForm[2.220446049250313*^-16]
INPUTFORM equal strings: False
```

