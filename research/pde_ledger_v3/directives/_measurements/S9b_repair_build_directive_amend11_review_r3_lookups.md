# Measurements — amendment 11 review, round 3 (generated 2026-10-16 05:48)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r3_lookups.sh` (sha256sum/grep/sed/cat only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md` are the reviewed round-3 text's (saved before launch); the engines' and harnesses' are the committed
files' at `9af9f705`. The legs' stdout is quoted from their evidence directories.

## The reviewed round-3 text, its scope, and the baseline after both legs
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r3.sha256
7d3283a488b81cc50f551411e5fad24e4a9e4b1f875cea5de5181489170f66b8  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md
7d3283a488b81cc50f551411e5fad24e4a9e4b1f875cea5de5181489170f66b8  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md
```

```
$ git diff -U0 4b4db191 -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep '^@@'
@@ -36,2 +36,2 @@ Dispositions: `directives/_measurements/S9b_SHARED_PHYSICS_amend2_review_disposi
@@ -51,3 +51,2 @@ predicate, and has each row that reads the local bulk density label where it is
@@ -57,3 +56,3 @@ clear. Dispositions:
@@ -63,4 +62,4 @@ this round. Reviewed until clear. Dispositions:
@@ -68,0 +68,7 @@ engines do not change in this round. Reviewed until clear. Dispositions:
@@ -274,27 +280,23 @@ comparator, the production runs and the record.
@@ -307,2 +309 @@ comparator, the production runs and the record.
@@ -342 +343 @@ comparator, the production runs and the record.
@@ -344 +345 @@ comparator, the production runs and the record.
@@ -347 +348 @@ comparator, the production runs and the record.
@@ -400,0 +402 @@ comparator, the production runs and the record.
@@ -422,0 +425 @@ comparator, the production runs and the record.
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r3.sha256 | grep -c ': OK$'
17
```

```
$ git diff --stat 9af9f705 -- research/pde_ledger_v3/scripts research/pde_ledger_v3/mathematica | wc -l
0
```

## Verdicts
```
$ grep -m1 'Verdict' _scratch/s9b_build/s9b_amend11_review_r3_claude.md
**Verdict:** NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r3_claude.md
6
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r3_grok.txt
Verdict: NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r3_grok.txt
3
```

## Claude 1: K11's difference and residual
```
$ sed -n '343,348p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md
      - both compact evaluations and their difference, with the structural mark on that pair;
      - the method's residual: the unmutated compact evaluation minus the full baseline payload reduced by the same
        configuration, printing both operands with their structural mark.
    - Each difference and residual is printed with its coverage: the retained grades it keeps, and its truncation
      order or sample domain. If the arithmetic is not exact, each also carries its precision and an error bound.
    - A tag the method does not reach is printed as not evaluated and not compared.
```

```
$ git show 4b4db191:research/pde_ledger_v3/directives/S9b_repair_build_directive.md | sed -n '307p'
      under item 13. This rule holds for K11's compact differences and method residual too. K11's sample domain is
```

```
$ sed -n '411,412p;321p;328p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
  If[SameQ[na,nb],Return[0]];
  If[TrueQ[change["CompleteRealDecision"]] && SameQ[change["ExistenceResult"],False],Return[0]];
    "Difference" -> difference[compactCorrupted[tag], compactBaseline[tag]],
    "MethodResidual" -> difference[compactBaseline[tag], baseline[tag]], "Coverage" -> coverage|>,
```

```
$ sed -n '219p;221p;1020p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_ablation.py
                delta = compact_difference(original[tag],corrupted[tag])
                residual = compact_difference(reduced,original[tag])
    if a==b or context is sp.false:return sp.S.Zero
```

```
$ grep '^C[123] ' _scratch/s9b_build/s9b_amend11_review_r3_claude_evidence/marks_and_k11_difference.stdout
C1 SameQ[corrupted, unmutated] False
C1 harness difference[corrupted, unmutated] <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"Domain\"}" -> <|"Status" -> "CHANGED", "Equality" -> (x > 0 && x < 0) || (x <= 0 && x >= 0), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{x}, Xor[x > 0, x < 0]], "ExistenceResult" -> Tru
C2 SameQ[ce2, ce1] False
C2 harness difference[ce2, ce1] <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"ConditionalExpression\"}" -> <|"Status" -> "CHANGED", "Equality" -> (x > 0 && x < 0) || (x <= 0 && x >= 0), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{x}, Xor[x > 0, x < 0]], "ExistenceResult" -> 
C3 harness difference on the line-411 call shape {x on x>0} vs {2x on x<0} (list) <|"ValueDifference" -> 0, "DomainDifference" -> <|"{1, \"Domain\"}" -> <|"Status" -> "CHANGED", "Equality" -> (x > 0 && x < 0) || (x <= 0 && x >= 0), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{x}, Xor
```

## Claude 2: what "identical" admits
```
$ sed -n '281,283p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md
      **identical** exactly when the payloads are structurally the same expression under a structural comparison
      the harness names in its output; otherwise mark it **not compared**. These marks claim only structural
      identity. Make no other claim of equality or change for a tag. Print no verdict.
```

```
$ sed -n '320,321p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
  na = canonical[a]; nb = canonical[b];
  If[SameQ[na,nb],Return[0]];
```

```
$ grep '^A1 \|^A4 .*k=1:\|^B1 ' _scratch/s9b_build/s9b_amend11_review_r3_claude_evidence/marks_and_k11_difference.stdout
A1 SameQ[1.0, 1.0+2^-51] False
A1 difference 4.440892098500626*^-16
A1 InputForm strings equal False
A1 strings {"1.", "1.0000000000000004"}
A4 machine 1.0 vs 1.0+k*2^-52, k=1: SameQ, strings equal, difference {True, False, 2.220446049250313*^-16}
B1 SameQ[e1,e2] False
B1 SameQ[Together[e1],e2] True
B1 e1 at x=1 Indeterminate
B1 e2 at x=1 2
```

```
$ grep '^A1 \|^A2 ' _scratch/s9b_build/s9b_amend11_review_r3_claude_evidence/witness_text_cases.stdout
A1 a = (x**2 - 1)/(x - 1) | b = x + 1
A1 a==b (sympy structural): False
A1 srepr(a)==srepr(b): False
A1 cancel(a)==cancel(b) (normalizing comparison): True
A1 a at x=1: nan | b at x=1: 2
A2 str(f1), str(f2): 1.00000000000000 1.00000000000000 | str equal: True
A2 srepr(f1), srepr(f2): Float('1.0', precision=53) Float('1.0000000000000002', precision=53)
A2 f1==f2: False | f2-f1: 2.22044604925031e-16
```

## Claude 3: a domain witness under a containing condition
```
$ sed -n '298,301p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md
      For a value witness, every attached condition is computed true; both values are finite in the ordinary sense,
      with the required integrals, limits and derivatives existing ordinarily, without regularization, principal
      value or analytic continuation. The values differ exactly, or by more than a printed bound covering both
      errors. For a domain witness, one copy's condition is computed true and the other's false. Unresolved
```

```
$ grep '^B[12] ' _scratch/s9b_build/s9b_amend11_review_r3_claude_evidence/witness_text_cases.stdout
B1 effective baseline set: Interval.open(0, oo) | effective corrupted set: Interval.open(0, oo)
B1 effective sets equal: True
B1 part sets: Db Interval.open(0, oo) | Dc Union(Interval.open(-oo, -5), Interval.open(0, oo)) | equal: False
B2 at x=-6: Db False | Dc True | C False | C&Db False | C&Dc False
```

## Claude 4: which entries are conditions; a set-valued domain; an identity's variable
```
$ sed -n '91,92p;966p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
                ('RHO_POSITIVE',sp.Gt(density,0)),
                ('RHO_NONPOSITIVE',sp.Le(density,0)),
            condition=cas(('RADIAL_DIFFERENTIAL_IDENTITY',r,attached_domain(sp.Interval.open(bmin,sp.oo)),attached_domain(domain),eq(GM,local_mass)))
```

```
$ grep '^C[12] \|^D[12] ' _scratch/s9b_build/s9b_amend11_review_r3_claude_evidence/witness_text_cases.stdout
C1 RHO_POSITIVE label: rho_0*(s9b_f(s9b_r) + 1) > 0
C1 RHO_NONPOSITIVE label: rho_0*(s9b_f(s9b_r) + 1) <= 0
C2 both labels true, as a set in y = f(r): EmptySet
D1 pointwise at r=2, GM=3/2: baseline True | corrupted False
D2 every-r reading, GM=3/2, baseline : counterexample set Union(Interval.open(1, 2), Interval.open(2, oo)) | statement true: False
D2 every-r reading, GM=3/2, corrupted : counterexample set Union(Interval.open(1, 4), Interval.open(4, oo)) | statement true: False
```

## Claude 5: j_n
```
$ sed -n '80p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  `ρ_br(x)` and `j_n(x)` are live radial profiles. `j_n` is the brane's normal exchange with the bulk and is
```

```
$ sed -n '1050p;1070p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    emit('F_MASS_SOLUTION',(attached_domain(eq(jn,0)),eq(Phi,flux),eq(liveV,forwardV)))
    labels={tag:('FORWARD_NO_FAR_ZONE_LOSS',attached_domain(eq(jn,0)),'Phi is the live outward coordinate-sphere mass flux')
```

```
$ sed -n '659,661p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
forwardPremise = <|"NormalExchange" -> (jn[r] == 0),"Flux" -> fluxInput,
  "Measure" -> "coordinate d3x"|>;
forwardWrap[z_] := <|"Premise" -> forwardPremise,"Object" -> z|>;
```

```
$ sed -n '303,304p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md

      An admissible profile assigns one specific radial function, with fixed coefficients, to each live profile
```

## Claude 6: a heuristic bound; a sampled quantifier
```
$ grep '^E2 \|^F[12] ' _scratch/s9b_build/s9b_amend11_review_r3_claude_evidence/witness_text_cases.stdout
E2 |v1-v2| = 0.00017725 | est1+est2 = 1.0e-36 | |v1-v2| > est1+est2: True
F1 sampled on 64 rational radii in [1, 10]: all true: True
F2 exact counterexample set on (1,10): {pi} | statement true: False
```

## Grok 1: one test for "computed true"
```
$ sed -n '259,265p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
  Inactive[ForAll][r, r >= rTurnFlyby, D[rayRadius, r] > 0] &&
  Inactive[ForAll][r, radarRadii, D[rayRadius, r] > 0];
traversedRadii = flybyRadii || radarRadii;
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
  physicalPropagation && speedChart && rhoBr[r] != 0 &&
  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];
traversalGate = rayConstruction && Inactive[ForAll][r, traversedRadii, physicalTraversal];
```

```
$ sed -n '81,83p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
def everywhere(variable, condition, domain):
    # An empty counterexample set is universal quantification, not sampling.
    return eq(sp.ConditionSet(variable, sp.Not(condition), domain), sp.EmptySet)
```

```
$ sed -n '400,404p;200,204p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
emit["LOCAL_NONRECIPROCAL_CONSTRUCTION", <|"OneForm" -> oddOneForm,
  "ExteriorDerivative" -> Simplify[oddExterior],
  "Closedness" -> oddClosedness,
  "Domain" -> (r > 0 && physicalTraversal && physicalPropagation && speedChart),
  "RadialPrimitive" -> Inactive[Integrate][oddOneForm[[1]], {r,rAnchor,rr}]|>];
regularityDomain = <|"Radius" -> r > bFar, "Anchor" -> rAnchor > bFar,
  "Profiles" -> {deltaProfile,velocityProfile,xiProfile,rhoBr,fProfile},
  "FunctionalDomain" -> {"required derivatives exist", "displayed improper integrals converge",
    "Abel inverse exists on the exterior half-line"},
  "OpticalCounting" -> {"delta O(1/r)", "V/c0 O(1/Sqrt[r])", "xiPrimeSquared O(1/r)"}|>;
```

```
$ cat _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/trueq_forall.wl.stdout
GATE InputForm[Inactive[ForAll][r, r > 1, r > 0]]
TRUEQ False
TRUEQ_ACTIVATE False
RESOLVE True
TRUEQ_RESOLVE True
```

```
$ sed -n '1,8p;12p' _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test_c.stdout
=== FALSE QUANTIFIED STATEMENT ===
STMT Eq(ConditionSet(s9b_r, s9b_r >= 10, Interval.open(1, oo)), EmptySet)
DOIT_EQ_TRUE False
DOIT_EQ_FALSE False
SOLVESET_ON_DOMAIN Interval(10, oo)
READING_B_HOLDS False
BINDER_SUBSTITUTED True
BINDER_DOIT True
=== CONTINUATION ===
```

```
$ grep 'FORMULA_TRUEQ\|WHOLE_ASSOC\|COUNTING_EQUALS\|CLOSEDNESS\|WITNESS_IF\|SECOND_FROM\|K4A_ZERO' _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test.wl.stdout
FORMULA_TRUEQ True
WHOLE_ASSOC_TRUEQ False
COUNTING_EQUALS_SPEC_FORMULA False
=== CLOSEDNESS SIBLING ===
WITNESS_IF_CLOSEDNESS_ATTACHED False
WITNESS_IF_CLOSEDNESS_IS_VALUE True
SECOND_FROM_RIGHT 1/5
SECOND_FROM_LEFT 0
K4A_ZERO_DERIVATIVE_EXISTS True
```

## Grok 2: the point and bound variables
```
$ cat _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/binders.wl.stdout
=== INTEGRATE BINDER ===
CASES {r, r}
SUBSTITUTED InputForm[Inactive[Integrate][1/4, {2, 1, Infinity}]]
ACTIVATE_INTACT 1
ACTIVATE_SUB FAILED
=== CONDITIONSET-LIKE FORALL RANGE SYMBOL ===
CASES {r, r, bFar, r}
RESOLVE_INTACT_BFAR True
RESOLVE_R_AND_BFAR True
TRUEQ_R_AND_BFAR False
```

```
$ grep 'CASES \|RESOLVE_ACTIVATE\|COLLAPSED' _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test.wl.stdout
CASES {r, r, rTurn, r, bFar, r}
RESOLVE_ACTIVATE False
COLLAPSED InputForm[Inactive[ForAll][2, True, True]]
TRUEQ_COLLAPSED False
RESOLVE_COLLAPSED True
```

```
$ grep 'BINDER' _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test.stdout _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test_c.stdout
_scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test.stdout:BINDER_ASSIGNED True
_scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test.stdout:BINDER_ASSIGNED_DOIT True
_scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test_c.stdout:BINDER_SUBSTITUTED True
_scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test_c.stdout:BINDER_DOIT True
_scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test_c.stdout:BINDER_EQ_TRUE True
```

## Grok 3: the profile sentence
```
$ sed -n '305,307p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r3.md
      the part carries (`δ`, `V`, `ξ_w`, and, where they enter, `ρ_br` and `f`). It is defined, with every
      derivative the part takes, at the witness point. It need not meet the condition being compared, nor the
      spec's order counting, which the computation does not test. A witness's point and profile are chosen for
```

```
$ grep -A4 'VALUE OUTSIDE ATTACHED DOMAIN' _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test.stdout
=== VALUE OUTSIDE ATTACHED DOMAIN ===
POINT 1 VALUES 1 1 FINITE True True DIFFER False
  DOMAIN_X_POSITIVE True
POINT -1 VALUES 1 -1 FINITE True True DIFFER True
  DOMAIN_X_POSITIVE False
```

```
$ grep -A1 'K10 GATE' _scratch/s9b_build/s9b_amend11_review_r3_grok_evidence/condition_test.wl.stdout
=== K10 GATE CANNOT HOLD FOR BOTH COPIES ===
BOTH False
```

