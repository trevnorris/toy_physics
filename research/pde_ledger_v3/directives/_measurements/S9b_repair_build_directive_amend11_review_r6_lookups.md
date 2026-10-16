# Measurements — amendment 11 review, round 6 (generated 2026-10-16 07:31)

Generator: `_scratch/s9b_build/gen/s9b_amend11_r6_lookups.sh` (sha256sum/grep/sed/cat/wc only). Line numbers in
`_scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r6.md` are the reviewed round-6 text's (saved before launch); the engines' are the committed files' at
`9af9f705`. The legs' stdout is quoted from their evidence directories.

## The reviewed round-6 text, its scope, and the baseline after both legs
```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r6.sha256
97de0851f7d2b33b099d6c624655a3d2f5e85c0321f5bbf0811dec05e710e13f  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r6.md
97de0851f7d2b33b099d6c624655a3d2f5e85c0321f5bbf0811dec05e710e13f  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r6.md
```

```
$ git diff -U0 dab6d4c5 -- research/pde_ledger_v3/directives/S9b_repair_build_directive.md | grep '^@@'
@@ -69,3 +69,4 @@ clear. Dispositions:
@@ -296,3 +297,13 @@ comparator, the production runs and the record.
@@ -301,2 +312 @@ comparator, the production runs and the record.
@@ -306,7 +316,4 @@ comparator, the production runs and the record.
@@ -314,2 +321,3 @@ comparator, the production runs and the record.
@@ -319,2 +327,2 @@ comparator, the production runs and the record.
@@ -421 +429,3 @@ comparator, the production runs and the record.
@@ -444 +454,3 @@ comparator, the production runs and the record.
```

```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend11_review_baseline_r6.sha256 | grep -c ': OK$'
20
```

```
$ git diff --stat 9af9f705 -- research/pde_ledger_v3/scripts research/pde_ledger_v3/mathematica | wc -l
0
```

## Verdicts
```
$ grep -m1 'Verdict' _scratch/s9b_build/s9b_amend11_review_r6_claude.md
**Verdict:** NEEDS REVISION
```

```
$ grep -c '^\*\*[0-9]\. ' _scratch/s9b_build/s9b_amend11_review_r6_claude.md
1
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend11_review_r6_grok.txt
Verdict: CLEAR
```

```
$ grep -c 'No findings' _scratch/s9b_build/s9b_amend11_review_r6_grok.txt
1
```

## Claude 1: pairing solutions across copies
```
$ sed -n '302,303p;321,322p' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r6.md
      conditions, identify its role: assigned at the point; ranged, with its range; bound; or fixed by that copy's
      construction. Print all equations fixing a quantity, each copy's solution set and the solution used, if any.
      For a value witness, both copies use the same independent inputs; each selected solution satisfies its copy's
      whole construction. Every attached condition is computed true; both values are ordinary finite values, with
```

```
$ sed -n '899,904p' research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
    polar,azimuth=sp.symbols('s9b_polar s9b_azimuth',real=True)
    coords=sp.Matrix([r*sp.sin(polar)*sp.cos(azimuth),r*sp.sin(polar)*sp.sin(azimuth),r*sp.cos(polar)])
    jacobian=sp.trigsimp(coords.jacobian((r,polar,azimuth)).det())
    area=sp.integrate(jacobian,(polar,0,sp.pi),(azimuth,0,2*sp.pi))
    mass_density = rho  # K8
    exchange=sp.simplify(-sp.diff(area*mass_density*velocity,r)/area)
```

```
$ grep 'EXCHANGE K8\|REBUILT_BASE\|ZERO_SOLUTION\|SQUARED_NONNEG\|same expression\|right (B)' _scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/join_limits_exact.stdout _scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/join_limits_exact.stdout:BRANCH 3 z 2 right (B) sqrt(2)/4 j_n BASE = -rho_br(z)*right = -3*sqrt(2)/4 | j_n K8 = -s9b_constant_density*right at 3 = -3*sqrt(2)/4
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/join_limits_exact.stdout:BRANCH 3 z 4 right (B) 1/16 j_n BASE = -rho_br(z)*right = -3/16 | j_n K8 = -s9b_constant_density*right at 3 = -3/16
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/join_limits_exact.stdout:BRANCH 4 z 2 right (B) -sqrt(2)/4 j_n BASE = -rho_br(z)*right = 3*sqrt(2)/4 | j_n K8 = -s9b_constant_density*right at 3 = 3*sqrt(2)/4
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/join_limits_exact.stdout:BRANCH 4 z 4 right (B) -1/16 j_n BASE = -rho_br(z)*right = 3/16 | j_n K8 = -s9b_constant_density*right at 3 = 3/16
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:EXCHANGE K8 s9b_constant_density*(-s9b_r*Derivative(s9b_V(s9b_r), s9b_r) - 2*s9b_V(s9b_r))/s9b_r
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 2 SIGNS (-1, -1) REBUILT_BASE_EQUALS_EMITTED True
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 3 SIGNS (-1, 1) REBUILT_BASE_EQUALS_EMITTED True
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 4 SIGNS (1, -1) REBUILT_BASE_EQUALS_EMITTED True
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 5 SIGNS (1, 1) REBUILT_BASE_EQUALS_EMITTED True
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:ZERO_SOLUTION_SET z>bmin {2, 4}
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:SQUARED_NONNEGATIVE_COUNTEREXAMPLES on (1,oo) EmptySet
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 2 j_n(z) BASE oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s9b_velocity_zero_radius**3 + 52*s9b_velocity_zero_radius**2 - 96*s9b_velocity_zero_radius + 64)) | K8 oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s9b
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 3 j_n(z) BASE -oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s9b_velocity_zero_radius**3 + 52*s9b_velocity_zero_radius**2 - 96*s9b_velocity_zero_radius + 64)) | K8 -oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 4 j_n(z) BASE oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s9b_velocity_zero_radius**3 + 52*s9b_velocity_zero_radius**2 - 96*s9b_velocity_zero_radius + 64)) | K8 oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s9b
_scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/solution_pairing.stdout:BRANCH 5 j_n(z) BASE -oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s9b_velocity_zero_radius**3 + 52*s9b_velocity_zero_radius**2 - 96*s9b_velocity_zero_radius + 64)) | K8 -oo*sign(sqrt(s9b_velocity_zero_radius**4 - 12*s
```

```
$ sed -n '572,574p' research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
      "AtComponentBoundary" -> <|"VelocityLimit" -> Inactive[Limit][value,r -> rBoundary],
        "ExchangeLimit" -> Inactive[Limit][jValue,r -> rBoundary],
        "Domain" -> "finite matching one-sided limits for a differentiable real velocity and flux"|>|>]]],solution["Cases"]];
```

## Observation for round 8: native one-sided limits at a double zero
```
$ grep 'native left' _scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/join_limits_exact.stdout
BRANCH 2 SIGNS (-1, -1) z 2 (A) native left,right (-sqrt(2)/4, -sqrt(2)/4) equal True | (B) explicit left,right (sqrt(2)/4, -sqrt(2)/4) equal False
BRANCH 2 SIGNS (-1, -1) z 4 (A) native left,right (-1/16, -1/16) equal True | (B) explicit left,right (1/16, -1/16) equal False
BRANCH 3 SIGNS (-1, 1) z 2 (A) native left,right (-sqrt(2)/4, sqrt(2)/4) equal False | (B) explicit left,right (sqrt(2)/4, sqrt(2)/4) equal True
BRANCH 3 SIGNS (-1, 1) z 4 (A) native left,right (-1/16, 1/16) equal False | (B) explicit left,right (1/16, 1/16) equal True
BRANCH 4 SIGNS (1, -1) z 2 (A) native left,right (sqrt(2)/4, -sqrt(2)/4) equal False | (B) explicit left,right (-sqrt(2)/4, -sqrt(2)/4) equal True
BRANCH 4 SIGNS (1, -1) z 4 (A) native left,right (1/16, -1/16) equal False | (B) explicit left,right (-1/16, -1/16) equal True
BRANCH 5 SIGNS (1, 1) z 2 (A) native left,right (sqrt(2)/4, sqrt(2)/4) equal True | (B) explicit left,right (-sqrt(2)/4, sqrt(2)/4) equal False
BRANCH 5 SIGNS (1, 1) z 4 (A) native left,right (1/16, 1/16) equal True | (B) explicit left,right (-1/16, 1/16) equal False
```

```
$ cat _scratch/s9b_build/s9b_amend11_review_r6_claude_evidence/left_limit_numeric.stdout
LEFT_LIMIT_DIRECTION - SIGNS (-1, -1)
q = 2 - 1e-3 operand 0.35417271718042048294
q = 2 - 1e-5 operand 0.35355957783837661688
q = 2 - 1e-7 operand 0.35355345246512319272
sqrt(2)/4 = 0.35355339059327376220
END
```

## Orchestrator item O1: review and record process in builder-facing text
```
$ grep -n 'round-8 review\|Round-8 review\|build review round 8' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend11_r6.md
71:amendment review round 5, printing each witness's readings for build review round 8 to check. Parts 2 and 3
308:      under this directive, the spec and the knife's definition. The record claims it only after round-8 review
429:- **Repair round 8 (amendment 11).** Item 14 as amended holds for the harness. Round-8 review checks each printed
454:- **Repair round 8 (amendment 11).** Item 14 as amended holds for the harness. Round-8 review checks each printed
```

