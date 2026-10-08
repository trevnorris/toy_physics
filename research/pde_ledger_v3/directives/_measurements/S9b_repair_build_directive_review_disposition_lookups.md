# Measurements — S9b repair build directive G2 dispositions (generated 2026-10-08 14:40)

Generator: `_scratch/s9b_build/gen/s9b_repair_build_directive_lookups.sh` (sed/grep/sha256sum only; lines cut at 330 characters). The reviewed version is the frozen copy `_scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md`.

```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
5b00f2c4c8a0f6b9723237575429c3459528fb5599f0dfe8d25c346c6e18aefe  _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_review_baseline.sha256
5b00f2c4c8a0f6b9723237575429c3459528fb5599f0dfe8d25c346c6e18aefe  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ grep -n -o 'Verdict: NEEDS REVISION' _scratch/s9b_build/s9b_repair_build_directive_review_codex_final.txt _scratch/s9b_build/s9b_repair_build_directive_review_grok.txt
_scratch/s9b_build/s9b_repair_build_directive_review_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_repair_build_directive_review_grok.txt:1:Verdict: NEEDS REVISION
```

## C1 — K6 changes a coefficient only
```
$ grep -n 'K6, radar source' _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
121:    - **K6, radar source:** the object from which the `ln(1/b²)` coefficient is extracted is replaced by one one-way
```

```
$ sed -n '183,190p' docs/development_pipeline.md
  tensor/sector/index/derivative/field dependence while preserving enough type, dimension, retained grade and
  domain to still run; a sign flip or scalar rescale *alone* is a coefficient test. A **coefficient** knife is
  allowed only when the list names it as a channel a FORM knife cannot see (the N6 Φ-coefficient precedent).
  Two-route objects: corrupt one route; the other must stay. [`per-tooth-ablation`, `xform-prefix-is-not-a-form-control`]
- **Scope = the claim-dependency frontier, not "every line."** Cover every independent physics-bearing
  construction whose output enters the comparator or supports a claim/control — each independent construction
  path gets a genuine FORM knife (N6 was five knives, not forty tags). Mechanical/arithmetic checks get their
  appropriate perturbation, never a relabeled FORM. Imported objects are marked explicitly: a consumer harness
```

```
$ grep -n 'K6_LOG_COEFFICIENT_RATIO\|DECAYING_ENDPOINT_LOG_LIMIT' _scratch/s9b_build/s9b_repair_build_directive_review_codex_evidence/physics_review_v2.stdout.txt _scratch/s9b_build/s9b_repair_build_directive_review_codex_evidence/knife_review.stdout.txt
_scratch/s9b_build/s9b_repair_build_directive_review_codex_evidence/physics_review_v2.stdout.txt:27:K6_LOG_COEFFICIENT_RATIO: 1/2
_scratch/s9b_build/s9b_repair_build_directive_review_codex_evidence/knife_review.stdout.txt:10:DECAYING_ENDPOINT_LOG_LIMIT: 0
```

```
$ grep -n 'K5, return leg' _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
120:    - **K5, return leg:** the round trip's return leg propagates in the outgoing leg's direction.
```

## C2 — no knife reaches the path-dependence construction
```
$ grep -n 'Path dependence (D6 B4)' _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
60:8. **Path dependence (D6 B4).** Whether the nonreciprocal part depends on the path or only on the endpoints is
```

```
$ grep -n '_CURL' _scratch/s9b_build/s9b_repair_build_directive_review_codex_evidence/knife_review.stdout.txt
2:SUPPLIED_RADIAL_CURL: [0, 0, 0]
3:K1_RADIAL_CURL: [0, 0, 0]
4:K2_RADIAL_CURL: [0, 0, 0]
5:K3_RADIAL_CURL: [0, 0, 0]
6:ARBITRARY_RADIAL_COEFFICIENT_CURL: [0, 0, 0]
8:ADVECTION_FIELD_FORM_CURL_GRADE: [-2*eta*x*z/(c0**2*(x**2 + y**2 + z**2)**2), -2*eta*y*z/(c0**2*(x**2 + y**2 + z**2)**2), -2*eta*z**2/(c0**2*(x**2 + y**2 + z**2)**2)]
```

## G1 — the forward case fixes the flux sign and omits the odd observables
```
$ grep -n 'steady inward brane mass flux' _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
69:        a new live symbol, the steady inward brane mass flux. Its relation to `GM` is open.
```

```
$ grep -n 'Print the deflection, the round-trip' _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
72:      - Print the deflection, the round-trip excess time and its `ln(1/b²)` coefficient, both effective `γ`s and
```

```
$ grep -n 'expected value, sign or limit' research/pde_ledger_v3/directives/S9b_repair_decision_list.md
114:- Any expected value, sign or limit. This includes the relation of `c_comp` to `c₀`, of `j_n` to zero, and the size
```

```
$ grep -n 'DIV_JN0_C\|DIV_JN0_MINUS_C\|NONRECIPROCAL_HALF_DIFF_DENSITY\|CHORD_INTEGRAL_C_PLUS_MINUS_C\|DEFLECTION_BORN_SYMMETRIC' _scratch/s9b_build/s9b_repair_build_directive_review_grok_evidence/item10_forward_case.stdout.txt
5:DIV_JN0_C 0
6:DIV_JN0_MINUS_C 0
9:NONRECIPROCAL_HALF_DIFF_DENSITY -C*z/(c_0**2*rho_0*(b**2 + z**2)**(3/2))
12:CHORD_INTEGRAL_C_PLUS_MINUS_C 0
16:DEFLECTION_BORN_SYMMETRIC 0
```

## G2 — item 12's stop conditions catch item 10
```
$ grep -n 'a sub-problem the spec does not name\|a premise the spec does not supply' _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
99:    - a sub-problem the spec does not name;
100:    - a premise the spec does not supply.
```

## G3 — K4 names `V_r`, a Part D symbol
```
$ grep -n 'K4a–c' _scratch/s9b_build/S9b_repair_build_directive_reviewed_v0.md
118:    - **K4a–c, profile freezes:** each of `δ`, `V_r` and `ξ_w` is replaced by a constant where it is differentiated
```

```
$ grep -n 'V_r' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
254:r ≡ |x| > 0 ,      V^i(x) = V_r(r) x^i/r ,
258:`V_r`, `ρ_br`, `j_n` and all unsupplied responses remain general live profiles/actions. Eulerian
462:The full inherited derivative keys are `ProfileDerivative` of `V_r`, `delta`, `f`, `h`, `j_n`,
464:`sqrt(x1²+x2²+x3²)` (O2-R §4). The names denote `V_r`, `δ`, `f`, `h`, `j_n`, `μ_⊥`, `ρ_br` and
```

```
$ sed -n '132p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
  Within these bounds, `δ`, `V` and `ξ_w` stay free radial profiles in Parts A–C. Part B's comparison
```

