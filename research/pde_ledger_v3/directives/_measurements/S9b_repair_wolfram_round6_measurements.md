# S9b Wolfram repair round 6 measurements

Directive amendment 9: existential value/domain differences; engine unchanged. No physics verdict or production run.

[Builder report](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/BUILD_REPORT.md); [final input hashes](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/final_source_hashes.json); [run/check records](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/final-checks.json); [validation](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/validation.json).

## Comparison method and scope

The prior polynomial/Boolean normalization remains a sufficient identity proof. A pointwise equality that does not simplify to False can still differ somewhere. The harness now forms an existential obligation using the complete two operands and their common attached conditions. Negative integer powers retain denominator exclusions; noninteger real powers use the nonnegative-base chart (strictly positive for negative exponents). Real algebraic formulas are passed to exact Resolve, with no samples, numerical tolerances, deadlines or size cutoffs. A True result prints CHANGED with the complete existential condition and proof result. A False result on a complete polynomial real decision prints zero; a False result on only a sufficient chart remains UNDECIDED. Formula support is structural, never a term/size selection. Every tag and every operand remains fully printed regardless of decidability.

Computed relations/predicates remain values; domain attachments remain a separate comparison. Value leaf contexts use the applicable ancestor attachment paths from both original payloads. Domain comparisons use symmetric difference; compatible ImplicitRegion objects use their conditions on their original variables. Real polynomial quantifiers are activated only for the exact solver. Definedness conditions for powers inside quantifiers are not moved across their binders: those complete obligations remain unresolved. General profile functions, derivatives, integral functionals and other unsupported analytic objects are not replaced by independent real values for change proofs; these print UNDECIDED unless the preceding exact identity normalization settles them. This introduces no profile family or physics premise. The engine’s supplied scalar assumptions are retained.

K11's exact extraction, configuration, coverage, payloads and residual calculation are unchanged. The amended comparison applies equally to K11's compact differences and method residuals. Coverage is not treated as an attached domain.

## Final guarded invocations

One kernel at a time, verified 16 GiB/zero swap/tasks 64 and host reserve; no time limits. Each job was polled until it ended. The engine source SHA stayed unchanged, and all final comparison output uses the final harness code. Earlier development measurements and the first full harness run are preserved in final-checks.json but are not used as final results below. That full run was allowed to finish with its source unchanged. A code inspection found that absence of a witness on a chart excluding poles was not a global equality proof; the final code leaves that case undecided. Numeric constants are excluded from free-variable quantification. A subsequent full-payload audit exposed the Element alternatives representation in live engine assumptions; the classifier now accepts that equivalent form, measured explicitly with the exact engine assumptions. These comparison-code corrections changed no physics construction.

### final2-finished-ablation

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-finished-ablation -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

Literal guard result:

```json
{
  "exitCode": 0,
  "wallSeconds": 970.8333426110912,
  "unit": "s11c-guard-e63108056bdb",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 970.7315944419242,
    "childPid": 407702
  }
}
```

Complete [stdout](/var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-finished-ablation/stdout), [stderr](/var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-finished-ablation/stderr), [resource samples](/var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-finished-ablation/resource-samples.jsonl).

### final2-source-comparison-measurements

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-source-comparison-measurements -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/comparison_measurements_final2.wl
```

Literal guard result:

```json
{
  "exitCode": 0,
  "wallSeconds": 6.555086299078539,
  "unit": "s11c-guard-88a378738a21",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 6.480053741019219,
    "childPid": 407415
  }
}
```

Complete [stdout](/var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-source-comparison-measurements/stdout), [stderr](/var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-source-comparison-measurements/stderr), [resource samples](/var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r6-final2-source-comparison-measurements/resource-samples.jsonl).

## Literal final comparison demonstrations

[Generating script](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/comparison_measurements_final2.wl) extracts the comparison definitions from the final harness; it adds no physics acceptance test.

```wl
LIVE_ASSUMPTIONS: Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar
LIVE_ASSUMPTIONS_SUPPORTED: True
POLYNOMIAL_IDENTITY: <|"ValueDifference" -> 0, "DomainDifference" -> 0|>
CROSSING_VALUES: <|"ValueDifference" -> <|"Status" -> "CHANGED", "Equality" -> x == 0, "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && 2*x != x], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> True|>, "Residual" -> Inactive[Subtract][2*x, x], "Corrupted" -> 2*x, "Baseline" -> x|>, "DomainDifference" -> 0|>
RATIONAL_VALUES: <|"ValueDifference" -> <|"Status" -> "CHANGED", "Equality" -> -x + x/(-1 + x) == 0, "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && x/(-1 + x) != x && -1 + x != 0], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> False|>, "Residual" -> Inactive[Subtract][x/(-1 + x), x], "Corrupted" -> x/(-1 + x), "Baseline" -> x|>, "DomainDifference" -> 0|>
CROSSING_PREDICATES: <|"ValueDifference" -> <|"Status" -> "CHANGED", "Equality" -> (-1 + x > 0 && x > 0) || (-1 + x <= 0 && x <= 0), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && (Xor[-1 + x > 0, x > 0])], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> True|>, "Residual" -> -Boole[x > 0] + Boole[x > 1], "Corrupted" -> x > 1, "Baseline" -> x > 0|>, "DomainDifference" -> 0|>
ATTACHED_COMMON_DOMAIN: <|"ValueDifference" -> <|"Value" -> <|"Status" -> "CHANGED", "Equality" -> x == 0, "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && x > 0 && x > 0 && 2*x != x], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> True|>, "Residual" -> Inactive[Subtract][2*x, x], "Corrupted" -> 2*x, "Baseline" -> x|>|>, "DomainDifference" -> 0|>
ATTACHED_EQUAL_STRATUM: <|"ValueDifference" -> 0, "DomainDifference" -> 0|>
DISJOINT_DOMAINS: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"Domain\"}" -> <|"Status" -> "CHANGED", "Equality" -> (x > 0 && x < 0) || (x <= 0 && x >= 0), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && (Xor[x > 0, x < 0])], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> True|>, "Residual" -> Boole[x > 0] - Boole[x < 0], "Corrupted" -> x > 0, "Baseline" -> x < 0|>|>|>
DOMAIN_CHANGE: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"Domain\"}" -> <|"Status" -> "CHANGED", "Equality" -> (-1 + x > 0 && x > 0) || (-1 + x <= 0 && x <= 0), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && (Xor[-1 + x > 0, x > 0])], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> True|>, "Residual" -> -Boole[x > 0] + Boole[x > 1], "Corrupted" -> x > 1, "Baseline" -> x > 0|>|>|>
QUANTIFIED_DOMAIN_CHANGE: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"Domain\"}" -> <|"Status" -> "CHANGED", "Equality" -> (Inactive[ForAll][y, -1 - 2*x + y^2 >= 0] && Inactive[ForAll][y, -x + y^2 >= 0]) || ( !Inactive[ForAll][y, -1 - 2*x + y^2 >= 0] &&  !Inactive[ForAll][y, -x + y^2 >= 0]), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, y, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && (Xor[Inactive[ForAll][y, -1 - 2*x + y^2 >= 0], Inactive[ForAll][y, -x + y^2 >= 0]])], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> True|>, "Residual" -> Boole[Inactive[ForAll][y, y^2 >= x]] - Boole[Inactive[ForAll][y, y^2 >= 1 + 2*x]], "Corrupted" -> Inactive[ForAll][y, y^2 >= x], "Baseline" -> Inactive[ForAll][y, y^2 >= 1 + 2*x]|>|>|>
SET_DOMAIN_CHANGE: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"Domain\"}" -> <|"Status" -> "CHANGED", "Equality" -> ImplicitRegion[-1 + x > 0, {x}] - ImplicitRegion[x > 0, {x}] == 0, "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && (Xor[-1 + x > 0, x > 0])], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> True|>, "Residual" -> Inactive[Subtract][ImplicitRegion[x > 1, {x}], ImplicitRegion[x > 0, {x}]], "Corrupted" -> ImplicitRegion[x > 1, {x}], "Baseline" -> ImplicitRegion[x > 0, {x}]|>|>|>
PROFILE_OBLIGATION: <|"ValueDifference" -> <|"Status" -> "UNDECIDED", "Equality" -> x == 0, "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && x + f[x] != f[x]], "ExistenceResult" -> Missing["UNDECIDED"], "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> False, "CompleteRealDecision" -> False|>, "Residual" -> Inactive[Subtract][x + f[x], f[x]], "Corrupted" -> x + f[x], "Baseline" -> f[x]|>, "DomainDifference" -> 0|>
COINCIDENT_POLES: <|"ValueDifference" -> <|"Status" -> "CHANGED", "Equality" -> x^(-1) == 0, "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && x^(-1) != 2/x && x != 0], "ExistenceResult" -> True, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> False|>, "Residual" -> Inactive[Subtract][x^(-1), 2/x], "Corrupted" -> x^(-1), "Baseline" -> 2/x|>, "DomainDifference" -> 0|>
DOMAIN_POLE_BOUNDARY: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"Domain\"}" -> <|"Status" -> "UNDECIDED", "Equality" -> (x^(-1) < 0 && x <= 0) || (x^(-1) >= 0 && x > 0), "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && (Xor[x^(-1) < 0, x <= 0]) && x != 0], "ExistenceResult" -> False, "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> True, "CompleteRealDecision" -> False|>, "Residual" -> Boole[x^(-1) < 0] - Boole[x <= 0], "Corrupted" -> x^(-1) < 0, "Baseline" -> x <= 0|>|>|>
NUMERIC_CONSTANT: <|"ValueDifference" -> <|"Status" -> "UNDECIDED", "Equality" -> -x + Pi*x == 0, "ChangeProof" -> <|"ExistenceCondition" -> Inactive[Exists][{c0, ell, b, bFar, ZE, ZR, r, rr, impact, GM, n, K, m, s, rho0, Phi, rAnchor, x}, Element[c0 | ell | b | bFar | ZE | ZR | r | rr | impact | GM | n | K | m | s | rho0 | Phi | rAnchor, Reals] && c0 > 0 && ell > 0 && b > bFar > 0 && ZE > 0 && ZR > 0 && r > 0 && rr > 0 && impact > 0 && rho0 > 0 && m > 0 && rAnchor > bFar && Pi*x != x], "ExistenceResult" -> Missing["UNDECIDED"], "Method" -> "exact real quantifier elimination on the common defined domain", "SupportedFormula" -> False, "CompleteRealDecision" -> False|>, "Residual" -> Inactive[Subtract][Pi*x, x], "Corrupted" -> Pi*x, "Baseline" -> x|>, "DomainDifference" -> 0|>
```

## Complete harness checks

[Full payload/mutation audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/payload_checks.json) checks every printed baseline and corrupted payload against the captured streams, every source copy against exactly one manifest mutation, all comparison shapes, all compact fields and full coverage. No tags, objects, terms or grades are truncated, sampled or omitted. Every ordinary knife and the dead-path control has 374 rows. K11 has 374 rows, of which 10 are reached and 364 explicitly NOT_EVALUATED; the full K11 copy is never executed. Both compact extracts are exact prefixes; every reached method residual is zero in both parts. The dead-path control has zero value/domain differences for all tags.

There are **73** rows containing changed leaves and **2000** containing undecided leaves; a structured row may contain both. All **219 CHANGED leaves** have **219 printed True existence results**. These are comparison measurements, not physics verdicts.

| Knife | Evaluated differences | Value zero | Domain zero | Both zero |
|---|---:|---:|---:|---:|
| K1 | 374 | 202 | 81 | 57 |
| K2 | 374 | 230 | 185 | 136 |
| K3 | 374 | 224 | 56 | 36 |
| K4a | 374 | 270 | 344 | 270 |
| K4b | 374 | 270 | 344 | 270 |
| K4c | 374 | 233 | 185 | 139 |
| K5 | 374 | 336 | 374 | 336 |
| K6 | 374 | 320 | 351 | 320 |
| K7 | 374 | 327 | 328 | 327 |
| K8 | 374 | 272 | 219 | 190 |
| K9a | 374 | 361 | 362 | 361 |
| K9b | 374 | 362 | 363 | 362 |
| K10 | 374 | 367 | 56 | 49 |
| K11 | 10 | 4 | 9 | 4 |
| HARNESS_DEAD_PATH | 374 | 374 | 374 | 374 |

Complete inventories: [engine tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/tags.txt), [conditional NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/not_established_tags.txt), [harness tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/harness_emitted_tags.txt), [harness NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/harness_not_established_tags.txt), [K11 reached](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/k11_evaluated_tags.txt), [K11 not evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/k11_not_evaluated_tags.txt), [undecided rows](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/undecided_rows.json), [changed rows](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/changed_rows.json), [component results](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/difference_rows.json).

[Full density-part audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/density_part_checks.json) verifies 24 positive/nonpositive pairs with identical Object fields in each of the fourteen full engine copies. [Baseline preservation](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/baseline_preservation.json) verifies the unchanged engine hash and byte-identical round-5/round-6 baseline streams. Earlier engine/demonstration/isolation checks remain in the unchanged-engine [round-5 measurements](/var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/directives/_measurements/S9b_wl_repair_round5_measurements.md); this round reruns the required full harness and its changed comparison demonstrations.

## Resources and events

| Final run | Peak bytes | Maximum sampled tasks | Full stdout bytes |
|---|---:|---:|---:|
| final2-finished-ablation | 3774926848 | 59 | 1190311085 |
| final2-source-comparison-measurements | 217907200 | 26 | 10616 |

[Resource summary](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair6/resource_summary.json). Zero final swap and memory-event counters. No refusal, guard kill, builder kill, task-limit or item 12 event in this round. No commits or reviews.
