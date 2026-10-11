# S9b Wolfram repair round 3 measurements

Directive amendment 6: separate computed-value changes from attached-domain changes in every evaluated harness difference and K11 method residual. The engine, knife construction sites and compact extraction boundary are unchanged. No physics verdict or production run is claimed.

[Builder handoff](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/BUILD_REPORT.md); [final input hashes](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/final_source_hashes.json); [all recorded invocations and executable checks](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/final-checks.json).

## Method

ConditionalExpression conditions and the explicit attachment keys Domain, Domains, BranchDomain, AlongRayDomain, ConstantDomain, GeneralConstantDomain, FunctionalDomain, Premise, Component and ZeroSetJoining are collected at their full structural paths. Bare computed predicates, including Closedness, ProfileEquation, Constraint, ConstraintOnOtherInputs and NoRealVelocityDomain, remain values. Conditions carried by nested Cases are visited without dropping any case. K11 coverage never enters this projection.

The CAS computes ValueDifference and DomainDifference independently. SameQ short-circuits identical parts to zero; arithmetic uses canonical subtraction; logical leaves use signed Boole subtraction. Changed categorical/structural metadata preserves both complete operands. No global Simplify/Reduce, profile sampling or further optical truncation is introduced. Complete original baseline/corrupted payloads are printed alongside the differences.

## Final guarded executions

The following runs were serial, with no wall-clock, CPU-time or inactivity deadline. Each exited 0 with zero stderr bytes and verified limits. The successful engine runs remain applicable because the engine did not change during this round. A development traversal warning affected only a superseded measurement attempt; its files are retained, and its outputs are not used below.

### finished-repository

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r3-final-finished-repository -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `22.87871941085905`. Live stdout growth observations: `32`.

### finished-isolated

Working directory: `/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/isolated-final`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r3-final-finished-isolated -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/isolated-final/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `23.30138709791936`. Live stdout growth observations: `31`.

### measurements

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r3-final2-measurements -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/measurements.wl
```

Guard wallSeconds: `23.78940407303162`. Live stdout growth observations: `34`.

### finished-ablation

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r3-final2-finished-ablation -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

Guard wallSeconds: `237.14456794899888`. Live stdout growth observations: `100`.

## Literal final comparison demonstrations

The demonstration script mechanically extracts the comparison definitions from the final harness. Its toy symbols exercise attachment/value handling without selecting physical profile values. It also computes the difference of every full engine payload with itself.

```wl
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_DOMAIN_ONLY_NUMBER: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"ConditionalExpression\"}" -> Boole[p] - Boole[q]|>|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_DOMAIN_ONLY_RELATION: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"ConditionalExpression\"}" -> Boole[p] - Boole[q]|>|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_VALUE_ONLY_PREDICATE: <|"ValueDifference" -> <|"Closedness" -> Boole[x == 0] - Boole[y == 0]|>, "DomainDifference" -> 0|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_COMPUTED_CONDITION: <|"ValueDifference" -> <|"ProfileEquation" -> Boole[x == 0] - Boole[y == 0]|>, "DomainDifference" -> 0|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_NESTED_CASE: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"Cases\", 1, \"Domain\"}" -> Boole[p] - Boole[q]|>|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_DOMAIN_ONLY_BOOLEAN: <|"ValueDifference" -> 0, "DomainDifference" -> <|"{\"ConditionalExpression\"}" -> Boole[p] - Boole[q]|>|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_VALUE_AND_DOMAIN: <|"ValueDifference" -> x - y, "DomainDifference" -> <|"{\"ConditionalExpression\"}" -> Boole[p] - Boole[q]|>|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_UNCHANGED: <|"ValueDifference" -> 0, "DomainDifference" -> 0|>
```

Identity measurements: **373**; each printed `<|"ValueDifference" -> 0, "DomainDifference" -> 0|>`. All rows, including the engine's physics demonstrations, are in the unabridged [measurement stream](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/measurements.txt) and the guarded stdout linked by the recorded invocation. Script: [measurements.wl](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/measurements.wl).

## Complete payload checks

[Mechanical audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/payload_checks.json) found no baseline/corrupted payload changes, missing fields, omitted tags or untagged diagnostic lines. Each of fourteen knives/subknives and the dead-path control printed **373 rows**. Every source copy has exactly its one manifest mutation, apart from Text export dropping its terminal newline. Both compact sources are exact prefix extractions, and the full K11 copy was not executed.

All reached K11 method residuals have `ValueDifference -> 0` and `DomainDifference -> 0`; all 363 unreached rows explicitly mark both evaluations, difference and residual NOT_EVALUATED. The common coverage keeps untruncated local dispersion/Fermat and all twelve retained grades in the nonreciprocal one-form, exact arithmetic and symbolic constructor domains, with no samples, seeds or extra truncation.

Difference-component observations (zero counts refer to the printed scalar zero, not an inferred physics verdict):

| Knife | Evaluated differences | Value part zero | Domain part zero | Both zero |
|---|---:|---:|---:|---:|
| K1 | 373 | 201 | 81 | 57 |
| K2 | 373 | 229 | 185 | 136 |
| K3 | 373 | 223 | 56 | 36 |
| K4a | 373 | 269 | 343 | 269 |
| K4b | 373 | 269 | 343 | 269 |
| K4c | 373 | 232 | 185 | 139 |
| K5 | 373 | 335 | 373 | 335 |
| K6 | 373 | 319 | 350 | 319 |
| K7 | 373 | 326 | 327 | 326 |
| K8 | 373 | 271 | 218 | 189 |
| K9a | 373 | 360 | 362 | 360 |
| K9b | 373 | 361 | 362 | 361 |
| K10 | 373 | 366 | 56 | 49 |
| K11 | 10 | 4 | 9 | 4 |
| HARNESS_DEAD_PATH | 373 | 373 | 373 | 373 |

Every row's component observations: [difference_rows.json](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/difference_rows.json). Complete inventories: [engine tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/tags.txt), [engine NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/not_established_tags.txt), [harness tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/harness_emitted_tags.txt), [harness NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/harness_not_established_tags.txt), [K11 evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/k11_evaluated_tags.txt), [K11 not evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/k11_not_evaluated_tags.txt), [measurement tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/measurement_tags.txt). There are **330 shared and 43 local engine tags**, no missing shared name and no duplicate. The **270 NOT_ESTABLISHED** engine occurrences are conditional outside-domain alternatives.

## Isolation and resources

The nonempty repository and isolated engine streams are byte-identical, SHA256 `ee0adba9726215abae9a606d794f07b3b5a0812d837c2ed373a04e523ecdabbd`. Only the copied source remains in isolation. Harness baseline matches that stream byte for byte. Redirected stdout grew while kernels were alive; `git status --porcelain research/pde_ledger_v3/mathematica/out/` produced no output. Final code hashes agree with the run inputs.

| Run | Reported memory peak (bytes) | Maximum sampled tasks |
|---|---:|---:|
| finished-repository | 429559808 | 42 |
| finished-isolated | 418246656 | 41 |
| measurements | 430960640 | 42 |
| finished-ablation | 2326315008 | 42 |

Final sampled swap is zero and all recorded final memory event counters are zero; [resource evidence](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair3/resource_summary.json). No refusal, kill, task-limit or item 12 event occurred this round. Prior history remains in the previous reports; the directive additionally records the earlier full-K11 kills at 8 and 16 GiB that motivated amendment 2. No commit, review, comparator, Part D, production output or interpretation was created.
