# S9b Wolfram repair round 4 measurements

Directive amendment 7: print the path-dependence predicate in the required sense and gate it on the full constructed along-ray branch domain. The harness continues to separate value and domain differences and to evaluate K11 by identical exact prefix extraction. No physics verdict, review or production run is claimed.

[Builder handoff](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/BUILD_REPORT.md); [final input hashes](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/final_source_hashes.json); [all recorded invocations and executable checks](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/final-checks.json).

## Construction and scope

The existing derived odd Fermat one-form supplies the exterior derivative and closedness test. The shared observable now uses the CAS-computed negation of that test, passed to the standard observable gate after rayConstruction, propagationGate and traversalGate are available. The flyby turning equation uses the flyby impact parameter; the radar turning equation and finite-endpoint angle equation fix the radar construction. The quantified radius domain is the union of the flyby exterior interval and the radar interval reaching both endpoints. Restrictions and forward copies retain their existing substitutions into these constructions.

The former mixed shared payload is retained as LOCAL_NONRECIPROCAL_CONSTRUCTION, with the complete one-form, exterior derivative, closedness, local construction domain and radial primitive. This is supporting construction evidence. The shared predicate itself carries the complete along-ray gate and its NOT_ESTABLISHED complement.

The compact boundary remains immediately after the local nonreciprocal construction. Consequently K11 still evaluates the same constructor scope on both copies. The gated shared predicate now lies after that boundary and is explicitly not evaluated by this compact method. The full baseline payload is nevertheless printed for it, as for every other unreached tag. No full K11 run or additional profile/grade narrowing is used.

Mechanical preservation evidence: [preservation_checks.json](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/preservation_checks.json). The previous mixed shared payload is byte-identical to the new local construction payload. Among existing engine tags, only the repaired shared predicate and local-name registry changed; no tag was removed.

## Final guarded executions

All runs below were serial, with no wall-clock, CPU-time or inactivity deadline. Each exited 0 with zero stderr bytes and verified limits. The first demonstration attempt used a raw-syntax comparison across ConditionalExpression evaluation and a replacement that did not match its intended gate; both measurement bugs were corrected. That superseded demonstration output is retained in scratch but is not used for the results below. The unchanged engine/harness from the successful runs were not replayed.

### finished-repository

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r4-final-finished-repository -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `20.454437725013122`. Live stdout growth observations: `29`.

### finished-isolated

Working directory: `/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/isolated-final`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r4-final-finished-isolated -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/isolated-final/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `22.256632823031396`. Live stdout growth observations: `33`.

### measurements

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r4-final2-measurements -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/measurements.wl
```

Guard wallSeconds: `23.28359479410574`. Live stdout growth observations: `34`.

### finished-ablation

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r4-final-finished-ablation -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

Guard wallSeconds: `237.97000257507898`. Live stdout growth observations: `102`.

## Literal final predicate and gate demonstrations

The identities compare complete evaluated conditional objects against the separately emitted branch conditions, and compare the outside condition with their complement. The four Boolean strata exercise the engine gate with its branch-condition outcomes bound to False/True; they select no physical radial profile.

```wl
WL_LOCAL_S9B_MEASUREMENT_PATH_DEPENDENCE_PREDICATE: False
WL_LOCAL_S9B_MEASUREMENT_PATH_DEPENDENCE_GATE_IDENTITIES: <|"OnDomainUsesBothConditions" -> True, "OutsideDomainIsComplement" -> True, "PredicateEqualsNotClosedness" -> True|>
WL_LOCAL_S9B_MEASUREMENT_PATH_DEPENDENCE_GATE_STRATA: {{<|"PropagationCondition" -> False, "TraversalCondition" -> False, "Object" -> <|"OnDomain" -> Undefined, "OutsideDomain" -> "NOT_ESTABLISHED", "Supplied" -> {"DISPERSION", "ADVECTION", "OPTICAL_RATIO", "EMBEDDING", "LAB_HELD"}|>|>, <|"PropagationCondition" -> False, "TraversalCondition" -> True, "Object" -> <|"OnDomain" -> Undefined, "OutsideDomain" -> "NOT_ESTABLISHED", "Supplied" -> {"DISPERSION", "ADVECTION", "OPTICAL_RATIO", "EMBEDDING", "LAB_HELD"}|>|>}, {<|"PropagationCondition" -> True, "TraversalCondition" -> False, "Object" -> <|"OnDomain" -> Undefined, "OutsideDomain" -> "NOT_ESTABLISHED", "Supplied" -> {"DISPERSION", "ADVECTION", "OPTICAL_RATIO", "EMBEDDING", "LAB_HELD"}|>|>, <|"PropagationCondition" -> True, "TraversalCondition" -> True, "Object" -> <|"OnDomain" -> False, "OutsideDomain" -> Undefined, "Supplied" -> {"DISPERSION", "ADVECTION", "OPTICAL_RATIO", "EMBEDDING", "LAB_HELD"}|>|>}}
```

The complete nonreciprocal construction, constructed ray domains and turning/endpoint equations are printed in the unabridged [measurement stream](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/measurements.txt) and guarded stdout. The same stream retains signed branch classification, response-domain, velocity and mass-balance demonstrations. Generating script: [measurements.wl](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/measurements.wl).

## Harness differences and complete payload checks

Every evaluated difference and K11 method residual continues to contain ValueDifference and DomainDifference. Attached domains are collected at their full structural paths, including nested cases; computed relations and predicates remain values. Identical parts compute to zero; changed logical leaves use exact signed indicator subtraction; categorical changes retain both full operands. K11 coverage is separate from attached domains.

Eight retained comparison-method demonstrations and **374** engine-payload identity measurements are in the full stream. Every identity printed `<|"ValueDifference" -> 0, "DomainDifference" -> 0|>`.

[Mechanical audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/payload_checks.json) found no baseline/corrupted payload changes, missing fields, omitted tags or untagged diagnostic lines. Fourteen knives/subknives and the dead-path control each printed **374 rows**. Each copy has exactly its manifest mutation, apart from Text export dropping the terminal newline. Both compact sources are exact prefix extractions; the full K11 copy was not executed.

K11 reaches **10** local tags and explicitly marks **364** tags not evaluated. Both parts of every reached method residual are zero. Its common coverage remains untruncated local dispersion/Fermat and all twelve retained grades of the nonreciprocal one-form, exact arithmetic and symbolic constructor domains, with no samples, seeds or extra truncation.

Difference-component observations (counts of printed scalar zeros, not inferred physics verdicts):

| Knife | Evaluated differences | Value part zero | Domain part zero | Both zero |
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
| K9a | 374 | 361 | 363 | 361 |
| K9b | 374 | 362 | 363 | 362 |
| K10 | 374 | 367 | 56 | 49 |
| K11 | 10 | 4 | 9 | 4 |
| HARNESS_DEAD_PATH | 374 | 374 | 374 | 374 |

Every row's component observations: [difference_rows.json](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/difference_rows.json). Complete inventories: [engine tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/tags.txt), [engine NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/not_established_tags.txt), [harness tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/harness_emitted_tags.txt), [harness NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/harness_not_established_tags.txt), [K11 evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/k11_evaluated_tags.txt), [K11 not evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/k11_not_evaluated_tags.txt), [measurement tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/measurement_tags.txt). There are **330 shared and 44 local engine tags**, no missing shared name and no duplicate. The **271 NOT_ESTABLISHED** engine occurrences are conditional outside-domain alternatives.

## Isolation and resources

Repository and isolated streams are nonempty and byte-identical, SHA256 `824917ebd0ecfe65f6c958b37b36e10e8eb73ea7fb1e0b758edd7558b8ad1685`. Only the copied source remains in isolation. Harness baseline matches that stream byte for byte. Redirected stdout grew while kernels were alive; `git status --porcelain research/pde_ledger_v3/mathematica/out/` produced no output. Final code hashes agree with run inputs.

| Run | Reported memory peak (bytes) | Maximum sampled tasks |
|---|---:|---:|
| finished-repository | 554364928 | 41 |
| finished-isolated | 419532800 | 42 |
| measurements | 420159488 | 41 |
| finished-ablation | 2326425600 | 42 |

Final sampled swap and all recorded final memory event counters are zero; [resource evidence](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair4/resource_summary.json). No refusal, kill, task-limit or item 12 event occurred this round. Prior authorized history remains in previous reports; the directive also records the full-K11 kills at 8 and 16 GiB that motivated amendment 2. No commit, review, comparator, Part D, production output or interpretation was created.
