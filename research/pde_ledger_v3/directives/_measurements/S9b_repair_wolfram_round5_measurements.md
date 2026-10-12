# S9b Wolfram repair round 5 measurements

Directive amendment 8: exact Part C response domains and density-sign labels; mathematical three-valued value/domain comparisons. No physics verdict or production run.

[Builder report](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/BUILD_REPORT.md); [final source hashes](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/final_source_hashes.json); [run/check records](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/final-checks.json); [validation](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/validation.json).

## Construction and comparison scope

`rho0 > 0` is supplied and printed. Each response uses a single exact object on both density-sign parts. Separate real nonnegative sound-speed roots define the fixed ratio, with a strictly positive asymptotic square; a zero local sound speed yields zero. Density powers use the supplied zero-exponent convention before pressure differentiation. CAS-derived derivatives retain zero-prefactor cases through Piecewise, preventing a later zero-density substitution from introducing an undefined product for a constant response. Taylor coefficients are computed only as the analytic germ at the supplied positive asymptotic anchor, not as a replacement response at a domain or knife site. Domains require the exact response to be real and positive; the density-sign labels add no predicate.

The full-copy label audit checks 24 pairs in the baseline and every K1–K10 copy, including all implied-exchange and forward objects. Each pair has both required labels and byte-identical Object fields. [Complete label audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/density_part_checks.json).

The harness projects attached domains at their complete structural paths, including conditional expressions and nested cases; computed equations/predicates remain values. Its proof method uses exact Expand-based polynomial normalization with shared indeterminates for complete non-polynomial operands, including functions, reciprocal powers and noninteger powers, numeric scaling of relations, and BooleanMinimize. Containers, quantified bodies and ImplicitRegion conditions are normalized recursively. No unknown-sign factor is divided from an inequality. A proof of equality gives zero. A proof of inequality gives CHANGED; an unresolved proof gives UNDECIDED with the full equality, residual, corrupted object and baseline object. Relations between complete non-polynomial operands are not assumed, and the procedure does not cancel unknown denominators. Rational/radical/function identities or quantifier/set equivalences beyond this proof procedure can remain undecided. This abstraction is used only for a sufficient proof of identity or universal inequivalence; a nonzero polynomial containing independent operands is not itself labelled CHANGED. The full objects participate: no object, tag, term, profile or grade is dropped, and there is no time or size cutoff. K11 uses the same procedure for differences and method residuals.

All **301 shared payloads outside Part C** are byte-identical to the prior round-4 baseline; no engine tag was added or removed. Only 28 Part C payloads and the two local supplied-input/response payloads changed. [Preservation comparison](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/preservation_checks.json) records every changed tag.

## Final serial guarded invocations

Every invocation below exited 0, with zero stderr bytes and verified limits. No wall-clock, CPU-time or inactivity deadline was added. Each job was polled until completion. All final reported output comes from the final source hashes.

### finished-repository

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r5-final4-finished-repository -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `23.511085689999163`. Live output-growth observations: `35`.

### finished-isolated

Working directory: `/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/isolated-final4`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r5-final4-finished-isolated -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/isolated-final4/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `23.369208425981924`. Live output-growth observations: `36`.

### measurements

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r5-final6-measurements -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/measurements.wl
```

Guard wallSeconds: `25.462137928931043`. Live output-growth observations: `37`.

### finished-ablation

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r5-final6-finished-ablation -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

Guard wallSeconds: `835.5647927629761`. Live output-growth observations: `583`.

## Literal final response and equality demonstrations

Full generating source: [measurements.wl](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/measurements.wl); [unabridged measurement stream](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/measurements.txt). Parameter substitutions below are demonstrations, not restrictions of shared-vocabulary profiles.

```wl
WL_LOCAL_S9B_MEASUREMENT_RESPONSE_DEFINED_VALUES: {<|"Inputs" -> {Global`rho0 -> 1, K -> 1, Global`m -> 1, Global`c0 -> 1, Global`n -> 1, Global`s -> 0, Global`fLocal -> -1}, "Responses" -> {1, 1, 1}, "Domains" -> {True, True, True}|>, <|"Inputs" -> {Global`rho0 -> 1, K -> 1, Global`m -> 1, Global`c0 -> 1, Global`n -> 3, Global`s -> 2, Global`fLocal -> -1}, "Responses" -> {1, 0, 0}, "Domains" -> {True, False, False}|>, <|"Inputs" -> {Global`rho0 -> 1, K -> 1, Global`m -> 1, Global`c0 -> 1, Global`n -> 3, Global`s -> 2, Global`fLocal -> -2}, "Responses" -> {1, 1, 1}, "Domains" -> {True, True, True}|>, <|"Inputs" -> {Global`rho0 -> 1, K -> -1, Global`m -> 1, Global`c0 -> 1, Global`n -> 3, Global`s -> 0, Global`fLocal -> -2}, "Responses" -> {1, Undefined, 1}, "Domains" -> {True, False, True}|>}
WL_LOCAL_S9B_MEASUREMENT_DENSITY_SOUND_ZERO: {<|"Exponent" -> 0, "SoundSquared" -> 0|>, <|"Exponent" -> 1, "SoundSquared" -> 1|>, <|"Exponent" -> 2, "SoundSquared" -> 0|>, <|"Exponent" -> 3, "SoundSquared" -> 0|>}
WL_LOCAL_S9B_MEASUREMENT_DENSITY_POWER_ZERO: {1, 1, 1}
WL_LOCAL_S9B_MEASUREMENT_DENSITY_POWER_DERIVATIVES: {{<|"Order" -> 1, "Exponent" -> 0, "AtZero" -> 0|>, <|"Order" -> 1, "Exponent" -> 1, "AtZero" -> 1|>, <|"Order" -> 1, "Exponent" -> 2, "AtZero" -> 0|>, <|"Order" -> 1, "Exponent" -> 3, "AtZero" -> 0|>}, {<|"Order" -> 2, "Exponent" -> 0, "AtZero" -> 0|>, <|"Order" -> 2, "Exponent" -> 1, "AtZero" -> 0|>, <|"Order" -> 2, "Exponent" -> 2, "AtZero" -> 2|>, <|"Order" -> 2, "Exponent" -> 3, "AtZero" -> 0|>}}
WL_LOCAL_S9B_MEASUREMENT_DENSITY_RESPONSE_DERIVATIVES_AT_ZERO: {{<|"ResponseIndex" -> 2, "Order" -> 1, "Value" -> 0|>, <|"ResponseIndex" -> 2, "Order" -> 2, "Value" -> 0|>}, {<|"ResponseIndex" -> 3, "Order" -> 1, "Value" -> 0|>, <|"ResponseIndex" -> 3, "Order" -> 2, "Value" -> 0|>}}
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_MATHEMATICAL_VALUE_IDENTITY: <|"ValueDifference" -> 0, "DomainDifference" -> 0|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_MATHEMATICAL_DOMAIN_IDENTITY: <|"ValueDifference" -> 0, "DomainDifference" -> 0|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_MATHEMATICAL_SET_IDENTITY: <|"ValueDifference" -> 0, "DomainDifference" -> 0|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_UNRESOLVED_FUNCTIONS: <|"ValueDifference" -> <|"Status" -> "UNDECIDED", "Equality" -> u[x] - v[x] == 0, "Residual" -> Inactive[Subtract][u[x], v[x]], "Corrupted" -> u[x], "Baseline" -> v[x]|>, "DomainDifference" -> 0|>
WL_LOCAL_S9B_MEASUREMENT_DIFFERENCE_DISTINCT_NUMBERS: <|"ValueDifference" -> <|"Status" -> "CHANGED", "Equality" -> False, "Residual" -> -1, "Corrupted" -> 1, "Baseline" -> 2|>, "DomainDifference" -> 0|>
```

The measurement stream also retains every velocity/balance, branch, gamma, ray-domain, path-dependence and payload-identity demonstration. All 374 engine-payload identities have zero value and domain differences. Density labels for all baseline ordinary, implied-exchange and forward objects are additionally measured there.

## Complete harness evidence

[Full payload/mutation audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/payload_checks.json) found no missing payloads, tags or fields, no unexpected untagged output, and exactly one manifest mutation per knife copy (apart from the Text exporter removing a terminal newline). Both compact sources are exact mechanical prefixes. The full K11 copy was not executed. Every tag has full baseline, both compact fields, difference, method residual and coverage; the 364 unreached tags explicitly print NOT_EVALUATED. The ten reached method residuals are zero in both parts. Compact coverage is unchanged: exact symbolic constructor domains, untruncated local dispersion/Fermat, all twelve one-form grades, no samples, seeds or extra truncation.

Counts below describe printed comparison results, not a physics verdict. An undecided row is separately inventoried and does not mean changed. There are **2005** rows with undecided difference/residual leaves and **144** with proved-inequality leaves; a structured row can contain both kinds.

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

Complete inventories: [engine tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/tags.txt), [conditional NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/not_established_tags.txt), [harness tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/harness_emitted_tags.txt), [harness NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/harness_not_established_tags.txt), [K11 reached](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/k11_evaluated_tags.txt), [K11 not evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/k11_not_evaluated_tags.txt), [undecided rows](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/undecided_rows.json), [proved-inequality rows](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/changed_rows.json), [component results](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/difference_rows.json), [measurement tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/measurement_tags.txt).

## Executable checks and resources

The repository and isolated streams are nonempty and byte-identical; only the copied source remains in isolation. The harness baseline is byte-identical too. The engine remains blind (no runtime file reads/imports/writes except stdout). Redirected stdout grew while each kernel was alive. `git status --porcelain research/pde_ledger_v3/mathematica/out/` is empty. Final run hashes match the code.

| Run | Reported memory peak (bytes) | Maximum sampled tasks |
|---|---:|---:|
| finished-repository | 434638848 | 57 |
| finished-isolated | 434991104 | 58 |
| measurements | 436744192 | 57 |
| finished-ablation | 3343511552 | 59 |

Final swap and memory-event counters are zero; [resource samples summary](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/resource_summary.json). No guard refusal, guard-triggered kill, task-limit or item 12 event occurred.

## Superseded development attempts

The first completed demonstrations exposed a measurement substitution using Global`K instead of System`K and an ImplicitRegion normalization defect that retained an unevaluated Part in a held variable argument. The builder explicitly stopped that defective development harness, preserving its sources and logs; exit 137 was an explicit SIGKILL after SIGTERM did not terminate the kernel, with guardReason null and zero memory events. It was not a guard-triggered resource kill. [Termination record](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/development_abort.json) and [preserved invocation/outcomes](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/development1_checks.json) record it. A second development harness was explicitly terminated to replace general radical algebra with rational algebra over complete radical operands; its [termination record](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/comparison_method_abort.json) and [preserved run](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/development4_checks.json) retain the evidence (exit 137, guardReason null, no memory event). A third was explicitly terminated to replace global fraction combination with the final polynomial identity procedure; [record](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/polynomial_method_abort.json) and [preserved run](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair5/development5_checks.json) give the same exit-137/guardReason-null/no-memory-event evidence. These were explicit comparison-method replacements, not timed cutoffs or resource-triggered kills. No output from any incomplete harness is used as a final result. Subsequent completed demonstrations exposed lost zero derivative coefficients at zero density; explicit CAS-derived zero cases corrected this and the corresponding pressure derivative. Final engine/isolation/demonstration runs and the full harness were then performed with the final sources. No automatic retries, raised resource limits, physics construction change, output narrowing, commit or review occurred.
