# S9b Wolfram repair round 2 measurements

Final source measurements for directive amendment 5 and spec amendment 2. These report computed objects and executable checks; they do not issue a physics verdict. Demonstration outputs remain under builder scratch paths as item 11 requires. No commits.

Builder handoff: [report](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/BUILD_REPORT.md). Source digests: [final inputs](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/final_source_hashes.json).

## Final guarded executions

All commands ran serially, without a time limit. Each exited 0 with zero stderr bytes and verified limits. Commands and working directories below are from the recorded executions.

### finished-repository

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r2-final-finished-repository -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `22.65895789093338`. Live stdout growth observations: `32`.

### finished-isolated

Working directory: `/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/isolated-final`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r2-final-finished-isolated -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/isolated-final/S9b_light_bending_mathematica_audit.wl
```

Guard wallSeconds: `22.561699615092948`. Live stdout growth observations: `32`.

### measurements

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r2-final-measurements -- math -script /var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/measurements.wl
```

Guard wallSeconds: `23.15520287491381`. Live stdout growth observations: `32`.

### finished-ablation

Working directory: `/var/projects/toy_physics-s9b-r2-wl`

```sh
python /var/projects/toy_physics-s9b-r2-wl/scripts/s11c_guarded_run.py --pool s9b --memory-gib 16 --tasks-max 64 --log-directory /var/projects/toy_physics-s9b-r2-wl/_scratch/s11c/s9b-wl-r2-final2-finished-ablation -- math -script /var/projects/toy_physics-s9b-r2-wl/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

Guard wallSeconds: `238.5078760450706`. Live stdout growth observations: `98`.

## Literal computed measurements

The signed-strata measurement uses symbolic positive `signedMagnitude`, positive radial metric/radius/radial wave covector, real local variables, and nonzero brane density. It substitutes the sign of the squared speed in the computed branch classifier; it chooses no radial profile family. The root-chart measurement uses positive squared speed and positive reference speed.

```wl
WL_LOCAL_S9B_MEASUREMENT_SIGNED_BRANCH_STRATA: <|"NegativeSquaredSpeed" -> {"decaying", "growing"}, "ZeroSquaredSpeed" -> "absent", "PositiveSquaredSpeed" -> Piecewise[{{"unable to traverse in a required direction", Global`aMetric*Global`vLocal^2 >= signedMagnitude}}, "real propagating"]|>
WL_LOCAL_S9B_MEASUREMENT_POSITIVE_ROOT_CHART: signedMagnitude > 0 && Global`c0 > 0 && deltaValue == (-Global`c0 + Sqrt[signedMagnitude])/Global`c0
```

Full, unabridged additional measurement stream: [measurements.txt](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/measurements.txt); generating script: [measurements.wl](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/measurements.wl). This includes all condition-velocity branches, integration constants/domains, response domains, substituted mass balances and the reduced gamma difference. All printed linearization, squared-velocity and substituted-balance residual values are `0`. Neither baseline nor measurements contain `Indeterminate`.

## Executable and payload checks

[Validation](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/validation.json) records nonempty byte-identical isolated/repository streams, isolation directory containing only its copied script, observed stdout growth and empty `mathematica/out/` status. Output SHA256: `ee0adba9726215abae9a606d794f07b3b5a0812d837c2ed373a04e523ecdabbd`. Final harness baseline is byte-identical to that stream.

[Payload audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/payload_checks.json) records all fourteen knife/subknife copies as exact single source mutations (apart from Text export removing the terminal newline), 373 payload rows each, and 373 dead-path rows. No baseline/corrupted payload omissions were found. Unequal list differences retain both full objects as inactive subtraction.

[K11 audit](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/k11_checks.json) records identical exact mechanical prefix extraction through the protected boundary: 10 reached tags, 363 explicitly not evaluated, no missing required fields, all reached method residuals zero, and no full-K11 execution. Coverage preserves all twelve retained grades of the nonreciprocal one-form, untruncated local dispersion/Fermat, exact arithmetic and symbolic domains; no samples, seeds or additional truncation.

Inventories: [engine tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/tags.txt), [engine NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/not_established_tags.txt), [harness tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/harness_emitted_tags.txt), [harness NOT_ESTABLISHED](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/harness_not_established_tags.txt), [K11 evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/k11_evaluated_tags.txt), [K11 not evaluated](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/k11_not_evaluated_tags.txt), [measurement tags](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/measurement_tags.txt). There are 373 engine tags (330 shared, 43 local), zero duplicates, and no missing shared name. The 270 engine NOT_ESTABLISHED occurrences are conditional outside-domain alternatives.

## Resource observations and scope

| Final run | Reported memory peak (bytes) | Maximum sampled tasks |
|---|---:|---:|
| s9b-wl-r2-final-finished-repository | 418275328 | 42 |
| s9b-wl-r2-final-finished-isolated | 417931264 | 42 |
| s9b-wl-r2-final-measurements | 418820096 | 42 |
| s9b-wl-r2-final2-finished-ablation | 2333646848 | 42 |

All recorded final memory event counters are zero; [resource observations](/var/projects/toy_physics-s9b-r2-wl/_scratch/s9b_wl_repair2/resource_summary.json). No guard refusal, guard kill, task-limit event or directive item 12 stop event occurred this round. Superseded development attempts are retained in scratch and are not used for the results above. Earlier authorized stop history is in the previous handoff. No production output, comparator, Part D, review or interpretation was created.
