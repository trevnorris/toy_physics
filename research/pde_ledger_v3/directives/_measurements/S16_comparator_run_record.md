# S16 comparator run record

Final code: `scripts/S16_cross_engine_comparator.py` and `scripts/test_S16_cross_engine_comparator.py` (paths relative to `research/pde_ledger_v3/`). Gate: `5994e9ac`; accepted engine sources/transcripts: `1d990ee6`; scientific decision list: `0464b80a`. No commit was made. No residual target was given.

Commands, from `/var/projects/toy_physics`, in execution order:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s16_cmp_guard/tests_deliverable_05 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s16/cmp_work/tests_deliverable_05.json
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s16_cmp_guard/full_05 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S16_moving_region.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S16_moving_region.out --output /var/projects/toy_physics/_scratch/s16/cmp_work/full_05.jsonl
```

| Run | Guard outcome | Guard wall seconds | Observed cgroup peak bytes |
|---|---|---:|---:|
| Tests `tests_deliverable_05` | Admitted; `exitCode: 0`, `limitsVerified: true`, child `exitCode: 0`, `guardReason: null` | 14.437464013986755 | 76738560 |
| Full comparison `full_05` | Admitted; `exitCode: 0`, `limitsVerified: true`, child `exitCode: 0`, `guardReason: null` | 7616.337004225 | 1114296320 |

Both runs used a 53687091200-byte memory cap, zero swap, 32-task cap, one logical CPU and no wall/CPU/inactivity deadline. Resource samples show maximum swap 0 and terminal memory-limit/OOM event counters 0. Receipts and samples are in `_scratch/s11c/s16_cmp_guard/{tests_deliverable_05,full_05}/`.

- Tests: [output](/var/projects/toy_physics/_scratch/s16/cmp_work/tests_deliverable_05.json), SHA256 `84611d3ff54b0415a522e6cf35ed4abf4a9884ca91e73349a44c80bbcbbd4210`. All 13 synthetic control groups completed. Script runtime 13.905590880953241 seconds; peak process RSS 77744 KiB.
- Full comparison: [output](/var/projects/toy_physics/_scratch/s16/cmp_work/full_05.jsonl), SHA256 `8e8a36ad85d6e8746f596252dad09577f55f4a7eed8cc5471e0c126257c5ec21`. Script runtime 7610.401615736 seconds; peak process RSS 980668 KiB. [Extracted accounting summary](/var/projects/toy_physics/_scratch/s16/cmp_work/full_05_summary.json).
- Comparator SHA256: `5a4812f631a8f1db747ba6a3fc6595efbb2f00013e22af3e207a16171c08f83e`; tests SHA256: `6d15ac27d6f88d169c3555ded0c8b3ef85915e2032d9b5ee4e010c4cd0965471`. Both and the accepted input hashes were rechecked after completion. SymPy 1.14.0; mpmath 1.3.0.

The first output record contains 47 binding rules, all 2028 generated bindings and all 91 name rows with citations. All 2028 bindings were compared; none was unformed or not compared. Accounting covers 26028 Py tags / 37578 WL tags, with parsed leaves 238488 / 580543 and compared leaves 161466 / 212319. Bound-end parsed leaves are 161466 / 213171. The 426 bindings below their parsed count are `coverage:1602` and every `comparison_notes:1603` through `comparison_notes:2027`: each leaves the WL `Missing` constructor/reason metadata unpaired while comparing its message.

The output lists every unjoined tag and unpaired component, with source-cited grouped reasons. It records 39663 evaluated samples and 162 normalized zero-count exclusions across three polynomial forms, three rational points and five strata. Numerical integral fallbacks: 0; unresolved coverage: 0. Numerical fallback controls exercised the 40/60-digit quadrature/error-estimate path synthetically. AST/InputForm readers and derivative-slot machinery were adapted from S14a_1 (`de81b9fc`); native reading, component/name maps, aliases, domains, dimensions, integral evaluation, mutation controls and accounting were reverified.

The earlier `full_04` output is partial because of the user-reported reboot and is not used as the completed comparison. Its output/logs remain intact; its old reservation was reconciled only after verifying terminal service state, no queued start and a dead owner, with the old registry preserved under `_scratch/s16/cmp_work/`. `full_05` is the completed final-code run.

Open questions: interpretation of the printed residuals and non-author build review belong to the S16 record; neither was performed here.
