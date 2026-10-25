# S15 comparator run record — complete comparison

The final comparator and synthetic tests completed through the guard. All 1616 bindings and 32040 samples were emitted, including all 810 formerly unresolved samples across 18 bindings and all samples of the three region-balance residuals. No binding is not compared, not formed, or below its independently parsed leaf count; no sample coverage remains unresolved. No residual target was given, and no interpretation of the residuals is made here.

Inputs: accepted engine sources/transcripts `3192123b`; S15 decision list revision 3 `91248b67`; comparator contract `e305f6d4`. Adapted S16 machinery at `f08db34a`: native readers, derivative-slot translation, leaf census and integral evaluation. Synthetic notation/mutation controls and complete parsing, binding and accounting on both accepted streams re-verified that machinery. No engine was edited or executed and no commit was made.

The final code evaluates sums by exact integral linearity, retains outer parameters during nested integration, and handles residual square roots with exact algebraic relations without choosing a root sign or changing derivative points. Sparse polynomial moments and per-binding cache clearing control memory. Eight workers preserve binding-table output order; synthetic tests compare mathematical output and accounting across one, two and eight workers and reversed submission. Original printed operands and as-computed exact residuals remain in the output. No supplied balance equation is imposed.

## Commands and outcomes

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15-comparator --all-cpus-exclusive --memory-gib 55 --log-directory /var/projects/toy_physics/_scratch/s11c/s15_cmp_guard/tests_20261025_12 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s15/cmp_work/tests_12.json
```

Tests: guard exit `0`, `limitsVerified: true`; child exit `0`, `guardReason: null`. Twenty-seven synthetic control groups completed and twenty-five comparison-disable groups rejected their reductions. Script runtime `464.8804181359883` s; guard wall `465.8017124659964` s; aggregate sampled peak memory `456359936` bytes. Output `_scratch/s15/cmp_work/tests_12.json`; sha256 `4da482fd3e66a8cd5567fb89052555253935bbc3c0af0b7534ae0efb642db415`.

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15-comparator --all-cpus-exclusive --memory-gib 55 --log-directory /var/projects/toy_physics/_scratch/s11c/s15_cmp_guard/full_20261025_07 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py --workers 8 --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_two_drain_region.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out --output /var/projects/toy_physics/_scratch/s15/cmp_work/full_07.jsonl
```

Full run: guard exit `0`, `limitsVerified: true`; child exit `0`, `guardReason: null`. Script runtime `8393.043212330987` s; guard wall `8395.158274214977` s; aggregate sampled peak memory `21581824000` bytes. Both jobs verified CPUs 0–7, the 55 GiB aggregate cap, zero swap and 32 tasks; no memory-limit or OOM events were recorded. Output `_scratch/s15/cmp_work/full_07.jsonl`; sha256 `46685b1326ec2e5d731c5518283ca845fec0154080358cf93d0579467f17ef50`.

## Tables, accounting and coverage

`RULE_NAMES`, `RULES` and `binding_contract` in `scripts/S15_cross_engine_comparator.py` define 52 rules and generate 1616 bindings. `NAME_TABLE` has 45 rows; `SLOT_TABLE` has four lexical-slot rows. The first output record prints all these tables and citations in full, before comparisons. Every binding was compared; no decision-list-named object is left unjoined. The following 45 records print all test functions and points.

Per-tag accounting includes all 2213 SymPy and 1353 Wolfram occurrences: 1352 SymPy tags bound, 861 unjoined; all 1353 Wolfram tags bound. The unjoined SymPy reasons are: 832 argument-dimension record fields check arguments against signatures/bounds, unlike Wolfram's field/object signatures; 27 empty-support coefficients have no counterpart among Wolfram's nonempty groups and their sum; one chart-domain prose line corresponds to internally used Wolfram assumptions with no printed prose; one free-parameter comparison note has no counterpart beside Wolfram's plain parameter emit. Every identity, raw payload and source-specific reason is printed in the per-tag records.

All-stream parsed/compared leaves are `5858040/5855423` (SymPy), `2649640/2649540` (Wolfram). Bound-end parsed/compared leaves are `5855423/5855423`, `2649540/2649540`; no binding has a shortfall. The remaining SymPy leaves are 2513 under unjoined tags and 104 explicitly unpaired constructor/key leaves under bound tags. Wolfram has 100 explicitly unpaired key leaves. Every such leaf and its reason appears in the output; none is silently omitted.

All 712 sampled bindings have the complete 3-form × 5-stratum × 3-point grid. Forms: mixed quadratic, cubic and quartic. Strata: generic, X2 velocity zero, both radii velocities zero, both of those zero, and X1 velocity zero, with higher derivatives retained. Points `(x_1,x_2,x_3,t)` are `(3/2,2/7,-1/4,1/6)`, `(7/4,3/8,4/11,-1/7)`, `(2,-1/3,2/5,2/9)`. Tests check the domain and shell at every point/form/stratum. The 363 integral-bearing bindings have 16335 evaluated sample pairs, all exact; no numerical quadrature or unresolved integral coverage appears in the full-run inventory. The numerical fallback remains covered by synthetic tests.

JSON/count/hash audits: `python3 _scratch/s15/cmp_work/collect_full07.py` and `python3 _scratch/s15/cmp_work/verify_final_accounting.py`; outputs `_scratch/s15/cmp_work/full_07.audit.json` and `full_07.coverage-audit.json`. They enumerate the repaired 810 samples, verify complete grids and table order, and retain the unpaired-leaf inventory. Final source pins and hashes are in `_scratch/s15/cmp_work/pins/full07/sha256.json`; the audited working files have those same tested/run hashes.

Earlier outputs and pins remain preserved. The orchestrator confirmed their stop of `full_05`; it was not a guard refusal or memory kill. Under the subsequent stop/preserve/test/restart authorization, Codex stopped `full_06` after detecting the remaining evaluator error. The corrected 21-binding diagnostic then evaluated all 945 samples before the final tests and complete run. No coverage was narrowed. The requested implementation, tests, complete run and run record are now complete.
