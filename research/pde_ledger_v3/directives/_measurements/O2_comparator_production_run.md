# Measurements — O2 comparator production run (generated 2026-10-07 21:45)

Generator: `_scratch/s9b_build/gen/o2_comparator_production_record.sh` (sha256sum/cat/wc/cmp only; this file is written by it). The comparator is the accepted r5 (`48c9bf33`). The run used the repository guard in pooled mode (pool `o2-comparator`, 6 GiB), from the main checkout, with no other guarded job live. A first launch (`_scratch/s11c/o2-production/comparator/`) exited 2 at start because its output directory did not exist; it produced no output.

## Commands
```
mkdir -p _scratch/s11c/o2-production/comparator-out
python3 scripts/s11c_guarded_run.py --pool o2-comparator --memory-gib 6 --log-directory _scratch/s11c/o2-production/comparator-02 -- python3 research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py --output _scratch/s11c/o2-production/comparator-out/comparison.jsonl --accounting _scratch/s11c/o2-production/comparator-out/accounting.jsonl --catalog _scratch/s11c/o2-production/comparator-out/catalog.json
cp _scratch/s11c/o2-production/comparator-out/comparison.jsonl research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out
cp _scratch/s11c/o2-production/comparator-out/accounting.jsonl research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator_accounting.jsonl
cp _scratch/s11c/o2-production/comparator-out/catalog.json research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator_catalog.json
```

## Guard outcome
```
$ cat _scratch/s11c/o2-production/comparator-02/child-outcome.json
{
  "exitCode": 0,
  "guardReason": null,
  "wallSeconds": 170.94895659200847,
  "childPid": 3727606
}
```

```
$ wc -c < _scratch/s11c/o2-production/comparator-02/stderr
0
```

```
$ cat _scratch/s11c/o2-production/comparator/child-outcome.json; cat _scratch/s11c/o2-production/comparator/stderr
{
  "exitCode": 2,
  "guardReason": null,
  "wallSeconds": 0.36689871083945036,
  "childPid": 3727379
}
operational_error: [Errno 2] No such file or directory: '_scratch/s11c/o2-production/comparator-out/catalog.json'
```

## Comparator source run (accepted r5 baseline)
```
$ sha256sum -c _scratch/s9b_build/o2_comparator_build_review_baseline_r5.sha256 2>&1 | head -3
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/ablate_O2_cross_engine_comparator.py: OK
```

## Outputs filed, and their identity with the reviewed r5 run
```
$ sha256sum _scratch/s11c/o2-production/comparator-out/comparison.jsonl research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out _scratch/s11c/o2-comparator-build-r5/comparison.jsonl
a1a736110d4884ba6ae14fe73e76bdfb3d7882c1f898909a96798764a44707b7  _scratch/s11c/o2-production/comparator-out/comparison.jsonl
a1a736110d4884ba6ae14fe73e76bdfb3d7882c1f898909a96798764a44707b7  research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator.out
a1a736110d4884ba6ae14fe73e76bdfb3d7882c1f898909a96798764a44707b7  _scratch/s11c/o2-comparator-build-r5/comparison.jsonl
```

```
$ sha256sum _scratch/s11c/o2-production/comparator-out/accounting.jsonl research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator_accounting.jsonl _scratch/s11c/o2-comparator-build-r5/accounting.jsonl
88d8d46248adeb1c1628311a35cc74f9f7d1b17392deef264392892216f1ccdb  _scratch/s11c/o2-production/comparator-out/accounting.jsonl
88d8d46248adeb1c1628311a35cc74f9f7d1b17392deef264392892216f1ccdb  research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator_accounting.jsonl
88d8d46248adeb1c1628311a35cc74f9f7d1b17392deef264392892216f1ccdb  _scratch/s11c/o2-comparator-build-r5/accounting.jsonl
```

```
$ sha256sum _scratch/s11c/o2-production/comparator-out/catalog.json research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator_catalog.json _scratch/s11c/o2-comparator-build-r5/catalog.json
ffa630b950378040bdd60e928aba6140502964d8f0cdc65d3285b8ebf465f22d  _scratch/s11c/o2-production/comparator-out/catalog.json
ffa630b950378040bdd60e928aba6140502964d8f0cdc65d3285b8ebf465f22d  research/pde_ledger_v3/scripts/out/O2_cross_engine_comparator_catalog.json
ffa630b950378040bdd60e928aba6140502964d8f0cdc65d3285b8ebf465f22d  _scratch/s11c/o2-comparator-build-r5/catalog.json
```

```
$ cmp _scratch/s11c/o2-production/comparator-out/comparison.jsonl _scratch/s11c/o2-comparator-build-r5/comparison.jsonl && echo COMPARISON_IDENTICAL_TO_REVIEWED_RUN
COMPARISON_IDENTICAL_TO_REVIEWED_RUN
```

```
$ cmp _scratch/s11c/o2-production/comparator-out/accounting.jsonl _scratch/s11c/o2-comparator-build-r5/accounting.jsonl && echo ACCOUNTING_IDENTICAL_TO_REVIEWED_RUN
ACCOUNTING_IDENTICAL_TO_REVIEWED_RUN
```

```
$ cmp _scratch/s11c/o2-production/comparator-out/catalog.json _scratch/s11c/o2-comparator-build-r5/catalog.json && echo CATALOG_IDENTICAL_TO_REVIEWED_RUN
CATALOG_IDENTICAL_TO_REVIEWED_RUN
```

