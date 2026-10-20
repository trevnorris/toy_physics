# Measurements: S12 comparator production output (generated 2026-10-20 17:05)

Generator: `_scratch/s12/gen/s12_comparator_production_record.sh` (sha256sum/cat/cmp/grep/head/git only; this file
is written by it). The comparator is the version reviewed in scoped round 2 and accepted by the user with the
engine-row note (`S12_comparator_build_review_r2_disposition.md`). The production output is that round's reviewed
run `r2_comparison_01` itself, not a rerun. Filing: `_scratch/s12/gen/s12_comparator_accept_file.sh` (cp only).

## The run

```
$ python3 -c "import json;d=json.load(open('_scratch/s11c/s12_cmp_guard/r2_comparison_01/invocation.json'));print(d['command'])"
['python3', '/var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py', '--py', '/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S12_order_source.out', '--wl', '/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S12_order_source.out', '--output', '/var/projects/toy_physics/_scratch/s12/cmp_work/r2_comparison_01.jsonl']
```

```
$ cat _scratch/s11c/s12_cmp_guard/r2_comparison_01/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 204.1622238169948,
  "unit": "s11c-guard-ea21413e930f",
  "stderrBytes": 26298,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 204.09423501399579,
    "childPid": 2082699
  }
}
```

```
$ python3 -c "import json;d=json.load(open('_scratch/s11c/s12_cmp_guard/r2_tests_01/invocation.json'));print(d['command'])"
['python3', '/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py', '--output', '/var/projects/toy_physics/_scratch/s12/cmp_work/r2_tests_01.out']
```

```
$ cat _scratch/s11c/s12_cmp_guard/r2_tests_01/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 24.5910155439924,
  "unit": "s11c-guard-ae42e8fba5d4",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 24.525772786990274,
    "childPid": 2082628
  }
}
```

```
$ grep -n '^Ran\|^OK\|^FAILED' _scratch/s12/cmp_work/r2_tests_01.out
41:Ran 38 tests in 24.102s
43:OK
```

## The committed comparator is the reviewed one

```
$ grep -F -e 'scripts/S12_cross_engine_comparator.py' -e 'scripts/test_S12_cross_engine_comparator.py' _scratch/s12/cmp/review/s12_cmp_build_baseline_r2.sha256 | sed 's|/var/projects/toy_physics/||'
4ed863f4b56ff865272c781015f86faf3d347b6d46ed592d661bb86b3512f462  research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py
b222dd0469f25da8adccfca581472fe958a3feb17df04e9c38da47c6069f602b  research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py
```

```
$ sha256sum research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py
4ed863f4b56ff865272c781015f86faf3d347b6d46ed592d661bb86b3512f462  research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py
b222dd0469f25da8adccfca581472fe958a3feb17df04e9c38da47c6069f602b  research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py
```

```
$ git diff --stat 72607874 -- research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py && echo 'no change since the reviewed commit 72607874'
no change since the reviewed commit 72607874
```

## Inputs

```
$ head -1 _scratch/s12/cmp_work/r2_comparison_01.jsonl | grep -o '"input_sha256"[^}]*}'
"input_sha256":{"py":"988149e6eb2162c78329eb6cbfa90762bf505fea448c716f7de2f33686036a97","wl":"ec0d1acdd0666456cb50ce785b0c5188208b55e88f19b339535a97bcaf7c9055"}
```

```
$ sha256sum research/pde_ledger_v3/scripts/out/S12_order_source.out research/pde_ledger_v3/mathematica/out/S12_order_source.out
988149e6eb2162c78329eb6cbfa90762bf505fea448c716f7de2f33686036a97  research/pde_ledger_v3/scripts/out/S12_order_source.out
ec0d1acdd0666456cb50ce785b0c5188208b55e88f19b339535a97bcaf7c9055  research/pde_ledger_v3/mathematica/out/S12_order_source.out
```

## Output filed, and its identity with the reviewed run

```
$ grep -F 'r2_comparison_01.jsonl' _scratch/s12/cmp/review/s12_cmp_build_baseline_r2.sha256 | sed 's|/var/projects/toy_physics/||'
91d4540952112a227a2b46cc8934c068f3e97941ae2e2c479c7f873b54326cf6  _scratch/s12/cmp_work/r2_comparison_01.jsonl
```

```
$ sha256sum _scratch/s12/cmp_work/r2_comparison_01.jsonl research/pde_ledger_v3/scripts/out/S12_cross_engine_comparator.out
91d4540952112a227a2b46cc8934c068f3e97941ae2e2c479c7f873b54326cf6  _scratch/s12/cmp_work/r2_comparison_01.jsonl
91d4540952112a227a2b46cc8934c068f3e97941ae2e2c479c7f873b54326cf6  research/pde_ledger_v3/scripts/out/S12_cross_engine_comparator.out
```

```
$ cmp _scratch/s12/cmp_work/r2_comparison_01.jsonl research/pde_ledger_v3/scripts/out/S12_cross_engine_comparator.out && echo COMPARISON_IDENTICAL_TO_REVIEWED_RUN
COMPARISON_IDENTICAL_TO_REVIEWED_RUN
```
