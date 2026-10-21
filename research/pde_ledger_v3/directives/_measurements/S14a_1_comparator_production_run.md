# Measurements: S14a part 1 comparator production output (generated 2026-10-21 14:44)

Generator: `_scratch/s14a/gen/s14a1_comparator_production_record.sh` (sha256sum/cat/cmp/grep/head/git only; this
file is written by it). The comparator is the version reviewed in build review round 1, where both legs were clear
(`S14a_1_comparator_build_review_r1_disposition.md`). The production output is that round's reviewed run
`r1_full_01` itself, not a rerun. Filing: `_scratch/s14a/gen/s14a1_comparator_accept_file.sh` (cp only).

## The run

```
$ python3 -c "import json;d=json.load(open('_scratch/s11c/s14a1_cmp_guard/r1_full_01/invocation.json'));print(d['command'])"
['python3', '/var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py', '--py', '/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out', '--wl', '/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out', '--output', '/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl']
```

```
$ cat _scratch/s11c/s14a1_cmp_guard/r1_full_01/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 106.72195675599505,
  "unit": "s11c-guard-f4b383ba54f3",
  "stderrBytes": 14173,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 106.65296896500513,
    "childPid": 2323430
  }
}
```

```
$ python3 -c "import json;d=json.load(open('_scratch/s11c/s14a1_cmp_guard/r1_tests_01/invocation.json'));print(d['command'])"
['python3', '/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py', '--output', '/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_tests_01.out']
```

```
$ cat _scratch/s11c/s14a1_cmp_guard/r1_tests_01/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 76.3144247730088,
  "unit": "s11c-guard-f7fd35b69ca6",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 76.24139743400156,
    "childPid": 2323269
  }
}
```

```
$ grep -n '^Ran\|^OK\|^FAILED' _scratch/s14a/cmp_work/r1_tests_01.out
35:Ran 32 tests in 75.627s
37:OK
```

## The committed comparator is the reviewed one

```
$ grep -F -e 'scripts/S14a_1_cross_engine_comparator.py' -e 'scripts/test_S14a_1_cross_engine_comparator.py' _scratch/s14a/cmp/review/s14a1_cmp_build_baseline_r1.sha256
d628d27a11eba4b67ac9ddcce7d73b6bd65f8325c30e7c8277f501b49da4e7c9  research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
65e96aeeab41aa362d2148e745e7cf7f2a94d5cf37846157597f7229b0217812  research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py
```

```
$ sha256sum research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py
d628d27a11eba4b67ac9ddcce7d73b6bd65f8325c30e7c8277f501b49da4e7c9  research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
65e96aeeab41aa362d2148e745e7cf7f2a94d5cf37846157597f7229b0217812  research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py
```

## Inputs

```
$ head -1 _scratch/s14a/cmp_work/r1_full_01.jsonl | grep -o '"input_sha256"[^}]*}'
"input_sha256":{"py":"3bfb17aa5205e7f0ddec1eea44c920d4f0bb4e3dffe4ccf6fcc74abe88fbe0db","wl":"7911ff19f667a388a2a54a3eb0d076ed06143d791ddb2b5bccfe397588652d5b"}
```

```
$ sha256sum research/pde_ledger_v3/scripts/out/S14a_1_projection.out research/pde_ledger_v3/mathematica/out/S14a_1_projection.out
3bfb17aa5205e7f0ddec1eea44c920d4f0bb4e3dffe4ccf6fcc74abe88fbe0db  research/pde_ledger_v3/scripts/out/S14a_1_projection.out
7911ff19f667a388a2a54a3eb0d076ed06143d791ddb2b5bccfe397588652d5b  research/pde_ledger_v3/mathematica/out/S14a_1_projection.out
```

```
$ git log -1 --format='%h %s' 94fddf2f
94fddf2f Accept the S14a part 1 build (cleared in scoped round 2)
```

## Output filed, and its identity with the reviewed run

```
$ grep -F 'r1_full_01.jsonl' _scratch/s14a/cmp/review/s14a1_cmp_build_baseline_r1.sha256
1db2009d9fa4b77e5d91fd138c43c6d6342585bb75ac56cd2d2f7c0d268c6ce2  _scratch/s14a/cmp_work/r1_full_01.jsonl
```

```
$ sha256sum _scratch/s14a/cmp_work/r1_full_01.jsonl research/pde_ledger_v3/scripts/out/S14a_1_cross_engine_comparator.out
1db2009d9fa4b77e5d91fd138c43c6d6342585bb75ac56cd2d2f7c0d268c6ce2  _scratch/s14a/cmp_work/r1_full_01.jsonl
1db2009d9fa4b77e5d91fd138c43c6d6342585bb75ac56cd2d2f7c0d268c6ce2  research/pde_ledger_v3/scripts/out/S14a_1_cross_engine_comparator.out
```

```
$ cmp _scratch/s14a/cmp_work/r1_full_01.jsonl research/pde_ledger_v3/scripts/out/S14a_1_cross_engine_comparator.out && echo COMPARISON_IDENTICAL_TO_REVIEWED_RUN
COMPARISON_IDENTICAL_TO_REVIEWED_RUN
```
