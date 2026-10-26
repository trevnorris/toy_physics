# Measurements: S15 comparator production output (generated 2026-10-26 00:21)

Generator: `_scratch/s15/gen/s15_comparator_production_record.sh` (sha256sum/cat/cmp/grep/head/jq/git only; this
file is written by it). The comparator is the version reviewed in rounds 0 and 1 and accepted in round 1
(`S15_comparator_build_review_r1_disposition.md`). The production output is the reviewed run `full_07` itself, not
a rerun. Filing: `_scratch/s15/gen/s15_comparator_accept_file.sh` (cp only).

## The full run

```
$ jq -c '.command' _scratch/s11c/s15_cmp_guard/full_20261025_07/invocation.json
["python3","-B","/var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py","--workers","8","--py","/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_two_drain_region.out","--wl","/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out","--output","/var/projects/toy_physics/_scratch/s15/cmp_work/full_07.jsonl"]
```

```
$ cat _scratch/s11c/s15_cmp_guard/full_20261025_07/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 8395.158274214977,
  "unit": "s11c-guard-45dc24f2a9a0",
  "stderrBytes": 3958799,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 8395.089541009016,
    "childPid": 319658
  }
}
```

## The tests (repair 1, the version accepted)

```
$ jq -c '.command' _scratch/s11c/s15_cmp_guard/tests_20261025_15/invocation.json
["python3","-B","/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py","--output","/var/projects/toy_physics/_scratch/s15/cmp_work/tests_15.json"]
```

```
$ cat _scratch/s11c/s15_cmp_guard/tests_20261025_15/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 817.9979908599926,
  "unit": "s11c-guard-ec972fcc4827",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 817.9349437160126,
    "childPid": 373210
  }
}
```

```
$ jq -c '.controls' _scratch/s15/cmp_work/tests_15.json
{"exclusive_guard_controls":"completed","nested_radial_derivative":"completed","nested_integral_parameters":"completed","split_integral_evaluation":"completed","sparse_integral_oracle":"completed","parallel_schedule":"completed","long_native_payload":"completed","basic_mutations":"completed","labels_and_conditions":"completed","equivalences":"completed","derivative_points":"completed","mixed_second_derivative_points":"completed","partial_time_semantics":"completed","integral_changes":"completed","binding_repoints":"completed","name_repoints":"completed","bound_name_repoints":"completed","pieces_groups_and_components":"completed","packaging_and_accounting":"completed","equation_relation_kind":"completed","repeated_records":"completed","supplied_equation_and_alias":"completed","coverage_controls":"completed","test_functions":"completed","duplicates":"completed","dag_sampling":"completed","quadrature":"completed"}
```

```
$ jq -c '[(.targeted_ablations|length), .synthetic_only, .runtime_seconds]' _scratch/s15/cmp_work/tests_15.json
[4800,true,816.9317231979803]
```

## The code: reviewed bytes

```
$ grep -F -e 'research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py' -e 'research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py' -e 'scripts/s11c_guarded_run.py' _scratch/s15/cmp_review/s15_cmp_build_baseline_r1.sha256
786ff802b0bcf54b1e8103f596bc7af3216eafc800a524106e23ee293fdd11aa  research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py
95cd49fd5127dcb3aa243e1793d23e6dcc5b2ff12b888bbe0d3c88a4b0eb88ae  research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py
86be9b55a56f7a3538cdf3852b74562bdc8c5bab19c9ad894cbc2c3b6bbe257b  scripts/s11c_guarded_run.py
```

```
$ sha256sum research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py scripts/s11c_guarded_run.py
786ff802b0bcf54b1e8103f596bc7af3216eafc800a524106e23ee293fdd11aa  research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py
95cd49fd5127dcb3aa243e1793d23e6dcc5b2ff12b888bbe0d3c88a4b0eb88ae  research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py
86be9b55a56f7a3538cdf3852b74562bdc8c5bab19c9ad894cbc2c3b6bbe257b  scripts/s11c_guarded_run.py
```

## The inputs: the accepted S15 transcripts

```
$ head -1 _scratch/s15/cmp_work/full_07.jsonl | jq -c '.input_sha256'
{"py":"4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c","wl":"811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6"}
```

```
$ sha256sum research/pde_ledger_v3/scripts/out/S15_two_drain_region.out research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out
4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c  research/pde_ledger_v3/scripts/out/S15_two_drain_region.out
811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6  research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out
```

```
$ git log -1 --format='%h %s' 3192123b
3192123b Accept the S15 build (review round 1 clear)
```

## The output: identical to the reviewed run

```
$ grep -F -e '_scratch/s15/cmp_work/full_07.jsonl' -e '_scratch/s15/cmp_work/full_07.audit.json' -e '_scratch/s15/cmp_work/full_07.coverage-audit.json' -e '_scratch/s15/cmp_work/tests_15.json' _scratch/s15/cmp_review/s15_cmp_build_baseline_r1.sha256
46685b1326ec2e5d731c5518283ca845fec0154080358cf93d0579467f17ef50  _scratch/s15/cmp_work/full_07.jsonl
09c022636cec61bf7f597a128c2edf8693863bcd4f31c5a22813579290f897f4  _scratch/s15/cmp_work/full_07.audit.json
b6444f9deeec7b623863cce1a2f10cf762708e4f04839a26c5769bc5b65a0c9c  _scratch/s15/cmp_work/full_07.coverage-audit.json
c0a390c8d6ede7f6c53bac2cb4103ba64773b2a93165967f6cf25dcdaf0d1077  _scratch/s15/cmp_work/tests_15.json
```

```
$ sha256sum _scratch/s15/cmp_work/full_07.jsonl research/pde_ledger_v3/scripts/out/S15_cross_engine_comparator.out _scratch/s15/cmp_work/full_07.audit.json _scratch/s15/cmp_work/full_07.coverage-audit.json _scratch/s15/cmp_work/tests_15.json
46685b1326ec2e5d731c5518283ca845fec0154080358cf93d0579467f17ef50  _scratch/s15/cmp_work/full_07.jsonl
46685b1326ec2e5d731c5518283ca845fec0154080358cf93d0579467f17ef50  research/pde_ledger_v3/scripts/out/S15_cross_engine_comparator.out
09c022636cec61bf7f597a128c2edf8693863bcd4f31c5a22813579290f897f4  _scratch/s15/cmp_work/full_07.audit.json
b6444f9deeec7b623863cce1a2f10cf762708e4f04839a26c5769bc5b65a0c9c  _scratch/s15/cmp_work/full_07.coverage-audit.json
c0a390c8d6ede7f6c53bac2cb4103ba64773b2a93165967f6cf25dcdaf0d1077  _scratch/s15/cmp_work/tests_15.json
```

```
$ cmp _scratch/s15/cmp_work/full_07.jsonl research/pde_ledger_v3/scripts/out/S15_cross_engine_comparator.out && echo COMPARISON_IDENTICAL_TO_REVIEWED_RUN
COMPARISON_IDENTICAL_TO_REVIEWED_RUN
```
