# Measurements: S16 comparator production output (generated 2026-10-24 10:51)

Generator: `_scratch/s16/gen/s16_comparator_production_record.sh` (sha256sum/cat/cmp/grep/head/jq/git only; this
file is written by it). The comparator is the version reviewed in the scoped round 2 and accepted
(`S16_comparator_review_r2_disposition.md`). The production output is that round's reviewed run `repair_r2_full_01`
itself, not a rerun. Filing: `_scratch/s16/gen/s16_comparator_accept_file.sh` (cp only).

## The full run

```
$ jq -c '.command' _scratch/s11c/s16_cmp_guard/repair_r2_full_01/invocation.json
["python3","-B","/var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py","--py","/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S16_moving_region.out","--wl","/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S16_moving_region.out","--output","/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_full_01.jsonl"]
```

```
$ cat _scratch/s11c/s16_cmp_guard/repair_r2_full_01/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 7705.871029508009,
  "unit": "s11c-guard-ecbad9f604ac",
  "stderrBytes": 179893,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 7705.8000143239915,
    "childPid": 153337
  }
}
```

## The tests

```
$ jq -c '.command' _scratch/s11c/s16_cmp_guard/repair_r2_tests_01/invocation.json
["python3","-B","/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py","--output","/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_tests_01.json"]
```

```
$ cat _scratch/s11c/s16_cmp_guard/repair_r2_tests_01/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 20.440217527997447,
  "unit": "s11c-guard-369bb95b82cb",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 20.373167177996947,
    "childPid": 153243
  }
}
```

```
$ cat _scratch/s16/cmp_work/repair_r2_tests_01.json
{
  "controls": {
    "mutation_control": "completed",
    "own_labels": "completed",
    "equivalents": "completed",
    "derivative_points": "completed",
    "integral_and_domains": "completed",
    "volume_partial_time": "completed",
    "count_zero_fields": "completed",
    "wrapper_resolution": "completed",
    "repair_ablations": "completed",
    "binding_repoints": "completed",
    "name_repoints": "completed",
    "piece_and_packaging": "completed",
    "domain_and_balance_factor": "completed",
    "coverage_controls": "completed",
    "alias_mutation": "completed",
    "duplicates": "completed",
    "quadrature_path": "completed"
  },
  "synthetic_only": true,
  "runtime_seconds": 19.873062699000002,
  "peak_rss_kib": 79972
}
```

## The code: reviewed bytes

```
$ grep -F -e 'research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py' -e 'research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py' _scratch/s16/cmp_review/s16_cmp_build_baseline_r2.sha256
f9c80a188588e56e155436345d8429d17299a1d0caac725be40cac525af1cd66  research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py
88df585fc55141a91e28016ba279e6835c587ee3176220cdf6005eb6526a6aff  research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py
```

```
$ sha256sum research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py
f9c80a188588e56e155436345d8429d17299a1d0caac725be40cac525af1cd66  research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py
88df585fc55141a91e28016ba279e6835c587ee3176220cdf6005eb6526a6aff  research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py
```

## The inputs: the accepted S16 transcripts

```
$ head -1 _scratch/s16/cmp_work/repair_r2_full_01.jsonl | jq -c '.input_sha256'
{"py":"da72cd98811b03a85d766743e0abad0df80b979f7b2d26cc0dfe230f8f5d1435","wl":"72d352df5f631ae5d71772ec992306fd5e1eb0c83843a69f0d65fae526699567"}
```

```
$ sha256sum research/pde_ledger_v3/scripts/out/S16_moving_region.out research/pde_ledger_v3/mathematica/out/S16_moving_region.out
da72cd98811b03a85d766743e0abad0df80b979f7b2d26cc0dfe230f8f5d1435  research/pde_ledger_v3/scripts/out/S16_moving_region.out
72d352df5f631ae5d71772ec992306fd5e1eb0c83843a69f0d65fae526699567  research/pde_ledger_v3/mathematica/out/S16_moving_region.out
```

```
$ git log -1 --format='%h %s' 1d990ee6
1d990ee6 Accept the S16 build (cleared in review round 1)
```

## The output: identical to the reviewed run

```
$ grep -F '_scratch/s16/cmp_work/repair_r2_full_01.jsonl' _scratch/s16/cmp_review/s16_cmp_build_baseline_r2.sha256
64be2c18529d781f8c6648dc5e0e19975b23d9167beb2fa6a6fc18a5750a7c1f  _scratch/s16/cmp_work/repair_r2_full_01.jsonl
```

```
$ sha256sum _scratch/s16/cmp_work/repair_r2_full_01.jsonl research/pde_ledger_v3/scripts/out/S16_cross_engine_comparator.out
64be2c18529d781f8c6648dc5e0e19975b23d9167beb2fa6a6fc18a5750a7c1f  _scratch/s16/cmp_work/repair_r2_full_01.jsonl
64be2c18529d781f8c6648dc5e0e19975b23d9167beb2fa6a6fc18a5750a7c1f  research/pde_ledger_v3/scripts/out/S16_cross_engine_comparator.out
```

```
$ cmp _scratch/s16/cmp_work/repair_r2_full_01.jsonl research/pde_ledger_v3/scripts/out/S16_cross_engine_comparator.out && echo COMPARISON_IDENTICAL_TO_REVIEWED_RUN
COMPARISON_IDENTICAL_TO_REVIEWED_RUN
```
