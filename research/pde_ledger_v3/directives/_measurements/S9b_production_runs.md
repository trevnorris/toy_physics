# Measurements — S9b production runs (generated 2026-10-12 08:00)

Generator: `_scratch/s9b_build/gen/s9b_production_record.sh` (sha256/cat/cmp/wc/git only; this file is written by it).
The engines are the reviewed round-5 engines (`bfa034c4`). Both build-review r5 legs found them sound; the open
finding of that review is in the two ablation harnesses, which repair round 6 changes without changing the engines.
By the user's decision (2026-10-12) the production runs and the comparator proceed on these engines meanwhile.
Both runs used the repository guard in pooled mode (pool `s9b-production`, 16 GiB; Wolfram also `--tasks-max 64`),
from the main checkout, run serially by `_scratch/s9b_build/gen/s9b_production_run.sh`.

## Commands
```
$ cat _scratch/s9b_build/gen/s9b_production_run.sh
#!/bin/bash
# Production runs of the S9b engines (the reviewed r5 engines, unchanged by round 6), serially through the guard.
cd /var/projects/toy_physics || exit 1
python3 scripts/s11c_guarded_run.py --pool s9b-production --memory-gib 16 --log-directory _scratch/s11c/s9b-production/py-audit -- /usr/bin/python -B research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py
echo "py exit $?"
python3 scripts/s11c_guarded_run.py --pool s9b-production --memory-gib 16 --tasks-max 64 --log-directory _scratch/s11c/s9b-production/wl-audit -- math -script research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
echo "wl exit $?"
```

```
cp _scratch/s11c/s9b-production/py-audit/stdout research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out
cp _scratch/s11c/s9b-production/wl-audit/stdout research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
```

## Guard outcomes
```
$ cat _scratch/s11c/s9b-production/py-audit/child-outcome.json _scratch/s11c/s9b-production/wl-audit/child-outcome.json
{
  "exitCode": 0,
  "guardReason": null,
  "wallSeconds": 423.009401806863,
  "childPid": 317726
}
{
  "exitCode": 0,
  "guardReason": null,
  "wallSeconds": 23.32588361017406,
  "childPid": 321455
}
```

```
$ wc -c < _scratch/s11c/s9b-production/py-audit/stderr; wc -c < _scratch/s11c/s9b-production/wl-audit/stderr
18452
0
```

```
$ head -3 _scratch/s11c/s9b-production/py-audit/stderr
S9B_WORK BRANCH_EXISTENCE/PATH_TRAVERSAL/BRANCH_TYPE: dispersion and signed speed
S9B_WORK RAY_DOMAIN FLYBY: turning and angular closure
S9B_WORK RAY_DOMAIN ER: turning and angular closure
```

## Engine sources run (the round-5 baseline; only the directive differs, under amendment 9 review)
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_review_baseline_r5.sha256 2>&1 | grep -E 'audit|exports|FAILED|WARNING'
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/scripts/S9b_exports.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: FAILED
sha256sum: WARNING: 1 computed checksum did NOT match
```

```
$ git diff --quiet bfa034c4 -- research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl research/pde_ledger_v3/scripts/S9b_exports.py && echo ENGINES_AND_EXPORT_UNCHANGED_SINCE_bfa034c4
ENGINES_AND_EXPORT_UNCHANGED_SINCE_bfa034c4
```

## Transcripts filed, and the same bytes as the builders' final baseline runs
```
$ sha256sum _scratch/s11c/s9b-production/py-audit/stdout research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out _scratch/s9b_build/r5_harness_output/py/baseline/stdout
b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142  _scratch/s11c/s9b-production/py-audit/stdout
b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142  research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out
b0b3e1aa2eb0928be33ebec553efc570e6a883394a7b484548e52cb4685cf142  _scratch/s9b_build/r5_harness_output/py/baseline/stdout
```

```
$ sha256sum _scratch/s11c/s9b-production/wl-audit/stdout research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out _scratch/s9b_build/r5_harness_output/wl/engine-repository/stdout
6519d10c86c9b16741bd32814a57216116c33717a0aec4150f61903617e9ffdb  _scratch/s11c/s9b-production/wl-audit/stdout
6519d10c86c9b16741bd32814a57216116c33717a0aec4150f61903617e9ffdb  research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out
6519d10c86c9b16741bd32814a57216116c33717a0aec4150f61903617e9ffdb  _scratch/s9b_build/r5_harness_output/wl/engine-repository/stdout
```

```
$ cmp _scratch/s11c/s9b-production/py-audit/stderr _scratch/s9b_build/r5_harness_output/py/baseline/stderr && echo PY_STDERR_IDENTICAL_TO_BUILDER
PY_STDERR_IDENTICAL_TO_BUILDER
```

