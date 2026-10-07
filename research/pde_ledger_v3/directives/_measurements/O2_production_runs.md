# Measurements — O2 production runs (generated 2026-10-07 13:20)

Generator: `_scratch/s9b_build/gen/o2_production_record.sh` (sha256/cat/sed/cmp only; this file is written by it). The engines are the accepted r3 (`0f2e0af8`). Both runs used the repository guard in pooled mode (pool `o2`, 6 GiB; Wolfram also `--tasks-max 64`, user-authorized for the O2 Wolfram job), from the main checkout, with no other guarded job live.

## Commands
```
python3 scripts/s11c_guarded_run.py --pool o2 --memory-gib 6 --log-directory _scratch/s11c/o2-production/py-audit -- python3 research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py
python3 scripts/s11c_guarded_run.py --pool o2 --memory-gib 6 --tasks-max 64 --log-directory _scratch/s11c/o2-production/wl-audit -- math -script research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl
cp _scratch/s11c/o2-production/py-audit/stdout research/pde_ledger_v3/scripts/out/O2_live_balance_sympy_audit.out
cp _scratch/s11c/o2-production/wl-audit/stdout research/pde_ledger_v3/mathematica/out/O2_live_balance_mathematica_audit.out
```

## Guard outcomes
```
$ cat _scratch/s11c/o2-production/py-audit/child-outcome.json _scratch/s11c/o2-production/wl-audit/child-outcome.json
{
  "exitCode": 0,
  "guardReason": null,
  "wallSeconds": 94.00317881791852,
  "childPid": 3039058
}
{
  "exitCode": 0,
  "guardReason": null,
  "wallSeconds": 6.2172754080966115,
  "childPid": 3039255
}
```

```
$ wc -c < _scratch/s11c/o2-production/py-audit/stderr; wc -c < _scratch/s11c/o2-production/wl-audit/stderr
0
0
```

## Engine sources run (accepted r3 baseline)
```
$ sha256sum -c _scratch/s9b_build/o2_build_review_baseline_r3.sha256 2>&1 | grep -v O2_exports
research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py: OK
research/pde_ledger_v3/scripts/O2_live_balance_sympy_ablation.py: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl: OK
research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_ablation.wl: OK
sha256sum: WARNING: 1 computed checksum did NOT match
```

## Transcripts filed
```
$ sha256sum _scratch/s11c/o2-production/py-audit/stdout research/pde_ledger_v3/scripts/out/O2_live_balance_sympy_audit.out _scratch/s11c/o2-production/wl-audit/stdout research/pde_ledger_v3/mathematica/out/O2_live_balance_mathematica_audit.out
b5c86fa328f6b60b39502d3e085af39e8e5f819ea2716f34057dcdb409e2de0d  _scratch/s11c/o2-production/py-audit/stdout
b5c86fa328f6b60b39502d3e085af39e8e5f819ea2716f34057dcdb409e2de0d  research/pde_ledger_v3/scripts/out/O2_live_balance_sympy_audit.out
a6a4fb6759b5ef100944c07f35e2ebb94e5f8084d883d67cf44b8aad10f51d36  _scratch/s11c/o2-production/wl-audit/stdout
a6a4fb6759b5ef100944c07f35e2ebb94e5f8084d883d67cf44b8aad10f51d36  research/pde_ledger_v3/mathematica/out/O2_live_balance_mathematica_audit.out
```

## Regenerated export: Dummy-index-only drift from the accepted r3 export
```
$ sha256sum _scratch/s9b_build/O2_exports_r3_accepted.py research/pde_ledger_v3/scripts/O2_exports.py
e1ec7512bdaa33570e374908be0fa7ff8013ddec697d4e4ccf016c07d37f8180  _scratch/s9b_build/O2_exports_r3_accepted.py
3c3b46a37d5918ba1edd51c8f8acffee97224bc2845ab1779c150d0093d3778b  research/pde_ledger_v3/scripts/O2_exports.py
```

```
$ grep -o 'dummy_index=-\?[0-9]*' _scratch/s9b_build/O2_exports_r3_accepted.py | sort | uniq -c; grep -o 'dummy_index=-\?[0-9]*' research/pde_ledger_v3/scripts/O2_exports.py | sort | uniq -c
  51561 dummy_index=8369554716719071657
  51561 dummy_index=-299293262346107954
```

```
$ cmp <(sed -E 's/dummy_index=-?[0-9]+/dummy_index=N/g' _scratch/s9b_build/O2_exports_r3_accepted.py) <(sed -E 's/dummy_index=-?[0-9]+/dummy_index=N/g' research/pde_ledger_v3/scripts/O2_exports.py) && echo IDENTICAL_AFTER_MASK
IDENTICAL_AFTER_MASK
```

