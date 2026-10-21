# S14a part 1 comparator run record

Contract: `43174e8c`; read-only engine inputs: `94fddf2f`. No commits were made. Both final jobs used the same comparator source; no engine was run.

Final source sha256s:
- `research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py`: `1e576f30d433dc4f36e67d1f916205c94cb97ee3427b65d4eed419501b119177`.
- `research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py`: `9dda4f5b082e15b4f67d79a4d72dd21af1343b507fec2a1bbe4709f2adc64be6`.

## Synthetic tests

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a1-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/tests_04 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s14a/cmp_work/tests_04.out
```

Guard admitted the job and verified its limits; guard exit `0`, child exit `0`, `guardReason: null`. No deadline or swap allowance. [Guard outcome](/var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/tests_04/outcome.json).
Output: [tests_04.out](/var/projects/toy_physics/_scratch/s14a/cmp_work/tests_04.out); sha256 `f3190a154f4e7a9b72bbf4a3059924ea50b2869c9e6e2920457275d925ee012b`.
Guard wall time: `53.10758459099452` seconds; cgroup peak `memory.peak`: `103006208` bytes.

## Full accepted-transcript comparison

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a1-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/full_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out --output /var/projects/toy_physics/_scratch/s14a/cmp_work/full_02.jsonl
```

Guard admitted the job and verified its limits; guard exit `0`, child exit `0`, `guardReason: null`. No deadline or swap allowance. [Guard outcome](/var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/full_02/outcome.json).
Output: [full_02.jsonl](/var/projects/toy_physics/_scratch/s14a/cmp_work/full_02.jsonl); sha256 `ff5e03e2ad708db7eb5e2c0b09c7aa445279374e1095f8df8cd713cb69580445`.
Guard wall time: `95.52405627901317` seconds; cgroup peak `memory.peak`: `196194304` bytes.

Tests: `27` run, `0` errors, `0` failures. Full-run accounting and coverage are inventoried in [full_02_inventory.json](/var/projects/toy_physics/_scratch/s14a/cmp_work/full_02_inventory.json).

Accepted transcript sha256s:
- `scripts/out/S14a_1_projection.out`: `3bfb17aa5205e7f0ddec1eea44c920d4f0bb4e3dffe4ccf6fcc74abe88fbe0db`.
- `mathematica/out/S14a_1_projection.out`: `7911ff19f667a388a2a54a3eb0d076ed06143d791ddb2b5bccfe397588652d5b`.

Development evidence is retained under `cmp_work/` and the guard directories: `tests_01` was manually interrupted after recorded test failures (guard exit `130`, child `-2`, `guardReason: null`; source snapshot `tests_01_sources/`); `tests_02` completed; `tests_03` stopped at a missing synthetic-fixture selector; `full_01` is superseded by the explicit dimension-domain handling. No guard refused admission or killed a job. These development runs are not the final evidence above.

No residual target was given. Residuals are not interpreted here.
