# S14a part 1 comparator run record

H1–H3 repair of reviewed baseline `88bf7bf9`, under contract `43174e8c`. Accepted read-only engine inputs: `94fddf2f`. Commands ran from `/var/projects/toy_physics`, one job at a time. No commits were made.

Final source sha256s:

- `research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py`: `d628d27a11eba4b67ac9ddcce7d73b6bd65f8325c30e7c8277f501b49da4e7c9`.
- `research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py`: `65e96aeeab41aa362d2148e745e7cf7f2a94d5cf37846157597f7229b0217812`.

H1: fresh integral binders followed by one capture-free lexical alpha renaming. H2: quadratic probes couple every pair of arguments; controls cover 1190 first-partial/evaluation-argument pairs and a shifted face derivative inside a shell integral. H3: the spaced comparison note uses the label path; active candidate accounting is computed and semantic rationale is explicitly curated.

## Synthetic tests

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a1-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/r1_tests_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_tests_01.out
```

Guard admitted the job and verified its limits: guard exit `0`, child exit `0`, `guardReason: null`. No deadline; zero swap. [Guard outcome](/var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/r1_tests_01/outcome.json).
Output: [r1_tests_01.out](/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_tests_01.out); sha256 `5d0738cc11a5126cb49380f249ebc2edc3ab11ee3fb38cfe8c5b1083d96520a5`.
Guard wall time: `76.3144247730088` seconds; cgroup peak `memory.peak`: `105504768` bytes.

## Full accepted-transcript comparison

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a1-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/r1_full_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out --output /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl
```

Guard admitted the job and verified its limits: guard exit `0`, child exit `0`, `guardReason: null`. No deadline; zero swap. [Guard outcome](/var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/r1_full_01/outcome.json).
Output: [r1_full_01.jsonl](/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl); sha256 `1db2009d9fa4b77e5d91fd138c43c6d6342585bb75ac56cd2d2f7c0d268c6ce2`.
Guard wall time: `106.72195675599505` seconds; cgroup peak `memory.peak`: `216891392` bytes.

Tests: `32` run, `0` errors, `0` failures; script runtime `75.62683781099622` seconds, peak RSS `106296` KiB. Full-comparison script runtime `105.5439142069954` seconds, peak RSS `200432` KiB.

[Accounting/coverage inventory](/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01_inventory.json); [partial-tag details](/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01_unpaired_components.json); [commands and resource provenance](/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_run_inventory.json).

Accepted transcript sha256s:

- `scripts/out/S14a_1_projection.out`: `3bfb17aa5205e7f0ddec1eea44c920d4f0bb4e3dffe4ccf6fcc74abe88fbe0db`.
- `mathematica/out/S14a_1_projection.out`: `7911ff19f667a388a2a54a3eb0d076ed06143d791ddb2b5bccfe397588652d5b`.

All `22` files previously under `cmp_work/` are unchanged, verified against `r1_existing_scratch_manifest.json`. Every repair output has a new `r1_` name. No admission refusal or guard kill occurred.

No residual target was given. No measured residual is interpreted. Completeness of comparison of named objects rests on binding-table curation; it is not computed by accounting.
