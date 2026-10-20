# S12 comparator run record — 2026-10-20

Contract: `S12_comparator_decision_list.md` at `3822b64e`; accepted engine sources/transcripts at `a3da43d9`; no commits made.
Final comparator SHA256: `9e1744595229168c2870224398c5d5bd504654d629e548e1b652882fa2566ce5`.
Final tests SHA256: `4c290d4db227406d5f3d3a9a5c657827913d2a3fd57ddd73d8007f771166f978`.
Inputs SHA256: {"py": "988149e6eb2162c78329eb6cbfa90762bf505fea448c716f7de2f33686036a97", "wl": "ec0d1acdd0666456cb50ce785b0c5188208b55e88f19b339535a97bcaf7c9055"}.

tests:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/tests_final_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s12/cmp_work/tests_final_02.out
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/tests_final_02/outcome.json`):
```json
{"exitCode":0,"wallSeconds":7.586903855000855,"unit":"s11c-guard-ce035d9d0887","stderrBytes":0,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":7.520178756996756,"childPid":2035278}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/tests_final_02.out`; SHA256 `96c7e323dc0c052a73852fa79698ce6963aeda544651ced5d1d7255020c566dc`.
Cgroup peak memory: `71036928` bytes.

full comparison:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/comparison_final_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S12_order_source.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S12_order_source.out --output /var/projects/toy_physics/_scratch/s12/cmp_work/comparison_final_02.jsonl
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/comparison_final_02/outcome.json`):
```json
{"exitCode":0,"wallSeconds":226.24581323101302,"unit":"s11c-guard-3a9d18babd7a","stderrBytes":26298,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":226.18029009798192,"childPid":2035341}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/comparison_final_02.jsonl`; SHA256 `8e5d01ce6736c9106d4e542edfd108cd0abf8842212adf6a66847e687104f3a9`.
Cgroup peak memory: `195858432` bytes.

Test result (literal last lines):
```text
----------------------------------------------------------------------
Ran 30 tests in 7.076s

OK
{"tests_run": 30, "errors": 0, "failures": 0, "runtime_seconds": 7.076304531015921, "peak_rss_kib": 72688}
```
Full-run resources and accounting (literal output records):
```json
{"kind":"resources","runtime_seconds":225.47244495898485,"peak_rss_kib":145548}
{"kind":"summary","counts":{"unjoined_py":84,"unjoined_wl":81,"joined":369,"joined_with_not_compared":46,"joined_with_not_formed":2,"joined_below_parsed_count":144},"exit_semantics":"zero for any findings; nonzero only for operational failure"}
```
No admission refusal or guard kill occurred. No time limit was supplied. No source engine was run.
Mechanical appendix extraction: `python3 /var/projects/toy_physics/_scratch/s12/cmp_work/assemble_report.py`; literal stdout is `_scratch/s12/cmp_work/assembly.out`.
Earlier scratch revisions remain in the same work/guard directories; their outputs are not the final run evidence.
Builder report: `S12_comparator_builder_report.md`; complete extracted tables and hashes: `_scratch/s12/cmp_work/run_manifest.jsonl`.
No residual target was given.
