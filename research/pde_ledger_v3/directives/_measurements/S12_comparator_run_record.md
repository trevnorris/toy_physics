# S12 comparator repair run record — 2026-10-20

Contract: `S12_comparator_decision_list.md` at `3822b64e`; reviewed baseline `e4dae246`; scope G1–G2, preserving F1–F4 and N1–N2. Accepted engine sources/transcripts `a3da43d9`. No commit.
Final comparator SHA256: `4ed863f4b56ff865272c781015f86faf3d347b6d46ed592d661bb86b3512f462`.
Final tests SHA256: `b222dd0469f25da8adccfca581472fe958a3feb17df04e9c38da47c6069f602b`.
Inputs SHA256: {"py": "988149e6eb2162c78329eb6cbfa90762bf505fea448c716f7de2f33686036a97", "wl": "ec0d1acdd0666456cb50ce785b0c5188208b55e88f19b339535a97bcaf7c9055"}.
Both runs were serial, pooled, with the requested 50 GiB cap and no time limit, using the final code. Earlier repair evidence is preserved in the run record at `e4dae246`.

final tests:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r2_tests_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s12/cmp_work/r2_tests_01.out
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r2_tests_01/outcome.json`):
```json
{"exitCode":0,"wallSeconds":24.5910155439924,"unit":"s11c-guard-ae42e8fba5d4","stderrBytes":0,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":24.525772786990274,"childPid":2082628}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/r2_tests_01.out`; SHA256 `adc4b8c1c71d0b22405dd6a44ec5eb5cdd933fbabaaf929123f790810249f1b6`.
Cgroup peak: `83120128` bytes; literal process resources: `{"tests_run":38,"errors":0,"failures":0,"runtime_seconds":24.10258426499786,"peak_rss_kib":84872}`.

final full comparison:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r2_comparison_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S12_order_source.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S12_order_source.out --output /var/projects/toy_physics/_scratch/s12/cmp_work/r2_comparison_01.jsonl
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r2_comparison_01/outcome.json`):
```json
{"exitCode":0,"wallSeconds":204.1622238169948,"unit":"s11c-guard-ea21413e930f","stderrBytes":26298,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":204.09423501399579,"childPid":2082699}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/r2_comparison_01.jsonl`; SHA256 `91d4540952112a227a2b46cc8934c068f3e97941ae2e2c479c7f873b54326cf6`.
Cgroup peak: `178499584` bytes; literal process resources: `{"kind":"resources","runtime_seconds":203.43066374000045,"peak_rss_kib":139264}`.

Final test result (literal last lines):
```text
----------------------------------------------------------------------
Ran 38 tests in 24.102s

OK
{"tests_run": 38, "errors": 0, "failures": 0, "runtime_seconds": 24.10258426499786, "peak_rss_kib": 84872}
```
Final full-run resources and accounting (literal output records):
```json
{"kind":"resources","runtime_seconds":203.43066374000045,"peak_rss_kib":139264}
{"kind":"summary","counts":{"unjoined_py":84,"unjoined_wl":81,"joined":369,"joined_with_not_compared":46,"joined_with_not_formed":2,"joined_below_parsed_count":144},"exit_semantics":"zero for any findings; nonzero only for operational failure"}
```
No admission refusal or guard kill occurred. No source engine was run.
Mechanical stored-output check: `python3 /var/projects/toy_physics/_scratch/s12/cmp_work/r2_validate_stored_output.py`; literal stdout is `_scratch/s12/cmp_work/r2_validation.out`.
Output: `_scratch/s12/cmp_work/r2_stored_output_validation.json`; SHA256 `1cfcdcb70792ddde4f074a5b25ca085f4b4e2b77e563552c892d75d4ea8771a5`.
Literal stdout: `{"old_nonlabel_record_count": 3536, "new_nonlabel_record_count": 3536, "nonlabel_records_identical": true, "configuration_except_label_layouts_identical": true, "changed_label_tags": ["calculus_scope", "completion", "dimension_basis", "engine"]}`.
Mechanical appendix extraction: `python3 /var/projects/toy_physics/_scratch/s12/cmp_work/assemble_r2_report.py`; literal stdout is `_scratch/s12/cmp_work/r2_assembly.out`.
Builder report: `S12_comparator_builder_report.md`; complete extracted tables and hashes: `_scratch/s12/cmp_work/r2_run_manifest.jsonl`.
No residual target was given.
