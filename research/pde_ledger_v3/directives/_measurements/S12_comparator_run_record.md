# S12 comparator repair run record — 2026-10-20

Contract: `S12_comparator_decision_list.md` at `3822b64e`; reviewed baseline `f1d075b8`; scope F1–F4 and N1–N2. Accepted engine sources/transcripts `a3da43d9`. No commit.
Final comparator SHA256: `70c8907b77200a27beb21f84eef838eba8566661b1b02faeeb125bbb46daf40b`.
Final tests SHA256: `99becfd23ae99e24955a8ed3f53e24784edc753daaf4e121c6f70150a4cb18e9`.
Inputs SHA256: {"py": "988149e6eb2162c78329eb6cbfa90762bf505fea448c716f7de2f33686036a97", "wl": "ec0d1acdd0666456cb50ce785b0c5188208b55e88f19b339535a97bcaf7c9055"}.
All runs were serial, pooled, with the requested 50 GiB cap and no time limit. Only the last two entries are evidence from the final code.

superseded test harness:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_tests_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s12/cmp_work/r1_tests_01.out
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_tests_01/outcome.json`):
```json
{"exitCode":1,"wallSeconds":23.546388293005293,"unit":"s11c-guard-3d29eef66a13","stderrBytes":0,"limitsVerified":true,"childOutcome":{"exitCode":1,"guardReason":null,"wallSeconds":23.47866639701533,"childPid":2054442}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/r1_tests_01.out`; SHA256 `2ab6f2597613d30dd44486132d24882f28f578b098e36d9788fca6dd15bb8d11`.
Cgroup peak: `82866176` bytes; literal process resources: `{"tests_run":36,"errors":0,"failures":3,"runtime_seconds":23.075704388000304,"peak_rss_kib":84168}`.

superseded tests before composite-argument edge repair:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_tests_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s12/cmp_work/r1_tests_02.out
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_tests_02/outcome.json`):
```json
{"exitCode":0,"wallSeconds":23.895485198008828,"unit":"s11c-guard-668eb5822ace","stderrBytes":0,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":23.824368446017615,"childPid":2054613}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/r1_tests_02.out`; SHA256 `d69d261fb03af96129a24ce942161bf42000e5787a9494cb6c867e58c2542da4`.
Cgroup peak: `82292736` bytes; literal process resources: `{"tests_run":36,"errors":0,"failures":0,"runtime_seconds":23.40565124701243,"peak_rss_kib":84288}`.

superseded comparison before composite-argument edge repair:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_comparison_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S12_order_source.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S12_order_source.out --output /var/projects/toy_physics/_scratch/s12/cmp_work/r1_comparison_01.jsonl
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_comparison_01/outcome.json`):
```json
{"exitCode":0,"wallSeconds":203.9009824960085,"unit":"s11c-guard-416376c01da1","stderrBytes":26298,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":203.8267485089891,"childPid":2054709}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/r1_comparison_01.jsonl`; SHA256 `cd78e066b666dc503a29cea5b1f00284e6b294a15e6e178cccbaa3a534736727`.
Cgroup peak: `178208768` bytes; literal process resources: `{"kind":"resources","runtime_seconds":203.1521755659778,"peak_rss_kib":139164}`.

final tests:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_tests_03 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S12_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s12/cmp_work/r1_tests_03.out
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_tests_03/outcome.json`):
```json
{"exitCode":0,"wallSeconds":24.494830347015522,"unit":"s11c-guard-3c0b1aaf978e","stderrBytes":0,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":24.426640743011376,"childPid":2055487}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/r1_tests_03.out`; SHA256 `f136e6e85f200a87fa53d6c9a297f0e94f9f2b867fc8bb4d731471a53ad239bb`.
Cgroup peak: `83128320` bytes; literal process resources: `{"tests_run":36,"errors":0,"failures":0,"runtime_seconds":23.941212727018865,"peak_rss_kib":84332}`.

final full comparison:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_comparison_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S12_order_source.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S12_order_source.out --output /var/projects/toy_physics/_scratch/s12/cmp_work/r1_comparison_02.jsonl
```
Guard outcome (literal JSON from `/var/projects/toy_physics/_scratch/s11c/s12_cmp_guard/r1_comparison_02/outcome.json`):
```json
{"exitCode":0,"wallSeconds":204.3378958920075,"unit":"s11c-guard-bae75c79f6a7","stderrBytes":26298,"limitsVerified":true,"childOutcome":{"exitCode":0,"guardReason":null,"wallSeconds":204.27504475598107,"childPid":2055567}}
```
Output: `/var/projects/toy_physics/_scratch/s12/cmp_work/r1_comparison_02.jsonl`; SHA256 `f04c95195312db4d2ed706596db87b740423366438f645bb8fc8183cde04778e`.
Cgroup peak: `178270208` bytes; literal process resources: `{"kind":"resources","runtime_seconds":203.59856851000222,"peak_rss_kib":139108}`.

The first test run exposed three assertions in two harness checks: aggregated sample sets hid the changed object, and Xor already handled the bare-inequality ablation. The final checks isolate the mutated object and exercise nested relation notation.
The next test/comparison pair precedes the final F1 edge repair: a variable entering one argument directly and another through a composite must retain both chain-rule contributions. Those outputs remain superseded.
Final test result (literal last lines):
```text
----------------------------------------------------------------------
Ran 36 tests in 23.941s

OK
{"tests_run": 36, "errors": 0, "failures": 0, "runtime_seconds": 23.941212727018865, "peak_rss_kib": 84332}
```
Final full-run resources and accounting (literal output records):
```json
{"kind":"resources","runtime_seconds":203.59856851000222,"peak_rss_kib":139108}
{"kind":"summary","counts":{"unjoined_py":84,"unjoined_wl":81,"joined":369,"joined_with_not_compared":46,"joined_with_not_formed":2,"joined_below_parsed_count":144},"exit_semantics":"zero for any findings; nonzero only for operational failure"}
```
No admission refusal or guard kill occurred. No source engine was run.
Mechanical appendix extraction: `python3 /var/projects/toy_physics/_scratch/s12/cmp_work/assemble_r1_report.py`; literal stdout is `_scratch/s12/cmp_work/r1_assembly.out`.
Builder report: `S12_comparator_builder_report.md`; complete extracted tables and hashes: `_scratch/s12/cmp_work/r1_run_manifest.jsonl`.
No residual target was given.
