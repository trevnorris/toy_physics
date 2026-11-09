# S15 separation-departure final run record

Current repair baseline: user-supplied reviewed commit `0d9675eb`; the preceding repair used `90441b89`. No git command or new commit was made. This repair follows the user decisions of 2026-11-06: Wolfram retains the six deferred mixed-drain region objects with decision-set not-obtained labels; SymPy attempts no item-4 object, including those independent of separation. Its earlier item-4 output is not a result of this build. Neither engine output was compared with the other.

## SymPy

Final source SHA256: `8b4ecfe996341f3fdeeb55725e219fa874eaaf5bad8ec02e15c4da100c7725c8`. Final run `repair4_decision_final_sympy_03`: guard and child exit 0, limits verified; runtime 1820.2957905728836 seconds; peak 388288512 bytes.
Transcript: `research/pde_ledger_v3/scripts/out/S15_separation_departure.out`; SHA256 `10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6`.

Exact argument vector:
```json
["python3", "/var/projects/toy_physics/scripts/s11c_guarded_run.py", "--pool", "s15", "--memory-gib", "2", "--log-directory", "/var/projects/toy_physics/_scratch/s15/guard/repair4_decision_final_sympy_03", "--", "env", "S15_PY_TRANSCRIPT=/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_separation_departure.out", "S15_PY_WORK_DIRECTORY=/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_sympy_03_work", "S15_PY_RUN_NAME=repair4_decision_final_sympy_03", "python3", "-u", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_sympy_03.py"]
```

All 336 named item-4 entries carry `user_decision_sympy_no_item4_2026_11_06`. No conditions or Q under conditions were attempted or printed in this final run. The complete named list, command, limits, source hash and transcript hash are in [the SymPy run manifest](S15_separation_departure_sympy_run_manifest.json). Items 1, 2, 3, 5 and 6 were computed. The same-engine baseline audit retained 18,587 non-item-4 records: 18 textual differences in pressure objects were checked by exact bound-variable canonicalization and held-integrand expansion, and all 18 residuals were zero in the guarded audit `repair4_decision_sympy_baseline_audit_01`. Its operands, residuals and outcome are recorded in the manifest; there was no cross-engine comparison.

History only: `repair4_decision_final_sympy_01` exited 1 before construction because the guard did not preserve a transcript-path environment setting. `repair4_decision_final_sympy_02` was stopped by the orchestrator at 2026-11-06T17:42:28Z with its launcher alive; guard exit 1, 7309.415392439114 seconds. Their recorded outcomes, peaks and source hashes are in the SymPy manifest. Neither run sets a label or supplies an item-4 result.

## Wolfram

Mathematical source SHA256: `366298d82e62d1f36a4ec715feeaa2a6a37903b69a9853d83a4fde50089f19d2`; current runner SHA256 `a5b6c68cc1c489912dfe37331eb275b64c4dd9cfcae6d990cd03e3f1018a7ba6`. Original segment-specific runner hashes are recorded in the manifest and the provenance table below. The original final mathematical source stayed unchanged. After a shared-cap memory failure, execution resumed only unfinished objects, then retried the three affected drain-1 region-rate objects alone at the largest admitted cap. Those continuations did not replay completed objects. The present label correction recomputed only the six sufficient-condition passes described below.
Transcript: `research/pde_ledger_v3/mathematica/out/S15_separation_departure.out`; SHA256 `7d8044652316b879e7d0c61495bcf6fd241ed243009b251687cf5353b2b710fb`.

| Guarded segment | Guard exit | Runtime seconds | Peak bytes |
|---|---:|---:|---:|
| repair4_decision_final_wolfram_01 | 1 | 10590.804926368175 | 8383946752 |
| repair4_decision_final_wolfram_remaining_01 | 0 | 78356.54220686108 | 8589934592 |
| repair4_decision_final_wolfram_alone_01 | 0 | 44205.496096736984 | 26377646080 |

Commands (argument vectors; the initial launch and resume records preserve the exact original paths):

```json
["python3", "/var/projects/toy_physics/scripts/s11c_guarded_run.py", "--pool", "s15", "--memory-gib", "8", "--log-directory", "/var/projects/toy_physics/_scratch/s15/guard/repair4_decision_final_wolfram_01", "--", "python3", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01_runner.py", "--source", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01.wl", "--transcript", "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_departure.out", "--work-directory", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01_work"]
```

```json
["python3", "/var/projects/toy_physics/scripts/s11c_guarded_run.py", "--pool", "s15", "--memory-gib", "8", "--log-directory", "/var/projects/toy_physics/_scratch/s15/guard/repair4_decision_final_wolfram_remaining_01", "--", "env", "S15_WL_RUN_NAME=repair4_decision_final_wolfram_remaining_01", "python3", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_remaining_01_runner.py", "--source", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01.wl", "--transcript", "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_departure.out", "--work-directory", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01_work", "--resume", "--exclude-object", "region/component_1/rate_labels_{drain1}"]
```

```json
["python3", "/var/projects/toy_physics/scripts/s11c_guarded_run.py", "--pool", "s15", "--memory-gib", "56", "--log-directory", "/var/projects/toy_physics/_scratch/s15/guard/repair4_decision_final_wolfram_alone_01", "--", "env", "S15_ALONE_MAXIMUM_CAP=1", "S15_WL_RUN_NAME=repair4_decision_final_wolfram_alone_01", "python3", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_alone_01_runner.py", "--source", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01.wl", "--transcript", "/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_departure.out", "--work-directory", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01_work", "--resume"]
```

The shared-cap failures are history, not out-of-memory labels. All three alone retries returned their conditions before six hours. No final Wolfram object carries a six-hour or memory-exhaustion label. The full per-kernel record and admission evidence are in [the Wolfram run manifest](S15_separation_departure_wolfram_run_manifest.json).

Wolfram labelled objects:

- `region/component_1/momentum_labels_{drain1, drain2}` — `user_decision_2026_11_06`.
- `region/component_1/rate_labels_{drain1, drain2}` — `user_decision_2026_11_06`.
- `region/component_2/momentum_labels_{drain1, drain2}` — `user_decision_2026_11_06`.
- `region/component_2/rate_labels_{drain1, drain2}` — `user_decision_2026_11_06`.
- `region/component_3/momentum_labels_{drain1, drain2}` — `user_decision_2026_11_06`.
- `region/component_3/rate_labels_{drain1, drain2}` — `user_decision_2026_11_06`.

For each of the six `user_decision_2026_11_06` objects, the original final run computed its printed sufficient conditions and Q under each. It made no exponent-sensitive attempt and did not obtain the complete conditions. The label-correction run independently recomputed those same sufficient conditions and Q under each, again without an exponent-sensitive attempt or complete conditions. Its six replacement transcript blocks carry the third-kind label set by the user decision. All 42 nondeferred separation-dependent item-4 traversals returned their computed conditions and restricted Q objects, with empty unresolved-branch lists. No complete condition set is claimed for any labelled object. The same-engine preservation audit retained 25,252 records outside the repaired item-4 entries. Its 270 text differences in retained item-4 intermediates were exact renamings of declared local parameter and direction-coefficient slots; no other retained record changed. All 210 prior numbered sufficient conditions and their 210 restricted Q objects were preserved. The declaration-derived maps and checks are in the manifest.

## Earlier attempts and larger-memory follow-up

[The attempt history](S15_separation_departure_attempt_history.md) and [its detailed inventory](S15_separation_departure_attempt_history.json) list the runs on each of the six mixed-drain region momentum/rate objects, for all three components, in each engine. They include trial designations, recorded results and later corrections (including the SymPy selection correction), start/end times, overlap with other guarded jobs, last steps, peaks and source provenance. These runs are history only and do not assign a six-hour or memory label to the six deferred objects.

Three early overwritten probe entry files have no recoverable historical hash: `repair3_90441_direction_coefficients_probe_01`, `_02`, and `repair3_90441_endpoint_basis_01`. The inventory explicitly records that limitation; hashes of current bytes are not presented as hashes of those runs. Preserved classifier-module hashes are supplied where available.

Larger-memory follow-up queue from the user decision: `region/component_i/momentum_labels_{drain1, drain2}` and `region/component_i/rate_labels_{drain1, drain2}`, for each `i=1,2,3`, in each engine. The current build makes no further exponent-sensitive attempt on those six. Wolfram computes and prints their sufficient conditions and Q under each; these do not constitute complete conditions. SymPy’s additional exclusion of all item 4 is recorded separately above.

The preceding unlimited-series trial report is `_scratch/s15/reading_build_work/repair3/series_unlimited/series_unlimited_report.md`, with inventory and source hashes beside it. Wolfram exhausted memory; the SymPy unlimited trial and size-rate measurement ended on the later user decision. No trial outcome is a final item-4 label.


## Round-5 label correction and final transcript provenance

The mathematical sources are unchanged. Only the Wolfram runner's six decision-label descriptions changed. The preserved r4 copies in both `out/preserved/` directories were read without modification. SymPy's final transcript is reused whole from `repair4_decision_final_sympy_03`; its SHA256 is unchanged, and its script is both source and direct runner (SHA256 `8b4ecfe996341f3fdeeb55725e219fa874eaaf5bad8ec02e15c4da100c7725c8`). No SymPy CAS run was required.

`repair5_wolfram_labels_01` ran the six Wolfram sufficient-condition passes with one fresh kernel per object. Guard and child exit 0, limits verified; runtime 254.77337615098804 seconds; peak 1176899584 bytes; cap 8 GiB. No exponent-sensitive classifier ran. Exact command:

```json
["python3", "/var/projects/toy_physics/scripts/s11c_guarded_run.py", "--pool", "s15", "--memory-gib", "8", "--log-directory", "/var/projects/toy_physics/_scratch/s15/guard/repair5_wolfram_labels_01", "--", "python3", "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/round5_labels/rerun_six.py"]
```

Source SHA256: `366298d82e62d1f36a4ec715feeaa2a6a37903b69a9853d83a4fde50089f19d2`. Partial-pass runner SHA256: `c6b6eba716c89afa5db44ec4eff7b02ebef7ba0f8f6af2c39f5a7d6b8916f504`; imported label-policy runner SHA256: `a5b6c68cc1c489912dfe37331eb275b64c4dd9cfcae6d990cd03e3f1018a7ba6`. The run manifest records each kernel's command, environment, start/end, result and part hash.

The final Wolfram transcript is assembled from preserved r4 bytes `[0,333550018)`, then the six newly generated blocks `part_1.out` through `part_6.out`, then preserved r4 bytes `[398077694,EOF)`. The three original guarded segments supplied the unreplaced ranges, by appending base and item-4 kernel output in the manifest's recorded execution order, including partial failed passes and successful retries. Original source/runner provenance:

| Run | Mathematical source SHA256 | Runner SHA256 |
|---|---|---|
| repair4_decision_final_wolfram_01 | `366298d82e62d1f36a4ec715feeaa2a6a37903b69a9853d83a4fde50089f19d2` | `ca79558aedf0fcfbb3c3c89e79fa48e5cb40767433947c121eed4d427bcce6b4` |
| repair4_decision_final_wolfram_remaining_01 | `366298d82e62d1f36a4ec715feeaa2a6a37903b69a9853d83a4fde50089f19d2` | `4aee6de98cafec2f5824a4536aacb9aec9efd7c0c32c771d2bd51c254cd6f61a` |
| repair4_decision_final_wolfram_alone_01 | `366298d82e62d1f36a4ec715feeaa2a6a37903b69a9853d83a4fde50089f19d2` | `4aee6de98cafec2f5824a4536aacb9aec9efd7c0c32c771d2bd51c254cd6f61a` |

Join and byte-audit command: `python3 /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/round5_labels/assemble.py`; script SHA256 `f04978d6fcd979938f940b1202bfab6a480f490df7c57e059e934cfa4de41378`. Every non-label object in the six regenerated blocks is byte-identical to its preserved counterpart, including 60 sufficient-condition records and 60 Q-under-condition records. The only changed lines are 18 descriptive label records and 300 progress diagnostics whose local Module serial numbers differ; their text is otherwise identical. Both unreplaced ranges were checked byte-for-byte by their hashes, and the joined file was reread and hashed. The preserved Wolfram SHA256 remains `5779d92eb7820b11fa51068526145c3f61e04b69d2f10f86d6b1af923531abce`; the preserved SymPy hash and current SymPy hash both remain the value above. No cross-engine output comparison was made.

The two run manifests contain the source and runner hashes, exact commands, original segment provenance, replacement part hashes, byte offsets, assembly script hash and complete preservation audit. The six current labels all say: “User decision 2026-11-06: this not-obtained label is set by that decision. This run computed the printed sufficient conditions and Q under each. It made no exponent-sensitive attempt and did not obtain the complete conditions. The run record lists the attempt history and the larger-memory follow-up.”

This round changes `mathematica/S15_separation_departure_run.py`, its final transcript, this record, and the two run manifests. Earlier artifact changes below are history of the preceding repair. No new commit was made.

## Changed artifacts in the preceding repair

- `scripts/S15_separation_departure.py`: removed item-4 attempts and their helper construction; prints the fourth decision label for every item-4 entry.
- `mathematica/S15_separation_departure.wl`: prints the six decision labels through the runner, preserves the independently constructed item-4 conditions for other objects, and does not turn returned unresolved branches into unauthorized final labels.
- `mathematica/S15_separation_departure_run.py`: preserves completed work across guarded segments, records shared-cap failures for alone retries, and applies only the authorized final resource-label rules.
- Both `*/out/S15_separation_departure.out` transcripts: regenerated plain files, with source and output hashes recorded above.
- This record, the two engine run manifests, and the Markdown/JSON attempt-history pair: commands, outcomes, per-object labels, preservation evidence and historical provenance.
- Scratch unlimited-trial report, inventory, manifest and `series_unlimited_files.sha256`: completed the preceding trial report from the files left by those runs.
