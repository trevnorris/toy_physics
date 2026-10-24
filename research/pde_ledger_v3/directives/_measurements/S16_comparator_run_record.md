# S16 comparator run record — K5 citation repair

Reviewed baseline: `de28e059`. This repair addresses only K5 in `S16_comparator_review_r1_disposition.md`, under `S16_comparator_decision_list.md`. Every printed SymPy citation in the binding rules/table, surface table, name table and reason strings was checked against the engine accepted at `1d990ee6`. The comparator changes only citation-string literals; the tests, Wolfram citations and all comparison behavior are unchanged. No commit was made. No residual target was given.

Commands, from `/var/projects/toy_physics`, in execution order:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s16_cmp_guard/repair_r2_tests_01 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_tests_01.json
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s16_cmp_guard/repair_r2_full_01 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S16_moving_region.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S16_moving_region.out --output /var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_full_01.jsonl
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s16_cmp_guard/repair_r2_audit_01 -- python3 -B /var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_audit.py
```

| Run | Guard outcome | Guard wall seconds | Observed cgroup peak bytes | Script runtime seconds | Peak process RSS KiB |
|---|---|---:|---:|---:|---:|
| `repair_r2_tests_01` | Admitted; `exitCode: 0`, `limitsVerified: true`, child `exitCode: 0`, `guardReason: null` | 20.440217527997447 | 78635008 | 19.873062699000002 | 79972 |
| `repair_r2_full_01` | Admitted; `exitCode: 0`, `limitsVerified: true`, child `exitCode: 0`, `guardReason: null` | 7705.871029508009 | 1152303104 | 7699.739048800999 | 1036104 |
| `repair_r2_audit_01` | Admitted; `exitCode: 0`, `limitsVerified: true`, child `exitCode: 0`, `guardReason: null` | 3.6937131169979693 | — | — | — |

Each job ran singly under the 53687091200-byte cap, zero swap, 32-task cap and one logical CPU, with no wall/CPU/inactivity deadline. For the tests and full run, sampled swap has maximum 0 and terminal memory-limit/OOM event counters are 0. Guard receipts and literal stdout/stderr are under `_scratch/s11c/s16_cmp_guard/` in the named directories. No accepted engine or Wolfram kernel was run.

- Tests: [output](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_tests_01.json), SHA256 `5805b68fc52fd434f20c07f789f462579b81aac87f0106f8e7db4d47fe37d44c`. All 17 synthetic control groups completed.
- Full comparison: [output](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_full_01.jsonl), SHA256 `64be2c18529d781f8c6648dc5e0e19975b23d9167beb2fa6a6fc18a5750a7c1f`.
- Citation/source/output audit: [output](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_audit_01.json), SHA256 `83c6b57ae5eb766cff2a67c7c3dfc1bb95608f16b07836637110b4961a0ae146`. It includes the accepted source excerpts for all 49 distinct citation specifications, their uses, the 38 changed string literals, full accounting and run receipts. Literal stdout is `_scratch/s11c/s16_cmp_guard/repair_r2_audit_01/stdout`.
- Final comparator SHA256 `f9c80a188588e56e155436345d8429d17299a1d0caac725be40cac525af1cd66`; unchanged tests SHA256 `88df585fc55141a91e28016ba279e6835c587ee3176220cdf6005eb6526a6aff`. Both were byte-compared after the run with the copies saved before launch in `_scratch/s16/cmp_work/repair_r2_final_code/`.
- Accepted transcript SHA256: SymPy `da72cd98811b03a85d766743e0abad0df80b979f7b2d26cc0dfe230f8f5d1435`; Wolfram `72d352df5f631ae5d71772ec992306fd5e1eb0c83843a69f0d65fae526699567`. Sources and transcript paths remain unchanged against `1d990ee6`.

The mechanical audit compares the comparator AST against `de28e059`, allowing only line-specification strings and SymPy citation spans in reason strings to change. It also checks the tests byte-for-byte against that commit and every citation's bounds against the 439-line accepted engine. Citation semantics were checked by reading the accepted source; the mechanical bounds check alone does not establish semantics.

The audit compares every new output record with `repair_r1_full_01.jsonl`, removing only explicit SymPy citation fields (`py_line`, the SymPy slot of each binding rule), SymPy citation spans in reason strings/summary keys, and `runtime_seconds`/`peak_rss_kib` in the resource record. It reports `non_citation_output_differences: 0`. Only the configuration record, 60779 tag-accounting records, the resource record and the summary differ before that removal. All operands, exact/resolved residuals, samples, labels, comparison coverage and leaf counts are unchanged. The round 1 output retains SHA256 `2cab8b06e674316e83049557191628eaf963eecfe8b9676ccb5558b950f80c78`.

The first full-output record prints all 47 binding rules (`RULES`, comparator line 292), all 2028 generated bindings (`BINDING_TABLE`, line 399), all 91 names (`name_contract`, line 239), and the two-row surface table. All 2028 bindings were compared; none was not formed or not compared. Per-tag accounting covers 26028 Py / 37578 WL tags, with 238488 / 580543 parsed leaves and 161466 / 212319 compared leaves. Bound-end parsed leaves are 161466 / 213171. The 426 bindings below their parsed counts are `coverage:1602` and `comparison_notes:1603` through `comparison_notes:2027`; the message is compared and WL `Missing` constructor/reason metadata remains unpaired.

Unjoined tags total 24322 Py / 36021 WL. Every tag and unpaired component remains printed, with these grouped reasons:

| Engine | Tags | Reason and corrected source citations |
|---|---:|---|
| Py / WL | 24264 / 36014 | Native dimension-sum audits enumerate different expression trees, without ordinal pairing (Py 409–434; WL 238–243). |
| Py / WL | 26 / 6 | Free-radius chart and split centre/radius velocities/speeds differ from surface total velocities/speeds (Py 51–54, 69–73, 113, 118–121; WL 56–58). Surface charts and rate pieces are bound; totals cannot bind individual parts. |
| Py | 6 | Dimensions of extra bulk-reference entries in `momentumSurfaceMotion` have no unused WL counterpart (Py 215–216, 359–407, 425–434; WL 130, 230–243). |
| Py | 21 | Count/dipole/momentum `weightedLocalRate`, `afterParts.bulk` and `afterParts.current` are unwindowed; WL's intermediates are window-weighted (Py 171–173; WL 87–95, 101). |
| Py | 5 | Engine-only Leibniz and generic pullback diagnostics; WL has different window intermediates (Py 60–63, 122–128; WL 74–95). |
| WL | 1 | Explicit moment-index labels occur only in WL; Py encodes indices in tags (Py 198–200; WL 119). |

The decision list's named geometry left unjoined is Py `geometry.chart`, `geometry.centerVelocity`, `geometry.radiusVelocity`, `geometry.centerNormalSpeed`, `geometry.radiusNormalSpeed` and WL `surfaces.velocity`, `surfaces.normalSpeed`, for the differing-object reason above. No shared requested balance, moment, rate or printed derivative piece is left unjoined.

Evaluation is unchanged: mixed quadratic/cubic/quartic live profiles, printed in 45 profile records, at three rational `(x1,x2,x3,t)` points `(1/5,2/7,-1/4,1/6)`, `(-2/9,3/8,4/11,-1/7)`, `(3/10,-1/3,2/5,2/9)`. Strata are generic, centre-velocity-zero, radius-velocity-zero, all-velocities-zero and count-zero; higher geometry derivatives stay live and `0 < r_in < R`. The count-zero density is its generic profile minus its exact spatial shell mean at the probe time; both the printed count and `N` vanish there. There are 39663 evaluated samples, 162 normalized-domain exclusions, no unresolved coverage and no numerical integral fallback. Synthetic controls exercise the 40/60-digit quadrature/error-estimate path.

Machinery remains the S14a_1 (`de81b9fc`) AST/InputForm readers and derivative-slot translation, with the existing S16 bindings, aliases, profiles, integral evaluation and accounting. The guarded reruns re-exercised all existing controls and the accepted streams; the additional mechanical audit verified unchanged computation/output apart from citations and resource measurements. Earlier outputs remain intact; previous code, tests and run record were copied to `_scratch/s16/cmp_work/repair_r2_baseline/`. Open questions: the non-author scoped citation/change check and interpretation of residuals belong to the orchestrator/S16 record; neither was performed here.
