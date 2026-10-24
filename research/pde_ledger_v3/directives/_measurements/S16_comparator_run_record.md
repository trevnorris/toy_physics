# S16 comparator run record — round 0 repair

The final comparator and synthetic tests repair K1–K4 in `S16_comparator_review_r0_disposition.md`, within `S16_comparator_decision_list.md`. Reviewed baseline: `98936790`; comparator gate: `5994e9ac`; accepted engine sources/transcripts: `1d990ee6`; scientific decision list: `0464b80a`. No commit was made. No residual target was given.

K1 preserves fixed-Cartesian partial time derivatives when evaluating volume integrands at the chart. K2 subtracts the generic density's exact shell mean at the probe time, a constant in space and time, leaving its derivatives live while the printed count and `N` vanish there. K3 prints `resolved_exact_residual` beside the unchanged as-computed field, preserving derivative slots, orders and points. K4 corrects the SymPy introduction citations to 101 (`N`), 102 (`D`) and 103 (`P`). The four added synthetic control groups cover these computational repairs and their individual ablations; no measured transcript is loaded by the tests.

Commands, from `/var/projects/toy_physics`, in execution order, with the final code:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s16_cmp_guard/repair_r1_tests_final_01 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s16/cmp_work/repair_r1_tests_final_01.json
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s16-comparator --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s16_cmp_guard/repair_r1_full_01 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S16_moving_region.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S16_moving_region.out --output /var/projects/toy_physics/_scratch/s16/cmp_work/repair_r1_full_01.jsonl
```

| Run | Guard outcome | Guard wall seconds | Observed cgroup peak bytes | Script runtime seconds | Peak process RSS KiB |
|---|---|---:|---:|---:|---:|
| `repair_r1_tests_final_01` | Admitted; `exitCode: 0`, `limitsVerified: true`, child `exitCode: 0`, `guardReason: null` | 20.492033766997338 | 78630912 | 19.954204236004443 | 79856 |
| `repair_r1_full_01` | Admitted; `exitCode: 0`, `limitsVerified: true`, child `exitCode: 0`, `guardReason: null` | 7726.7241978039965 | 1152303104 | 7720.600736742999 | 1035920 |

Both ran singly with a 53687091200-byte memory cap, zero swap, 32-task cap, one logical CPU and no wall/CPU/inactivity deadline. Terminal resource samples record swap 0 and memory-limit/OOM event counters 0. Guard outcomes, limits and resource samples remain in `_scratch/s11c/s16_cmp_guard/{repair_r1_tests_final_01,repair_r1_full_01}/`. No Wolfram kernel or accepted engine was run.

- Tests: [output](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r1_tests_final_01.json), SHA256 `e27719e5b4a79b307d17ac9efbbfb223a4e10e220f9704245d801120847d3cbe`; all 17 synthetic control groups completed, including individual K1–K3 ablations that each trigger their control's assertion.
- Full comparison: [output](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r1_full_01.jsonl), SHA256 `2cab8b06e674316e83049557191628eaf963eecfe8b9676ccb5558b950f80c78`.
- Final comparator SHA256 `6a55ae5e31e0b284b2beb7d94babc3759637819c802465fd775c21fd9195b5cc`; tests SHA256 `88df585fc55141a91e28016ba279e6835c587ee3176220cdf6005eb6526a6aff`. Both were byte-compared after the run with the copies saved before launch in `_scratch/s16/cmp_work/repair_r1_final_code/`.
- Accepted transcript SHA256: SymPy `da72cd98811b03a85d766743e0abad0df80b979f7b2d26cc0dfe230f8f5d1435`; Wolfram `72d352df5f631ae5d71772ec992306fd5e1eb0c83843a69f0d65fae526699567`. Sources and transcript paths were checked unchanged against `1d990ee6`.

The full output's first record prints the complete source-cited tables: 47 rules (`RULES`), 2028 generated bindings (`BINDING_TABLE`) and 91 names (`name_contract`) in `scripts/S16_cross_engine_comparator.py`. All 2028 bindings were compared; none was not formed or not compared. Per-tag accounting covers 26028 Py / 37578 WL tags, with 238488 / 580543 parsed leaves and 161466 / 212319 compared leaves. Bound-end parsed leaves are 161466 / 213171. All 426 bindings below their parsed count are `coverage:1602` and `comparison_notes:1603` through `comparison_notes:2027`: their message is compared and the WL `Missing` constructor/reason metadata stays unpaired.

Unjoined tags total 24322 Py / 36021 WL. Every member and unpaired component remains in the output. Grouped reasons (counts include associated native notes/dimensions where applicable):

| Engine | Tags | Reason |
|---|---:|---|
| Py / WL | 24264 / 36014 | Engine-local dimension-sum audits enumerate different expression trees; no term/reference pairing by ordinal (Py 418–445; WL 238–243). |
| Py / WL | 26 / 6 | Py free-radius chart and split centre/radius velocities and normal speeds differ from WL surface total velocities/speeds (Py 114–122; WL 56–58). Surface charts and rate pieces are bound; totals cannot bind their parts. |
| Py | 6 | Dimensions of extra bulk-reference entries in `momentumSurfaceMotion` have no unused WL counterpart (Py 223–224, 398–445; WL 130, 230–243). |
| Py | 21 | `weightedLocalRate`, `afterParts.bulk`, `afterParts.current` for count, three dipoles and three momenta are unwindowed intermediates; WL's counterparts are window-weighted objects (Py 183–185; WL 87–95, 101). |
| Py | 5 | Leibniz operand/derivative/orientation and generic pullback operand/derivative are Py-only diagnostics (Py 123–130; WL 74–95 has different window intermediates). |
| WL | 1 | Explicit `moment.indices` labels occur only in WL; Py puts indices in tags (WL 119; Py 205–207). |

For the decision list's named geometry, the unjoined objects are precisely Py `geometry.chart`, `geometry.centerVelocity`, `geometry.radiusVelocity`, `geometry.centerNormalSpeed`, `geometry.radiusNormalSpeed` and WL `surfaces.velocity`, `surfaces.normalSpeed`, for the differing-object reason above. No shared requested balance, moment, rate or printed derivative piece was left unjoined.

The output prints 45 full test-function profiles: mixed quadratic/cubic/quartic forms, each at the three rational `(x1,x2,x3,t)` points `(1/5,2/7,-1/4,1/6)`, `(-2/9,3/8,4/11,-1/7)`, `(3/10,-1/3,2/5,2/9)`, in the generic, centre-velocity-zero, radius-velocity-zero, all-velocities-zero and count-zero strata. Higher geometry derivatives remain live; `0 < r_in < R`. The count-zero density construction is printed in the configuration and its polynomials in the profiles. `definitions:37` (`N`) and `definitions:38` (the printed count integral) each have nine count-zero samples, all printing both operand values as `0`. Unnormalized comparisons include this stratum; 162 normalized samples are outside their domain. There are 39663 evaluated samples, no unresolved coverage and no numerical integral fallback in the full run. Synthetic controls exercise the 40/60-digit quadrature/error-estimate path.

All 885 exact comparisons retain the as-computed residual and add its resolved form; 228 resolved strings differ from the as-computed strings. The 114 as-computed residuals containing literal `Subs(0,` remain printed; none of the resolved residuals contains that literal. These are output counts, not a residual target or an interpretation of the physics.

Machinery remains adapted from S14a_1 (`de81b9fc`): native AST/InputForm readers and partial-derivative slot translation. The final run re-exercised native reading, bindings, name maps, aliases, domains, dimensions, integral evaluation and per-leaf accounting; the tests re-exercised the existing controls plus K1–K3 controls and individual ablations. Mechanical extraction command: `python3 _scratch/s16/cmp_work/repair_r1_summarize.py > _scratch/s16/cmp_work/repair_r1_summary_stdout.json`; its [summary](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r1_full_01_summary.json) and literal stdout retain the counts, resources, hashes and count-zero sample records.

Round 0 outputs, including `full_05.jsonl`, remain intact. The previous run record and code were copied to `_scratch/s16/cmp_work/repair_r1/` before editing. Open questions: non-author round 1 review and interpretation of residuals belong to the orchestrator/S16 record; neither was performed here.
