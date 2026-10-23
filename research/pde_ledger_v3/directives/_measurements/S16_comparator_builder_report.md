# S16 comparator builder report

The builder's final message, copied verbatim from its transcript (`_scratch/s16/cmp_build/build_resume_2.log`, lines 1707–1746).

Completed the [comparator](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py), [synthetic tests](/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py), and [run record](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S16_comparator_run_record.md).

`binding_contract` defines 47 rules generating 2028 bindings; `name_contract` defines 91 name rows. All tables and citations appear in the first [full-output](/var/projects/toy_physics/_scratch/s16/cmp_work/full_05.jsonl) record.

All 2028 bindings were compared; none was unformed or not compared. The 426 below-count bindings are `coverage:1602` and `comparison_notes:1603`–`comparison_notes:2027`: WL `Missing` metadata remains unpaired while its message is compared.

| Engine | Tags | Unjoined tags | Parsed leaves | Compared leaves |
|---|---:|---:|---:|---:|
| Py | 26028 | 24322 | 238488 | 161466 |
| WL | 37578 | 36021 | 580543 | 212319 |

Bound-end parsed leaves: Py 161466; WL 213171. Every tag and unpaired component is enumerated in the output.

| Unjoined reason | Py | WL |
|---|---:|---:|
| Dimension audits use different expression-tree terms/reference addends | 24264 | 36014 |
| Free-radial chart or split velocities/speeds versus surface totals | 26 | 6 |
| Unwindowed IBP versus window-weighted intermediates | 21 | 0 |
| Engine-specific Leibniz/pullback diagnostics | 5 | 0 |
| Py-only bulk reference-motion dimensions | 6 | 0 |
| Explicit moment-index labels versus indices encoded in tags | 0 | 1 |

Unpaired bound content includes Py bulk reference-motion entries and alias dimensions, WL solved-law/windowed intermediates, and constructor/key/`Missing` metadata. Source citations and complete membership are retained.

Evaluation used three polynomial forms with live geometry derivatives, at `(1/5,2/7,-1/4,1/6)`, `(-2/9,3/8,4/11,-1/7)`, and `(3/10,-1/3,2/5,2/9)`. Strata: generic, center velocity zero, radius velocities zero, both zero, and count zero. There were 39663 evaluated samples and 162 normalized zero-count exclusions; numerical integral fallbacks and unresolved coverage were both 0.

Reused S14a_1 AST/InputForm readers and derivative-slot machinery (`de81b9fc`); reverified native reading, mappings, aliases, domains, dimensions, integrals, mutation controls, and accounting.

All 13 synthetic control groups completed. Tests and final full run both returned exit code 0, with verified limits and `guardReason: null`. Exact commands, output SHA256 hashes, and receipts are in the run record.

| Run | Guard wall seconds | Observed cgroup peak bytes |
|---|---:|---:|
| Tests | 14.437464013986755 | 76738560 |
| Full | 7616.337004225 | 1114296320 |

The reboot-interrupted output is preserved separately. Final accounting and hashes were validated; no reservation remains. No commit was made. Inputs remain at `1d990ee6`, under gate `5994e9ac`.

No residual target was given. Scientific interpretation and non-author review remain for the S16 record.
hook: Stop
hook: Stop Completed
