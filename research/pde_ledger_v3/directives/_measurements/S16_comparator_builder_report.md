# S16 comparator builder report

The builder's final message from repair 1 (round 0 findings K1–K4), copied verbatim from its transcript (`_scratch/s16/cmp_repair1/repair.log`, lines 38117–38147). The round 0 build's report is preserved at `98936790`.

Repaired K1–K4 in the [comparator](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py) and [tests](/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py): fixed partial-time reading, live zero-count fields, additional wrapper-resolved residuals, and corrected citations.

The [rewritten run record](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S16_comparator_run_record.md) contains the exact guarded commands, outcomes, output paths and SHA256 hashes.

- All **17 synthetic control groups completed**, including individual K1–K3 ablations.
- Tests and full comparison exited **0**, with verified limits and no guard kill.
- All **2028 bindings** were compared; none was unformed or not compared.
- The first [output](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r1_full_01.jsonl) record prints all **47 rules, 2028 bindings and 91 name rows**, defined in `RULES`, `BINDING_TABLE` and `name_contract`.
- Accounting: **26028 / 37578** Py/WL tags; **238488 / 580543** parsed leaves; **161466 / 212319** compared leaves.
- The **426** bindings below their parsed counts are `coverage:1602` and `comparison_notes:1603`–`2027`; only WL `Missing` constructor/reason metadata remains unpaired.

Unjoined tags, grouped by reason:

| Reason | Py | WL |
|---|---:|---:|
| Different native dimension-sum audits | 24264 | 36014 |
| Free-radius chart and split velocities/speeds versus surface totals | 26 | 6 |
| Extra bulk-reference dimensions without counterparts | 6 | 0 |
| Unwindowed versus window-weighted intermediates | 21 | 0 |
| Engine-only calculus diagnostics | 5 | 0 |
| Explicit moment-index labels | 0 | 1 |

Every member and reason remains printed. No shared requested balance, moment, rate or derivative piece was left unjoined.

Evaluation used three printed polynomial forms, three rational points and all five strata. Both `N` and the printed count vanish at every count-zero probe. There are **39663 evaluated samples**, **162 normalized-domain exclusions**, **0 unresolved coverage** and **0 numerical integral fallbacks**. All **885** exact comparisons retain both residual forms.

Reused S14a_1 readers and derivative-slot machinery; re-exercised mappings, aliases, domains, dimensions, integrals and accounting.

Script runtimes: tests **19.954204236004443 s**; full comparison **7720.600736742999 s**. Peak RSS: **79856 / 1035920 KiB**.

No residual target was given. No commit was made; reviewed baseline remains `98936790`, and round 0 outputs remain intact. Open work belongs to the orchestrator: round 1 review and interpretation of the residuals.
