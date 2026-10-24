# S16 comparator builder report

The builder's final message from repair 2 (round 1 finding K5), copied verbatim from its transcript (`_scratch/s16/cmp_repair2/repair2.log`, lines 43004–43033). The earlier reports are preserved at `98936790` (build) and `de28e059` (repair 1).

Repaired K5 in the [comparator](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py): **38 citation-string literals changed**, checked against `1d990ee6`. Tests are unchanged. The [audit](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_audit_01.json) found **no output differences beyond citations and resource measurements**.

The [rewritten run record](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S16_comparator_run_record.md) contains the exact guarded commands, outcomes, output paths and SHA256 hashes.

- All **17 test groups completed**. Tests, full comparison and audit exited **0**, with verified guard limits and no kills.
- The first [full-output record](/var/projects/toy_physics/_scratch/s16/cmp_work/repair_r2_full_01.jsonl) prints **47 rules, 2028 bindings, 91 names and two surfaces**, from `RULES`, `BINDING_TABLE`, `name_contract` and `SURFACE_TABLE`.
- All **2028 bindings** were compared; none was unformed or not compared.
- Py/WL accounting: **26028 / 37578 tags**, **238488 / 580543 parsed leaves**, **161466 / 212319 compared leaves**.
- The **426** bindings below parsed counts remain `coverage:1602` and `comparison_notes:1603`–`2027`; WL `Missing` constructor/reason metadata stays unpaired.

Unjoined tags remain:

| Reason | Py | WL |
|---|---:|---:|
| Different native dimension-sum audits | 24264 | 36014 |
| Free-radius chart and split velocities/speeds versus surface totals | 26 | 6 |
| Extra bulk-reference dimensions without counterparts | 6 | 0 |
| Unwindowed versus window-weighted intermediates | 21 | 0 |
| Engine-only calculus diagnostics | 5 | 0 |
| Explicit moment-index labels | 0 | 1 |

Every member and reason is printed. The named geometry above remains unjoined because its printed objects differ; no shared requested balance, moment, rate or derivative piece is unjoined.

Evaluation retains mixed quadratic/cubic/quartic profiles at `(1/5,2/7,-1/4,1/6)`, `(-2/9,3/8,4/11,-1/7)` and `(3/10,-1/3,2/5,2/9)`. All five strata remain covered: generic, centre velocity zero, radii velocities zero, both zero, and count zero. There are **39663 evaluated samples**, **162 normalized-domain exclusions**, and **no unresolved coverage or numerical integral fallbacks**.

Existing S14a_1 reader/derivative machinery and S16 comparison machinery were retained and re-exercised.

Script runtimes: tests **19.873062699000002 s**; full comparison **7699.739048800999 s**. Peak RSS: **79972 / 1036104 KiB**.

No residual target was given. No commit was made; reviewed baseline remains `de28e059`. Earlier outputs remain intact. The remaining scoped review belongs to the orchestrator.
