O2 comparator builder r3 — bounded builder report

Deliverable: repair OPEN object deltas and result-path controls, preserve items 1–8, and supply guarded synthetic/ablation/production evidence.
Baseline commit: `f0be35c5e3b9b64549b9c5402420469f60b87d5d`; no new commits, engine executions, spawned agents or review rounds.
Files: `../scripts/O2_cross_engine_comparator.py` (C), `../scripts/test_O2_cross_engine_comparator.py` (T), `../scripts/ablate_O2_cross_engine_comparator.py` (A); [literal measurements](O2_comparator_builder_r3_measurements.txt).
Run directory B: `_scratch/s11c/o2-comparator-build-r3/`; earlier round directories/reports were not changed.
Changed: canonical profile/derivative keys, alpha-bound placeholders and object inventories; lossless raw trees/censuses remain separate. No physics or closure was added; no residual target was given.
Reused: HEAD's lossless readers, cited join/name/action/balance tables, subtraction, signed traversal, layout and accounting; the reader framing originally adapts S11c-b. Synthetic reader/subtraction/orientation/layout controls and production coverage were re-verified; no other comparator residual machinery was imported.
Tables: B/`catalog.json` and the measurements reproduce all 159 join, 47 name, 45 action and 4 balance rows, each with both engine construction-line citations and spec object.
Accounting: B/`accounting.jsonl` and the measurements list every row's parsed/compared leaves and all unjoined paths/reasons: {'joined': 159, 'unjoined': 202}.
Stream census (parsed/accounted/compared; arithmetic coverage, not agreement): PY 1857811/1857811/4751; WL 319397/319397/3959. Unaccounted PY/WL: 0/0; generic unjoined reasons: 0.
Commands: G = `python scripts/s11c_guarded_run.py --pool o2-comparator --memory-gib 6`; S = `research/pde_ledger_v3/scripts/`; paths B/S below expand to those directories. No time limits.
Tests: `G --log-directory B/tests-02 -- /usr/bin/time -f 'elapsed_seconds=%e peak_rss_kib=%M' -o B/tests-02-resources.txt python S/test_O2_cross_engine_comparator.py`; 54 tests completed, exit 0; elapsed_seconds=1.37 peak_rss_kib=55472.
Ablations: `G --log-directory B/ablations-01 -- python S/ablate_O2_cross_engine_comparator.py --output B/ablation-copies`; 20 isolated copies, each complete suite, all required assertion failures observed, no test errors; guard exit 0, child wall seconds 29.831767752068117, maximum per-copy peak RSS 55464 KiB. Each copy's exact failures/runtime/memory are in measurements and B/`ablation-copies/summary.json`.
Production: `G --log-directory B/full-01 -- python S/O2_cross_engine_comparator.py --output B/comparison.jsonl --accounting B/accounting.jsonl --catalog B/catalog.json`; exit 0, runtime 313.1616699050646 seconds, peak RSS 1332660 KiB. All comparison products remain in B/`comparison.jsonl`; this report assigns no interpretation.
No unintended test failure, kill or admission refusal occurred; no assertion was weakened. Initial tests-01 had 53 tests; the final suite added spelling/repetition coverage. Controlled ablation failures are intentional.

Item → implementation (C/T mean the filenames above, under S); ablation → failing test (full failure lists in measurements):

1. Same-role joins/accounting — C:433,584,1198,1271; `same_role_join` → T:296 `test_mass_loss_join_keeps_rhs_separate`.
2. Signed held/body traversal — C:777,1062; `plain_sign` → T:625 `test_plain_sum_orientation`; `held_aggregate_sign` → T:256 `test_held_aggregate_orientation`; `sympy_derivative_body_sign` → T:394 `test_all_held_linear_orientations`.
3. Computed controls/subtraction — T:22,38,235; C:630,1150; `applied_functions_as_symbols` → T:223 `test_exact_argument_content`; `addition_instead_of_subtraction` → T:229 `test_common_translation_preserves_residual`; `constant_residual` → T:70 `test_one_sided_operand`.
4. OPEN name consistency — C:282,612,984; `operand_binding_removed` → T:327 `test_open_operand_binding_consistency`.
5. Every balance entry/closed sum — C:1044,1084,1123; `closed_part_omitted` → T:349 `test_closed_balance_entry`.
6. Derivative points — C:630,701; `evaluation_point_collapsed` → T:379 `test_derivative_evaluation_points`.
7. Computed cited-table OPEN fields/repoints — C:526,555,860,984,1018,1084; `open_field_deltas_omitted` → T:434 `test_action_field_differences`; `action_repoint_ignored` → T:476 `test_each_action_binding_repoint` (all rows); balance repoints T:496.
8. Declared transpose — C:1319; `transpose_dropped` → T:246 `test_transpose_layout_matches_explicit_transpose`.
9. Canonical live objects/binders/operand sets — C:864,869,984,1004; `canonical_derivative_key_removed` → T:564 `test_profile_derivative_is_one_object_key`; `canonical_sqrt_removed` → T:537 `test_identical_open_objects_have_empty_deltas`; `object_sets_as_occurrence_counts` → T:575 `test_object_inventory_ignores_repeated_spelling`. Bilateral freeze/order/repoint and alpha-binder controls: T:582,612.
10. Both relation operands/all closed terms — C:1175,1114; `relation_left_only`, `relation_right_only` → T:634 `test_each_relation_operand_both_engines`; `closed_part_first_term` → T:652 `test_each_closed_term_in_multiterm_balance`.
Lossless-reader preservation also ablated: `function_metadata_lost` → T:170 `test_lossless_function_metadata_survives`; all fixtures enter `run()` with native serialized streams, including malformed-input and table-rejection tests.

Open questions: interpretation of measured differences belongs to sub-step 7; non-author build review remains with the orchestrator. No new implementation prerequisite or physical premise was opened. Stop.
