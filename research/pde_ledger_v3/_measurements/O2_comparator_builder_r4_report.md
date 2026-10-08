O2 comparator builder r4 — amendment 1

Baseline/amended directive: `44f78617552aa6e1b7ea9ba180485cf9791a1cb6`. No commits, engine executions, agents or reviews launched; earlier outputs/reports unchanged.
Files C/T/A: `research/pde_ledger_v3/scripts/{O2_cross_engine_comparator.py,test_O2_cross_engine_comparator.py,ablate_O2_cross_engine_comparator.py}` respectively. [Complete measurements](O2_comparator_builder_r4_measurements.txt). B = `_scratch/s11c/o2-comparator-build-r4/`.

1. O4: C:413 separates operand/count-label joins; C:1322,1339 records the unmatched repeated identity, relation constructor and hold occurrence. T:836 detects restoration of the wrong container grain; O4 rows now contain no hold-action differences.
2. Scope/limit: C:988,1002,1094,1164 emits exactly three inventories and their limit; repetitions counted on raw local occurrences, outside leaves counted disjointly. Constructor trees/censuses/hashes are labelled serialization diagnostics. OPEN argument positions and multiplicities are not compared.
3. Labels: C:332 binds both carried-momentum labels to their Wolfram counterparts (PY 251–252 / WL 204, spec §5); C:1002 includes text and unbound names. The O4 count label is also bound (PY 396 / WL 340).
4. Controls: T:679,688,718 adds the missing finite-contract controls; A:33 lists one ablation per bullet. Each copy runs its designated native-serialized control(s) through run(); the unmodified full suite runs separately.

Reused HEAD readers, scalar algebra/subtraction, orientation traversal, transpose and accounting; re-verified by synthetic controls and full-stream census. Reader framing remains the S11c-b adaptation. Superseded binder/tree assertions now address live inventories and printed limits under amendment 1.
Tables: B/catalog.json and measurements contain all 160 joins, 50 names, 45 actions and 4 balances, each with both construction citations and spec object. B/accounting.jsonl and measurements contain every parsed/compared count and unjoined path/reason.
Accounting: 160 joined; 204 unjoined. Parsed/accounted/compared leaves: PY 1857811/1857811/4751; WL 319397/319397/3959. Unaccounted leaves: 0/0; generic reasons: 0.
Commands: G = `python scripts/s11c_guarded_run.py --pool o2-comparator --memory-gib 6`; S = `research/pde_ledger_v3/scripts/`. No time limits. Fully expanded commands and literal outputs are in measurements.
Tests: `G --log-directory B/tests-02 -- /usr/bin/time -f 'elapsed_seconds=%e peak_rss_kib=%M' -o B/tests-02-resources.txt python S/test_O2_cross_engine_comparator.py`.
Ablations: `G --log-directory B/ablations-01 -- python S/ablate_O2_cross_engine_comparator.py --output B/ablation-copies`.
Full: `G --log-directory B/full-01 -- python S/O2_cross_engine_comparator.py --output B/comparison.jsonl --accounting B/accounting.jsonl --catalog B/catalog.json`.
Outcomes/resources: 70 tests, exit 0 (elapsed_seconds=1.65 peak_rss_kib=56208); 33 ablations, every designated assertion failure observed, no test errors (child 11.622178299119696 s; maximum copy RSS 55904 KiB); full exit 0 (170.71045624720864 s; peak RSS 1333064 KiB).
No unintended code failure, kill or admission refusal occurred; no assertion was relaxed to repair a failure. No residual target was given. Production differences are left uninterpreted.

Item 6 map, in directive order: R=residuals/operands; J=joins/names; P=profiles/derivatives; O=OPEN; B=balances; L=layout/translation; I=limit; T=output. Each entry gives ablation → failing test (test_ prefix omitted); full names, file:lines and all failures are in measurements.
R1 `subtraction_to_addition` → `residual_reconstructs_left_operand`; R2 `operands_after_residual` → `operands_precede_residual`; R3 `constant_residual` → `one_sided_operand`; R4 `form_ignored` → `form_change`.  
J1 `repoints_ignored` → `each_name_binding_repoint`, `each_join_binding_repoint`; J2 `join_requires_same_tag` → `moved_tag_and_key_stays_joined`; J3 `unsupported_as_zero` → `computed_and_open_same_role`; J4 `named_inventory_declared_only` → `unbound_open_name_is_in_difference`.  
P1 `stripped_argument_restored` → `applied_argument_stripped`; P2 `live_profile_frozen` → `live_profile_frozen`; P3 `derivative_order_ignored` → `derivative_order`; P4 `evaluation_point_collapsed` → `derivative_evaluation_points`; P5 `binder_scope_erased` → `contract_binder_structure`.  
O1 `head_difference_omitted` → `contract_changed_open_head`; O2 `named_difference_omitted` → `contract_changed_named_operand`; O3 `labels_omitted` → `contract_changed_label`; O4 `live_difference_omitted` → `contract_changed_live_object`; O5 `orientation_difference_omitted` → `contract_changed_orientation`; O6 `canonical_key_removed` → `identical_open_objects_have_empty_deltas`; O7 `missing_sibling_as_exact` → `nested_sibling_removed`.  
B1 `plain_sign` → `plain_sum_orientation`; B2 `held_aggregate_sign` → `held_aggregate_orientation`; B3 `sympy_body_sign` → `all_held_linear_orientations`; B4 `closed_part_first_term` → `each_closed_term_in_multiterm_balance`.  
L1 `transpose_dropped` → `transpose_layout_matches_explicit_transpose`; L2 `relation_left_only` → `each_relation_operand_both_engines`; L3 `sqrt_to_cbrt` → `wolfram_sqrt_translation`.  
I1 `repeated_live_omitted` → `limit_repeated_live_object`; I2 `outside_count_omitted` → `limit_outside_leaves`.  
T1 `verdict_token_inserted` → `output_has_no_verdict_tokens`.  
Further controls: `duplicate_tables_admitted` → `duplicate_join_rejected`, `duplicate_name_rejected_both_sides`; `boolean_as_algebra` → `boolean_does_not_hide_algebraic_sibling`; `o4_container_grain` → `o4_components_have_same_role_joins`.

Open questions: interpretation belongs to sub-step 7 and build review remains with the orchestrator. The printed inventories do not establish argument placement, multiplicity or equality of content outside them. Stop.
