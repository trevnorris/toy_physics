# S15 comparator run record — targeted control revision

The control revision completes 27 synthetic control groups and 4800 targeted ablations through the guard. This replaces the earlier global-disable evidence, which did not establish that each targeted comparison step was checked. The reviewed baseline is `0ad8773e`, **reviewed, not accepted**. Only the tests and this record changed in this revision; the comparator and guard are byte-identical to the complete `full_07` run. Accordingly, the complete comparison was not rerun. All 1616 bindings and 32040 samples in that retained run remain accounted for; no binding is not formed, not compared, or below its parsed leaf count, and no sample coverage is unresolved. No residual target was given.

Inputs: accepted engine sources/transcripts `3192123b`; S15 decision list revision 3 `91248b67`; comparator contract `e305f6d4`. Adapted S16 machinery at `f08db34a`: native readers, derivative-slot translation, leaf census and integral evaluation. Synthetic notation/mutation controls and complete parsing, binding and accounting on both accepted streams re-verified that machinery. No engine was edited or executed and no commit was made.

The unchanged comparator evaluates sums by exact integral linearity, retains outer parameters during nested integration, and handles residual square roots with exact algebraic relations without choosing a root sign or changing derivative points. Sparse polynomial moments and per-binding cache clearing control memory. Eight workers preserve binding-table output order; synthetic tests compare mathematical output and accounting across one, two and eight workers and reversed submission. Original printed operands and as-computed exact residuals remain in the output. No supplied balance equation is imposed.

## Commands and outcomes

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15-comparator --all-cpus-exclusive --memory-gib 55 --log-directory /var/projects/toy_physics/_scratch/s11c/s15_cmp_guard/tests_20261025_15 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py --output /var/projects/toy_physics/_scratch/s15/cmp_work/tests_15.json
```

Tests: guard exit `0`, `limitsVerified: true`; child exit `0`, `guardReason: null`, empty stderr. Twenty-seven synthetic groups completed and 4800 individual targeted ablations reached their specified assertion. Script-reported runtime `816.9317231979803` s; guard wall `817.9979908599926` s; aggregate sampled peak memory `521035776` bytes. Output `_scratch/s15/cmp_work/tests_15.json`; sha256 `c0a390c8d6ede7f6c53bac2cb4103ba64773b2a93165967f6cf25dcdaf0d1077`. The independent groups use eight workers; the scheduling control then checks one, two and eight comparator workers and reverse submission.

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15-comparator --all-cpus-exclusive --memory-gib 55 --log-directory /var/projects/toy_physics/_scratch/s11c/s15_cmp_guard/full_20261025_07 -- python3 -B /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py --workers 8 --py /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_two_drain_region.out --wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out --output /var/projects/toy_physics/_scratch/s15/cmp_work/full_07.jsonl
```

Retained complete run with the unchanged final comparator: guard exit `0`, `limitsVerified: true`; child exit `0`, `guardReason: null`. Script runtime `8393.043212330987` s; guard wall `8395.158274214977` s; aggregate sampled peak memory `21581824000` bytes. Both jobs verified CPUs 0–7, the 55 GiB aggregate cap, zero swap and 32 tasks; no memory-limit or OOM events were recorded. Output `_scratch/s15/cmp_work/full_07.jsonl`; sha256 `46685b1326ec2e5d731c5518283ca845fec0154080358cf93d0579467f17ef50`.

## Tables, accounting and coverage

`RULE_NAMES`, `RULES` and `binding_contract` in `scripts/S15_cross_engine_comparator.py` define 52 rules and generate 1616 bindings. `NAME_TABLE` has 45 rows; `SLOT_TABLE` has four lexical-slot rows. The first output record prints all these tables and citations in full, before comparisons. Every binding was compared; no decision-list-named object is left unjoined. The following 45 records print all test functions and points.

Per-tag accounting includes all 2213 SymPy and 1353 Wolfram occurrences: 1352 SymPy tags bound, 861 unjoined; all 1353 Wolfram tags bound. The unjoined SymPy reasons are: 832 argument-dimension record fields check arguments against signatures/bounds, unlike Wolfram's field/object signatures; 27 empty-support coefficients have no counterpart among Wolfram's nonempty groups and their sum; one chart-domain prose line corresponds to internally used Wolfram assumptions with no printed prose; one free-parameter comparison note has no counterpart beside Wolfram's plain parameter emit. Every identity, raw payload and source-specific reason is printed in the per-tag records.

All-stream parsed/compared leaves are `5858040/5855423` (SymPy), `2649640/2649540` (Wolfram). Bound-end parsed/compared leaves are `5855423/5855423`, `2649540/2649540`; no binding has a shortfall. The remaining SymPy leaves are 2513 under unjoined tags and 104 explicitly unpaired constructor/key leaves under bound tags. Wolfram has 100 explicitly unpaired key leaves. Every such leaf and its reason appears in the output; none is silently omitted.

All 712 sampled bindings have the complete 3-form × 5-stratum × 3-point grid. Forms: mixed quadratic, cubic and quartic. Strata: generic, X2 velocity zero, both radii velocities zero, both of those zero, and X1 velocity zero, with higher derivatives retained. Points `(x_1,x_2,x_3,t)` are `(3/2,2/7,-1/4,1/6)`, `(7/4,3/8,4/11,-1/7)`, `(2,-1/3,2/5,2/9)`. Tests check the domain and shell at every point/form/stratum. The 363 integral-bearing bindings have 16335 evaluated sample pairs, all exact; no numerical quadrature or unresolved integral coverage appears in the full-run inventory. The numerical fallback remains covered by synthetic tests.

JSON/count/hash audits of the complete output remain `_scratch/s15/cmp_work/full_07.audit.json` and `full_07.coverage-audit.json`. They enumerate the repaired 810 samples, verify complete grids and table order, and retain the unpaired-leaf inventory. The new standard-library audit is `python3 _scratch/s15/cmp_work/review_controls/audit_final.py`, output `review_controls/final_audit.json`. It verifies the final test pin, unchanged comparator and guard hashes, complete-output hash, 4800 consecutive ablation identities, all 52 binding rules and all 45 name rows, and 4012 content-addressed native fixture files. It performs no CAS computation.

## Targeted controls

Each numeric mutation is checked separately for symbolic residual and sample changes at its binding. The test-only function replacements alter one production operation and one binding, while the native reader, extraction, binding, counterpart and remaining comparison channels continue. Exact-residual faults retain samples; sampling faults retain exact residuals and the other operand's samples. Other bindings' comparisons are checked unchanged. Only the intended assertion text counts as rejection: an inactive baseline, unrelated assertion, grammar error or missing key cannot earn ablation credit. Missing/repeated objects and malformed structures specifically target formation diagnostics. No global `disable_comparison` option is used.

The sparse oracle expands the synthetic polynomial densely, expresses sine/cosine monomials as finite Fourier sums, and integrates each exponential mode at the interval endpoints. It never calls the comparator's moment routine or uses its beta/gamma formula. For each of the pressure, volume and rational-cancellation fixtures, the production `angular_moment(..., full=True)` is then halved with evaluator caches cleared. The native comparison fails its independent-oracle equality assertion, while the oracle operand and symbolic residual remain unchanged.

The following fault names specify the altered step. The JSON output's `targeted_ablations` array contains one record for **every** case: control group, consecutive instance number, binding/rule and semantic paths where applicable, native fixture paths, configuration/name/slot changes, exact source replacement and the assertion reached. Thus loop instances (including every derivative argument/context, binding rule and name row) remain separately countable. The per-group index below points into that exhaustive record.

| Fault | Targeted step and alteration, with other comparisons retained |
| --- | --- |
| `exact` | Symbolic subtraction in `compare_value`: subtract the SymPy operand from itself for the selected binding. |
| `sample_py`, `sample_wl` | Evaluate the counterpart in place of the mutated operand at sampling only; exact subtraction stays active. |
| `sample_offset` | Alter one computed sample operand by adding one (Boolean samples are negated); sample residual is recomputed. Checks notation/oracle equivalence and scheduling invariance. |
| `text` | Suppress only the label/basis text difference; operands remain printed. |
| `text_equality` | Emit the side-by-side difference for equal text too; the equal-label comparison changes. |
| `dimension` | Use the same dimension multiset on both sides of the residual; dimension evaluation remains active. |
| `dimension_sample` | Evaluate the Wolfram dimensions in the SymPy sample slot; dimension multiset residual remains active. |
| `dimension_offset` | Double only the SymPy dimension sample values; retain the dimension residual. |
| `formation` | Omit only the `not_formed` diagnostic in the unsupported-object branch; well-formed comparisons continue. |
| `domain` | Replace evaluated domain conditions by true only in the sampling guard; exact comparison and printed domain context remain. |
| `leaf_credit` | Suppress crediting compared leaf paths while retaining parsing, residuals and samples. |
| `unpaired_leaf` | Suppress unpaired-component enumeration; keep the parsed census and comparisons. |
| `parsed_census` | Report one leaf per parsed subtree; retain comparison and unpaired-component enumeration. |
| `record_ordinal` | Assign repeated records to sorted identities in their stream order, instead of their content/source context; both native formats still pass through the production reader and comparator. |
| `duplicate_validation` | Remove only the relevant duplicate end/name/slot/rule rejection in `validate_tables`. |
| `coverage_reduction` | Run each witness with its own point/form/argument-coupling/vanishing/stratum/function-collision reduction alone; the full-coverage witness assertion fails and symbolic subtraction remains unchanged. |
| `full_range_moment_halved` | Halve only production full-range angular moments, clear caches and compare against the independent serialized Fourier oracle. |

| Control | Ablation instances | Faults (case counts) |
| --- | --- | --- |
| `exclusive_guard_controls` | Ancillary resource control | Direct affinity/memory/swap predicate mutations and shared/exclusive lock attempts; enclosing guard verifies the actual allocation. Does not target comparator operations. |
| `nested_radial_derivative` | 1–135 | `sample_offset` (45); `exact` (45); `sample_py` (45) |
| `nested_integral_parameters` | 136–146 | `sample_offset` (3); `exact` (4); `sample_py` (4) |
| `split_integral_evaluation` | 147–155 | `sample_offset` (3); `exact` (3); `sample_py` (3) |
| `sparse_integral_oracle` | 156–167 | `sample_offset` (3); `full_range_moment_halved` (3); `exact` (3); `sample_py` (3) |
| `parallel_schedule` | 168–176 | `sample_offset` (3); `exact` (3); `sample_py` (3) |
| `long_native_payload` | 177–179 | `sample_offset` (1); `exact` (1); `sample_py` (1) |
| `basic_mutations` | 180–190 | `exact` (5); `sample_py` (5); `formation` (1) |
| `labels_and_conditions` | 191–206 | `text` (10); `formation` (1); `exact` (2); `sample_py` (2); `domain` (1) |
| `equivalences` | 207–226 | `sample_offset` (18); `dimension_offset` (1); `text_equality` (1) |
| `derivative_points` | 227–2862 | `exact` (1318); `sample_py` (1318) |
| `mixed_second_derivative_points` | 2863–3502 | `exact` (320); `sample_py` (320) |
| `partial_time_semantics` | 3503–3511 | `sample_offset` (3); `exact` (3); `sample_py` (3) |
| `integral_changes` | 3512–3519 | `exact` (4); `sample_py` (4) |
| `binding_repoints` | 3520–3622 | `text` (7); `exact` (45); `sample_py` (45); `dimension` (3); `dimension_sample` (3) |
| `name_repoints` | 3623–3712 | `exact` (45); `sample_wl` (45) |
| `bound_name_repoints` | 3713–3721 | `exact` (4); `sample_wl` (4); `duplicate_validation` (1) |
| `pieces_groups_and_components` | 3722–4349 | `exact` (314); `sample_py` (314) |
| `packaging_and_accounting` | 4350–4364 | `exact` (5); `sample_py` (3); `sample_wl` (2); `formation` (2); `unpaired_leaf` (1); `parsed_census` (1); `leaf_credit` (1) |
| `equation_relation_kind` | 4365–4365 | `formation` (1) |
| `repeated_records` | 4366–4373 | `dimension` (1); `dimension_sample` (1); `formation` (2); `record_ordinal` (2); `exact` (1); `sample_py` (1) |
| `supplied_equation_and_alias` | 4374–4377 | `exact` (2); `sample_py` (2) |
| `coverage_controls` | 4378–4388 | `coverage_reduction` (11) |
| `test_functions` | Ancillary profile invariants | Every declared profile/form/stratum/point is checked for distinctness, argument dependence, derivative separation and geometry. Corresponding comparison witnesses are ablated individually under `coverage_controls`. |
| `duplicates` | 4389–4393 | `duplicate_validation` (4); `formation` (1) |
| `dag_sampling` | 4394–4798 | `sample_offset` (135); `exact` (135); `sample_py` (135) |
| `quadrature` | 4799–4800 | `exact` (1); `sample_py` (1) |

## Provenance and changed files

Final test source sha256 `95cd49fd5127dcb3aa243e1793d23e6dcc5b2ff12b888bbe0d3c88a4b0eb88ae`; pinned copy `_scratch/s15/cmp_work/review_controls/tests_15.pinned.py`. Comparator sha256 `786ff802b0bcf54b1e8103f596bc7af3216eafc800a524106e23ee293fdd11aa` and guard sha256 `86be9b55a56f7a3538cdf3852b74562bdc8c5bab19c9ad894cbc2c3b6bbe257b` equal their `full07` pins and the reviewed baseline. The complete output hash above was recomputed.

Changed repository files in this revision: `scripts/test_S15_cross_engine_comparator.py` and this record (both paths here relative to `research/pde_ledger_v3/`). No other repository file was changed. In particular, the earlier modification to repository-root `scripts/s11c_guarded_run.py` is already in `0ad8773e`; it is explicitly disclosed here and was not changed again. The pre-existing `.claude/skills/prep-for-compact/SKILL.md` modification was left untouched. No engine was edited or run; no commit was made.

All prior outputs remain preserved. `tests_13` and `tests_14` exited 1 on test-harness bugs (wrong dimension fault operand; an intentionally malformed scheduling sibling treated as an active mutation baseline), not guard refusals or kills. The corrected scheduling diagnostic exited 0, followed by the complete final `tests_15` run above. Their outputs/logs and the final working-file inventory are retained under `_scratch/s15/cmp_work/review_controls/` and `_scratch/s11c/s15_cmp_guard/`; the inventory is `review_controls/files_written.json`. No points, forms, strata, test-function class or measured comparison coverage changed.
