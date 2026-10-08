O2 comparator builder r5 — amendment 1

Baseline: `f8e06242099e0737ce493a6b1fc3a04bb74f387f`; amended directive: `44f78617`. No commits, agents, engine executions or reviews. Earlier scratch directories/reports unchanged.
Files C/T/A: `research/pde_ledger_v3/scripts/{O2_cross_engine_comparator.py,test_O2_cross_engine_comparator.py,ablate_O2_cross_engine_comparator.py}` respectively. [Complete measurements](O2_comparator_builder_r5_measurements.txt). B = `_scratch/s11c/o2-comparator-build-r5/`.

1. Named inventory/leaf policy: C:1055 handles OPEN descriptors before structural roles, excludes provenance and inactive operator spelling, inventories the native descriptor operand, and consumes action-head leaves at any depth. C:995 states the disjoint leaf policy. T:881,897,915 checks exact named inventories and exact limit counts through run().
2. Balance orientations: T:834 asserts each balance entry’s printed sign, inventory sign and difference; T:853–878 covers plain sums, Derivative, Lambda, native sums, variation, Function and inactive Total/Map/D. A:81,107 supplies separate term_orientation knives; the existing action-path knives and closed-part control remain.
3. Limit: C:988 explicitly names held derivative/variation variables, aggregate sets/binders, coefficient magnitudes, provenance, options and other syntax outside the three inventories. T:932 confirms the wrapper-coordinate mutation leaves those differences unchanged and that the limit names the excluded content. No new comparison was added.

Reused the HEAD readers, declared tables, canonical live keys, scalar algebra/subtraction, orientation calculation, transpose and accounting. Re-verified with native-serialized controls and the full-stream accounting/limit census; the readers retain the S11c-b framing adaptation.
Tables: B/catalog.json and measurements give every one of the 160 joins, 50 names, 45 actions and 4 balances with both construction citations and spec objects. B/accounting.jsonl and measurements give every parsed/compared count and every unjoined path/reason.
Accounting: 160 joined; 204 unjoined. Parsed/accounted/compared leaves: PY 1857811/1857811/4751; WL 319397/319397/3959. Unaccounted: 0/0; generic unjoined reasons: 0.
Output census: no `wl::6`, `wl::3,8`, `wl::D` or `wl::Map` in named inventories; all occurrence and balance-entry limit counts conserve parsed leaves and are nonnegative. Diagnostic raw syntax remains printed.
Commands: G = `python scripts/s11c_guarded_run.py --pool o2-comparator --memory-gib 6`; S = `research/pde_ledger_v3/scripts/`. No time limits; expanded commands, guard receipts and literal outputs are in measurements.
Tests: `G --log-directory B/tests-02 -- /usr/bin/time -f 'elapsed_seconds=%e peak_rss_kib=%M' -o B/tests-02-resources.txt python S/test_O2_cross_engine_comparator.py`.
Ablations: `G --log-directory B/ablations-01 -- python S/ablate_O2_cross_engine_comparator.py --output B/ablation-copies`.
Full: `G --log-directory B/full-01 -- python S/O2_cross_engine_comparator.py --output B/comparison.jsonl --accounting B/accounting.jsonl --catalog B/catalog.json` (after datalad get of both inputs).
Outcomes/resources: 83 tests, exit 0 (elapsed_seconds=2.22 peak_rss_kib=56704); 47 ablations, every required assertion failure observed, no test errors (child 22.275778204901144 s; maximum copy RSS 56428 KiB); full exit 0 (176.45849719410762 s; peak RSS 1333076 KiB).
Own test correction: tests-01 had two mistaken leaf-total expectations counting association keys as terminal leaves. Corrected 4→2 and 3+extra→2+extra, retaining exact equalities and all exclusion/sign assertions; tests-02 passed. No kill or admission refusal. No residual target was given.

New balance-entry ablations → observed failing tests (T file above; complete 47-variant map and literal failures in measurements):

| Ablation | Failing test(s) |
|---|---|
| `balance_negative_number` | `test_balance_plain_sum_orientation` (T:853) |
| `balance_inactive_body` | `test_balance_inactive_total_orientation` (T:872), `test_balance_inactive_map_orientation` (T:875), `test_balance_inactive_d_orientation` (T:878) |
| `balance_derivative_body` | `test_balance_derivative_body_orientation` (T:857) |
| `balance_product_sign` | `test_balance_plain_sum_orientation` (T:853) |
| `balance_lambda_body` | `test_balance_lambda_body_orientation` (T:860) |
| `balance_variation_body` | `test_balance_variation_body_orientation` (T:866) |
| `balance_function_body` | `test_balance_function_body_orientation` (T:869) |
| `balance_native_sum_body` | `test_balance_native_sum_body_orientation` (T:863) |
| `balance_inactive_total` | `test_balance_inactive_total_orientation` (T:872) |
| `balance_inactive_map` | `test_balance_inactive_map_orientation` (T:875) |
| `balance_inactive_d` | `test_balance_inactive_d_orientation` (T:878) |

Open questions: production interpretation belongs to sub-step 7; review remains with the orchestrator. Empty OPEN differences cover only the three declared inventories. Stopped.
