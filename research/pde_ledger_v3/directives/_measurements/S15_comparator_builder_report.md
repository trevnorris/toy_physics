# S15 comparator builder report

The builder's final message from repair 1 (round 0 findings C1 and C2), copied verbatim from its transcript (`_scratch/s15/cmp_build/repair1.log`, lines 23976–23992). The build's report is preserved at `0ad8773e`.

Completed the control revision. The [run record](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md) gives exact commands, hashes, and each control’s targeted operation, alteration and detecting assertion.

- **27 synthetic groups and 4,800 targeted ablations completed.** Exact-residual faults retain sampling; sampling faults retain exact comparisons. Unrelated exceptions and inactive baselines earn no credit.
- The sparse-moment oracle now uses independent Fourier integration. It detects halved production full-range moments in all three oracle fixtures.
- Final guarded tests exited **0**: wall time **817.9979908599926 seconds**, sampled aggregate peak **521,035,776 bytes**.
- The comparator source is unchanged. Its hash and the complete `full_07` output hash were verified, so no complete rerun was required.
- That retained run compared **1,616 bindings and 32,040 samples**. No binding was unformed, uncompared or below its parsed leaf count; no sample remained unresolved.
- Its first output record contains **52 binding rules, 1,616 bindings, 45 name rows and four slot rows**, with citations.
- Accounting remains: SymPy **1,352 bound / 861 unjoined**; Wolfram **1,353 bound / 0 unjoined**. The unjoined reasons are 832 different dimension-check counterparts, 27 absent empty-support counterparts, and two unmatched prose/note lines. No decision-list-named object is unjoined.
- Bound parsed/compared leaves remain **5,855,423/5,855,423** and **2,649,540/2,649,540**. Remaining leaves are explicitly accounted for in the output and record.
- All three forms, three rational points and five strata remain unchanged. All **16,335 integral sample pairs** evaluated exactly; none required numerical quadrature or remained unresolved.
- Retained full-run wall time: **8395.158274214977 seconds**; sampled aggregate peak: **21,581,824,000 bytes**. Both final jobs used the verified 55 GiB guard allocation.
- Reused S16 readers, derivative translation, leaf accounting and integral machinery remain covered by the synthetic controls and retained complete-stream accounting.

Only the [tests](/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py) and run record changed in this revision. The earlier change to [scripts/s11c_guarded_run.py](/var/projects/toy_physics/scripts/s11c_guarded_run.py) is in reviewed-but-not-accepted commit `0ad8773e`; it was not changed again. Retained scratch files are enumerated in the [inventory](/var/projects/toy_physics/_scratch/s15/cmp_work/review_controls/files_written.json).

No engine was edited or run. No commit was made. No residual target was given.
