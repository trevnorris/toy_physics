# S15 comparator builder report

The builder's final message, copied verbatim from its transcript (`_scratch/s15/cmp_build/resume4.log`, lines 1357–1373). Its earlier interim reports are in `resume1.log`, `resume2.log` and `resume3.log` beside it.

Completed the comparator and synthetic tests. Final guarded tests and the complete comparison both exited 0. The [run record](/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md) contains both exact commands, guard outcomes, output paths and SHA256 hashes.

- The [comparator](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py:287) contains 52 binding rules, 45 name-table rows and 4 lexical-slot rows. The first record of [full output](/var/projects/toy_physics/_scratch/s15/cmp_work/full_07.jsonl) prints these and all 1,616 generated bindings with citations, before comparisons.
- All 1,616 bindings were compared. None was not formed, not compared, or compared below its independently parsed leaf count. No decision-list-named object remains unjoined.
- Per-tag accounting covers 2,213 SymPy occurrences: 1,352 bound and 861 unjoined; all 1,353 Wolfram occurrences are bound.
- Unjoined SymPy occurrences comprise 832 dimension-record fields checking different counterparts, 27 empty-support coefficients without printed Wolfram counterparts, one chart-domain prose line and one free-parameter comparison note without printed counterparts. Every occurrence and its reason is listed.
- Bound parsed/compared leaves are 5,855,423/5,855,423 for SymPy and 2,649,540/2,649,540 for Wolfram. Total parsed leaves are 5,858,040 and 2,649,640 respectively. The remaining leaves are explicitly accounted for: 2,513 under unjoined SymPy tags, 104 unpaired SymPy constructor/key leaves and 100 unpaired Wolfram key leaves.
- All 32,040 samples evaluated, including the previously unresolved 810 samples across 18 bindings. The 363 integral-bearing bindings produced 16,335 integral sample pairs, all evaluated exactly; none required numerical quadrature or remained unresolved.
- Evaluation retains mixed quadratic, cubic and quartic forms with distinct functions and coupled arguments. Points `(x₁,x₂,x₃,t)` are `(3/2,2/7,−1/4,1/6)`, `(7/4,3/8,4/11,−1/7)` and `(2,−1/3,2/5,2/9)`.
- All five strata were covered: generic, stationary `X₂`, stationary radii, both stationary, and stationary `X₁`. Higher time derivatives remain generic; geometry checks cover every sample configuration. Full test-function definitions are printed.
- Exact integral linearity makes the sums tractable. Nested integration now preserves free outer angular parameters and handles algebraic roots exactly. Original operands, samples and as-computed residuals remain printed.
- The [synthetic tests](/var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py) completed 27 control groups and 25 comparison-disable groups. Final sources were hash-verified unchanged between testing and the complete run.
- Reused S16 machinery includes native readers, derivative-slot translation, leaf accounting and integral machinery. Synthetic native-format controls and the complete accepted streams re-verified their application here.
- Eight workers ran with canonical output ordering and per-object progress. The guard verified CPUs 0–7, the 55 GiB cap, zero swap and the task cap.
- Tests: guard wall time **465.8017124659964 seconds**; sampled aggregate peak memory **456,359,936 bytes**. Full comparison: **8395.158274214977 seconds**; sampled aggregate peak **21,581,824,000 bytes**. Neither final run recorded memory-limit or OOM events.

Accepted inputs remain at `3192123b`, the decision list at `91248b67`, and the comparator contract at `e305f6d4`. Prior outputs were preserved. No engine was edited or run, and no commit was made. No residual target was given.
