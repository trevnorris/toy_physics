# Independent review — S11 engine fix round 2 DECISION LIST, before a builder runs

Read-only. ⛔ Change no file. `CLAUDE.md` rule 7's trigger: no builder launches until its decision list has
had two legs. The list is the one artifact the builder trusts and is otherwise checked zero times.

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_engine_fix_round2_brief.md`
Evidence: `research/pde_ledger_v3/directives/_measurements/S11_engine_fix_round2_brief.md`

## Read the authorities FIRST, then the list

`CLAUDE.md` (rules 2, 3, 5, 7, 10, 11, 14); `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`;
the engine at `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` (baseline `19d4a44a`).

## Context

Round 1 cleared a `MemoryError` wall: `run_cell('MAIN',5)` now completes, the rank construction was
verified by two legs to be the exact algebraic rank of the live matrix, and a FORM ablation showed nothing
is typed. It also introduced one regression and exposed several pre-existing defects. This list is round 2.

## What to check

1. **ARE ALL EIGHT ITEMS REAL, AND CORRECTLY LOCATED?** Verify each at its cited line. Report any that is
   misdiagnosed, already fixed, or whose stated cause is wrong. Round 1's list had its central item wrong,
   and two legs caught it — that is what this gate is for.
2. **ITEM 3 (the regression).** The list says an undecided zero test must not silently discard a candidate,
   and forbids reverting to a structural test. Is that property implementable without inventing the routine
   (rule 3), and is it tight enough to exclude a wrong object?
3. **ITEM 5 (spec compliance).** The list says a hard-coded size gate refuses a computation that completes,
   and cites `SHARED_PHYSICS:245` and `:285`. Verify against the spec. Is removing the gate the right call,
   or is this a spec question rather than an engine question?
4. **ITEM 4 vs ITEM 2.** Both concern the mid-loop publish. Do they conflict, and would applying both leave
   a coherent publish path?
5. **ACCEPTANCE.** Round 1's acceptance was unreachable by its own item 1. Is this one reachable by these
   eight items? Check each acceptance clause against what the items actually change.
6. **RULE 5.** Does the list, or its measurements file, leak an expected value or a prior result a builder
   could converge on? Say so precisely if it does.
7. **WHAT IS MISSING OR WRONG**, including the premise that these are engine fixes rather than spec
   questions.

## Method

- **Every claim carries the command that produced it and its literal output** (rule 2). A prose
  re-derivation is discarded.
- ⛔ Do not run the full package loop — it has OOM-killed this machine. Use
  `/home/trevnorris/.s11_build/repro_d5.py` (one cell, minutes). Wrap runs in `timeout 600`.
- ⛔ Do not edit any file. A watchdog kills any S11 python over 6 GB RSS.
- Report a finding only if it bears on whether the physics could come out wrong, or a wrong result survive.

## Deliverable

A verdict per numbered item, every additional defect with file:line and the command that found it, and a
plain statement of anything the list gets wrong.
