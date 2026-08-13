# Independent review — the S11 SymPy engine after fix round 1

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`
(Codex-written, uncommitted). The brief it was built to:
`research/pde_ledger_v3/directives/_legs/S11_engine_fix_round1_brief.md`.
Physics authority: `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`. `CLAUDE.md` binds.

Background: five runs failed to publish `S11_exports.py` because `run_cell('MAIN',5)` raised `MemoryError`
inside expression-tree row reduction. This round replaced the rank and Q8-minor constructions.

## ⛔⛔ A FORM ABLATION IS MANDATORY. ⛔ NOT OPTIONAL.

Change the **structure** of a load-bearing object — flip a sign *and* an off-diagonal, collapse two
independent symbols into one — re-run, and report the **literal** diff. ⭐ A **coefficient** rescale tests
arithmetic; only a **FORM** change tests physics, because scaling never leaves the family.
⚠ Measured 2026-08-04: a script emitted a physics conclusion as a typed sentence with no computation behind
it, and the ablation produced **byte-identical output**. ⛔ **Eight fidelity legs missed it by reading.**
⇒ ⭐ Ask of every claim: **which line computed this?** Give the line number or report it uncomputed.
⚠ An `assert` before an emit hides this — report any that precedes the value it guards.

## What to check

1. ⛔⛔ **IS THE NEW RANK THE SAME ALGEBRAIC OBJECT?** The brief named Q4's exact `rank`/`stacked_rank` of
   the live `M_r` and deliberately did **not** name a routine. The builder chose `DomainMatrix` over a
   rational-function domain. **Verify it returns the exact algebraic rank of the same matrix**, not a
   generic rank at a specialization, not a rank over a different domain, not any cheaper object that
   returns an integer. This is the defect the round was most at risk of.
2. ⛔ **VERIFY `D2/D3/D4 changed 0` YOURSELF.** The builder claims no value moved at those dimensions.
   Recompute and say. ⚠ If a value *did* move, that is a **finding**, ⛔ not a regression — do not treat
   prior output as a reference.
3. ⛔⛔ **THE MID-LOOP PUBLISH IS A SUSPECTED BRIEF VIOLATION — measure it.** There are now two
   `write_exports` sites (`:1731` mid-loop at MAIN completion, `:1772` end of run), while `:1713` still
   computes `skipped_pairs_list = declared_pairs - completed_pairs`. The brief required that *no row claims
   a run record it does not have*, and said to report the conflict under §10 rather than publish. **Does
   the mid-loop export contain a `skipped_pairs` row asserting later packages were skipped when they had
   not yet been attempted?** Show the row.
4. ⭐ **THE THREE-VALUED ZERO PREDICATE.** It now returns `None` when undecided. Does SymPy's row reduction
   handle `None` the way the builder assumes, and can an undecided entry still be taken as a pivot?
5. ⭐ **THE §10 ENTRIES.** The D5 run records *"solve unavailable for 30-equation exact system"* and
   *"canonical locus unavailable"*. Are those honest unavailable constructions, or a computation quietly
   degraded to a token that a consumer would read as a result?
6. ⭐ **THE `out/` TEE** (`:1807-1809`). The builder smoke-tested it against `/tmp`, **not** the repo path.
   Verify it writes the real file without altering stdout content or ordering.
7. **Anything else**, including whether the fix is sound at all.

## Method

- **Every claim carries the command that produced it and its literal output.** A prose re-derivation is
  discarded (`CLAUDE.md` rule 2). Write your own probe scripts and report their absolute paths.
- Report a finding only if it catches a way the physics could be wrong, or a wrong result could survive.

## Operational constraints — these bind both legs identically

- ⛔ **Do not run the full package loop.** It takes hours and has OOM-killed this machine. Use
  `/home/trevnorris/.s11_build/repro_d5.py`, which runs one cell in ~minutes.
- ⛔ Copy the script to `/tmp/<yourname>/` and ablate the **copy**. ⛔ Never modify the working tree.
- ⛔ Wrap every run in `timeout 600`. A 600 s hit is a failed ablation — report it and move on.
- ⚠ A memory watchdog kills any S11 python over 6 GB RSS and logs to `~/.s11_build/watchdog.log`.

## Deliverable

A verdict per numbered item, the literal output of your form ablation, every defect with file:line and the
command that found it.
