# Review a large change to a cross-engine checking harness

⛔ **READ-ONLY on the repository.** Do not edit, create, or commit any file under `research/`. Scratch
files go under `/tmp` only; copy files there if you need to modify them.
⛔ **Do not run Mathematica, `math -script`, or `wolframscript`.**
⚠ Wrap every command in `timeout 900`. On a timeout, make it cheaper — ⛔ never raise the limit.

## The change

Committed as `92461853`:

```
git -C /var/projects/toy_physics show 92461853 -- \
  research/pde_ledger_v3/reduction/engine_output_checks.py \
  research/pde_ledger_v3/reduction/checks_S10.yaml \
  research/pde_ledger_v3/reduction/test_engine_output_checks.py
```

**The directive it was built from:** `research/pde_ledger_v3/directives/S10_harness_repair_directive.md`

## Background

The harness compares the tagged stdout of two independently-written symbolic-algebra engines — a
Mathematica `.wl` and a SymPy `.py` — whose committed output sits in
`research/pde_ledger_v3/mathematica/out/` and `research/pde_ledger_v3/scripts/out/`.

**The defect this change fixes:** the layer that partitions main-package tags into RESPONSIVE and
INVARIANT against control packages — the harness's **ablation** test, the one thing distinguishing a
computed output from a typed one — matched **zero** tags, because it keyed on a naming convention from an
earlier step. It reported `compared=0 responsive=0 invariant=0`, ⛔ **which reads as health.**

## The questions, in priority order

**Q1 — ⭐ Find any remaining check in this file that can report success without having compared
anything.** This is the defect class the change exists to fix. A previous review found several instances.
⭐ **Assume the fix did not catch them all, and assume the fix introduced new ones.** For every layer the
harness reports, answer: *what does its output look like when its input is empty, mis-keyed, or filtered
to nothing, and is that distinguishable from a clean pass?* Name each with line numbers. ⛔ Do not stop at
the first.

**Q2 — Is the new coverage guard honest, or can it be satisfied vacuously?** The change reports
`CONTROL_COVERAGE: main=836 compared=314 uncovered=522 uncovered_fraction=0.62 no_ablation_D=[D2,D5]`.
⭐ **Verify those numbers yourself** from the output files. Then ask: can a layer pass this guard while
still being inert over a region it claims to cover? ⭐ Construct the case if so.

**Q3 — The dimension counter was split into `compared` and `vacuous`.** ⭐ Is the split drawn in the right
place? Does `compared` count only tags where a real dimensional comparison happened? Is anything counted
as compared that was not? ⭐ Verify against actual payloads, not against the code's intent.

**Q4 — Regressions, and ⛔ verify the attributions rather than accepting them.** The S9 configuration must
still work; its recorded baseline is `CONTROL_RESPONSE: compared=170 responsive=150 invariant=20
unparsed=1 unpaired=8`. The unit suite was `2 failed, 15 passed` before this change and is now
`5 failed, 12 passed`. ⭐ **Establish the cause of each new failure yourself.**

**Q5 — Did the change weaken anything?** Guards removed, error paths softened to warnings, a report line
that now says less than it did, a test whose assertion was loosened.

**Q6 — Anything in the change you believe is wrong.**

## ⚠ One known issue that is NOT the subject of this review

The S10 acceptance run currently exits 2 before printing counters, on a nested dimension payload a
follow-up change handles. ⭐ To see the other layers, copy `checks_S10.yaml` to `/tmp` and blank its
`derived_dimensions` mapping.

## Output format

```
VERDICT: <one line>
BLOCKING: <numbered; each with file:line and why it blocks>
Q1 SILENT-ZERO LAYERS: <table: layer, line, what empty input looks like, distinguishable?>
Q2 / Q3 / Q4 / Q5 / Q6: <one short block each>
SOLID: <what you checked and found correct — with what you actually ran>
```

⛔ Do not report the physics conclusions of the step; report on the instrument.
⛔ "Looks fine" without saying what you ran is not a review.
