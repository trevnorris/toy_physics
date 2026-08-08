# Independent review — S9 export chain, fix round 1

## Artifact

Uncommitted changes in the working tree at `/var/projects/toy_physics`, on top of commit `753ae7b1`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
research/pde_ledger_v3/scripts/extract_knob_inventory.py
research/pde_ledger_v3/scripts/S9_exports.py                     (regenerated)
research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out  (regenerated)
```

The pre-fix state of each is `git show 753ae7b1:<path>`. `git diff` shows exactly this round's work.
The Wolfram engine was **not** changed in this round.

The fix directive is `research/pde_ledger_v3/directives/S9_export_chain_fix_round1.md`; the build directive
is `directives/S9_export_chain_rebuild_directive.md`; the governing plan is `S9_REWRITE_PLAN.md`.
**Read the directives last, not first.**

## Context

The first build (`753ae7b1`) was reviewed by two legs. Its **physics held** — both legs derived S9
independently in CAS and reproduced every value, and form ablation confirmed the twelve standard-name
emissions are read out of the live computation. Its **export machinery did not**. This round fixes six of
the seven findings.

⛔ **F-2 was deliberately NOT fixed and is not in scope for this review.** It is with the user. For the
record: each engine binds a standard name to an object through its own name→object table, and re-pointing
one entry made the light cone an `omega^2` rather than a speed squared while every check passed. **Do not
report F-2 again, and do not propose a fix for it.**

## What the fix round claims

```
EXPORT_COMPUTED_DIMENSION_COUNT: 3    DIMENSIONED: rho_br, mu_R, transverse_speed_squared
CLASS_TALLY: DERIVED 134, KNOB 2, PREMISE 3   PREMISES: lagrangian, plane_wave_ansatz, assumptions
PY_MAPPED_DIFF_EXIT: 0, and each of the seven formerly-excluded payloads corrupts to exit=1
ANNOTATION_DRIFT: removing mu_R's KNOB annotation moves the export (KNOB 1, LEDGER 138)
CROSS_ENGINE_INFORMATIVE_AGREEMENT: 10/10, excluding the two self-zero rows
DEGENERATE case prints CANDIDATES_FOUND: [] before the AssertionError
WL_BYTE_IDENTICAL_EXIT: 0
```

## Required method

**Derive independently where a value is in question, and do it in CAS.** A prose re-derivation is worth
nothing here. Write your script before opening the artifact, save it and its literal stdout to named
absolute paths, and report those paths.

**Ablate. Copy to a scratch directory and ablate the copy — never the working tree.** ⛔ Never `git
checkout` anything: the artifact under review is **uncommitted** and a checkout destroys it.

Ask of every claim: **which line computed this?** Give the line number or report it as uncomputed.

## What to attack

**① Did F-1 actually close, or did the count just get smaller?** Three entries now carry a `dim`. Verify
each is a dimension the solve produced **for that object**, not a copy of its value and not inferred.
Then check the converse: is there an object whose dimension the run genuinely computes but which now has
**no** `dim`? A silently-dropped dimension is the same defect wearing the other face. In particular the run
still emits per-coefficient dimensions for the control packages — decide whether their absence from the
export is correct or a loss.

**② Is the `dim` on `transverse_speed_squared` the right vector, and is it live?** This is the light cone
and the step rests on it. Derive its dimension yourself. Then mutate the source it is read from and confirm
the exported `dim` moves.

**③ Does class really come from the annotation now, in both directions?** The drift probe removes an
annotation. Also try: change a class from `DERIVED` to `PREMISE` and back; add a new annotated declaration;
give a declaration a tag outside the closed set; leave a declaration untagged. Does the export follow, and
does the extractor report the untagged one as a finding?

**④ Are the three `PREMISE` entries the right three, and are they the only three?** The claim is
`lagrangian`, `plane_wave_ansatz`, `assumptions`. Read the engine and the step record. Is anything else in
the 134 `DERIVED` actually posited rather than derived? Consider the control-package actions and anything
whose form is chosen rather than solved for. **A wrong classification here corrupts what the knit phase
will believe S9 established.**

**⑤ Is the mapped diff now able to fail on everything it claims to cover?** The report shows exit=1 for
seven corrupted payloads. Re-run that yourself. Then check what is still excluded — five telemetry rows —
and confirm each genuinely has no baseline predecessor. Try corrupting an excluded row and say what is lost.
Check the map is not so permissive that it maps two different baseline objects onto one new name.

**⑥ Is `10/10` honest?** Two rows were dropped as self-zero. Verify that both really are each engine's own
internal residual, already zero when emitted, and check whether any of the remaining ten is *also*
uninformative for the same reason — a row that would read zero whatever the physics did.

**⑦ F-7's guard.** Confirm the operands are printed before the guard for **all four** gated objects, not
just the one demonstrated. Report any other assert in the file that precedes the value it guards.

**⑧ Regression: is the physics still what it was?** Both `.out` files and the export changed. Confirm no
computed physics value moved between `753ae7b1` and now, other than the intended `dim`/`class` fields.
Re-run a **form** ablation to confirm the twelve emissions are still live readout and were not accidentally
frozen by this round's restructuring.

## Physics filter

> Report a finding only if it catches a way the physics could be wrong.

**This project's failure mode is adding machinery, not omitting it.** A finding that a piece of this is
unnecessary is worth as much as one that it is insufficient. Do not propose a new check unless its absence
lets a specific wrong answer through — and name that answer.

A leg returning "nothing survives the filter" is weak evidence alone: state what you checked, what you
could not, and what would have had to be true for you to find something.

## Operational

- One Mathematica kernel at a time — the licence has **two seats** and another leg may be running.
  `timeout 600` on every kernel run; a 600 s hit is a failed ablation, report it and move on.
- ⛔ Never modify the working tree. ⛔ Never `git checkout`. ⛔ Do not commit. ⛔ Do not restore `reduction/`.
- Save every ablation script **and** its literal stdout to named absolute paths; report them.
- Absolute paths throughout.

## Output

Rank most-severe first: **claim · evidence (a literal command and its output, or `file:line`) · what must
change.**
