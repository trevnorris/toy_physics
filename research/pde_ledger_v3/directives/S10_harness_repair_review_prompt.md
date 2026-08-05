# Review a repair directive for a cross-engine checking harness

⛔ **READ-ONLY.** Do not edit, create, or commit any file in the repository. You may run `python3` and
write scratch files under `/tmp` only.
⛔ **Do not run Mathematica, `math -script`, or `wolframscript`** under any circumstances.
⚠ Any script you run: wrap in `timeout 600`. If it times out, make it cheaper — never raise the limit.

## What you are reviewing

`research/pde_ledger_v3/directives/S10_harness_repair_directive.md` — a list of repairs to be applied to
`research/pde_ledger_v3/reduction/engine_output_checks.py` and `reduction/checks_S10.yaml`.

## Context

Two independently-written symbolic-algebra engines — `mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
and `scripts/S10_brane_mode_spectrum_sympy_audit.py` — computed the same physics from one shared written
specification, `directives/S10_SHARED_PHYSICS.md`. Their committed outputs are under
`mathematica/out/` and `scripts/out/`. The harness is the instrument that compares them.

The engines are **finished and frozen**. The directive is meant to change **only the instrument** (with one
narrow exception it names). Its author is the same person who wrote the specification, so it inherits that
author's assumptions.

## The questions, in priority order

**Q1 — ⭐ Find another check in `engine_output_checks.py` that can report success without having compared
anything.** The directive's item H0 exists because one layer's tag-matching pattern silently matched zero
tags and the layer then reported `compared=0` on a line that reads like health. ⭐ **Assume that is not
the only one.** Go through every layer the harness reports and, for each, answer: *what does its output
look like when its input is empty, mis-keyed, or filtered to nothing, and is that distinguishable from a
clean pass?* Name every one you find, with line numbers. ⛔ Do not stop at the first.

**Q2 — ⭐ Item H3 asserts that a dimension table built from ONE package, with the brane dimension left
symbolic, is correct for every package in a sweep over `D = 2,3,4,5`. DERIVE THIS YOURSELF AND CHECK IT.**
⛔ Do not verify it against the directive's own reasoning, and ⛔ do not verify it against the engines'
`_SPECIALIZED` tags, which are that same claim already evaluated.
⭐ Work from the action written in `S10_SHARED_PHYSICS.md` §2: read off what the coefficients' dimensions
must be for the action to be dimensionally homogeneous at general `D`, then determine whether every
homogeneity verdict the harness would issue is independent of `D`. ⭐ Then say whether a *numeric*
per-package table and the *symbolic* table can ever disagree, and construct a case if they can.

**Q3 — item H0b makes "a configured layer compared nothing" an operational failure.** ⭐ Work out what it
would do to the existing `checks_S9.yaml` run and to any legitimate case where a package genuinely has no
counterpart. Does it convert a real signal into noise, or noise into a real signal? Be concrete.

**Q4 — item H2 collapses a list of per-coefficient dimension vectors into the single one they share, and
fails if they differ.** ⭐ Is failing right? Consider what it means for two coefficients in the same family
to carry different dimensions, and whether the physics allows it. ⚠ The directive also flags one known
weakness in its own prototype; ⭐ say whether the fix it asks for is sufficient and find any other input
shape that would be misread.

**Q5 — items H5 and H6 set things aside.** H5 declares a list of sentinel names non-dimensional; H6 sets
aside nine dimensional-homogeneity failures as artefacts of substituting a bare number for a dimensionful
symbol. ⭐ For each, answer: **what real defect would this hide?** ⛔ Be specific and adversarial — name the
concrete failure that would now go unreported, or say plainly that you could not construct one.

**Q6 — the acceptance section demands each repaired layer be shown able to fail by corrupting a copy of an
engine output.** ⭐ Is that demonstration strong enough to establish the layer works, or can a layer pass
it while still being broken? ⭐ If you can construct such a layer, do.

**Q7 — anything the directive misses**, any instruction in it that is wrong or self-contradictory, and
anything it asks for that would damage the harness or the engines.

## Output format

```
VERDICT: <one line>
BLOCKING: <numbered; each with file:line and why it blocks>
Q1 SILENT-ZERO LAYERS: <table: layer, line, what empty input looks like>
Q2 DERIVATION: <your own derivation and its conclusion>
Q3 / Q4 / Q5 / Q6: <one short block each>
NITS: <numbered>
SOLID: <what you checked and found correct — with what you actually ran>
```

⛔ Do not tell me the physics conclusion of the step. ⭐ Report on the instrument and the directive only.
⛔ "Looks fine" without saying what you ran is not a review.
