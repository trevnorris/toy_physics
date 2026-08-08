# Independent review — S9 export chain, fix round 2

## Artifact

Uncommitted changes in `/var/projects/toy_physics` on top of commit `4d12d0b6`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
research/pde_ledger_v3/scripts/extract_knob_inventory.py
research/pde_ledger_v3/scripts/S9_exports.py                                (regenerated)
research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out  (regenerated)
```

Pre-fix state: `git show 4d12d0b6:<path>`. `git diff` is exactly this round. The Wolfram engine was not
changed. Directives: `directives/S9_export_chain_fix_round2.md` (this round),
`directives/S9_export_chain_rebuild_directive.md` (the build). Plan: `S9_REWRITE_PLAN.md`.
**Read the directives last.**

## Context

Round 1 fixed six findings; two legs then found four more, all in the class system and the agreement
count. This round fixes those four. ⛔ **F-2 (the standard-name→object binding) is with the user, is out of
scope, and must not be reported again or worked around.**

## What this round claims

```
mangled separator on an annotation      -> extractor exit 1 with a finding   (was: silent, .out byte-identical)
out-of-set tag                          -> extractor exit 1 with a finding
PY_S9_EXPORT_CLASS_TALLY emitted into the .out, and it moves when a class moves
PREMISE entries now 5: lagrangian, plane_wave_ansatz, assumptions, dim_energy_density, dim_squared_velocity
explicit named fallback replaces the comment on the dict accumulator; 125 entries take it, enumerated
INDEPENDENT_CROSS_ENGINE_AGREEMENT: 9/9  -- SPEED_DIMENSION_DIFFERENCE restored;
  COEFFICIENT_DIMENSION_DIFFERENCE, IMPLIED_SPEED_DIMENSION, DYNAMICAL_MATRIX_ROUTE_RESIDUAL excluded as entailed
17 tags have no fb29bba2 predecessor, enumerated
regression: no exported value moved, no dim moved, only the two intended class moves; WL byte-identical;
  form ablation still moves 86 exported values
```

## Required method

Derive independently where a value is in question, **in CAS**; a prose re-derivation is worth nothing.
Write your script before opening the artifact; save it and its literal stdout to named absolute paths and
report them. **Ablate copies in a scratch directory — never the working tree, and never `git
checkout`/`restore`/`stash`: the artifact is uncommitted and that would destroy it.** For every claim ask
*which line computed this?*

## What to attack

**① Is the structural declaration scan actually structural, or a wider regex?** The round claims a
malformed tag is now caught because declaration sites are enumerated structurally. Find the ways round: an
annotation on a continuation line; a multiple-assignment target; an augmented or annotated assignment; a
declaration inside an `if`/`for`/function that still produces an exported object; a tag that is valid but
on the wrong object; a Unicode look-alike separator. **Does any of them still land an object in the default
class silently?** That is G-1 reopening.

**② Is `PY_S9_EXPORT_CLASS_TALLY` computed, or assembled?** It must be read out of the collection that was
actually exported. Mutate a class and confirm the emitted tally moves. Then check it cannot be satisfied by
a literal beside the computation — this project's signature defect.

**③ The 125-entry default.** The round replaced an unowned comment with a named fallback and enumerated who
takes it. Is the enumeration correct and complete? More importantly: **is `DERIVED` the honest class for all
125?** Look for hand-typed reference values, chosen forms, and anything a control would move. The builder
explicitly flagged one it would not decide — `root_scaling_quadratic`, constructed as
`lambda_scale**2 * root`, a mixed live-root/hand-selected-reference object. **Settle it with evidence:
posited, derived, or genuinely mixed.**

**④ Is `9/9` right?** Verify each of the three exclusions is genuinely *entailed* by rows that are counted —
not merely correlated, and not excluded for reading zero. Then attack the nine: is any of them **also**
entailed by the others, or unable to disagree for a reason not yet named? The rule is that a row is
excluded for entailment and never for reading zero; check it was applied in both directions.

**⑤ Are the two new `PREMISE` entries right, and are five the right number?** `dim_energy_density` and
`dim_squared_velocity` are now premises. Confirm both are hand-typed rather than solved. Then ask whether
anything else among the 132 `DERIVED` is posited — two were found by reading last round, so read.

**⑥ Regression.** Confirm no exported value moved against `4d12d0b6` beyond the two intended class changes,
that the Wolfram engine still re-runs byte-identical, and that a **form** ablation still moves the exported
values — 86 is fewer than round 1's 119 and that difference needs explaining, not assuming.

## Physics filter

> Report a finding only if it catches a way the physics could be wrong, or changes what may be claimed.

**This project's failure mode is adding machinery, not omitting it.** A finding that a piece of this is
unnecessary is worth as much as one that it is insufficient. Propose a new check only if its absence lets a
specific wrong answer through — and name that answer.

A leg returning "nothing survives the filter" is weak evidence alone: state what you checked, what you
could not, and what would have had to be true for you to find something. ⚠ Last round one leg cleared this
work outright while the other found four real defects.

## Operational

One Mathematica kernel at a time (two seats); `timeout 600` per kernel run; a 600 s hit is a failed
ablation — report and move on. Do not modify the working tree. Do not commit. Do not restore `reduction/`.
Absolute paths. Save every ablation script and its literal stdout; report the paths.

## Output

Rank most-severe first: **claim · evidence (literal command and output, or `file:line`) · what must change.**
