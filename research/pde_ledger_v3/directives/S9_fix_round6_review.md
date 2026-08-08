# Independent review — S9 export chain, fix round 6 (small)

## Artifact

Uncommitted in `/var/projects/toy_physics`, on top of `71f33901`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
research/pde_ledger_v3/scripts/S9_exports.py                                (regenerated)
research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out  (regenerated)
```

`git diff` is exactly this round. Pre-round state: `git show 71f33901:<path>`. The Wolfram engine is
untouched. Directive: `directives/S9_export_chain_fix_round6.md`. Plan: `S9_REWRITE_PLAN.md`.
**Read directives last.** ⛔ **F-2 is with the user, out of scope, not to be re-reported.**

## Context

The export is the single source of truth a later step — ultimately the knit phase — looks things up in.
Three of this round's four items serve that; one is a classification correction.

## What this round claims

```
U-1  LEDGER field order now display, value, dim, class, step; all 139 pre-existing values byte-identical
U-2  assert set(posited_output_classes.values()) <= set(CLASS_TAGS)
     'PREMISEE' -> exit 1;  '' -> exit 1;  unablated run byte-identical
     omission and wrong-but-valid tags REMAIN OPEN (need an independent oracle)
U-3  field_dimension = Matrix([[1],[0],[0]]) and q_dimension = Matrix([[-2],[0],[0]]) enter LEDGER as PREMISE
     new count 141;  KNOB 2, PREMISE 7, DERIVED 132;  round-trip residual 0
U-4  root_scaling_quadratic reclassified PREMISE -> DERIVED
regression: no pre-existing exported value, display, dim or step moved; WL re-run byte-identical;
     gradient FORM ablation moves 86 exported values
```

## What to attack

**① U-3 is the item with physics behind it — check the values, not just their presence.** `field_dimension`
is the step's `[u] = L` posit and `q_dimension` is the wavevector-norm posit. Derive both yourself and say
whether the exported vectors are right. Then check the claim that motivated exporting them: doubling
`field_dimension` moves `INERTIA_`, `STIFFNESS_` and `BARE_FIELD_COEFFICIENT_DIMENSION` in **both** engines
to the same wrong values with both exiting 0. Does exporting them change what a consumer can detect, or is
it only visibility? Say which — the round claims only visibility.

**② Did adding two entries disturb anything?** 139 → 141. Verify no pre-existing entry moved in any field,
that the two new ones carry no `dim` for the right reason (the solve produced none), that the round-trip
still covers all 141, and that the class tally now sums to the reconstructed count — the previous round left
a tally that could sum to 138 against 139.

**③ U-2's assert, from both sides.** Confirm it fires on out-of-vocabulary and empty strings, and confirm it
does **not** fire on anything legitimate. Then confirm the round is honest that omission and a
wrong-but-valid tag still pass silently — construct both and report what the export says.

**④ U-4 is a judgement call the user may reverse — give them the measurement.** `root_scaling_quadratic` was
`PREMISE`, now `DERIVED`. It is `lambda_scale**2 × a derived root`; its exponent is typed but not exported;
its identically-built sibling `root_scaling_scaled` is `DERIVED`; it moves under a form ablation. Is
`DERIVED` right, is `PREMISE` right, or does the class vocabulary genuinely lack a term for it? Argue from
measurement, not from the directive's reasoning.

**⑤ U-1 is ordering only.** Confirm every value is byte-identical and that nothing about reconstruction
changed.

**⑥ Regression.** No pre-existing exported value moves against `71f33901`; Wolfram still re-runs
byte-identical; a form ablation still moves exported values.

## Required method

Derive independently in CAS where a value is in question; a prose re-derivation is worth nothing. Write your
script before opening the artifact; save it and its literal stdout to named absolute paths and report them.
Ablate copies in a scratch directory — **never the working tree, and never `git checkout`/`restore`/`stash`:
the artifact is uncommitted.** For every claim ask *which line computed this?*

## Physics filter

> Report a finding only if it catches a way the physics could be wrong, or changes what may be claimed.

**This project's failure mode is adding machinery, not omitting it.** A finding that a piece of this is
unnecessary is worth as much as one that it is insufficient. Propose a new check only if its absence lets a
specific wrong answer through, and name it.

A leg returning "nothing survives" is weak evidence alone: say what you checked, what you could not, and
what would have had to be true to find something.

## Operational

One Mathematica kernel at a time (two seats); `timeout 600` per run. Do not modify the working tree. Do not
commit. Do not restore `reduction/`. Absolute paths. Save every ablation script and its stdout; report paths.

## Output

Rank most-severe first: **claim · evidence (literal command and output, or `file:line`) · what must change.**
