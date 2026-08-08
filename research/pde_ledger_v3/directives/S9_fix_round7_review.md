# Independent review — S9 export chain, fix round 7

## Artifact

Uncommitted in `/var/projects/toy_physics`, on top of `9c4ad8fc`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
research/pde_ledger_v3/scripts/S9_exports.py                                (regenerated)
research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out  (regenerated)
```

`git diff` is exactly this round. Pre-round: `git show 9c4ad8fc:<path>`. The Wolfram engine is untouched.
Directive: `directives/S9_export_chain_fix_round7.md`. Plan: `S9_REWRITE_PLAN.md`. **Read directives last.**
⛔ **F-2 is with the user, out of scope, not to be re-reported.**

## Context

Two of round 7's three items repair defects **round 6 introduced**. The third is a residual D10 already
required. Nothing here is new capability.

## What this round claims

```
V-1  LEDGER key q_dimension -> wavevector_norm_dimension, the name the .wl already emits
     classified as a premise of the SymPy engine ONLY; the .wl derives that vector
     whether py should derive it too is named as an OWED SPEC QUESTION, not done
V-2  the class assert moved from posited_output_classes onto the reconstructed ledger
     all four mutations exit 1 (posited OOV, empty, hardcoded-input OOV, default OOV)
     relocation-only run byte-identical to 9c4ad8fc; engine line count unchanged at 850
     caveat recorded: it still fires AFTER write_text, so a failed run leaves a bad export on disk
V-3  EXPORT_ROUNDTRIP_COUNT_RESIDUAL and EXPORT_CLASS_TALLY_RESIDUAL emitted, both 0
     round-trip residual now covers class and step as well as value and dim
     duplicate-name case: exit 1, counts 142/141, count residual 1
regression: 141 keys before and after, only the rename; after mapping it, no value or record moved;
     WL byte-identical; gradient FORM ablation moves 86 exported values
```

## What to attack

**① Is the rename right, and is it complete?** Verify the `.wl` really emits that object under
`WAVEVECTOR_NORM_DIMENSION` and really computes rather than types it. Then check the rename did not leave a
stale `q_dimension` reference anywhere, and that the exported **value** is unchanged. Ask the harder
question: is `wavevector_norm_dimension` the right standard name, or is there a better one both engines
could adopt? ⚠ And check whether any **other** LEDGER key is named after an engine-local placeholder — this
round fixed one instance; the class may be larger.

**② Is the relocated assert genuinely a control outside the thing it polices?** Re-run all four mutations
yourself. Then try to defeat it: a class tag that is in `CLASS_TAGS` but wrong for the object; a class
injected after reconstruction; a fifth source of class tags. State plainly what it does and does not cover
— round 6's assert was described as covering more than it did, and that is the failure being repaired.

**③ Do the two new residuals actually catch what they claim?** Construct the duplicate-name case and any
other way to make live and reconstructed counts diverge, or the tally diverge from the count. Check the
residuals are computed from live operands rather than restated. ⚠ And check the **new asserts sit after
their emits**, per this project's print-then-guard rule.

**④ The classification claim.** `wavevector_norm_dimension` is now described as a premise of the SymPy
engine only. Is that right — and does the `LEDGER` say so, or does a consumer just see `PREMISE` with no
indication that the other engine derives it? If the export cannot express "typed here, derived there," say
so; that is a finding about the schema, not a request to build one.

**⑤ Regression.** 141 keys, the rename only, no value moved after mapping it, Wolfram byte-identical, form
ablation still 86.

## Required method

Derive independently in CAS where a value is in question; a prose re-derivation is worth nothing. Write your
script before opening the artifact; save it and its literal stdout to named absolute paths and report them.
Ablate copies in a scratch directory — **never the working tree, and never `git checkout`/`restore`/`stash`:
the artifact is uncommitted.** For every claim ask *which line computed this?*

## Physics filter

> Report a finding only if it catches a way the physics could be wrong, or changes what may be claimed.

**This project's failure mode is adding machinery, not omitting it** — and rounds 6 and 7 are evidence:
round 6's additions created the defects round 7 repairs. A finding that something here should be **removed**
is worth more than one that something should be added.

A leg returning "nothing survives" is weak evidence alone: say what you checked, what you could not, and
what would have had to be true to find something.

## Operational

One Mathematica kernel at a time (two seats); `timeout 600` per run. Do not modify the working tree. Do not
commit. Do not restore `reduction/`. Absolute paths. Save every ablation script and its stdout; report paths.

## Output

Rank most-severe first: **claim · evidence (literal command and output, or `file:line`) · what must change.**
