# Independent review — S9 export chain, fix round 4

## Artifact

One uncommitted file in `/var/projects/toy_physics`, on top of commit `ef0943ec`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
```

`git diff` is exactly this round. The extractor, the Wolfram engine, `S9_exports.py` and both `.out` files
are unchanged — the round reports the export and stdout regenerate **byte-identical**, which is itself a
claim to check. Pre-round state: `git show ef0943ec:<path>`.

Directive: `directives/S9_export_chain_fix_round4.md`. Plan: `S9_REWRITE_PLAN.md`. **Read directives last.**
⛔ **F-2 is with the user, out of scope, and must not be re-reported.**

## What this round claims

```
posited_output_classes = {"ROOT_SCALING_QUADRATIC": "PREMISE"}   at module scope
class_tag = posited_output_classes.get(suffix) or output_class_for(value)
assert set(posited_output_classes) <= set(main_outputs)

renamed key   -> exit 1 at that assert
unrenamed     -> stdout and exports byte-identical to ef0943ec
class provenance over 139 exports: 14 declaration-identity, 1 posited map, 124 default
```

And, from a mutation run on **each** of the twelve cross-engine rows:

```
10 counted; 2 excluded under rule 2 (FULL_ROOT_MULTISET, DYNAMICAL_MATRIX_ROUTE_RESIDUAL);
ZERO rule-1 (entailment) exclusions
```

⚠ That last line reverses three earlier positions of mine, including that
`COEFFICIENT_DIMENSION_DIFFERENCE` is "literally stiffness minus inertia so it cannot disagree" — the round
says the **subtraction is separately coded in each engine** and a sign reversal moves only that row.

## What to attack

**① Is the class fix real, or a rename of the same lie?** The previous mechanism was a hardcoded string
compared against a tag name. This is a module-scope dict keyed on the same tag name, plus an assert.
**Is that materially different, or the same fragility with a guard bolted on?** Test it: rename the dict
key; rename the `main_outputs` key; rename both consistently; delete the assert and rename. **Does any
route silently reclassify?** State plainly whether the guard fails on every weaker state or only the one
demonstrated.

**② Is `ZERO rule-1 exclusions` honest?** This is the round's most surprising claim and it inverts earlier
reasoning. For each of the ten counted rows, verify the mutation really moves that row **and leaves the rows
it was previously said to follow from unmoved** — in particular `COEFFICIENT_DIMENSION_DIFFERENCE`,
`IMPLIED_SPEED_DIMENSION`, `TRANSVERSE_SPEED_SQUARED` and `FULL_ROOT_MULTISET`. **Then attack from the other
side:** find a row where the two engines *cannot* differ because both compute it the same way from the same
emitted operands. A wrongly-counted row inflates the agreement claim exactly as a wrongly-excluded one
deflated it.

**③ The two rule-2 exclusions.** Both are said to print their value and then exit non-zero. Confirm the
value is emitted before the guard in both cases, and that the exit really is non-zero. Then ask whether
"a disagreement becomes a build failure" is the right reason to stop counting a row — or whether such a row
should be counted with its guard noted.

**④ Regression.** The round changed a source file and claims the generated export and stdout are
byte-identical to `ef0943ec`. Verify. Confirm the Wolfram engine still re-runs byte-identical, and that a
form ablation still moves exported values — the round reports 86 on the gradient form.

**⑤ Anything the round left.** The extractor was in scope and untouched, deliberately. The `field_dimension`
/ `length_dimension` alias is reported latent. The `.premises` file is reported orphaned. Judge whether
leaving each is right.

## Required method

Derive independently in CAS where a value is in question; a prose re-derivation is worth nothing here. Write
your script before opening the artifact and save it and its literal stdout to named absolute paths. Ablate
copies in a scratch directory — **never the working tree, and never `git checkout`/`restore`/`stash`: the
artifact is uncommitted.** For every claim ask *which line computed this?*

## Physics filter

> Report a finding only if it catches a way the physics could be wrong, or changes what may be claimed.

**This project's failure mode is adding machinery, not omitting it.** A finding that a piece of this is
unnecessary is worth as much as one that it is insufficient. Propose a new check only if its absence lets a
specific wrong answer through — and name it.

A leg returning "nothing survives" is weak evidence alone: say what you checked, what you could not, and
what would have had to be true to find something.

## Operational

One Mathematica kernel at a time (two seats); `timeout 600` per run. Do not modify the working tree. Do not
commit. Do not restore `reduction/`. Absolute paths. Save every ablation script and its stdout; report paths.

## Output

Rank most-severe first: **claim · evidence (literal command and output, or `file:line`) · what must change.**
