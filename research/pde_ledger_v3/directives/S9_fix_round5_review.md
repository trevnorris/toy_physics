# Independent review — S9 export chain, fix round 5

## Artifact

Uncommitted in `/var/projects/toy_physics`, on top of `4ac4f1fc`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py   (modified; 851 -> 846 lines)
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.premises   (DELETED)
```

`git diff` is exactly this round. `S9_exports.py`, both `.out` files and the Wolfram engine are unchanged —
the round claims the export and stdout regenerate **byte-identical**, which is itself a claim to check.
Pre-round state: `git show 4ac4f1fc:<path>`. Directive: `directives/S9_export_chain_fix_round5.md`.
Plan: `S9_REWRITE_PLAN.md`. **Read directives last.**
⛔ **F-2 is with the user, out of scope, and must not be re-reported.**

## What this round claims

It **deleted** the class identity-matching machinery (`annotated_objects`, `output_class_for`,
`output_without_declaration_site_class`) and folded all six `PREMISE` names into the guarded
`posited_output_classes` map.

```
stdout md5 bef65b4b…  export md5 71407fad…   both match 4ac4f1fc
846 lines
GRADIENT witness:  exit 0, lagrangian = PREMISE          (previously silently DERIVED)
COLLISION witness: exit 0, transverse_multiplicity = DERIVED  (previously silently CONTROL)
11 counted, 1 excluded (DYNAMICAL_MATRIX_ROUTE_RESIDUAL, exact zero-matrix guard)
.premises deleted after grep showed no reader
```

## What to attack

**① Did deleting the identity route close the holes, or move them?** Re-run both witnesses yourself. Then
look for new ones: a `PREMISE` name that no longer resolves; an object whose class was previously right *by
identity* and is now wrong or defaulted; an entry in `posited_output_classes` naming something that is not
emitted; a falsy class string. The round also swapped `... or ...` for `.get(suffix, "DERIVED")` — check
that changed nothing it should not have.

**② Is the guard sound now, or still partial?** A leg previously established that the assert fails closed on
key-divergence but **not** on omission of an entry, nor on a wrong class string. With six entries rather
than one, that exposure is six times larger. State precisely which weaker states fail closed and which pass
silently — and whether the emitted class tally moves in the ones that pass.

**③ The 11 / 1 table, attacked from both sides.** For each counted row: does a mutation move it while
leaving the rows it is said to follow from unmoved, **and** can its value be recomputed from other counted
rows plus premises both engines type? ⚠ Note `INERTIA_COEFFICIENT_DIMENSION`,
`STIFFNESS_COEFFICIENT_DIMENSION` and `COEFFICIENT_DIMENSION_DIFFERENCE` are marked mutually recomputable —
**any two determine the third.** Is counting all three an overstatement of the agreement, or does each carry
something the others do not? Give the answer per row, not a rule.

**④ The single exclusion.** `DYNAMICAL_MATRIX_ROUTE_RESIDUAL` is excluded on the claim that its guard is
exact equality against the zero matrix, so **any** non-zero value exits 1. Verify that this is an invariant
of the guard rather than an outcome of the mutations tried — that is exactly the distinction that broke the
previous round's `FULL_ROOT_MULTISET` exclusion.

**⑤ The deletion.** Confirm nothing read `.premises`. Then judge: does deleting it lose anything the export
does not now carry?

**⑥ Regression.** Verify the two md5s against `4ac4f1fc`, the line count, that the Wolfram engine still
re-runs byte-identical, and that the gradient form ablation still moves 86 exported values.

## Required method

Derive independently in CAS where a value is in question; a prose re-derivation is worth nothing. Write your
script before opening the artifact and save it and its literal stdout to named absolute paths. Ablate copies
in a scratch directory — **never the working tree, and never `git checkout`/`restore`/`stash`: the artifact
is uncommitted and one file is a deletion that a restore would resurrect.** For every claim ask *which line
computed this?*

## Physics filter

> Report a finding only if it catches a way the physics could be wrong, or changes what may be claimed.

**This project's failure mode is adding machinery, not omitting it** — this round removed some, and a
finding that it removed too little is as valuable as one that it removed too much. Propose a new check only
if its absence lets a specific wrong answer through, and name that answer.

A leg returning "nothing survives" is weak evidence alone: say what you checked, what you could not, and
what would have had to be true to find something.

## Operational

One Mathematica kernel at a time (two seats); `timeout 600` per run. Do not modify the working tree. Do not
commit. Do not restore `reduction/`. Absolute paths. Save every ablation script and its stdout; report paths.

## Output

Rank most-severe first: **claim · evidence (literal command and output, or `file:line`) · what must change.**
