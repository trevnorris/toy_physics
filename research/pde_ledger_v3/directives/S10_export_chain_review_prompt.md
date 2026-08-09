# Independent physics review — S10 onto the export chain

## Artifact

Uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py
git diff -- research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out
```

plus the new untracked `research/pde_ledger_v3/scripts/S10_exports.py` (**generated**) and
`research/pde_ledger_v3/directives/S10_export_chain_directive.md`.

Written by Codex in one pass. You are one of two independent legs with this identical prompt; you will not
see the other's output.

## What this is

`scripts/S9_exports.py` is a generated flat `LEDGER` — one record per object S9 concluded,
`{display, value, dim?, class, step}`, keys carrying the component count of the construction that produced
them. **S10 is the first thing to import it.** It binds S9's objects instead of re-declaring them, adds its
own derivations across several component counts, overwrites what it re-derives, and writes the merged
`S10_exports.py` for S11.

S10 could not run at all before this change — it imported a module from a deleted directory.

## What to check

By mutation, not by reading. **Ranked: the first three are what this build exists to establish.**

1. **Does the import actually bind S9's objects?** A consumer that re-declares a symbol gets a *different*
   object — same printed name, different assumptions, and a difference that does not simplify to zero but
   **prints as though it does**. Check every symbol S10 uses that S9 also exports. Any S10 declaration
   shadowing an imported object is the defect this question exists to catch.
2. **Does the overwrite work, and can it fail?** Where S10 writes a key that already exists it is supposed
   to emit the upstream value, its own value, and the residual, then guard. Make S10 write a *wrong* value
   to an existing key and show what happens. Then ask the harder question: **is the upstream operand
   genuinely independent of the one S10 computed**, or does it descend from the same source? A residual
   whose operands share a source is zero by construction — three such guards have already been built and
   deleted on S9.
3. **Are the exported keys fixed, or minted from the answer?** S10's solver returns a number of roots and
   discovers a number of strata; before this change its tags were labelled by those indices. A ledger key
   whose *name* depends on the result cannot survive a form control. Change the action's form so the
   cardinality changes, and report what happens to the keys.
4. **Is the D-key partition right, object by object?** A key carries the component count of the
   construction that produced the object. ⛔ Free-symbol inspection gives the wrong answer in both
   directions and is not the test.
5. **Did the derivation move?** The engine now imports upstream objects with stronger assumptions. The
   directive claims some values moved as a result and names them. **Verify that list is complete** — a
   value that moved and is not on it is the finding.
6. **Is anything typed that the computation could produce?** Counts, dimensions, ranks, multiplicities, the
   variable part of any key. ⚠ The authored *name* of an object is exempt; a CAS expression does not carry
   its own name. ⚠ Supplied **premises** are exempt; a premise is an input.
7. **Did this build add an in-run check on the export writer that cannot fail?** Beyond the overwrite
   residual in item 2, it was told not to. Verify.

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY** — a change to the *structure* of a load-bearing object, not a rescale.
Run a control at the unchanged structure first and show your harness reproduces the artifact.

⛔⛔ **A test that passes on a weaker fix is not a test.** For each guard, construct the weakest change it
should reject and show whether it does.

Report with line numbers: anything typed where it could be read; a check whose expected value lives inside
the artifact it checks; an `assert` preceding the value it guards; an emitted object redundant with
another; a key or count assembled beside the computation rather than read out of it.

Ask of every claim in the directive: **which line computed this?** Give the line number or report it
uncomputed.

### Write your own script first

⛔ Write your own derivation/verification script **before** opening the artifact; save it and its literal
stdout to named absolute paths and report them. **Without them your claims will be discarded** — a prose
re-derivation is the defect class this rebuild exists to remove, relocated into the review.

## Do not read

- Anything under `directives/` beginning `S10_export_chain_decisions`, `S10_decisions_review_prompt`, or
  `S9_export_` / `S9_structural_` — those are the orchestrator's decision lists and earlier review briefs.
- Anything under `/tmp/s10_export_chain_evidence/` or `/tmp/S10_*` — that is the builder's own evidence.
  ⭐ Produce your own.
- Any Codex or review transcript; none is in the repository.

## Operational constraints

- Ablate **copies** under `/tmp`. ⛔ Never modify the working tree — it holds the uncommitted change and
  there is no committed baseline for the new file at all.
- ⚠ **This engine is slow** — it sweeps several component counts with a per-operation timeout. Wrap every
  run in `timeout 1800`. A hit is a failed ablation: report it and move on. ⛔ Never raise it further.
- ⛔ Do not start `wolframscript`; no Wolfram engine is in scope.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

Carried limits, out of scope unless this change altered one: S9's assumption set hand-types its predicates
so the assumptions it exports are invisible in its own record; S9's export boundary is a hand-typed class
comment; S9's unit-basis convention is unpoliced inside S9; `field_dimension` is an alias of
`length_dimension`; cross-engine (py↔wl) naming and the comparator are a separate pass.
