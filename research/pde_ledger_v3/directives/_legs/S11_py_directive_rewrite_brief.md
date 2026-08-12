# Rewrite brief — the S11 SymPy build directive

## What to produce

Rewrite `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_build_directive.md`
(new file; the existing `S11_sympy_no_ledger_build_directive.md` is superseded — ⛔ delete it).

It is the build directive for the S11 SymPy engine. ⛔ It is **not** the physics: the physics is
`research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`, which binds and wins every conflict.
`CLAUDE.md` binds.

## ⭐⭐ THE ONE RULE THAT SHAPES THIS DOCUMENT

⛔⛔ **POINT AT AN OBLIGATION THE SHARED SPEC ALREADY CARRIES. ⛔ NEVER RE-WORD IT.**
⚠⚠ Measured on the previous version: **every** restatement came out **weaker** than the original, and all
three of its blocking findings were the same defect. ⇒ ⭐ **Before writing any rule, search the spec for the
obligation it duplicates.** If one exists, cite it and move on.

The three that must not recur:

| the previous directive wrote | the spec already said it |
|---|---|
| a per-cell **completion registry** the engine maintains | `§7` — *the run record must be **observed**, not declared*; and `§7`'s no-dropping-a-cell-for-cost rule |
| *"absent interface content does not turn `§Q11`'s closure inventory into an operational failure"* | **corollary 4** — ⛔ emission is never conditional on a payload's value. ⚠ Naming `§Q11` as an exception re-opened the identical defect on `§Q10`'s pinned failure object ⇒ **state the property, ⛔ never a list of exceptions** |
| *"nonzero exit is reserved for an operational failure"* | `§7` — ⚠ *"operational failure"* proved reachable by **route choice**, not impossibility |

## ⭐⭐ WHAT IS NEW SINCE THE PREVIOUS VERSION — ⛔ its scope was reversed

⛔⛔ **THE ENGINE WRITES `S11_exports.py`.** ⚠ The previous directive said it writes no ledger; ⭐ the user
reversed that — **every step writes its output so the next step can import it**, and that accumulated
record is the point of the chain: the list of everything the model defines and the true list of knobs.

⭐ **The write rule is `F9`** in `directives/S11_export_chain_decisions_v2.md`. ⛔ Do not restate `F9`;
point at it. ⭐ Its shape, for orientation only: three outcomes per key against the imported `LEDGER` —
absent ⇒ bare, proved-equal ⇒ bare, not-proved ⇒ `s11_`-prefixed with the imported row untouched.

⚠ **`F9` is the write rule, ⛔ not an acceptance harness.** ⛔ Do not add guards, registries or invariants
to the engine to police it ⇒ ⭐ the same file's **build-review checklist** is what polices it, and it is run
by the **review legs against the built script**, ⛔ not by the script against itself.
⚠ ⭐ **A check cannot audit its own input** — ⛔ do not build one.

## ⚠ Two measured traps in the code this engine inherits from

⛔ Neither is a rule to restate — ⭐ they are facts about `S10_brane_mode_spectrum_sympy_audit.py` that a
writer reusing its export path will hit:

1. ⛔ `exact_reconstruction_match` (`:1906-1910`) decides equality by `type` **and `srepr`** — a pure
   **spelling** test. ⚠ Algebraically equal objects in different normal forms read as different.
2. ⛔ The provenance check (`:2087`) hard-codes one predecessor step by name.

## Boundaries

- ⛔ **Do not compute, quote or imply any S11 result** — no spectrum, determinant, root, count or sign. The
  engine has not been run and those answers must not exist in any readable artifact.
- ⛔ Do not state what any object equals, is expected to be, or was measured (rule 5).
- ⛔ Do not restate `§8`'s tag grammar, `§5`'s corollaries, or `§7`'s run-record and cost rules. Cite them.
- ⭐ **Shorter is better.** ⚠ The previous version was 81 lines and three of them were the blockers.
  ⛔ Length in this workstream has tracked defect count, not coverage.

## Deliverable

The rewritten directive, and a short note listing every spec obligation you **pointed at** rather than
re-worded, with its section number.
