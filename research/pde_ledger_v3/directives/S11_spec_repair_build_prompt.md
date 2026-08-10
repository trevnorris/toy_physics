# Apply the S11 spec repair

## What you are editing

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **914 lines**, at
HEAD. ⛔ This is the **only** file you edit.

## The decision list — apply it exactly

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_spec_repair_decisions_v2.md`

It has been through two independent review legs and is folded. ⭐ **Apply items 1–5 and nothing else.**

- ⛔ Do **not** implement anything under *"What this round does not do"* or *"Registered separately"*.
  Those are deliberate exclusions, not oversights. `C17` and `C18` in `DEFECT_REGISTER.md` are **open on
  purpose** and are ⛔ not repaired here.
- ⛔ Do **not** touch either engine, any `out/` file, or any other document.
- ⭐ Everything the list does not name stays **byte-identical**.

## What this specification is

Two engines — one SymPy, one Wolfram — are rebuilt from this file **independently**, and their outputs are
compared. ⇒ ⭐ **Any freedom you leave in the wording becomes a divergence between two engines that both
followed it.** That is the failure mode every item on the list exists to remove, so where the list pins a
name, a payload shape, or a term decomposition, **transcribe it exactly** rather than paraphrasing.

⛔⛔ **THIS FILE MAY NEVER STATE WHAT ANYTHING COMES OUT TO BE.** No value, no count, no sign, no rank, no
dimension, no expected effect of any package. It says **what to compute**, as equations and definitions.

⚠⚠ **`DEFECT_REGISTER.md`'s `C16`, `C17` and `C18` entries contain MEASURED VALUES** — root lists, mode
counts, stratum equations. You will read those entries because the decision list cites them.
⛔⛔ **NOT ONE OF THOSE NUMBERS MAY APPEAR IN THE REPAIRED SPECIFICATION**, and no sentence may be written
from which one could be inferred. ⭐ The register explains *why* a repair is needed; the spec says only
*what to compute*.

## Style — match what is already there

⭐ Read enough of the file first to match its existing voice, structure and markup conventions. Repairs
should be indistinguishable from the surrounding text in form. ⛔ Do not restructure sections the list does
not name, do not renumber, and do not "improve" adjacent prose.

⚠ Several items change a **cross-reference target** as well as a definition — for instance a rule that
reads *"from the emitted stiffness functional"* must end up consistent with a package now carrying two
members. ⭐ After each item, re-read the sections that cite what you changed and make them agree. ⛔ Leaving
a stale cross-reference is the same defect as not applying the item.

## Report back — ⛔ under 25 lines

1. The file's new line count, and one line per item saying **where** you applied it.
2. ⭐ Anything in the decision list you found **ambiguous, contradictory, or impossible to apply** without
   choosing something it did not specify — **name the choice you made.** ⭐ This is wanted and is more
   valuable than a clean application.
3. ⭐ Any place the repair made an **existing** passage wrong, inconsistent, or unreachable.
4. ⛔ Do **not** report what any computation will produce, and ⛔ do not evaluate whether the repair is a
   good idea.
