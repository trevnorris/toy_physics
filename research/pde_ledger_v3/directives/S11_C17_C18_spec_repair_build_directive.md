# Build directive — repair `C17` and `C18` in `S11_SHARED_PHYSICS.md`

⭐ **Target file:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`
(1005 lines). ⛔ **Edit ONLY that file.**

⭐ **What must become true** is the seven decisions in
`directives/S11_C17_C18_spec_repair_decisions_v2.md`. ⭐ Read it first, in full. ⛔ Do not re-litigate it.
⚠ Its decision list had two independent review legs; ⛔ this build does not reopen them.

⭐ **Scope — `T1`–`T6` only.** ⛔ **`T7` is a COMPARATOR-CONTRACT obligation and is OUT OF SCOPE here**, with
one exception: if `T7` requires the **spec** to constrain how a boolean-valued object is *emitted*, make
that constraint and ⛔ nothing more.

## ⛔⛔⛔ THE ONE RULE THAT OVERRIDES EVERY OTHER INSTRUCTION

⛔⛔ **THE SPEC SAYS WHAT TO COMPUTE. ⛔ IT NEVER SAYS WHAT ANYTHING EQUALS, IS EXPECTED TO BE, OR WAS
MEASURED.**

⚠ The decision list and `DEFECT_REGISTER.md#C17`/`#C18` contain **measured results of this step's own
system** — root counts, a named coincidence locus, which engine emitted which stratum. ⛔⛔ **NONE of that
may appear in the repaired spec, in any form from which an engine builder could derive it:** ⛔ no count,
⛔ no value, ⛔ no dimension, ⛔ no sign, ⛔ no named locus, ⛔ no example that happens to be the answer.
⇒ ⭐ **You may say "emit the sub-locus where the count changes." ⛔ You may NOT say which locus that is, or
that there is one.**

⚠ A builder who can *derive* a value from your prose has been given the value. ⭐ **The test is
derivability, ⛔ not literal presence.**

## Where the defects live — ⭐ loci verified against the file today

| | section | lines |
|---|---|---|
| `C18` | `§5`'s locus protocol — the five suffixes | **`:230-251`** |
| `C17` | `§Q8b` — `STRATUM<s>_POINT` and the `Q3`/`Q4` rerun | **`:608-641`** |
| both | `§Q3` (`:307`), `§8` tag grammar (`:917`), `§9` (`:965`) | ⚠ change **only** as `T1`–`T6` require |

⛔ Do not renumber sections. ⛔ Do not reflow untouched prose. ⭐ Keep the file's existing voice and marker
conventions.

## ⛔ Constraints

1. ⭐ **Extend the `§5` protocol; ⛔ do not delete any of its five objects.** They are load-bearing and were
   written against a measured CAS limitation stated at `:232-237`.
2. ⛔ **Do not specify an admissibility ALGORITHM**, and ⛔ do not require any real-domain capability
   symmetrically — `:232-237` already measures that the two CASes differ there. ⭐ Name the **object** and
   the **status**; let each engine reach it however it can.
3. ⛔ **Do not introduce recursive stratification** or any canonicalisation of components.
4. ⭐ Every new tag must fit `§8`'s existing grammar, and ⛔ must not collide with a name already in use.
5. ⚠ `T3`'s witness exchange makes one engine read a value the other emitted. ⭐ **State it as an obligation
   on the ENGINE'S OWN computation** — each engine evaluates the received point against **its own**
   `_EQUATIONS` and **its own** `M`. ⛔ It must not become an instruction to import, agree with, or adjust
   toward the other engine.
6. ⛔ **A fold leaves stale prose above the fix.** ⭐ After editing a section, re-read that section **and
   its introduction** for a sentence that the edit made false. ⚠ This has bitten this file before.

## Deliverables

1. The edited file.
2. ⭐ A summary listing, per decision `T1`–`T6`: **which lines you changed, and the sentence that now makes
   it true.** ⛔ If a decision required no change, say so and why.
3. ⭐ Any place where two decisions conflict, or a decision cannot be met without stating a value ⇒
   ⛔ **STOP and report it. ⛔ Do not invent a compromise.**

⛔ Do not commit. ⛔ Do not edit any other file. ⛔ Do not run either engine.
