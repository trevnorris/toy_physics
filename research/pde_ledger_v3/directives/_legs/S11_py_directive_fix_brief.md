# Fix brief — the S11 SymPy build directive. ⛔ The stdout half is CLEAR; ⛔ every item below is the EXPORT.

## Authority

Edit `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_build_directive.md` only.
⛔ Change no other file. `CLAUDE.md` binds; `S11_SHARED_PHYSICS.md` is the physics authority.

⭐ **Two legs reviewed it.** ⭐ Item 2 (pointer fidelity), item 4 (leakage) and item 5 (43 lines cut, **35
obligations classified, nothing dropped except what is below**) all came back clean, and the **stdout tag
stream is buildable from the `§§1-10` pointer alone.** ⛔ Do not rewrite what is working; ⛔ do not lengthen
the file for its own sake.

⛔⛔ **Everything below sits in the EXPORT.** ⭐ The directive names the export's **rule** (`F9`) and names
neither its **population**, its **key spelling**, nor what **publishes** it.

## ⭐⭐ 1 · THE POPULATION — decided, ⛔ not open

⭐ **The export carries every tag of the run's primary package, as S9 and S10 both do.** ⛔ Not a knob-only
register, ⛔ not every package.
⚠ **Why it is not a design choice:** it is the existing chain's shape, and ⛔ the narrow reading makes `F9`
a **no-op** — the keys that collide are the primary package's, and if they are never candidates the rule
the user decided has nothing to act on.
⛔ Say which package that is by pointing at `§7`, ⛔ not by naming one here.

## ⭐⭐ 2 · THE TAG→KEY TRANSFORM — ⛔ cannot be specified, ⭐ MUST be visible

⛔ **Do not supply a transform.** ⚠ `S11_SHARED_PHYSICS.md:1084` lets a builder name any object the file
does not name ⇒ ⛔ no rule can be total, and a mechanical transform is **design 1 of four blocked**
(`DEFECT_REGISTER.md#c20`).

⛔⛔ **But an engine in which `F9` NEVER FIRES must not be indistinguishable from one in which it fired and
found nothing.** ⭐ A builder keying rows differently from the imported ledger's own spelling produces
**zero** collisions, takes the absent branch on every row, and the run looks identical.

⇒ ⭐⭐ **What must become true: the engine EMITS, as computed objects, the key it derived for each row and
which of those the imported `LEDGER` already held.** ⭐ Operands, ⛔ not a claim and ⛔ not a count — a
review leg reads them and sees the collision surface for itself (rule 2). ⛔ Do not add a guard, a
threshold or an assertion over them.

## ⭐⭐ 3 · WHAT PUBLISHES THE EXPORT — `F6` requires the choice BEFORE the writer exists

⚠ `F6` offers two: publish only if every declared primary-package cell completed, **or** carry
machine-readable per-cell completeness every consumer checks. ⛔ The directive currently bans the mechanism
for the second and never picks the first.
⭐ **Choose the first**, and ⛔ do not build a completeness field: ⚠ `§7`'s observed run record already
carries which cells completed, and ⛔ a second one authored beside it is the registry that was already cut.

## ⭐ 4 · Three smaller ones

- ⭐ **`class` on a written row.** ⚠ The predecessor step chose it from an authored per-quantity map whose
  keys are all dead here. ⛔ Do not port that map. ⭐ Name what the class must be a property **of**; `§8`'s
  declaration classes and `§9` are the material.
- ⛔ **Restore the boundary the rewrite dropped: the engine modifies no file but its two products.**
  ⚠ `S10_exports.py` is a committed artifact, is `§Q6r`'s own operand, and `F4` is recorded **still open**
  — ⭐ a builder has a live reason to reach for it.
- ⭐ **`F9`'s equality and residual obligations moved** into `F9`'s rule section (they were misfiled under
  its review checklist, where a builder never reads them). ⛔ Do not restate them; ⭐ the pointer to `F9`
  now carries them.

## Boundaries

- ⛔ **No S11 result** — no spectrum, determinant, root, count, sign or dimension. None exists yet.
- ⛔ No value, expectation or measurement (rule 5). ⚠ Items 1-3 are all statable without one.
- ⛔ Do not restate `§8`, `§5`'s corollaries, `§7`, or `F9`. Cite them.
- ⭐ **Shorter is better.** ⚠ Length in this workstream has tracked defect count, ⛔ not coverage.

## Deliverable

The edited directive, and a note listing what you added for items 1-4 and which section you pointed at for
each.
