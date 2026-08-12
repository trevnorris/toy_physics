# Fix brief 4 — the S11 SymPy build directive. ⛔ ONE root cause, ⭐ three symptoms fall out of it.

## Authority

Edit `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_sympy_build_directive.md` only.
⛔ Change no other file. `CLAUDE.md` binds; `S11_SHARED_PHYSICS.md` is the physics authority.

⚠⚠ **The population sentence in this directive is FACTUALLY WRONG, and it is the orchestrator's error, ⛔ not
yours.** ⭐ It was decided from precedent reasoning and ⛔ the census was never run. ⇒ ⭐ **Three separate
findings across two review rounds are symptoms of it.** ⛔ Fix the root; ⛔ do not patch the symptoms
individually.

⭐ **Cleared and ⛔ not to be touched:** the `F4`/non-regeneration sentence; the `F1`+`F9` naming sentence;
the `F6` first-branch choice; the read-only boundary; the key/imported-key emission; the `F2`/`F9c`
supersession property; the `S10_…py` line citations; leakage.

## ⛔⛔ 1 · THE POPULATION — ⭐ measured, ⛔ and the current sentence captures half of it

> `Its S11 candidate population is every tag emitted for the package §7 identifies as primary.`

⭐ **Measured on the committed predecessor export** (evidence in
`directives/_measurements/S11_sympy_build_directive.md`):

```
total rows      : 617
tag-shaped keys : 0
KNOB rows       : ['rho_br', 'mu_R']
STRUCTURAL rows : ['D', 'length_dimension', 'time_dimension', 'mass_dimension']
```

⇒ ⭐⭐ **The predecessor's export has TWO sources**, and its own code says so at
`scripts/S10_brane_mode_spectrum_sympy_audit.py:2005` (`add_symbol_records`): ⭐ **tag-derived rows, keyed by
the object's name** — ⛔ not by the tag — ⭐ **plus one row for every free symbol occurring in those values**
that does not already have one.

⚠⚠ **Every `KNOB` and every `STRUCTURAL` row in the whole chain comes from the second source.** ⇒ ⛔ a
tag-only population cannot carry a knob, and ⭐ *"the true list of knobs"* is the stated purpose of the
export.

⇒ ⭐ **State the population to match what is measured.** ⛔ Do not name any symbol; ⛔ do not state which
symbols are or are not already present in the import; ⛔ do not specify how symbols are discovered — ⭐ name
the property: the export carries the primary package's computed objects **and** the free symbols they
contain.

## ⭐ 2 · THE THREE SYMPTOMS — ⛔ verify each is closed by item 1, ⛔ do not fix them separately

- ⭐ **The `class` rule's only operative clause** (*"a declared symbol carries its declared class"*) had **no
  referent** in a tag-only population. ⇒ ⭐ with symbols in the population it does. ⛔ Leave the rule's
  wording alone unless item 1 fails to give it a referent.
- ⭐ **A coefficient of the primary package that the import lacks** could never get a row. ⇒ ⭐ it is a free
  symbol of that package's action, so item 1 admits it.
- ⭐ **`§Q6r`'s two-step lookup** needs a coefficient row carrying a `dimension_key`. ⇒ ⭐ item 1 supplies the
  coefficient row. ⛔ Do **not** add a rule for how a `dimension_key` is chosen — ⚠ re-pointing one is
  `OWED` item 4, the failure where every repo check passes.

## ⛔ 3 · *"§Q6r's published-export lookup"* IS AMBIGUOUS — ⭐ one clause

⚠ `§Q6r`'s lookup is defined **into the imported ledger**, so the phrase attributes to `§Q6r` a term it does
not have, and splits two ways: ⛔ (a) the lookup `§Q6r` already performs — ⚠ then the object **duplicates
`§Q6r`'s own emitted objects**, measures nothing new, and splits one named object across several payloads,
which `§8` forbids; ⭐ (b) that lookup's **shape applied to the export this step publishes** — the only
reading that measures anything.
⇒ ⭐ **Say (b) explicitly**, in terms of the export this engine writes.

## ⛔ 4 · `F3` — ⭐ POINT AT IT, ⛔ do not transcribe it

> `Meet **F3**: the merged export alone must let its consumer recompute what each row claims; do not
> prescribe its row shape.`

⛔ **This transcribes `F3`'s weaker clause and drops its headline scope — *in the row*.** ⚠ Demonstrated
satisfiable with **zero change** on the committed predecessor: a row claiming agreement carries no operands,
while a row that happens to contain them sits a few lines above it in the same file, and the file-scoped
reading calls that met. ⇒ ⭐ `F3`'s stated defect survives intact.
⚠ **And as written it collides with the cleared `F4` sentence:** it quantifies over *every* row of the merged
export, including carried-forward rows this build must not touch and rows `F9c` orders left exactly as they
stand.
⇒ ⭐ **Point at `F3` the way the adjacent sentence points at `F1` and `F9`** — *apply it; do not restate it*
— ⭐ and bound it to the rows **this step writes**.

## Boundaries

- ⛔ **No S11 result** — no spectrum, determinant, root, count, sign or dimension.
- ⛔ No value or expectation (rule 5). ⛔ Do not name a coefficient or state any row's membership.
- ⛔ Do not restate `§4`, `§8`, `§5`'s corollaries, `§Q6r`, `F1`, `F3` or `F9`. Cite them.
- ⭐ Keep it short; ⛔ items 3 and 4 should each shrink the file.

## Deliverable

The edited directive, and a note saying what you changed for each item and — for item 2 — whether item 1
gave each symptom a referent, or whether one still needs its own sentence.
