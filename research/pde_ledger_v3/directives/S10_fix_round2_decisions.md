# S9/S10 export chain — fix round 2. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

In scope — **all of these, because items below change each**:
`scripts/S9_light_requires_shear_sympy_audit.py` · `scripts/S10_brane_mode_spectrum_sympy_audit.py` ·
`scripts/extract_knob_inventory.py` · both generated `scripts/S9_exports.py` and `scripts/S10_exports.py` ·
both regenerated `scripts/out/S*_audit.out` · a new directive under `directives/`.
⛔ Out of scope: every `.wl`, every step record, every `.tex`.

⚠ **S9 is reopened deliberately** — F5 is a defect at S9's construction site that S10 inherits and cannot
repair without breaking the untouched-entries-are-identical property. Regenerate S9 first, then S10.
⚠ S10 is slow: wrap runs in `timeout 1800`.

⭐ The physics is independently confirmed (a leg reproduced 26/26 tags from its own derivation). Every item
below is the export chain. ⛔ The derivation, action, ansatz and every computed physics value stay untouched.

---

## F1 · The ledger must carry the objects its own values are written in.

Measured on the committed artifact: **158** distinct symbol names appear across S10's exported values and
dimensions, and **3** are ledger keys (`D`, `mu_R`, `rho_br`). S9: 23 names, the same 3. So `omegaSquared`
(in 60 exported values), `kx`/`ky`/`kz`, `k1…k5`, `a1…a5`, `t`, `x`/`y`/`z` are objects a consumer can only
re-declare — and this step's own record already measures what re-declaring costs: substituting a look-alike
into the determinant **left it unchanged and printed nothing**.

**Decision: a consumer binds every symbol it needs; it never re-declares one.**

⚠ **The genuine difficulty, and it is yours to resolve:** some of those 158 names are not objects at all —
they are **authored field names** inside the `indexed_derivations` rows, built with `sp.Symbol` because a
struct key was wanted. ⭐ **The separation must be structural, ⛔ never a hand-maintained list of which
names are labels** — a list is the same hand-typing defect relocated. ⛔ If you conclude the two populations
cannot be separated structurally, **name that and stop.**

⭐ `declaration_classes(...)` is called at `S10:1872` and `del`'d at `:1873` — the classification is computed
and thrown away, which is why the ledger tally reads `COORDINATE 0` while `omegaSquared` is annotated
`# COORDINATE ·`. ⭐ S9 already solved the analogous problem for `KNOB`/`STRUCTURAL`; do what it does.

## F2 · An imported entry's provenance cannot change silently.

`S10:1888-1903` compares `value` only, then `:1908-1917` writes `'class': record.class_tag, 'step': 'S10'`
unconditionally. ⇒ a key collision on a same-valued upstream record **relabels its class and step with
residual 0 and every guard green.**

⚠ Measured, so state it accurately: **145 of 145 S9 classes survive into the committed S10 ledger and none
changed.** This is a **missing guard**, ⛔ not a present corruption.

**Decision: what the overwrite changes is emitted and guarded — not only the value.**

## F3 · One object, one record.

**48 of the 50** records carrying a `dim` export that identical object **again** as a standalone record;
the only two that do not are `rho_br` and `mu_R`. ⇒ a downstream step differencing a record's `dim` against
its sibling gets **zero by construction** — the recurring defect, pre-seeded into S11's inputs.

**Decision: the ledger does not carry the same object under two keys.** ⭐ You choose which representation
survives; **report which, and what a consumer loses.**

## F4 · Every symbol name in the written ledger maps to exactly one object.

⭐⭐ **This is the one in-run check sanctioned on the export writer** — it walks the *written module* and
asserts a property *of* it, ⛔ it is not another residual whose operands descend from one source.

Measured today: **158 names, 0 bound to more than one object** — so ⛔ a passing run proves nothing by
itself. **Demonstrate it able to fail**: introduce a same-named symbol with different assumptions and show
the run rejects it.

## F5 · An exported value is a well-formed SymPy object all the way down.

`S10`'s `q6_dimension_solution` and `S9`'s `dim_solution` nest a **raw Python `list`** inside the stored
object. Measured: `.free_symbols` and `.atoms()` both raise `AttributeError`; `srepr`/round-trip pass, which
is why every existing guard is green. ⇒ a consumer cannot traverse the record, **and F4's scan is blind
exactly here.**

**Decision: every stored value is traversable.** ⛔ Do not special-case the two records — the construction
that produced a raw list is the thing to fix, and it is in both engines.

---

## Constraints

- ⭐ Every S9 record S10 does not overwrite stays identical between the two generated modules; the three
  overwrites stay.
- ⭐ Report both complete generated-module diffs and both complete `.out` diffs against the current working
  tree, and how many lines this round added and removed.
- ⛔ Do not add any in-run check on the export writer beyond F4. **Four** have been built and deleted on this
  project; each compared two operands descending from one source, and one passed while the change it
  policed was reverted entirely.
- ⛔ Do not state anywhere what a count, tally or partition size comes out as. Emit it and let it be read.

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** For each item, construct the weakest change that
should be rejected and show whether it is.

⛔⛔ **A FORM ablation is mandatory** — change the structure of a load-bearing object, not a coefficient, and
report the literal diff. A coefficient rescale tests arithmetic; only a form change tests physics.

## Deliverables

The changed files, both literal stdouts, all four diffs, and every ablation script with its stdout at named
absolute paths **outside the repository**.
