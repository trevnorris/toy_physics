# Decision list v2 — the export chain

**Orchestrator-written, 2026-08-10, folding the ten findings that blocked v1.**
⚠ v1 (`36589024`) is a **record**; ⛔ do not read its `E1`–`E7` as decisions.

⛔ **This list is not itself a build directive.** It settles the chain questions the **S11 PY decision
list** must answer; ⭐ that list gets the two legs, before any builder launches (rule 7).

⭐⭐ **What v1 got wrong, in one line:** it made cross-step collisions impossible by construction — ⛔ and
the collision **is** the measurement. Two steps deriving one object must be able to **meet**.

---

## The decisions

**F1 · Storage keys stay FLAT. `D5` is unchanged; ⛔ there is no producer prefix.**
⭐ A step writes the object's name. A later step re-deriving that object writes **the same key**.

**F2 · Before writing a key that exists in the imported `LEDGER`, the writer compares the OBJECT.**
⭐ This is `DEFECT_REGISTER.md:675`, which already prescribed it.
- ⭐ **Same object** ⇒ a **re-derivation**: emit both operands and the residual, then guard.
- ⛔ **Different object** ⇒ a **finding that fails loudly**, naming both. ⛔ Never a silent overwrite.
⚠ Deciding "same object" is the load-bearing part and belongs to the S11 PY list.

**F3 · A re-derived row carries its own evidence, in the row.**
⛔ `corroborated_steps` alone is an agreement claim with no operands in the file that carries it, and
`S11:527-529` **specifies `Q6r` to propagate it** ⇒ a consumer forwards a claim it cannot check.
⭐ A consumer reading **only** the merged export must be able to recompute the residual.

**F4 · S10's export is REGENERATED under `F3`. ⛔ It is NOT re-keyed and its tag names do not move.**
⭐ Both legs agree by measurement that S10 need not be re-keyed before S11 is built.
⚠ ⭐ Every value must be **byte-identical** apart from the fields `F3` adds — ⛔ if any value moves, it was
not a regeneration.

**F5 · `C19` is a REAL deviation, ⛔ not a traversal gap.** `S10_SHARED_PHYSICS.md:197` orders
*"Emit `M_A`, `M_B`, and `M_A − M_B`"*; both engines emit `Q2_MATRIX_A/B`.
⭐ S10's step record must **disclose** it. ⛔ The rename itself is its own gated workstream — ⚠ injectivity
across the worklist first — and ⭐ **S11's build does not depend on it**, because S11's spec fixes its own
quantity names.

**F6 · An export is published only if every declared `MAIN` cell completed** — ⭐ or it carries
machine-readable per-cell completeness that every consumer checks. ⛔ One or the other, chosen before the
writer is built (`S11:895`).

**F7 · A `Q6r` lookup that cannot resolve is reported as unresolved, generically.**
⛔ No placeholder, ⛔ no exception, ⛔ no membership stated anywhere a builder can read.

---

## ⛔ What this list does not decide

- ⭐ The **"same object" predicate** of `F2`, and the row shape of `F3` ⇒ the S11 PY decision list.
- The `C19` rename worklist ⇒ its own gate, per `F5`.
- `C17`, `C18`, S10's requirements registers.
