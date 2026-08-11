# S11 PY engine — decision list v3

**Orchestrator-written, 2026-08-11.** ⛔ v1 (`11bf8e05`) and v2 are **records of blocked gates**; ⛔ do not
build from either.

**Deliverable:** `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`, a **full rewrite**,
which on running writes `research/pde_ledger_v3/scripts/S11_exports.py`.

**The physics is `directives/S11_SHARED_PHYSICS.md` at `cf4a21a4`.** ⛔ This list adds no physics, names no
object the spec does not name, and states nothing about what any computation produces.
⭐ Where this list and the spec appear to conflict, **the spec wins and the conflict is reported**, ⛔ not
resolved by the builder.

**The chain rules are `F1`–`F8` in `directives/S11_export_chain_decisions_v2.md`.** ⛔ A builder does not
re-choose them. ⭐ In particular `F1` — **storage keys are flat, and a later step re-deriving one object
writes the same key, because the collision is the measurement.**

⛔ The existing file is replaced, ⛔ not repaired: `:21` imports from a `reduction/` directory that no longer
exists.

---

## The decisions

**P1 · What the export carries.**
⭐ Every object the `MAIN` package derives, at every `D` in its sweep; ⭐ the premises those results hold
under; ⭐ the symbol rows the exported values reference; ⭐ the sweep's own observed run record.
⛔ The seven control packages are ablation **evidence**, ⛔ not exports.
⛔ The previous-step comparison is **not** exported as a standalone object — ⭐ its operands and residual
ride on the rows they concern (`F3`), and ⛔ a solver-condition tag is engine internals and is not exported.

⚠ ⭐ **The boundary must be decidable for the symbol population too, ⛔ not only for tagged objects** — a
symbol row is manufactured by traversing exported values, so ⛔ a rule stated only over tags leaves it
undecided.

**P2 · The row schema is READ FROM THE IMPORTED ARTIFACT AT RUN TIME, ⛔ never typed into the engine.**
⭐ The writer discovers the fields the imported rows carry and writes rows of that shape. ⛔ Do not type a
field list, ⛔ do not add a field the imported artifact does not carry except where `F3` requires one, and
⛔ do not drop a field from an imported row.
⚠ ⭐ **A field list typed beside the writer is a value the engine did not read** — the defect class this
rebuild exists to remove.

**P3 · A row this step re-derives carries its evidence IN THE ROW** (`F3`).
⭐ Both operands and the residual, as CAS objects in the row, so a consumer reading **only** the merged
export can recompute the comparison. ⛔ A corroboration claim with no operands in the file that carries it
may not be written, and ⛔ may not be propagated from an imported row that lacks them.

**P4 · "Same object" is decided by TWO operands, ⛔ never by one predicate.**
⭐ Emit, for every key that exists in both the import and this step's output:
- the **structural identity** of the two reconstructions, assumptions included;
- the **residual**, wherever subtraction is defined for that object;
- and ⭐ **which of the two distinguished them**, when they disagree.

⛔ A single "residual reduces to zero" test is not admissible: it is undefined for most of the objects the
chain carries — tuples, relations, boolean and predicate payloads, strings — and it loses domain
information, calling `x/x` and `1` the same object. ⛔ A zero test that reads `residual == 0` is also wrong
for matrix payloads.
⭐ **The predicate must be total over every value type the export admits**, and where it cannot decide, it
⭐ reports `UNDECIDED` with both operands. ⛔ It may not throw, and ⛔ it may not silently skip a row.

⭐ Under `F2`: same object ⇒ a re-derivation; ⛔ different object ⇒ a finding that fails loudly, naming both.
⛔ Never a silent overwrite, and ⛔ the builder does not resolve a collision by renaming.

**P5 · Chain integrity compares WHOLE IMPORTED ROWS against a snapshot, ⛔ not values, and ⛔ not against
the engine's record of what it wrote.**
⭐ Snapshot the imported `LEDGER` at import as an **independent reconstruction**, ⛔ not a reference.
⭐ Fix the set of keys this step **declares** it re-derives **before the write loop opens**; ⛔ it may not be
computed from what was written.
⭐ Emit three operands: the imported rows that **differ in any field** between snapshot and output; the
**declared** set; and ⭐ **the first minus the second** — a row that changed without being declared.
⚠ ⛔ Not the symmetric difference: a declared re-derivation that reproduces its row exactly is correct, and
a guard that flags it is wrong.
⚠ ⭐ Any field may carry physics: a `dimension_key` rewritten on an imported row redirects a dimension
lookup with **no value anywhere changing**.

**P6 · Publication is atomic with respect to completeness.**
⭐ An export is published only if every declared `MAIN` cell completed, read from the sweep's **observed**
run record, ⛔ never from the declared sweep.
⛔ **And an incomplete or failed run must leave no previously published artifact importable as though it
were this run.** ⚠ Refusing to write is not enough: the previous file remains on disk, reads as complete,
and hands a later step physics from a run that did not happen.

**P7 · Every eligible emitted object corresponds to exactly one written row.**
⭐ Emit the eligible emitted set, the written key set, and their difference, then guard.
⚠ ⭐ Without this, a completed computation can print its tag and be **absent from the export** while the
round trip, chain integrity and completeness checks all pass — ⛔ the absent-computation class, which is
the most expensive defect this project has had.

**P8 · Reconstruction is verified by a round trip in the run** (`D3`): re-read what was written, and emit
both operands and the residual.

**P9 · Import `S10_exports.LEDGER`, bind what this step consumes, export the merged dict.**
⛔ Do not fence, filter or scope the import. ⚠ The downstream step seeing the upstream one is the point.
⭐ A `Q6r` lookup that cannot resolve is reported as **unresolved, generically** (`F7`) — ⛔ no placeholder,
⛔ no exception, ⛔ no membership of either object stated anywhere a builder can read it.

**P10 · Every declared symbol carries a class tag and an English description, read live.**
⛔ No `VERDICT`, no `PASS`, no `FAIL`, no summary judgement (spec §5). ⛔ A physics finding exits `0`;
non-zero is for operational failure only.

**P11 · ⭐⭐ REPORTING IS SUCCESS.** ⭐ If an item here is ambiguous, under-determined, ill-posed,
tautological, or cannot be built from what the spec supplies, ⭐ **report it and build the rest.** ⛔ Do not
invent a mechanism to cover the gap, and ⛔ do not resolve a conflict with the spec yourself.

---

## ⛔ What this list does not decide

- Anything the spec decides, and anything `F1`–`F8` decide.
- The Wolfram engine, the comparator contract, the step record, the paper card, any register.
- The S10 regeneration that `F8` requires; ⭐ this engine is built against whatever key namespace the
  imported artifact presents at run time, ⛔ and hardcodes no key from it.
