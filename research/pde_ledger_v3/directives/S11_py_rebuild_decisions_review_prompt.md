# Independent review — the S11 PY engine rebuild decision list, before the builder launches

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_py_rebuild_decisions.md`

It is **orchestrator-written**. It is the list a builder will rewrite a SymPy physics engine from. Nothing
has been built yet. You are one of two independent legs; the other is not visible to you.

⚠ The engine's **physics** comes from `directives/S11_SHARED_PHYSICS.md`, which has already passed its own
two review legs and is **not** under review here. ⭐ This list governs everything else: what the engine
exports, what it imports, how keys are named, and what guards it carries.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the decision list first

For a document, blindness comes from **reading order**.

1. `research/pde_ledger_v3/S9_REWRITE_PLAN.md` — the export-chain pattern and its settled decisions
   `D1`–`D11`. **This is the contract the list inherits.**
2. `research/pde_ledger_v3/scripts/S10_exports.py` — the previous step's generated output, which this
   engine imports. ⭐ It is large; **import it and inspect it programmatically** rather than reading it.
3. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — skim `§7` (packages), `§8` (tag grammar),
   `§Q6r`, and `§5`. ⭐ You need what the engine will EMIT to judge what the list says to EXPORT.
4. **Write down, before step 5:** what *you* would have to settle before handing a builder this rewrite.
   ⭐ Keep that list.
5. **Only now** read the decision list.

## What to check

### 1. ⭐⭐ IS THE EXPORT BOUNDARY RIGHT, AND IS IT STATED PRECISELY ENOUGH TO APPLY?

The list says: export every object `MAIN` emits at every `D`; control packages are evidence, not exports.
⭐ **Verify the claim it rests on** — that `S10_exports.py` carries no control-package key. Run it.
⚠ Then the harder question: is that boundary **correct for S11**? S11 has eight packages and a `_LOCAL_`
convention S9 lacked. ⛔ Does excluding all seven controls lose something a later step or the comparator
needs? ⭐ Say so concretely if it does.

### 2. ⭐⭐ IS ITEM 2's `_LOCAL_` RULE COHERENT?

It splits `Q6r`'s outputs: derived vectors export, imported/comparison objects do not.
⭐ Check that split against what `§Q6r` actually emits in the repaired spec. ⚠ Is the line drawn in a place
a builder can apply without judgement, and does it hold for **every** `_LOCAL_` object the engine will
produce, not just `Q6r`'s? ⛔ Name any object it fails to classify.

### 3. ⭐⭐ IS THE KEY-NAMING TRANSFORM ACTUALLY MECHANICAL?

The list claims `<QUANTITY>` lowercased plus `_d<n>` is a mechanical transform matching the previous step.
⭐ **Test it**: take real tag names from `§6`/`§8` of the spec and real keys from `S10_exports.py`, apply
the rule, and report collisions, ambiguities, or mismatches with the existing convention.
⚠ Specifically: do `ROOT<r>` and `STRATUM<s>` scopes survive it? Can two distinct emitted objects map to
one key? ⛔ A collision silently overwrites a chain entry.

### 4. ⭐ IS THE IMPORT/OVERWRITE RULE SAFE?

⭐ Does "overwritten in place, nothing deleted, untouched entries identical" hold as an invariant a run can
actually check? ⚠ What happens when this step derives an object the previous step exported **under a
different key** — is that an overwrite, a new entry, or a silent divergence? ⛔ That case is the one the
naming rule makes possible.

### 5. ⭐ ARE THE GUARDS REAL, OR DO ANY READ ZERO BY CONSTRUCTION?

⚠ A residual whose two operands come from the same route is structurally zero and tests nothing.
⭐ For each of the three guards, ask: **what would have to be broken for it to move?** ⛔ Name any guard
that cannot fail, and say what would make it able to.

### 6. ⛔ DOES THE LIST LEAK AN ANSWER, OR COMMISSION DAMAGE?

⛔ Does it state, imply, or supply a reason from which a builder could derive what any computation
produces? ⚠ **The test is DERIVABILITY, not literal presence** — a justification can fix a value without
printing it.
⚠ Is any item broader than the problem it names, or does it override something the spec already settles?
⭐ The list says the spec wins on conflict — check that no item quietly contradicts the spec anyway.

### 7. ⭐ WHAT IS MISSING?

⭐ Compare your step-4 list. What must be settled before a builder starts that this list leaves open?
⚠ A rewrite of this size fails on the thing nobody decided.

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, or a way this list would cause a
**wrong engine, or an engine that cannot be compared to its Wolfram sibling**, to be built. ⛔ Not style.

## Method

- ⭐⭐ **Quote both sides for every finding**: the list's text, and the source it fails against.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
  ⛔ A path that resolves is not a source.
- ⭐ Where a question is settled by computation — the export boundary, the naming transform, key collisions
  — **write a script, run it, and paste its literal stdout.** ⛔ Prose is not accepted for a contested claim.
- ⛔ Do **not** edit anything, and ⛔ do not write the engine. Read-only.
- ⭐ End with which of your step-4 items the list handles and which it misses, and state explicitly **what
  you checked that could have failed and did not.**
