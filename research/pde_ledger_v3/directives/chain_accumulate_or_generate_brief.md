# ⭐⭐ AUTHOR one decision — is the merged list ACCUMULATED at each write, or GENERATED from the steps?

⚠ **You are the author.** ⛔ Not a reviewer. A fresh Claude agent and Grok will review what you write.
⛔ **One decision.** ⛔ Do not redesign the export chain, ⛔ do not solve the naming problem, ⛔ do not write
a build directive. ⚠ Four attempts to design more than one thing at once were all blocked.

## The requirement — ⛔ NOT in question, ⛔ do not weaken it

⭐⭐ **The chain exists to end with (a) one authoritative list of everything defined across the whole model,
and (b) the true knob frontier** — what remains undetermined after every step has had its say. ⚠ *(user,
2026-08-11.)* ⛔ **Any answer that does not deliver both is not an answer.**
⭐ `S9_REWRITE_PLAN.md#D7` states the frontier as *"the union across steps minus anything a later step
derives — and the overwrite is what marks that."*

## The question

⭐ Today, `S9_REWRITE_PLAN.md#D5`: each step **imports the previous step's `LEDGER`, adds its own, overwrites
what it derives, and exports the MERGED flat dict.** ⇒ each artifact physically **copies** every upstream
row.

⭐⭐ **Should it stay that way, or should each step export only what it derives, with the merged list and the
knob frontier GENERATED from the per-step exports on demand?**

## ⭐ What the answer must weigh — ⛔ measure, ⛔ do not assume

- ⭐ **Measure the copying.** How many rows of each committed export are copies of an upstream step's? The
  files are `scripts/S9_exports.py` and `scripts/S10_exports.py`. ⭐ Project the trend forward over the
  steps still to be built (`REBUILD_HANDOFF.md` lists them).
- ⭐⭐ **The distinction that motivates the question:** a change to an upstream **value** genuinely couples
  the steps — a later step computed *with* that value, so its physics must be recomputed. ⛔ A change to an
  upstream **name or reference** moves no physics anywhere. ⚠ Under the merge, both cost the same: every
  downstream artifact is invalidated. ⭐ **Is that distinction real, and does it survive measurement?**
- ⚠ ⭐ **The requirement's own reliability.** An accumulated list is correct only if every intermediate copy
  is correct. A generated one is correct whenever it is generated. ⭐ Which better delivers **(a)** and
  **(b)** above, ⛔ and say what each loses.
- ⭐ **What the merge buys that a generated view does not.** ⛔ Argue this side properly: one import giving a
  consumer everything; the overwrite being visible as a diff between two committed files; knob retirement
  (`KNOB → DERIVED`, `step` moving) marked at the write. ⛔ Do not dismiss these.
- ⚠ **What breaks either way.** ⭐ Name the consumers that exist today and what each would have to change.
  ⛔ Do not assume a consumer exists; check.

## ⛔ Constraints

1. **`CLAUDE.md`** in full. ⛔ **Rule 5** — ⛔ no measured physics value, count, membership or expected
   outcome in a builder-facing artifact.
2. ⭐ **PHYSICS, ⛔ NOT PROCESS.** ⛔ Machinery that catches no way the physics could be wrong must be cut.
   ⚠ *"Simple here is usually better than more process."* (user)
3. ⭐ If this **supersedes** `D5`, `D7`, or an `F`-item of `S11_export_chain_decisions_v2.md`, ⭐ say so
   **explicitly with the reason recorded**. ⛔ Never silently. ⚠ A settled decision may be superseded; ⛔ a
   contradicted-but-still-standing one is worse than either.
4. ⭐⭐ **REPORTING IS SUCCESS.** ⛔ If the question cannot be settled from what exists, ⭐ say what is missing
   and what would settle it. ⛔ Do not invent a third mechanism to escape the choice — ⚠ every mechanism
   invented to close a gap in this workstream has been blocked by two legs.

## ⭐ Context you should read, ⛔ not inherit

`DEFECT_REGISTER.md#c20` records four blocked export-chain designs and the measurement that killed each.
⭐ Read it so you do not re-propose one. ⛔ It does **not** answer this question; ⚠ the question is upstream
of all four.

## Report back — ⛔ under 15 lines

1. The file you wrote and its line count.
2. ⭐ **The decision, in one sentence.**
3. ⭐ What it supersedes, and what that costs.
4. ⭐ What you could not settle and what would settle it.

⛔ Read-only except the file you create. ⛔ `timeout 600` on any run; ⭐ scripts under `/tmp` with their
literal stdout at named absolute paths.
