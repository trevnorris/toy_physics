# Independent review — the direct-import architecture, ROUND 2

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifact — orchestrator-written, second version

`/var/projects/toy_physics/research/pde_ledger_v3/EXPORT_CHAIN_ARCHITECTURE.md`

⛔ **Nothing has been built or deleted.** ⭐ This review decides whether it gets built.

## ⭐ What changed since round 1, and what you should NOT re-litigate

⭐ Round 1 was reviewed by two legs and **rejected**. Its findings are accepted, not re-argued:

- ⛔ no engine can emit the export record v1 specified (`named c_gamma/c_L: 0`, `class field: 0`,
  `source_locus: 0` in five of six outputs) ⇒ **the export record is gone**; steps now chain by **direct
  Python import** and the objects never leave SymPy
- ⛔ the bridging adapter would be *"the unacknowledged third list"* ⇒ **there is no adapter**
- ⛔ `import ⊄ export` is a spelling check ⇒ **replaced** by a provenance requirement
- ⛔ the comparison list does not get short; S10's size is structural ⇒ **conceded in the cost section**
- ⛔ premise-only objects have no cross-engine contention ⇒ **re-framed**: they are free knobs, and the
  control on them is a register, not contention

⭐ **Attack the new version.** ⛔ Do not spend the review restating round 1.

## ⭐⭐ What to attack

**① The free-knob register — is the diagnosis right?** ⭐ The proposal claims v3 has **none**, citing
`ANSATZ_LEDGER.md` living only under the old ledger, `kind: parameter` being a type rather than a
provenance, and `provenance_status` describing relations rather than knobs. ⭐ **Verify that**, and then
answer the harder question: **can the set of free knobs be COMPUTED** — every free symbol in a step that is
neither derived in-step nor imported — ⛔ or must it be hand-declared? ⚠ If it can be computed, the register
is a report and this is cheap; ⭐ if it must be hand-written, it is another maintained list and the proposal
should say so.

**② Direct import — does it actually work on these engines?** ⭐ Look at the real scripts. Can the
derivation be wrapped in a function without restructuring the physics? ⭐ What are the import-time side
effects — module-level computation, global state, `emit` calls at import? ⭐ Does `S10-py`'s unconditional
registry reload survive, and what breaks when it does not? ⭐ Give the concrete blocker list per engine.

**③ The provenance check is the load-bearing new item and is unscoped.** ⭐ *"An object claimed `derived`
traces to this step's own action / EL path."* ⭐ **Is that mechanically checkable at all** in SymPy — and in
Wolfram? ⭐ If yes, sketch what it measures. ⛔ If it is not checkable, say so plainly: the proposal then
rests on an unbuildable check and the blind-derivation guarantee is prose.

**④ Does narrowing `relations.yaml` to knob-constraints hold up?** ⭐ The claim is that `R4`/`R5` become
dataflow while `R1`/`R2` stay, because **no v3 step derives any medium quantity**. ⭐ Verify that from the
repository. ⭐ Then judge: is a relation among undetermined quantities a **meaningful, checkable** statement,
or bookkeeping? ⛔ And does the homogeneity gate retain any discriminating power over just those three?

**⑤ Re-run cost.** ⭐ The proposal claims S11 → S10 → S9 is ~37 minutes and calls it tolerable. ⭐ Check the
arithmetic and the growth: what does this become at 20+ steps, and does the chain fan out or stay linear?
⛔ Is "cacheable later" a real option here, or does caching reintroduce a stored intermediate — the exact
thing this design removes?

**⑥ What does this still fail to catch?** ⭐ The proposal explicitly disclaims common-mode errors and says
contention exists only for objects both engines derive. ⭐ **Enumerate, from the committed outputs, which
cross-step objects each engine actually derives** — and therefore which are genuinely under contention and
which are not. ⚠ That list is the honest scope of the whole design and it does not exist yet.

**⑦ Migration order.** ⭐ It proposes: register first, then fix `R4`/`R5`'s binding to the squared speed both
engines already emit, then pilot import on one step pair. ⭐ Is that the right order, and is step 2 really
available today without engine changes? ⛔ If `R4` can be closed by a config change alone, say so — that is
the shortest path to the project's stated requirement and it should not wait behind an architecture.

## Method

⭐ Read the repository first and form your own view. ⭐ **Then** the proposal. ⛔ Not in the other order.

- ⛔ Do not modify the repository. ⛔ Do not build anything.
- ⭐ **Absolute paths outside the repository.** ⚠ A `cd` into a temp directory has failed four times in this
  session, once editing live files.
- ⭐ Literal command and output for every claim about what is or is not emitted.
- ⛔ One Mathematica kernel at a time (two seats); ⛔ `timeout 900` each.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it adopts.

⭐ Worth most: *"the provenance check is not buildable"* · *"the free-knob set cannot be computed"* ·
*"direct import breaks on engine X for reason Y"* · *"R4 can/cannot be closed by config alone"* ·
*"object Z crosses steps and neither engine independently derives it."*

⚠ **This project's failure mode is adding machinery, not omitting it.** ⭐ A finding that a piece of this is
unnecessary is worth as much as one that it is insufficient.
⚠ A leg returning *"nothing survives the filter"* is weak evidence alone — state what you checked, what you
could not, and what would have had to be true to find something.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
