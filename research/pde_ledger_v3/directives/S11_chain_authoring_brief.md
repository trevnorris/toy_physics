# ⭐⭐ YOU ARE THE AUTHOR NOW — write the decision list the S11 SymPy build runs on

⚠⚠ **This is not a review.** `CLAUDE.md` **rule 15** — *if successive revisions keep breeding defects in the
material just changed, change the author* — has been invoked by the user. ⭐ **Three orchestrator-written
attempts have been blocked by two independent legs each, with the same root cause every time.** ⛔ Do not
write a fourth version of any of them. ⭐ **Author the thing yourself, from the constraints.**

You will not review this artifact. A fresh Claude agent and Grok will.

---

## The deliverable

⭐ **A decision list a builder can build from**, covering everything the S11 SymPy engine must settle to
write its export, plus whatever chain-level rule the shared key namespace needs.

⭐ **You choose the decomposition** — one document or two, and what belongs at chain level versus engine
level. ⚠ The previous split (a chain rule `F8` plus an engine list) was the orchestrator's and it produced
cross-references that both legs found under-determined. ⛔ Do not inherit it out of politeness.

⭐ **Write it to** `research/pde_ledger_v3/directives/` under names you choose. ⛔ Do not edit
`S11_SHARED_PHYSICS.md`, ⛔ any engine, ⛔ any committed export, and ⛔ any of the three blocked artifacts.

⚠⚠ **LENGTH IS THE DEFECT RATE, measured on this project:** a six-page decision list drew ~30 findings and
0 lines built; the one-page rewrite of the same material reproduced the physics 26/26 on the first build.
⇒ ⭐ **Write what must become TRUE. ⛔ Not how.**

---

## The problem, stated without a proposed solution

A chain of steps each import the previous step's flat `LEDGER` and export a merged one. Two steps analyse
**different actions** and compute **the same kinds of object** — a root ordering, an expanded Lagrangian, a
dynamical matrix. The existing keys name the kind of object and not which system it belongs to, so one key
is claimed by two different objects.

⭐ **What must remain true, and it is the reason the chain exists:** when two steps derive **one object**,
they must **meet and be compared**. ⛔ Machinery that prevents all collisions destroys the check.
⭐ **What must become true:** an object must not be mistaken for a different object of the same kind.

⛔ **Nothing more is decided.** ⭐ Whether identity belongs in the key, in the row, in a separate object, or
somewhere none of the three attempts considered, ⭐ **is yours to decide.**

---

## ⛔ The three attempts and why each died — ⛔ do not repeat any of them

Full leg reports, ⭐ read them:

| | |
|---|---|
| `directives/_legs/round2_codex_leg.txt` · `round2_grok_leg.txt` | scoped keys — **blocked** |
| `directives/_legs/round3_codex_leg.txt` · `round3_grok_leg.txt` | slug + action descriptor — **blocked** |

The blocked artifacts, ⭐ as records: `S11_py_rebuild_decisions.md` (v1, `11bf8e05`),
`S11_py_decisions_v2.md`, `S11_py_decisions_v3.md`, and `F8` in `S11_export_chain_decisions_v2.md`.
⭐ v1's block message is `git log -1 --format=%B 11bf8e05`.

⚠⚠ **The invariant across all three: each tried to make a KEY STRING carry object identity.** Your own
round-3 measurement is the sharpest statement of it — an already-wrong frozen map passed **all six** rename
controls, and `QUANTITY_NAME_HAS_NON_PAYLOAD_BINDING_EVIDENCE=False`.

⭐ **One candidate you should ATTACK, ⛔ not adopt:** put identity in the **row** as data rather than in the
key, so a key is only an index. ⛔ It is the orchestrator's idea, it has had no legs, and the orchestrator's
last three ideas were all blocked. ⭐ **If it fails, say why and write something else.**

---

## ⛔ Constraints — an artifact that violates one of these is not a deliverable

1. **`CLAUDE.md`** in full. ⭐ Rules 2, 3, 5, 6, 11 and 14 have each been broken by a previous attempt.
2. **`F1`–`F7`** of `S11_export_chain_decisions_v2.md` are settled and two-legged. ⭐ You may supersede one
   **explicitly, with the reason recorded** — ⛔ never silently. ⚠ `F8` is blocked; ⛔ treat it as a record.
   ⚠ `F4`'s "not re-keyed" clause rested on a measurement taken before the namespace was censused.
3. **The physics is `git show cf4a21a4:…/S11_SHARED_PHYSICS.md`.** ⛔ Add no physics; ⛔ name no object it
   does not name; ⛔ state nothing about what any computation produces.
4. ⛔⛔ **RULE 5 — ⛔ NO MEASURED PHYSICS IN A BUILDER-FACING ARTIFACT.** ⛔ No value, count, membership,
   sign, spectrum, movement, or expected outcome a builder could converge on. ⚠ **The last two attempts both
   leaked, the second one round after being caught.** ⭐ A property and its required operands are fine; ⛔ a
   measured result is not.
5. ⭐ **The user's standing filter: PHYSICS, ⛔ NOT PROCESS.** ⛔ An item that adds machinery without catching
   a way the physics could be wrong must be cut. ⚠ The user's words this round: *"simple here is usually
   better than more process"* and *"aim to reduce unnecessary naming collisions."*
6. ⭐ Every guard **emits both operands and the residual, then guards** — ⛔ never an assertion whose only
   outcome is `0`, ⛔ never a check that compares an artifact against its own input.
7. ⭐ **`REPORTING IS SUCCESS`** must appear in what you write, and you must obey it yourself: ⛔ if a
   question cannot be settled from what exists, ⭐ **say so and leave it open** rather than inventing a
   mechanism. ⚠ Every mechanism invented to close such a gap in this workstream has been blocked.

## ⭐ What the list must settle, at minimum

⭐ These are the questions a builder cannot proceed without. ⛔ How each is answered is yours.

- What the export carries, decidable for **tag-derived rows and for the symbol rows manufactured by
  traversal** — ⚠ a rule stated only over tags leaves the second population undecided.
- The row schema, and what happens when an imported row's shape is not what was expected. ⚠ Both a typed
  contract and run-time discovery have been measured to fail differently; ⭐ neither is presumed.
- What makes two objects **the same object**, ⭐ **total** over every value type the export admits — ⚠ most
  admitted values do not support subtraction — and ⭐ what happens when it cannot decide. ⛔ `F2` currently
  has no branch for that.
- What a re-derived row must carry so a consumer reading **only** the merged export can recompute the
  comparison (`F3`).
- How an imported row is protected from change this step did not intend, ⭐ given that any field can carry
  physics and that a guard must not audit its own input.
- What must be true before an export may be published, and what must be true of a **previously** published
  artifact after a run that did not complete.
- How a computed object that never reached the export is detected.
- Whether the committed upstream export must change, ⭐ and if so, what makes that safe — ⚠ re-pointing a
  name is the failure mode where every check in the repository passes and the physics silently moves.

## Report back — ⛔ under 20 lines

1. The files you wrote and their line counts.
2. ⭐ Which of the minimum questions you **left open**, and why it cannot be settled from what exists.
3. ⭐ Anything in the constraints you found contradictory, under-determined, or unsatisfiable.
4. ⛔ Do not summarise the physics, ⛔ do not state what anything equals, and ⛔ do not claim the design is
   correct.

⛔ Read-only on everything except the files you create. ⛔ `timeout 600` on any run; scripts under `/tmp`,
with their literal stdout at named absolute paths.
