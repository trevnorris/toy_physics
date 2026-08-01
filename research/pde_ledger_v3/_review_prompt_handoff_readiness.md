# REVIEW — handoff readiness. Read-only, one pass, deliberately narrow.

Repository root `/var/projects/toy_physics`, branch `ledger-v2-rebuild`.
**Run `git log --oneline -8` and `git status` first.** ⛔ Do not trust a hash written in any document,
including this one.

⛔⛔ **SCOPE — this is NOT a content audit.** Three prior reviews of this work died having emitted nothing
because they tried to verify every claim in the corpus. ⭐ **This one has a single question:**

> **Can a fresh session, reading only the orientation budget, pick this up and do the right thing?**

⇒ Judge the **handoff**, ⛔ not the physics. You are simulating the next session.

## How to run it

1. Read `research/pde_ledger_v3/NEXT_SESSION.md` **first and in full** — it is the handoff.
2. Read the four other documents its `ORIENTATION BUDGET` names, ⛔ and nothing else except what steps 3–5
   below require.
3. **Try to execute `S0.5`** from `research/pde_ledger_v3/V3_STEP_PLAN.md` as written — ⛔ do NOT actually
   edit anything; walk it. Open the files and code it names. Report where you would be **blocked, guessing,
   or about to do the wrong thing**.
4. **Check the handoff's factual claims against the tree**, since a stale handoff is the failure mode that
   cost this project a session already:
   - does `## ▶ WHERE WE ARE` match `git log`?
   - does each of the six residuals' stated status (`DONE` / `PARTLY` / `OPEN` / deferred) match reality?
   - do the `PENDING WORK` items and `DECISIONS OWED` describe things that are actually outstanding?
5. **Contradictions between the five budget documents.** ⭐ Especially: does any of them still say something
   another one now retracts? This corpus's recurring defect is a correction landing in one place while its
   siblings survive — it has recurred **six times** in one session, and four of the six were found only
   after a reviewer had already passed the file.

## What to report

⭐ **Rank by "would this make the next session do the wrong thing?"** ⛔ Cosmetic issues are not wanted.

Specific things worth checking, since each was decided late and may not have propagated:
- `S0.5` says the two quantities and `R3` are **DELETED, not retired**. ⛔ Does anything anywhere still say
  retire, or assume a retired-but-present row?
- `C-M1` has no mutation left once `R3` is deleted. Is that recorded where someone executing `S0.5` will
  see it, and is it clear that **choosing the replacement is the user's call**?
- `CHARTER` §1.3 (`{#falsification-standard}`) states what counts as failure, including that producing
  **more** than an incomplete reference theory is expected. Does anything else contradict it?
- The acceptance payload must be **independently re-derived, never preserved**. Is that instruction intact
  and unambiguous?

## Operating constraints

- **READ ONLY.** ⛔ Do not modify, stage, or commit. ⛔ One pass, no clarifying questions.
- ⭐ **Ship the report.** If you approach any limit, emit what you have and list what you did not reach.
  ⛔ A partial report that ships beats a thorough one that does not — three prior runs proved it.
- ⛔ For any claim of the form "there is no X", `wc -l` the artifact and read it before asserting.

## Output

```
# Handoff readiness — <your name/model>

## VERDICT
<READY / READY WITH GAPS / NOT READY> + 2-3 sentences

## WOULD CAUSE WRONG ACTION
<numbered; each: what a fresh session would do, why, file:line, what to change.
 EMPTY IS A REAL VERDICT.>

## STALE OR WRONG IN THE HANDOFF
<claims in NEXT_SESSION.md that do not match the tree>

## CONTRADICTIONS BETWEEN BUDGET DOCUMENTS

## WALKING S0.5 — where I would be blocked or guessing

## NOT REACHED
```

⛔ Do not be agreeable — ⛔ but do not invent problems. If a fresh session could pick this up cleanly, say so.
