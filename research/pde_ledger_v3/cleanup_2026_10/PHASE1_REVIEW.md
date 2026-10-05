# Phase 1 review (Claude) and Phase 2 instructions

**Reviewed:** `FINDINGS.md` at `174f4fa9` and `INVENTORY.tsv` at `91f560e1`. **Reviewer:** Claude (orchestrator).
**Verdict:**
- Accept Part B as a cited reading guide, and C as the conflict list.
- Part A fails directive rule 4, and so does the document's goal of giving the user a clear picture. Fix it first,
  as step 0 of Phase 2.
- Part D is approved with the amendments below.
- Part E is superseded by the tiered keep rules now in directive Phase 3a.

## 0. Fix the opening of FINDINGS first

The user is a programmer. Today's Part A reads as qualifiers stacked on jargon: "selected face drives",
"LAB_HELD/RHO4_CONSTANT", "normalized/strong-edge observable". Replace the top of FINDINGS with a section
**"Where the light sector stands"** of 15 lines or fewer, in plain words. It must answer five questions:

1. What did S11c ask? In one sentence.
2. What is now known about light on a **uniform** brane? Include how strong the evidence is: dual-engine or single
   engine, and reviewed how.
3. What is known about light where the brane **changes**? Is there a mechanism? Is there a size?
4. What is wrong with, or open in, the equations themselves? That is, the repairs and the c2 term.
5. Which later step owns each open question, and what is next?

Jargon may appear only after the plain sentence that says what it means. Keep every qualification, in plain words:
"one engine, not cross-checked" rather than "single-engine CONDITIONAL". Leave Part B as the detailed version
underneath.

**Register example.** This shows the register only. Don't copy it; get the content from the records:
> On a perfectly uniform brane, light does not couple to the bulk at all, at linear order. Both engines derived this
> independently (S11b). A later one-engine check found it still holds when the light and sound speeds are equal.

## Answers to the Phase 1 conflicts (C1–C8)

- **C1 (status versus queue).**
  - Phase 2 replaces `STATUS.md` entirely with the short front door. No old clause survives; they are in the tag.
  - In `V3_STEP_PLAN.md`, rewrite the S11b/S11c table and the S11c-a…e sub-step statuses to the current state.
    Don't leave a "NEXT" or "in flight" anywhere in S11.
- **C2 (old b/c2 closes versus repaired operators).**
  - Leave each record's original verification text as it is.
  - Add a section called **"Changed after close"** that lists each repair (commit, what changed, scope of its
    review, literal verdicts).
  - State plainly that the original verification describes the pre-repair version. The repaired version carries
    only the scoped review. For c2 this includes "VALUES are unaffected": that statement was about the pre-repair
    operator.
- **C3 (mixed-term wording).**
  - Use the confirmed source fact, not the stronger closeout wording: the c2 fold sets the direct three-leg
    `[0,2]` entry to zero while keeping the iterated product
    (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:400–418`, confirmed in the joint disposition).
  - The kept grades include `eta*sigma_W`, so order counting does not justify that zero.
  - **Whether the complete closed response has a nonzero direct term is UNRESOLVED.**
  - Correct the closeout's "It omits a direct mixed height–slope contribution", which assumes such a term exists.
- **C4 (no-leak scope versus live conversion).** Nothing transfers. Record the clean-condition result as a measured
  symmetry selection rule on one engine. Its no-leak zero needs the drain frozen, so it is OPEN and goes to S12.
  Don't extend it.
- **C5 (MacCullagh framing).** The ledger records hold the claim boundary:
  - **Uniform:** confinement is derived for the stated linear model.
  - **Non-uniform:** not established.
  - **Material admissibility:** open.

  In `docs/s11_maccullagh_differentiation.md`, **rewrite** only the sentences that claim more than that, in place.
  ⛔ Don't add a banner or annotation.
- **C6 (exploratory review disagreements).** Don't adjudicate. The throat/EM notes stay exploratory and paused, and
  are listed as input to S22/Q2. Each assessment's own "claims not accepted" section is enough.
- **C7 (paper).**
  - Add an S11c section to Part I, after S11b. It needs a visible PARTIAL statement and its limits in the main text,
    not in the suppressed Verification field.
  - Build the PDF and report the build command with its literal tail output.
- **C8 (Lean stale prose).** Rewrite the stale "pending" and "in progress" text to match the closure reviews.
  Application debts, which are uses of the theorems in S11c-d, are listed separately from theorem review status.

## Part D: approved, with these amendments

- **D1, the S11c-d step record.**
  - Approved. Keep it to about 150 lines or fewer, and give it the same plain opening as above.
  - It must also cover workstream 06, the clean-condition packet, as the C4 statement.
  - Link evidence rather than restating the chronology. It contains no process history (launches, rounds, resumes).
- **D2, the closeout.** Approved. The closeout is the S11c umbrella: about 60 lines or fewer, pointing at the a, b,
  c1, c2 and d records. Apply C3.
- **D3, the b and c2 records.** Approved, in the form given in C2.
- **D4, the plan and STATUS.**
  - Approved. Keep STATUS to about 80 lines or fewer.
  - Keep the plan's real dependencies (D4's list). Name S12 as next.
  - Mention the cleanup itself in one line: the history lives in the archive tag, and the cleanup record is under
    `cleanup_2026_10/`.
- **D5, S10/S11 versus Lean.** Approved, but only for stale or contradicting statements. No expansion of the
  records.
- **D6, the paper.** Approved, as in C7.
- **D7.** Approved, as in C5 and C6. No muonium edits beyond keeping its ownership link.

## Phase 2 STOP report

The report should list:
- the files changed, each with its line count before and after;
- the PDF build result;
- any place where applying these answers needed a judgement call.

Ground rules 1–10 still apply. Then STOP for Claude's review, with Grok as the second leg.
