# S12 record decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_record_decision_list_review_p1_lookups.md`, beside this file
(generator `_scratch/s12/gen/recdl_p1_lookups.sh`). "K-n" cites a section of that file.

## The pass

- **Artifact.** `directives/S12_record_decision_list.md`, as first written, uncommitted at review on top of
  `fe485b73`. It is preserved unchanged in the commit that adds this file. Its hash matched the baseline after the
  pass, and the working tree was unchanged apart from this file and the lookups (K1).
- **Prompt and legs.** `_scratch/s12/rec/recdl_review_prompt_p1.md`, sha256 `80ec80a3…` (K1). Both legs received it
  identically, in parallel.
  - **Codex:** `gpt-6.1-sol`, effort xhigh. Report `/var/projects/_review_scratch/s12-recdl-p1-codex/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-recdl-p1-grok/REVIEW.md`.
- **Verdicts.** Codex: CLEAR. Grok: FINDINGS, 2 (K2).
- **Gate.** This is a decision list under G2: one pass, findings verified, folded once, then the record author starts.

## Findings

**Q1. Item 5 gives every premise one status, which the step's decision list does not state.** (Grok 1, High.)
**CONFIRMED.**
- **The list.** Item 5 asks, for each of U1–U4, the "Fields at a face" rule and the fast-ordering assumption, for "its
  status (adopted, conditional, re-tested at the knit)" (K3).
- **The gap.**
  - The step's decision list gives these different standings (K3). U1–U4 are each "an adopted substrate input to a
    conditional model". How U1 enters is "an identification, flagged". U1's re-ordering clause is "not encoded" and
    stays an OPEN admissibility requirement. The face rule is a statement of the fields' function class. Item 2 is "a
    labelled conditional specialization".
  - The step's decision list never mentions a knit (K3).
  - Read as one template, item 5 would have the record give all of them one status, the unencoded re-ordering clause
    included.
- **What must be true.** Each of these is reported with the standing the step's decision list gives it, quoted from
  that list, and with nothing it does not give.

**Q2. Item 9 misses a status line that says S12 waits on the closeout's review.** (Grok 2, Medium.) **CONFIRMED.**
- **The list.** Item 9 updates lines that call the closeout "review pending" or "under review" (K4).
- **The gap.** `STATUS.md:88` says "S12 starts once the light closeout's review is clear." It uses neither phrase, so
  item 9 does not reach it (K4). The acceptance commit `e070c7f3` did not edit `STATUS.md` (K4).
- **What must be true.** Every line in those three files that describes the light-sector closeout as awaiting review,
  in whatever words, states its acceptance at `e070c7f3` instead.

## On the CLEAR leg

Codex found nothing under the filter. Both findings rest on the list's text, the step's decision list and the status
files themselves (K3, K4).

## Outcome

Two findings, confirmed. The list as reviewed is preserved, and the fold makes each "What must be true" hold. That is
the one fold G2 allows. The record author starts after it.
