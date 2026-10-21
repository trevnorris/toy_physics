# S14a part 1 decision list, narrow pass 3: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_decision_list_review_p3_lookups.md`, beside this file
(generator `_scratch/s14a/gen/dl_p3_lookups.sh`). "K-n" cites a section of that file.

## The pass

- **Artifact.** Revision 3 of `directives/S14a_1_decision_list.md`, sha256 `74e74a37…`. It is the repair of pass 2's
  Q1 under the user's "Fix, narrow check" (`S14a_1_decision_list_review_p2_disposition.md`, committed `b2c1956d`).
  It matched the baseline after the pass, and the working tree was unchanged (K1).
- **Prompt and legs.** `_scratch/s14a/dl_review_prompt_p3.md`, sha256 `8d69e39a…` (K1). It was scoped to three
  checks:
  - Q1 is repaired;
  - nothing but the Q1 line and the revision header moved since revision 2 (`ccfb8e8e`);
  - the edit introduces no new defect.

  Both legs received it identically, in parallel. Codex (xhigh) wrote
  `/var/projects/_review_scratch/s14a1-dl-p3-codex/REVIEW.md`, and Grok (`grok-4.7`, effort high) wrote
  `/var/projects/_review_scratch/s14a1-dl-p3-grok/REVIEW.md`. Both exited 0 (K1).
- **Verdicts.** Codex: CLEAR. Grok: CLEAR (K2).

## What the legs established

- **Q1 is repaired.** The list prints "the face channel minus the order lost" once. The old wording occurs zero times
  (K3). Both legs show that the sentence has one reading, and that the order lost keeps P3's signed definition.
- **Nothing else moved.** The diff from revision 2 is 5 insertions and 3 deletions (K3), in two hunks: the header and
  the Q1 phrase. Codex's byte comparison prints that revision 3 equals revision 2 with exactly those two replacements.
  Grok's prints that the remainders are equal, that the order lost's definition is unchanged, and that the intended
  claim is unchanged (K3).

## Adjudication

No finding. Nothing outstanding changes what is computed or what may be claimed. That is the substantive exit
(`CLAUDE.md` G4), and it clears the list's spec content. G2 is met: one two-leg pass (pass 1), folded once (revision
2). Pass 2 then reviewed the fold until clear, and this narrow pass checked its one repair. Revision 3 is the list the
build directive is derived from.
