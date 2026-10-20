# S12 decision list, review pass 5 (scoped, revision 4): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_decision_list_review_p5_lookups.md`, beside this file
(generator `_scratch/s12/gen/dl_p5_lookups.sh`). "R-n" cites a section of that file.

## The pass

- **Artifact.** Revision 4 of `directives/S12_decision_list.md`, uncommitted at review on top of `45e42474`
  (revision 3, preserved). It is committed unchanged with this file. Its hash matched the baseline after the pass, and
  the working tree was unchanged apart from this file and the lookups (R1).
- **Prompt and legs.** `_scratch/s12/dl_review_prompt_p5.md`, sha256 `8bc4e5dd…` (R1). It asked whether S5-1 to S5-3
  hold, and whether the face rule, as now scoped, is consistent with the committed law and U1's encoding and leaves
  out nothing the law produces at a face. Both legs received it identically.
  - **Codex:** `gpt-6.1-sol`, effort xhigh. Report `/var/projects/_review_scratch/s12-dl-p5-codex/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-dl-p5-grok/grok_stdout.txt`.
- **Verdicts.** Both legs exited 0, and both report CLEAR, with each of S5-1 to S5-3 holding (R1, R2).

## What the legs checked

- **S5-1 (Q1).** The face rule and the OPEN item now name `n`, `u`, `J_χ` and the interior order. Nothing in the
  list can be read to forbid a face part of `n Γ_B`, or to make the order balance's face content a face condition.
- **S5-2.** Item 1 still includes every face part of `n Γ_B` and what the total-material balance requires at each
  face. The list gives no value, sign or form for either.
- **S5-3.** Every changed block belongs to S5-1, S5-2 or the gate paragraph.
- **Consistency.** Each leg derived, in its own script, what the committed law produces at a moving, tilted face for
  fields in the stated class (R3). Each finds the scoped rule consistent with the law and U1. Each finds that nothing
  the law produces from that class is left out, and that a flux concentrated on a face is an extra input, which the
  rule and the OPEN item leave out as the user chose.

## Outcome

No findings. The S12 decision list is cleared as revision 4. The builder may repair the S12 build against it. The
list is committed exactly as reviewed.
