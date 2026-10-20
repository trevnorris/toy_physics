# S12 decision list, review pass 3 (scoped): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_decision_list_review_p3_lookups.md`, beside this file
(generator `_scratch/s12/gen/dl_p3_lookups.sh`). "Q-n" cites a section of that file.

## The pass

- **Scope.** This is the third round, so it was raised with the user before it ran. The user chose "Scoped review":
  check the two object repairs and the pointer-only inventories, and confirm nothing else changed.
- **Artifact.** Revision 2 of `directives/S12_decision_list.md`, uncommitted at review on top of `ea711e24`. It is
  committed unchanged with this file. Its hash matched the baseline after the pass, and the repository was unchanged
  (Q1, Q2).
- **Prompt and legs.** `_scratch/s12/dl_review_prompt_p3.md`, sha256 `c426ea06…` (Q1). It asked whether four
  properties hold, S3-1 to S3-4, and for anything in the changed text that would change the computed object or enlarge
  a claim. Both legs received it identically.
  - **Codex:** `gpt-6.1-sol`, effort xhigh. Report `/var/projects/_review_scratch/s12-dl-p3-codex/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-dl-p3-grok/grok_review.txt`.
- **Verdicts.** Both legs exited 0. Both report CLEAR, and both report that each of S3-1 to S3-4 holds (Q1).

## What the legs checked

- **S3-1 (N1).** The computed object is named for U1's support condition. The re-ordering clause is an OPEN
  admissibility requirement with no rate and no local sign. No sentence claims the object covers all of U1.
- **S3-2 (N2).** The specialization's two face-offset arguments, with the function general in both, build in no
  single one-argument choice. It remains a labelled specialization, not a limit.
- **S3-3 (N3–N6).** Every inventory pointer resolves to the row it names. The S14a item quotes the plan exactly. The
  short sentences kept beside the pointers close, drop or reassign none of the cited obligations.
- **S3-4.** Every changed block of the diff belongs to S3-1 to S3-3 or to the gate paragraph.

## Outcome

No findings. The S12 decision list is cleared as revision 2: the decision-list gate (G2) and the spec row are both
met. The builder may start from it. The list is committed exactly as reviewed. Its gate paragraph still reads as
awaiting review, and this disposition and its commit record the clearance.
