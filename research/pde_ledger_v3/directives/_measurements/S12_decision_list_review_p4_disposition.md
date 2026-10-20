# S12 decision list, review pass 4 (scoped, revision 3): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_decision_list_review_p4_lookups.md`, beside this file
(generator `_scratch/s12/gen/dl_p4_lookups.sh`). "P-n" cites a section of that file.

## The pass

- **Why.** The S12 build's review rounds 0 and 1 found three defects with one cause: the list never said what kind of
  function each field is at a face (`S12_build_review_r{0,1}_disposition.md`). The user chose to state that in the
  list, and to leave order carried along a face out for now (2026-10-20). Revision 3 does that.
- **Artifact.** Revision 3 of `directives/S12_decision_list.md`, uncommitted at review on top of `04a24c36`. Revision
  2 is unchanged at that commit. Revision 3 is preserved unchanged in the commit that adds this file. Its hash matched
  the baseline after the pass, and the working tree was unchanged apart from this file and the lookups (P1).
- **Prompt and legs.** `_scratch/s12/dl_review_prompt_p4.md`, sha256 `c207d863…` (P1). It asked whether S4-1 to S4-5
  hold, and whether the stated class is consistent with the committed law and U1. Both legs received it identically.
  - **Codex:** `gpt-6.1-sol`, effort xhigh. Report `/var/projects/_review_scratch/s12-dl-p4-codex/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-dl-p4-grok/grok_stdout.txt`.
- **Verdicts.** Codex: CLEAR, S4-1 to S4-5 hold. Grok: one finding, High; S4-3 fails and the other four hold (P2).

## Findings

**Q1. The face rule can be read to cover the source `n Γ_B` itself.** (Grok 1, High.) **CONFIRMED, as a defect in
the wording.**
- **The list.**
  - Item 1 covers "every part concentrated on a face" (line 83).
  - The face rule ends "no field has a part concentrated on a face itself" (line 91).
  - The new OPEN item excludes "any part of a field concentrated on a face" (line 128) (P5).
- **What the legs agree on.** Each leg computed the order balance's coefficient at a moving face from fields in the
  stated class, and the two agree on it. It carries the jumps of `J_χ` and the interior traces of the order, density
  and velocity (P4).
- **Where they differ.** They differ on where that coefficient goes (P3).
  - **Codex** reads the face rule as applying to the four named fields only. The coefficient is then the face part of
    `n Γ_B`, which item 1 includes.
  - **Grok** reads "field" as covering `n Γ_B` as well. The coefficient must then vanish, as a second face condition
    beside the material one. A build following item 1 would then put a face part in `n Γ_B` that the OPEN item
    excludes, or drop the condition.
- **Which reading the list intends.** Codex's, and it is the one consistent with U1.
  - Revision 2, cleared at `a48ce830`, already included in `n Γ_B` "any part concentrated on a face if the order
    jumps there" (P5).
  - U1 says material carried past the slab's thickness de-orders, and the order may jump at a face (P5). With a jump,
    that loss happens at the face, and it is the face part of `n Γ_B`.
  - Requiring the coefficient to vanish would be a choice the list does not make: no conversion at a face. Wherever
    the order jumps, it would forbid ordered material from crossing a face unless `J_χ` cancels it.
- **Why it is still a defect.** Revision 3's text admits Grok's reading, and a builder could take it. That would
  change what is computed at each face.
- **What must be true.** The face rule applies to the fields the source is built from: `n`, `u`, `J_χ` and the
  interior order. It does not apply to `n Γ_B`, which the build computes and whose face parts are in item 1. The
  OPEN item names the same four fields. Nothing else changes.

## Outcome

One finding, confirmed. It changes what a builder could compute, so revision 3 does not clear. Revision 3 is
preserved as is. Revision 4 makes the repair above, and gets its own scoped review before the builder repairs.
