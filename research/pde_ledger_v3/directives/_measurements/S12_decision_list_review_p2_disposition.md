# S12 decision list, review pass 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S12_decision_list_review_p2_lookups.md`, beside this file
(generator `_scratch/s12/gen/dl_p2_lookups.sh`). "M-n" cites a section of that file.

## The pass

- **Artifact.** Revision 1 of `directives/S12_decision_list.md`, the pass-1 fold, uncommitted at review on top of
  `cd505a87`. It is preserved unchanged in the commit that adds this file. Its hash matched the baseline after the
  pass, and the repository was unchanged (M1).
- **Prompt and legs.** `_scratch/s12/dl_review_prompt_p2.md`, sha256 `b35adc41…` (M1). It asked for a review of the
  whole revised list and a check of each pass-1 repair. Both legs received it identically.
  - **Codex:** `gpt-6.1-sol`, effort xhigh. Report `/var/projects/_review_scratch/s12-dl-p2-codex/REVIEW.md`.
  - **Grok:** `grok-4.7`, effort high. Report `/var/projects/_review_scratch/s12-dl-p2-grok/grok_review.txt`.
- **Verdicts.** Codex: FINDINGS, 1. Grok: FINDINGS, 5 (M1).
- **What holds.** Both legs report the balances faithful and homogeneous. They find that the target-plus-deviation
  encoding, the earned-limit wording, the resting-brane wording and the galactic edge are gone, and that the two power
  obligations stay separate. Codex reports each of K2–K6 repaired at its site, finds no answer in the object
  definition and no further freeze, and finds the object able to disagree across engines.

## Findings

**N1. The support condition leaves out U1's re-ordering clause.** (Codex 1, Medium.) **CONFIRMED.**
- **The list.** Outside the region `χ_B` vanishes; inside it is a live field. The object is the source "on order
  fields U1 allows" (M2).
- **The gap.** U1 also says "Bulk material that settles into the slab re-orders" (M2), and return's committed meaning
  is bulk material re-ordering into the brane state (M2). The support condition admits material that settles inside
  and never re-orders. U4's global sign does not exclude that.
- **Repair.** Name the computed object for what it uses: U1's support condition. Keep the re-ordering clause as an
  OPEN admissibility requirement. Choose no rate and no local sign.

**N2. The fast-ordering specialization does not fix its arguments.** (Grok 1, Medium.) **CONFIRMED.**
- **The list.** "A given function of position relative to the faces", with the faces "equivalently a live centre and
  thickness" (M3).
- **The gap.** Different face-relative arguments, such as a fraction of the thickness or an offset from the centre,
  give different sources. A builder choosing one would add an unflagged assumption.
- **Repair.** State the specialization's arguments so that no single choice is built in.

**N3. "The return law" widens S14a's item.** (Grok 2, Medium.) **CONFIRMED.**
- **The plan.** S14a derives "the controlled `ω → 0` return law" (`V3_STEP_PLAN.md:636–637`); O2's S14a row says
  "controlled DC return bridge" (M4). The return controllers stay in S12's source inventory.
- **Repair.** Use the plan's words.

**N4. The IR/return gloss drops part of the cited row.** (Grok 3, Medium.) **CONFIRMED.**
- **The source.** The cited row also has extrapolation, the treatment of zero modes and moduli, and the
  outgoing/absorbing prescription (`docs/two_throat_simulation_handoff_spec.md:340`). Stage 045 places the return
  controller in the zero-mode complement (`:50`) (M5).
- **Repair.** Point at the row; do not paraphrase it.

**N5. The bulk-direction carry is not S12's alone.** (Grok 4, Medium.) **CONFIRMED.**
- **The sources.** It is an OPEN reduction through `𝒥_map` and the normal material response `𝒩_br^live`
  (`O2_SHARED_PHYSICS.md:246, 584`), whose owners are listed at `:134`. The O2 record gives S12 the additional
  partners and reaction system, and gives the native reduction to O6 and the normal response (`O2_steady_brane_balance.md:396`)
  (M6).
- **Repair.** Point at those rows.

**N6. The S9b/O2 domain sentence drops restrictions.** (Grok 5, Low.) **CONFIRMED.**
- **The source.** The sentence also has a spherically symmetric mass, the far field, linear optical waves and the
  leading eikonal (`O2_steady_brane_balance.md:74–76`) (M7).
- **Repair.** Point at the sentence.

## The pattern, and the change of method

All six findings are in material the pass-1 fold changed. N3–N6 are paraphrases of sources that dropped or widened
their content. The inventories are a record-stage deliverable that the builder does not use. In revision 2 the list
therefore points at each source row instead of restating it, and the record writes the inventories under their own
review.

## Outcome

Six findings, all confirmed. Revision 1 is preserved as is. N1 and N2 change what the build computes, so revision 2
needs review before any builder. That would be the third review round, so it is raised with the user before it runs.
