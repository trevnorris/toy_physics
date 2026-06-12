# Phase C — Batch 1 adjudication record (Claude+Codex consult, 2026-06-11)

Signed-off adjudication shape (execution plan §5): structured Claude+Codex consult, `codex exec -s read-only`, user final gate on FIND_STANDS + conceptual. Consult log: `codex_logs/step5_batch_c1_adjudication_consult_codex.txt` (session `019eba37-fc91-7731-a57b-005573a2499a`); prompt: `codex_logs/step5_batch_c1_adjudication_consult_prompt.md`. Codex independently verified every cited file/line before adjudicating.

Batch 1 = 28 parameter-value families (published_target ∪ HIGH). **27 families NO. Two items routed to adjudication. NO FATAL FLAW SURVIVES** — both outcomes are non-fatal close-out fixes. Both concur with Claude's preliminary position.

## Item 1 — `fit_stage192_chi_q_unit_relations` — adversarial YES (paper_card_overclaim)
- **Defense: PARTIAL. Adjudication: PARTIAL — finding STANDS, non-fatal.**
- The card `paper/stages/stage_192.tex:13/20/21` asserts Stage 192 "collapses Packet A to Δ_Q=χ_Q−1 / N_Q=χ_Q⁻¹", but stage-192 notes are orbit/quotient projector calculus only (χ_{0,*} in the matrices) and explicitly do NOT yet compose the projector with the defect packet (note:376). χ_Q=1 is fixed at stage 194 (notes:162/184); Δ_Q:=χ_Q−1 and N_Q=χ_Q⁻¹ at stage 195 (notes:185/215). The identities are genuine downstream — just misattributed to 192 via byte-identical shared boilerplate across cards 192/193/194.
- **Disposition (Step 6 close-out, content fix → Codex applies, Claude reviews):** correct stage_192.tex (and the shared 192/193 boilerplate block) so it does not attribute χ_Q/Δ_Q/N_Q to stage 192; attribute to 194/195. Not fatal — no derivation dependency cone to hard-stop. The `fam_0058_delta_q` "benign boilerplate" counter-reading is acknowledged; adjudication: it is a real card overclaim AND non-fatal (both true).

## Item 2 — `chain_aspect_37_20` — central classification (mandatory user gate, PC3)
- **Defense: PARTIAL. Adjudication: FIND_STANDS on the classification correction — non-fatal, non-conceptual.**
- **L/a = 37/20 (and derived r_F1) is `free_choice`, NOT `published_target`.** Stage 073:54 introduces it as "a reference-branch numerical freeze, not a new theorem"; stage 121:58 derives r_F1 = √((12/π²)(37/20)²−1) deterministically from the carried value; `benchmarks.yaml` has NO 37/20 / 1.85 / 1.778 entry. The genuine external class is 54/5 (P0_target/N_Q → GR 2G/5c⁵, with benchmark IDs). Project memory's grouping of 37/20 with the 54/5 class was an error.
- **Disposition (Step 6 close-out, provenance-label cleanup):** relabel `constraint_kind` `published_target → free_choice` in the stale r_F1 bundles `fit_stage156_carried_constants_block__r_f1.yaml:191`, `fit_stage158_canonical_quartet_literals__r_f1.yaml:193`, `fit_stage161_rf1_family1_radius__r_f1.yaml:43`. Bundles 162/163 already say free_choice. `fit_stage163_Rstar_quarter` already self-flags its scanner mislabel (1/4 = R_q, not r_F1).
- **USER GATE:** surfaced to the user as the central-adjudication resolution (37/20 = posited free_choice, disclosure adequate). The "FIND_STANDS" here = the classification-correction finding stands; it is NOT a surviving fatal flaw.

## Other batch-1 close-out fixes (NO verdicts, non-fatal — Step 6)
- stage247 script wording ":239 independently derived / :248 falsifiable" overstates the tautological λ_L closure.
- Pe_fail_chi (089) false "082" upstream anchor (script 61-62/92 + CHECKPOINT_CONSTANT_PROVENANCE.md:111); the 089-notes "086" anchor is correct for the rho window.
- fam_0351 stale "Stage-18" anchor → Stage 035 (CHECKPOINT_CONSTANT_PROVENANCE.md:376).
- fam_0212 V_known bundle `published_target → free_choice` (illustrative).
- stage_189.tex bare StatusExactClosure w/o card-level GR-target disclosure; g_UV "carried from 245" stale symbol-anchor.

**Overall: layer-2 batch-1 confirms the program's external-fit surface is honestly disclosed; the single overclaim (stage_192 card) and the classification/label corrections are all non-fatal close-out items. No program revision required.**
