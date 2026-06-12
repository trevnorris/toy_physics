You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u07_gain_closure_forcing`, stages 141, 142, 144, 148, 151, 156) in the INTERNAL-CONSISTENCY COMPLETENESS-CRITIC spot-check. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for internal_consistency values (THE INVERTED TEST — READ CAREFULLY)
These values are classified `constraint_kind: internal_consistency` — i.e. the ledger CLAIMS each is a DERIVED consequence (an algebraic identity / Solve / closure / extremum / limit / normalization) of upstream POSITED inputs, NOT itself an inserted choice. Unlike the free_choice sweep, the fatal flaw here is the OPPOSITE mis-classification: **a value that is genuinely a POSITED free_choice (an inserted number / a tuned knob) OR a BACK-SOLVE from an EXTERNAL published number, but is mis-classified `internal_consistency` AND labelled "derived"/"forced"/"exact"/"fixed by"/"a theorem"/"not free" in a paper CARD or per-stage NOTES line** — a posit-as-derived (or fit-as-derived) that ESCAPED the free_choice + published_target sweeps via mis-classification. These candidates were selected BECAUSE their slug carries an overclaim-flavored word (forced/fixed_exactly/backsolve/injection/not_free/_req/unique/...). The slug is the PROMPT to scrutinize, NOT evidence.
VERDICT RULE: **NO** = the value is a GENUINE derivation — open the notes/card and find the actual upstream step that PRODUCES it (a pole condition, a Solve, a closure relation, an extremum, a uniqueness criterion, a parity/symmetry forcing zero, an exact-algebra consequence of a NAMED posited input); the slug's "forced/exact" correctly names that consequence. **YES** = the value is actually inserted/tuned/back-solved (no upstream step produces it; it is chosen to hit a downstream or external target and then relabelled) — give the specific value + the card/notes file+line that calls it derived/forced/exact. Also flag (non-fatal) any stale upstream stage-anchor, value-divergence between members, or a slug that overstates vs the actual (honest) card/notes wording.

## Focus
Mouth-gain / fixed-point-closure band. M_q/R_q "no_longer_free" (141/142): verify each is fixed by the named fixed-point/closure, not a posit (these are the gains the free_choice batch-4 u-units found honestly conditioned on "once the mixed tube is identified"). r_F1 "radical_constant_forced_by_rf1" (148): r_F1 = sqrt(4107-100pi^2)/(10pi) closed form -- confirm the radical is the derived constant, not the decimal hardcode (punch-list F2 is a script-form cosmetic on 162, not here). Sigma0_can "unique_traction_renorm" (156): confirm derived renormalization. delta_Pi_act "fixed_by_covariances" (151) + Pi_star "unique_regular_canonical_branch" (144): confirm closure.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage141_gain_pair_no_longer_free (param=M_q): redteam_adversarial/provenance/fit_stage141_gain_pair_no_longer_free__m_q.yaml
- fit_stage142_three_quantities_no_longer_free (param=R_q): redteam_adversarial/provenance/fit_stage142_three_quantities_no_longer_free__r_q.yaml
- fit_stage144_unique_regular_canonical_branch (param=Pi_star): redteam_adversarial/provenance/fit_stage144_unique_regular_canonical_branch__pi_star.yaml
- fit_stage148_radical_constant_forced_by_rf1 (param=r_F1): redteam_adversarial/provenance/fit_stage148_radical_constant_forced_by_rf1__r_f1.yaml
- fit_stage151_correction_fixed_by_covariances (param=delta_Pi_act): redteam_adversarial/provenance/fit_stage151_correction_fixed_by_covariances__delta_pi_act.yaml
- fit_stage156_unique_traction_renormalization (param=Sigma0_can): redteam_adversarial/provenance/fit_stage156_unique_traction_renormalization__sigma0_can.yaml

## Sources to open as needed
Each provenance bundle carries a `constraints[].constraint_kind` + `evidence_citation` and often a `provenance_findings` block stating the Phase-B synthesis agent's RATIONALE for the internal_consistency call (with cited notes/card lines). Treat that rationale as a CLAIM TO ATTACK, not ground truth — verify it against the actual files. Paper cards: research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes: research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (EM-projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage; read how the value is LABELED and find the line that DERIVES it. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (mis-classified posit/fit dressed as derived?) | one-clause reason citing BOTH the label line AND the derivation-step line (or its absence)`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific value + the overclaiming card/notes file+line + why no upstream step derives it.
3. Classification assessment (the crux): for each member, is `internal_consistency` CORRECT (a genuine derivation exists) or should it be free_choice / published_target (a posit / external back-solve)? Name any mis-classification stage+line, or state the classification is sound.
4. Proof you looked: which card+notes files you opened for each anchor stage, the exact label-word, and the exact derivation-step line you found (or confirmed absent).
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
