You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u02_early_thresholds_B`, stages 038, 041, 044, 045, 048, 049, 052, 062, 067) in the INTERNAL-CONSISTENCY COMPLETENESS-CRITIC spot-check. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for internal_consistency values (THE INVERTED TEST — READ CAREFULLY)
These values are classified `constraint_kind: internal_consistency` — i.e. the ledger CLAIMS each is a DERIVED consequence (an algebraic identity / Solve / closure / extremum / limit / normalization) of upstream POSITED inputs, NOT itself an inserted choice. Unlike the free_choice sweep, the fatal flaw here is the OPPOSITE mis-classification: **a value that is genuinely a POSITED free_choice (an inserted number / a tuned knob) OR a BACK-SOLVE from an EXTERNAL published number, but is mis-classified `internal_consistency` AND labelled "derived"/"forced"/"exact"/"fixed by"/"a theorem"/"not free" in a paper CARD or per-stage NOTES line** — a posit-as-derived (or fit-as-derived) that ESCAPED the free_choice + published_target sweeps via mis-classification. These candidates were selected BECAUSE their slug carries an overclaim-flavored word (forced/fixed_exactly/backsolve/injection/not_free/_req/unique/...). The slug is the PROMPT to scrutinize, NOT evidence.
VERDICT RULE: **NO** = the value is a GENUINE derivation — open the notes/card and find the actual upstream step that PRODUCES it (a pole condition, a Solve, a closure relation, an extremum, a uniqueness criterion, a parity/symmetry forcing zero, an exact-algebra consequence of a NAMED posited input); the slug's "forced/exact" correctly names that consequence. **YES** = the value is actually inserted/tuned/back-solved (no upstream step produces it; it is chosen to hit a downstream or external target and then relabelled) — give the specific value + the card/notes file+line that calls it derived/forced/exact. Also flag (non-fatal) any stale upstream stage-anchor, value-divergence between members, or a slug that overstates vs the actual (honest) card/notes wording.

## Focus
Threshold / no-longer-free / self-duality band. "no_longer_free" (044 xi, 049 zeta_n_phys) must be conditioned on a NAMED upstream determination (verify the notes show what fixes it). zeta_req appears twice at 048 (threshold + unique) -- same value, confirm both are the SAME derived threshold, not two posits. Theta "not_free_entropy_constant" (062): verify it is the derived entropy constant, not a tuned value. 067 self-duality w_g/w_f "forced" + C_res "unique_global_max": confirm extremum/duality derivation, not insertion.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage038_delta_fixed_by_eps_eta (param=delta): redteam_adversarial/provenance/fit_stage038_delta_fixed_by_eps_eta__delta.yaml
- fit_stage041_n_req_exact_support_loading (param=n_req): redteam_adversarial/provenance/fit_stage041_n_req_exact_support_loading__n_req.yaml
- fit_stage044_xi_no_longer_free (param=xi): redteam_adversarial/provenance/fit_stage044_xi_no_longer_free__xi.yaml
- fit_stage045_tracking_forced (param=R_U): redteam_adversarial/provenance/fit_stage045_tracking_forced__r_u.yaml
- fit_stage048_zeta_req_threshold (param=zeta_req): redteam_adversarial/provenance/fit_stage048_zeta_req_threshold__zeta_req.yaml
- fit_stage048_zeta_req_unique (param=zeta_req): redteam_adversarial/provenance/fit_stage048_zeta_req_unique__zeta_req.yaml
- fit_stage049_zeta_phys_no_longer_free (param=zeta_n_phys): redteam_adversarial/provenance/fit_stage049_zeta_phys_no_longer_free__zeta_n_phys.yaml
- fit_stage052_asymmetry_forced (param=Omega_asym): redteam_adversarial/provenance/fit_stage052_asymmetry_forced__omega_asym.yaml
- fit_stage062_theta_not_free_entropy_constant (param=Theta): redteam_adversarial/provenance/fit_stage062_theta_not_free_entropy_constant__theta.yaml
- fit_stage067_C_res_unique_global_max (param=C_res_squared): redteam_adversarial/provenance/fit_stage067_c_res_unique_global_max__c_res_squared.yaml
- fit_stage067_w_ratio_self_duality_forced (param=w_g_over_w_f): redteam_adversarial/provenance/fit_stage067_w_ratio_self_duality_forced__w_g_over_w_f.yaml

## Sources to open as needed
Each provenance bundle carries a `constraints[].constraint_kind` + `evidence_citation` and often a `provenance_findings` block stating the Phase-B synthesis agent's RATIONALE for the internal_consistency call (with cited notes/card lines). Treat that rationale as a CLAIM TO ATTACK, not ground truth — verify it against the actual files. Paper cards: research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes: research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (EM-projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage; read how the value is LABELED and find the line that DERIVES it. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (mis-classified posit/fit dressed as derived?) | one-clause reason citing BOTH the label line AND the derivation-step line (or its absence)`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific value + the overclaiming card/notes file+line + why no upstream step derives it.
3. Classification assessment (the crux): for each member, is `internal_consistency` CORRECT (a genuine derivation exists) or should it be free_choice / published_target (a posit / external back-solve)? Name any mis-classification stage+line, or state the classification is sound.
4. Proof you looked: which card+notes files you opened for each anchor stage, the exact label-word, and the exact derivation-step line you found (or confirmed absent).
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
