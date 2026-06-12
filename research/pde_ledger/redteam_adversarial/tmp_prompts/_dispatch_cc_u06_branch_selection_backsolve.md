You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u06_branch_selection_backsolve`, stages 125, 126, 127, 128, 130) in the INTERNAL-CONSISTENCY COMPLETENESS-CRITIC spot-check. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for internal_consistency values (THE INVERTED TEST — READ CAREFULLY)
These values are classified `constraint_kind: internal_consistency` — i.e. the ledger CLAIMS each is a DERIVED consequence (an algebraic identity / Solve / closure / extremum / limit / normalization) of upstream POSITED inputs, NOT itself an inserted choice. Unlike the free_choice sweep, the fatal flaw here is the OPPOSITE mis-classification: **a value that is genuinely a POSITED free_choice (an inserted number / a tuned knob) OR a BACK-SOLVE from an EXTERNAL published number, but is mis-classified `internal_consistency` AND labelled "derived"/"forced"/"exact"/"fixed by"/"a theorem"/"not free" in a paper CARD or per-stage NOTES line** — a posit-as-derived (or fit-as-derived) that ESCAPED the free_choice + published_target sweeps via mis-classification. These candidates were selected BECAUSE their slug carries an overclaim-flavored word (forced/fixed_exactly/backsolve/injection/not_free/_req/unique/...). The slug is the PROMPT to scrutinize, NOT evidence.
VERDICT RULE: **NO** = the value is a GENUINE derivation — open the notes/card and find the actual upstream step that PRODUCES it (a pole condition, a Solve, a closure relation, an extremum, a uniqueness criterion, a parity/symmetry forcing zero, an exact-algebra consequence of a NAMED posited input); the slug's "forced/exact" correctly names that consequence. **YES** = the value is actually inserted/tuned/back-solved (no upstream step produces it; it is chosen to hit a downstream or external target and then relabelled) — give the specific value + the card/notes file+line that calls it derived/forced/exact. Also flag (non-fatal) any stale upstream stage-anchor, value-divergence between members, or a slug that overstates vs the actual (honest) card/notes wording.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It sits on prior-overclaim / mis-classification-risk ground (see Focus: e.g. a literal "backsolve"/"injection"/"5pn data"/"calibration_match" slug, or the stage-192 / chi_Q neighborhoods). Be especially literal: quote the EXACT card/notes line and its label-word, and the EXACT line that shows the derivation step (or its absence). The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
THE genuine posit-as-derived RISK surface: slugs literally say "backsolve". xi_star "admixture_backsolve" (126), x_star_exp/x_star_slab "backsolve" (127): a back-solve from a target is EXACTLY the fit-as-derived failure mode -- scrutinize HARD. Open notes127/126: is x_star/xi_star back-solved from a DOWNSTREAM target value (-> should be free_choice or published_target, a mis-classification YES), or is it the closed-form solution of an upstream equation (-> genuine internal_consistency, NO)? "lower_branch_unique/ uniquely_selected" (125/128 g_star): a branch SELECTION by a uniqueness criterion is a derivation IF the criterion is upstream; confirm it is not a hand-pick. r_F1 "fixed_by_geometry" (128) ties to the 37/20 geometry. Pi_star "exact_unique" (130).

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage125_lower_branch_unique (param=branch_sign): redteam_adversarial/provenance/fit_stage125_lower_branch_unique__branch_sign.yaml
- fit_stage126_xi_star_admixture_backsolve (param=xi_star): redteam_adversarial/provenance/fit_stage126_xi_star_admixture_backsolve__xi_star.yaml
- fit_stage127_slab_depth_unique (param=x_star_slab): redteam_adversarial/provenance/fit_stage127_slab_depth_unique__x_star_slab.yaml
- fit_stage127_x_star_exp_backsolve (param=x_star_exp): redteam_adversarial/provenance/fit_stage127_x_star_exp_backsolve__x_star_exp.yaml
- fit_stage127_x_star_slab_backsolve (param=x_star_slab): redteam_adversarial/provenance/fit_stage127_x_star_slab_backsolve__x_star_slab.yaml
- fit_stage128_lower_branch_uniquely_selected (param=g_star): redteam_adversarial/provenance/fit_stage128_lower_branch_uniquely_selected__g_star.yaml
- fit_stage128_rF1_fixed_by_geometry (param=r_F1): redteam_adversarial/provenance/fit_stage128_rf1_fixed_by_geometry__r_f1.yaml
- fit_stage130_pi_star_exact_unique (param=Pi_star): redteam_adversarial/provenance/fit_stage130_pi_star_exact_unique__pi_star.yaml

## Sources to open as needed
Each provenance bundle carries a `constraints[].constraint_kind` + `evidence_citation` and often a `provenance_findings` block stating the Phase-B synthesis agent's RATIONALE for the internal_consistency call (with cited notes/card lines). Treat that rationale as a CLAIM TO ATTACK, not ground truth — verify it against the actual files. Paper cards: research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes: research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (EM-projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage; read how the value is LABELED and find the line that DERIVES it. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (mis-classified posit/fit dressed as derived?) | one-clause reason citing BOTH the label line AND the derivation-step line (or its absence)`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific value + the overclaiming card/notes file+line + why no upstream step derives it.
3. Classification assessment (the crux): for each member, is `internal_consistency` CORRECT (a genuine derivation exists) or should it be free_choice / published_target (a posit / external back-solve)? Name any mis-classification stage+line, or state the classification is sound.
4. Proof you looked: which card+notes files you opened for each anchor stage, the exact label-word, and the exact derivation-step line you found (or confirmed absent).
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
