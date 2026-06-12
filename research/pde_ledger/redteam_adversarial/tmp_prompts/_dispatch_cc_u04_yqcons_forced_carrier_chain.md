You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u04_yqcons_forced_carrier_chain`, stages 091, 092, 093, 095, 097, 098, 099, 100) in the INTERNAL-CONSISTENCY COMPLETENESS-CRITIC spot-check. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for internal_consistency values (THE INVERTED TEST — READ CAREFULLY)
These values are classified `constraint_kind: internal_consistency` — i.e. the ledger CLAIMS each is a DERIVED consequence (an algebraic identity / Solve / closure / extremum / limit / normalization) of upstream POSITED inputs, NOT itself an inserted choice. Unlike the free_choice sweep, the fatal flaw here is the OPPOSITE mis-classification: **a value that is genuinely a POSITED free_choice (an inserted number / a tuned knob) OR a BACK-SOLVE from an EXTERNAL published number, but is mis-classified `internal_consistency` AND labelled "derived"/"forced"/"exact"/"fixed by"/"a theorem"/"not free" in a paper CARD or per-stage NOTES line** — a posit-as-derived (or fit-as-derived) that ESCAPED the free_choice + published_target sweeps via mis-classification. These candidates were selected BECAUSE their slug carries an overclaim-flavored word (forced/fixed_exactly/backsolve/injection/not_free/_req/unique/...). The slug is the PROMPT to scrutinize, NOT evidence.
VERDICT RULE: **NO** = the value is a GENUINE derivation — open the notes/card and find the actual upstream step that PRODUCES it (a pole condition, a Solve, a closure relation, an extremum, a uniqueness criterion, a parity/symmetry forcing zero, an exact-algebra consequence of a NAMED posited input); the slug's "forced/exact" correctly names that consequence. **YES** = the value is actually inserted/tuned/back-solved (no upstream step produces it; it is chosen to hit a downstream or external target and then relabelled) — give the specific value + the card/notes file+line that calls it derived/forced/exact. Also flag (non-fatal) any stale upstream stage-anchor, value-divergence between members, or a slug that overstates vs the actual (honest) card/notes wording.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It sits on prior-overclaim / mis-classification-risk ground (see Focus: e.g. a literal "backsolve"/"injection"/"5pn data"/"calibration_match" slug, or the stage-192 / chi_Q neighborhoods). Be especially literal: quote the EXACT card/notes line and its label-word, and the EXACT line that shows the derivation step (or its absence). The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The Y_Q_cons / Yhat_Q_cons "forced carrier" module chain (091-100) feeding the chi_Q quotient. Heavy "forced"/"module_forced"/"forced_carrier" slug density. This is the b052471 chi_Q-provenance neighborhood (chi_Q=1 fixed DOWNSTREAM at 105, not here). Two risks: (1) K_geom "forced_by_identity" (091) -- verify it is forced by a genuine algebraic identity, NOT the K_geom=3K_pole tension that 092 FREES (punch-list E2: 091 says "forced", 092 says "no longer forced unless"; that tension is OPENLY logged, not a flaw, but confirm); (2) the carrier-ledger-label "forced" wording (092/095/ 097/098) is a LABEL on a carried module, not a posited number -- verify the module value is derived/carried, not inserted. N_Q "factorization_forced" (100) + Gamma_5 "forced_by_yhatQcons_script" (099): confirm script-forcing is an algebraic factorization, not a back-solve.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage091_K_geom_forced_by_identity (param=K_geom): redteam_adversarial/provenance/fit_stage091_k_geom_forced_by_identity__k_geom.yaml
- fit_stage091_module_forced_paper (param=Y_Q_cons): redteam_adversarial/provenance/fit_stage091_module_forced_paper__y_q_cons.yaml
- fit_stage092_forced_carrier_ledger_label (param=Y_Q_cons): redteam_adversarial/provenance/fit_stage092_forced_carrier_ledger_label__y_q_cons.yaml
- fit_stage093_module_forced_statement (param=Y_Q_cons): redteam_adversarial/provenance/fit_stage093_module_forced_statement__y_q_cons.yaml
- fit_stage095_forced_carrier_ledger_label (param=Y_Q_cons): redteam_adversarial/provenance/fit_stage095_forced_carrier_ledger_label__y_q_cons.yaml
- fit_stage097_yhatQcons_forced_carrier (param=Yhat_Q_cons): redteam_adversarial/provenance/fit_stage097_yhatqcons_forced_carrier__yhat_q_cons.yaml
- fit_stage098_forced_carrier_card (param=Yhat_Q_cons): redteam_adversarial/provenance/fit_stage098_forced_carrier_card__yhat_q_cons.yaml
- fit_stage098_support_demand_exact (param=support_ratio): redteam_adversarial/provenance/fit_stage098_support_demand_exact__support_ratio.yaml
- fit_stage099_forced_by_yhatQcons_script (param=Gamma_5): redteam_adversarial/provenance/fit_stage099_forced_by_yhatqcons_script__gamma_5.yaml
- fit_stage100_factorization_forced_script (param=N_Q): redteam_adversarial/provenance/fit_stage100_factorization_forced_script__n_q.yaml

## Sources to open as needed
Each provenance bundle carries a `constraints[].constraint_kind` + `evidence_citation` and often a `provenance_findings` block stating the Phase-B synthesis agent's RATIONALE for the internal_consistency call (with cited notes/card lines). Treat that rationale as a CLAIM TO ATTACK, not ground truth — verify it against the actual files. Paper cards: research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes: research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (EM-projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage; read how the value is LABELED and find the line that DERIVES it. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (mis-classified posit/fit dressed as derived?) | one-clause reason citing BOTH the label line AND the derivation-step line (or its absence)`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific value + the overclaiming card/notes file+line + why no upstream step derives it.
3. Classification assessment (the crux): for each member, is `internal_consistency` CORRECT (a genuine derivation exists) or should it be free_choice / published_target (a posit / external back-solve)? Name any mis-classification stage+line, or state the classification is sound.
4. Proof you looked: which card+notes files you opened for each anchor stage, the exact label-word, and the exact derivation-step line you found (or confirmed absent).
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
