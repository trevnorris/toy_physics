You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u01_early_algebraic_forcing_A`, stages 006, 019, 024, 027, 029, 032, 035, 036) in the INTERNAL-CONSISTENCY COMPLETENESS-CRITIC spot-check. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for internal_consistency values (THE INVERTED TEST — READ CAREFULLY)
These values are classified `constraint_kind: internal_consistency` — i.e. the ledger CLAIMS each is a DERIVED consequence (an algebraic identity / Solve / closure / extremum / limit / normalization) of upstream POSITED inputs, NOT itself an inserted choice. Unlike the free_choice sweep, the fatal flaw here is the OPPOSITE mis-classification: **a value that is genuinely a POSITED free_choice (an inserted number / a tuned knob) OR a BACK-SOLVE from an EXTERNAL published number, but is mis-classified `internal_consistency` AND labelled "derived"/"forced"/"exact"/"fixed by"/"a theorem"/"not free" in a paper CARD or per-stage NOTES line** — a posit-as-derived (or fit-as-derived) that ESCAPED the free_choice + published_target sweeps via mis-classification. These candidates were selected BECAUSE their slug carries an overclaim-flavored word (forced/fixed_exactly/backsolve/injection/not_free/_req/unique/...). The slug is the PROMPT to scrutinize, NOT evidence.
VERDICT RULE: **NO** = the value is a GENUINE derivation — open the notes/card and find the actual upstream step that PRODUCES it (a pole condition, a Solve, a closure relation, an extremum, a uniqueness criterion, a parity/symmetry forcing zero, an exact-algebra consequence of a NAMED posited input); the slug's "forced/exact" correctly names that consequence. **YES** = the value is actually inserted/tuned/back-solved (no upstream step produces it; it is chosen to hit a downstream or external target and then relabelled) — give the specific value + the card/notes file+line that calls it derived/forced/exact. Also flag (non-fatal) any stale upstream stage-anchor, value-divergence between members, or a slug that overstates vs the actual (honest) card/notes wording.

## Focus
Earliest band. Each slug claims an algebraic FORCING of a value by a pole condition / isotropy / eigenvector / locus. Verify each value is a genuine upstream algebraic consequence (the notes show a Solve / pole / locus step producing it), NOT an inserted number relabeled "forced"/"fixed_by"/"unique". ampere_sign (006) is an EM-projected stage -> notes live at repo-root notes/em_projected/ (sign fixed by an audit/convention, verify it is a derived sign not a tuned value). mhat_0 "completely_fixed" (032) and xi_req/g_B_req "unique/req" (035/036): confirm these are determined by the stated constraint, not back-solved.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage006_projected_vector_signs_fixed_by_audit (param=ampere_sign): redteam_adversarial/provenance/fit_stage006_projected_vector_signs_fixed_by_audit__ampere_sign.yaml
- fit_stage019_K_Sigma_fixed_by_one_pole (param=K_Sigma): redteam_adversarial/provenance/fit_stage019_k_sigma_fixed_by_one_pole__k_sigma.yaml
- fit_stage024_grouped_bundle_forced_isotropic (param=grouped_lane_coefficients_20_21_22): redteam_adversarial/provenance/fit_stage024_grouped_bundle_forced_isotropic__grouped_lane_coefficients_20_21_22.yaml
- fit_stage027_max_coupling_branch_not_free (param=kappa_theta): redteam_adversarial/provenance/fit_stage027_max_coupling_branch_not_free__kappa_theta.yaml
- fit_stage029_deltaK_algebraically_forced (param=DeltaK_ax): redteam_adversarial/provenance/fit_stage029_deltak_algebraically_forced__deltak_ax.yaml
- fit_stage032_mhat_minus_completely_fixed (param=mhat_0): redteam_adversarial/provenance/fit_stage032_mhat_minus_completely_fixed__mhat_0.yaml
- fit_stage035_xi_req_unique_locus (param=xi_req): redteam_adversarial/provenance/fit_stage035_xi_req_unique_locus__xi_req.yaml
- fit_stage036_no_free_hidden_eigenvector_param (param=g_B_req): redteam_adversarial/provenance/fit_stage036_no_free_hidden_eigenvector_param__g_b_req.yaml

## Sources to open as needed
Each provenance bundle carries a `constraints[].constraint_kind` + `evidence_citation` and often a `provenance_findings` block stating the Phase-B synthesis agent's RATIONALE for the internal_consistency call (with cited notes/card lines). Treat that rationale as a CLAIM TO ATTACK, not ground truth — verify it against the actual files. Paper cards: research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes: research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (EM-projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage; read how the value is LABELED and find the line that DERIVES it. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (mis-classified posit/fit dressed as derived?) | one-clause reason citing BOTH the label line AND the derivation-step line (or its absence)`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific value + the overclaiming card/notes file+line + why no upstream step derives it.
3. Classification assessment (the crux): for each member, is `internal_consistency` CORRECT (a genuine derivation exists) or should it be free_choice / published_target (a posit / external back-solve)? Name any mis-classification stage+line, or state the classification is sound.
4. Proof you looked: which card+notes files you opened for each anchor stage, the exact label-word, and the exact derivation-step line you found (or confirmed absent).
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
