You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u12_tierB_crossfamily_sample`, stages 002, 026, 031, 040, 054, 077, 094, 115, 138, 139, 148, 187, 194, 222, 224, 228, 230, 251) in the INTERNAL-CONSISTENCY COMPLETENESS-CRITIC spot-check. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for internal_consistency values (THE INVERTED TEST — READ CAREFULLY)
These values are classified `constraint_kind: internal_consistency` — i.e. the ledger CLAIMS each is a DERIVED consequence (an algebraic identity / Solve / closure / extremum / limit / normalization) of upstream POSITED inputs, NOT itself an inserted choice. Unlike the free_choice sweep, the fatal flaw here is the OPPOSITE mis-classification: **a value that is genuinely a POSITED free_choice (an inserted number / a tuned knob) OR a BACK-SOLVE from an EXTERNAL published number, but is mis-classified `internal_consistency` AND labelled "derived"/"forced"/"exact"/"fixed by"/"a theorem"/"not free" in a paper CARD or per-stage NOTES line** — a posit-as-derived (or fit-as-derived) that ESCAPED the free_choice + published_target sweeps via mis-classification. These candidates were selected BECAUSE their slug carries an overclaim-flavored word (forced/fixed_exactly/backsolve/injection/not_free/_req/unique/...). The slug is the PROMPT to scrutinize, NOT evidence.
VERDICT RULE: **NO** = the value is a GENUINE derivation — open the notes/card and find the actual upstream step that PRODUCES it (a pole condition, a Solve, a closure relation, an extremum, a uniqueness criterion, a parity/symmetry forcing zero, an exact-algebra consequence of a NAMED posited input); the slug's "forced/exact" correctly names that consequence. **YES** = the value is actually inserted/tuned/back-solved (no upstream step produces it; it is chosen to hit a downstream or external target and then relabelled) — give the specific value + the card/notes file+line that calls it derived/forced/exact. Also flag (non-fatal) any stale upstream stage-anchor, value-divergence between members, or a slug that overstates vs the actual (honest) card/notes wording.

## NOTE: this is a TIER-B SAMPLE unit
These slugs carry WEAKER/structural words (exact/closure/canonical/derived/carried/fixed/determined). Across all 7 free_choice batches the repeatedly-confirmed finding was that such words attach to an algebraic IDENTITY / closure / law among symbols, never to a posited number. You are SPOT-CHECKING that the same holds for internal_consistency. Same inverted test; same YES/NO rule.

## Focus
Tier-B SAMPLE (the "samples" half of the mandate). Tier-B = 135 internal_ consistency bundles whose slugs carry WEAKER/structural overclaim words (exact/closure/canonical/derived/carried/fixed/determined). Across all 7 free_choice batches the decisive, repeatedly-confirmed finding was that these words attach to an algebraic IDENTITY / closure / law among symbols, never to a posited number. This unit SPOT-CHECKS that the SAME holds for internal_consistency by sampling one representative per family-word across the bands. For each: confirm the "exact/derived/canonical/closure" attaches to a genuine derivation step (identity/Solve/normalization/extremum), NOT to an inserted value relabeled. Specifically named seed-list Tier-B items to INCLUDE: 224 grouped_signature_exact (b_P0_slope), 228 k_compat_exact_ stiffness (k_compat), 230 stage228_rigid_slope_carry (s_minus_den). These are representative; the full Tier-B population (135) is logged as scanned-not-individually-deep-audited with this sampling rationale.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage002_overlap_factor_4pi_derived (param=(sample)): redteam_adversarial/provenance/fit_stage002_overlap_factor_4pi_derived__overlap_factor_4pi.yaml
- fit_stage026_kappa_exact_overlap (param=(sample)): redteam_adversarial/provenance/fit_stage026_kappa_exact_overlap__kappa.yaml
- fit_stage031_hf_identity_derived_not_assumed (param=(sample)): redteam_adversarial/provenance/fit_stage031_hf_identity_derived_not_assumed__dlambda_minus_dalpha.yaml
- fit_stage040_lambda_exact_closed_form (param=(sample)): redteam_adversarial/provenance/fit_stage040_lambda_exact_closed_form__lambda_minus.yaml
- fit_stage054_robin_branch_determined (param=(sample)): redteam_adversarial/provenance/fit_stage054_robin_branch_determined__a_k.yaml
- fit_stage077_If_exact_canonical (param=(sample)): redteam_adversarial/provenance/fit_stage077_if_exact_canonical__i_f.yaml
- fit_stage094_static_split_assigned_not_derived (param=(sample)): redteam_adversarial/provenance/fit_stage094_static_split_assigned_not_derived__c_geom.yaml
- fit_stage115_gq_exact_two_branch_law (param=(sample)): redteam_adversarial/provenance/fit_stage115_gq_exact_two_branch_law__g_q.yaml
- fit_stage138_rq_quarter_derived_not_assumed (param=(sample)): redteam_adversarial/provenance/fit_stage138_rq_quarter_derived_not_assumed__r_q.yaml
- fit_stage139_gains_land_exactly_on_canonical (param=(sample)): redteam_adversarial/provenance/fit_stage139_gains_land_exactly_on_canonical__m_q_comp_star.yaml
- fit_stage148_canonical_first_order_shifts_exact (param=(sample)): redteam_adversarial/provenance/fit_stage148_canonical_first_order_shifts_exact__delta_pi_lambda.yaml
- fit_stage187_exact_level_set_orbit_identity (param=(sample)): redteam_adversarial/provenance/fit_stage187_exact_level_set_orbit_identity__c_nt_star.yaml
- fit_stage194_chi_q_fixed_on_canonical_branch (param=(sample)): redteam_adversarial/provenance/fit_stage194_chi_q_fixed_on_canonical_branch__chi_q.yaml
- fit_stage222_exact_overlap_constant_kappa (param=(sample)): redteam_adversarial/provenance/fit_stage222_exact_overlap_constant_kappa__kappa.yaml
- fit_stage224_grouped_signature_exact (param=b_P0_slope): redteam_adversarial/provenance/fit_stage224_grouped_signature_exact__b_p0_slope.yaml
- fit_stage228_k_compat_exact_stiffness (param=k_compat): redteam_adversarial/provenance/fit_stage228_k_compat_exact_stiffness__k_compat.yaml
- fit_stage230_stage228_rigid_slope_carry (param=s_minus_den): redteam_adversarial/provenance/fit_stage230_stage228_rigid_slope_carry__s_minus_den.yaml
- fit_stage251_gamma_safe_exact_replacement (param=(sample)): redteam_adversarial/provenance/fit_stage251_gamma_safe_exact_replacement__gamma_safe.yaml

## Sources to open as needed
Each provenance bundle carries a `constraints[].constraint_kind` + `evidence_citation` and often a `provenance_findings` block stating the Phase-B synthesis agent's RATIONALE for the internal_consistency call (with cited notes/card lines). Treat that rationale as a CLAIM TO ATTACK, not ground truth — verify it against the actual files. Paper cards: research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes: research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (EM-projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage; read how the value is LABELED and find the line that DERIVES it. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (mis-classified posit/fit dressed as derived?) | one-clause reason citing BOTH the label line AND the derivation-step line (or its absence)`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific value + the overclaiming card/notes file+line + why no upstream step derives it.
3. Classification assessment (the crux): for each member, is `internal_consistency` CORRECT (a genuine derivation exists) or should it be free_choice / published_target (a posit / external back-solve)? Name any mis-classification stage+line, or state the classification is sound.
4. Proof you looked: which card+notes files you opened for each anchor stage, the exact label-word, and the exact derivation-step line you found (or confirmed absent).
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
