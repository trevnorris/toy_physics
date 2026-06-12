You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u04_session2_3_helicity`, stages 248, 249) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
Session-II/III launch scales and helicity/peak readbacks. r_contact=0.18 (248
chosen_contact_scale_benchmark) is a DECLARED "chosen contact radius" fixed by the Session-II benchmark
specialization -- a posited benchmark scale, explicitly not derived and not external. E_sub=2.5 (248
launch_parameters + readback, same value reused) is a posited sub-barrier launch energy ("not a new theorem
by itself"). The Stage-249 helicity quantities -- H_sub_minus_final (helicity_ratio_matched), H_int_aligned
(session2_helicity_readbacks), hint_aligned (session3_peak_readback) -- are UN-DERIVED Session-II/III run
OUTPUTS read back as diagnostics; the Stage-249 notes are explicit this is "diagnostic rather than a spin
theorem" / "does NOT prove intrinsic spin-1/2 ... not yet a theorem of microscopic spin support". C_sigma
(249 c_sigma) and R_int (249 r_int) are likewise reported run readbacks (R_int computed from un-derived
Session-II outputs R_pk/R_int). Vocabulary "matched"/"aligned"/"chosen benchmark" describes a
FIT/READBACK, not a derivation. Test = every value is honestly a chosen benchmark scale / posited launch
energy / un-derived run readback, and NO card or notes line upgrades any of them to a derived/forced
result or a "theorem". Check stage_248/249 cards + notes for the label-words.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage248_chosen_contact_scale_benchmark (param=r_contact): redteam_adversarial/provenance/fit_stage248_chosen_contact_scale_benchmark__r_contact.yaml
- fit_stage248_session2_launch_parameters (param=E_sub): redteam_adversarial/provenance/fit_stage248_session2_launch_parameters__e_sub.yaml
- fit_stage248_session2_readback (param=E_sub): redteam_adversarial/provenance/fit_stage248_session2_readback__e_sub.yaml
- fit_stage249_helicity_ratio_matched (param=H_sub_minus_final): redteam_adversarial/provenance/fit_stage249_helicity_ratio_matched__h_sub_minus_final.yaml
- fit_stage249_session2_helicity_readbacks (param=H_int_aligned): redteam_adversarial/provenance/fit_stage249_session2_helicity_readbacks__h_int_aligned.yaml
- fit_stage249_session3_peak_readback (param=hint_aligned): redteam_adversarial/provenance/fit_stage249_session3_peak_readback__hint_aligned.yaml
- fit_stage_249_c_sigma (param=C_sigma): redteam_adversarial/provenance/fit_stage_249_c_sigma__c_sigma.yaml
- fit_stage_249_r_int (param=R_int): redteam_adversarial/provenance/fit_stage_249_r_int__r_int.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
