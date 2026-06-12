You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u05_session4_collapse`, stages 250, 251) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
Session-IV collapse rate/time constants and one Stage-250 free_choice readback. chi_rm (250 chi_rm) is the
ONLY free_choice candidate left at Stage 250 (the benchmark_masses_and_window / goldilocks_inputs candidates
-- chi_peak,chi_raw,E_max,g_UV,mu_eta,t_transit_*,m_s -- were already audited in batch 1 as the
published_target/external cluster); chi_rm is a reported readback / posited scan quantity. gamma_crit (251
session4_rate_constants + session4_time_inputs, same value) is a posited critical collapse rate constant.
t_collapse_0 and t_cross (251 session4_time_inputs) are Session-IV time-input sample points (NOTE: the
separate candidate fit_stage_250_t_cross is "derived by construction" internal_consistency and is NOT in
this unit; here t_cross is the 251 free_choice time-input). cref (251 cref) is a posited reference speed/
coefficient. The s_c entry in session4_time_inputs is internal_consistency -> EXCLUDED. Test = each is
honestly a reported readback / posited rate or time constant / reference coefficient, NONE dressed as
derived/forced. Check stage_250/251 cards + notes for the label-words.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage_250_chi_rm (param=chi_rm): redteam_adversarial/provenance/fit_stage_250_chi_rm__chi_rm.yaml
- fit_stage251_session4_rate_constants (param=gamma_crit): redteam_adversarial/provenance/fit_stage251_session4_rate_constants__gamma_crit.yaml
- fit_stage251_session4_time_inputs (param=gamma_crit, t_collapse_0, t_cross): redteam_adversarial/provenance/fit_stage251_session4_time_inputs__gamma_crit.yaml
- fit_stage251_session4_time_inputs (param=gamma_crit, t_collapse_0, t_cross): redteam_adversarial/provenance/fit_stage251_session4_time_inputs__t_collapse_0.yaml
- fit_stage251_session4_time_inputs (param=gamma_crit, t_collapse_0, t_cross): redteam_adversarial/provenance/fit_stage251_session4_time_inputs__t_cross.yaml
- fit_stage_251_cref (param=cref): redteam_adversarial/provenance/fit_stage_251_cref__cref.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
