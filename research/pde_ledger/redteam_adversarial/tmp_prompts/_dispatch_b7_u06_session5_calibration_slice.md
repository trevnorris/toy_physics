You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u06_session5_calibration_slice`, stages 252, 253) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The Stage-252 three-to-one split and the Stage-253 PHYSICAL-CALIBRATION SLICE -- the headline benchmark
stage (the 119.2336... figure). This is the campaign's final and most calibration-heavy unit; vocabulary
traps abound ("K_turn_force_match", "calibration_slice_recovery", "calibration"). f_lat=3/4 (252
three_to_one_split; also 253 nominal_constants, stage252_slice_inputs, _253_f_lat -- duplicate scanner
slugs for the SAME literal) is OPENLY "a posited (by-hand) split, not derived ... a calibration consistency
check, not a theorem that the realized branch must use the speed-independent split" (Phase-B). mu_eta (252
_252_mu_eta; also 253 nominal / session5_recovery / stage252_slice_inputs -- 4x scanner duplication) is a
posited material/effective coupling for the calibration slice. K_turn (253 k_turn_force_match,
nominal_constants, _253_k_turn -- triple scanner duplication; "force_match" in the slug) is force-MATCHED
to reproduce the benchmark turning point -- a CALIBRATION FIT (free_choice), fatal ONLY if a card/notes line
claims K_turn is derived/forced-by-first-principles rather than calibration-matched. gamma_lattice_legacy /
gamma_lattice_red (253 nominal / gamma_lattice_red / session5_recovery) are the legacy vs reduced lattice
damping calibration constants. Upsilon_lat (253 upsilon_lat_calibration + _calibration_factor, 2x) is the
lattice calibration factor. K_corr (253 _253_k_corr) is a posited correction coefficient. The two
..._sympy_audit / ..._mathematica_audit MEGA-SLUGS at 252/253 (parameter_name = the whole audit-script
filename) are pure scanner artifacts -> not in this unit; LOG only. Test = every calibration-slice constant
is honestly a posited by-hand split / calibration fit / nominal slice constant / calibration factor, and
NO card or notes line dresses any of them as derived/forced/uniquely-determined/a-theorem. "force_match"
and "calibration" must denote a FIT, not a derivation. HIGH duplication -> note the triple/quad scanner
dups (dedup item, not a content fault). Orchestrator reads stage_252/253 cards + notes.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage252_three_to_one_split (param=f_lat): redteam_adversarial/provenance/fit_stage252_three_to_one_split__f_lat.yaml
- fit_stage_252_mu_eta (param=mu_eta): redteam_adversarial/provenance/fit_stage_252_mu_eta__mu_eta.yaml
- fit_stage253_k_turn_force_match (param=K_turn): redteam_adversarial/provenance/fit_stage253_k_turn_force_match__k_turn.yaml
- fit_stage253_calibration_slice_nominal_constants (param=f_lat, gamma_lattice_legacy, K_turn, mu_eta): redteam_adversarial/provenance/fit_stage253_calibration_slice_nominal_constants__f_lat.yaml
- fit_stage253_calibration_slice_nominal_constants (param=f_lat, gamma_lattice_legacy, K_turn, mu_eta): redteam_adversarial/provenance/fit_stage253_calibration_slice_nominal_constants__gamma_lattice_legacy.yaml
- fit_stage253_calibration_slice_nominal_constants (param=f_lat, gamma_lattice_legacy, K_turn, mu_eta): redteam_adversarial/provenance/fit_stage253_calibration_slice_nominal_constants__k_turn.yaml
- fit_stage253_calibration_slice_nominal_constants (param=f_lat, gamma_lattice_legacy, K_turn, mu_eta): redteam_adversarial/provenance/fit_stage253_calibration_slice_nominal_constants__mu_eta.yaml
- fit_stage253_gamma_lattice_red (param=gamma_lattice_red): redteam_adversarial/provenance/fit_stage253_gamma_lattice_red__gamma_lattice_red.yaml
- fit_stage253_session5_calibration_slice_recovery (param=gamma_lattice_red, mu_eta): redteam_adversarial/provenance/fit_stage253_session5_calibration_slice_recovery__gamma_lattice_red.yaml
- fit_stage253_session5_calibration_slice_recovery (param=gamma_lattice_red, mu_eta): redteam_adversarial/provenance/fit_stage253_session5_calibration_slice_recovery__mu_eta.yaml
- fit_stage253_stage252_slice_inputs (param=f_lat, mu_eta): redteam_adversarial/provenance/fit_stage253_stage252_slice_inputs__f_lat.yaml
- fit_stage253_stage252_slice_inputs (param=f_lat, mu_eta): redteam_adversarial/provenance/fit_stage253_stage252_slice_inputs__mu_eta.yaml
- fit_stage253_upsilon_lat_calibration (param=Upsilon_lat): redteam_adversarial/provenance/fit_stage253_upsilon_lat_calibration__upsilon_lat.yaml
- fit_stage253_upsilon_lat_calibration_factor (param=Upsilon_lat): redteam_adversarial/provenance/fit_stage253_upsilon_lat_calibration_factor__upsilon_lat.yaml
- fit_stage_253_f_lat (param=f_lat): redteam_adversarial/provenance/fit_stage_253_f_lat__f_lat.yaml
- fit_stage_253_k_corr (param=K_corr): redteam_adversarial/provenance/fit_stage_253_k_corr__k_corr.yaml
- fit_stage_253_k_turn (param=K_turn): redteam_adversarial/provenance/fit_stage_253_k_turn__k_turn.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
