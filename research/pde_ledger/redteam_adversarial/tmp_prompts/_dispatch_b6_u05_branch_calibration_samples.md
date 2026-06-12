You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u05_branch_calibration_samples`, stages 224, 233, 234, 236, 237) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
Branch-calibration assumptions, held-fixed axes, prescribed free targets, lane-weight conventions, a
script-only sample point, and one scanner artifact. Members: lambda_20=1 (224
weak_axisymmetric_lane_signature) = one entry of the imported grouped (l=2) lane signature
(lambda_20,lambda_21,lambda_22)=(1,1/2,-1), a posited convention. Delta_rm (224 delta_rm) = a scanner
truncation of Delta_norm; value 0 is the CONDITIONAL calibration assumption Delta_norm=0 ("if the real
branch already hits the universal quadrupole normalization"). `hat` (224) is NOT a real parameter -- it is
a LaTeX-command fragment scanner-extracted from \hat m_0 (mhat_0 is an undetermined branch quantity with
only a lower bound at 224); flag as a scanner artifact (benign). Delta_norm (233
calibrated_branch_delta_norm) value 0 is again a CALIBRATION / branch selection ("On a calibrated branch
with Delta_norm=0..."), explicitly NOT derived. Xi_1 (234 canonical_least_deformation_family) is NOT a
fixed number -- it is a PRESCRIBED first-order defect target ("Xi_1 = xi") that a "canonical
least-deformation family realizes any chosen Xi_1"; "canonical" attaches to the family construction, not
to forcing Xi_1. T_U (236) is one axis of a dependent triple, explicitly HELD FIXED ("at fixed T_U"), not
derived. Rratio_base=7/6 (237 branch_sample_point) is a SCRIPT-ONLY arbitrary substitution sample used to
numerically verify symbolic identities. Test = every member is honestly a convention / conditional
calibration / held-fixed axis / prescribed-free target / script sample / scanner artifact, and NONE is
dressed as derived/forced. "canonical_least_deformation" must label the family, not force Xi_1. Check
notes224/233/234/236 + script237 for the label-words.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage224_weak_axisymmetric_lane_signature (param=lambda_20): redteam_adversarial/provenance/fit_stage224_weak_axisymmetric_lane_signature__lambda_20.yaml
- fit_stage_224_delta_rm (param=Delta_rm): redteam_adversarial/provenance/fit_stage_224_delta_rm__delta_rm.yaml
- fit_stage_224_hat (param=hat): redteam_adversarial/provenance/fit_stage_224_hat__hat.yaml
- fit_stage233_calibrated_branch_delta_norm (param=Delta_norm): redteam_adversarial/provenance/fit_stage233_calibrated_branch_delta_norm__delta_norm.yaml
- fit_stage234_canonical_least_deformation_family (param=Xi_1): redteam_adversarial/provenance/fit_stage234_canonical_least_deformation_family__xi_1.yaml
- fit_stage_236_t_u (param=T_U): redteam_adversarial/provenance/fit_stage_236_t_u__t_u.yaml
- fit_stage237_branch_sample_point (param=Rratio_base): redteam_adversarial/provenance/fit_stage237_branch_sample_point__rratio_base.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
