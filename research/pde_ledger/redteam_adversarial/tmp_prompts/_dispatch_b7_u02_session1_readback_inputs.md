You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u02_session1_readback_inputs`, stages 245, 246) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
Session-I readback sample points, probe stiffnesses, and the logarithmic dressing-chart reference. ALL are
recorded interactive-run values or posited sample/probe inputs, fed back as consistency checks into the
EXACT Stage-245 finite (U,V) compiler (the compiler is internal_consistency-derived; the readback inputs
are free_choice). epsilon_eta(r_eval)~0.2893 (245 eps_eta_session1_match) is eps_eta=eps_eta_ref e^V
evaluated at a probe/sample point, not first-principles-derived. U_obs/V_obs (245 readback / readback_inputs)
are the recorded matched dressing readback (U_obs~0.1431 "relaxed-rigid-mouth point", script Float), not
ledger-derived. a_U,a_V (>0 transfer-shape / dressing-leg stiffnesses), eps_eta_ref (chart reference point
of V=ln(eps_eta/eps_eta_ref)), g_UV (posited closure parameter, no derived numeric value) -- all posited,
low-severity (no derived numeric value to classify). a_0,b_0,r_eval,r_sigma (246 session1_inputs) and a0
(246 readback_point) are representative Session-I mouth/source sample values reproducing the recorded
diagnostics. Test = every value is honestly an interactive readback / posited probe sample / chart
reference, and NO card or notes line claims any of them derived/forced/uniquely-determined. WATCH: do not
mis-read the script comments "F5 pin independently derived benchmark quantities" / "F4 lowering theorem" --
those are about OTHER (lambda_L / potential-lowering) quantities, not these readback inputs. Check
notes245/246 for the label-words.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage245_eps_eta_session1_match (param=epsilon_eta): redteam_adversarial/provenance/fit_stage245_eps_eta_session1_match__epsilon_eta.yaml
- fit_stage245_session1_readback (param=U_obs): redteam_adversarial/provenance/fit_stage245_session1_readback__u_obs.yaml
- fit_stage245_session1_readback_inputs (param=U_obs, V_obs, a_U, a_V, eps_eta_ref, g_UV): redteam_adversarial/provenance/fit_stage245_session1_readback_inputs__u_obs.yaml
- fit_stage245_session1_readback_inputs (param=U_obs, V_obs, a_U, a_V, eps_eta_ref, g_UV): redteam_adversarial/provenance/fit_stage245_session1_readback_inputs__v_obs.yaml
- fit_stage245_session1_readback_inputs (param=U_obs, V_obs, a_U, a_V, eps_eta_ref, g_UV): redteam_adversarial/provenance/fit_stage245_session1_readback_inputs__a_u.yaml
- fit_stage245_session1_readback_inputs (param=U_obs, V_obs, a_U, a_V, eps_eta_ref, g_UV): redteam_adversarial/provenance/fit_stage245_session1_readback_inputs__a_v.yaml
- fit_stage245_session1_readback_inputs (param=U_obs, V_obs, a_U, a_V, eps_eta_ref, g_UV): redteam_adversarial/provenance/fit_stage245_session1_readback_inputs__eps_eta_ref.yaml
- fit_stage245_session1_readback_inputs (param=U_obs, V_obs, a_U, a_V, eps_eta_ref, g_UV): redteam_adversarial/provenance/fit_stage245_session1_readback_inputs__g_uv.yaml
- fit_stage246_session1_inputs (param=a_0, b_0, r_eval, r_sigma): redteam_adversarial/provenance/fit_stage246_session1_inputs__a_0.yaml
- fit_stage246_session1_inputs (param=a_0, b_0, r_eval, r_sigma): redteam_adversarial/provenance/fit_stage246_session1_inputs__b_0.yaml
- fit_stage246_session1_inputs (param=a_0, b_0, r_eval, r_sigma): redteam_adversarial/provenance/fit_stage246_session1_inputs__r_eval.yaml
- fit_stage246_session1_inputs (param=a_0, b_0, r_eval, r_sigma): redteam_adversarial/provenance/fit_stage246_session1_inputs__r_sigma.yaml
- fit_stage246_session1_readback_point (param=a0): redteam_adversarial/provenance/fit_stage246_session1_readback_point__a0.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
