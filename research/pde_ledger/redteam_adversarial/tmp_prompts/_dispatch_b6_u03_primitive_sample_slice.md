You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u03_primitive_sample_slice`, stages 222, 225, 228) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The posited finite-throat ADMISSIBLE SAMPLE BRANCH -- a tuple of ILLUSTRATIVE numerical placements carried
222 -> 225 -> 228. Two compound bundles at 222 (numerical_primitive_slice: K,M,Omega_U,Omega_W,
lambda_B,lambda_R,lambda_U,lambda_W,varpi; primitive_sample_slice: a,c_s,K,lam_B,lam_R,lam_U,lam_W,M,
Omega_U,Omega_W,varpi -- the lam_* are the same couplings under a script alias) + the standalone lambda_B
+ the carried sample points OmU=Omega_U=1 (225) and M=1 (228 primitive_parameter_tuple). The tuple is
posited at notes 222 line 289-293 (a=c_s=1 at line 293) and stage_222.tex line 15 states "the displayed
slice is a numerical placement"; notes line 7/434 flag the whole slice as illustrative, not branch-derived.
Test = every sample-slice value is honestly an illustrative placement / carried sample point, and NO card
or notes line claims any of them derived/forced/uniquely-determined. HIGH PARAM COUNT -> orchestrator
confirms the stage_222.tex:15 placement disclosure + notes illustrative flag cover every param.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__k.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__m.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__omega_u.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__omega_w.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__lambda_b.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__lambda_r.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__lambda_u.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__lambda_w.yaml
- fit_stage222_numerical_primitive_slice (param=K, M, Omega_U, Omega_W, lambda_B, lambda_R, lambda_U, lambda_W, varpi): redteam_adversarial/provenance/fit_stage222_numerical_primitive_slice__varpi.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__a.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__c_s.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__k.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__lam_b.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__lam_r.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__lam_u.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__lam_w.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__m.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__omega_u.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__omega_w.yaml
- fit_stage222_primitive_sample_slice (param=a, c_s, K, lam_B, lam_R, lam_U, lam_W, M, Omega_U, Omega_W, varpi): redteam_adversarial/provenance/fit_stage222_primitive_sample_slice__varpi.yaml
- fit_stage_222_lambda_b (param=lambda_B): redteam_adversarial/provenance/fit_stage_222_lambda_b__lambda_b.yaml
- fit_stage225_stage223_sample_point (param=OmU): redteam_adversarial/provenance/fit_stage225_stage223_sample_point__omu.yaml
- fit_stage228_primitive_parameter_tuple (param=M): redteam_adversarial/provenance/fit_stage228_primitive_parameter_tuple__m.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
