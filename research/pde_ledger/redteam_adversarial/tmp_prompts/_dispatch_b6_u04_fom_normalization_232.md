You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u04_fom_normalization_232`, stages 232) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The Stage-232 figure-of-merit / normalization literals. Members: L_over_a (value 20) appears in two
DUPLICATE scanner slugs (family1_geometry_normalization + lambda_ell_geometry_literals) -- it is the
inverse thin-wall lock a/ell=20 inside Lambda_ell=20*sqrt(2)*pi/x01, originating at the Stage-073 posit
epsilon_r=1/20; it is NOT the batch-1 aspect ratio L/a=37/20 (stale_provenance_anchor mis-name). DO NOT
confuse the two; DO NOT re-litigate the 37/20 adjudication. Xi_prefactor_100 = the literal prefactor 100
in the FoM Xi=100*Theta_w*Lambda_ell^2 (no in-stage derivation; structure traces to Stage 066). N_Q=1 =
the canonical outgoing-mode normalization CONVENTION (touches the batch-1 N_Q quadrupole family but here
is plainly =1, a convention). lambda_mu=1 = a benchmark evaluation point ("for the benchmark", notes 232
line 149). Test = each is a posited convention / FoM prefactor / normalization choice / geometric literal
traceable to an upstream posit, and NONE is dressed as derived/forced. Orchestrator reads stage_232 card
+ notes; confirm L_over_a=20 is the thin-wall lock (not 37/20) and N_Q=1 is labeled a normalization
condition. NOTE the 232 notes file is the "known_5pn_data_injection" notes -- the 5PN reference is the
INTERNAL EM-projected chain, not an external publication (memory audit correction).

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage232_family1_geometry_normalization (param=L_over_a): redteam_adversarial/provenance/fit_stage232_family1_geometry_normalization__l_over_a.yaml
- fit_stage232_lambda_ell_geometry_literals (param=L_over_a): redteam_adversarial/provenance/fit_stage232_lambda_ell_geometry_literals__l_over_a.yaml
- fit_stage232_fom_prefactor_100 (param=Xi_prefactor_100): redteam_adversarial/provenance/fit_stage232_fom_prefactor_100__xi_prefactor_100.yaml
- fit_stage232_nq_canonical_normalization (param=N_Q): redteam_adversarial/provenance/fit_stage232_nq_canonical_normalization__n_q.yaml
- fit_stage_232_lambda_mu_1 (param=lambda_mu_1): redteam_adversarial/provenance/fit_stage_232_lambda_mu_1__lambda_mu_1.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
