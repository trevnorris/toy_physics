You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u03_session1_couplings`, stages 247) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
Stage-247 Session-I softening/closure couplings. lambda_W (247 lambda_w_natural_specialization) -- the slug
says "natural_specialization"; VERIFY "natural" attaches to a specialization CHOICE of the carried lambda_W
coupling (a posited convention / natural value for the benchmark), NOT a claim that lambda_W is
derived/forced. a0 (247 session1_softening_point) is the strongest-softening sample evaluation point (do not
mis-attribute the nearby "F5 pin"/"F4 lowering theorem" script comments to a0; those concern lambda_L).
alpha_2 (247 alpha_2) is a posited closure/shape coefficient. NOTE the lambda_L back-solve candidates at 247
(benchmark_pinned_to_paper_figures, lambda_l_backsolve, lambda_l_fixed_by_recorded_point, lambda_l_closure)
are NOT in this unit -- they are published_target / internal back-solves to the paper figures (out of scope;
neighbors), and the Stage-247 card is itself cautious about them ("useful as an audit point, not as a
general theorem"). Test = lambda_W / a0 / alpha_2 are each honestly a posited specialization / sample point /
closure coefficient, NONE dressed as derived/forced. Check stage_247 card + notes for the label-words.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage247_lambda_w_natural_specialization (param=lambda_W): redteam_adversarial/provenance/fit_stage247_lambda_w_natural_specialization__lambda_w.yaml
- fit_stage247_session1_softening_point (param=a0): redteam_adversarial/provenance/fit_stage247_session1_softening_point__a0.yaml
- fit_stage_247_alpha_2 (param=alpha_2): redteam_adversarial/provenance/fit_stage_247_alpha_2__alpha_2.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
