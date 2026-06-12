You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u01_precursor_loading_ratio`, stages 240, 243) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The Stage-240 minimal-isotropic-precursor loading module and the Stage-243 lift profile amplitude.
c0_min=3/4 (240 minimal_isotropic_precursor; with c_1=1/4) is THE LOAD-BEARING POSITED INPUT: the entire
downstream "Exact selected loading ratio" chain rho_alpha=1/c0=4/3, zeta_req=c_1/c0=1/3, Pi_tr=4 C_mix/3 is
DERIVED from it by the inverse contact-plus-pole compiler. The Stage-240 card line 15 OUTPUT field labels
that chain "Exact selected loading ratio". CRUX: "exact" here is an EXACT ALGEBRAIC CONSEQUENCE of the
posited c0_min module (honest-conditional, the gamma_0 batch-5 shape) -- it is FATAL only if the card or
notes dress c0_min ITSELF as derived/forced/non-tunable/uniquely-determined. The card ALSO hedges
explicitly ("should not be read as a second independent proof of actual nonlinear branch realization";
Checks: "Check whether the row is algebraic, numerical placement, or open branch-realization"). Verify the
notes label c0_min as a posited "minimal isotropic module" / "canonical module value" (canonical=module
construction, NOT a forced number). If c0_min itself is honestly posited, verdict NO (or a gamma_0-style
non-fatal PARTIAL if you judge the bare "Exact selected loading ratio" wording overstates depth by lacking
an inline "conditional on the posited c0 module" qualifier -- report that as PARTIAL, not YES). Note c0_min
is an ANSATZ_LEDGER section-5 borderline (c0/c0_min dedup) already at USER gate for its layer-3 disposition;
do NOT decide the dedup, only the framing. E_w (243 e_w) is the posited weak-lift profile amplitude used to
evaluate the exact work integral W_w=sqrt(2 pi)/8 ell_w j0 E0; per Phase-B it "carries no inserted value to
classify as a fit ... the profile is posited for the lift, not derived" -> verify honest-posited. Read
stage_240 + stage_243 cards and notes for the exact label-words.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage240_minimal_isotropic_precursor (param=c0_min): redteam_adversarial/provenance/fit_stage240_minimal_isotropic_precursor__c0_min.yaml
- fit_stage_243_e_w (param=E_w): redteam_adversarial/provenance/fit_stage_243_e_w__e_w.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
