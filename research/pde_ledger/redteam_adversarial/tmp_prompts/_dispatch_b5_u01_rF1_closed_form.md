You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u01_rF1_closed_form`, stages 162, 163, 165, 168) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The family-1 radius r_F1 closure carried 162 -> 163 -> 165 -> 168 (re-fit drift check). This is a direct
CONTINUATION of batch-4 u03 (r_F1 at 122/123/124), now at its "closed-form / canonical-point" anchors.
Naming is overclaim-adjacent: "r_f1_family1_value", "rF1_family1_radius", "family1_canonical_point_block",
"rstar_closed_form_4107", "rexact_closed_form_4107". Test = is r_F1 disclosed as a geometry-selected /
carried family-1 branch choice (the sqrt(4107-100 pi^2)/(10 pi) closed form evaluated at L/a=37/20), or is
it claimed forced/fixed/exact/canonical-as-unique by the geometry? "closed_form"/"rexact"/"4107" must label
the FORMULA/algebraic closure, not assert the VALUE is derived. "canonical_point" must mean a chosen
canonical reference point, not a uniqueness theorem. Check all anchors carry a NUMERICALLY IDENTICAL r_F1
(1.77799353547498 = sqrt(4107-100 pi^2)/(10 pi); a fit would drift). Note 163 bundle carries r_F1 only here
(its g_star/s_can/sigma0_can params are internal_consistency, out of scope). Orchestrator will spot-check
stage_162/163/165/168 cards + notes regardless of verdict.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage162_r_f1_family1_value (param=r_F1): redteam_adversarial/provenance/fit_stage162_r_f1_family1_value__r_f1.yaml
- fit_stage162_rF1_family1_radius (param=r_F1): redteam_adversarial/provenance/fit_stage162_rf1_family1_radius__r_f1.yaml
- fit_stage163_family1_canonical_point_block (param=r_F1): redteam_adversarial/provenance/fit_stage163_family1_canonical_point_block__r_f1.yaml
- fit_stage165_rstar_closed_form_4107 (param=r_F1): redteam_adversarial/provenance/fit_stage165_rstar_closed_form_4107__r_f1.yaml
- fit_stage168_rexact_closed_form_4107 (param=r_F1): redteam_adversarial/provenance/fit_stage168_rexact_closed_form_4107__r_f1.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
