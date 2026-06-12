You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u01_aspect_ratio_37_20`, stages 121, 139, 146) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The L/a = 37/20 aspect-ratio geometry family carried 121 -> 139 -> 146 (re-fit drift check)
plus the 121 geometry-selection siblings L_W, L, mathfrak_r. SENSITIVE: "r_not_tunable" is the
single strongest overclaim-loaded name in the batch -- a value literally named "not tunable" must
be checked against its card/notes label: is mathfrak_r (and the 37/20 ratio) disclosed as a
posited reference freeze / geometric selection, or claimed forced/exact/non-tunable/derived?
37/20 is the same value batch-2 flagged at stage 073 (punch-list A5/A6). "reference_freeze",
"exact_branch", "literal", "primitive" naming -> verify each is honestly a posited/selected
geometric choice, not a derived consequence. Check L_W=L identification is an identification, not
a derivation. Orchestrator will spot-check stage_121/139/146 cards + notes regardless of verdict.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage121_aspect_ratio_37_20_reference_freeze (param=L_over_a): redteam_adversarial/provenance/fit_stage121_aspect_ratio_37_20_reference_freeze__l_over_a.yaml
- fit_stage121_lw_equals_l_identification (param=L_W): redteam_adversarial/provenance/fit_stage121_lw_equals_l_identification__l_w.yaml
- fit_stage121_r_geom_exact_branch (param=L): redteam_adversarial/provenance/fit_stage121_r_geom_exact_branch__l.yaml
- fit_stage121_r_not_tunable (param=mathfrak_r): redteam_adversarial/provenance/fit_stage121_r_not_tunable__mathfrak_r.yaml
- fit_stage139_aspect_ratio_literal (param=L_over_a): redteam_adversarial/provenance/fit_stage139_aspect_ratio_literal__l_over_a.yaml
- fit_stage146_aspect_ratio_37_20_primitive (param=L_over_a): redteam_adversarial/provenance/fit_stage146_aspect_ratio_37_20_primitive__l_over_a.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
