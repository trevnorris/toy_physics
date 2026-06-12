You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u11_aspect_37_20_origin`, stages 073) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
L/a = 37/20 ORIGIN bundles. NOTE: classification ALREADY adjudicated in batch 1 = free_choice
(user-accepted 2026-06-11). Focus here is PURELY the card/notes overclaim test (stage_192 failure mode):
does any stage-073 card or notes line call 37/20 derived/forced/exact/a-theorem rather than a posited
reference-branch freeze? Classification is settled; report only disclosure adequacy.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage073_Lambda_ell_fixed_37 (param=L_over_a): redteam_adversarial/provenance/fit_stage073_lambda_ell_fixed_37__l_over_a.yaml
- fit_stage073_aspect_ratio_frozen (param=L_over_a): redteam_adversarial/provenance/fit_stage073_aspect_ratio_frozen__l_over_a.yaml
- fit_stage073_aspect_ratio_L_over_a (param=L_over_a): redteam_adversarial/provenance/fit_stage073_aspect_ratio_l_over_a__l_over_a.yaml
- fit_stage073_paper_fixes_geometry_ratios (param=L_over_a): redteam_adversarial/provenance/fit_stage073_paper_fixes_geometry_ratios__l_over_a.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
