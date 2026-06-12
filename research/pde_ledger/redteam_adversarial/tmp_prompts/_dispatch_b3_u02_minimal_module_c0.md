You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u02_minimal_module_c0`, stages 088, 090) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
Minimal-module coefficients at 088 + downstream 090. c0/c_0 (3 bundles), c_contact (088+090),
Omega_Q pole frequency. Naming is overclaim-loaded: "c0_c1_derived_not_defined",
"natural_identification_fixes_exactly", "minimal_module_coefficients_paper_quoted",
"upstream_contradiction" -> the posit-dressed-as-derived test bites hardest here. c0=3/4 is an
ANSATZ_LEDGER branch-determinable value; classification settled free_choice, test = disclosure only.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage088_c0_c1_derived_not_defined (param=c0): redteam_adversarial/provenance/fit_stage088_c0_c1_derived_not_defined__c0.yaml
- fit_stage088_minimal_module_c0_c1 (param=c_0): redteam_adversarial/provenance/fit_stage088_minimal_module_c0_c1__c_0.yaml
- fit_stage088_minimal_module_coefficients_paper_quoted (param=c0): redteam_adversarial/provenance/fit_stage088_minimal_module_coefficients_paper_quoted__c0.yaml
- fit_stage088_natural_identification_fixes_exactly (param=c_contact): redteam_adversarial/provenance/fit_stage088_natural_identification_fixes_exactly__c_contact.yaml
- fit_stage088_pole_frequency_OmegaQ (param=Omega_Q): redteam_adversarial/provenance/fit_stage088_pole_frequency_omegaq__omega_q.yaml
- fit_stage090_minimal_module_upstream_contradiction (param=c_contact): redteam_adversarial/provenance/fit_stage090_minimal_module_upstream_contradiction__c_contact.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
