You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage156_carried_constants_block__phase_c_adversarial.md

This is a FAMILY audit (parameter-value family `chain_aspect_37_20`) — cover EVERY member with its OWN sub-verdict.


## Focus
THE CENTRAL ADJUDICATION of the whole audit. The Family-1 canonical block (r_F1 and the carried constants at stages 156/158/161) bottoms out on ONE posited input L/a = 37/20 (origin stage 073/121). Is 37/20 a free_choice (posited toy-geometry aspect ratio: notes call it 'a reference-branch numerical freeze, not a new theorem'; no external number named, no back-solve) OR a published_target (external fit)? Project memory once grouped it with the external 54/5 class; in-band evidence says posit. ALSO: member fit_stage163_Rstar_quarter is a HIGH self_contradictory_origin (value 1/4 = R_q=R_*, NOT r_F1 — scanner mislabel). Give per-member verdicts AND a clear free_choice-vs-published_target determination with your reasoning.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage156_carried_constants_block (param=g_star): redteam_adversarial/provenance/fit_stage156_carried_constants_block__g_star.yaml
- fit_stage156_carried_constants_block (param=Pi_star): redteam_adversarial/provenance/fit_stage156_carried_constants_block__pi_star.yaml
- fit_stage156_carried_constants_block (param=r_F1): redteam_adversarial/provenance/fit_stage156_carried_constants_block__r_f1.yaml
- fit_stage156_carried_constants_block (param=Sigma0_star): redteam_adversarial/provenance/fit_stage156_carried_constants_block__sigma0_star.yaml
- fit_stage156_carried_constants_block (param=T_hat_star): redteam_adversarial/provenance/fit_stage156_carried_constants_block__t_hat_star.yaml
- fit_stage158_canonical_quartet_literals (param=r_F1): redteam_adversarial/provenance/fit_stage158_canonical_quartet_literals__r_f1.yaml
- fit_stage158_canonical_quartet_literals (param=S_can): redteam_adversarial/provenance/fit_stage158_canonical_quartet_literals__s_can.yaml
- fit_stage158_canonical_quartet_literals (param=Sigma0_can): redteam_adversarial/provenance/fit_stage158_canonical_quartet_literals__sigma0_can.yaml
- fit_stage158_canonical_quartet_literals (param=T_can): redteam_adversarial/provenance/fit_stage158_canonical_quartet_literals__t_can.yaml
- fit_stage161_rF1_family1_radius (param=r_F1): redteam_adversarial/provenance/fit_stage161_rf1_family1_radius__r_f1.yaml
- fit_stage163_Rstar_quarter (param=r_F1): redteam_adversarial/provenance/fit_stage163_rstar_quarter__r_f1.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open what you need to ground each finding.

## Benchmarks (use sourced entries, NEVER memory)
If this family claims an external match: GR mass-quadrupole 2G/5c^5 -> benchmarks.yaml family_ids fam_0187_p0_target / fam_0093_gamma_5 / fam_0174_n_q; CODATA proton/electron mass ratio 1836.15267343 -> fam_0457_m_s / fam_0458_m_s_ratio_1836. Read research/pde_ledger/redteam_adversarial/benchmarks.yaml. The rendered prompt already embeds the relevant entry if applicable.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id | YES/NO (fatal flaw?) | one-clause reason` (one row per member).
2. Family verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific external value / broken step / overclaiming card+line.
3. Disclosure assessment (layer-2 crux): honestly labeled target/benchmark/open, or overclaimed derived/forced/exact? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which benchmark entries (if any external-match claim) checked digit-by-digit, and which fit-vs-derive tests you ran.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts — read and reason.
