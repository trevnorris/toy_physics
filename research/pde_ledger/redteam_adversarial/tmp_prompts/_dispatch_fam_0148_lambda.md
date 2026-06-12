You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage038_lambda_external_target_coefficient__phase_c_adversarial.md

This is a FAMILY audit (parameter-value family `fam_0148_lambda`) — cover EVERY member with its OWN sub-verdict.


## Focus
Lambda (stage 038) is the radiative demand scale that PACKAGES the external 2.5PN GR quadrupole target (HIGH published_target). Disclosed honestly?

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage038_lambda_external_target_coefficient (param=Lambda): redteam_adversarial/provenance/fit_stage038_lambda_external_target_coefficient__lambda.yaml

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
