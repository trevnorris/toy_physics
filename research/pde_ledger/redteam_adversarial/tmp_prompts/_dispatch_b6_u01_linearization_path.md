You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u01_linearization_path`, stages 200, 202, 203) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## Focus
Early-band symbolic / script-only linearization & path coordinates plus one declared target-orbit value.
Members: eps_beta (200) is a SYMBOLIC linearization (perturbation) direction beta=1+eps*eps_beta, not a
fitted number -> verify it is a symbol, never a posited-as-derived value. C_tr_target (202,
"split_u_ratio_uniquely_fixed") is a DECLARED target-orbit specification ("Fix the target orbit by the
exact target values"); "uniquely_fixed" is overclaim-ADJACENT -> verify "uniquely_fixed" means the orbit
is fixed BY the chosen target (a parameterization choice), NOT that C_tr_target itself is derived/forced.
beta_path=2^(2*tau-1) (203) and KW_t=KW0*exp(t*kapW) (203) are SCRIPT-ONLY illustrative one-parameter path
coordinates driving the Stage-203 graph-crossing / intermediate-value demonstration -> verify they are
explicitly path/illustrative, not fitted/derived. Test = none of the four is dressed as derived; all are
symbolic / declared-target / script-illustrative. Check notes200/202/203 for the label-words.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage_200_eps_beta (param=eps_beta): redteam_adversarial/provenance/fit_stage_200_eps_beta__eps_beta.yaml
- fit_stage202_split_u_ratio_uniquely_fixed (param=C_tr_target): redteam_adversarial/provenance/fit_stage202_split_u_ratio_uniquely_fixed__c_tr_target.yaml
- fit_stage_203_beta_path (param=beta_path): redteam_adversarial/provenance/fit_stage_203_beta_path__beta_path.yaml
- fit_stage_203_kw_t (param=KW_t): redteam_adversarial/provenance/fit_stage_203_kw_t__kw_t.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
