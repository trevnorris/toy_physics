You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u06_P2_chain`, stages 191, 192, 193, 194, 195) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
MAXIMALLY SENSITIVE. The P_2 multipole-order convention carried 192 -> 193 -> 194 -> 195, plus the 191
carried packet coefficients chi_0_star + F_star. Stages 192-195 are the literal neighborhood of the
canonical stage_192 posit-dressed-as-derived failure, and punch-list A1 (chi_Q/Delta_Q/N_Q attributed to
192 when they originate downstream at 194/195) is an OPEN PARTIAL/verdict_logged finding here.
SCOPING: P_2 is the MULTIPOLE-ORDER CONVENTION (the '2' = l=2 quadrupole sector) per its own constraint
notes -- free_choice (sector/branch convention), the most benign subtype. The fit-vs-derive test: is the
P_2='2' convention honestly a sector/branch choice across 192-195, or is any card/notes line dressing it
(or chi_0_star / F_star) as derived/forced? chi_0_star + F_star are carried packet coefficients (191
b_star is internal_consistency, out of scope) -> verify they are honestly carried/posited. This unit is
DISTINCT from the A1 chi_Q overclaim (do not re-litigate A1; just confirm no NEW posit-as-derived on
P_2/chi_0_star/F_star). Orchestrator will independently read stage_191/192/193/194/195 cards + notes.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage191_carried_coefficients_in_packet_interconversion (param=chi_0_star, F_star): redteam_adversarial/provenance/fit_stage191_carried_coefficients_in_packet_interconversion__chi_0_star.yaml
- fit_stage191_carried_coefficients_in_packet_interconversion (param=chi_0_star, F_star): redteam_adversarial/provenance/fit_stage191_carried_coefficients_in_packet_interconversion__f_star.yaml
- fit_stage_192_p_2 (param=P_2): redteam_adversarial/provenance/fit_stage_192_p_2__p_2.yaml
- fit_stage_193_p_2 (param=P_2): redteam_adversarial/provenance/fit_stage_193_p_2__p_2.yaml
- fit_stage_194_p_2 (param=P_2): redteam_adversarial/provenance/fit_stage_194_p_2__p_2.yaml
- fit_stage_195_p_2 (param=P_2): redteam_adversarial/provenance/fit_stage_195_p_2__p_2.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
