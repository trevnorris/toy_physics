You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u02_barrier_benchmark_Vknown`, stages 222, 223) of **free_choice** ansaetze. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for free_choice values (READ CAREFULLY)
These values are classified `free_choice`: posited modeling choices, NOT claimed external matches. So the GR/CODATA benchmark test does NOT bite here and `Benchmarks: []` is expected and correct. The fatal-flaw you are hunting is the **posit-dressed-as-derived overclaim** (the exact failure mode found at stage_192): a value that is genuinely a CHOICE but is labeled "derived" / "forced" / "exact" / "non-tunable" / "a theorem" / "fixed by" in a paper CARD or a per-stage NOTES line, or is silently propagated as if it were a consequence. A free_choice value that is HONESTLY disclosed as posited/carried/a reference freeze/an open choice is NOT a flaw. Also flag (non-fatal) any stale upstream stage-anchor or value-divergence between members.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It has prior overclaim history / borderline-classification status (see Focus). Be especially literal: quote the EXACT card/notes line and its label-word for every value. The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The V_known illustrative barrier benchmark (=1.181909222592) and its derived requirement DeltaV_req.
THE AUDIT-CORRECTION SITE: the scanner candidate-hint flagged V_known "likely published_target (external
barrier import)", but the audit correction RESOLVED it free_choice -- stage-222/223 notes REPEATEDLY and
EXPLICITLY disclaim external provenance, calling it an ILLUSTRATIVE reduced-barrier benchmark with no
external source named and no back-solve shown. Test here is NOT "is it external" (settled NO) but: does
ANY stage_222/223 card or notes line dress the illustrative V_known value as derived/forced/external/
published? A provenance-gap (the float's genealogy below 222 is undocumented) is NOT an overclaim as long
as it is disclosed illustrative. DeltaV_req(1)=V_known(1)-epsilon=1.081909222592 is SYMBOLICALLY derived
from V_known (the survival inequality is internally derived) -> verify DeltaV_req only inherits V_known's
illustrative status, no independent overclaim. Orchestrator will read stage_222/223 cards + notes.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage222_illustrative_barrier_benchmark (param=V_known): redteam_adversarial/provenance/fit_stage222_illustrative_barrier_benchmark__v_known.yaml
- fit_stage_222_deltav_req (param=DeltaV_req): redteam_adversarial/provenance/fit_stage_222_deltav_req__deltav_req.yaml
- fit_stage223_barrier_benchmark_and_eta_window (param=V_known): redteam_adversarial/provenance/fit_stage223_barrier_benchmark_and_eta_window__v_known.yaml

## Sources to open as needed
Paper cards research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (em_projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage and read how the value is LABELED. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (fatal overclaim?) | one-clause reason citing the card/notes line that labels it`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific posited value + the overclaiming card/notes file+line that calls it derived/forced/exact.
3. Disclosure assessment (layer-2 crux): for each member, is the free_choice value honestly labeled (posited/carried/reference-freeze/open) or overclaimed (derived/forced/exact/non-tunable)? Name any overclaiming stage+line, or state disclosure is adequate.
4. Proof you looked: which card+notes files you opened for each anchor stage and which exact label-words you checked.
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
