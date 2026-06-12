You are an ADVERSARIAL AUDITOR in a layer-2 fit-vs-derive audit of a toy-physics derivation ledger. You are an adversary, not a collaborator: assume the work is wrong; find the fatal flaw. INTERNAL CONSISTENCY IS NOT EVIDENCE OF CORRECTNESS.

FIRST read this rendered prompt IN FULL and ADOPT the FROZEN DIRECTIVE it embeds (your role + fit-vs-derive test + red-flag scan + required-output rules):
  /var/projects/toy_physics/research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage073_aspect_ratio_L_over_a__phase_c_adversarial.md

This is a FAMILY audit (unit `u03_pe_req_chain`, stages 059, 061, 075, 078, 079, 082, 089, 090) in the INTERNAL-CONSISTENCY COMPLETENESS-CRITIC spot-check. Cover EVERY member with its OWN sub-verdict (and per-parameter where a member carries several).

## The test for internal_consistency values (THE INVERTED TEST — READ CAREFULLY)
These values are classified `constraint_kind: internal_consistency` — i.e. the ledger CLAIMS each is a DERIVED consequence (an algebraic identity / Solve / closure / extremum / limit / normalization) of upstream POSITED inputs, NOT itself an inserted choice. Unlike the free_choice sweep, the fatal flaw here is the OPPOSITE mis-classification: **a value that is genuinely a POSITED free_choice (an inserted number / a tuned knob) OR a BACK-SOLVE from an EXTERNAL published number, but is mis-classified `internal_consistency` AND labelled "derived"/"forced"/"exact"/"fixed by"/"a theorem"/"not free" in a paper CARD or per-stage NOTES line** — a posit-as-derived (or fit-as-derived) that ESCAPED the free_choice + published_target sweeps via mis-classification. These candidates were selected BECAUSE their slug carries an overclaim-flavored word (forced/fixed_exactly/backsolve/injection/not_free/_req/unique/...). The slug is the PROMPT to scrutinize, NOT evidence.
VERDICT RULE: **NO** = the value is a GENUINE derivation — open the notes/card and find the actual upstream step that PRODUCES it (a pole condition, a Solve, a closure relation, an extremum, a uniqueness criterion, a parity/symmetry forcing zero, an exact-algebra consequence of a NAMED posited input); the slug's "forced/exact" correctly names that consequence. **YES** = the value is actually inserted/tuned/back-solved (no upstream step produces it; it is chosen to hit a downstream or external target and then relabelled) — give the specific value + the card/notes file+line that calls it derived/forced/exact. Also flag (non-fatal) any stale upstream stage-anchor, value-divergence between members, or a slug that overstates vs the actual (honest) card/notes wording.

## NOTE: this unit is orchestrator-flagged SENSITIVE
It sits on prior-overclaim / mis-classification-risk ground (see Focus: e.g. a literal "backsolve"/"injection"/"5pn data"/"calibration_match" slug, or the stage-192 / chi_Q neighborhoods). Be especially literal: quote the EXACT card/notes line and its label-word, and the EXACT line that shows the derivation step (or its absence). The orchestrator will independently spot-check your verdict against the actual files regardless of what you return.

## Focus
The Pe_req (Peclet-required threshold) chain, carried across 059-090. SLUGS include "forced_not_assumed" (089), "thresholds_exact" (059), "locked_triple" (090), "unique_constructive" (079). KNOWN non-fatal item C1 (punch-list): Pe_fail_chi at 089 carries a FALSE "082" upstream anchor in the script (61-62/92) + CHECKPOINT_CONSTANT_PROVENANCE.md:111; the correct rho-window anchor is "086" (089-notes). That is a stale-anchor metadata item, NOT a mis-classification -- but verify Pe_req itself is a DERIVED threshold (a blocked/unblocked Peclet condition), the same value carried forward, not a tuned number re-inserted per stage. Confirm 059 misattributed_carry is the stale-anchor (metadata) not a value change.

## Member bundles (read each; relative to /var/projects/toy_physics/)
- fit_stage059_pe_req_misattributed_carry (param=Pe_req): redteam_adversarial/provenance/fit_stage059_pe_req_misattributed_carry__pe_req.yaml
- fit_stage059_thresholds_exact (param=Pe_req): redteam_adversarial/provenance/fit_stage059_thresholds_exact__pe_req.yaml
- fit_stage_061_pe_req (param=Pe_req): redteam_adversarial/provenance/fit_stage_061_pe_req__pe_req.yaml
- fit_stage_075_pe_req (param=Pe_req): redteam_adversarial/provenance/fit_stage_075_pe_req__pe_req.yaml
- fit_stage_078_pe_req (param=Pe_req): redteam_adversarial/provenance/fit_stage_078_pe_req__pe_req.yaml
- fit_stage079_unique_constructive_Pe_req (param=Pe_req): redteam_adversarial/provenance/fit_stage079_unique_constructive_pe_req__pe_req.yaml
- fit_stage_082_pe_req (param=Pe_req): redteam_adversarial/provenance/fit_stage_082_pe_req__pe_req.yaml
- fit_stage089_Pe_req_forced_not_assumed (param=Pe_req): redteam_adversarial/provenance/fit_stage089_pe_req_forced_not_assumed__pe_req.yaml
- fit_stage090_locked_triple_Pe_req (param=Pe_req): redteam_adversarial/provenance/fit_stage090_locked_triple_pe_req__pe_req.yaml

## Sources to open as needed
Each provenance bundle carries a `constraints[].constraint_kind` + `evidence_citation` and often a `provenance_findings` block stating the Phase-B synthesis agent's RATIONALE for the internal_consistency call (with cited notes/card lines). Treat that rationale as a CLAIM TO ATTACK, not ground truth — verify it against the actual files. Paper cards: research/pde_ledger/paper/stages/stage_<NNN>.tex ; per-stage notes: research/pde_ledger/notes/stages/moving_throat_pde_stage<NNN>_*.md (EM-projected stages 004-020 notes live at repo-root notes/em_projected/). Open the card AND the notes for each anchor stage; read how the value is LABELED and find the line that DERIVES it. Open what you need to ground each finding.

## REQUIRED OUTPUT (compact, under ~350 words; raw data for an orchestrator, NOT a human-facing message)
1. Per-member verdict table: `candidate_id (param) | YES/NO (mis-classified posit/fit dressed as derived?) | one-clause reason citing BOTH the label line AND the derivation-step line (or its absence)`.
2. Unit verdict: fatal flaw? YES/NO. If YES, the single most important one in one sentence with the specific value + the overclaiming card/notes file+line + why no upstream step derives it.
3. Classification assessment (the crux): for each member, is `internal_consistency` CORRECT (a genuine derivation exists) or should it be free_choice / published_target (a posit / external back-solve)? Name any mis-classification stage+line, or state the classification is sound.
4. Proof you looked: which card+notes files you opened for each anchor stage, the exact label-word, and the exact derivation-step line you found (or confirmed absent).
Do NOT propose fixes or next steps. Pure falsification pass. DO NOT run scripts - read and reason.
