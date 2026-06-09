---
unit_id: 206
batch: VI.1
created_at: 2026-06-09T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 206

Two of the three findings are `paper_misalignment` items requiring user resolution
(F1, F2). The orchestrator is HOLDING for user resolution — do NOT edit paper.tex,
notes/, or the script's stage-number labels until the user has chosen a direction.
F3 is a `stale_output` informational item handled by the verifier's re-run, not a
Codex edit.

Apply nothing on this unit until a follow-up directive authorizes a direction.

## F1 — paper_misalignment (card Verification field is stale)

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_206.tex:11` quote: "SymPy audit: \StageFile{scripts/...stage206...sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl` — a complete, fully-passing Mathematica audit exists (all M1-M7, F3a, F3b PASS). The appendix `\claimstatus` at `paper/appendices/stage_appendix_part06.tex:918` also omits the second engine.

## Resolve before fix_loop

The card says "Mathematica audit: none yet" but a passing Mathematica `.wl` was added in the pass-1 dual-engine retrofit. Should the card be updated to cite it?

Possible directions (the user picks one):
- (a) Card is stale → update `stage_206.tex:11` to cite `\StageFile{mathematica/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_mathematica_audit.wl}` in place of "none yet" (paper-side edit; Codex does not auto-apply paper edits without authorization).
- (b) The `.wl` is considered provisional / not yet to be advertised → leave the card and record the rationale.

## F2 — paper_misalignment (Stage 239 → Stage 205 label)

**Subtype:** notes_contradicts_script

**Notes side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem.md:200-208` quote: "If the curvature envelope collapses to a point, \(\underline K_1=\overline K_1=L_1\) … which is exactly the **Stage 205** quadratic logarithmic predictor written in the oriented variables."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py:131` comment "V. Collapse to the Stage 239 quadratic log predictor"
- `:133` `subbanner("V. Collapse to the Stage 239 quadratic log predictor")`
- `:145` `expect_zero("Stage 206/239 log-predictor collapse", ...)`

Stage 239 is "Rigid-mouth physical normal form" — unrelated. The verified identity is correct; only the stage label is wrong (numbering-drift artifact). Canonical target per the notes is Stage 205.

## Resolve before fix_loop

Section V of the SymPy script labels the collapse target "Stage 239"; the notes (the canonical carrier) say the collapse is to the Stage 205 quadratic log predictor. Confirm the correct upstream stage so the label can be corrected.

Possible directions (the user picks one):
- (a) Notes are correct (expected) → change the three script labels at `:131`, `:133`, `:145` from "Stage 239"/"206/239" to "Stage 205"/"206/205"; re-run SymPy. No math change.
- (b) The script is correct and the notes/narrative are wrong → flag for deeper review (would contradict the entire Stage 204/205/206 narrative).

The orchestrator will not invoke Codex on this unit until the user has chosen a direction for F1 and F2.

## F3 — stale_output (informational; verifier handles)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.txt`

**Issue:** The committed SymPy output (mtime 2026-06-02 11:32) predates the `.py` (2026-06-03 15:59) and still banners "STAGE 189 — …" / "STAGE 189 SYMPY AUDIT PASSED" while the current `.py` banners "STAGE 206". All numeric/symbolic results in the stale output still agree with the current script (every check `= 0`), so this is informational. The verifier's independent re-run will refresh it; the refresh must occur AFTER any F2 label fix so the section-V labels in the `.txt` also update.

**Required change:** none by Codex. The orchestrator/verifier re-runs `python3 …sympy_audit.py` and recommits the refreshed `.txt`.

**Verification:** refreshed `.txt` banner reads "STAGE 206 …"; mtime newer than the `.py`.
