---
unit_id: 210
batch: VI.1
created_at: 2026-06-09T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 210

This directive contains NO Codex-applied script edits. One finding is a
`paper_misalignment` held for user resolution (F1); the other is an
informational `stale_output` resolved by a re-run, not a script edit (F2).
Codex must apply nothing until the user resolves direction on F1.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_210.tex:11` quote:
  "SymPy audit: \StageFile{scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl:1-317` — a Mathematica audit exists, runs, and passes; the committed output ends "STAGE 210 MATHEMATICA AUDIT PASSED".

## Resolve before fix_loop

The card's Verification field says "Mathematica audit: none yet", but a passing,
independent Mathematica audit `.wl` is present (added in the pass-1 dual-engine
retrofit). Which way should this be reconciled?

Possible directions (the user picks one):
- (a) Paper is stale → update card line 11 to cite
  `mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl`
  (e.g. "Mathematica audit: \StageFile{mathematica/...mathematica_audit.wl}."). Paper-side edit, no script change. (Expected resolution — the `.wl` is real and independent.)
- (b) The `.wl` should not count → leave the card, but then the dual-engine
  requirement for this stage is unmet (contradicts feedback-dual-engine-required). Unlikely.

The orchestrator will not invoke Codex on this unit until the user has chosen a
direction. Codex must not edit paper.tex/notes/ for F1.

## F2 — stale_output (informational; no script edit)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.txt`

**Issue:** The committed SymPy output (mtime 2026-05-11) predates the SymPy
script (mtime 2026-06-03) and prints the stale banner `STAGE 193` (lines 11,
173) whereas the current `.py` prints `STAGE 210` (lines 35, 256). This is the
known P4-52 stale-banner set: the output predates the banner fix. The check
substance (all residuals `= 0`) is unaffected. The Mathematica output is fresh.

**Required change:** None to the script. Re-run
`python3 scripts/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_sympy_audit.py`
and recommit the refreshed `.txt`; the orchestrator's independent re-run does
this. Do NOT edit the `.py`.

**Verification command:** After the orchestrator re-runs `redteam exec-sympy 210`,
confirm the refreshed `.txt` shows the `STAGE 210` banner at top and bottom and
exits 0.
