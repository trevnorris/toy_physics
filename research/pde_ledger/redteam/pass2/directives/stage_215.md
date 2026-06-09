---
unit_id: 215
batch: VI.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 215

There is exactly one finding, and it is a `paper_misalignment` pending user resolution. Codex applies NOTHING on this unit until the user chooses a direction. Do not edit `paper/`, `notes/`, or the scripts to "fix" this. The orchestrator holds before invoking Codex.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_215.tex:11` quote: "\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_sympy_audit.py}.  Mathematica audit: none yet.}"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl:1-256` — a complete, passing Mathematica audit (blocks M1-M7) exists; saved output ends "Stage 215 Mathematica audit passed."

## Resolve before fix_loop

The card's `\stagefield{Verification}` line still says "Mathematica audit: none yet," but a passing dual-engine Mathematica audit (`mathematica/...stage215..._mathematica_audit.wl`) was added in the pass-1 retrofit and is verified INDEPENDENT (real quantifier elimination vs. the SymPy discrete enumeration). The card understates verification coverage.

Possible directions (the user picks one):
- (a) Paper card is stale → update `stage_215.tex:11` so the Verification line cites the `.wl` audit (e.g. "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl}.") and drop "none yet". No script change. (This is a paper-side edit; per file-ownership it is applied by Codex only after the user authorizes a follow-up directive, Claude reviews.)
- (b) The `.wl` is intentionally not yet "blessed" for the card → leave the card as-is and record the deferral; no edit. (Inconsistent with the rest of the pass-2 dual-engine retrofit reconciliation, but the user's call.)

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
