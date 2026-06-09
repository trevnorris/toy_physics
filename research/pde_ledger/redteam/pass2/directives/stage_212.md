---
unit_id: 212
batch: VI.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 212

The only finding is a `paper_misalignment`. Codex applies NOTHING here — the
orchestrator holds for user resolution. Do not edit `paper/`, `notes/`, or scripts
until the user picks a direction in a follow-up directive.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale verification-status line)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_212.tex:11` quote:
  "\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py}.  Mathematica audit: none yet.}"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl:1-206` — a full Mathematica audit (M1-M4) that PASSES.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.txt` — fresh transcript ending "Stage 212 Mathematica audit passed."

## Resolve before fix_loop

The card's `\stagefield{Verification}` line says "Mathematica audit: none yet," but
a passing dual-engine Mathematica audit `.wl` is present (added in the pass-1
retrofit) with a fresh passing output. The card status line is stale.

Possible directions (the user picks one):
- (a) Card is stale → update `paper/stages/stage_212.tex:11` to cite
  `scripts/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_sympy_audit.py`
  AND `mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl`
  (drop "none yet"). No script change. (This is the expected resolution.)
- (b) The `.wl` was not meant to ship for this stage → remove/quarantine it and keep
  the card as-is. (Unlikely; the dual-engine policy requires a `.wl` wherever
  Mathematica can independently verify, which it can here.)

The orchestrator will not invoke Codex on this unit until the user has chosen a
direction. Note: this is a prose-status edit only — no math change either way; the
engines already agree and all checks pass.
