---
unit_id: 216
batch: VI.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 216

Apply each non-paper_misalignment finding below in order. There are NONE — the sole
finding is a `paper_misalignment` pending user resolution. Codex applies nothing on
this unit until the user chooses a direction in a follow-up directive. Do not edit
paper.tex, notes/, or scripts.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale verification metadata)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_216.tex:11` quote:
  "SymPy audit: \StageFile{scripts/...stage216..._sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl` — a complete Mathematica audit that DERIVES the stage's load-bearing objects (Lagrange-Solve gradient optimum, matrix-quadratic-form leverage, spectral `w_Σ≤4` bound, quadratic-Solve bracket) and PASSES all six check families (saved output: "All Stage 216 Mathematica audit checks passed.").

## Resolve before fix_loop

The stage-216 card's `\stagefield{Verification}` line states "Mathematica audit: none yet," but a passing Mathematica audit `.wl` now exists (added in the pass-1 dual-engine retrofit, dated 2026-06-02). The sibling stage 218 card already cites its `.wl` in this same field. This is a stale documentation-metadata mismatch, not a math disagreement.

Possible directions (the user picks one):
- (a) Card is stale → update `paper/stages/stage_216.tex:11` to name the `.wl` exactly as stage 218 does, e.g. append "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl}." No script change. (This is a paper-side edit Codex performs only under an explicit follow-up directive.)
- (b) The `.wl` is considered out of scope for this card and should not be cited → leave the card; no change.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
