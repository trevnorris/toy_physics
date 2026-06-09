---
unit_id: 173
batch: V.1
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T16:49:03-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 173

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what the finding names. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all `expectZero` checks passing. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl:26-99`

**Issue:** The `.wl` is a line-by-line port of the SymPy script (`scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py`): identical lane-expansion variable names (`d0A`, `d2A`, `d4A`, `n0A`), identical `Series`→order-1-`Coefficient` extraction mirroring SymPy's `series(...).diff(eps).subs(eps,0)`, identical canonical substitutions `d2 -> -(1/9) d0` and `d4 -> -(1/27) d0`, identical `Solve[...]/First` choreography for the even-preserving collapse, and a byte-for-byte identical carry-forward `Print` block (wl L87-93 vs py L113-119). This violates the second-engine independent-derivation policy: a shared conceptual error would pass both engines unchanged.

**Required change:**
Re-author the derivation portion of the `.wl` (lines 32-84) so Mathematica reaches the slope laws by an independent route instead of echoing the SymPy `Series`/`Coefficient` choreography. Keep the six `expectZero` checks and their target right-hand sides exactly as they are now (those are the paper's deliverables and must continue to pass). Recommended independent route:

1. Replace the `d0A/d2A/d4A/n0A` named-ratio + `Series`/`Coefficient` pipeline with direct first-order implicit differentiation of the defining response relations evaluated at `eps = 0`. For example:
   - Define the responses as functions of `eps` only at the point of differentiation, e.g.
     `u21 = FullSimplify[(D[-(d2 + eps*lam*d21)/(d0 + eps*lam*d01), eps] /. eps -> 0)/lam];`
     `u41 = FullSimplify[(D[((d2 + eps*lam*d21)^2 - (d0 + eps*lam*d01)*(d4 + eps*lam*d41))/(d0 + eps*lam*d01)^2, eps] /. eps -> 0)/lam];`
     `p1  = FullSimplify[(D[(n0 + eps*lam*n01)/(d0 + eps*lam*d01), eps] /. eps -> 0)/lam];`
   - This computes the same first-order slopes by analytic `D[...]`-at-`eps=0` rather than truncated-series-coefficient extraction, and does NOT reuse the SymPy intermediate naming `d0A/d2A/...`.
2. Use distinct local naming for any helper intermediates (do not mirror the `.py`'s `u21_can`, `u41_can` step-for-step beyond what the assertions require; e.g. inline the canonical substitution into the `expectZero` argument).
3. Remove the verbatim carry-forward `Print` block (wl L87-93) or replace it with a Mathematica-native one-line summary; it must not be a byte-copy of the `.py` block.
4. Leave all six `expectZero` calls (`u2 slope identity`, `u4 canonical formula`, `P1/P0 formula`, `hidden-even residual`, `D21 + D01/9`, `D41 + D01/27`) and the `Xi_load` / lane prints intact, with their current target expressions, so the verified claims are unchanged.

Self-test confirmation (already checked by auditor): the `D[ratio, eps]/.eps->0` route depends genuinely on `eps` through the `eps*lam*d01` terms, so each first-order slope is non-trivial (matches the current nonzero general-form outputs); no new constant is introduced, so the paper round-trip is preserved.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 173` and confirms (a) the script exits 0 with all six `expectZero` checks PASS and the same residual-0 outputs, and (b) the `.wl` no longer reproduces the `.py`'s named `d0A/d2A/d4A/n0A` choreography or the byte-identical carry-forward print block.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl`
- summary: Replaced the transliterated series/coefficient derivation with direct Mathematica differentiation at eps=0, renamed helper intermediates, derived the even-preserving collapse by linear coefficient extraction, and replaced the copied carry-forward block with a one-line summary.
- deviation: none
