# Batch I.1 (stages 001-012) — pde_ledger SECOND PASS

Date: 2026-06-04

## Method

- v2 paper-grounded auditor PLUS a new **exhaustive script→doc
  value-reconciliation augmentation** (see
  `redteam/pass2/RECONCILIATION_AUGMENTATION.md`).
- One clean per-stage audit agent per stage (no context contamination).
- Independent exec reliability gate (re-run both engines; refresh committed
  outputs).
- Codex as the fix-applier (orchestrator does not hand-apply directive fixes).
- Clean per-stage verify agents.

## Result

All 12 stages reached `verified` at v2 depth + value-reconciliation
augmentation. **12/12 verified, `material_change=false` on all 12. No
stop-cold.**

| Outcome | Stages |
| --- | --- |
| Clean (9) | 001, 002, 004, 005, 007, 008, 010, 011, 012 |
| Findings → resolved → verified (3) | 003, 006, 009 |

## The three resolutions

- **006** — card Ampère sign reconciled to the engine form:
  eq:stage006-ampere `∇×H − ∂_tD = μ₀J + L_mix`
  → `−∇×H − ∂_tD + L_mix = μ₀J` (documentation sign typo; contradicted the
  card's own Gauss law + both engines; nothing downstream consumes the vector
  Ampère sign). Codex-applied, Claude-reviewed.
- **003** — notes grouped-P2 index garble `d_{237/238/239}` → `d_{2,20/21/22}`
  (fixed), plus the `.wl` `lRed` consolidated to one parenthesized assignment
  (latent construction fragility; math no-op, output byte-identical). The
  apparent "Mathematica Lagrangian doubling" alarm was an **exec-DISPROVEN
  false positive**: the multi-line `lRed` silently drops 4 lines via WL
  newline-before-minus parsing, and the parenthesized re-add compensates, so
  net `lRed` is correct.
- **009** — generic-μ₁ coverage EXTENDED: a generic-kernel first-moment (μ₁)
  check added in BOTH engines (Gamma kernel w(u)=u·e^{-u}, μ₁=2; O(ℓ) leading
  term asserted = q0+ℓ·μ₁·q1), so the card's general μ₁ formula is now
  exercised, not just the exponential special case (μ₁=1).

## Checkpoints

001, 002, 003 re-verified (with value-reconciliation augmentation). 003's
"Lagrangian doubling" was an exec-disproven false positive; no trust impact on
the carried results.

## Status

12/12 verified, all `material_change=false`; no stop-cold. Paper/notes edits
(stage 006 `.tex`, stage 003 notes) Codex-applied + Claude-reviewed and closed
this batch (not deferred). See `notes/PAPER_CLEANUP_TRACKER.md` (P5-02).
