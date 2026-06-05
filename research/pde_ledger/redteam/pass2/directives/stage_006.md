---
unit_id: 006
batch: I.1
created_at: 2026-06-04T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-04T18:20:30-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 006 (RESOLVED; user-approved 2026-06-04)

The single finding is a `paper_misalignment` the user has resolved. This directive AUTHORIZES a paper-side (`.tex`) edit. No script change. After editing, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation`. There is no script to run for this finding — do NOT modify or run the scripts (both engines are already correct).

## F1 — paper_misalignment (projected Ampère sign) → RESOLVED direction (b): fix the CARD

**Decision evidence (orchestrator cross-stage investigation, 2026-06-04):** both independent engines verify `−∇×H_flux − ∂_t D_flux + L_mix = μ₀ J_proj`. This engine sign is the honest consequence of the script's component map `G^{i0}=D_i`, `G^{ij}=ε^{ijk}H_k`, which ALSO yields the card's own Gauss law `+∇·D = +μ₀ρ`; the stage-006 script even defines `curlH_i := −(∇×H)_i` (sympy/wl). The disputed `∇×H` sign and the H↔F/D map live ENTIRELY in stage 006 — stages 005/008 stay covariant (`F^{μν}`), 010/011/012/023 are sign-agnostic P2/bundle moments, and the EM-extension endgame (243/244/247) rides a scalar (continuity-leakage / Poynting-work) branch that does not consume the vector Ampère sign. Nothing downstream relies on `+∇×H`. The card's displayed `+∇×H` is the lone inconsistency — it contradicts its own script, its own Gauss law, and both engines. Therefore the CARD is corrected to the engine form (NOT the engines).

**Target:** `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_006.tex:38-41` (eq:stage006-ampere)

**Required change:** change the displayed projected Ampère law from
`\nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux} = \mu_0\mathbf J_{\rm proj}+\mathbf L_{\rm mix}`
to the engine-verified form
`-\nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux}+\mathbf L_{\rm mix} = \mu_0\mathbf J_{\rm proj}`
i.e. flip the sign of the `\nabla\times\mathbf H_{\rm flux}` term and move `\mathbf L_{\rm mix}` to the LHS. Then SCAN the rest of `stage_006.tex` (especially any `\stagefield{Output}` sentence or surrounding prose that states the Ampère/curl-H sign convention) and make it consistent with the corrected equation. Do NOT change the Gauss law (`+∇·D=+μ₀ρ` is already correct and must stay). Quote each edited line in the `## Applied` block.

**Also check (low-risk):** the part appendix row `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` — find the stage-006 row; if it restates the Ampère sign, make it consistent; if it is only a status/summary line, leave it. Report which.

**Verification:** the card's eq:stage006-ampere now reads `−∇×H_flux − ∂_t D_flux + L_mix = μ₀ J_proj`, matching both engines' load-bearing assertion and the card's own Gauss law. No script changed; `material_change: false` (no downstream re-stale). A clean verify agent will confirm card↔engine alignment.

## Applied: F1

- files_changed:
  - `paper/stages/stage_006.tex`
  - `redteam/pass2/directives/stage_006.md`
- summary: Corrected the projected Ampere law to the engine-verified curl-H sign with the leakage term on the LHS; the appendix stage-006 row is only a status/summary line and was left unchanged.
- deviation: none
- edited_lines:
  - `-\nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux}`
  - `+\mathbf L_{\rm mix}=\mu_0\mathbf J_{\rm proj}.`
