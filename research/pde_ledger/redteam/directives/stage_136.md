---
unit_id: 136
batch: IV.4
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 136

This directive contains a single `paper_misalignment` finding. Codex must NOT auto-apply any edit. The orchestrator is to hold this unit until the user resolves the question below. Codex never edits `paper/`, `notes/`, or other prose documents.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script (closest fit; the notes' internal provenance is contradictory and the paper card does not anchor the carry-forward constants to a named upstream stage)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_136.tex:9` quote: "This stage imports the shell/mixed core, the mouth source law, outlet consistency, core-to-mouth gain maps, and self-matched susceptibility closure."
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_136.tex:11` quote: "SymPy audit: none yet.  Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:1` quote: "# Moving-Throat PDE — Stage 238: Mouth-Layer Fixed-Point Status After the Coupled Solve"
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:7` quote: "After Stages 184–186, the mouth-source selection problem has narrowed again."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:30` quote: "M_s \approx 1.50882951349316 - 0.658075937605429\,M_q."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md:40` quote: "\Sigma_m^* \approx 0.451485277739090."

**Script side:**
- (no script exists; this is `is_status_only_candidate: True`)

## Resolve before fix_loop

The notes file backing stage 136 has internally inconsistent provenance: the H1 header (L1) labels the file "Stage 238", the body (L7) cites "Stages 184–186" as the upstream narrowing, and the filename / paper card label the unit as stage 136. Two load-bearing numerical constants (`M_s ≈ 1.50882951349316 − 0.658075937605429 M_q` at L30, and `Sigma_m^* ≈ 0.451485277739090` at L40) are quoted as carry-forward from this upstream chain, but the paper card's `\stagefield{Inputs}` (L9) does not name a specific upstream stage that derived them. As a status-only ledger entry, stage 136 must point to a stage where these were actually located; right now nothing in the paper card or notes does so unambiguously.

The user must pick one of the following directions:

Possible directions:
- (a) The notes file was originally written for a different stage (e.g., 238) and renamed without rewriting its body. Action: rewrite the notes H1 and the "Stages 184–186" reference to point to the actual upstream stages whose scripts numerically located `(M_s, M_q)` and `Sigma_m^*`, and add those stage indices to the card's `\stagefield{Inputs}` line.
- (b) The status entry at stage 136 is forward-looking and the constants are not yet derived anywhere — they are conjectural carry-forward placeholders. Action: mark them as such in the card (`\claimstatus` already flags `\StatusOpen{}` but the body block does not distinguish derived from placeholder values), and downstream consumers (stages 146–153 per card L27) must inherit the status tag, not the numeric values, until a stage actually derives them.
- (c) The constants are derived in an upstream stage whose audit is not yet under is_status_only_candidate coverage — in which case the carry-forward chain is broken at the level of the auditor's manifest, not the paper. Action: identify the upstream unit whose script located these constants and either add a `\stagefield{Inputs}` cross-reference here or convert this unit to a non-status-only audit that re-locates the fixed points in-script.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No script edits are prescribed; no `paper/` or `notes/` edits are authorized for Codex.
