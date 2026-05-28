---
unit_id: 132
batch: IV.4
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 132

This directive contains only a `paper_misalignment` finding. Codex must NOT apply any edits to paper/, notes/, or scripts/ for this unit until the user has resolved the question below in a follow-up directive. The orchestrator should halt before invoking Codex.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script (notes' carry-forward chain is internally inconsistent with this unit's stage number; no script exists to compare against, but the carry-forward chain is the verification surface for a status-only unit)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_132.tex:1` quote: `\section[Stage 132]{Stage 132: Mouth Boundary-Layer Status After Explicit Source-Law Extraction}`
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_132.tex:11` quote: `\stagefield{Verification}{SymPy audit: none yet.  Mathematica audit: none yet.}`
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:1298` quote: `\input{stages/stage_132}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md:2` quote: `# Moving-Throat PDE — Stage 234: Mouth Boundary-Layer Status After Explicit Source-Law Extraction`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md:4` quote: `After Stages 180–182, the mouth-source side is no longer an abstract profile problem.`

**Script side:**
- (none; no script exists for unit 132)

## Resolve before fix_loop

The notes file `notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md` is identified by filename and slug with Stage 132, but its H1 says "Stage 234" and its body attributes the load-bearing source-law results (explicit GNLS/localized-Maxwell source family, Family-1 bias \(\mathfrak g_\Pi\), monotonicity, \(\Pi_*\approx 1.50882951349316\), and the parent threshold identity) to "Stages 180–182." Stages 180–182 lie *after* 132 in the ordering and cannot be upstream of this status card. Which is correct?

Possible directions (the user picks one):

- (a) **Notes were misfiled** — the notes really belong to Stage 234, not 132. Move the file to `notes/stages/moving_throat_pde_stage234_mouth_boundary_layer_status.md`, fix any cross-references, and supply (or commission) a different notes file for stage 132 that points to genuine upstream sources strictly less than 132. No script change.
- (b) **Notes were stale, stage 132 is the correct anchor** — the carry-forward chain referenced as "Stages 180–182" should actually be lower-numbered upstream stages (the auditor cannot enumerate them because it cannot read other units' files). Update the notes to (i) correct the H1 to "Stage 132" and (ii) re-attribute the upstream chain to the actual lower-numbered units that perform the GNLS/localized-Maxwell source-law extraction. No script change yet, but a follow-up `missing_sympy`/`missing_mathematica` audit may be warranted once the corrected upstream chain is identified.
- (c) **The stage is genuinely conjectural at this position in the ledger** — the notes' results are forward-looking and stage 132 is a placeholder. Update `paper/stages/stage_132.tex` `\claimstatus` and `Inputs` lines to reflect that the listed checks (positivity, normalizations, Family-1 compensation) are *targets*, not validated, and update the boxed prose claim accordingly. No script change.
- (d) **Renumber the stage** to its intended ledger position (likely 234 given the notes H1) and propagate the rename through `paper/stages/`, `paper/appendices/stage_appendix_part04.tex`, the notes filename, and any downstream `\stagefield{Downstream use}` cross-references that mention this card. No script change.

Until the user picks a direction, Codex must not edit `paper/stages/stage_132.tex`, the notes file, or create scripts under the `stage132` slug.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
