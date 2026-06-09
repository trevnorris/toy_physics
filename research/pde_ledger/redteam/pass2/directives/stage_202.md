---
unit_id: 202
batch: VI.1
created_at: 2026-06-09T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 202

The only finding is a `paper_misalignment` (F1). Codex does nothing — the orchestrator is holding for user resolution. Do NOT edit `paper/`, `notes/`, or any script to "fix" this until the user picks a direction in a follow-up directive.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale verification line)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_202.tex:11` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage202_free_quintuple_target_graph_sympy_audit.py}.  Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage202_free_quintuple_target_graph.md:549-552` quote (Supporting files): lists only `…_sympy_audit.py` and `…_sympy_audit_output.txt`; no `.wl` named.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl:1-252` — a complete Mathematica audit exists and passes (output `…/mathematica/output/…_mathematica_audit.txt` shows M0–M7 all PASS, `Exit[0]`).

## Resolve before fix_loop

The card's `\stagefield{Verification}` line and the notes' Supporting-files list say there is no Mathematica audit for stage 202, but a full passing `.wl` is present (added in the pass-1 dual-engine retrofit). Should the prose be updated to record the Mathematica audit?

Possible directions (the user picks one):
- (a) Prose is stale → update `stage_202.tex:11` to read "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage202_free_quintuple_target_graph_mathematica_audit.wl}." and add the `.wl` to the notes Supporting-files list (lines 549-552). No script change. (Expected direction.)
- (b) The `.wl` should not be claimed by this stage → leave prose as-is and remove/relocate the `.wl` (a script-side decision, unlikely given it passes and matches the paper claims exactly).

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. There are NO non-paper_misalignment findings, so there is nothing for Codex to apply regardless.
