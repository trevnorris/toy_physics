---
unit_id: 003
batch: I.1
created_at: 2026-06-04T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T00:24:19Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 003 (RESOLVED; user-approved 2026-06-04)

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with `files_changed`, `summary` (one sentence), and `deviation` (or "none"). Do NOT introduce features, refactors, or stylistic changes beyond what each finding specifies. After any SCRIPT edit, RUN the affected script under `timeout 600 math -script <path>` and iterate until it exits 0 with all in-file checks passing.

> **NOTE — the original F1 (Mathematica "Lagrangian doubling") and F2 (stale output) are WITHDRAWN.** The orchestrator confirmed by direct execution + instrumentation that the current `.wl` runs (exit 0), passes all 19 checks, and is byte-identical to the committed transcript. The `lRed = lRed + (...)` block is NOT a doubling: the original multi-line `lRed` (lines 54-59) captures only the *kinetic* terms because lines 56-59 each begin with `-` and the WL script parser treats them as separate, discarded statements; the parenthesized re-add then supplies the stiffness/BdG/coupling exactly once. Net `lRed` is the correct single-coefficient Lagrangian. There is no math defect. F4 below hardens this fragile construction; F3 fixes a genuine notes garble.

## F3 — paper_misalignment (notes index garble) → RESOLVED direction (a) [AUTHORIZED notes edit]

**Subtype:** notes_contradicts_script. The user chose direction (a): the notes indices are garbled; correct them to match the `.tex` card and both scripts (already correct). This directive AUTHORIZES this notes edit.

**Target:** `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage003_bdg_coupling.md:367-372`

**Required change:** In the grouped trace / anisotropy formulas, replace the garbled superscript indices `d_{237}^{\rm eff}`, `d_{238}^{\rm eff}`, `d_{239}^{\rm eff}` with the correct grouped-P2 channel indices `d_{2,20}^{\rm eff}`, `d_{2,21}^{\rm eff}`, `d_{2,22}^{\rm eff}` (channels A ∈ {20,21,22}), matching `paper/stages/stage_003.tex:131-140` and `scripts/...stage003_bdg_sympy_audit.py:321-323` (`d220/d221/d222`). Change ONLY the three index tokens in those formulas; do not alter any other notes text, equation structure, or values. No script change. No script run for this finding.

**Verification:** the three formulas in `md:367-372` now read `d_{2,20}^{eff}/d_{2,21}^{eff}/d_{2,22}^{eff}`; no other byte of the notes file changed.

## Applied: F3

- files_changed:
  - `notes/stages/moving_throat_pde_stage003_bdg_coupling.md`
- summary: Corrected the grouped-P2 effective coefficient indices from 237/238/239 to 2,20/2,21/2,22 in the three formulas.
- deviation: none

## F4 — robustness consolidation of `lRed` (user-approved) [SCRIPT edit, re-run required]

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl:54-67`

**Issue (latent fragility, not a current math error):** the multi-line assignment `lRed = <line55> <line56> ... <line59>;` (lines 54-59) silently drops lines 56-59 — each starts with `-`, so WL parses them as separate discarded statements; the original `lRed` captures only the kinetic terms. Correctness is rescued ONLY by the subsequent `lRed = lRed + ( ... );` re-add (lines 62-67), which a future "remove the redundant duplicate" cleanup would naturally delete — silently breaking the stage.

**Requirement (Codex designs the exact edit):** rewrite the `lRed` definition as a SINGLE, robust, fully-parenthesized assignment that captures ALL terms exactly once — the complete Lagrangian: qa/qL inertia (`maa,maL,mLL`), the K-potential (`-1/2 kaa qa^2 - kaL qa qL - 1/2 kLL qL^2`), both BdG kinetic+potential terms (`xa`,`xb`), and all four wall–matter couplings (`c1a,c1b,c2a,c2b`) — and remove the now-unnecessary `lRed = lRed + (...)` re-add block. KEEP the `Clear[qa0, qL0, xa0, xb0, vqa, vqL, vxa, vxb, aqa, aqL, axa, axb];` statement (its locals are used downstream by `lAlg`/`timeD`/`backEL`).

**Acceptance criteria (MUST hold — this is a no-op on the math):**
- The net `lRed` is algebraically identical to the current (post-re-add) value: single-coefficient full Lagrangian.
- A fresh run `timeout 600 math -script <path>` prints the SAME 19 PASS lines (qa/qL/xa/xb equations = 0, derived D0_eff vs Deff = {{0,0},{0,0}}, series match, dispersion/roots, Sections II–IV) and exits 0.
- The regenerated output `.txt` is byte-identical to the current committed transcript (the orchestrator will diff it). If it is NOT byte-identical, the consolidation changed the math — STOP and report rather than "fixing" the assertions.

**Verification:** orchestrator re-runs `redteam exec-mathematica 003`, diffs the new transcript against the prior committed one (must be identical), and a clean verify agent confirms the single-assignment `lRed` captures every term once.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- summary: Consolidated `lRed` into one parenthesized full-Lagrangian assignment and removed the compensating re-add block.
- deviation: none
