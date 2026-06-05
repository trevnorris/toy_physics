---
unit_id: 080
batch: III.4
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 080

## ⛔ ORCHESTRATOR DECISION — ENTIRE DIRECTIVE DEFERRED, NO CODEX INVOCATION

The single finding below flags five comment/docstring labels (`:5,6,27,35,65`) that all cite `Stage-61`/`Stage 62`. **These are CROSS-stage references to OTHER stages (078 and 079), not stage-080 self-labels.** Stage 080's own self-labels are already canonical (docstring line 3 reads `SymPy audit for Stage 080.`; banner line 21 reads `STAGE 080`). Per the settled Reading-2 in-loop policy, the red-team loop fixes only a stage's own UNAMBIGUOUS self-labels; cross-references are owned by the dedicated `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (content-keyed, never offset-swept) and are LEFT UNTOUCHED here. The auditor over-flagged (the same docstring-cross-ref class was correctly marked CLEAN on sibling stage 078 this batch). **Codex is NOT invoked on 080.** The committed transcripts are refreshed by the orchestrator's independent exec re-run; the verify pass treats 080 as deferred-clean (no source change). The original finding is retained below for the audit trail.

## F1 — stale cross-stage label fix (label-only, comments/docstring) — DEFERRED, NOT APPLIED

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:5,6,27,35,65`

**Issue:** The sympy docstring and four inline comments cite the upstream sources as "Stage-61"/"Stage 62", which are stale pre-renumber labels from the +17 EM-extension realignment (61→078, 62→079). The canonical upstream sources, per `paper/stages/stage_080.tex:7` (`the Stage~078 Peclet windows`) and the Stage-079 demand map / `kappa_F1 = 12321/5`, `eta_F1 = 37` constants, are **Stage 078** (the four Pe thresholds) and **Stage 079** (the `A_F1 Omega_Pe^2` demand map). No math changes; the constants and formulas in the script are already correct.

**Required change:**
Edit only these comment/docstring lines:
- line 5: `Converts the Stage-61 Family-1 Pe_req thresholds into explicit quadrupole-demand` → `Converts the Stage-078 Family-1 Pe_req thresholds into explicit quadrupole-demand`
- line 6: `thresholds zeta_req using the Stage-62 Family-1 demand map.` → `thresholds zeta_req using the Stage-079 Family-1 demand map.`
- line 27: `# Family-1 constants from Stage 62.` → `# Family-1 constants from Stage 079.`
- line 35: `# Stage-61 explicit Pe thresholds.` → `# Stage-078 explicit Pe thresholds.`
- line 65: `# Stage-61 numerical Pe thresholds (96.5285..., 11220.5..., 22.0062...,` → `# Stage-078 numerical Pe thresholds (96.5285..., 11220.5..., 22.0062...,`

Leave line 21 (`banner("STAGE 080 — FAMILY-1 ZETA THRESHOLDS")`) unchanged — it is already correct.

**Verification command:**
`grep -n "Stage.6[12]\|Stage 6[12]\|Stage-6[12]" scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py` returns no matches, and `python3 scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py` exits 0 with the same printed values as the committed transcript.
