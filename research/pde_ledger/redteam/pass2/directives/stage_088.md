---
unit_id: 088
batch: III.5
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T00:06:58Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 088

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`python3 scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`) and confirm it exits 0 with the identical transcript (the change is docstring-only, non-functional).

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (stale self-label, numbering)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py:3` and `:5`

**Issue:** The SymPy docstring carries the pre-renumber OWN-number self-label. Line 3 reads `moving_throat_pde_stage71_loading_ratio_from_minimal_module_sympy_audit.py` and line 5 reads `SymPy audit for Stage 71.`, but the file's canonical stage number is 088 (its filename, the paper card `stage:088`, the manifest, and the in-script banner at line 29 all say 088). "Stage 71" / "stage71" is the known +17 EM-extension renumber drift (088 − 17 = 071). This is the file labeling itself; cross-references elsewhere (`stage085` at lines 112 and in the `.wl`) are correctly numbered and out of scope.

**Required change:**
- At line 3, change `moving_throat_pde_stage71_loading_ratio_from_minimal_module_sympy_audit.py` to `moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`.
- At line 5, change `SymPy audit for Stage 71.` to `SymPy audit for Stage 088.`
- Do NOT alter `stage085` references (lines 112 etc.) — those are correct cross-references.

**Verification command:**
After Codex applies, the verifier will confirm `grep -n "stage71\|Stage 71" scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py` returns nothing, that the banner (line 29) and filename still read 088, and that `redteam exec-sympy 088` exits 0 with the unchanged transcript.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py`
- summary: Updated the SymPy audit docstring self-label from stage71/Stage 71 to stage088/Stage 088.
- deviation: none
