---
unit_id: 072
batch: III.3
created_at: 2026-05-26T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T13:08:49-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 072

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, normally the orchestrator holds for user resolution. F1 here is `paper_misalignment` subtype `notes_contradicts_script` but with a clearly unambiguous direction (the paper card / filename / notes / appendix row all label this unit "Stage 072"; the script's internal strings are the stale artifact of a prior renumbering). The user is not required to resolve a question of fact — only to authorize the labelling sync. Codex is authorized to apply F1 as a script-side string replacement; do NOT touch paper.tex, notes/, or any prose document.

If F1's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_072.tex:1` quote: `\section[Stage 072]{Stage 072: Explicit Branch Placement Map and Threshold Surfaces}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage072_explicit_branch_thresholds.md` (filename uses "stage072"; body header uses "Stage 072")
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex:122` quote: `072 & Explicit branch threshold surfaces & \StatusExactClosure{} & Fail/succeed surfaces in \((\chi_s,\Lambda_\ell,\Upsilon_w)\). \\`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py:3` quote: `moving_throat_pde_stage55_explicit_branch_thresholds_sympy_audit.py`
- `:5` quote: `SymPy audit for Stage 55:`
- `:24` quote: `banner("STAGE 55 — EXPLICIT BRANCH THRESHOLD SURFACES")`
- `:110` quote: `banner("STAGE 55 THEOREM LEDGER")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl:26` quote: `banner["STAGE 055 — EXPLICIT BRANCH THRESHOLD SURFACES"];`

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py:3,5,24,110`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl:26`

**Issue:**
Script docstrings and banner strings still reference the pre-renumbering identifier "Stage 55" / "Stage 055". The filename, paper card, notes file, and appendix row all use "Stage 072". The discrepancy is purely cosmetic but bleeds into saved transcripts (`STAGE 55 — EXPLICIT BRANCH THRESHOLD SURFACES`, `STAGE 55 THEOREM LEDGER`, `STAGE 055 — EXPLICIT BRANCH THRESHOLD SURFACES`).

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py`:

- Line 3: replace `moving_throat_pde_stage55_explicit_branch_thresholds_sympy_audit.py` with `moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py`.
- Line 5: replace `SymPy audit for Stage 55:` with `SymPy audit for Stage 072:`.
- Line 24: replace the banner argument string `STAGE 55 — EXPLICIT BRANCH THRESHOLD SURFACES` with `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES`.
- Line 110: replace the banner argument string `STAGE 55 THEOREM LEDGER` with `STAGE 072 THEOREM LEDGER`.

Do not modify any other content of the SymPy script. Do not change comments, code, variable names, or assertions.

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl`:

- Line 26: replace the banner argument string `STAGE 055 — EXPLICIT BRANCH THRESHOLD SURFACES` with `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES`.

Do not modify any other content of the Mathematica script.

**Self-test note:**
- Edit is pure string replacement; no assertion or code-flow change. The eight `expect_zero`/`expectZero` checks continue to operate on the same symbolic expressions and continue to evaluate to 0.
- Verified that the new strings match the paper card identifier "Stage 072" and the filename, so no new `paper_misalignment` is introduced.
- No variable-independence trap (no diff added).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 072` and `redteam exec-mathematica 072` and confirm:
1. The saved `.txt` transcripts' headers now read `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES` and `STAGE 072 THEOREM LEDGER` (SymPy) and `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES` (Mathematica).
2. The SymPy script's module docstring now references the correct filename `…stage072…_sympy_audit.py` and the correct claim `SymPy audit for Stage 072:`.
3. `grep -in "stage 55\|stage 055\|stage55\|stage055" scripts/moving_throat_pde_stage072_*.py mathematica/moving_throat_pde_stage072_*.wl` returns no matches.
4. Both scripts exit 0, and all eight pre-existing PASS lines (`shell fail asymptotic`, `shell suff asymptotic`, `compression fail asymptotic`, `compression suff asymptotic`, plus the four `Delta_0`/`Delta_inf` shell/comp leading-order ratio anchors) remain.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl`
- summary: Updated the stale Stage 55/055 docstring and banner labels to Stage 072 in the SymPy and Mathematica audit scripts.
- deviation: none
