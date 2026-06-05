---
unit_id: 074
batch: III.4
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T21:38:08Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 074

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Both findings are low-severity script/output-band self-label drift (remnants of the EM-extension renumber). No assertion, constant, or asserted value changes. Do NOT touch any math, any `expect_zero`/`expectZero` call, any symbol declaration, paper.tex, or notes/.

**ORCHESTRATOR SCOPE (seat policy): this is a `.py`-only session.** Edit ONLY the SymPy `.py` (F2 docstring self-label) and RUN ONLY `python3 <sympy path>` to refresh the SymPy transcript. Do NOT run `math -script` and do NOT edit the `.wl` — the `.wl` source carries no stale self-label (its banner already reads `STAGE 074`); the orchestrator's independent `exec-mathematica` re-run refreshes the stale `.wl` committed transcript afterward. (F1 below names both transcripts for completeness, but Codex only refreshes the SymPy one.)

## F2 — paper_misalignment (self-label only; mechanical, no user resolution needed)

**Subtype:** notes_contradicts_script (stale self-label, not a value disagreement)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:3`

**Issue:** The module docstring's first line carries a stale stage number `moving_throat_pde_stage57_family1_healing_lock_sympy_audit.py`. The file's actual name and the docstring's own next line (`SymPy audit for Stage 074:`) are 074. This is an unambiguous self-label; correct it mechanically.

**Required change:**
At `..._sympy_audit.py:3`:
- before: `moving_throat_pde_stage57_family1_healing_lock_sympy_audit.py`
- after:  `moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`

Make no other change to the file. Do not alter the banner (line 27 already reads `STAGE 074`), the assertions, or the symbol setup.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`
- summary: Corrected the SymPy audit module docstring self-label from stage57 to stage074.
- deviation: none

## F1 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage074_family1_healing_lock_mathematica_audit.txt`

**Issue:** Both committed transcripts predate their scripts and still carry the pre-renumber banner (`STAGE 57` in the SymPy `.txt`, `STAGE 057` in the Mathematica `.txt`). The current scripts emit `STAGE 074`. All numeric/symbolic result lines already match; only the banner is stale.

**Required change:**
No source edit beyond F2. **Codex re-runs ONLY the SymPy script** to regenerate its committed transcript:
- `python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`

Capture stdout to `output/...sympy_audit.txt`. The Mathematica committed transcript (stale `STAGE 057` banner) is refreshed by the orchestrator's independent `exec-mathematica` re-run — Codex does NOT run `math -script` (seat policy).

**Verification command:**
After Codex applies, the verifier re-runs both engines and confirms: SymPy output line 3 and Mathematica output line 3 both read `STAGE 074 — HEALING-LENGTH LOCK AND SUPPORT SCALE`; result lines unchanged (`chi_s = 37/2`, `kappa = 12321/5`, `alpha (numeric) = 49.640709100495331260`); both scripts exit 0 with all PASS lines present; `grep -n "stage57\|STAGE 57\|STAGE 057" scripts/...sympy_audit.py scripts/output/...sympy_audit.txt mathematica/output/...mathematica_audit.txt` returns nothing.

## Applied: F1

- files_changed:
  - `scripts/output/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.txt`
- summary: Regenerated the SymPy audit transcript from the current Stage 074 Python script.
- deviation: none
