---
unit_id: 168
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 168

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (cosmetic banner mislabel)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl:26`

**Issue:** Both scripts open with a banner that prints the wrong unit number, `STAGE 151`, even though this is unit 168 (filename, sympy docstring line 3, and Mathematica closing line 124 all say 168). The wrong label is propagated into both saved transcripts (line 11 of each `.txt`). It is a label-only defect — no assertion, symbol, or numeric value is affected.

**Required change:**
1. In the SymPy script, line 30, change:
   - before: `banner("STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION")`
   - after:  `banner("STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION")`
2. In the Mathematica script, line 26, change:
   - before: `banner["STAGE 151 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION"];`
   - after:  `banner["STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION"];`

Change only the literal `151` to `168` inside each banner string. Do not touch any other line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 168` and `redteam exec-mathematica 168` and confirm: (a) line 11 of both saved `.txt` transcripts now reads `STAGE 168 — EXACT OFF-BUNDLE SLIPPAGE DECOMPOSITION`, and (b) both scripts still exit 0 with all five zero-residual / `PASS` lines unchanged.

## Applied: F1

files_changed:
- scripts/moving_throat_pde_stage168_off_bundle_slippage_sympy_audit.py
- mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl

summary: Changed the literal `151` to `168` inside the opening banner string of both the SymPy and Mathematica audit scripts.

deviation: none
