---
unit_id: 034
batch: II.1
created_at: 2026-06-04T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T05:07:59Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 034

Apply the finding below. After applying, append an `## Applied: F<n>` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Edit exactly the file:line strings named. These are LABEL-ONLY edits to comment/print strings — do NOT touch any equation, value, variable, assertion, or `expect_zero`/`expectZero` target.

After editing, RUN the affected script (`python3 <py>`) and confirm it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (residual self-labels + stale transcripts)

**Context:** The committed `.txt` transcripts predate the June-3 numbering commit (`e2a4780`), which fixed the scripts' main `banner(...)` headers to the canonical `STAGE 34.x` / `STAGE 034` but left residual stale "Stage 17" self-labels in the SymPy module docstring and the closing print-summary. The print-summary prints INTO the transcript. The Mathematica `.wl` carries no residual self-label. All math is unchanged (every residual `= 0` / `PASS`); these are label-only.

**Required change** — replace the stale stage number `17` with the canonical stage number, formatted to match the `.py`'s OWN main banner (2-digit `STAGE 34.x`). Change nothing but the number:
- `scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py:3` — `Moving-throat PDE Stage 17 SymPy audit.` → `Stage 34`.
- `scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py:94` — `print("All Stage 17 checks passed.")` → `Stage 34`.

If you judge any of these strings load-bearing for a downstream parser, append `## Blocked: F1` with the concern instead of editing.

**Verification:** After edit, `python3 <py>` exits 0; the regenerated SymPy transcript closing line reads `All Stage 34 checks passed.`; every residual line remains `0`. (The Mathematica output staleness resolves on the orchestrator's standard re-run; no `.wl` edit needed.)

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`
- summary: Replaced the stale Stage 17 labels in the SymPy audit docstring and closing print with Stage 34.
- deviation: none
