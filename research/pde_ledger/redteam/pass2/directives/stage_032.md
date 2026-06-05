---
unit_id: 032
batch: II.1
created_at: 2026-06-04T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-04T22:52:20-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 032

Apply the finding below. After applying, append an `## Applied: F<n>` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Edit exactly the file:line strings named. These are LABEL-ONLY edits to comment/print strings — do NOT touch any equation, value, variable, assertion, or `expect_zero`/`expectZero` target.

After editing, RUN the affected scripts (`python3 <py>` and `math -script <wl>`) and confirm each exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (residual self-labels + stale transcripts)

**Context:** The committed `.txt` transcripts predate the June-3 numbering commit (`e2a4780`), which fixed the scripts' main `banner(...)` headers to the canonical `STAGE 32.x` / `STAGE 032` but left residual stale "Stage 15" self-labels in the module docstring and the closing print-summary. The print-summary prints INTO the transcript, so a plain re-run would still emit "All Stage 15 checks passed." All math is unchanged (every residual `= 0` / `PASS`); these are label-only.

**Required change** — replace the stale stage number `15` with the canonical stage number, formatted to match each file's OWN main banner (the `.py` uses 2-digit `STAGE 32.x`; the `.wl` uses its banner's format). Change nothing but the number:
- `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:3` — `Moving-throat PDE Stage 15 SymPy audit.` → canonical (`Stage 32`).
- `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py:217` — `print("All Stage 15 checks passed.")` → `Stage 32`.
- `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl:198` — `Print["All Stage 15 checks passed."];` → canonical (match the `.wl` banner's stage-number format).

If you judge any of these strings load-bearing for a downstream parser, append `## Blocked: F1` with the concern instead of editing.

**Verification:** After edit, `python3 <py>` and `math -script <wl>` exit 0; the regenerated transcripts' closing line reads `All Stage 32 checks passed.` (matching the `STAGE 32.x` banners); every residual line remains `0` / `PASS`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- summary: Replaced the stale Stage 15 self-labels in the SymPy docstring and both closing pass summaries with Stage 32.
- deviation: none
