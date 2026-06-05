---
unit_id: 033
batch: II.1
created_at: 2026-06-04T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-04T23:04:26-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 033

Apply the finding below. After applying, append an `## Applied: F<n>` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Edit exactly the file:line strings named. These are LABEL-ONLY edits to comment/print strings — do NOT touch any equation, value, variable, assertion, or `expect_zero`/`expectZero` target.

After editing, RUN the affected scripts (`python3 <py>` and `math -script <wl>`) and confirm each exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (residual self-labels + stale transcripts)

**Context:** The committed `.txt` transcripts predate the June-3 numbering commit (`e2a4780`), which fixed the scripts' main `banner(...)` headers to the canonical `STAGE 33.x` / `STAGE 033` but left residual stale "Stage 16" self-labels in the module docstring, the closing print-summary, and a `.wl` comment. The print-summary prints INTO the transcript. All math is unchanged (every residual `= 0` / `PASS`); these are label-only.

**Required change** — replace the stale stage number `16` with the canonical stage number, formatted to match each file's OWN main banner (the `.py` uses 2-digit `STAGE 33.x`; the `.wl` uses its banner's format). Change nothing but the number:
- `scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:3` — `Moving-throat PDE Stage 16 SymPy audit.` → `Stage 33`.
- `scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:133` — `print("All Stage 16 checks passed.")` → `Stage 33`.
- `mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl:134` — the comment `the Stage 16.1 monotonicity identity and the Stage 16.6 gate identity at` → `Stage 33.1` / `Stage 33.6` (sub-stage refs; keep the `.k` suffix, change only the base number).

If you judge any of these strings load-bearing for a downstream parser, append `## Blocked: F1` with the concern instead of editing.

**Verification:** After edit, `python3 <py>` and `math -script <wl>` exit 0; the regenerated transcript closing line reads `All Stage 33 checks passed.`; every residual line remains `0` / `PASS`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl`
  - `scripts/output/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.txt`
- summary: Replaced stale Stage 16 labels with the canonical Stage 33/033 labels and regenerated the affected transcripts.
- deviation: none
