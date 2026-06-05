---
unit_id: 036
batch: II.1
created_at: 2026-06-04T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-04T23:06:07-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 036

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts and their saved outputs.

## F1 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt`

**Issue:** Both committed transcripts predate the current scripts. The scripts' main `banner(...)` headers were already fixed to canonical `STAGE 36.x` / `STAGE 036.x` (sympy L63,77,88,157; wl L31,98,156), but the committed outputs still show the old `STAGE 19.x` labels (captured pre-fix). **In addition, two residual stale "Stage 19" self-labels remain in the SymPy source** that the banner-relabel pass missed — the module docstring and the closing print-summary (the latter prints INTO the transcript, so a plain re-run alone would still emit "All Stage 19 checks passed."). The math content of every assertion is unchanged (all `= 0` / PASS); these are label-only.

**Required change** — (a) fix the two residual SymPy self-labels (label-only; match the `.py`'s own 2-digit `STAGE 36.x` banner format; change only the number `19` → `36`; touch no equation/value/assertion):
- `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:3` — `Moving-throat PDE Stage 19 SymPy audit.` → `Stage 36`.
- `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:163` — `print("All Stage 19 checks passed.")` → `Stage 36`.

(b) Make NO change to `.wl` logic or to any math; then re-run both engines to regenerate the saved outputs:
- `python3 .../scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py > .../scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`
- `math -script .../mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl > .../mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt`

If you judge either string load-bearing for a downstream parser, append `## Blocked: F1` with the concern instead of editing that string.

**Verification command:**
After edit + regeneration the SymPy output's first banner reads `STAGE 36.1 — EXACT SUPPORT-FEASIBILITY FUNCTION`, its closing line reads `All Stage 36 checks passed.`, and the Mathematica output's first banner reads `STAGE 036 — SUPPORT-FEASIBILITY FRONTIER`; all residual lines remain `= 0` and all `PASS:` lines remain. Both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt`
- summary: Updated the two residual SymPy Stage 19 self-labels to Stage 36 and regenerated both audit transcripts from the current scripts.
- deviation: none
