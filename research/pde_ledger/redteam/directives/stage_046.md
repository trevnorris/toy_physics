---
unit_id: 046
batch: III.1
created_at: 2026-05-26T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 046

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

Do NOT introduce new features, refactors, or stylistic changes.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side (notes — authoritative derivation source per audit protocol):**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:87-91` quote (P_R):
  ```
  P_R
  = 4 R^4 xi^3
    + 54 R^2 delta^2 xi + 90 R^2 delta xi^2 + 36 R^2 xi^3
    + 230 R delta^3 + 324 R delta^2 xi + 230 R delta xi^2
    + 81 delta^3 + 243 delta^2 xi + 243 delta xi^2 + 81 xi^3.
  ```
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:131-133` quote (P_1):
  ```
  P_1
  = 18 R^2 delta^2 xi + 36 R^2 delta xi^2 + 22 R^2 xi^3
    + 81 R delta^3 + 248 R delta^2 xi + 99 R delta xi^2
    + 230 delta^3 + 423 delta^2 xi + 360 delta xi^2 + 99 xi^3,
  ```
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:136-139` quote (P_2 fragment): `... + 237 R^2 xi^4 + ...`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:67-79` (P_R) uses `162 R delta^3` and `162 R delta xi^2`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:98-109` (P_1) uses `180 R delta^2 xi` and `162 delta^3`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:111-128` (P_2) uses `220 R^2 xi^4`.
- Mathematica `D[fTr, r]` and `Factor[fFlat - fTr]` (`/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.txt:14, 20`) independently confirm the script-side values `162`, `162`, `180`, `162`, `220`.

## Resolve before fix_loop

The notes give five coefficient values (230, 230 in `P_R`; 248, 230 in `P_1`; 237 in `P_2`) that disagree with what both engines independently compute (162, 162, 180, 162, 220 respectively). The paper card itself does not quote these coefficients, so only the notes are affected. Which side is correct?

Possible directions (the user picks one):
- (a) **Scripts are correct; notes are arithmetic typos.** This is the math-evident direction: both engines agree on `162/162/180/162/220` from the same definition of `F_tr` the notes use. Action: edit the notes at the five quoted lines to:
  - P_R: replace `230 R delta^3` with `162 R delta^3` and `230 R delta xi^2` with `162 R delta xi^2`.
  - P_1: replace `248 R delta^2 xi` with `180 R delta^2 xi` and `230 delta^3` with `162 delta^3`.
  - P_2: replace `237 R^2 xi^4` with `220 R^2 xi^4`.
  No script change; no re-run required.
- (b) **Notes are correct; scripts use a different `F_tr` definition than the notes intend.** Action: revisit `F_tr` at `scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:44-48` and `mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl:40-44`. Identify the alternative form that would give 230/230/248/230/237 and re-run both engines. (This is unlikely — both engines agree on the closed-form factorization given the `F_tr` definition both share.)
- (c) Both sides have independent typos that happen to differ — flag for deeper hand-derivation review of `dF_tr/dR` and `F_flat - F_tr`.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. If the user picks (a), the resolution edits a notes file only — outside Codex's normal red-team scope — and so a follow-up directive must explicitly authorize the notes-side edit. If the user picks (b), a follow-up directive will specify the script-side edit and the verifier will re-run both engines.

## Applied: F1

files_changed: notes/stages/moving_throat_pde_stage046_tracking_branch_bounds.md:90,132-133,137
summary: Fixed 5 coefficient typos in notes auxiliary polynomials P_R (230→162 twice), P_1 (248→180, 230→162), P_2 (237→220). Scripts unchanged (already correct per both engines' independent derivation). Per user-approved Q4 (a) in batch III.1 v2.
deviation: none — notes-only fix; no script re-run needed.
