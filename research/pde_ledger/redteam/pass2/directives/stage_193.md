---
unit_id: 193
batch: V.3
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 193

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file ranges named.

After editing, RUN the affected scripts and iterate until they exit 0 with all in-file checks passing. Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.txt`

**Issue:** The committed SymPy output `.txt` predates the current `.py` (output mtime 2026-06-01 11:43:38; script mtime 2026-06-03 15:59:11). The stale transcript still carries the old `STAGE 176` banner (lines 3, 127) whereas the current script emits `STAGE 193` (`.py` L49, L164). The math content is correct and unchanged; only the captured banner/ledger label is stale.

**Required change:**
No source edit. Re-run the SymPy script and overwrite its saved output:
`python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.py`
and capture stdout to
`/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage193_isotropic_grouped_p2_target_surface_sympy_audit.txt`.
(The Mathematica `.wl`/`.txt` pair is already fresh — `.txt` newer than `.wl` — so no Mathematica action is required, but a re-run is harmless.)

**Verification command:**
The verifier runs `redteam exec-sympy 193`; the regenerated `.txt` banner must read `STAGE 193 ...` and `STAGE 193 LEDGER`, every `expect_zero` line must print `= 0`, and the script must exit 0.
