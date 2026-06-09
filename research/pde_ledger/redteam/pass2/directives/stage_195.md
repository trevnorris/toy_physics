---
unit_id: 195
batch: V.3
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 195

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts and their saved outputs.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.txt`

**Issue:** The committed SymPy saved output is stale. Its mtime (2026-06-01 11:53:34) predates the current `.py` (2026-06-03 15:59:11), and its content reflects an earlier revision: the banner reads `STAGE 178 — EXACT SOURCE-MAP REDUCTION…` (line 3) and `STAGE 178 LEDGER` (line 122), whereas the current `scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py` prints `STAGE 195 …` (line 35) and `STAGE 195 LEDGER` (line 189). No source edit is required — the `.py` math is correct (all `expect_zero` checks reduce to 0); only the committed transcript is out of date.

**Required change:**
1. Do NOT edit the `.py` source.
2. Re-run the SymPy script and overwrite the committed saved output:
   `python3 scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py` and capture stdout to `scripts/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.txt` (matching the repo's existing capture convention for this stage).
3. Confirm the refreshed `.txt` shows the `STAGE 195` banner (not `STAGE 178`) and every `expect_zero(...)` line ends in `= 0`.

**Verification command:**
After Codex applies, the verifier runs the SymPy engine for unit 195 and confirms: the saved `.txt` mtime ≥ the `.py` mtime; line 3 reads `STAGE 195 — …`; the ledger banner reads `STAGE 195 LEDGER`; and the script exits 0 with all in-file checks `= 0`.
