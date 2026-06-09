---
unit_id: 192
batch: V.3
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 192

Only one finding, and it is `stale_output`. No source edit is required — the SymPy
`.py` is mathematically current; only its committed `.txt` transcript is stale
(carries the pre-renumber `STAGE 175` banner from 2026-05-11, predating the
2026-06-03 script). Codex must NOT edit the `.py` or `.wl`.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.txt`

**Issue:** The saved SymPy output (mtime 2026-05-11 12:48) predates the current
script (mtime 2026-06-03 15:59) and still prints `STAGE 175 — …` (txt line 11) and
`STAGE 175 LEDGER` (txt line 634), whereas the current script prints `STAGE 192`
(`scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py:35`)
and `STAGE 192 LEDGER` (py:196). The math content is unchanged; only the stage
banner differs. The Mathematica output is already fresh.

**Required change:** None to source. Re-run the SymPy script
(`python3 scripts/moving_throat_pde_stage192_orbit_quotient_projectors_sympy_audit.py`,
exit 0, all `expect_zero` residuals zero) so the committed output is refreshed with
the `STAGE 192` banner. The orchestrator's independent re-run already accomplishes
this.

**Verification command:** After the re-run, the SymPy output `.txt` line 11 reads
`STAGE 192 — EXACT ORBIT/QUOTIENT PROJECTORS …` and the ledger header reads
`STAGE 192 LEDGER`. No assertion content changes.
