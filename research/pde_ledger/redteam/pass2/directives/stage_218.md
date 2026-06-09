---
unit_id: 218
batch: VI.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 218

There is exactly one finding (F1, `stale_output`, low severity) and it requires NO Codex edit.

## F1 — stale_output (no Codex action)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.txt`

**Issue:**
The SymPy `.py` (mtime 2026-06-03 15:59:11, commit `e2a4780` "numbering reconciliation Phase 1:
doc-only stage-label fixes") was modified after its `.txt` output was captured (mtime
2026-06-02 12:33:13, commit `a12029e`). The only change was the banner string `STAGE 201`→`STAGE 218`;
the captured `.txt` still prints `STAGE 201` on line 3. Every numeric/logical line in the stale `.txt`
matches what the current script emits — there is NO content disagreement, only the stale banner label.

**Required change:**
None. Do NOT edit the script or the output by hand. The orchestrator's standard independent re-run
(`redteam exec-sympy 218`) regenerates the `.txt` with the corrected `STAGE 218` banner. The
Mathematica output is already fresh.

**Verification command:**
After the orchestrator re-runs `python3` on the SymPy script, confirm `.txt` line 3 reads
`STAGE 218 — FULL SUPPORT-<=5 ...` and all numeric/PASS lines are unchanged.
