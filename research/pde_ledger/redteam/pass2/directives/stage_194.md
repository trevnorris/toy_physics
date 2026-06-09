---
unit_id: 194
batch: V.3
created_at: 2026-06-09T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 194

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — paper_misalignment (subtype: paper_missing_script_claim — card-text lag)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_194.tex:11` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_mathematica_audit.wl:1` — a complete Mathematica audit exists and passes (saved output all `PASS`, banner `STAGE 194`).

### Resolve before fix_loop

The card says "Mathematica audit: none yet," but a passing `.wl` audit for this stage exists. Should the card be updated to reference the existing `.wl`?

Possible directions (the user picks one):
- (a) Card text is stale → update card line 11 to name the existing `.wl` audit (paper-side prose edit; not Codex/script — follow-up directive after user OK).
- (b) The `.wl` is provisional/not yet blessed → leave the card as-is until the second engine is formally accepted.

The orchestrator will not invoke Codex on F1 until the user has chosen a direction. (No script change in either direction.)

## F2 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.txt`

**Issue:** The committed SymPy output (`.txt`, mtime 2026-06-01 11:43:38) is older than its script (`.py`, mtime 2026-06-03 15:59:11) and its banner still reads `STAGE 177` (lines 3, 143: `STAGE 177 — …`, `STAGE 177 LEDGER`), whereas the current `.py` prints `STAGE 194` (lines 35, 208). The result values in the stale `.txt` are otherwise unchanged from the current source; only the banner is out of date. No source edit is required — the fix is to regenerate the output.

**Required change:**
- Do NOT edit the `.py` source.
- Re-run the SymPy script and recommit its refreshed output:
  `python3 scripts/moving_throat_pde_stage194_outgoing_l2_fingerprint_and_deformation_algebra_sympy_audit.py` (redirecting to the committed `.txt` path per the project's normal output-capture convention).
- Confirm the script exits 0 with all `expect_zero` lines printing `= 0`.

**Verification command:**
After Codex re-runs, the verifier confirms the refreshed `.txt` line 3 reads `STAGE 194 — …` (not `177`), its `STAGE … LEDGER` banner reads `STAGE 194 LEDGER`, the `.txt` mtime is ≥ the `.py` mtime, and the script exited 0 with every check at 0.
