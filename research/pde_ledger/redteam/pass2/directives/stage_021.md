---
unit_id: 021
batch: I.2
created_at: 2026-06-04T20:40:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-04T21:15:34-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 021

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt` (regenerate) and, optionally, the leftover label literal at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:195`

**Issue:** The Mathematica saved output (`…_mathematica_audit.txt`, mtime 2026-05-25 22:24:38) predates the `.wl` script (mtime 2026-06-03 15:59:11). The only intervening `.wl` change (commit `e2a4780`) was a banner-label fix `STAGE 004` → `STAGE 021` at line 35. The captured output still shows the stale banner `STAGE 004 — MAXWELL + MIXED-SECTOR REDUCTION` (output line 3) and `Stage 004 Mathematica audit passed.` (output line 82). No math result is affected — every residual in the output is `= 0` / `PASS` and matches the current script. This is informational only.

**Required change:**
1. (Optional, cosmetic) In `…_mathematica_audit.wl:195`, change the leftover stale closing label
   `Print["Stage 004 Mathematica audit passed."];`
   to
   `Print["Stage 021 Mathematica audit passed."];`
   This is a non-load-bearing print string (not part of the earlier line-35 banner fix); no equation, variable, or assertion changes.
2. Re-run `math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl` and recapture stdout into `mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt` so the saved banner reads `STAGE 021` (and the closing line reads `Stage 021 …` if step 1 was applied).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 021` and confirms: the script exits 0, output line 3 reads `STAGE 021 — MAXWELL + MIXED-SECTOR REDUCTION`, every `PASS:` line and `= 0` residual is unchanged from the prior transcript, and (if step 1 applied) the closing line reads `Stage 021 Mathematica audit passed.`

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt`
- summary: Updated the stale closing stage label and regenerated the Mathematica audit transcript from the current script.
- deviation: none
