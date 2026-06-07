---
unit_id: 135
batch: IV.4
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T23:14:23Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 135

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line range named.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl:78`

**Issue:** Line 78 asserts `expectApprox["closure residual", residual, 0, 10^-14];`, but `residual` (line 64) is `N[piStar - sigmaStar*(4 - sStar), 30]` and `sigmaStar` (line 60) was solved from exactly the equation `piStar == sigmaVar*(4 - sStar)`, i.e. `sigmaStar = piStar/(4 - sStar)`. Therefore `piStar - sigmaStar*(4 - sStar) ≡ 0` by construction, independent of the kernel value or any physics — the assertion can never fail and verifies nothing. The SymPy counterpart already recognized this and demoted it from an assert to a bare print (sympy:84-86, comment "was the original tautological closure residual"). The Mathematica side is out of parity and emits a spurious `PASS: closure residual` line.

**Required change:**
Delete the single assertion line at `.wl:78`:
```
expectApprox["closure residual", residual, 0, 10^-14];
```
Do NOT replace it with another assertion. Leave the residual PRINT at line 72 (`Print["Pi_star - Sigma_star*(4 - S_star) = ", residual];`) exactly as-is so the value is still reported as a sanity probe, matching the SymPy script (sympy:85-86). The kernel value is already independently exercised by the numeric anchors at lines 74-77 and the range checks at lines 67 and 79; no compensating check is needed.

Before/after (only line 78):
- before: `expectApprox["closure residual", residual, 0, 10^-14];`
- after:  (line removed)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 135` and confirm: (a) `.wl` no longer contains an `expectApprox`/`expect*` assertion whose name is "closure residual"; (b) the saved Mathematica transcript no longer contains the `PASS: closure residual` line; (c) the script still `Exit[0]` with every remaining `expect*` check PASSing; (d) the print at line 72 (`Pi_star - Sigma_star*(4 - S_star) = 0...`) is still present.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl`
- summary: Removed the tautological closure residual assertion while leaving the residual print intact.
- deviation: none
