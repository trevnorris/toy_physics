---
unit_id: 156
batch: IV.6
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 156

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py:126`

**Issue:**
The paper card and notes commit six load-bearing numerical outputs for this stage: \(\Sigma_0^{\rm can}\approx 4.651033550168867\), \(\widehat T_{m,\rm can}\approx 1.4467083664567613\), \(\mathfrak g_{\rm can}\approx 0.758035078944663\), \(\mathcal R_{\rm can}=1/4\), \(\mathcal S_{\rm can}\approx 0.6703621156734617\), \(\Pi_{\rm can}\approx 3.871564377479002\). The SymPy script computes and prints all six (lines 103–116) but its existing `assert` block (lines 123–126) only checks `|g_can - g_star| <= 1e-10` and `|R_can - 0.25| <= 1e-10`. Those two checks are not independent verifications of the paper claim: the first is the bisection's own drive target (`bisect` at line 103 was specifically constructed to make `g_fp(Sigma0_can) = g_star`), and the second is algebraically forced once the first passes — with the hardcoded `rF1 = 1.77799353547498` and `g_star = 0.758035078944663`, `R = (g_star - rF1)^2/(1 + rF1^2) ≈ 0.2499999...` independent of the rest of the solve. The Mathematica counterpart asserts all six via `expectApprox` (lines 147–152). Add the four missing asserts to the SymPy side so both engines guard the same load-bearing claim.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py`, replace lines 123–126 (the existing two-assert block) with the following six-assert block. Lines 117–122 and lines 127–131 are unchanged.

Before (lines 123–126):
```python
if abs(g_can - g_star) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore g = g_*.")
if abs(R_can - 0.25) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore R = 1/4.")
```

After (lines 123–134, replacing the four-line block above):
```python
if abs(g_can - g_star) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore g = g_*.")
if abs(R_can - 0.25) > 1e-10:
    raise AssertionError("Renormalized canonical branch did not restore R = 1/4.")
if abs(Sigma0_can - 4.651033550168867) > 1e-10:
    raise AssertionError("Sigma0_can deviates from notes value 4.651033550168867.")
if abs(S_can - 0.6703621156734617) > 1e-11:
    raise AssertionError("S_can deviates from notes value 0.6703621156734617.")
if abs(Pi_can - 3.871564377479002) > 1e-10:
    raise AssertionError("Pi_can deviates from notes value 3.871564377479002.")
if abs(T_hat_can - 1.4467083664567613) > 1e-10:
    raise AssertionError("T_hat_can deviates from notes value 1.4467083664567613.")
```

Rationale for tolerances: the SymPy transcript shows printed agreement to 12+ digits with the notes values, so 1e-10 / 1e-11 are comfortable margins; they are one decade looser than the Mathematica `expectApprox` tolerances (1e-11 / 1e-12 in lines 147–152 of the `.wl`) to account for the SymPy script's absence of the `phShift = ph - Min[ph]` overflow guard that the Mathematica script uses. This keeps both engines green at current numerical agreement while still flagging any future drift exceeding five significant figures.

Do not edit any other lines. Do not add comments. Do not touch the Mathematica script.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 156`. The script must (a) still exit 0, (b) still emit the existing printed lines for `Sigma0_can`, `g_can`, `S_can`, `R_can`, `Pi_can`, `T_hat_can`, and (c) the saved transcript's `# Exit code: 0` / `# Status: PASS` header must remain. The four new assert lines are silent on success (they only raise on failure); the verifier will inspect the source file to confirm the new asserts are present at the indicated locations, and may optionally run a perturbation check (e.g., change `N` to `241` and confirm one of the new asserts now fails) to confirm the new asserts are not trivially-true.
