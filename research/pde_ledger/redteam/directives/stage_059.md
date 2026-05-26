---
unit_id: 059
batch: III.2
created_at: 2026-05-26
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T11:04:18-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 059

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:72`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:72`

**Issue:**
Both scripts set up positivity-ordered scaffolding `Xi_fail_ordered = Pe_req/(Delta0 + delta_gap)` and `Xi_suff_ordered = Pe_req/Delta0` with `delta_gap > 0`, but never assert the ordering inequality `Xi_fail_ordered <= Xi_suff_ordered` that the paper notes section 4 states explicitly ("Since `Delta_inf >= Delta_0 > 0`, one has `Xi_fail <= Xi_suff`"). The variables are computed and then abandoned. Add the missing positivity check using the engines' existing helpers.

**Required change:**

SymPy: insert one new assertion line directly after the current line 72 (the `zeta_req_branch = ...` line). Before/after of the relevant block:

Before (lines 68-72):
```python
delta_gap = sp.symbols("delta_gap", positive=True, real=True)
DeltaInf_ordered = Delta0 + delta_gap
Xi_fail_ordered = sp.simplify(Pe_req / DeltaInf_ordered)
Xi_suff_ordered = sp.simplify(Pe_req / Delta0)
zeta_req_branch = sp.simplify(A_K * Omega(Pe_req) ** 2)
```

After (lines 68-73):
```python
delta_gap = sp.symbols("delta_gap", positive=True, real=True)
DeltaInf_ordered = Delta0 + delta_gap
Xi_fail_ordered = sp.simplify(Pe_req / DeltaInf_ordered)
Xi_suff_ordered = sp.simplify(Pe_req / Delta0)
zeta_req_branch = sp.simplify(A_K * Omega(Pe_req) ** 2)
expect_positive("Xi_suff - Xi_fail (ordered)", Xi_suff_ordered - Xi_fail_ordered)
```

Mathematica: insert one new assertion line directly after the current line 72 (the `xiSuffOrdered = ...` line). Before/after of the relevant block:

Before (lines 68-72):
```mathematica
Clear[deltaGap];
$Assumptions = $Assumptions && deltaGap > 0;
deltaInfOrdered = FullSimplify[delta0 + deltaGap, Assumptions -> $Assumptions];
xiFailOrdered = FullSimplify[peReq/deltaInfOrdered, Assumptions -> $Assumptions];
xiSuffOrdered = FullSimplify[peReq/delta0, Assumptions -> $Assumptions];
```

After (lines 68-73):
```mathematica
Clear[deltaGap];
$Assumptions = $Assumptions && deltaGap > 0;
deltaInfOrdered = FullSimplify[delta0 + deltaGap, Assumptions -> $Assumptions];
xiFailOrdered = FullSimplify[peReq/deltaInfOrdered, Assumptions -> $Assumptions];
xiSuffOrdered = FullSimplify[peReq/delta0, Assumptions -> $Assumptions];
expectPositive["Xi_suff - Xi_fail (ordered)", xiSuffOrdered - xiFailOrdered];
```

Do not touch any other line. The helpers `expect_positive` (SymPy line 31) and `expectPositive` (Mathematica line 32) are already defined and operate under the positive-real assumption blocks already in scope (`Pe_req > 0`, `Delta0 > 0`, `delta_gap > 0` / `peReq > 0`, `delta0 > 0`, `deltaGap > 0`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 059` and `redteam exec-mathematica 059`. The new SymPy transcript should include a line of the form `Xi_suff - Xi_fail (ordered) = Pe_req*delta_gap/(Delta0*(Delta0 + delta_gap))` (or an algebraically equivalent factored form). The new Mathematica transcript should include `PASS: Xi_suff - Xi_fail (ordered)`. Both scripts must still exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl`
- summary: Added the missing positivity assertions that prove `Xi_suff_ordered - Xi_fail_ordered` is positive under the ordered gap assumptions.
- deviation: none
