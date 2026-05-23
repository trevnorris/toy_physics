---
unit_id: 053
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
applied_at: 2026-05-22T23:36:02Z
findings_applied: 1
findings_blocked: 0
---

# Codex directive — unit 053

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py:72-76`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl:63-72`

**Issue:**
The "small-alpha linear coefficient" check in both engines never reads the actual series expansion of `Omega_alpha`. In SymPy, `expected_linear = sp.simplify(2 / pi - sp.Rational(1, 2))` is then asserted to equal `(4 - pi) / (2 * pi)` — a pure algebraic identity that holds regardless of the integration. In Mathematica, `linearCoeff = FullSimplify[2/Pi - 1/2]` plays the same role. The series object (`series_small` / `seriesSmall`) is computed and printed but never enters the assertion. If the integrand or closed form were wrong, this check would still pass. The fix is to extract the linear coefficient from the series itself and assert it equals `(4-pi)/(2pi)`.

**Required change:**

SymPy file `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py`, lines 72-76. The current block is:

```python
series_small = sp.series(Omega_alpha, alpha, 0, 3).removeO()
print("small-alpha series =", series_small)
expected_linear = sp.simplify(2 / pi - sp.Rational(1, 2))
print("linear coefficient =", expected_linear)
expect_zero("linear coefficient - (4-pi)/(2pi)", expected_linear - (4 - pi) / (2 * pi))
```

Replace it with:

```python
series_small = sp.series(Omega_alpha, alpha, 0, 3).removeO()
print("small-alpha series =", series_small)
linear_coeff = sp.simplify(series_small.coeff(alpha, 1))
print("linear coefficient =", linear_coeff)
expect_zero("linear coefficient - (4-pi)/(2pi)", linear_coeff - (4 - pi) / (2 * pi))
```

Mathematica file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl`, line 64. The current line is:

```mathematica
linearCoeff = FullSimplify[2/Pi - 1/2];
```

Replace it with:

```mathematica
linearCoeff = FullSimplify[Coefficient[seriesSmall, alpha, 1], Assumptions -> ell > 0];
```

Lines 63 (`seriesSmall = ...`), 69 (`Print["linear coefficient = ", fmt[linearCoeff]];`), and 72 (`expectZero["linear coefficient - (4-Pi)/(2Pi)", linearCoeff - (4 - Pi)/(2 Pi)];`) are unchanged. The other assertion lines on 70-71 are unrelated and must not be touched.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 053` and `redteam exec-mathematica 053` and confirm:
1. The new lines `linear_coeff = sp.simplify(series_small.coeff(alpha, 1))` and `Coefficient[seriesSmall, alpha, 1]` appear in their respective files.
2. The literals `expected_linear = sp.simplify(2 / pi - sp.Rational(1, 2))` and `linearCoeff = FullSimplify[2/Pi - 1/2]` are gone.
3. Both scripts still exit 0 (the actual linear coefficient is `(4-pi)/(2pi)`, so the substantive assertion still passes).

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl`
- summary: Replaced the tautological linear coefficient checks with coefficients extracted from the computed small-alpha series in both engines.
- deviation: none
