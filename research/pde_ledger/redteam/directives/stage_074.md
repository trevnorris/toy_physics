---
unit_id: 074
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-23T05:10:54Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 074

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:29-40`

**Issue:**
The SymPy script defines `chi_lock = sp.simplify(Lambda_ell / 2)` (line 33) and then asserts `chi_lock - Lambda_ell/2 == 0` (line 39). This is tautological: the subtraction equals zero by construction, regardless of whether the GNLS healing-length relation `ell = hbar/(2 m c_s)` actually produces `chi_s = Lambda_ell/2`. The Mathematica counterpart (lines 33-37 of the `.wl`) derives `chi_s` from `m c_s L / hbar` with the substitutions `c_s -> hbar/(2 m ell)` then `L/ell -> Lambda_ell`. The SymPy script must mirror that derivation chain so the assertion becomes non-tautological.

**Required change:**

Edit the block from line 29 through line 40 of `scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`.

Replace the current block:
```python
Lambda_ell = sp.symbols("Lambda_ell", positive=True)
chi_s = sp.symbols("chi_s", positive=True)
kappa = sp.symbols("kappa", positive=True)

chi_lock = sp.simplify(Lambda_ell / 2)
kappa_lock = sp.simplify(4 * chi_lock**2 + sp.Rational(4, 5) * Lambda_ell**2)

print("chi_s (locked) =", chi_lock)
print("kappa(Lambda_ell) =", kappa_lock)

expect_zero("chi_s - Lambda_ell/2", chi_lock - Lambda_ell / 2)
expect_zero("kappa - (9/5) Lambda_ell^2", kappa_lock - sp.Rational(9, 5) * Lambda_ell**2)
```

with:
```python
Lambda_ell = sp.symbols("Lambda_ell", positive=True)
hbar, m_psi, c_s, ell, L = sp.symbols("hbar m_psi c_s ell L", positive=True)

# Physical definition of the dimensionless support scale chi_s = m c_s L / hbar.
chi_def = m_psi * c_s * L / sp.symbols("hbar", positive=True)
chi_def = m_psi * c_s * L / hbar

# Apply the GNLS healing/compliance width: ell = hbar / (2 m c_s),
# equivalently c_s = hbar / (2 m_psi ell).
chi_in_ell = sp.simplify(chi_def.subs(c_s, hbar / (2 * m_psi * ell)))
print("chi (after healing-length substitution) =", chi_in_ell)

# Re-express L/ell as the dimensionless ratio Lambda_ell.
chi_lock = sp.simplify(chi_in_ell.subs(L, Lambda_ell * ell))

# Family-1 branch coefficient: kappa = 4 chi_s^2 + (4/5) Lambda_ell^2.
# Coefficients 4 and 4/5 come from the Family-1 Euler-Lagrange branch
# (carried forward from the earlier Family-1 stages); this stage only
# verifies that, with chi_s locked to Lambda_ell/2, kappa reduces to
# (9/5) Lambda_ell^2.
kappa_lock = sp.simplify(4 * chi_lock**2 + sp.Rational(4, 5) * Lambda_ell**2)

print("chi_s (locked) =", chi_lock)
print("kappa(Lambda_ell) =", kappa_lock)

expect_zero("chi_s - Lambda_ell/2", chi_lock - Lambda_ell / 2)
expect_zero("kappa - (9/5) Lambda_ell^2", kappa_lock - sp.Rational(9, 5) * Lambda_ell**2)
```

Remove the redundant first assignment of `chi_def` (the line that re-uses `sp.symbols("hbar", positive=True)`) — the two-line pair above is shown for clarity; use only the second `chi_def = m_psi * c_s * L / hbar` line and delete the first. After the edit, `chi_def` should be defined exactly once, immediately after the symbol declarations.

Do not touch the reference-branch block at lines 42-55 or the final ledger print at lines 57-59. Those remain valid once `chi_lock` is derived non-tautologically.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 074` and confirm:
- The script source contains the symbols `hbar`, `m_psi`, `c_s`, `ell`, `L`.
- The intermediate print `chi (after healing-length substitution) = L/(2*ell)` appears in the output.
- The assertion `chi_s - Lambda_ell/2 = 0` still passes.
- The script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`
- summary: Replaced the direct chi lock assignment with the healing-length substitution chain deriving `chi_s = Lambda_ell/2`.
- deviation: none

## F2 — insufficient_verification

**Target:** `scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py:33-34`

**Issue:**
The kappa formula `4 chi_s^2 + (4/5) Lambda_ell^2` is introduced on line 34 without comment or provenance. A reader cannot tell from the script alone whether the `4` and `4/5` are Family-1 branch coefficients or an unjustified hard-coded form. Combined with F1, the resulting check then reduces to the arithmetic identity `(1/2)^2 * 4 + 4/5 = 9/5`, which tests no physics. Adding the substitution chain from F1 fixes the chi_s side; this finding requires an inline comment anchoring the kappa coefficients.

**Required change:**

The replacement block in F1 already includes the required comment immediately above the `kappa_lock = ...` line:

```
# Family-1 branch coefficient: kappa = 4 chi_s^2 + (4/5) Lambda_ell^2.
# Coefficients 4 and 4/5 come from the Family-1 Euler-Lagrange branch
# (carried forward from the earlier Family-1 stages); this stage only
# verifies that, with chi_s locked to Lambda_ell/2, kappa reduces to
# (9/5) Lambda_ell^2.
```

If F1 is applied as specified above, F2 is satisfied by the same edit. If F1 is blocked, apply only the comment block: insert the four-line comment immediately above the existing `kappa_lock = sp.simplify(4 * chi_lock**2 + sp.Rational(4, 5) * Lambda_ell**2)` line (currently line 34) and make no other change.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 074` and confirm:
- The four-line comment beginning `# Family-1 branch coefficient` appears immediately above the `kappa_lock = ...` definition.
- The assertion `kappa - (9/5) Lambda_ell^2 = 0` still passes.
- The script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage074_family1_healing_lock_sympy_audit.py`
- summary: Added the required inline provenance comment for the Family-1 kappa coefficients above the `kappa_lock` definition.
- deviation: none
