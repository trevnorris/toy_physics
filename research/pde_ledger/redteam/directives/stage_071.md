---
unit_id: 071
batch: III.3
created_at: 2026-05-22T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-22T20:09:59-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 071

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py:65-69`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl:77-81`

**Issue:**

In both scripts, `K_m` is defined as `T_X / ell` and then `eta = K_m * L / T_X` is checked against `L/ell`. Because `eta = (T_X/ell) * L / T_X = L/ell` is an algebraic identity in the symbol `T_X`, the check passes regardless of the actual closed form of `T_X` and provides no physics verification. The fix pins `K_m` to the closed form derived from `T_X = pi a^2 ell hbar^2 / (3 m rho_w)`, replacing the tautology with a substantive equality that fails if `I_f` or any factor in `T_X` is wrong.

**Required change:**

### SymPy script

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`, replace the existing block at lines 65-69:

Before (lines 65-69):
```python
Km = sp.simplify(Tx / ell)
eta = sp.simplify(Km * L / Tx)
print("K_m =", Km)
print("eta =", eta)
expect_zero("eta - L/ell", eta - L / ell)
```

After:
```python
Km = sp.simplify(Tx / ell)
eta = sp.simplify(Km * L / Tx)
print("K_m =", Km)
print("eta =", eta)
Km_expected = sp.pi * a**2 * hbar**2 / (3 * m * rho_w)
expect_zero("K_m - pi a^2 hbar^2 / (3 m rho_w)", Km - Km_expected)
expect_zero("eta - L/ell (from closed-form K_m)", (Km_expected * L / Tx) - L / ell)
```

This adds two assertions: one pins `K_m` to its closed form (failure modes: missing factor of `ell`, wrong `I_f`, sign error); the other reconstructs `eta` from the closed-form `K_m` and `T_X` (failure modes: any algebraic error in `T_X` or its `ell`-dependence). The original tautological `expect_zero("eta - L/ell", ...)` is removed.

### Mathematica script

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl`, replace the existing block at lines 77-81:

Before (lines 77-81):
```mathematica
km = FullSimplify[tx/ell, Assumptions -> $Assumptions];
eta = FullSimplify[km*L/tx, Assumptions -> $Assumptions];
Print["K_m = ", fmt[km]];
Print["eta = ", fmt[eta]];
expectZero["eta - L/ell", eta - L/ell];
```

After:
```mathematica
km = FullSimplify[tx/ell, Assumptions -> $Assumptions];
eta = FullSimplify[km*L/tx, Assumptions -> $Assumptions];
Print["K_m = ", fmt[km]];
Print["eta = ", fmt[eta]];
kmExpected = Pi*a^2*hbar^2/(3*m*rhoW);
expectZero["K_m - pi a^2 hbar^2 / (3 m rhoW)", km - kmExpected];
expectZero["eta - L/ell (from closed-form K_m)", (kmExpected*L/tx) - L/ell];
```

The original `expectZero["eta - L/ell", eta - L/ell];` is removed and replaced by the two substantive checks.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 071` and `redteam exec-mathematica 071`. Both should exit 0. The new transcripts must contain:

- SymPy: a line `K_m - pi a^2 hbar^2 / (3 m rho_w) = 0` and a line `eta - L/ell (from closed-form K_m) = 0`.
- Mathematica: `PASS: K_m - pi a^2 hbar^2 / (3 m rhoW)` and `PASS: eta - L/ell (from closed-form K_m)`.

The old `eta - L/ell = 0` line (without the `(from closed-form K_m)` qualifier) must NOT appear in either transcript, since it was the tautology being removed.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl`
- summary: Replaced the tautological eta-only checks with closed-form K_m checks and eta reconstruction from closed-form K_m.
- deviation: none
