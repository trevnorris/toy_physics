---
unit_id: 049
batch: III.2
created_at: 2026-05-22T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
verification_status: pending
applied_at: 2026-05-22T16:54:53-06:00
findings_applied: 2
findings_blocked: 0
---

# Codex directive — unit 049

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py:65-68`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:45`

**Issue:**
In both engines, `k_n` (resp. `kN`) is defined as `simplify((n+1/2)*pi/L)` and then asserted equal to `(n+1/2)*pi/L`. This is `X − X = 0` by construction; the assertion cannot fail and tests no physics. The docstring's first claim ("Exact D/N half-wave momentum") is therefore not actually verified. Replace each tautology with a substantive boundary-condition check: χ_n(s) = √(2/L) sin(k_n s) must satisfy the Neumann condition at s = L, i.e. `cos(k_n · L) = 0` (Dirichlet at s = 0 is automatic from the sine). For integer n this reduces to `cos((n+1/2)π) = 0`, which is non-trivially zero.

**Required change:**

SymPy file `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py`, replace the block at lines 65-68

Before:
```python
    expect_zero(
        "k_n - (n+1/2) pi / L",
        k_n - (n + sp.Rational(1, 2)) * sp.pi / L,
    )
```
After:
```python
    expect_zero(
        "k_n satisfies D/N Neumann boundary",
        sp.cos(k_n * L),
    )
```

Mathematica file `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`, replace line 45

Before:
```wolfram
expectZero["k_n - (n+1/2) pi / L", kN - (n + 1/2) Pi/l];
```
After:
```wolfram
expectZero["k_n satisfies D/N Neumann boundary", Cos[kN l]];
```

Leave all other lines, including the surrounding `Print` statements and the `chi_n` / `chiN` definitions, untouched.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 049` and `redteam exec-mathematica 049` and confirm the new check `k_n satisfies D/N Neumann boundary` appears in both saved outputs with residual 0, and that no assertion of the form `k_n - (n+1/2) pi / L` remains.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`
- summary: Replaced tautological half-wave momentum equality checks with D/N Neumann boundary residual checks.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:27,48,53`

**Issue:**
The `.wl` script duplicates the SymPy choreography line for line (same helper-function set, same variable choreography, same assertion sequence, same banner text). To break the mirror, remove the `uniformDnOverlap` helper and derive `overlapFormula` directly from `Integrate[chiN, {s, 0, l}]`, with `i0` obtained from `overlapFormula` by substitution `n -> 0` rather than from the helper. This routes the Mathematica overlap derivation through the integrator (an independent engine implementation), instead of echoing the SymPy closed form via a parallel helper.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`:

1. Delete line 27 (the `uniformDnOverlap` helper definition):
```wolfram
uniformDnOverlap[n_, l_] := FullSimplify[Sqrt[2 l]/((n + 1/2) Pi)];
```

2. Replace line 48
Before:
```wolfram
overlapFormula = uniformDnOverlap[n, l];
```
After:
```wolfram
overlapFormula = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];
```

3. Replace line 53
Before:
```wolfram
i0 = uniformDnOverlap[0, l];
```
After:
```wolfram
i0 = FullSimplify[overlapFormula /. n -> 0];
```

Leave the SymPy script and all other lines of the Mathematica script unchanged. In particular, do not alter `overlapRatio`, `twinSupportRatio`, the zeta-twin assertions, or the carry-forward print block.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 049` and confirm: (a) the saved output still shows `uniform overlap integral = 0` with `PASS`, (b) the closed-form line `I_n closed form` matches the integrator result (e.g. `(2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` or algebraic equivalent), and (c) the script no longer contains `uniformDnOverlap`. The script must still exit 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`
- summary: Removed the Mathematica overlap helper and derived the overlap formula and zeroth overlap from direct integration.
- deviation: none
