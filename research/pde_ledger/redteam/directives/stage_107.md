---
unit_id: 107
batch: IV.2
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 107

Apply the finding below. After applying, append an `## Applied: F1` block under the finding
with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line
ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py:58-60`

**Issue:**

The SymPy script computes `Sigma2_evenmatch = sp.simplify(sol[Sigma2])` and
`Sigma4_evenmatch = sp.simplify(sol[Sigma4])` (lines 58-59) but only prints them. The
paper's notes (`notes/stages/moving_throat_pde_stage107_general_dtn_deformation.md`) state
two boxed closed forms:
```
Sigma_2 = -(3 S beta^2 - 3 S + Sigma_0)/9,
Sigma_4 = -(3 S beta^4 - 3 S + Sigma_0)/27.
```
The Mathematica twin asserts both at `mathematica/...stage107..._audit.wl:68-69` via
`expectZero["Sigma2 exact formula", ...]` and `expectZero["Sigma4 exact formula", ...]`. The
SymPy script has no analogous assertion, so two of the three boxed Stage 107 deliverables
are unguarded in SymPy.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py`,
immediately after the existing `print('Sigma4_evenmatch = ...')` line (currently line 59) and
before the `chi_even = sp.simplify(chiQ.subs(sol))` line (currently line 61), insert two new
`expect_zero` checks that lock in the paper's boxed closed forms.

Before (lines 58-60):
```python
print('Sigma2_evenmatch =', sp.simplify(sol[Sigma2]))
print('Sigma4_evenmatch =', sp.simplify(sol[Sigma4]))

```

After:
```python
print('Sigma2_evenmatch =', sp.simplify(sol[Sigma2]))
print('Sigma4_evenmatch =', sp.simplify(sol[Sigma4]))
expect_zero(
    'Sigma2 exact formula',
    sol[Sigma2] - (-(3*S*beta**2 - 3*S + Sigma0)) / sp.Integer(9),
)
expect_zero(
    'Sigma4 exact formula',
    sol[Sigma4] - (-(3*S*beta**4 - 3*S + Sigma0)) / sp.Integer(27),
)

```

Notes for Codex:
- Use the existing `expect_zero` helper defined at the top of the file (lines 13-17); do not
  introduce a new helper.
- Use `sp.Integer(9)` / `sp.Integer(27)` so the division stays symbolic (the existing code
  in this script uses `sp.Integer(...)` for the same purpose; see lines 27 and 30).
- The `S`, `beta`, `Sigma0`, `Sigma2`, `Sigma4`, and `sol` names are already in scope (see
  lines 21-22 and 54).
- Do not alter any other line in the script.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 107` and confirm that:
1. The transcript now contains two new lines, `Sigma2 exact formula = 0` and
   `Sigma4 exact formula = 0`, both before the existing `chi_Q under canonical-even matching`
   line.
2. The script still exits 0.
3. The existing `normalized expansion direct-formula = 0` and `chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0) = 0` lines still appear.
