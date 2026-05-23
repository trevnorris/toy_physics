---
unit_id: 057
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-22T23:42:23Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 057

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:107-110`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:90-93`

**Issue:**
In both scripts, the `y_req identity` assertion subtracts the literal RHS of `y_req_sq`'s definition from `y_req_sq` itself, producing an algebraic tautology that any expression for y_req² would pass. The check fails to verify that y_req² actually solves the defining equation ζ_req = Ω_Pe²(κ+π²/4)/(κ+y²). Replace the tautological self-subtraction with a round-trip substitution into the defining equation (analogous to the `kappa_req defining equation` check the Mathematica script already performs for κ_req at lines 82-85).

**Required change:**

SymPy patch. The current block at lines 107-110 reads:

```python
expect_zero(
    "y_req identity",
    y_req_sq - ((Omega_Pe**2 / zeta_req) * (kappa + sp.pi**2 / 4) - kappa),
)
```

Replace it with:

```python
expect_zero(
    "y_req defining equation",
    zeta_req - sp.simplify(
        Omega_Pe**2 * (kappa + sp.pi**2 / 4) / (kappa + y_req_sq)
    ),
)
```

Do not modify the definition of `y_req_sq` at line 92-93 or the surrounding `print` statement at line 101. Leave the rest of the script untouched.

Mathematica patch. The current block at lines 90-93 reads:

```mathematica
expectZero[
  "y_req identity",
  yReqSq - ((omegaPe^2/zetaReq) (kappa + Pi^2/4) - kappa)
];
```

Replace it with:

```mathematica
expectZero[
  "y_req defining equation",
  zetaReq - FullSimplify[
    omegaPe^2 (kappa + Pi^2/4)/(kappa + yReqSq),
    Assumptions -> $Assumptions
  ]
];
```

Do not modify the definition of `yReqSq` at line 76 or the surrounding `Print` statement at line 80. Leave the rest of the script untouched.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 057` and `redteam exec-mathematica 057` and confirm:
1. Both scripts exit 0.
2. The SymPy saved output now contains a line `y_req defining equation = 0` in place of the prior `y_req identity = 0`.
3. The Mathematica saved output now contains `PASS: y_req defining equation` in place of `PASS: y_req identity`.
4. The new assertion is non-tautological: substituting y² = (Ω²/ζ_req)(κ+π²/4) − κ into Ω²(κ+π²/4)/(κ+y²) collapses to ζ_req only because the y_req² formula is correctly derived; any other expression for y_req_sq would leave a nonzero residual.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl`
- summary: Replaced the tautological y_req identity checks with round-trip defining-equation residual checks.
- deviation: none
