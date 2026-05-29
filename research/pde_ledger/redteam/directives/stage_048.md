---
unit_id: 048
batch: III.1
created_at: 2026-05-26T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T07:48:27Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 048

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py:58`

**Issue:**
The SymPy script computes and prints the right-endpoint limit `lim_{xi -> 1^-} F_tr` at line 58 but never asserts it. The notes (`notes/stages/moving_throat_pde_stage048_support_compensation_theorem.md`, section 3, lines 124-126) declare this divergence as load-bearing for the IVT existence of `xi_req in (0,1)`. The Mathematica counterpart (`mathematica/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.wl`, lines 55-64) asserts a strengthened softening-coefficient identity that captures the divergence and its leading coefficient. The SymPy script must be brought into parity, asserting the same strengthened identity.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py`, immediately after the existing line 58:

```python
print("limit xi->1^- of F_tr =", sp.limit(F_tr, xi, 1, dir="-"))
```

insert the following lines (do not delete the print at line 58 — keep it for context):

```python
soft_coeff = sp.simplify(sp.limit((1 - xi) * F_tr, xi, 1, dir="-"))
print("softening coefficient for F_tr =", soft_coeff)
expect_zero(
    "(1-xi) F_tr softening coefficient",
    soft_coeff
    - (9 * delta + 9 + 2 * R ** 2) ** 2 * (9 * delta + 9 + 2 * R) ** 2
    / (81 * (9 * delta ** 2 + 18 * delta + 9 + 2 * R ** 2) ** 2),
)
```

The closed-form RHS is taken from the Mathematica script's `softCoeffExpected` at lines 59-62 of `mathematica/moving_throat_pde_stage048_support_compensation_theorem_mathematica_audit.wl`:

```mathematica
softCoeffExpected = FullSimplify[
  (9*delta + 9 + 2*r^2)^2*(9*delta + 9 + 2*r)^2/(81*(9*delta^2 + 18*delta + 9 + 2*r^2)^2),
  Assumptions -> delta > 0 && r > 0
];
```

No other lines change. The existing assertion at line 57 (`F_tr(xi=0)-1`) is unchanged. The existing assertion at lines 79-82 (`if limit_phys != sp.oo: raise AssertionError(...)`) is unchanged. The "M_crit - G_tr" block (lines 59-66) is unchanged.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 048`. Confirm:
1. The script exits 0.
2. A new line `softening coefficient for F_tr = ...` appears in the output transcript between the existing `limit xi->1^- of F_tr = oo` line and the existing `M_crit - G_tr = ...` line.
3. A new line `(1-xi) F_tr softening coefficient = 0` appears immediately after.
4. The SymPy output for `soft_coeff` textually agrees (modulo ordering/sign convention) with the Mathematica output line 20: `softening coefficient for F_tr = ((9 + 9*delta + 2*r)^2*(9 + 9*delta + 2*r^2)^2)/(81*(9*(1 + delta)^2 + 2*r^2)^2)`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py`
- summary: Added the required SymPy softening-coefficient computation and zero-residual assertion immediately after the existing right-endpoint limit print.
- deviation: none
