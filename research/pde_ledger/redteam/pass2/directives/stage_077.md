---
unit_id: 077
batch: III.4
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T21:37:59Z
verification_status: pending
needs_user_resolution: false
findings_applied: 1
findings_blocked: 0
---

# Codex directive — unit 077

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line range named.

After editing, RUN the affected script (`python3 scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`) and iterate until it exits 0 with all in-file checks passing.

Do NOT touch paper.tex, notes/, or any prose documents. Do NOT touch the Mathematica script (it is already correct).

## F1 — symbol_assumption_error

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:33`

**Issue:** The integration variable `xi` is declared `positive=True`, but it is the variable of the symmetric full-line integral `If = sp.integrate(sp.simplify(chi**2), (xi, -sp.oo, sp.oo))` (line 36) and the physically load-bearing support cut point is `xi_* ≈ -0.3856 < 0`. The positivity assumption contradicts the setup and disagrees with the Mathematica mirror (which correctly declares `Element[{xi, alphaR}, Reals]`). It is currently dormant (the explicit `(-oo, oo)` integral bounds, not the symbol assumption, set the domain, so `I_f = 1/3` is unaffected), but it is a latent trap for any future symbolic step added to this block.

**Required change:**
Replace line 33:

```python
xi, alpha_r, lambda_mu = sp.symbols("xi alpha_r lambda_mu", positive=True, real=True)
```

with (keep `alpha_r`, `lambda_mu` positive — they genuinely are; make `xi` real-only):

```python
xi = sp.symbols("xi", real=True)
alpha_r, lambda_mu = sp.symbols("alpha_r lambda_mu", positive=True, real=True)
```

Do not change any other line. The symbolic and numeric results must remain identical (`I_f = 1/3`, `xi_* ≈ -0.38558106921542562404`, `<rho>_chi = 0.19261900555649309…`, `<rho^2>_chi = 0.16274529400326462…`, `Theta_w^(chi) = 4.0686323500816155…`, `Theta_w^(J) = 0.92755203253930797…`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 077` and confirm the script exits 0 and the saved output still prints `I_f - 1/3 = 0`, the cut-point line, and the four numeric values above unchanged (the positivity removal cannot alter an explicit-bounds definite integral or the float-based numeric block).

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`
- summary: Split the symbol declaration so `xi` is real-only while `alpha_r` and `lambda_mu` remain positive real parameters.
- deviation: none
