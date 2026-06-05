---
unit_id: 077
batch: III.4
verdict: verified
overall: verified
material_change: false
verified_at: 2026-06-05
---

# Verification — unit 077 (pass 2)

**Verdict: verified.** (Independent read+reason verification by a clean verify agent;
this record persists the agent's returned verdict — the agent completed the analysis
but did not write the file itself.)

## Finding F1 — symbol_assumption_error — RESOLVED

The fix is exactly the prescribed declaration split at
`scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:33`:

```
-xi, alpha_r, lambda_mu = sp.symbols("xi alpha_r lambda_mu", positive=True, real=True)
+xi = sp.symbols("xi", real=True)
+alpha_r, lambda_mu = sp.symbols("alpha_r lambda_mu", positive=True, real=True)
```

`xi` now carries `real=True` only; `alpha_r`, `lambda_mu` remain `positive=True, real=True`
(they genuinely are positive parameters). Nothing else in the file changed.

### Adversarial checks (no surviving assertion relies on `xi > 0`)
- **A1** (`I_f - 1/3`, line ~40): the integral at line ~37 has explicit `(-oo, oo)` bounds,
  so the domain is set by the bounds, not by the symbol assumption. Result `I_f = 1/3` is
  independently reproduced by the Mathematica engine (which declares `xi ∈ Reals`, never
  positive). Not reliant on `xi > 0`.
- **A2** (back-substitution `1 - alpha_r·S(xi_*)**2`): uses `xi_star` (an expression in
  `alpha_r`), not the free symbol `xi`. Independent of `xi`'s assumptions.
- **A3–A7** (numeric block): all use mpmath floats (`R1`, `R2`, `Theta_chi`, `Theta_J`),
  never the SymPy `xi`. Independent.

There is no second `xi` declaration anywhere; the positivity is gone everywhere.

## Output freshness / no result change
The orchestrator independently re-ran SymPy (exit 0). The committed transcript is
byte-identical to HEAD: `I_f - 1/3 = 0`, `xi_*(alpha_r=10) = -0.38558106921542562404`,
and the four numeric values (`<rho>_chi`, `<rho^2>_chi`, `Theta_w^(chi)`, `Theta_w^(J)`)
unchanged — confirming a safe latent-trap removal that moves no result. Mathematica
(already correct domain) untouched and agrees to ~1e-50.

`material_change: false`.
