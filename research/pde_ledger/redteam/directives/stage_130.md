---
unit_id: 130
batch: IV.4
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 130

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:33` (insert after this line)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:57` (insert after this line)

**Issue:**
The paper notes (`notes/stages/moving_throat_pde_stage130_mouth_bias_map.md:78-82`) box the strict inequality `dg_Π/dΠ > 0`, which underpins the "unique Π_*" half of the paper card's `\stagefield{Output}`. Both scripts verify only the covariance *identity* and a single-point derivative value at Π_*. Add a multi-point positivity check that exercises the strict-sign claim on a range that brackets Π_*.

**Required change (SymPy):**
After line 33 of `scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py` (immediately after the `raise AssertionError("Incorrect point-source limit.")` line and before the comment `# Unique Family-1 compensation point`), insert:

```python
# Strict monotonicity sweep: dg/dPi > 0 for Pi > 0 (notes boxed result)
dgPi = sp.diff(gPi, Pi)
for val in (sp.Rational(1, 10), sp.Rational(1, 2), sp.Integer(1), sp.Rational(15088, 10000), sp.Integer(3), sp.Integer(10)):
    deriv_val = sp.N(dgPi.subs(Pi, val), 30)
    print(f"dg/dPi at Pi={val} = {deriv_val}")
    if deriv_val <= 0:
        raise AssertionError(f"Strict monotonicity dg/dPi > 0 failed at Pi={val}.")
```

**Required change (Mathematica):**
After line 57 of `mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl` (immediately after the `expectZero["point-source limit", gInf - 1];` line and before the `piStar = ...` line), insert:

```mathematica
(* Strict monotonicity sweep: dg/dPiM > 0 for piM > 0 (notes boxed result) *)
dgPi = D[gPi, piM];
Module[{vals = {1/10, 1/2, 1, 15088/10000, 3, 10}, dv},
  Do[
    dv = N[dgPi /. piM -> v, 40];
    Print["dg/dpiM at piM=", fmt[v], " = ", fmt[dv]];
    If[TrueQ[dv > 0], pass["dg/dpiM > 0 at piM=" <> ToString[v]],
      fail["dg/dpiM > 0 at piM=" <> ToString[v], dv]],
    {v, vals}]
];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 130` and `redteam exec-mathematica 130` and confirm:
- SymPy transcript contains six lines beginning `dg/dPi at Pi=` and no `AssertionError`.
- Mathematica transcript contains six `PASS: dg/dpiM > 0 at piM=...` lines.
- Both scripts exit 0.

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:15` (insert after this line)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:42` (insert after this line)

**Issue:**
The paper notes (`notes/stages/moving_throat_pde_stage130_mouth_bias_map.md:34-40`) box the closed form `g_Π = 2Π(2Π e^Π + π) / ((4Π² + π²)(e^Π − 1))`. Neither script asserts equality between the integral-evaluated `g_Π` and this boxed form; agreement is only visual. Add a direct equality assertion in both engines.

**Required change (SymPy):**
After line 15 of `scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py` (immediately after `gPi = sp.simplify(sp.integrate(sigma*f, (z, 0, L)))`), insert:

```python
gPi_boxed = 2*Pi*(2*Pi*sp.exp(Pi) + sp.pi) / ((4*Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1))
if sp.simplify(gPi - gPi_boxed) != 0:
    raise AssertionError("g_Pi does not match paper boxed closed form.")
```

**Required change (Mathematica):**
After line 42 of `mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl` (immediately after `gPi = FullSimplify[Integrate[sigma*f, {z, 0, lM}], Assumptions -> $Assumptions];`), insert:

```mathematica
gPiBoxed = 2*piM*(2*piM*Exp[piM] + Pi) / ((4*piM^2 + Pi^2) * (Exp[piM] - 1));
expectZero["g_Pi matches paper boxed closed form", gPi - gPiBoxed];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 130` and `redteam exec-mathematica 130` and confirm:
- SymPy transcript shows no `AssertionError`; the printed `g_Pi = ...` line remains unchanged.
- Mathematica transcript contains a `PASS: g_Pi matches paper boxed closed form` line.
- Both scripts exit 0.
