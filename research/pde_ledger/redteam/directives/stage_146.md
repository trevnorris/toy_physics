---
unit_id: 146
batch: IV.5
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 146

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:71-75`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:85-90`

**Issue:**
The "affine law" assertions are restatements of distributive arithmetic on the script's own definitions, not tests of the paper's affine identities for the source moments. `g_eps = (1-eps)*gminus + eps*gbar` already equals `gminus + eps*(gbar - gminus)` by algebra, so `expect_zero` cannot fail. Replace these with assertions that compare the *integral* form of \(\bar g_\epsilon\) and \(\bar S_\epsilon\) (linear combination of integrals against the kernels) against the affine right-hand side.

**Required change:**

(SymPy) In `scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py`, replace lines 71-75:

Before:
```python
# Exact affine law for convex family moments.
g_eps = (1-eps)*gminus + eps*gbar
S_eps = (1-eps)*Sformula.subs(Pi, Pi_star) + eps*Sbar
expect_zero("g_eps affine law", sp.expand(g_eps - (gminus + eps*(gbar-gminus))))
expect_zero("S_eps affine law", sp.expand(S_eps - (Sformula.subs(Pi, Pi_star) + eps*(Sbar-Sformula.subs(Pi, Pi_star)))))
```

After:
```python
# Exact affine law for convex family moments — tested via integration,
# not via redefining the affine form. Use a concrete positive normalized
# test profile varsigma_test(x) = 6 x (1 - x) (a polynomial bump,
# positive on (0,1), with integral 1 on [0,1]).
varsigma_test = 6*x*(1-x)
Sigma_eps = (1-eps)*Sigma.subs(Pi, Pi_star) + eps*varsigma_test
gbar_phys = sp.integrate(Sigma_eps*sp.cos(sp.pi*x/2), (x, 0, 1))
Sbar_phys = sp.integrate(Sigma_eps*Kq, (x, 0, 1))
gbar_v    = sp.integrate(varsigma_test*sp.cos(sp.pi*x/2), (x, 0, 1))
Sbar_v    = sp.integrate(varsigma_test*Kq, (x, 0, 1))
expect_zero(
    "g_eps affine law (integral form)",
    sp.simplify(gbar_phys - (gminus + eps*(gbar_v - gminus))),
)
expect_zero(
    "S_eps affine law (integral form)",
    sp.simplify(Sbar_phys - (Sformula.subs(Pi, Pi_star) + eps*(Sbar_v - Sformula.subs(Pi, Pi_star)))),
)
```

Notes for Codex:
- `Sigma` and `Kq` are already defined above (lines 21-22). `Sigma.subs(Pi, Pi_star)` substitutes the numeric `Pi_*` so the integrals become numeric in `eps` only.
- `gminus` (line 47) and `Sformula.subs(Pi, Pi_star)` (already used on line 73 above) carry the numeric branch values; the residuals should reduce to 0 because both `gbar_phys` and the affine RHS share `gminus = integral of Sigma_*(x) cos(pi x/2)` and `Sformula.subs(Pi, Pi_star) = integral of Sigma_*(x) Kq(x)` (the latter only holds if F2 is also applied — see F2).
- `sp.simplify` may be slow if `sp.integrate` returns symbolic combinations of erf-like terms; if so, fall back to a numeric check at, e.g., `eps -> sp.Rational(1, 10)` and `eps -> sp.Rational(1, 2)` with `sp.N(..., 30)` and `assert abs(...) < 1e-25`.

(Mathematica) In `mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`, replace lines 85-90:

Before:
```mathematica
gEps = (1 - eps)*gMinus + eps*gBar;
sEps = (1 - eps)*sStar + eps*sBar;
expectZero["g_eps affine law", Expand[gEps - (gMinus + eps*(gBar - gMinus))]];
resSEps = Chop[Expand[sEps - (sStar + eps*(sBar - sStar))]];
Print["S_eps affine law = ", fmt[resSEps]];
If[TrueQ[resSEps === 0], pass["S_eps affine law"], fail["S_eps affine law", resSEps]];
```

After:
```mathematica
(* Affine laws tested via integral form, not via algebraic restatement. *)
varsigmaTest = 6*x*(1 - x);
sigmaEps = (1 - eps)*(sigma /. p -> pStar) + eps*varsigmaTest;
gBarPhys = Integrate[sigmaEps*Cos[Pi*x/2], {x, 0, 1}];
sBarPhys = Integrate[sigmaEps*kq, {x, 0, 1}];
gBarV    = Integrate[varsigmaTest*Cos[Pi*x/2], {x, 0, 1}];
sBarV    = Integrate[varsigmaTest*kq, {x, 0, 1}];
expectZero[
  "g_eps affine law (integral form)",
  FullSimplify[gBarPhys - (gMinus + eps*(gBarV - gMinus))]
];
resSEps = Chop[N[FullSimplify[sBarPhys - (sStar + eps*(sBarV - sStar))], 40]];
Print["S_eps affine law (integral form) = ", fmt[resSEps]];
If[TrueQ[resSEps === 0] || (NumericQ[resSEps] && Abs[resSEps] < 10^-25),
  pass["S_eps affine law (integral form)"],
  fail["S_eps affine law (integral form)", resSEps]
];
```

Notes for Codex:
- The Mathematica `S_eps` check uses `sStar` (numeric, 40-digit) rather than the symbolic `Sformula /. p -> pStar`. Because `sStar` and `Integrate[sigma*kq /. p -> pStar, {x,0,1}]` are equal but only to numeric precision, the residual is allowed to be `< 10^-25` rather than literal `=== 0`. This mirrors how the existing `kernel check at Pi=...` already accepts a `10^-12` tolerance for numeric comparisons.

**Claim manifest:** n/a (assertion-rewrite, not a missing-script finding).

**Verification command:**
The verifier will run `redteam exec-sympy 146` and `redteam exec-mathematica 146`. Both should still exit 0. The new transcripts must contain the line `g_eps affine law (integral form) = 0` and either `S_eps affine law (integral form) = 0` or a residual under `10^-25`.

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:28` (insertion point: after line 28, before line 30)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:48` (insertion point: after line 48)

**Issue:**
The closed-form moments `gPi`/`gFormula` and `Sformula`/`sFormula` are the load-bearing inputs to every numeric value in this unit. The Mathematica script symbolically checks `gFormula` against direct integration of `sigma*Cos[Pi*x/2]` (line 48), but the SymPy script never performs this symbolic check; \(\mathcal S_q(\Pi)\) is verified only at three numeric samples in both engines.

**Required change:**

(SymPy) Insert after line 28 in `scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py`:

```python
# Symbolic verification that the closed-form moments equal the canonical
# integrals against the kernel functions. This anchors the formulas before
# the numeric Pi_* root-finder runs.
gPi_direct = sp.integrate(Sigma*sp.cos(sp.pi*x/2), (x, 0, 1))
expect_zero("g(Pi) direct-formula", sp.simplify(gPi_direct - gPi))
Sq_direct = sp.integrate(Sigma*Kq, (x, 0, 1))
expect_zero("S_q(Pi) direct-formula", sp.simplify(Sq_direct - Sformula))
```

Notes for Codex:
- If `sp.integrate` returns an unevaluated `Integral` for one of these, fall back to:
  ```python
  for pp in [sp.Rational(7,10), sp.Rational(11,10), sp.Rational(17,10), sp.Rational(23,10)]:
      diff = sp.N(sp.integrate((Sigma*sp.cos(sp.pi*x/2)).subs(Pi, pp), (x,0,1)) - gPi.subs(Pi, pp), 30)
      print(f"g(Pi) numeric sample Pi={pp}: diff={diff}")
      if abs(float(diff)) > 1e-25:
          raise AssertionError("g(Pi) numeric mismatch")
  ```
  and analogously for `Sq_direct`. Prefer symbolic; fall back only if symbolic returns `Integral`.

(Mathematica) Insert after line 48 in `mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`:

```mathematica
sDirect = FullSimplify[Integrate[sigma*kq, {x, 0, 1}], Assumptions -> p > 0];
expectZero["S_q(Pi) direct-formula", sDirect - sFormula];
```

Notes for Codex:
- `sigma`, `kq`, `sFormula` are all defined on lines 38-42 above.
- If `Integrate` returns unevaluated under `Assumptions -> p > 0`, retry with `Assumptions -> p > 0 && p < Pi/2` (the kernel denominator `(kap^2 - p^2)` is the only thing that could cause a branch issue and the small-`p` regime suffices for the formula identity). If still unresolved, fall back to a 4-point numeric check at `p -> {7/10, 11/10, 17/10, 23/10}` with WorkingPrecision 50 and tolerance 10^-25.

**Claim manifest:** n/a.

**Verification command:**
The verifier will run `redteam exec-sympy 146` and `redteam exec-mathematica 146`. Both should still exit 0. The SymPy transcript must contain `g(Pi) direct-formula = 0` and `S_q(Pi) direct-formula = 0`; the Mathematica transcript must contain `S_q(Pi) direct-formula = 0`.

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:32`

**Issue:**
The Mathematica banner literal says `STAGE 129 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS`. This is the stage-146 audit script. The wrong stage number is stale scaffolding.

**Required change:**

Edit line 32 from:
```mathematica
banner["STAGE 129 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];
```
to:
```mathematica
banner["STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];
```

**Claim manifest:** n/a.

**Verification command:**
The verifier will confirm the new Mathematica transcript contains the literal string `STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS` in the banner block.
