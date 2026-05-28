---
unit_id: 143
batch: IV.5
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 143

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:39-57`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:53-69`

**Issue:** Neither script asserts the paper's bottom-line deliverables. The endpoint limits ($g_0=2/\pi$, $g_\infty=1$), the constants ($R_\infty\approx 0.1454544523$, $S_\infty=1$, $\widehat T_m/\sqrt\Pi\approx 0.7256691307$), and the positivity of the three decomposition pieces ($e^\Pi-1-\Pi-\Pi^2/2>0$, $\pi^2-2\pi>0$, $\pi^2/2-4>0$) are all printed but never gated by an assertion. A regression that broke any would still PASS.

**Required change:**

In `scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py`:

1. Add a helper near the top of the script (after the existing `expect_zero` definition, around line 17):
   ```python
   def expect_equal(name, lhs, rhs):
       diff = sp.simplify(lhs - rhs)
       print(f"{name}: lhs - rhs = {diff}")
       if diff != 0:
           raise AssertionError(f"{name} mismatch: {lhs} vs {rhs}")

   def expect_positive(name, expr):
       val = sp.simplify(expr)
       print(f"{name}: {val}")
       if not (val.is_positive is True or (val.is_number and float(val) > 0)):
           raise AssertionError(f"{name} not positive: {val}")
   ```

2. After the existing `print("  quadratic coeff =", sp.simplify(pi**2/2-4))` line (currently line 37), add three positivity assertions for the pieces that don't depend on Π:
   ```python
   expect_positive("pi**2 - 2*pi > 0", pi**2 - 2*pi)
   expect_positive("pi**2/2 - 4 > 0", pi**2/2 - 4)
   # exp-remainder positivity via Taylor coefficient
   exp_rem_series = sp.series(sp.exp(Pi) - 1 - Pi - Pi**2/2, Pi, 0, 5).removeO()
   expect_equal("exp remainder leading term is Pi**3/6", exp_rem_series.coeff(Pi, 3), sp.Rational(1, 6))
   ```

3. After the existing `print("lim_{Pi->oo} g_Pi =", ginf)` line (currently line 43), add:
   ```python
   expect_equal("lim_{Pi->0+} g_Pi == 2/pi", g0, 2/pi)
   expect_equal("lim_{Pi->oo} g_Pi == 1", ginf, sp.Integer(1))
   ```

4. After the existing `print("lim That/sqrt(Pi) =", that_ratio)` line (currently line 57), add:
   ```python
   expect_equal("R_infty == (1-r)**2/(1+r**2)", Rinf, (1-r)**2/(1+r**2))
   expect_equal("S_infty == 1", Sinf, sp.Integer(1))
   expect_equal("lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty))", that_ratio, sp.sqrt(sp.Rational(9,20)/(1-Rinf)))
   ```

In `mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl`:

5. Add a helper after the existing `expectZero` definition (around line 31):
   ```mathematica
   expectEqual[name_String, lhs_, rhs_] := Module[{res},
     res = FullSimplify[lhs - rhs, Assumptions -> $Assumptions];
     Print[name, ": lhs - rhs = ", fmt[res]];
     If[TrueQ[res === 0], pass[name], fail[name, res]];
   ];

   expectPositive[name_String, expr_] := Module[{val},
     val = FullSimplify[expr, Assumptions -> $Assumptions];
     Print[name, ": ", fmt[val]];
     If[TrueQ[Simplify[val > 0]], pass[name], fail[name, val]];
   ];
   ```

6. After the existing `Print["  quadratic coeff = ", fmt[FullSimplify[Pi^2/2 - 4]]];` line (currently line 51), add:
   ```mathematica
   expectPositive["Pi^2 - 2*Pi > 0", Pi^2 - 2*Pi];
   expectPositive["Pi^2/2 - 4 > 0", Pi^2/2 - 4];
   expRemSeries = Normal[Series[Exp[piM] - 1 - piM - piM^2/2, {piM, 0, 4}]];
   expectEqual["exp remainder leading term is piM^3/6", Coefficient[expRemSeries, piM, 3], 1/6];
   ```

7. After the existing `Print["lim_{Pi->oo} g_Pi = ", fmt[gInf]];` line (currently line 58), add:
   ```mathematica
   expectEqual["lim_{piM->0+} g_Pi == 2/Pi", g0, 2/Pi];
   expectEqual["lim_{piM->oo} g_Pi == 1", gInf, 1];
   ```

8. After the existing `Print["lim That/sqrt(Pi) = ", fmt[tHatRatio]];` line (currently line 69), add:
   ```mathematica
   expectEqual["R_infty == (1-r)^2/(1+r^2)", rInf, (1 - r)^2/(1 + r^2)];
   expectEqual["S_infty == 1", sInf, 1];
   expectEqual["lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty))", tHatRatio, Sqrt[(9/20)/(1 - rInf)]];
   ```
   (Note: the assertion on `sInf` here uses the variable `sInf` as it stands; F2 below will replace its hardcoded definition with an actual limit, after which this assertion gates that limit.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 143` and `redteam exec-mathematica 143` and confirm:
- sympy transcript shows new lines for `pi**2 - 2*pi > 0`, `pi**2/2 - 4 > 0`, `exp remainder leading term is Pi**3/6`, `lim_{Pi->0+} g_Pi == 2/pi`, `lim_{Pi->oo} g_Pi == 1`, `R_infty == (1-r)**2/(1+r**2)`, `S_infty == 1`, `lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty))` and exits 0.
- mathematica transcript shows matching `PASS:` lines and exits 0.

## F2 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:60-64`

**Issue:** Line 62 reads `sInf = 1;` — hardcoded. The full $S_q$ formula is defined on line 60 but its limit is never taken. `Sigma0` (the dynamical $\Pi/(1-R_qS_q)$) and `That` (the dynamical $\widehat T_m=\sqrt{9\Pi/(20[1-R_qS_q])}$) are never defined; `sigmaRatio` and `tHatRatio` are computed by algebraic substitution into the limiting form rather than as limits of dynamical objects. The corresponding SymPy lines 45-57 actually take limits via `sp.limit`.

**Required change:**

Replace the block on lines 60-69 of `mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl` with:

```mathematica
sQ = piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2));
rQ = (gPi - r)^2/(1 + r^2);
sigma0 = piM/(1 - rQ*sQ);
that = Sqrt[(9/20)*sigma0];

Clear[piInf2, piInf3, piInf4, piInf5];
rInf = FullSimplify[Limit[rQ /. piM -> piInf2, piInf2 -> Infinity], Assumptions -> piInf2 > 0];
sInf = FullSimplify[Limit[sQ /. piM -> piInf3, piInf3 -> Infinity], Assumptions -> piInf3 > 0];
sigmaRatio = FullSimplify[Limit[sigma0/piM /. piM -> piInf4, piInf4 -> Infinity], Assumptions -> piInf4 > 0];
tHatRatio = FullSimplify[Limit[that/Sqrt[piM] /. piM -> piInf5, piInf5 -> Infinity], Assumptions -> piInf5 > 0];

Print["R_infty = ", fmt[rInf]];
Print["S_infty = ", fmt[sInf]];
Print["lim Sigma0/Pi = ", fmt[sigmaRatio]];
Print["lim That/sqrt(Pi) = ", fmt[tHatRatio]];
```

Concretely: (a) replace the hardcoded `sInf = 1;` with a `Limit[sQ /. piM -> piInf3, piInf3 -> Infinity]`; (b) define `rQ`, `sigma0`, `that` as dynamical objects (in `piM`) before the limits; (c) compute `rInf` and `sigmaRatio` and `tHatRatio` as limits of those dynamical objects rather than as algebraic substitutions.

**Verification command:**
After Codex applies, `redteam exec-mathematica 143` should still print:
- `R_infty = 0.1454544522604201261...` (same numeric value as before, but now via a `Limit` of `rQ`).
- `S_infty = 1` (now derived as a limit of the full `sQ` formula).
- `lim Sigma0/Pi = 4107/(20*Pi*Sqrt[4107 - 100*Pi^2])` (now a limit of `sigma0/piM`).
- `lim That/sqrt(Pi) = (111*Sqrt[3/Pi])/(20*(4107 - 100*Pi^2)^(1/4))` (now a limit of `that/Sqrt[piM]`).

The F1 assertions on `R_infty`, `S_infty`, and `tHatRatio` should still pass after this refactor.

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:42-47`

**Issue:** The Mathematica block defining `num`, `decomp`, and calling `expectZero["numerator - exact positive decomposition", num - decomp]` is a token-for-token transliteration of the SymPy block at lines 28-32 of the `.py`. Both scripts subtract the same hand-rolled three-piece decomposition from the same intermediate `num`. Second-engine policy requires independent derivation.

**Required change:**

After the existing line 46 (`expectZero["numerator - exact positive decomposition", num - decomp];`), add an independent positivity check that does NOT mirror the SymPy choreography. Use Mathematica's `Reduce`:

```mathematica
(* Independent positivity verification: prove num > 0 for piM > 0 directly *)
numPositiveCheck = Reduce[num > 0, piM, Reals] /. {(Element[piM, Reals] && piM > 0) -> True, (piM > 0) -> True};
Print["Reduce[num > 0, piM, Reals] = ", fmt[numPositiveCheck]];
If[TrueQ[Simplify[numPositiveCheck === True || numPositiveCheck === (piM > 0)]],
  pass["num > 0 for piM > 0 via Reduce"],
  fail["num > 0 for piM > 0 via Reduce", numPositiveCheck]];
```

This proves $1-\mathfrak g_\Pi>0$ (equivalently $\mathfrak g_\Pi<1$) by a structurally different mechanism: Mathematica's `Reduce` over the reals, rather than subtraction against a hand-built decomposition. The two checks together (decomposition identity + Reduce positivity) cover the paper claim by independent paths.

**Verification command:**
After Codex applies, `redteam exec-mathematica 143` should print a new line `Reduce[num > 0, piM, Reals] = ...` followed by `PASS: num > 0 for piM > 0 via Reduce` (or, if Mathematica's `Reduce` returns `piM > 0` instead of `True`, the conditional accepts that form). Script must still exit 0.
