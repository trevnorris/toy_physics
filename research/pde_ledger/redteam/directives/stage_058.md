---
unit_id: 058
batch: III.2
created_at: 2026-05-22T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-22T17:44:52-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 058

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:49-62`

**Issue:**
The Mathematica script imports the closed-form antiderivatives `fc`, `fs` and the combined drop `delta` from the SymPy script rather than deriving them independently. Lines 49-50 declare `fc`, `fs` as the same character-for-character ansatz as SymPy lines 56-57; lines 59-62 transcribe SymPy's `Delta` expression. Mathematica is fully capable of computing these via `Integrate[]` symbolically; the second-engine policy requires it to do so.

**Required change:**
1. Replace lines 49-50 with independent symbolic integration. Before:
   ```
   fc = Exp[Pe*x]*(Pe*Cosh[alpha*x] - alpha*Sinh[alpha*x])/(Pe^2 - alpha^2);
   fs = Exp[Pe*x]*(Pe*Sinh[alpha*x] - alpha*Cosh[alpha*x])/(Pe^2 - alpha^2);
   ```
   After:
   ```
   fc = FullSimplify[Integrate[Exp[Pe*x]*Cosh[alpha*x], x], Assumptions -> $Assumptions && Pe != alpha];
   fs = FullSimplify[Integrate[Exp[Pe*x]*Sinh[alpha*x], x], Assumptions -> $Assumptions && Pe != alpha];
   ```
   Keep the existing `expectZero` derivative cross-checks at lines 51-52 as regression tests (they should still pass since differentiating the integrated antiderivative reproduces the integrand). Relabel them to "Ic antiderivative regression (Mma re-derivation)" and "Is antiderivative regression (Mma re-derivation)".

2. Replace lines 59-62 with the independent integral form. Before:
   ```
   delta = FullSimplify[
     Pe/(Exp[Pe] - 1)*((1 - Cosh[alpha])*ic + (eta/alpha + Sinh[alpha])*is)/w,
     Assumptions -> $Assumptions
   ];
   ```
   After:
   ```
   delta = FullSimplify[
     Integrate[kernel*sigmaPe, {x, 0, 1}, Assumptions -> $Assumptions && Pe != alpha],
     Assumptions -> $Assumptions && Pe != alpha
   ];
   ```
   This rebuilds Delta from the physical definition (Delta = integral of K * Sigma) rather than from the SymPy-imported combination. Add immediately after this line a regression check that the new `delta` agrees with the old combination form, labelled `"delta agrees with (1-cosh)Ic + (eta/alpha+sinh)Is combination"`:
   ```
   deltaCombination = FullSimplify[
     Pe/(Exp[Pe] - 1)*((1 - Cosh[alpha])*ic + (eta/alpha + Sinh[alpha])*is)/w,
     Assumptions -> $Assumptions && Pe != alpha
   ];
   expectZero["delta independent integral matches combination form", delta - deltaCombination];
   ```

3. Relabel the remaining cross-engine-shared labels to mark them as Mathematica re-derivations. At line 41 change `"Kprime identity"` to `"Kprime identity (Mma re-derivation)"`. At line 47 change `"Sigma normalization"` to `"Sigma normalization (Mma re-derivation)"`. At line 71 change `"Delta0 formula"` to `"Delta0 formula (Mma re-derivation)"`. At line 74 change `"Delta0 integral identity"` to `"Delta0 integral identity (Mma re-derivation)"`. At line 83 change `"Delta_inf formula"` to `"Delta_inf formula (Mma re-derivation)"`. At line 92 change `"weak-coupling constant term"` to `"weak-coupling constant term (Mma re-derivation)"`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 058` and confirm:
(a) the script source contains `Integrate[Exp[Pe*x]*Cosh[alpha*x], x` and `Integrate[Exp[Pe*x]*Sinh[alpha*x], x`,
(b) the script source contains `Integrate[kernel*sigmaPe, {x, 0, 1}`,
(c) the new `"delta independent integral matches combination form"` check appears and passes,
(d) the relabelled assertions appear and pass,
(e) the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`
- summary: Replaced copied antiderivative and Delta formulas with Mathematica symbolic integrations and relabelled cross-engine checks as Mathematica re-derivations.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:86-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:88-92`

**Issue:**
Docstring items (5) "exact sharp-bottom endpoint Delta_inf" and (6) "fixed-point branch bracket Pe_* in [Xi Delta_0, Xi Delta_inf]" are not exercised. The script only prints Pe_lo, Pe_hi (no monotonicity / non-emptiness assertion), and Delta_inf is defined as K(1) without linking it to lim_{Pe -> infinity} Delta(Pe), which is the physical meaning of "sharp-bottom endpoint".

**Required change:**

1. In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`, insert the following block between line 86 (after `print("Pe_hi = Xi Delta_inf =", Pe_hi)`) and line 88 (`# Weak-coupling branch law.`):
   ```python
   # Bracket non-emptiness: Delta_inf >= Delta_0 for all alpha, eta > 0.
   bracket_gap = sp.simplify(sp.together(Delta_inf - Delta0))
   bracket_gap_expected = sp.simplify(
       ((alpha**2 - eta) * (sp.cosh(alpha) - 1) + alpha * eta * sp.sinh(alpha))
       / (alpha**2 * W)
   )
   expect_zero("bracket gap closed form", bracket_gap - bracket_gap_expected)
   for a_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(3)]:
       for e_val in [sp.Rational(1, 10), sp.Integer(1), sp.Integer(10)]:
           val = float(sp.N(bracket_gap.subs({alpha: a_val, eta: e_val})))
           if val <= 0:
               raise AssertionError(
                   f"bracket gap non-positive at alpha={a_val}, eta={e_val}: {val}"
               )
   print("bracket gap positivity sweep = PASS")

   # Delta_inf is the sharp-bottom (Pe -> oo) endpoint of Delta(Pe).
   Delta_inf_limit = sp.simplify(sp.limit(Delta, Pe, sp.oo))
   print("Delta(Pe -> oo) =", Delta_inf_limit)
   expect_zero("Delta_inf as Pe -> oo limit", Delta_inf_limit - Delta_inf_expected)
   ```

2. In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`, insert the following block between line 88 (after `Print["Pe_hi = Xi Delta_inf = ", fmt[peHi]];`) and line 90 (the `deltaSeries = ...` line):
   ```
   bracketGap = FullSimplify[deltaInfExpected - delta0Expected, Assumptions -> alpha > 0 && eta > 0];
   bracketGapExpected = FullSimplify[
     ((alpha^2 - eta)*(Cosh[alpha] - 1) + alpha*eta*Sinh[alpha])/(alpha^2*w),
     Assumptions -> alpha > 0 && eta > 0
   ];
   expectZero["bracket gap closed form", bracketGap - bracketGapExpected];
   bracketGapValues = Flatten[Table[
     N[bracketGap /. {alpha -> aV, eta -> eV}],
     {aV, {1/10, 1, 3}}, {eV, {1/10, 1, 10}}
   ]];
   If[AnyTrue[bracketGapValues, # <= 0 &],
     fail["bracket gap positivity sweep", bracketGapValues],
     pass["bracket gap positivity sweep"]
   ];

   deltaInfLimit = FullSimplify[
     Limit[delta, Pe -> Infinity, Assumptions -> alpha > 0 && eta > 0],
     Assumptions -> alpha > 0 && eta > 0
   ];
   Print["Delta(Pe -> oo) = ", fmt[deltaInfLimit]];
   expectZero["Delta_inf as Pe -> oo limit", deltaInfLimit - deltaInfExpected];
   ```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 058` and `redteam exec-mathematica 058` and confirm:
(a) new assertion labels "bracket gap closed form", "bracket gap positivity sweep", "Delta_inf as Pe -> oo limit" appear and pass in both engines,
(b) printed line `Delta(Pe -> oo) = ...` matches the existing `Delta_inf` printed value symbolically,
(c) both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`
- summary: Added bracket-gap closed-form, positivity sweep, and sharp-bottom Pe-to-infinity endpoint checks in both engines.
- deviation: SymPy rewrites the Pe-to-infinity residual in exponential form before `expect_zero` so the equivalent limit expression simplifies to zero.

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:80`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:83`

**Issue:**
`expect_zero("Delta_inf formula", Delta_inf - Delta_inf_expected)` and its Mathematica mirror are algebraic substitution: `Delta_inf := simplify(K.subs(x, 1))` equals `(cosh(alpha)+(eta/alpha)*sinh(alpha)-1)/W` by direct evaluation (using cosh(0)=1), and `Delta_inf_expected` is that same expression. The check cannot fail and adds no physics content. After F2 adds the genuine "Delta_inf as Pe -> oo limit" check, this line should be relabelled to make clear it is a substitution sanity print, not the docstring's "exact sharp-bottom endpoint" verification.

**Required change:**

1. In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:80`, change:
   ```python
   expect_zero("Delta_inf formula", Delta_inf - Delta_inf_expected)
   ```
   to:
   ```python
   expect_zero("Delta_inf direct substitution (sanity)", Delta_inf - Delta_inf_expected)
   ```

2. In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:83`, change (after F1 relabel takes effect — apply this F3 edit AFTER F1):
   ```
   expectZero["Delta_inf formula (Mma re-derivation)", deltaInf - deltaInfExpected];
   ```
   to:
   ```
   expectZero["Delta_inf direct substitution (sanity, Mma re-derivation)", deltaInf - deltaInfExpected];
   ```

**Verification command:**
After Codex applies, the verifier confirms the assertion labels at the cited line are renamed exactly as specified and the assertions still pass (the underlying check is unchanged — only the label changes, since the substantive replacement lives in F2).

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`
- summary: Relabelled the Delta_inf direct substitution checks as sanity checks in both engines.
- deviation: none

## F4 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:89-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:90-92`

**Issue:**
`expect_zero("weak-coupling constant term", Delta_series.coeff(Pe, 0) - Delta0)` is mathematically guaranteed for any function analytic at Pe = 0 (the zeroth Taylor coefficient equals the limit at 0). Since Delta is analytic at Pe = 0 (the Pe/(exp(Pe)-1) factor extends smoothly), this assertion cannot fail. It tests the series engine, not the physics. The docstring's weak-coupling claim should constrain the *first-order* coefficient.

**Required change:**

1. In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`, replace lines 89-91:
   ```python
   Delta_series = sp.series(Delta, Pe, 0, 2).removeO()
   print("Delta(Pe) small-Pe series =", Delta_series)
   expect_zero("weak-coupling constant term", sp.expand(Delta_series).coeff(Pe, 0) - Delta0)
   ```
   with:
   ```python
   Delta_series = sp.series(Delta, Pe, 0, 2).removeO()
   print("Delta(Pe) small-Pe series =", Delta_series)
   expect_zero("weak-coupling constant term", sp.expand(Delta_series).coeff(Pe, 0) - Delta0)
   Pe1_coeff = sp.simplify(sp.expand(Delta_series).coeff(Pe, 1))
   print("Delta(Pe) first-order coefficient =", Pe1_coeff)
   Pe1_val = sp.N(Pe1_coeff.subs({alpha: sp.Integer(1), eta: sp.Integer(1)}))
   if sp.Abs(Pe1_val) < sp.Rational(1, 10**8):
       raise AssertionError(
           f"weak-coupling first-order coefficient vanishes at alpha=eta=1: {Pe1_val}"
       )
   print("weak-coupling first-order coefficient nonvanishing at alpha=eta=1: PASS")
   ```

2. In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`, after line 92 (after the F1-relabelled `"weak-coupling constant term (Mma re-derivation)"` `expectZero` line), insert:
   ```
   pe1Coeff = FullSimplify[SeriesCoefficient[deltaSeries, {Pe, 0, 1}], Assumptions -> alpha > 0 && eta > 0];
   Print["Delta(Pe) first-order coefficient = ", fmt[pe1Coeff]];
   pe1Val = N[pe1Coeff /. {alpha -> 1, eta -> 1}];
   If[Chop[pe1Val] === 0,
     fail["weak-coupling first-order coefficient vanishes at alpha=eta=1", pe1Val],
     pass["weak-coupling first-order coefficient nonvanishing at alpha=eta=1"]
   ];
   ```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 058` and `redteam exec-mathematica 058` and confirm:
(a) both scripts print `Delta(Pe) first-order coefficient = ...` with a non-zero closed-form expression,
(b) the new non-vanishing assertion appears and passes in both engines,
(c) both scripts exit 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`
- summary: Added first-order weak-coupling coefficient printing and a nonvanishing check in both engines.
- deviation: none
