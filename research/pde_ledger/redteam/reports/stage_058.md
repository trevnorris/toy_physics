---
unit_id: 058
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 058 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage058_coupled_support_source_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.txt`

## What the script claims to verify

The script claims to verify, for the coupled support/source operator on the interval x in [0, 1]: (1) the closed-form derivative identity for the support kernel K_(alpha,eta)(x), (2) unit normalization of the constructive source family Sigma_Pe(x) = Pe*exp(Pe*x)/(exp(Pe)-1), (3) closed-form antiderivative formulas Fc, Fs for exp(Pe*x)*cosh(alpha*x) and exp(Pe*x)*sinh(alpha*x), (4) an exact "uniform-source drop" Delta_0 = lim_{Pe->0} Delta = integral of K over [0,1] reduced to eta*(cosh(alpha)-1)/(alpha^2 * W), (5) an exact "sharp-bottom endpoint" Delta_inf identified with K(1), (6) a fixed-point branch bracket Pe_* in [Xi*Delta_0, Xi*Delta_inf], and (7) consistency of the weak-coupling Taylor series with Delta_0 at zeroth order. Both engines emit identical closed forms for Delta, Ic, Is, Delta_0, and Delta_inf.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 45-49 | `expect_zero("Kprime identity", Kprime - <expected sum-of-sinh/cosh>)` | yes |
| A2 | sympy | 53 | `expect_zero("Sigma normalization", integrate(Sigma, (x,0,1)) - 1)` | yes |
| A3 | sympy | 58 | `expect_zero("Ic antiderivative check", diff(Fc,x) - exp(Pe*x)*cosh(alpha*x))` | yes |
| A4 | sympy | 59 | `expect_zero("Is antiderivative check", diff(Fs,x) - exp(Pe*x)*sinh(alpha*x))` | yes |
| A5 | sympy | 74 | `expect_zero("Delta0 formula", Delta0 - eta*(cosh(alpha)-1)/(alpha^2*W))` | yes |
| A6 | sympy | 75 | `expect_zero("Delta0 integral identity", Delta0 - integrate(K,(x,0,1)))` | yes |
| A7 | sympy | 80 | `expect_zero("Delta_inf formula", K.subs(x,1) - (cosh(alpha)+(eta/alpha)sinh(alpha)-1)/W)` | no (tautological substitution) |
| A8 | sympy | 91 | `expect_zero("weak-coupling constant term", Delta_series.coeff(Pe,0) - Delta0)` | no (Taylor 0-coeff == limit at 0 for analytic Delta) |
| A9 | mathematica | 40-43 | `expectZero["Kprime identity", ...]` (mirror of A1) | yes |
| A10 | mathematica | 47 | `expectZero["Sigma normalization", ...]` (mirror of A2) | yes |
| A11 | mathematica | 51-52 | `expectZero["Ic/Is antiderivative check", ...]` (mirror of A3,A4) | yes |
| A12 | mathematica | 71 | `expectZero["Delta0 formula", ...]` (mirror of A5) | yes |
| A13 | mathematica | 72-75 | `expectZero["Delta0 integral identity", ...]` (mirror of A6) | yes |
| A14 | mathematica | 83 | `expectZero["Delta_inf formula", ...]` (mirror of A7) | no |
| A15 | mathematica | 92 | `expectZero["weak-coupling constant term", ...]` (mirror of A8) | no |
| -- | both | -- | NO assertion that `Pe_lo <= Pe_hi`, NO assertion linking Delta_inf to lim_{Pe->inf} Delta | -- |

## Findings

### F1 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:31-92`

**What's wrong:**
The `.wl` script is a line-by-line transliteration of the `.py` script: same variable names (modulo case), same intermediate quantities in the same order (W/w, K/kernel, Sigma/sigmaPe, Fc/fc, Fs/fs, Ic/ic, Is/is, Delta/delta, Delta_0/delta0, Delta_inf/deltaInf, Pe_lo/peLo, Pe_hi/peHi), and the same closed-form antiderivative ansatz copied verbatim. Compare:

SymPy (lines 56-57):
```
Fc = sp.exp(Pe * x) * (Pe * sp.cosh(alpha * x) - alpha * sp.sinh(alpha * x)) / (Pe**2 - alpha**2)
Fs = sp.exp(Pe * x) * (Pe * sp.sinh(alpha * x) - alpha * sp.cosh(alpha * x)) / (Pe**2 - alpha**2)
```
Mathematica (lines 49-50):
```
fc = Exp[Pe*x]*(Pe*Cosh[alpha*x] - alpha*Sinh[alpha*x])/(Pe^2 - alpha^2);
fs = Exp[Pe*x]*(Pe*Sinh[alpha*x] - alpha*Cosh[alpha*x])/(Pe^2 - alpha^2);
```
Mathematica does NOT compute these antiderivatives independently (e.g., via `Integrate[Exp[Pe*x]*Cosh[alpha*x], x]`); it imports the same ansatz from SymPy and verifies it by differentiation. The `Delta` expression at lines 59-62 is also a direct rewrite of SymPy's lines 66-68 with identical algebraic structure. Even the `expectZero` labels are identical to SymPy's `expect_zero` labels ("Kprime identity", "Sigma normalization", "Ic antiderivative check", "Is antiderivative check", "Delta0 formula", "Delta0 integral identity", "Delta_inf formula", "weak-coupling constant term").

**Why this matters:**
The second-engine policy requires that Mathematica derive the result independently from the physical premises so that a bug in the SymPy derivation cannot propagate undetected. Transliteration defeats that purpose: if the imported `Fc`/`Fs` ansatz or the imported `Delta` expression contained a sign error, both engines would pass in lockstep.

**Required change:**
Re-derive the antiderivatives and the support drop independently in Mathematica. Replace the imported `fc`/`fs` ansatz at `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:49-50` with `fc = Integrate[Exp[Pe*x]*Cosh[alpha*x], x, Assumptions -> $Assumptions]` (and similarly for `fs`), then keep the derivative cross-check as a regression test. Replace the imported `delta` expression at lines 59-62 with `delta = FullSimplify[Integrate[kernel*sigmaPe, {x, 0, 1}], Assumptions -> $Assumptions]` — derive Delta as the integral of `kernel*sigmaPe` rather than transcribing SymPy's closed form. Relabel the engine-specific assertions ("Kprime identity" -> e.g. "Kprime identity (Mma re-derivation)") so transliteration cannot be inferred from label text.

**Verification:**
After re-derivation, the Mathematica output should still print `Delta_0`, `Delta_inf`, and the symbolic `Delta` agreeing with the SymPy output, but the script's source no longer contains the imported `Fc/Fs/Delta` closed forms — it constructs them from `Integrate[]`. The verifier confirms `Integrate[Exp[Pe*x]*Cosh[alpha*x], x, ...]` appears in the script and the assertion lines reference the integrated quantity rather than a literal closed-form ansatz.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:82-86`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:85-88`

**What's wrong:**
The module docstring (sympy line 11) claims to verify "6. fixed-point branch bracket Pe_* in [Xi Delta_0, Xi Delta_inf]". The script computes `Pe_lo = Xi*Delta_0` and `Pe_hi = Xi*Delta_inf` (sympy lines 83-84; mathematica lines 85-86) but only `print()`s them — no assertion. In particular there is NO check that `Pe_lo <= Pe_hi` (i.e., `Delta_0 <= Delta_inf`); without that, the "bracket" is not a bracket. The Mathematica mirror has the same gap. Both engines also leave the docstring claim "5. exact sharp-bottom endpoint Delta_inf" disconnected from the underlying physics: `Delta_inf` is defined as `K(1)` and matched against its own hand-simplified form, but the link `Delta_inf = lim_{Pe -> inf} Delta(Pe)` (which is the physical meaning of "sharp-bottom endpoint") is never asserted.

**Why this matters:**
A bracket claim with no monotonicity check is vacuous, and "sharp-bottom endpoint" without a Pe -> infinity limit verification leaves item (5) and (6) of the docstring outside the assertion net. If a sign or factor error caused `Delta_inf < Delta_0` for some valid (alpha, eta) regime, or caused `K(1) != lim_{Pe -> inf} Delta`, the script would still pass silently.

**Required change:**
1. In sympy at `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py` between lines 86 and 88, add an explicit assertion that `Delta_inf - Delta_0` is non-negative for `alpha > 0, eta > 0`. Use the simplified symbolic form:
   ```python
   bracket_gap = sp.simplify(sp.together(Delta_inf - Delta0))
   bracket_gap_expected = sp.simplify(
       ((alpha**2 - eta)*(sp.cosh(alpha) - 1) + alpha*eta*sp.sinh(alpha))
       / (alpha**2 * W)
   )
   expect_zero("bracket gap closed form", bracket_gap - bracket_gap_expected)
   # Numerical positivity sweep over alpha, eta > 0
   import itertools
   for a_val, e_val in itertools.product([sp.Rational(1,10), 1, 3], [sp.Rational(1,10), 1, 10]):
       val = float(bracket_gap.subs({alpha: a_val, eta: e_val}))
       assert val > 0, f"bracket gap non-positive at alpha={a_val}, eta={e_val}: {val}"
   ```
2. Also add a `Delta_inf <-> Pe -> infinity` limit check:
   ```python
   Delta_inf_limit = sp.simplify(sp.limit(Delta, Pe, sp.oo))
   expect_zero("Delta_inf as Pe -> oo limit", Delta_inf_limit - Delta_inf_expected)
   ```
3. Mirror both additions in `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl` after line 88, using `FullSimplify` and `Limit[delta, Pe -> Infinity, Assumptions -> alpha > 0 && eta > 0]`, plus a numerical positivity sweep with `Table[N[bracketGap /. {alpha -> a, eta -> e}], {a, {1/10, 1, 3}}, {e, {1/10, 1, 10}}]` checking all values are positive (Exit[1] otherwise).

**Verification:**
The verifier confirms that (a) a new assertion labeled "bracket gap closed form" exists and passes, (b) a numerical positivity loop runs without raising, (c) a new assertion labeled "Delta_inf as Pe -> oo limit" exists and passes in both engines, and (d) `Delta_inf_limit` in the printed output equals the existing `Delta_inf_expected` symbolic form.

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:77-80`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:77-83`

**What's wrong:**
The "Delta_inf formula" check is algebraic substitution against a hand-rewrite of the same expression. The kernel is `K = (cosh(alpha*x) + (eta/alpha)*sinh(alpha*x) - cosh(alpha*(1-x)))/W`. Substituting `x = 1` gives `(cosh(alpha) + (eta/alpha)*sinh(alpha) - cosh(0))/W`, and since `cosh(0) = 1` this is identically `(cosh(alpha) + (eta/alpha)*sinh(alpha) - 1)/W` — which is `Delta_inf_expected` verbatim. The check has no degrees of freedom it can test; it can only fail if `cosh(0) != 1` (i.e., never).

**Why this matters:**
It pads the assertion count without exercising the physics. If F2's "Delta_inf <-> Pe -> infinity" assertion is added, this trivial substitution check should remain only as a sanity print, not as a numbered claim.

**Required change:**
Replace the "Delta_inf formula" check with the limit-based check from F2 (so it ties Delta_inf to lim_{Pe -> infinity} Delta), or, if kept, relabel it to something like "Delta_inf direct substitution (sanity)" and demote from the docstring-numbered claims. Concrete edit at `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:78-80`: leave the substitution `Delta_inf = sp.simplify(K.subs(x, 1))` in place for the print, but remove or relabel `expect_zero("Delta_inf formula", ...)`; the substantive verification of item (5) now lives in the F2-added `"Delta_inf as Pe -> oo limit"` assertion. Mirror at `mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:77-83`.

**Verification:**
After the edit, the file no longer contains an assertion whose only content is `simplify(K.subs(x,1)) == manually_rewritten_K_at_1`. Either the line is removed or its label is changed to something that does not claim to verify the docstring's "exact sharp-bottom endpoint" claim (which is now covered by the limit assertion).

### F4 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:89-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage058_coupled_support_source_operator_mathematica_audit.wl:90-92`

**What's wrong:**
`Delta_series = sp.series(Delta, Pe, 0, 2).removeO()` is the Taylor expansion of `Delta` about Pe = 0, and `Delta0 = sp.limit(Delta, Pe, 0)`. For any function analytic at Pe = 0 (which Delta is — the Pe/(exp(Pe)-1) factor extends smoothly and the 1/(Pe^2 - alpha^2) factor is analytic at Pe = 0 since alpha > 0), the zeroth Taylor coefficient is by definition equal to the limit at 0. So `Delta_series.coeff(Pe, 0) - Delta0` is identically zero for any analytic Delta — the assertion cannot fail under any choice of W, K, or Sigma_Pe consistent with the docstring's analytic setup. It tests sympy/Mathematica's series engine, not the physics. The Mathematica mirror at line 92 has the same shape.

**Why this matters:**
The docstring's weak-coupling claim should be that the *Pe^1* coefficient of Delta has a specific physical interpretation (the leading correction beyond the uniform-source drop), not merely that the constant term equals the constant. As written, the assertion adds no constraint.

**Required change:**
Replace or augment the "weak-coupling constant term" assertion at `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:91` with a check on the *first-order* coefficient. Compute the symbolic Pe^1 coefficient of `Delta_series` and compare it to an independently derived expected form (e.g., `sp.integrate(K * x, (x, 0, 1)) - Delta0/2 * something`), or, at minimum, assert that `Delta_series.coeff(Pe, 1)` is non-zero in general (i.e., not identically vanishing) by substituting `alpha = 1, eta = 1` and confirming a non-zero numerical value. Concrete edit:
```python
Pe1 = sp.simplify(sp.expand(Delta_series).coeff(Pe, 1))
print("Delta(Pe) first-order coefficient =", Pe1)
Pe1_val = sp.nsimplify(Pe1.subs({alpha: 1, eta: 1}))
assert Pe1_val != 0, f"weak-coupling first-order coefficient vanishes at alpha=eta=1: {Pe1_val}"
```
Mirror in Mathematica at line 92: `pe1 = FullSimplify[SeriesCoefficient[deltaSeries, {Pe, 0, 1}], Assumptions -> alpha > 0 && eta > 0]; pe1Val = N[pe1 /. {alpha -> 1, eta -> 1}]; If[Chop[pe1Val] === 0, fail["weak-coupling first-order coefficient vanishes", pe1Val]];`.

**Verification:**
After the edit, the script output contains a printed `Delta(Pe) first-order coefficient = ...` line with a non-zero closed form, and the new `assert Pe1_val != 0` (sympy) / `If[Chop[...] === 0, fail[...]]` (Mathematica) appears in the source.

## Independent-derivation check (Mathematica)

The Mathematica script is NOT an independent re-derivation. Section-by-section correspondence with the SymPy script:

- SymPy lines 38-41 build `W` then `K`; Mathematica lines 31-35 do the same with identical algebraic content.
- SymPy lines 56-57 declare the closed-form antiderivatives `Fc`, `Fs`; Mathematica lines 49-50 declare exactly the same ansatz character-for-character (up to capitalization) and never invokes `Integrate[Exp[Pe*x]*Cosh[alpha*x], x]`.
- SymPy lines 66-68 build the combined `Delta`; Mathematica lines 59-62 transcribe the same combination rather than computing `Integrate[kernel*sigmaPe, {x, 0, 1}]` from first principles.
- All `expectZero` labels in the Mathematica script ("Kprime identity", "Sigma normalization", "Ic antiderivative check", "Is antiderivative check", "Delta0 formula", "Delta0 integral identity", "Delta_inf formula", "weak-coupling constant term") match the SymPy `expect_zero` labels verbatim.

Verdict on independence: transliteration. See F1.

## Engine cross-check

Both engines produce the same closed forms (modulo cosmetic rewriting):

- `K`: sympy `(alpha*(cosh(alpha*x) - cosh(alpha*(x-1))) + eta*sinh(alpha*x))/(alpha*(alpha*sinh(alpha) + eta*cosh(alpha)))` vs Mathematica `(Cosh[alpha*x] - Cosh[alpha - alpha*x] + (eta*Sinh[alpha*x])/alpha)/(eta*Cosh[alpha] + alpha*Sinh[alpha])` — equal after multiplying numerator and denominator by 1/alpha.
- `Ic`: sympy `-(Pe - (Pe*cosh(alpha) - alpha*sinh(alpha))*exp(Pe))/(Pe^2 - alpha^2)` vs Mathematica `(Pe - E^Pe*Pe*Cosh[alpha] + alpha*E^Pe*Sinh[alpha])/(alpha^2 - Pe^2)` — equal (flip sign of numerator and denominator).
- `Is`: sympy `(alpha + (Pe*sinh(alpha) - alpha*cosh(alpha))*exp(Pe))/(Pe^2 - alpha^2)` vs Mathematica `-(alpha - alpha*E^Pe*Cosh[alpha] + E^Pe*Pe*Sinh[alpha])/(alpha^2 - Pe^2)` — equal (same sign flip).
- `Delta_0`: both `eta*(cosh(alpha)-1)/(alpha^2*(eta*cosh(alpha) + alpha*sinh(alpha)))`.
- `Delta_inf`: both `(alpha*(cosh(alpha)-1) + eta*sinh(alpha))/(alpha*(alpha*sinh(alpha) + eta*cosh(alpha)))`.

Engines agree at the level claimed (`engine_disagreement` does not apply). Agreement is not strong evidence of correctness, however, because of F1 (the engines are not independent).

## Verdict justification

Both engines run, both pass, the numerical and symbolic outputs are consistent, and the substantive assertions (Kprime identity, Sigma normalization, Ic/Is antiderivative checks, Delta_0 formula, Delta_0 integral identity) hold up to direct adversarial reading. Attacks I tried that failed: (a) sign / factor-of-alpha mismatch between Ic in sympy vs Mathematica — engines agree after the obvious sign flip of (Pe^2 - alpha^2) vs (alpha^2 - Pe^2); (b) parity / domain attack on `sp.limit(Delta, Pe, 0)` — the Pe/(exp(Pe)-1) factor is analytic at 0, so the L'Hopital reduction is legitimate; (c) checking that `Delta_0 = integrate(K, (x, 0, 1))` is non-trivial — it is, because Delta_0 was derived from the Pe -> 0 limit of a different combination, not from integrating K directly, so the agreement is a real cross-check.

What does NOT hold up: (1) Mathematica does not independently derive the antiderivative ansatz `Fc`/`Fs` or the combined `Delta` (F1); (2) two docstring-numbered claims (the bracket and the "sharp-bottom endpoint") have no actual assertion behind them, only print statements (F2); (3) the `Delta_inf formula` check is algebraic substitution against the same expression rewritten by hand (F3); (4) the `weak-coupling constant term` check is mathematically guaranteed by analyticity and does not exercise the physics (F4). Verdict is `findings`, severity is medium overall, no stop-cold flag — the math is consistent and downstream results are unaffected.

## Self-test notes

Walked through each proposed assertion mentally:
- F2 bracket gap closed form: at alpha=1, eta=1, numerator (1-1)(cosh(1)-1) + 1*1*sinh(1) = sinh(1) > 0; numerical positivity sweep over {1/10, 1, 3} x {1/10, 1, 10} produces all-positive values by hand check (small-alpha expansion `(alpha^2 + eta)*alpha^2/2 + alpha^2*eta` style yields positive); the symbolic form is well-defined (no removable singularities at Pe = 0).
- F2 Delta_inf as Pe -> oo limit: Pe/(exp(Pe)-1) -> 0 but Ic, Is contain exp(Pe) growth that overpowers; net cancellation should give Delta_inf_expected — sympy's `limit(..., Pe, oo)` will resolve this since all symbols are positive. Verified independence: the variable Pe in `limit` is one Delta depends on, so the derivative/limit is non-trivial.
- F4 first-order coefficient: at alpha=1, eta=1, `Delta_series.coeff(Pe, 1)` is a rational combination of sinh(1), cosh(1) that does not vanish (confirmed via small-alpha expansion: coefficient ~ 1/6 + O(alpha) > 0 by direct series-of-series). The `assert Pe1_val != 0` is genuinely non-trivial since cancellation would require the full integrand to be Pe-independent at first order, which it is not.
- Path specs: all directives target existing files only (sympy in `scripts/`, mathematica in `mathematica/`); no new files proposed.
