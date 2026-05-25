---
unit_id: 079
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 079 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.txt`

## What the script claims to verify

The scripts construct the Family-1 amplitude `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)` from `kappa_F1 = 12321/5`, `eta_F1 = 37`, and the Robin root `y_F1` of `y*tan(y) = eta_F1` near `y=1.53`. They then define the closed-form transport-bias map `Omega(Pe) = pi*Pe*(2*Pe*exp(Pe) + pi) / ((4*Pe^2 + pi^2)*(exp(Pe) - 1))` and the demand variable `zeta_F1(Pe) = A_F1 * Omega(Pe)^2`. The substantive claims verified are: (i) `Omega(0+) = 1`, (ii) `Omega(inf) = pi/2`, (iii) the small-Pe linearization `zeta_F1(Pe) ≈ A_F1 * (1 + ((4-pi)/pi)*Pe) + O(Pe^2)`, and (iv) the resulting bracket `A_F1 <= zeta_req <= A_F1 * pi^2/4` for which a constructive `Pe_req` exists. The Mathematica script asserts the same Omega limits and series coefficients, but reduces the bracket-endpoint checks to numeric comparisons against hardcoded targets copied from the SymPy output.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 52 | `expect_small("Omega(0+) - 1", Omega0 - 1)` (symbolic limit) | yes |
| A2 | sympy | 53 | `expect_small("Omega(inf) - pi/2", Omega_inf - pi/2)` (symbolic limit) | yes |
| A3 | sympy | 62 | `expect_small("zeta_F1(0+) - A_F1", zeta0 - A_F1)` (symbolic limit of A_F1*Omega^2) | yes |
| A4 | sympy | 63 | `expect_small("zeta_F1(inf) - A_F1 pi^2/4", zeta_inf - A_F1 * pi^2/4)` | yes |
| A5 | sympy | 69 | `expect_small("small-Pe coefficient check", zeta_series - A_F1*(1+((4-pi)/pi)*Pe))` | yes |
| A6 | mathematica | 46 | `expectApprox["Robin residual", yF1*Tan[yF1] - etaF1, 0, 1e-28]` | yes (residual of nonlinear root, substantive) |
| A7 | mathematica | 58 | `expectZero["Omega(0+) - 1", omega0 - 1]` (symbolic limit) | yes |
| A8 | mathematica | 59 | `expectZero["Omega(inf) - Pi/2", omegaInf - Pi/2]` (symbolic limit) | yes |
| A9 | mathematica | 67 | `expectApprox["zeta_F1(0+) numeric check", zeta0, 1.00005192880219492404747131934, 1e-14]` — `zeta0 := N[aF1*omega0^2, 50]` with `omega0=1`, so this reduces to `aF1 == hardcoded_target`. Target is copied from the SymPy `.txt` (line 24: `zeta_F1(0+) = 1.00005192880219492404747131934`). | no — tautological in zeta_F1 structure, hardcoded target |
| A10 | mathematica | 68 | `expectApprox["zeta_max^(F1) numeric check", zetaInf, 2.46752922945601123498982913352, 1e-14]` — `zetaInf := N[aF1*omegaInf^2, 50]` with `omegaInf=Pi/2`, so this reduces to `aF1*Pi^2/4 == hardcoded_target`. Target is copied from the SymPy `.txt` (line 25: `zeta_F1(inf) = 2.46752922945601123498982913352`). | no — hardcoded target from other engine; check does not independently exercise A_F1 derivation |
| A11 | mathematica | 73 | `expectApprox["small-Pe constant coefficient check", Coeff[seriesDiff,pe,0], 0, 1e-28]` | yes |
| A12 | mathematica | 74 | `expectApprox["small-Pe linear coefficient check", Coeff[seriesDiff,pe,1], 0, 1e-28]` | yes |

## Findings

### F1 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl:67-68`

**What's wrong:**
The two bracket-endpoint checks for `zeta_F1` are written as approximate numeric comparisons against decimal constants:
```
expectApprox["zeta_F1(0+) numeric check", zeta0, ToExpression["1.00005192880219492404747131934`30"], 10^-14];
expectApprox["zeta_max^(F1) numeric check", zetaInf, ToExpression["2.46752922945601123498982913352`30"], 10^-14];
```
These constants are not derived anywhere in the `.wl` script. They match byte-for-byte the SymPy script's own runtime output (`scripts/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.txt` lines 24-25). Mathematica is therefore "verifying" its own `zeta0`/`zetaInf` against numbers that were copied from the SymPy engine — this is not an independent check, it is round-tripping the SymPy result through Mathematica.

**Why this matters:**
The audit's stated purpose is dual-engine cross-verification. If the SymPy engine were wrong by, say, 5e-15 at digit 16 (well inside the `1e-14` tolerance), the Mathematica "check" would silently confirm the wrong number rather than catch it. The check provides no information beyond "Mathematica's `aF1` agrees with SymPy's `A_F1` to 14 digits", which the SymPy script's existence already guarantees by construction.

**Required change:**
Replace lines 67-68 with symbolic identity checks anchored in the engine's own derivation. The substantive claims are `zeta_F1(0+) = A_F1` (the lower bracket endpoint) and `zeta_F1(inf) = A_F1 * Pi^2/4` (the upper bracket endpoint). Use symbolic `Limit` of `aF1*omega^2` (rather than `aF1*omega0^2` with `omega0` already pre-evaluated) to exercise the same derivation path SymPy uses in A3/A4.

**Verification:**
After fix, lines 67-68 should contain `expectZero` calls referencing `Limit[aF1*omega^2, pe -> 0, Direction -> "FromAbove"] - aF1` and `Limit[aF1*omega^2, pe -> Infinity] - aF1*Pi^2/4`, both with no decimal constants. Output should contain `PASS: zeta_F1(0+) - A_F1` and `PASS: zeta_F1(inf) - A_F1 Pi^2/4` (or similarly named) and no decimal targets near 1.00005... or 2.46752....

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl:61-62,67-68`

**What's wrong:**
Lines 61-62 define `zeta0 := N[aF1*omega0^2, 50]` and `zetaInf := N[aF1*omegaInf^2, 50]`. By the previous lines (58-59), `omega0 = 1` and `omegaInf = Pi/2` are already established as exact symbolic values. So `zeta0` is identically `aF1` and `zetaInf` is identically `aF1*Pi^2/4` by construction at the line they are defined. The checks at lines 67-68 then compare these to numbers; but the only thing being tested is that `aF1` (which is itself a fixed numeric value from line 40) equals a pre-baked target. There is no independent derivation of the limit of `aF1*omega^2` — the limit was implicitly substituted before the assertion, so the assertion cannot detect an error in the limit calculation itself.

**Why this matters:**
The intent of A3/A4 (in SymPy) is to verify that `sp.limit(A_F1*Omega^2, Pe, 0)` symbolically reduces to `A_F1` and `sp.limit(A_F1*Omega^2, Pe, oo)` reduces to `A_F1*pi^2/4`. The Mathematica script bypasses the symbolic limit by substituting the pre-computed `omega0` and `omegaInf` before the multiplication, so the check no longer exercises the symbolic-limit machinery at all. A bug in Mathematica's `Limit` applied to `aF1*omega^2` (e.g., a wrong branch at infinity) would be invisible.

**Required change:**
Same as F1: replace the two checks with `expectZero` of `Limit[aF1*omega^2, pe -> 0, Direction -> "FromAbove"] - aF1` and `Limit[aF1*omega^2, pe -> Infinity] - aF1*Pi^2/4`. This forces Mathematica to take the symbolic limit of the *product* (not of `omega` alone), matching what SymPy's A3/A4 do. The intermediate `zeta0`/`zetaInf` numeric variables on lines 61-62 can be retained only for the `Print` statements at 65-66, or deleted entirely.

**Verification:**
The new assertions must contain `Limit[aF1*omega^2, ...]` (not `aF1*omega0^2` or `aF1*omegaInf^2`) in the residual expression. Output should show the residuals reduce to exact zero (or sub-1e-40 numerically), not just 1e-15.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl:32-74`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.py:28-69`

**What's wrong:**
The Mathematica script reproduces the SymPy script's exact algebraic choreography step-for-step, with identical numeric constants and identical closed-form expressions, rather than independently re-deriving `Omega(Pe)` or `A_F1` from a physical premise. Corresponding sections:

SymPy lines 33-37:
```python
kappa_F1 = sp.Rational(12321, 5)
eta_F1 = sp.Integer(37)
y_F1 = sp.nsolve(y * sp.tan(y) - eta_F1, 1.53, tol=1e-34, maxsteps=100)
A_F1 = (kappa_F1 + sp.pi**2 / 4) / (kappa_F1 + y_F1**2)
```
Mathematica lines 37-40:
```mathematica
kappaF1 = 12321/5;
etaF1 = 37;
yF1 = y /. FindRoot[y*Tan[y] == etaF1, {y, 1.53}, WorkingPrecision -> 80, AccuracyGoal -> 30, PrecisionGoal -> 30];
aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yF1^2), 50];
```

SymPy line 45:
```python
Omega = sp.simplify(sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi) / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1)))
```
Mathematica lines 48-51:
```mathematica
omega = FullSimplify[
  Pi*pe*(2*pe*Exp[pe] + Pi)/((4*pe^2 + Pi^2)*(Exp[pe] - 1)),
  Assumptions -> $Assumptions
];
```

SymPy lines 67-69 (expected series):
```python
expected_series = sp.expand(A_F1 * (1 + ((4 - sp.pi) / sp.pi) * Pe))
```
Mathematica line 71:
```mathematica
expectedSeries = FullSimplify[aF1*(1 + ((4 - Pi)/Pi)*pe), Assumptions -> $Assumptions];
```

These are not independent derivations; the Mathematica file is a syntactic port of the SymPy file. Neither script derives `Omega(Pe)` from an underlying PDE / boundary-value problem — both take the closed-form expression as given. A second-engine audit whose only job is to compute `Limit` of the same closed-form expression at the same two points adds only that Mathematica's `Limit` and `Series` primitives agree with mpmath's, not that the closed form itself is correct.

**Why this matters:**
The second-engine policy exists so that an algebraic mistake in transcribing the underlying physics into a closed form is caught by independent re-derivation. Here, if the closed form for `Omega(Pe)` were wrong (e.g., a factor of 2 missing in the numerator, or a swap of `pi` and `Pe`), both engines would compute the same limits and series of the same wrong expression and both would "PASS". The same is true for the choice of `(kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)` for `A_F1`.

**Required change:**
This is the hardest finding to resolve mechanically because the script does not contain any earlier-stage derivation of `Omega(Pe)` or `A_F1` from a PDE. A minimal independent corroboration the Mathematica script could perform without inventing new scope:

1. **Independent integral form for `Omega(Pe)`**: in the convection-diffusion / Robin family, the cell-average bias factor admits an integral representation. A check that `Omega(Pe)` equals `(Pi^2/2) * Integrate[expression in pe, ...]` over a stated interval would corroborate the closed form from a different angle. However, the script gives no such integral, so this would amount to scope expansion and is out-of-bounds for the directive.

2. **Independent computation of `Omega` from a Fourier series / eigenfunction expansion**: same caveat.

Given the rule "do not propose new features or scope extensions", the actionable correction is narrower: add a single substantive *independent* arithmetic identity for `Omega(Pe)` that Mathematica can check without echoing the SymPy steps. For example, the Robin / convection-diffusion `Omega` satisfies the algebraic functional identity
`Omega(Pe) * (4*Pe^2 + Pi^2) * (Exp[Pe] - 1) == Pi*Pe*(2*Pe*Exp[Pe] + Pi)`
which is the closed form *re-stated as a product equation*. A `FullSimplify` of that residual being zero is mechanical and not independent. A more substantive option: verify the **derivative identity** `D[Omega[Pe], Pe]` evaluated at `Pe -> 0` equals `(4-Pi)/(2*Pi)` (the slope implied by the small-Pe series). Computing `Limit[D[omega, pe], pe -> 0, Direction -> "FromAbove"]` independently in Mathematica and comparing to `(4-Pi)/(2*Pi)` is a non-trivial check on the closed form that does not appear in the SymPy script and forces Mathematica's symbolic differentiation engine to corroborate the SymPy series result through a different code path. (The SymPy series in A5 is computed by `sp.series`; an explicit symbolic derivative in Mathematica is a different computation that lands on the same number iff the closed form is right.)

Therefore: add an `expectZero["Omega'(0+) - (4-Pi)/(2 Pi)", Limit[D[omega, pe], pe -> 0, Direction -> "FromAbove"] - (4 - Pi)/(2*Pi)]` after the existing limit checks.

**Verification:**
After fix, the Mathematica `.wl` should contain a `D[omega, pe]` call (not present today) and an `expectZero` line with target `(4-Pi)/(2*Pi)` evaluated as a Limit at `pe -> 0`. Output should show `PASS: Omega'(0+) - (4-Pi)/(2 Pi)`. This forces Mathematica's symbolic differentiation to corroborate the closed form independently of SymPy's `series` engine.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script. See F3 for the corresponding-section quotes. The closed form for `Omega(Pe)` and the structure of `A_F1` are introduced as fiat in both scripts without independent derivation; the only difference between the two files is engine syntax and the substitution of hardcoded numeric targets for SymPy's symbolic-limit checks (F1, F2).

## Engine cross-check

Both engines agree at the level claimed:
- `Omega(0+) = 1` — SymPy: exact zero; Mathematica: exact zero.
- `Omega(inf) = Pi/2` — SymPy: exact zero; Mathematica: exact zero.
- `A_F1` — SymPy at 30 digits: `1.00005192880219520364707725466`. Mathematica at 50 digits: `1.00005192880219532865933408371...`. These agree to ~15 digits, which is consistent with SymPy's `nsolve(... tol=1e-34)` actually delivering only ~15-digit precision (mpmath default precision was not raised; the SymPy `Robin residual = -1.42e-14` confirms this). Within the assertion tolerances (`1e-12`), this is fine.
- Series coefficients — both reproduce `(4-Pi)/Pi` linear coefficient to high precision.

No `engine_disagreement` finding is warranted.

## Verdict justification

The findings are real but bounded. The SymPy script's substantive checks (A1-A5) all hold up: they exercise symbolic limits of the `Omega(Pe)` closed form, compare zeta_F1's endpoints to symbolic targets derived in-script, and verify the small-Pe linearization. The Mathematica script's A6-A8 and A11-A12 are also substantive. However, A9 and A10 (lines 67-68) are tautological numeric round-trips of SymPy-produced constants, and the entire `.wl` is a syntactic port of the `.py` rather than an independent derivation. F1 and F2 are easy mechanical fixes (swap two `expectApprox` for `expectZero`); F3 requires adding one substantive independent check that the closed form for `Omega(Pe)` corresponds to the claimed small-Pe slope, computed through Mathematica's symbolic-differentiation path rather than SymPy's series path.

Attacks tried that did not break the math:
- Verified the small-Pe series by hand: `Omega(Pe) ≈ 1 + ((4-pi)/(2*pi))*Pe + O(Pe^2)`, so `Omega^2 ≈ 1 + ((4-pi)/pi)*Pe + O(Pe^2)`. The script's `expected_series = A_F1*(1 + ((4-pi)/pi)*Pe)` is correct.
- Verified the limit at Pe → infinity: numerator `~2*pi*Pe^2*exp(Pe)` for large Pe; denominator `~4*Pe^2*exp(Pe)`. Ratio `→ pi/2`. Correct.
- Verified the limit at Pe → 0 by series: numerator `~pi^2*Pe`, denominator `~pi^2*Pe` (leading order). Ratio `→ 1`. Correct.
- The `Limit` warnings in the Mathematica output (`Assumptions that involve the limit variable are ignored`) are benign — they refer to the `pe > 0` assumption being dropped when `pe` is the limit variable, which does not invalidate the limit.
- Robin residual: SymPy `~1.4e-14` vs Mathematica `~2.3e-55`. Both well inside their tolerances; the SymPy precision shortfall is a non-blocking note (the assertion tolerance is `1e-12`, so it passes; the `tol=1e-34` in `sp.nsolve` is misleading but does not cause a failed claim within this stage).

No `UNFIXABLE` or `CRITICAL_DOWNSTREAM` verdict: F1/F2 affect only how the Mathematica script tests, not the values themselves; F3 asks Mathematica to add one corroborating check, not to change any value. No downstream unit's quoted value would change.

## Self-test notes

- Variable independence: the proposed `D[omega, pe]` in F3 — `omega` is a non-trivial function of `pe` (it depends on `pe` through `Exp[pe]`, `pe^2`, and a linear `pe` factor), so the derivative is non-trivially nonzero, and `Limit[D[omega, pe], pe -> 0, Direction -> "FromAbove"]` is a real symbolic operation, not a trivial zero. By the manual series expansion `Omega(Pe) = 1 + ((4-pi)/(2*pi))*Pe + O(Pe^2)`, the slope at 0+ is `(4-pi)/(2*pi)`, which is nonzero (`≈ 0.1366`).
- Trivial-case pre-check for F1/F2: substituting the symbolic limits, `Limit[aF1*omega^2, pe -> 0]` should give `aF1*1 = aF1`, so the residual `... - aF1` is exactly 0. `Limit[aF1*omega^2, pe -> Infinity]` gives `aF1*(Pi/2)^2 = aF1*Pi^2/4`, residual exactly 0. Both `expectZero` will pass.
- Path specifications: the Mathematica script lives at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl`. All directive targets use this absolute path.
