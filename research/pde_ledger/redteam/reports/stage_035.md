---
unit_id: 035
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T21:52:09Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 035 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.txt`

## What the script claims to verify

The scripts claim to verify, for the moving-throat PDE D/N overlap constants `kappa0^2 = 8/pi^2`, `kappa1^2 = 16/(9*pi^2)`: (1) the dimensionless shape function `F(xi,delta) = N(x)*A/(beta0*kappa0^2)` equals the closed form `(9*delta + 11*xi)^4 / (81*(1 - xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2)`; (2) `dF/dxi` equals a closed form with the manifestly-positive numerator polynomial `(9*delta+11*xi)^3*(81*delta^3 + 72*delta^2 + 189*delta^2*xi + 297*delta*xi^2 + 121*xi^3)`, and the onset/softening limits `F(0,delta) = 1` and `lim_{xi->1^-} (1-xi)*F = (9*delta+11)^4 / (81*(9*delta^2+18*delta+11)^2)`; (3) the required total loading `alpha_req(xi,delta) = 9*pi^2*A*xi*(xi+delta) / (8*(9*delta+11*xi))` and its critical limit `alpha_crit = 9*pi^2*A*(1+delta)/(8*(11+9*delta))`; (4) near-onset series expansions of `F` and `alpha_req` through O(xi^2). The Mathematica script additionally checks `R_target == pi^2*A*NQ/(8*beta0)`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 58 | `simplify(F - F_target) == 0` | yes |
| A2 | sympy | 73 | `simplify(dF - dF_target) == 0` | partial (verifies rewritten form; does not assert positivity) |
| A3 | sympy | 74 | `simplify(F_target.subs(xi,0) - 1) == 0` | yes |
| A4 | sympy | 78 | `simplify(soft_const - soft_const_target) == 0` | yes |
| A5 | sympy | 82 | `simplify(soft_scaled_series - soft_const_target) == 0` | yes |
| A6 | sympy | 89 | `simplify(alpha_req - alpha_req_target) == 0` | yes |
| A7 | sympy | 94 | `simplify(alpha_crit - alpha_crit_target) == 0` | yes |
| A8 | sympy | 105 | `simplify(F_series - F_series_target) == 0` | yes |
| A9 | sympy | 110 | `simplify(alpha_series - alpha_series_target) == 0` | yes |
| A10 | mathematica | 61 | `expectZero[F - fTarget]` | yes |
| A11 | mathematica | 62 | `expectZero[rTarget - rTargetClosed]` | yes (algebraic identity, but uses two independently typed expressions for `R_target`) |
| A12 | mathematica | 71 | `expectZero[dF - dFTarget]` | partial (same partial as A2) |
| A13 | mathematica | 72 | `expectZero[(fTarget /. xi -> 0) - 1]` | yes |
| A14 | mathematica | 83 | `expectZero[softConst - softConstTarget]` | yes |
| A15 | mathematica | 90 | `expectZero[softScaledSeries - softConstTarget]` | yes |
| A16 | mathematica | 105 | `expectZero[alphaReq - alphaReqTarget]` | yes |
| A17 | mathematica | 106 | `expectZero[alphaCrit - alphaCritTarget]` | yes |
| A18 | mathematica | 111 | `expectZero[gBReqSqOverVarpi2 - (alphaReqTarget - alphaMix)]` | **no — tautological** |
| A19 | mathematica | 123 | `expectZero[fSeries - fSeriesTarget]` | yes |
| A20 | mathematica | 124 | `expectZero[alphaSeries - alphaSeriesTarget]` | yes |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:108-111`

**What's wrong:**

The Mathematica script defines, on lines 108-109:

```
alphaMix = FullSimplify[Chi^2/(OmegaU^2*Delta0), Assumptions -> $Assumptions];
gBReqSqOverVarpi2 = FullSimplify[alphaReqTarget - alphaMix, Assumptions -> $Assumptions];
```

Then on line 111 asserts:

```
expectZero["support split residual", gBReqSqOverVarpi2 - (alphaReqTarget - alphaMix)];
```

The argument to `expectZero` is `FullSimplify[alphaReqTarget - alphaMix] - (alphaReqTarget - alphaMix)`. After `FullSimplify` inside `expectZero` is applied, this is mathematically guaranteed to equal zero by construction: it is `E - E` where `E = alphaReqTarget - alphaMix`. The check cannot fail regardless of what `alphaReqTarget` or `alphaMix` actually equal — so it verifies nothing about the support-coupling split.

The corresponding SymPy script (line 96-99) computes the same quantity but does not include this assertion at all (it only prints it), confirming the assertion was added in the Mathematica port to balance the output but did not encode a substantive check.

**Why this matters:**

`expectZero` in the report counts as a verification step but adds zero information. If the intent was to verify that `alpha_req_target - alpha_mix` matches a particular closed form (e.g., a `g_B,req^2 / varpi^2` formula stated elsewhere), that target is not present in either script. As written, the line decorates the output with a false PASS that gives the reader misplaced confidence.

**Required change:**

Remove the tautological assertion at line 111. Keep the surrounding print so the value is still emitted. Delete only the `expectZero` line; do not add a substitute claim because no `g_B,req^2/varpi^2` closed-form target exists in the SymPy script either.

**Verification:**

After the patch, the Mathematica output file should no longer contain a `support split residual = 0 / PASS: support split residual` line pair (currently lines 37-38 of `moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.txt`). The `Print["g_B,req^2 / varpi^2 = ...", ...]` line above remains. Script must still exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:34-54,64-69`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py:33-54,65-70`

**What's wrong:**

The Mathematica script is a line-by-line port of the SymPy script's algebraic choreography, not an independent re-derivation. Corresponding excerpts:

SymPy (lines 33-50):
```
kappa0_sq = sp.Rational(8) / sp.pi**2
kappa1_sq = sp.Rational(16) / (9 * sp.pi**2)
...
x = sp.simplify(A * xi)
DeltaK_sub = sp.simplify(A * delta)
alpha_x = sp.simplify(x * (x + DeltaK_sub) / (kappa0_sq * (x + DeltaK_sub) + kappa1_sq * x))
N_x = sp.simplify(
    beta0 * (kappa0_sq * (x + DeltaK_sub) + kappa1_sq * x) ** 4
    / (kappa0_sq * (A - x) * (kappa0_sq * (x + DeltaK_sub) ** 2 + kappa1_sq * x**2) ** 2)
)
```

Mathematica (lines 34-48):
```
kappa0Sq = 8/Pi^2;
kappa1Sq = 16/(9*Pi^2);
...
x = FullSimplify[A*xi, ...];
deltaKSub = FullSimplify[A*delta, ...];
alphaX = FullSimplify[x*(x + deltaKSub)/(kappa0Sq*(x + deltaKSub) + kappa1Sq*x), ...];
nX = FullSimplify[beta0*(kappa0Sq*(x + deltaKSub) + kappa1Sq*x)^4/
    (kappa0Sq*(A - x)*(kappa0Sq*(x + deltaKSub)^2 + kappa1Sq*x^2)^2), ...];
```

This is the same expression written in the other engine's syntax — identical variable choreography, identical bracketing, identical intermediate symbols. The same pattern repeats for the closed-form target polynomials:

SymPy (lines 66-70):
```
dF_target = sp.simplify(
    (9 * delta + 11 * xi) ** 3
    * (81 * delta**3 + 72 * delta**2 + 189 * delta**2 * xi + 297 * delta * xi**2 + 121 * xi**3)
    / (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)
)
```

Mathematica (lines 64-69):
```
dFTarget = FullSimplify[
  (9*delta + 11*xi)^3*(81*delta^3 + 72*delta^2 + 189*delta^2*xi + 297*delta*xi^2 + 121*xi^3)/
    (81*(1 - xi)^2*(9*delta^2 + 18*delta*xi + 11*xi^2)^3),
  ...
];
```

Identical literal polynomial coefficients (81, 72, 189, 297, 121 in the numerator; 81, 18, 11 in the denominator) appear in both. Because the same numeric target is hard-coded as `dFTarget` in both engines, the second engine cannot catch an error in those literal coefficients — both would simply confirm `D[fTarget, xi]` equals the hard-coded polynomial. The independent-derivation requirement is that the two engines start from the same physical premises (the D/N overlap constants and the definition of `N_x`) and **independently** arrive at the closed-form polynomial. Here the closed-form polynomials are dropped into both scripts as the same literal.

The same observation applies to `fTarget`, `softConstTarget`, `alphaReqTarget`, `alphaCritTarget`, `fSeriesTarget`, `alphaSeriesTarget` — all use identical literal expansions in both scripts.

**Why this matters:**

If the closed-form polynomial coefficients in the SymPy target (e.g., the `121*xi^3` coefficient in `dFTarget`) were wrong, the Mathematica script would not detect it: the Mathematica script uses the same literal `121*xi^3` and would `FullSimplify` to the same wrong polynomial. The second engine provides no cross-check on the closed forms themselves — only that two algebra systems agree about a tautology. The transliteration violates the second-engine policy: both engines must derive the result independently from the physical premises.

**Required change:**

In the Mathematica script, replace the hard-coded `dFTarget`, `softConstTarget`, `alphaCritTarget`, `fSeriesTarget`, `alphaSeriesTarget`, and `fTarget` definitions with **independent** derivations from the engine's own algebra, then compare the engine-derived closed form against the engine-derived `D[f, xi]`/limit/series. Specifically:

1. Lines 51-54 (`fTarget`): Replace with `fTarget = Together[FullSimplify[f, Assumptions -> $Assumptions]];` — i.e., let Mathematica produce its own canonical form of `f` (no literal polynomial drop-in). Then the assertion on line 61 becomes `f - fTarget` which is still trivially zero — so additionally introduce an independent rational-function canonicalization, e.g.:
   ```
   fIndep = Together[(9*delta + 11*xi)^4/(81*(1 - xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2)];
   expectZero["F - independent rational form", Cancel[Together[f - fIndep]]];
   ```
   where `fIndep` is the literal closed form being tested. The point is to keep one literal closed-form for the test (that's the claim under verification), but ensure the path from the physical premises to `f` is independent in each engine.

2. Lines 64-69 (`dFTarget`): Replace with engine-independent derivation:
   ```
   dF = FullSimplify[D[f, xi], Assumptions -> $Assumptions];
   dFIndep = (9*delta + 11*xi)^3*(81*delta^3 + 72*delta^2 + 189*delta^2*xi + 297*delta*xi^2 + 121*xi^3)/
             (81*(1 - xi)^2*(9*delta^2 + 18*delta*xi + 11*xi^2)^3);
   expectZero["dF/dxi - manifestly positive form", Cancel[Together[dF - dFIndep]]];
   ```
   That is, compute `dF` from `f` (Mathematica's own `f`, not the literal `fTarget`), and compare against the claimed closed form `dFIndep`. This way, if SymPy's `dF_target` had wrong coefficients, the Mathematica script would compute `D[f, xi]` independently and the residual would be nonzero.

3. Apply the same pattern to `softConst`, `alphaReqTarget`, `alphaCritTarget`, `fSeriesTarget`, `alphaSeriesTarget`: in each case, compute the LHS (limit, series, etc.) from the Mathematica `f` or `alphaX` directly, and compare against the literal claimed closed form. Do not derive the LHS by `FullSimplify`-ing the same literal closed form.

The principle: the closed-form target may be a literal (it's the claim under test), but the LHS to which it is compared must come from independent engine computation, not from the same literal.

**Verification:**

After the patch, `dF` (or equivalent) must be defined as `D[f, xi]` where `f` is computed via lines 50 of the script (from `nX`), **not** by simplifying the same literal `dFTarget`. The Mathematica output should still show `dF/dxi - manifestly positive form = 0` (the math is correct), but the residual is now computed by Mathematica's independent `D[f, xi]` minus the literal closed form, which would surface any literal-coefficient error. Run the script; it must exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script does NOT independently derive the closed forms; it ports the SymPy script's algebra (same intermediate expressions, same hard-coded polynomial coefficients, same series targets). See F2 for quoted side-by-side excerpts. This is the textbook `mathematica_transliteration` pattern: same variable names (modulo capitalization), same arithmetic structure, identical literal targets.

## Engine cross-check

Both engines produce matching final residuals (all `= 0`) and matching output forms:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `F(xi,delta)` | `-(9δ + 11ξ)^4 / (81(ξ-1)(2ξ² + 9(δ+ξ)²)²)` | `-1/81*(9*delta + 11*xi)^4/((-1 + xi)*(9*delta^2 + 18*delta*xi + 11*xi^2)^2)` (algebraically equal) |
| `dF/dxi residual` | 0 | 0 |
| `softening constant` | `(9δ+11)^4/(81(9δ²+18δ+11)²)` | `(11+9δ)^4/(81(11+9δ(2+δ))²)` (algebraically equal) |
| `alpha_req` | `9π²A ξ(δ+ξ)/(8(9δ+11ξ))` | `(9 A Pi² xi(δ+xi))/(72δ+88xi)` (algebraically equal) |
| `alpha_crit` | `9π²A(δ+1)/(8(9δ+11))` | `9 A(1+δ)Pi²/(88+72δ)` (algebraically equal) |
| series F O(xi²) coeff residual | 0 | 0 |
| series alpha O(xi²) coeff residual | 0 | 0 |

Numerical agreement is consistent. The concern is not numerical disagreement — it is that the agreement is structural (both engines fed the same literal targets), so this cross-check has reduced epistemic weight. See F2.

## Verdict justification

The math itself holds up: every substantive assertion (A1-A17, A19, A20) verifies a claim that is non-trivial and the residuals correctly evaluate to zero. The closed-form polynomials for `F`, `dF`, the softening constant, `alpha_req`, `alpha_crit`, and the near-onset series are internally consistent. However, two issues warrant findings: (1) one Mathematica assertion (A18) is tautological and adds no verification value; (2) the Mathematica script's closed-form targets are hard-coded as the same literals used in the SymPy script, so the two engines do not provide independent cross-checks of those coefficients. I attempted to break the math by: testing parity/sign on the limit at `xi -> 1`, the series expansion through `O(xi^2)` (computed by hand for delta=1 to verify the `28/(27*delta^2)` coefficient — yes, comes from the binomial expansion of the squared denominator at xi=0), the positivity of the `dF` numerator polynomial (coefficients `81, 72, 189, 297, 121` all positive — manifest), and the algebraic factorization of `alpha_x` denominator `(8(9δ+11xi))/(9π²)` after substituting `kappa0^2, kappa1^2`. All hold up. Verdict: `findings`, two non-fatal items.

## Self-test notes

- Variable independence: not applicable (no `sp.diff` over a variable absent from the integrand in this script; all derivatives are `d/dxi` of expressions explicitly containing `xi`).
- Symmetry/parity: not applicable (no unbounded-domain integrals; all checks are algebraic identities or limits/series).
- Trivial-case pre-check: substituted `xi=0` mentally — `F(0,delta) = (9δ)^4/(81*1*(9δ²)^2) = 6561δ^4/(6561δ^4) = 1` ✓. Substituted `delta=1, xi=0` into `dF` numerator polynomial: `81+72+0+0+0 = 153`; matches the explicit expansion of `4*81*(-1) + 9*9 - 44*(-1)*9 = -324+81+396 = 153` ✓ (sanity check on F2's claim that the closed-form polynomial agrees with `D[f, xi]`).
- Path specifications: the directive's targets all refer to existing files in `mathematica/` and `scripts/`; no missing-script findings, so no path-only edits.
