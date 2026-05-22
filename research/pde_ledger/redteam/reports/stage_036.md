---
unit_id: 036
batch: II.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 6
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 036 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt`

## What the script claims to verify

The scripts claim to verify (per docstring): (1) the exact dimensionless support-feasibility function `G(xi,delta) = 9*xi*(xi+delta)/(9*delta+11*xi)`; (2) exact monotonicity of G and its endpoint values `G(0,delta)=0` and `G_max = lim_{xi->1^-} G = 9*(1+delta)/(9*delta+11)`; (3) the exact factorization of the required support loading as `g_{B,req}^2/varpi^2 = (pi^2 A / 8)*(G - M_mix)`; (4) the near-onset asymptotics `G ~ xi - 2*xi^2/(9*delta)`; and they also exhibit witness samples for admissible/inadmissible parameter choices that anchor the support inequality `M_mix <= G` and the R-target test `F(xi,delta) >= 1`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 68-71 | `simplify(gBreq_sq_over_varpi2 - (pi^2 A/8)(G - Mmix_expr)) == 0` | no (tautological by construction) |
| A2 | sympy | 72 | `simplify(R_target - pi^2 A NQ/(8 beta0)) == 0` | no (tautological — R_target was just defined as that) |
| A3 | sympy | 78 | `simplify(dG/dxi - dG_target) == 0` | yes |
| A4 | sympy | 79 | `simplify(G.subs(xi,0)) == 0` | partial (trivial substitution) |
| A5 | sympy | 83 | `simplify(Gmax - 9*(1+delta)/(9*delta+11)) == 0` | yes (exercises symbolic limit) |
| A6 | sympy | 88-91 | `(pi^2 A/8)*(G(xi_req)-Mmix) - gBreq(xi_req) == 0` | no (same identity as A1 with renamed variable) |
| A7 | sympy | 99 | `bool(F_sample >= 1)` | yes (numeric witness) |
| A8 | sympy | 117-120 | `F(1/2,1) - R_target_host == 0` (kappa-based) | yes (non-trivial — uses kappa expansion) |
| A9 | sympy | 122-125 | `bool(Mmix_good <= G_sample)` numeric witness | partial (witness only) |
| A10 | sympy | 126 | `bool(Rational(9,10) < 1)` | no (purely "9/10 < 1") |
| A11 | sympy | 128-131 | `bool(Mmix_bad > G_sample)` numeric witness | partial (witness only) |
| A12 | sympy | 137 | `simplify(G_series - (xi - 2*xi^2/(9*delta))) == 0` | yes (series expansion) |
| B1 | mathematica | 57 | `g - 8*alphaReq/(Pi^2*A) == 0` | no (algebraic identity by construction) |
| B2 | mathematica | 58 | `g - gTarget == 0` | no (g and gTarget both built from the same alphaReq closed form) |
| B3 | mathematica | 59 | `rTarget - Pi^2*A*NQ/(8*beta0) == 0` | no (tautological) |
| B4 | mathematica | 60-63 | `gBReqSqOverVarpi2 - (Pi^2 A/8)(gTarget - mMix) == 0` | no (tautological mirror of A1) |
| B5 | mathematica | 71 | `dG - dGTarget == 0` | yes |
| B6 | mathematica | 72 | `gTarget /. xi->0 == 0` | partial |
| B7 | mathematica | 82 | `gMax - gMaxTarget == 0` | yes |
| B8 | mathematica | 83 | `(Pi^2 A/8)*gMaxTarget - alphaCrit == 0` | no (alphaCrit defined as exactly that) |
| B9 | mathematica | 90-93 | `(Pi^2 A/8)((gTarget /. xi->xiReq) - mMix) - (gBReqSqOverVarpi2/.xi->xiReq) == 0` | no (mirror of A6) |
| B10 | mathematica | 101 | `fSample >= 1` | yes |
| B11 | mathematica | 111-114 | `(fTarget /. ...) - rTargetHost == 0` | yes (numeric kappa witness) |
| B12 | mathematica | 115-119 | `mMixGood <= gSample` numeric | partial |
| B13 | mathematica | 120 | `9/10 < 1` | no |
| B14 | mathematica | 121-125 | `mMixBad > gSample` numeric | partial |
| B15 | mathematica | 131 | `gSeries - (xi - 2*xi^2/(9*delta)) == 0` | yes |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:57,72`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:50,59`

**What's wrong:**
In the SymPy script line 57 defines
`R_target = sp.simplify(NQ * A / (beta0 * (sp.Rational(8) / sp.pi**2)))`
which simplifies to `pi**2 * A * NQ / (8 * beta0)`. Line 72 then asserts
`expect_zero("R_target - pi^2 A NQ/(8 beta0)", R_target - sp.pi**2 * A * NQ / (8 * beta0))`.
This is `X - X = 0` for the same expression. The assertion cannot fail no matter what the physics is.

The Mathematica script mirrors this: line 50 defines `rTarget = FullSimplify[NQ*A/(beta0*(8/Pi^2)), ...]`, line 59 asserts `rTarget - Pi^2*A*NQ/(8*beta0) == 0`.

**Why this matters:**
The script's docstring identifies the `R_target` relation as one of the structural quantities being audited. A line that re-states a definition as a "check" provides zero verification of the physics. It is misleading evidence of coverage.

**Required change:**
Remove the tautological line 72 (SymPy) and line 59 (Mathematica). The `R_target` symbol still appears in the printed output for context. If the script wants a non-trivial test of R_target, that role is already filled by the kappa-based derivation at A8/B11 (lines 104-120 / 102-114).

**Verification:**
On re-run, the lines `R_target - pi^2 A NQ/(8 beta0) = 0` and `PASS: R_target - Pi^2 A NQ/(8 beta0)` must be absent from the saved outputs. All remaining assertions must still pass.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:54,56,60-61,68-71,88-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:41-43,49,51,60-63,90-93`

**What's wrong:**
The SymPy script defines, at lines 53-61,
`F = (9*delta+11*xi)**4 / (81*(1-xi)*(9*delta**2+18*delta*xi+11*xi**2)**2)`,
`G = 9*xi*(xi+delta)/(9*delta+11*xi)`,
`alpha_mix = Chi**2/(OmegaU**2*Delta0)`,
`Mmix_expr = 8*alpha_mix/(pi**2*A)`,
`a_req = 9*pi**2*A*xi*(xi+delta)/(8*(9*delta+11*xi))`,
`gBreq_sq_over_varpi2 = a_req - alpha_mix`.

Inspection: `(pi**2*A/8)*G = (pi**2*A/8)*9*xi*(xi+delta)/(9*delta+11*xi) = a_req` and `(pi**2*A/8)*Mmix_expr = alpha_mix` follow from the literal closed forms used, with no derivation step in between. Therefore the assertion at lines 68-71,
`gBreq_sq_over_varpi2 - (pi**2*A/8)*(G - Mmix_expr) = (a_req - alpha_mix) - (a_req - alpha_mix) = 0`,
is algebraically guaranteed by the way the symbols were introduced; it does not test the claimed "exact factorization of the required support loading through G - M_mix" (docstring item 3) — it restates how the symbols were typed in.

Lines 88-91 re-apply the same identity at `xi -> xi_req`; this is a syntactic substitution and provides no additional information.

The Mathematica script mirrors this construction at lines 41-43 (`alphaReq` written explicitly as `9*Pi^2*A*xi*(xi+delta)/(8*(9*delta+11*xi))`, `gTarget` as `9*xi*(xi+delta)/(9*delta+11*xi)`), and the corresponding "checks" at lines 57, 60-63, 90-93 are tautological in exactly the same way.

**Why this matters:**
The "exact factorization through G - M_mix" is the central support-inequality content of this checkpoint stage. If it is only `X - X = 0` after both sides were written in the same form, then the script is not actually verifying that factorization — it is asserting that a definition equals itself. A future change to either `a_req` or `G`'s closed form would not be caught by this assertion unless the change broke the textual match; small algebraic errors propagated into the closed-form template would pass silently.

**Required change:**
In the SymPy script, derive `a_req` (or equivalently `(pi**2*A/8)*G`) symbolically from a different starting point so the equality is no longer by construction. The script already has the kappa-based derivation machinery for the host sample (lines 104-116); generalize it symbolically. Concretely, introduce symbolic `A_sym, beta0_sym, xi, delta` and define
`x_sym = A_sym*xi`, `deltaK_sym = A_sym*delta`,
`N_sym = beta0_sym*(KAPPA0_SQ*(x_sym+deltaK_sym) + KAPPA1_SQ*x_sym)**4 / (KAPPA0_SQ*(A_sym-x_sym)*(KAPPA0_SQ*(x_sym+deltaK_sym)**2 + KAPPA1_SQ*x_sym**2)**2)`
`R_target_sym = N_sym*A_sym/(beta0_sym*KAPPA0_SQ)`
and then assert `simplify(R_target_sym - F) == 0`. This is a non-trivial identity that depends on the specific values `KAPPA0_SQ = 8/pi**2`, `KAPPA1_SQ = 16/(9*pi**2)` being correct. With that anchor in place, the lines 68-71 and 88-91 assertions, which mechanically follow from the closed forms once F (and hence G) is anchored to the kappa expansion, can remain as bookkeeping confirmations of the factorization.

In the Mathematica script, add the analogous symbolic derivation (do not copy the SymPy variable names — derive it independently from the same physical premises).

**Verification:**
On re-run, the SymPy output must contain a new assertion line of the form `R_target_sym - F = 0` that passes, and the Mathematica output must contain a corresponding `RTargetSym - fTarget == 0` PASS line.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:79,83`

**What's wrong:**
Line 79 defines
`alphaCrit = FullSimplify[9*Pi^2*A*(1 + delta)/(8*(11 + 9*delta)), ...]`
and line 83 asserts
`expectZero["(Pi^2 A / 8) G_max - alpha_crit", (Pi^2*A/8)*gMaxTarget - alphaCrit]`.
Since `gMaxTarget = 9*(1+delta)/(9*delta+11)` (line 78), `(Pi^2*A/8)*gMaxTarget = 9*Pi^2*A*(1+delta)/(8*(9*delta+11)) = alphaCrit` is `X - X = 0`. The assertion cannot fail.

**Why this matters:**
Same pattern as F1: a defined-by-definition relation is being presented as a verification. The Mathematica script claims an extra check that the SymPy script doesn't have, but the extra check is tautological.

**Required change:**
Delete lines 79 and 83 of the Mathematica script. The `G_max - closed form` assertion at line 82 already pins down `gMaxTarget` symbolically; there is no additional content in restating `(Pi^2*A/8)*gMaxTarget == alphaCrit` when `alphaCrit` is defined as exactly that.

**Verification:**
On re-run, the lines `(Pi^2 A / 8) G_max - alpha_crit = 0` and `PASS: (Pi^2 A / 8) G_max - alpha_crit` must be absent from the saved Mathematica output.

### F4 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:126`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:120`

**What's wrong:**
SymPy line 126:
`expect_true("inadmissible sample: R_target < 1 blocks the branch", bool(sp.Rational(9, 10) < 1), "R_target=9/10")`.
The boolean argument is `bool(Rational(9,10) < 1)`, i.e., the literal `True`. No `R_target` is computed from any parameter setting; the script simply asserts `9/10 < 1`. The label "inadmissible sample" implies a parameter choice was made and produced an inadmissible R_target, but no such parameter choice exists in the script — only the bare number `9/10` is compared to `1`.

Mathematica line 120 mirrors the issue: `expectTrue["inadmissible sample: R_target < 1 blocks the branch", 9/10 < 1, "R_target=9/10"]`.

**Why this matters:**
This assertion provides no test of the inadmissibility branch. A future change that broke R_target's actual formula would not flip this assertion. The label is misleading evidence of coverage.

**Required change:**
Either (a) remove the assertion (preferred, since the admissibility branch is already exercised at line 99 / 101 with a real F_sample value), or (b) construct an actual sample where `F(xi_bad, delta_bad) < 1` is computed from the closed form for some choice of `xi_bad, delta_bad` for which the inequality holds, and assert that result.

Apply option (a): delete line 126 in the SymPy script and line 120 in the Mathematica script.

**Verification:**
On re-run, the lines `inadmissible sample: R_target < 1 blocks the branch = R_target=9/10` and the corresponding `PASS:` line must be absent from both saved outputs.

### F5 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl` (entire script)

**What's wrong:**
The Mathematica script is structurally a port of the SymPy script: same variables in the same order, same closed-form targets typed in directly, same numeric host sample, same `M_mix_good = G - 1/10` / `M_mix_bad = G + 1/10` witness strategy. Examples of corresponding sections:

(a) SymPy line 54: `G = sp.simplify(9 * xi * (xi + delta) / (9 * delta + 11 * xi))`.
    Mathematica line 43: `gTarget = FullSimplify[9*xi*(xi + delta)/(9*delta + 11*xi), ...]`.
    The Mathematica closed form is written down identically; it is not derived from any prior expression.

(b) SymPy lines 75-76:
    `dG_target = sp.simplify(9 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) / (9 * delta + 11 * xi) ** 2)`.
    Mathematica lines 66-69:
    `dGTarget = FullSimplify[9*(9*delta^2 + 18*delta*xi + 11*xi^2)/(9*delta + 11*xi)^2, ...]`.
    Same closed form copied directly.

(c) SymPy lines 104-115 (host kappa derivation: A_host=3, beta0_host=5, x_host=A*xi_sample, deltaK_host=A*delta_sample, N_host built from KAPPA0_SQ, KAPPA1_SQ).
    Mathematica lines 102-110: aHost=3, beta0Host=5, xHost=aHost*xiSample, deltaKHost=aHost*deltaSample, nHost built from kappa0Sq, kappa1Sq. Identical variable choreography with renamed identifiers.

(d) SymPy line 134: `G_series = sp.series(G, xi, 0, 3).removeO()`. Mathematica line 128: `gSeries = FullSimplify[Normal[Series[gTarget, {xi, 0, 2}]], ...]`. Both compute the same Taylor expansion and compare to the same hand-written target `xi - 2*xi^2/(9*delta)`.

None of the Mathematica "closed forms" (`gTarget`, `dGTarget`, `gMaxTarget`, `gSeriesTarget`, `alphaCrit`) are derived inside the Mathematica script; they are written down by hand, matching SymPy's hand-written targets. The second-engine policy requires the engines to derive results independently from the physical premises so that an algebraic error or convention mismatch in one engine is exposed by the other.

**Why this matters:**
A transliteration cannot catch errors in the closed forms themselves. If the SymPy script's `dG_target` were wrong (e.g., a coefficient flipped), the Mathematica script would write down the same wrong target and "verify" it. The two engines must arrive at the same closed form via different algebraic routes, not by typing the same target into both.

**Required change:**
Restructure the Mathematica script so each closed-form target is *derived* in Mathematica, not declared. Specifically, in the .wl script:

1. Replace the explicit declaration of `dGTarget` at lines 66-69 with `dGTarget = FullSimplify[D[gTarget, xi]*((9*delta + 11*xi)^2/9), Assumptions -> $Assumptions]` (or another Mathematica-native simplification path), and assert that `9*dGTarget/(9*delta+11*xi)^2 == dG` after the derivative step. This forces Mathematica to perform the differentiation and algebraic factorization on its own.

2. Replace the explicit declaration of `gMaxTarget` at line 78 with the Mathematica limit form already computed at line 76, i.e., set `gMaxTarget = gMax` (where `gMax` is the symbolic Limit), then verify `gMax == 9*(1 + delta)/(9*delta + 11)` by `FullSimplify[gMax - 9*(1+delta)/(9*delta+11), Assumptions -> delta > 0] == 0`. This way the closed form is the side being tested, not the side being declared.

3. Replace `gSeriesTarget = FullSimplify[xi - 2*xi^2/(9*delta), ...]` at line 129 with a derivation: expand `gTarget` symbolically by series first and assert the result equals `xi - 2*xi^2/(9*delta)` via `Coefficient[gSeries, xi, k]` checks for k=0,1,2 against the expected coefficients `{0, 1, -2/(9*delta)}` extracted via `Solve` or by direct substitution into a polynomial-fitting expression. The "target" must not be declared up front; the coefficients must be read out of `gSeries`.

4. For the kappa-based numeric witness (lines 102-114), substitute the Mathematica-native `Resolve` / `Reduce` framework or rely on `FullSimplify` over `Rationals` so that the witness computation is structured differently from SymPy's `sp.simplify(...)` chain. At minimum, do not rename SymPy identifiers character-for-character.

**Verification:**
After Codex's edits, the Mathematica script's output must still pass, but the structural diff must show derived (not declared) targets. The verifier will spot-check that `dGTarget`, `gMaxTarget`, and `gSeriesTarget` no longer appear as up-front closed-form declarations matching SymPy's.

### F6 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:104-120`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:102-114`

**What's wrong:**
The kappa-based derivation of R_target (the "stronger middle-conjunct witness" at SymPy lines 104-120 and Mathematica lines 102-114) only checks the identity at a single numeric sample (`A_host=3, beta0_host=5, xi_sample=1/2, delta_sample=1`). The relation `F(xi,delta) == R_target_sym(xi,delta,A,beta0)` holds symbolically for all positive A, beta0, all delta>0, all 0<=xi<1; restricting the check to one sample means an error in any of the coefficients `11, 9, 18` inside F or N that happens to cancel at (1/2, 1) would not be caught.

For a checkpoint stage with stated "exact" claims, the witness should be symbolic across (xi, delta, A, beta0). The infrastructure (symbol declarations, kappa values) is already in place — only the symbolic instantiation step is missing.

**Why this matters:**
F is the load-line function whose ratio against `1` controls admissibility at the support-feasibility frontier. A symbolic identity `F = R_target_sym` would close the docstring's "Exact dimensionless support-feasibility function" claim. As written, the script verifies only the closed form `9*xi*(xi+delta)/(9*delta+11*xi)` of G — F is taken on faith except at one sample.

**Required change:**
Add a symbolic version of the host computation, using already-declared symbolic `A` (line 46) and a newly declared symbolic `beta0`. Concretely, after line 120 of the SymPy script, add:

```
A_sym = sp.symbols("A_sym", positive=True, real=True)
beta0_sym = sp.symbols("beta0_sym", positive=True, real=True)
x_sym = A_sym * xi
deltaK_sym = A_sym * delta
N_sym = (
    beta0_sym
    * (KAPPA0_SQ * (x_sym + deltaK_sym) + KAPPA1_SQ * x_sym) ** 4
    / (
        KAPPA0_SQ
        * (A_sym - x_sym)
        * (KAPPA0_SQ * (x_sym + deltaK_sym) ** 2 + KAPPA1_SQ * x_sym ** 2) ** 2
    )
)
R_target_sym = sp.simplify(N_sym * A_sym / (beta0_sym * KAPPA0_SQ))
expect_zero("symbolic kappa derivation: F(xi,delta) - R_target_sym", sp.simplify(F - R_target_sym))
```

Add the analogous block in the Mathematica script after line 114, using independent Mathematica simplification (`FullSimplify[Together[Expand[fTarget - rTargetSym]], Assumptions -> ...]`).

**Verification:**
On re-run, both outputs must contain a new PASS line of the form `symbolic kappa derivation: F(xi,delta) - R_target_sym = 0` (SymPy) / `PASS: symbolic kappa derivation: F(xi,delta) - R_target_sym` (Mathematica).

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script. Closed-form targets (`gTarget`, `dGTarget`, `gMaxTarget`, `gSeriesTarget`, `alphaCrit`) are written down by hand in identical form to the SymPy script rather than being derived inside Mathematica. The host sample (A=3, beta0=5, xi=1/2, delta=1) is numerically identical with renamed identifiers. See F5 for quoted excerpts. The Mathematica script's unique line (the `alphaCrit` check) is itself tautological; see F3.

## Engine cross-check

Both engines produce identical printed output:

- `G(xi,delta) = 9*xi*(delta + xi)/(9*delta + 11*xi)` (SymPy) / `(9*xi*(delta + xi))/(9*delta + 11*xi)` (Mathematica) — agree.
- `F(xi,delta) = -(9*delta + 11*xi)**4 / (81*(xi-1)*(...)^2)` (SymPy) / `-1/81*(9*delta + 11*xi)^4/((-1+xi)*(...)^2)` (Mathematica) — algebraically identical (factor of `-1/(xi-1)` vs `-1/(-1+xi)`).
- Sample values `M_mix=53/145, G=27/58, R_target_host=1414562/558009` agree.
- All assertions pass in both engines.

No engine disagreement. However, agreement here means little, given that the Mathematica script is a transliteration (F5) — the two engines were not run independently against the same physical premises.

## Verdict justification

The script holds up algebraically: the differentiation, limit, and series checks (A3, A5, A12; B5, B7, B15) genuinely exercise computer-algebra machinery, and the kappa-derived host sample (A8/B11) is a real numeric witness for `F = R_target` at a non-trivial parameter point. But several of the headline assertions — `R_target - pi^2 A NQ/(8 beta0)`, `g_B,req^2/varpi^2 - (pi^2 A/8)(G - M_mix)`, `(Pi^2 A/8) G_max - alphaCrit` in Mathematica, and the "inadmissible sample" `9/10 < 1` — are guaranteed by how the symbols were defined and provide no test. The Mathematica script is structurally a port of the SymPy script and so cannot catch errors in the closed-form targets. Findings F1-F4 each remove a hollow assertion; F5 forces the second engine to be a genuine cross-check; F6 closes the docstring's "exact" claim with a symbolic kappa derivation. No `UNFIXABLE` or `CRITICAL_DOWNSTREAM` is warranted — the fixes are local edits to this stage's scripts that strengthen rather than invalidate the verified content.

## Self-test notes

Walked through F2/F6's proposed symbolic kappa derivation by hand: with `KAPPA0_SQ = 8/pi^2`, `KAPPA1_SQ = 16/(9*pi^2) = (2/9)*KAPPA0_SQ`, `x = A*xi`, `deltaK = A*delta`, the inner-numerator expression `KAPPA0_SQ*(x+deltaK) + KAPPA1_SQ*x = KAPPA0_SQ*((11/9)*x + deltaK) = (KAPPA0_SQ*A/9)*(11*xi + 9*delta)`, and the inner-denominator expression `KAPPA0_SQ*(x+deltaK)^2 + KAPPA1_SQ*x^2 = (KAPPA0_SQ*A^2/9)*(9*delta^2 + 18*delta*xi + 11*xi^2)`. Substituting into `N*A/(beta0*KAPPA0_SQ)` collapses to `(9*delta+11*xi)^4 / (81*(1-xi)*(9*delta^2+18*delta*xi+11*xi^2)^2) = F`, confirming the assertion `F - R_target_sym == 0` is a non-trivial identity that will simplify to 0 only if the kappa coefficients are correct. Verified the parity/symmetry of the series check (F's denominator is symmetric in delta, even under `delta -> delta`; the series expansion is in `xi`, not on a symmetric domain — no parity trap). Confirmed `A_sym` and `beta0_sym` actually appear in N_sym (not just as wrappers around an identity), so they will cancel in `R_target_sym = N_sym*A_sym/(beta0_sym*KAPPA0_SQ)` only because they are factored out at construction — verified by hand that `A_sym` and `beta0_sym` enter and leave the simplified expression. Path specifications: directive targets the existing `.py` in `scripts/` and `.wl` in `mathematica/`; no new files required.
