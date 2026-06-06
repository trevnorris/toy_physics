---
unit_id: 110
batch: IV.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 110

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl:33-54`. The
direct ratio-series route was removed and replaced by two undetermined-coefficient
`Solve` systems (order-by-order polynomial inversion):

1. **z-jet (coefficients route), wl:33-41.** Introduces a generic degree-5 polynomial
   `yRJet = Sum[yCoeff[k] z^k, {k,0,5}]`, forms the defining identity residual
   `lambdaR*yRJet - (-3 + rho)`, extracts its `z^k` coefficient equations for `k=0..5`
   via `Table[Coefficient[Expand[...], z, k] == 0]`, and solves the resulting LINEAR
   system `First[Solve[jetEquations, jetCoeffs]]`. `yRSeries` is the reconstructed jet.
   This is an order-by-order inversion of `lambdaR * Y_R = (-3 + rho)`, exactly the
   directive's first sanctioned route — NOT a series of the `(-3+rho)/lambdaR` ratio.
2. **rho-jet (linearization route), wl:47-54.** Independently, the linearization is
   obtained by the same undetermined-coefficient method on `(3 - rho)*rhoJet - 3`:
   a degree-2 generic polynomial `rhoJet = Sum[rhoCoeff[k] rho^k]`, the `rho^k`
   coefficient equations, and `First[Solve[rhoEquations, rhoCoeffs]]`. This inverts
   `(3 - rho)*chi_Q^R = 3` order-by-order — NOT `Series[chiR, {rho,0,2}]`.

Coefficient extraction from the reconstructed `yRSeries` (wl:43-46) still uses
`Coefficient[...,z,k]` with the `/I` and `/(1/27)` normalizations, but those operate on
a polynomial the `.wl` itself built from a linear solve, not on a SymPy-mirrored series.

**Assessment:**
Genuine independence — resolved.

- **Series-on-ratio gone:** `grep "Series["` over the current `.wl` returns nothing
  (exit 1). Both `Series[yR,{z,0,5}]` (the shared black-box of the `(-3+rho)/lambdaR`
  ratio) and `Series[chiR,{rho,0,2}]` are eliminated. Acceptance criterion 1 met.
- **Structurally distinct primitive:** the `.py` (`scripts/...sympy_audit.py:11,17`) reaches
  the coefficients by `sp.series(Y_R, z, 0, 6)` on the explicit ratio and the
  linearization by `sp.series(chi_R, rho, 0, 3)`. The new `.wl` reaches BOTH by an
  undetermined-coefficient linear `Solve` against the defining identity — a different
  native primitive (`Solve` of a coefficient-matched linear system vs. a series
  expansion of a closed-form ratio) and a different intermediate structure (generic
  jet + residual equations vs. a directly-expanded fraction). It is not a renamed
  transliteration. Acceptance criterion 3 met.
- **No regression:** `c2`, `c4`, `c5`, `chi_Q^R = -3/(-3+rho)` (= `3/(3-rho)`), and the
  linearization `1 + rho/3 + rho^2/9` are all still computed, all five `expectZero`
  asserts remain (wl:63-67) with the SAME RHS literals, and the Mathematica log shows
  all five `PASS`. Emitted values unchanged (`c2 = (9-3 rho)^-1`, `c4 = -1/9 (-4+rho)/(-3+rho)^2`,
  `c5 = (27-9 rho)^-1`, `chi_Q^R = -3/(-3+rho)`). Acceptance criterion 2 met.
- **Assertions non-tautological:** the LHS coefficients now come from solving the linear
  jet system, the RHS are the paper's boxed literals; the `expectZero`/`FullSimplify`
  gate fails if the solve produces anything else. No hardcoded LHS, no tautology.
- **`rho != 3` assumption preserved:** `$Assumptions = Element[{z, rho}, Reals] && rho != 3;`
  retained (wl:29).
- **Collateral edit (benign):** wl:41 applies a single rewrite rule
  `(-z^2/(3*(-3+rho))) -> z^2/(9-3*rho)` to `yRSeries`. This is a cosmetic re-normalization
  of one algebraically-equal term, affecting only the printed `Y_R(z)` display so it
  byte-matches the prior committed output; it does NOT manufacture `c2`, since `c2` is
  re-extracted and the assertion `c2 - 1/(9-3 rho) == 0` independently re-validates via
  `FullSimplify`. Non-load-bearing for correctness or independence.

## Exec log assessment

**SymPy:** exit=0. `.py` UNCHANGED (reference engine, per directive). Notable lines:
`chi_Q^R = -3/(rho - 3)`, `chi_Q^R linearized = rho**2/9 + rho/3 + 1`, `stage110: PASS`.

**Mathematica:** exit=0. Notable lines:
`chi_Q^R = -3/(-3 + rho)`, `chi_Q^R linearized = 1 + rho/3 + rho^2/9`,
`PASS: c2 - 1/(9 - 3 rho)`, `PASS: c4 - (4 - rho)/(9 (3 - rho)^2)`,
`PASS: c5 - 1/(27 - 9 rho)`, `PASS: chi_Q^R - 3/(3 - rho)`,
`PASS: chi_Q^R linearized - (1 + rho/3 + rho^2/9)`, `Stage 110 Mathematica audit passed.`

**Output freshness:** confirmed. The committed `.txt` (mtime 01:57:20) is newer than the
`.wl` script (mtime 01:52:59) — regenerated post-fix. `git status` shows no uncommitted
change to the output, and `git diff HEAD` on the committed `.txt` returns no diff:
the deliverable output is **byte-identical** to HEAD (method-only change). Acceptance
criterion 4 met. SymPy `.txt` (01:57:20) likewise fresh.

## Material-change assessment

`material_change`: false. This is a method-only re-author of the second engine. No
emitted value changed; the committed Mathematica output is byte-identical to HEAD and the
SymPy reference engine was not touched. Downstream units depending on stage 110's
deliverables (`c2/c4/c5`, `chi_Q^R = 3/(3-rho)`, linearization) see no numeric change.

## Side observations (non-blocking)

- The wl:41 normalization rewrite rule is a slightly brittle idiom (it pattern-matches one
  specific surface form to preserve byte-identical display). It is harmless here because
  the assertions re-validate `c2` independently, but a future refactor that changes the
  `Solve` output's surface form would silently no-op the rule and only alter the printed
  `Y_R(z)` string — never a correctness regression, since the `expectZero` gates are the
  real check. Non-blocking.

## Verdict justification

The user-authorized re-author closes F1. The `.wl` now reaches the Robin coefficients and
`chi_Q^R = 3/(3-rho)` (plus its `rho`-linearization) via undetermined-coefficient linear
`Solve` against the defining identities `lambdaR*Y_R = (-3+rho)` and `(3-rho)*chiR = 3`
— a structurally distinct native primitive from the SymPy `Series`-on-ratio route. Every
`Series[]` call is gone (grep exit 1), all five assertions remain non-tautological and
PASS, the `rho != 3` domain is preserved, and the committed Mathematica output is
byte-identical to HEAD. Independence defect resolved, no correctness or numeric change.
Verdict: verified.
