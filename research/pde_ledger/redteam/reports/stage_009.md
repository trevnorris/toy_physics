---
unit_id: 009
batch: I.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-20T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 009 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.txt`
- mathematica output: `(missing)`

## What the script claims to verify

The script claims to extend the projection-first Maxwell analysis to a near-throat domain, producing five auditable results: (i) the exact projected inhomogeneous Maxwell law on a throat interval with surviving mixed-sector flux derivative `<∂_w(Z F^{wν})>`; (ii) a symmetric narrow-kernel Taylor expansion whose first nontrivial width correction is `O(σ²)` (m₂ moment); (iii) a one-sided mouth kernel `W_ell(w) = exp(-w/ell)/ell` for which the apparent `1/ell` boundary singularity cancels and the projected derivative expands as `q1 + ell·q2 + ell²·q3 + …`; (iv) zero-mode effective parameters `μ_eff = μ₀<S>/<Z>`, `ξ_eff = ξ<Z>/<H>`, with exact cancellation when `H=Z` or `S=CZ`; (v) concrete Gaussian-localizer examples giving series `1 − σ²/λ² + 3σ⁴/(2λ⁴) + O(σ⁶)` (symmetric) and `1 − 2ell²/λ² + 12ell⁴/λ⁴ − 120ell⁶/λ⁶ + O(ell⁸)` (mouth). The bottom line is delivered by SymPy assertions in sections 2, 3, 4, and 5; sections 1 and 6 are prose only.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 86 | `assert_zero("Gaussian even-kernel Q moments", avg_Q_even_gauss - avg_Q_even.subs({sigma:1, m2:1, m4:3}))` | partial |
| A2 | sympy | 87 | `assert_zero("Gaussian even-kernel derivative moments", avg_dQ_even_gauss - avg_dQ_even.subs({sigma:1, m2:1, m4:3}))` | partial |
| A3 | sympy | 115 | `assert_zero("half-line boundary recombination", recombined_half - avg_dQ_half)` | yes |
| A4 | sympy | 116 | `assert_zero("half-line derivative expansion", avg_dQ_half - (q1 + ell*q2 + ell**2*q3 + ell**3*q4))` | yes |
| A5 | sympy | 192 | `assert_zero("H=Z before effective-gauge cancellation", IH_conc.subs(H_sym, Z_conc) - IZ_conc)` | no (tautology) |
| A6 | sympy | 193 | `assert_zero("S=CZ before effective-coupling cancellation", IS_conc.subs(S_sym, C*Z_conc) - C*IZ_conc)` | no (tautology) |
| A7 | sympy | 194 | `assert_zero("H=Z effective gauge", sp.simplify((xi*IZ_conc/IH_conc).subs(H_sym, Z_conc)) - xi)` | no (tautology) |
| A8 | sympy | 195 | `assert_zero("S=CZ effective coupling", sp.simplify((mu0*IS_conc/IZ_conc).subs(S_sym, C*Z_conc)) - C*mu0)` | no (tautology) |
| A9 | sympy | 222 | `assert_zero("mouth Gaussian asymptotic from erfc closed form", IWZ_mouth_gauss_series - (1 - 2*ell**2/lam**2 + 12*ell**4/lam**4 - 120*ell**6/lam**6))` | partial (literal vs literal) |
| A10 | sympy | 223 | `assert_zero("mouth Gaussian asymptotic from Taylor integration", IWZ_mouth_gauss_series - IWZ_mouth_taylor_integral)` | yes |
| — | — | — | symmetric Gaussian series in section 5 (`<Z> = 1 - σ²/λ² + 3σ⁴/(2λ⁴) + …`) | NO ASSERTION |
| — | — | — | symmetric μ_eff^(sym), ξ_eff^(sym) series in section 4 | NO ASSERTION |
| — | — | — | mouth μ_eff^(mouth), ξ_eff^(mouth) series in section 4 | NO ASSERTION |

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `(missing)` — no Mathematica script for unit 009

**What's wrong:**
The manifest entry at `redteam/MANIFEST.yaml:235–262` records `is_status_only_candidate: false`, `is_checkpoint: false`, and `mathematica: path: null, exists: false`. Per the audit prompt's status rules, a non-checkpoint, non-status-only unit requires both SymPy and Mathematica independent verification scripts. Only the SymPy script (`moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`) exists. There is no `.wl`/`.m` companion that re-derives the near-throat projected-Maxwell results (sections 2, 3, 4, 5) from independent algebra.

**Why this matters:**
Without a second engine, every assertion in the SymPy script is single-engine. A SymPy-specific simplification glitch (e.g., a wrong `series` truncation at `r → ∞`, or `integrate` returning a branch that hides a sign) cannot be cross-checked. The most exposed pieces are (a) the closed-form hand-typed `erfc` form on line 214 — which feeds the entire mouth-Gaussian asymptotic in section 5 — and (b) the half-line integration-by-parts recombination at lines 110–115 where SymPy must compute the `1/ell` boundary piece and its W'-integral consistently.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl` (and corresponding output `.../scripts/output/...near_throat_mathematica_audit.txt`) that independently re-derives the same claims through different Mathematica idioms (e.g., `Integrate` with explicit `Assumptions`, `Series[…, {ell, 0, …}]`, `FullSimplify`, `Print[…]` plus `If[FullSimplify[diff] =!= 0, Exit[1]]` guards). It must not transliterate the SymPy code; in particular, derive the erfc closed form by direct Mathematica integration rather than retyping the literal.

**Verification:**
After Codex creates the file, the verifier runs `redteam exec-mathematica 009` and confirms (a) the script exits 0, (b) it contains assertion guards (`Exit[1]`) for at least the claims M1–M5 enumerated in the directive, (c) its printed series for the symmetric and mouth Gaussian cases match the SymPy output literals `1 - sigma^2/lambda^2 + 3*sigma^4/(2*lambda^4)` and `1 - 2 ell^2/lambda^2 + 12 ell^4/lambda^4 - 120 ell^6/lambda^6`.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:188–195`

**What's wrong:**
The four "exact profile-alignment" assertions on lines 192–195 are tautological. Inspect the construction:

```
H_sym = sp.Function("H")(w)
S_sym = sp.Function("S")(w)
IH_conc = sp.Integral(W_conc * H_sym, (w, 0, sp.oo))   # UNEVALUATED Integral
IS_conc = sp.Integral(W_conc * S_sym, (w, 0, sp.oo))   # UNEVALUATED Integral
IZ_conc = sp.integrate(W_conc * Z_conc, (w, 0, sp.oo)) # evaluated
assert_zero("H=Z before effective-gauge cancellation", IH_conc.subs(H_sym, Z_conc) - IZ_conc)
assert_zero("S=CZ before effective-coupling cancellation", IS_conc.subs(S_sym, C*Z_conc) - C*IZ_conc)
assert_zero("H=Z effective gauge", sp.simplify((xi*IZ_conc/IH_conc).subs(H_sym, Z_conc)) - xi)
assert_zero("S=CZ effective coupling", sp.simplify((mu0*IS_conc/IZ_conc).subs(S_sym, C*Z_conc)) - C*mu0)
```

After `IH_conc.subs(H_sym, Z_conc)`, the unevaluated `IH_conc` literally becomes `sp.Integral(W_conc*Z_conc, (w, 0, sp.oo))` — the same integrand `IZ_conc` integrates. Both evaluate to the same closed form under `simplify`, so A5 reduces to `IZ_conc - IZ_conc == 0`. Likewise A6 reduces to `C*IZ_conc - C*IZ_conc == 0`. A7 then computes `xi * IZ_conc / IZ_conc = xi`, and A8 computes `mu0 * C * IZ_conc / IZ_conc = C*mu0`. These checks are guaranteed by symbolic substitution + linearity of the integral; they verify nothing about the physics claim that "when the gauge weight profile coincides with the kinetic weight profile, the effective gauge parameter is exact." That claim has actual content only when `H` and `Z` enter the *zero-mode-reduced* expression through different routes (e.g., from different sectors of the bulk equation), or when the equality is tested under at least one nontrivial perturbation (e.g., `H = Z + epsilon h(w)` with `epsilon → 0`).

**Why this matters:**
The script asserts `STATUS: PASS` partly on the strength of these four checks. They give a false signal of verification: any change to the underlying derivation that breaks the `H=Z ⟹ ξ_eff=ξ` reasoning would still leave these checks passing because they exercise only the symbol-substitution rule, not the derivation path.

**Required change:**
Replace the four tautological assertions on lines 192–195 with one or both of the following substantive checks (keep them adjacent to line 195):

1. Perturbative test: define `epsilon = sp.symbols("epsilon", real=True)`, set `H_pert = Z_conc + epsilon*(h0 + h1*w + h2*w**2/2)`, compute `xi_eff_pert = sp.simplify(xi * IZ_conc / sp.integrate(W_conc * H_pert, (w, 0, sp.oo)))`, expand in `epsilon` to first order, assert that the `epsilon^0` coefficient equals `xi` (recovers the exact-cancellation claim) AND that the `epsilon^1` coefficient is `-xi * ∫W·(h0+h1*w+h2*w²/2) / IZ_conc` (recovers the leading correction). Do the analogous test for `S_pert = C*Z_conc + epsilon*(s0 + s1*w + s2*w**2/2)`.

2. Drop the four current assertions entirely and instead re-derive `ξ_eff = ξ` from the symbolic expansions `xi_eff_half` (line 154) by substituting `h0 = z0, h1 = z1, h2 = z2` and asserting `sp.simplify(xi_eff_half - xi) == 0`. This makes the cancellation a property of the actual expanded ratio, not of the symbol-substitution layer.

Prefer option (1): it adds discriminating power without removing assertions.

**Verification:**
After Codex applies, the new lines must contain `assert_zero(...)` calls whose left-hand sides are *not* of the form `expr.subs(F(w), φ(w)) - expr_with_φ`. The verifier reruns the SymPy script and confirms exit 0; the audit output must include at least one new line label distinct from "H=Z …" / "S=CZ …".

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:147–156` (effective-parameter expansions)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:205–208` (symmetric Gaussian)

**What's wrong:**
Three quantitative claims that the script prints — and that section 4 / section 5 prose calls out as bottom-line results — are never asserted:

- Line 149: `mu_eff_sym_series = sp.series(mu_eff_sym, sigma, 0, 3).removeO()` produces the symmetric μ_eff expansion. The output (line 87 of the .txt) shows `mu0*s0/z0 + sigma**2*(-m2*mu0*s0*z2/(2*z0**2) + m2*mu0*s2/(2*z0))`. No assertion confirms this expansion equals the hand-derivable formula `μ₀ s₀/z₀ · [1 + (m₂σ²/2)(s₂/s₀ − z₂/z₀)] + O(σ⁴)`.
- Line 150: `xi_eff_sym_series` analogously unasserted; the printed form `xi*z0/h0 + sigma**2*(m2*xi*z2/(2*h0) - h2*m2*xi*z0/(2*h0**2))` should be checked against `ξ z₀/h₀ · [1 + (m₂σ²/2)(z₂/z₀ − h₂/h₀)] + O(σ⁴)`.
- Lines 155–156: analogous mouth expansions `mu_eff_half_series`, `xi_eff_half_series` are likewise only printed.
- Line 207–208: `IWZ_sym_gauss_series = sp.series(IWZ_sym_gauss, sigma, 0, 5).removeO()` is the symmetric Gaussian sampling. The text on line 233 claims `<Z> = 1 - σ²/λ² + 3σ⁴/(2λ⁴) + O(σ⁶)`, but no `assert_zero` ties the series to this literal (compare the mouth case on line 222, which *does* assert against a literal — even if hardcoded).

Section 3's analogous claim about `<Q>_ell = q0 + ell*q1 + ell²*q2 + ell³*q3 + ell⁴*q4` (line 110, printed on output line 55) is also unasserted, though A4 covers `<∂_w Q>_ell`.

**Why this matters:**
The "first mismatch measures local shape-slippage" corollaries (output lines 94–96) and the claim that "the first projection-width corrections are even, O(σ²)" (output line 46) ride entirely on these unasserted series. A SymPy `series` boundary change or a typo in `mu_eff_sym`/`xi_eff_sym` definitions would not be caught.

**Required change:**
Add four `assert_zero(...)` lines immediately after the four printed series, anchoring each to its expected canonical form:

1. After line 149, assert:
   `assert_zero("symmetric μ_eff series", mu_eff_sym_series - (mu0*s0/z0 + (m2*sigma**2/2)*(mu0*s2/z0 - mu0*s0*z2/z0**2)))`
2. After line 150, assert:
   `assert_zero("symmetric ξ_eff series", xi_eff_sym_series - (xi*z0/h0 + (m2*sigma**2/2)*(xi*z2/h0 - xi*z0*h2/h0**2)))`
3. After line 155, assert:
   `assert_zero("mouth μ_eff series", mu_eff_half_series - (mu0*s0/z0 + ell*(mu0*s1/z0 - mu0*s0*z1/z0**2)))`
4. After line 156, assert:
   `assert_zero("mouth ξ_eff series", xi_eff_half_series - (xi*z0/h0 + ell*(xi*z1/h0 - xi*z0*h1/h0**2)))`

Also add, after line 208, the symmetric Gaussian series literal check:
`assert_zero("symmetric Gaussian asymptotic literal", IWZ_sym_gauss_series - (1 - sigma**2/lam**2 + 3*sigma**4/(2*lam**4)))`.

Finally, after line 110, add the `<Q>_ell` literal check:
`assert_zero("half-line Q expansion", avg_Q_half - (q0 + ell*q1 + ell**2*q2 + ell**3*q3 + ell**4*q4))`.

**Verification:**
The verifier runs SymPy; the rerun output must contain six new label lines: "symmetric μ_eff series", "symmetric ξ_eff series", "mouth μ_eff series", "mouth ξ_eff series", "symmetric Gaussian asymptotic literal", "half-line Q expansion". Exit code 0.

### F4 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py:210–222`

**What's wrong:**
Line 214 hardcodes a closed-form `erfc` expression rather than deriving it from the integral:

```
IWZ_mouth_gauss = sp.simplify(sp.integrate(Wm * Zg, (w, 0, sp.oo)))   # line 212: the real integral
r = sp.symbols("r", positive=True)
IWZ_mouth_gauss_r = sp.sqrt(sp.pi) * lam * r * sp.erfc(lam * r / 2) * sp.exp(lam**2 * r**2 / 4) / 2   # line 214: hardcoded
IWZ_mouth_gauss_series = sp.simplify(
    sp.series(IWZ_mouth_gauss_r, r, sp.oo, 8).removeO().subs(r, 1 / ell)
)
```

The asymptotic series `IWZ_mouth_gauss_series` is taken from the hand-typed `IWZ_mouth_gauss_r`, NOT from the SymPy-evaluated integral `IWZ_mouth_gauss` computed on line 212. The actual SymPy integration result on line 212 is computed, printed, and then ignored for verification. Assertion A9 on line 222 then compares the hand-typed-derived series to another hardcoded literal `(1 - 2*ell**2/lam**2 + 12*ell**4/lam**4 - 120*ell**6/lam**6)`. Only A10 (line 223) provides an independent cross-check (Taylor-expand `Zg` then integrate against `Wm`); this is what actually validates the result, but the hand-typed `IWZ_mouth_gauss_r` is structurally a pre-baked answer rather than a derivation.

A typo in line 214 (e.g., `lam * r / 2` written as `lam * r`) would propagate into both `IWZ_mouth_gauss_series` and the comparison at A9 in a way that A10 *would* catch — but only because A10 happens to provide an independent path. If A10 were removed or weakened, the script would degenerate into "literal == literal."

**Why this matters:**
The unit's section 5 mouth-Gaussian numerical fingerprint (output line 122) is the only concrete numeric the script produces; relying on a hardcoded closed form to derive it means the SymPy `integrate` call on line 212 is decorative. The reader/auditor cannot tell from the script alone whether SymPy actually agrees that the integral equals `sqrt(pi)*lambda*(1 - erf(lambda/(2*ell)))*exp(lambda**2/(4*ell**2))/(2*ell)`.

**Required change:**
Replace lines 213–217 so the asymptotic series is derived from the SymPy-evaluated integral, not from a hand-typed closed form. Specifically:

Before (lines 213–217):
```
r = sp.symbols("r", positive=True)
IWZ_mouth_gauss_r = sp.sqrt(sp.pi) * lam * r * sp.erfc(lam * r / 2) * sp.exp(lam**2 * r**2 / 4) / 2
IWZ_mouth_gauss_series = sp.simplify(
    sp.series(IWZ_mouth_gauss_r, r, sp.oo, 8).removeO().subs(r, 1 / ell)
)
```

After:
```
r = sp.symbols("r", positive=True)
# substitute ell = 1/r in the SymPy-evaluated integral so we can series-expand at r → ∞
IWZ_mouth_gauss_r = sp.simplify(IWZ_mouth_gauss.subs(ell, 1/r))
IWZ_mouth_gauss_series = sp.simplify(
    sp.series(IWZ_mouth_gauss_r, r, sp.oo, 8).removeO().subs(r, 1 / ell)
)
# guard: the SymPy integral really equals the erfc closed form (otherwise the substitution above is meaningless)
IWZ_mouth_gauss_erfc = sp.sqrt(sp.pi) * lam * sp.erfc(lam / (2*ell)) * sp.exp(lam**2 / (4*ell**2)) / (2*ell)
assert_zero("mouth Gaussian integral = erfc closed form", IWZ_mouth_gauss - IWZ_mouth_gauss_erfc)
```

This makes the closed-form a tested *consequence* of `sp.integrate`, not a pre-baked answer. The existing A9 / A10 assertions on lines 222–223 then remain valid.

**Verification:**
Verifier reruns SymPy; new output must contain the line label "mouth Gaussian integral = erfc closed form". Exit 0. The printed `IWZ_mouth_gauss` and its series must match the existing output literals (lines 121–122 of the .txt) up to algebraic equivalence.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` script exists. See F1.

## Engine cross-check

Not applicable — only the SymPy engine is present. See F1.

## Verdict justification

The script's structurally good assertions are A3, A4, A10 — these are non-trivial cross-checks (integration-by-parts recombination, explicit `<∂_w Q>_ell` expansion, Taylor-vs-closed-form asymptotic match) that would catch real algebra errors. A1, A2 are weak but valid moment-matching checks. The four "exact profile-alignment" assertions A5–A8 are tautological (symbol substitution + integral linearity), the mouth-Gaussian section uses a hand-typed closed form (F4), and several headline series in sections 4 and 5 print without an assertion (F3). The Mathematica engine is absent (F1), violating the non-status-only, non-checkpoint policy. None of the issues are mathematically irreconcilable (no `UNFIXABLE`); none would invalidate quoted constants in later units (the unit's outputs are symbolic expansions and one numerical fingerprint, both of which would survive the proposed substantive replacements), so no `CRITICAL_DOWNSTREAM`. Attacks tried: looking for sign errors in integration-by-parts (correctly handled: the `-q0/ell` boundary piece cancels the `+q0/ell` in `-<W'Q>`); looking for a missing branch in the erfc series at `r → ∞` (the series is real-coefficient and the asymptotic regime is `r → ∞` ⇔ `ell → 0`, which is the stated near-mouth limit; no branch ambiguity); looking for assumption errors in `sigma, ell, lam, xi, mu0 > 0` (consistent with the physical setup of widths and couplings; `m2, m4` declared `real` with no positivity, which is correct since they are kernel moments and `m4` could be small relative to `m2²` without breaking anything).
