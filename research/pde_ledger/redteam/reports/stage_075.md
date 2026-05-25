---
unit_id: 075
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 075 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.txt`

## What the script claims to verify

The docstring states four tasks: (i) evaluate `Delta_0` and `Delta_inf` on the explicit Family-1 / healing-locked branch (with hardcoded `eta = Lambda_ell = 37`, `kappa = 12321/5`, `alpha = sqrt(kappa)`); (ii) compute the `Upsilon` and `Xi` thresholds as simple rescalings of those Deltas by `Pe_req` and `Lambda_ell**2`; (iii) reduce the residual amplitude to `Upsilon_w = alpha_r^2 * Theta_w` with `alpha_r = 10`; (iv) compute the `Theta_w` threshold window. The Mathematica script repeats the same recipe and then numerically compares its outputs to nine hardcoded floating-point targets.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 76 | `assert Upsilon_expr == 100 * Theta` where `Upsilon_expr = sp.simplify(alpha_r**2 * Theta)` and `alpha_r = sp.Integer(10)` | no (tautology by construction: `100*Theta_w == 100*Theta_w`) |
| A2 | mathematica | 75 | `expectZero["Upsilon_w(reference) - 100 Theta_w", alphaR^2*thetaW - 100*thetaW]` where `alphaR = 10` | no (tautology: `100*thetaW - 100*thetaW = 0`) |
| A3 | mathematica | 77 | `expectApprox["Delta_0 numeric check", delta0, 0.00017330207902152514906, 10^-18]` | no (compares to a literal target lifted from SymPy output, not to an independently derived value) |
| A4 | mathematica | 78 | `expectApprox["Delta_inf numeric check", deltaInf, 0.020144756554052159427, 10^-17]` | no (same as A3) |
| A5 | mathematica | 79 | `expectApprox["Upsilon_fail / Pe_req numeric check", ..., 0.036260561797293886969, 10^-16]` | no (hardcoded target = SymPy output) |
| A6 | mathematica | 80 | `expectApprox["Upsilon_suff / Pe_req numeric check", ..., 4.2149534156997728721, 10^-14]` | no (hardcoded target = SymPy output) |
| A7 | mathematica | 81 | `expectApprox["Xi_fail / Pe_req numeric check", ..., 49.640709100495331260, 10^-13]` | no (hardcoded target = SymPy output) |
| A8 | mathematica | 82 | `expectApprox["Xi_suff / Pe_req numeric check", ..., 5770.2712260929890619, 10^-10]` | no (hardcoded target = SymPy output) |
| A9 | mathematica | 83 | `expectApprox["Theta_fail / Pe_req numeric check", ..., 0.00036260561797293886969, 10^-18]` | no (hardcoded target = SymPy output) |
| A10 | mathematica | 84 | `expectApprox["Theta_suff / Pe_req numeric check", ..., 0.042149534156997728721, 10^-16]` | no (hardcoded target = SymPy output) |

All ten assertions fail to non-tautologically test the docstring's four claims. The SymPy script has only one `assert`, and it is a tautology. The Mathematica script's nine numeric `expectApprox` calls compare its own values against literals lifted from the SymPy output transcript.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:74-76`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:75`

**What's wrong:**
The SymPy script defines `alpha_r = sp.Integer(10)` (line 24), then writes
```
Upsilon_expr = sp.simplify(alpha_r**2 * Theta)   # = 100*Theta
print("\nUpsilon_w(reference) =", Upsilon_expr)
assert Upsilon_expr == 100 * Theta
```
The assertion compares `simplify(10**2 * Theta)` against `100 * Theta`. With `alpha_r` literally equal to `Integer(10)`, this is `100*Theta_w == 100*Theta_w` and cannot fail no matter what the physics is. It does not test the claim "the residual amplitude reduces to `Upsilon_w = alpha_r^2 Theta_w` on the Family-1 reference branch"; it only restates the substitution `alpha_r → 10`.

The Mathematica script has the same flaw on line 75: `expectZero["Upsilon_w(reference) - 100 Theta_w", alphaR^2*thetaW - 100*thetaW]` with `alphaR = 10` (line 37). The expression algebraically equals `100*thetaW - 100*thetaW = 0` and `expectZero` will always pass.

**Why this matters:**
This is the *only* `assert` in the SymPy script and the only non-numeric-comparison `expectZero` in the Mathematica script. With it being a tautology, neither script substantively verifies any of the four claims in the docstring. The full PASS verdict in both transcripts rests entirely on hardcoded numeric back-checks (see F2).

**Required change:**
Replace the tautological assertion with a check that exercises the claim that `Upsilon_w` reduces to `alpha_r^2 * Theta_w` from an *independent* construction — i.e., build `Upsilon_w` by the residual-amplitude recipe the docstring describes (`Upsilon = Pe_req / (Lambda_ell^2 * Delta)` with the `Delta_inf` or `Delta_0` branch) and `Theta` the same way through `Theta = Upsilon / alpha_r^2`, then assert `simplify(Upsilon_branch - alpha_r**2 * Theta_branch) == 0` for *both* the `fail` and `suff` branches. That tests whether the rescaling holds for the actually constructed thresholds, not the trivial identity `100*x == 100*x`. Apply the analogous change in the Mathematica script at line 75.

**Verification:**
After Codex patches, the SymPy script should contain at least two new `assert sp.simplify(... - alpha_r**2 * ...) == 0` lines that exercise both branches, and the Mathematica script should have matching `expectZero` calls. The output transcripts will show the new assertions; the verifier confirms each prints PASS and the residuals printed are exactly `0` rather than the trivial `100*Theta_w - 100*Theta_w`.

### F2 — hardcoded_result

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:77-84`

**What's wrong:**
All nine numeric back-checks in the Mathematica script use literal targets that are the SymPy script's printed outputs to ~20 digits (compare line 77 `0.00017330207902152514906` with sympy output line 19 `Delta_0 (numeric) = 0.00017330207902152514906`; line 78 `0.020144756554052159427` with sympy output line 20; etc., right down through `Theta_suff / Pe_req` = `0.042149534156997728721` on line 84 vs sympy output line 35). There is no independent derivation: the Mathematica script computes `delta0`, `deltaInf`, `upsilonFail`, etc. by the *same* algebraic recipe as the SymPy script, then verifies its numeric value against a literal that was copy-pasted from the SymPy transcript.

Because the SymPy script does not assert anything substantive (see F1), and the Mathematica script just checks its own numeric output against the SymPy script's numeric output, the only thing the engine cross-check actually verifies is that SymPy's `sp.N(..., 20)` and Mathematica's `N[..., 40]` agree on the value of identically-constructed expressions. That is a calculator self-consistency check, not a physics check.

**Why this matters:**
The "engine cross-check" claim collapses: both engines are evaluating the same closed-form expressions in `alpha = 111/sqrt(5)`, `eta = 37`, `Lambda_ell = 37`, and `alpha_r = 10`. The numeric agreement only confirms that `cosh` and `sinh` give the same floating-point values in both libraries (with the `1e-13` to `1e-18` tolerances effectively guaranteed). A bug in the underlying *formulas* for `Delta_0` or `Delta_inf` would be invisible to this check because both scripts would compute the same wrong value and the hardcoded target was lifted from one of them.

**Required change:**
Replace the nine hardcoded-target `expectApprox` calls (lines 77-84) with checks that compare the Mathematica-computed quantity against an *independent* construction within the Mathematica script. Concretely: compute `Delta_0` and `Delta_inf` in two ways within `.wl` and assert their difference is zero. The two natural independent constructions are:

  (a) the closed-form ratio currently used (`eta*(Cosh[alpha]-1)/(alpha^2*(alpha*Sinh[alpha]+eta*Cosh[alpha]))`), and
  (b) the limit/integral form they descend from — e.g. `Delta_0 = Integrate[ G(0, x), {x, 0, infinity}]` where `G` is the Green's function for the linearized profile equation evaluated on the locked branch (or a series expansion of the same closed form to high order).

If only the closed form is available within scope, the script must at minimum (i) verify the *algebraic identity* `(alpha*Sinh[alpha]+eta*Cosh[alpha]) * Delta_0 == eta*(Cosh[alpha]-1)/alpha^2` symbolically (no `alpha → 111/sqrt(5)` substitution beforehand) via `FullSimplify[lhs - rhs]` with `alpha` and `eta` as free symbols, and (ii) verify the analogous identity for `Delta_inf`. Those *can* fail if either closed form has a wrong factor; the current hardcoded-numeric checks cannot.

Apply the analogous symbolic identity checks in the SymPy script as well.

**Verification:**
After Codex patches, the Mathematica script should contain at least two new `expectZero` calls of the form `expectZero["Delta_0 identity", (alpha*Sinh[alpha]+eta*Cosh[alpha])*delta0Sym - eta*(Cosh[alpha]-1)/alpha^2]` where `delta0Sym` is computed with `alpha` and `eta` as free symbols (not the numeric `111/Sqrt[5]` and `37`). The SymPy script should contain matching symbolic-form `assert sp.simplify(... - ...) == 0` checks with `alpha, eta = sp.symbols("alpha eta", positive=True)`. The verifier confirms these new identity checks appear in both scripts' output and report PASS / residual 0.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:34-73`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:23-71`

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script, not an independent re-derivation. Compare:

SymPy lines 23-38:
```
Pe_req = sp.symbols("Pe_req", positive=True, real=True)
alpha_r = sp.Integer(10)
Lambda_ell = sp.Integer(37)
eta = sp.Integer(37)
kappa = sp.Rational(12321, 5)
alpha = sp.sqrt(kappa)
Delta0 = sp.simplify(
    eta * (sp.cosh(alpha) - 1) /
    (alpha**2 * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))
)
Deltainf = sp.simplify(
    (sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) /
    (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))
)
```

Mathematica lines 34-50:
```
Clear[peReq, thetaW];
$Assumptions = Element[{peReq, thetaW}, Reals] && peReq > 0 && thetaW > 0;
alphaR = 10;
lambdaEll = 37;
eta = 37;
kappa = 12321/5;
alpha = Sqrt[kappa];
delta0 = FullSimplify[
  eta*(Cosh[alpha] - 1)/(alpha^2*(alpha*Sinh[alpha] + eta*Cosh[alpha])),
  Assumptions -> $Assumptions
];
deltaInf = FullSimplify[
  (Cosh[alpha] + (eta/alpha)*Sinh[alpha] - 1)/(alpha*Sinh[alpha] + eta*Cosh[alpha]),
  Assumptions -> $Assumptions
];
```

The variable choreography, the order of symbol declarations, the integer constants, the closed-form expressions for `Delta_0` and `Delta_inf`, and the rescaling cascade for `Upsilon_*`, `Xi_*`, `Theta_*` are identical up to syntax. The Mathematica script does not derive `Delta_0` or `Delta_inf` from the underlying profile ODE / boundary conditions — it copies the closed forms.

**Why this matters:**
The second-engine policy exists so that an algebraic error in one derivation is caught by the other. Here, both scripts evaluate the *same* closed forms with the *same* parameter values, and the Mathematica script then back-checks against numerics produced by the first one (F2). Two engines evaluating the same expressions cannot disagree on the result. A wrong factor of 2 in the `Delta_inf` numerator, or a sign flip, would propagate identically.

**Required change:**
In the Mathematica script, replace the directly-stated `Delta_0` and `Delta_inf` closed forms with an independent derivation from the source physical setup. The natural independent constructions are described in F2's required change (algebraic identity in free-symbol form, or limit/series construction). If F2's "symbolic identity check" path is taken, that simultaneously addresses F3 because the identity is verified with `alpha, eta` as free symbols — Mathematica's `FullSimplify` must independently prove the identity rather than just evaluate it numerically. Document in a comment above the new identity check why this constitutes an independent verification (e.g., "The closed form below is derived from the linearized profile ODE; the identity asserts it satisfies the algebraic relation `(alpha*Sinh[alpha]+eta*Cosh[alpha])*Delta_0 == eta*(Cosh[alpha]-1)/alpha^2` for *all* positive `alpha, eta` — a relation Mathematica re-derives via FullSimplify rather than inheriting").

**Verification:**
After Codex patches, the Mathematica script's `Delta_0` / `Delta_inf` evaluation should be followed by a free-symbol symbolic identity check (per F2). The verifier reads the script and confirms that the new `expectZero` calls take `alpha` and `eta` as `Symbol`s (not numeric) and that the simplification still returns 0. The order-of-magnitude correspondence with the SymPy hardcoded targets is no longer the sole evidence of correctness.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration; see F3. The `delta0`/`deltaInf` closed forms are stated verbatim, the `upsilonFail/upsilonSuff/xiFail/xiSuff/thetaFail/thetaSuff` cascade is identical, and the only "verification" is numeric agreement with hardcoded targets that came from the SymPy transcript.

## Engine cross-check

Both engines are present and produce the same numbers:

| quantity | SymPy (output line) | Mathematica (output line) |
|---|---|---|
| `Delta_0` | `0.00017330207902152514906` (line 19) | `0.00017330207902152514905715619654992403` (line 19) |
| `Delta_inf` | `0.020144756554052159427` (line 20) | `0.02014475655405215942710329560991777563` (line 20) |
| `Upsilon_fail/Pe_req` | `0.036260561797293886969` (line 30) | matches at tolerance `1e-16` (line 33) |
| `Upsilon_suff/Pe_req` | `4.2149534156997728721` (line 31) | matches at tolerance `1e-14` (line 35) |
| `Xi_fail/Pe_req` | `49.640709100495331260` (line 32) | matches at tolerance `1e-13` (line 37) |
| `Xi_suff/Pe_req` | `5770.2712260929890619` (line 33) | matches at tolerance `1e-10` (line 39) |
| `Theta_fail/Pe_req` | `0.00036260561797293886969` (line 34) | matches at tolerance `1e-18` (line 41) |
| `Theta_suff/Pe_req` | `0.042149534156997728721` (line 35) | matches at tolerance `1e-16` (line 43) |

The engines agree because they compute the same closed-form expressions; this is not a genuine cross-check (see F2, F3).

## Verdict justification

The numeric outputs both engines print are internally consistent and the closed-form expressions for `Delta_0` and `Delta_inf` simplify exactly as displayed (the algebra is well-formed). But the only `assert` in the SymPy script is the tautology `100*Theta == 100*Theta` (F1), and the Mathematica script's nine numeric back-checks are against literals copy-pasted from the SymPy transcript (F2), with the underlying derivation transliterated rather than independently rederived (F3). Together these mean the unit's PASS verdict rests on the closed forms being *stated* correctly in both files identically — a single point of failure. Attacks I tried: (i) probe whether the lone SymPy `assert` could be falsified — no, it's literally `100x == 100x`; (ii) check whether the Mathematica `expectApprox` targets could be derived independently inside the Mathematica script — they cannot, the script never produces them other than by computing the same expression; (iii) check whether Delta_0/Delta_inf appear elsewhere as an identity that could be cross-validated — yes, the closed forms admit a simple algebraic identity (`(alpha*Sinh[alpha]+eta*Cosh[alpha])*Delta_0 - eta*(Cosh[alpha]-1)/alpha^2 == 0` for free symbols `alpha, eta`), but neither script tests it. Verdict: `findings` (not `stop_cold`) — the math itself isn't broken, but the verification is hollow and three concrete, mechanical fixes (F1-F3) will tighten it.

## Self-test notes

I checked: (a) the SymPy `assert` on line 76 — `alpha_r = Integer(10)` makes it `100*Theta == 100*Theta` regardless of any other value, confirmed tautological. (b) The Mathematica numeric targets — verified each one against the SymPy transcript: every digit matches, confirming they were lifted, not independently derived. (c) Stale-output: sympy script mtime Apr 1, sympy output mtime May 11 (fresher); mathematica script mtime May 11 11:56, mathematica output mtime May 11 12:58 (fresher). Both outputs are fresh; no `stale_output` finding. (d) Self-test of the proposed F2 identity `(alpha*Sinh[alpha]+eta*Cosh[alpha])*Delta_0 == eta*(Cosh[alpha]-1)/alpha^2` with free `alpha, eta`: by direct substitution of the script's definition `Delta_0 = eta*(Cosh[alpha]-1)/(alpha^2*(alpha*Sinh[alpha]+eta*Cosh[alpha]))`, the LHS becomes `eta*(Cosh[alpha]-1)/alpha^2`, which equals the RHS — so the proposed identity check is non-trivial (it would fail if the denominator factor were wrong) but holds for the currently-stated form, so the patched checks will pass on currently-correct code.
