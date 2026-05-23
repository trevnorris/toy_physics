---
unit_id: 072
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 072 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.txt`

## What the script claims to verify

The script declares (by hand) two boundary-layer functions `Delta_0 = eta (cosh(alpha)-1) / (alpha^2 (alpha sinh(alpha) + eta cosh(alpha)))` and `Delta_inf = (cosh(alpha) + (eta/alpha) sinh(alpha) - 1) / (alpha sinh(alpha) + eta cosh(alpha))` with `eta = Lambda_ell` and `alpha = sqrt(4 chi_s^2 + (4/5) Lambda_ell^2)`. It then defines two threshold surfaces `Upsilon_fail = Pe_req/(Lambda_ell^2 Delta_inf)` and `Upsilon_suff = Pe_req/(Lambda_ell^2 Delta_0)` and prints corresponding `V_0^2` forms. Four assertions verify two "asymptotic" identities in each of two regimes: shell-gradient-dominated (`kappa ~ (4/5) Lambda_ell^2`, so `alpha ~ (2/sqrt(5)) Lambda_ell`) and compression-dominated (`kappa ~ 4 chi_s^2`, so `alpha ~ 2 chi_s`). Each assertion compares the value of `Pe_req/(Lambda_ell^2 * Delta_*_limit)` to a separately-hand-written target constant.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 65 | `expect_zero(..., Pe_req/(Lambda_ell^2 * DeltaInf_shell) - Upsilon_fail_shell)` | partial |
| A2 | sympy | 66 | `expect_zero(..., Pe_req/(Lambda_ell^2 * Delta0_shell) - Upsilon_suff_shell)` | partial |
| A3 | sympy | 78 | `expect_zero(..., Upsilon_fail_comp - 2 Pe_req chi_s / Lambda_ell^2)` | partial |
| A4 | sympy | 79-80 | `expect_zero(..., Upsilon_suff_comp - 4 Pe_req chi_s^2 (Lambda_ell + 2 chi_s)/Lambda_ell^3)` | partial |
| A5 | mathematica | 71 | `expectZero[..., peReq/(lambdaEll^2 deltaInfShell) - upsilonFailShell]` | partial |
| A6 | mathematica | 72 | `expectZero[..., peReq/(lambdaEll^2 delta0Shell) - upsilonSuffShell]` | partial |
| A7 | mathematica | 83 | `expectZero[..., upsilonFailComp - 2 peReq chiS/lambdaEll^2]` | partial |
| A8 | mathematica | 84-87 | `expectZero[..., upsilonSuffComp - 4 peReq chiS^2 (lambdaEll + 2 chiS)/lambdaEll^3]` | partial |

All eight rows are "partial": each assertion verifies an arithmetic identity between two hand-built algebraic expressions (a "limit form" of `Delta_*` and a separately-typed "shell/comp target"), but neither side of the comparison is derived from the full `Delta_0`/`Delta_inf` expressions that the script's docstring/comments claim to be taking the asymptotic of. The `Delta_*_shell` and `Delta_*_comp` blocks are large-alpha asymptotic forms typed by hand, never obtained via `sp.series`/`sp.limit` or Mathematica's `Series`/`Asymptotic`. The check therefore never exercises the connection between `Delta_0(kappa, eta)` (lines 33-40) and its asymptotic in either physical regime.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py:57-80`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl:61-87`

**What's wrong:**
Both engines claim to verify the shell-gradient-dominated and compression-dominated asymptotics of `Delta_0`/`Delta_inf`, but the four `expect_zero` calls never reference the full closed forms `Delta_0`/`Delta_inf` defined on lines 33-40 (sympy) / 37-44 (wl). Instead, both engines construct hand-written large-alpha forms

```python
# sympy lines 63-64
Delta0_shell = sp.simplify(Lambda_ell / ((c * Lambda_ell)**2 * (c * Lambda_ell + Lambda_ell)))
DeltaInf_shell = sp.simplify((1 + Lambda_ell / (c * Lambda_ell)) / (c * Lambda_ell + Lambda_ell))
```

and

```python
# sympy lines 69-70
Delta0_comp = sp.simplify(Lambda_ell / ((2 * chi_s)**2 * (2 * chi_s + Lambda_ell)))
DeltaInf_comp = sp.simplify((1 + Lambda_ell / (2 * chi_s)) / (2 * chi_s + Lambda_ell))
```

These are leading-order asymptotic forms (obtained by replacing `cosh(alpha) - 1 -> cosh(alpha)` and `sinh(alpha) -> cosh(alpha)` at large `alpha`, then dividing) — but the script never tests that they are in fact the asymptotics of the symbolic `Delta_0`, `Delta_inf` defined earlier. Each asserted residual is a closed-form algebraic identity in the substituted `c, chi_s, Lambda_ell`:

- A1/A5 verifies the arithmetic `(2 + 1/c)/((c+1) Lambda_ell) * 1/Lambda_ell = 2/(sqrt(5) Lambda_ell)` with `c = 2/sqrt(5)`.
- A2/A6 verifies `c^2 (c+1) = (4/5)(1 + 2/sqrt(5))`.
- A3/A7 verifies that `(1 + L/(2 chi_s))/(2 chi_s + L) = 1/(2 chi_s)` (a trivial algebraic identity, independent of any physics).
- A4/A8 verifies that `Pe_req chi_s^2 (2 chi_s + L) 4 / L^3` matches its own factored form.

In all four cases the assertion is reducible to a numeric/symbolic arithmetic identity that does not depend on the values of `Delta_0`, `Delta_inf` ever being correct. The Mathematica engine reproduces the same pattern. The docstring/comments ("Explicit branch placement map and threshold surfaces" and "Shell-gradient dominated asymptotics: kappa ~ (4/5) Lambda_ell^2"; "Compression dominated asymptotics: kappa ~ 4 chi_s^2") claim the asymptotics of `Delta_0`/`Delta_inf` are being verified; they are not.

A particularly suspicious sub-case: in A3/A7, `DeltaInf_comp = (1 + Lambda_ell/(2 chi_s))/(2 chi_s + Lambda_ell)` simplifies algebraically to `1/(2 chi_s)` for *any* `Lambda_ell, chi_s`, with no asymptotic regime invoked at all. So `Pe_req/(Lambda_ell^2 * DeltaInf_comp) - 2 Pe_req chi_s/Lambda_ell^2` is identically zero, regardless of whether `Delta_inf` actually reduces to that form in the `chi_s -> infty` limit.

**Why this matters:**
The unit's headline physics — that the fail/suff threshold surfaces have specific limiting forms in shell-dominated and compression-dominated regimes — is unverified. A sign error, missing factor, or wrong functional form in `Delta_0`/`Delta_inf` would not be caught by any of the four assertions, because none of them references those full expressions. The agreement seen in the saved outputs (`shell fail asymptotic = 0`, etc.) is the kind of trivial-pass case the audit policy explicitly warns against (output transcripts that just say PASS while the underlying check is a tautological algebraic identity).

**Required change:**
Replace each of the four "asymptotic" residual checks in both engines with one that actually takes the asymptotic of the full `Delta_0`/`Delta_inf` symbolically. Concretely, for each engine, introduce a series/limit step that extracts the leading term of `Delta_0(kappa, eta)` and `Delta_inf(kappa, eta)` in the appropriate regime and compare *that* to the hand-written `Delta_*_shell`/`Delta_*_comp`. See directive F1 for the precise edits and parameterization.

**Verification:**
After Codex applies, the verifier will see (a) a new printed line in each engine reading e.g. `Delta_0 shell leading order = ...` and `Delta_inf shell leading order = ...`; (b) the four `expect_zero` / `expectZero` calls now reference the leading-order extracted from the full `Delta_0`/`Delta_inf`, not from the hand-built `Delta_*_shell`/`Delta_*_comp`. Both scripts must still exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl:28-87`

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script's algebra: identical variable choreography (`chiS` / `chi_s`, `lambdaEll` / `Lambda_ell`, `kappa`, `eta`, `alpha`, `delta0`, `deltaInf`, `upsilonFail`, `upsilonSuff`, `v0FailSq`, `v0SuffSq`, `cShell`, `delta0Shell`, `deltaInfShell`, `delta0Comp`, `deltaInfComp`, `upsilonFailComp`, `upsilonSuffComp`), identical right-hand sides, identical hand-typed asymptotic targets, identical banner/print layout, identical assertion order. There is no independent derivation pathway in the `.wl`.

Corresponding sections:

SymPy (lines 29-43):
```
kappa = sp.simplify(4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2)
eta = Lambda_ell
alpha = sp.sqrt(kappa)
Delta0 = sp.simplify(
    eta * (sp.cosh(alpha) - 1) /
    (alpha**2 * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))
)
DeltaInf = sp.simplify(
    (sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) /
    (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))
)
Upsilon_fail = sp.simplify(Pe_req / (Lambda_ell**2 * DeltaInf))
Upsilon_suff = sp.simplify(Pe_req / (Lambda_ell**2 * Delta0))
```

Mathematica (lines 33-47):
```
kappa = FullSimplify[4*chiS^2 + (4/5)*lambdaEll^2, ...];
eta = lambdaEll;
alpha = FullSimplify[Sqrt[kappa], ...];
delta0 = FullSimplify[
  eta*(Cosh[alpha] - 1)/(alpha^2*(alpha*Sinh[alpha] + eta*Cosh[alpha])), ...];
deltaInf = FullSimplify[
  (Cosh[alpha] + (eta/alpha)*Sinh[alpha] - 1)/(alpha*Sinh[alpha] + eta*Cosh[alpha]), ...];
upsilonFail = FullSimplify[peReq/(lambdaEll^2*deltaInf), ...];
upsilonSuff = FullSimplify[peReq/(lambdaEll^2*delta0), ...];
```

And, more critically, the asymptotic blocks at sympy lines 58-66 and wl lines 63-72 are identical letter-for-letter — same `c = 2/sqrt(5)`, same `Delta0_shell`, `DeltaInf_shell` definitions, same target constants. Similarly for the compression block (sympy 69-80 vs wl 76-87).

**Why this matters:**
The two-engine policy exists so that one CAS's simplification path (e.g., Mathematica's automatic `Sinh[]` normal-form choice vs SymPy's `Together`) cannot mask a hidden branch or sign error. When the second script merely re-evaluates the first script's algebra in different syntax, any mistake — including the leading-order replacement ambiguity flagged in F1 — passes through both engines identically. The reported agreement is not independent.

**Required change:**
After F1 lands, the `.wl` script's asymptotic blocks will already differ from the SymPy script's structurally (because each engine extracts its asymptotic via its own native machinery — `sp.series`/`sp.limit` vs `Series`/`Limit`). The remaining transliteration concern is the upfront block defining `Delta_0`/`Delta_inf` and `Upsilon_fail`/`Upsilon_suff`. Refactor the `.wl` so that it derives `Delta_0` and `Delta_inf` from a different starting point than the SymPy script: solve the ODE that has these as Green's-function-style support functions on `[0,1]` with Robin BCs of the form `u'(0) = eta u(0), u'(1) = 0`, then read off `Delta_0 = u(0) / (forcing)` and `Delta_inf = u(1)`. See directive F2 for the precise edits.

**Verification:**
After Codex applies, the verifier will see in the `.wl` output a printed line corresponding to a DSolveValue boundary-value-problem solution, with a `expectZero[..., delta0Derived - delta0]` and `expectZero[..., deltaInfDerived - deltaInf]` check. The `.wl` script must still exit 0, and its variable names for the BVP intermediates must differ from the SymPy script (which won't have a BVP block).

## Independent-derivation check (Mathematica)

Not independent. As documented in F2, the `.wl` script reproduces the SymPy script's variable choreography step-for-step, including the closed forms of `Delta_0`, `Delta_inf`, the threshold surfaces, and the hand-typed asymptotic comparison targets. The Mathematica output's letter-by-letter agreement with the SymPy output is guaranteed by construction.

## Engine cross-check

Both outputs agree at the level claimed:

| Quantity | SymPy (.txt) | Mathematica (.txt) |
|---|---|---|
| `kappa` | `4*Lambda_ell**2/5 + 4*chi_s**2` | `(4*(5*chiS^2 + lambdaEll^2))/5` |
| `Delta_0` | hyperbolic form involving `sqrt(5)*sqrt(Lambda_ell**2+5*chi_s**2)/5` and `cosh`/`sinh` | algebraically equivalent hyperbolic form (presented via `Sinh^2` and `Coth`) |
| `Upsilon_fail_shell` | `2*sqrt(5)*Pe_req/(5*Lambda_ell)` | `(2*peReq)/(Sqrt[5]*lambdaEll)` |
| `Upsilon_suff_shell` | `4*Pe_req*(2*sqrt(5) + 5)/25` | `(4*(5 + 2*Sqrt[5])*peReq)/25` |
| `Upsilon_fail_comp` | `2*Pe_req*chi_s/Lambda_ell**2` | `(2*chiS*peReq)/lambdaEll^2` |
| `Upsilon_suff_comp` | `4*Pe_req*chi_s**2*(Lambda_ell + 2*chi_s)/Lambda_ell**3` | `(4*chiS^2*(2*chiS + lambdaEll)*peReq)/lambdaEll^3` |

Agreement is exact (modulo CAS notation). However, given the transliteration in F2, this is not independent evidence of correctness — both engines arrive at identical answers because they execute identical algebra in different syntax.

## Verdict justification

The script's stated physical claims — shell-gradient-dominated and compression-dominated asymptotics of the threshold surfaces — are not actually exercised by any of the four `expect_zero` calls in either engine. The assertions reduce to algebraic identities between hand-built target expressions and hand-built "limit forms" of `Delta_*`, with no symbolic asymptotic step ever touching the full `Delta_0`/`Delta_inf` expressions. One sub-case (A3/A7, compression `Delta_inf`) is even an unconditional algebraic identity, holding regardless of any asymptotic regime. Independent of that, the `.wl` script is a transliteration of the `.py` script and provides no second-engine check. Both findings are medium severity; the underlying `Delta_0`/`Delta_inf` formulas appear self-consistent (the hyperbolic forms reduce to the leading-order limit forms in the expected regimes when one does take the limit by hand), so neither finding is `UNFIXABLE` or downstream-critical. Outputs are fresh (sympy: script Apr 1, output May 11; wl: script May 11 11:56, output May 11 12:57). Verdict: `findings`, count 2.

## Self-test notes

I checked: (1) **arithmetic of the hand-built limit forms** — manually verified that with `c = 2/sqrt(5)`, `DeltaInf_shell = (1+1/c)/((c+1) L) = sqrt(5)/(2 L)` and `Pe_req/(L^2 * DeltaInf_shell) = 2 Pe_req/(sqrt(5) L)`, matching `Upsilon_fail_shell`; similarly `c^2(c+1) = (4/5)(1+2/sqrt(5))` gives the suff target; the compression `DeltaInf_comp = 1/(2 chi_s)` identically. So the existing assertions are real algebraic identities — they just don't reference `Delta_0`/`Delta_inf`. (2) **leading-order asymptotic of the full `Delta_0`** — for alpha large, `cosh(alpha)-1 ~ cosh(alpha)`, `sinh(alpha) ~ cosh(alpha)`, so `Delta_0 ~ eta cosh(alpha)/(alpha^2 (alpha+eta) cosh(alpha)) = eta/(alpha^2(alpha+eta))`, which matches the hand-built `Delta0_shell` and `Delta0_comp` when `alpha = c*Lambda_ell` or `alpha = 2 chi_s` and `eta = Lambda_ell`. So the directive's proposed `sp.series`/`Series` extraction should reproduce the existing target — meaning the fix won't change pass/fail of the four checks, but it will make them meaningful. (3) **variable-independence trap** — directive F1 takes `series(Delta_0, chi_s, 0, 1)` etc.; verified `Delta_0` depends on `chi_s` via `alpha = sqrt(4 chi_s^2 + 4 Lambda_ell^2 / 5)` so the series in `chi_s` is non-trivial. (4) **alternative-regime parameter trap** — "shell-dominated" means `chi_s -> 0` at fixed large `Lambda_ell`, but the script labels the regime by `kappa ~ (4/5) Lambda_ell^2`, which is the `chi_s -> 0` limit only when `Lambda_ell` is large enough that `alpha` itself is large (otherwise the small-alpha series of `Delta_0` is `~ eta/(2(alpha^2 + eta) + ...) != eta/(alpha^2(alpha+eta))`). The directive's series expansion must be in *both* `chi_s -> 0` AND `Lambda_ell -> infinity`, or equivalently in `alpha -> infinity`. The cleanest formulation is `sp.series` in `1/alpha` around `0` after substituting `alpha` as a free symbol. The directive uses that approach.
