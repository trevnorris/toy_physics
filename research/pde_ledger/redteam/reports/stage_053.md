---
unit_id: 053
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 053 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage053_overlap_boost_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage053_overlap_boost_mathematica_audit.txt`

## What the script claims to verify

The scripts compute the "exact overlap-boost window" for the lowest sine support lane `chi0(s) = sqrt(2/L) sin(pi s/(2L))` on `[0, L]`. They assert: (i) the Dirac-delta upper-bound value `Omega_max = L * sup(chi0) / I_W` equals `pi/2`, so `A_I,max = pi^2/4`; (ii) the explicit bottom-biased exponential family `sigma_alpha(s) = alpha exp(alpha s/L)/(exp(alpha)-1)` has total mass `L` (matching the uniform branch) and the inner product `Omega_alpha = <sigma_alpha, chi0>/I_W` matches the typed-in closed form `pi*alpha*(2*alpha*e^alpha + pi)/((4 alpha^2 + pi^2)(e^alpha - 1))`; (iii) `Omega_alpha` interpolates from `1` (alpha->0) to `pi/2` (alpha->infinity); (iv) the small-alpha linear coefficient equals `(4-pi)/(2pi)`. Together they bound the pure-overlap rescue threshold by `zeta_req <= pi^2/4`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero("Omega_max - pi/2", Omega_max - pi/2)` | yes |
| A2 | sympy | 43 | `expect_zero("A_I,max - pi^2/4", A_I_max - pi**2/4)` | yes |
| A3 | sympy | 57 | `expect_zero("same total source strength", Sigma_total - L)` | yes |
| A4 | sympy | 63 | `expect_zero("Omega_alpha closed form", Omega_alpha - Omega_alpha_simpler)` | yes |
| A5 | sympy | 69 | `expect_zero("uniform branch value", Omega0 - 1)` | yes |
| A6 | sympy | 70 | `expect_zero("antinode concentration limit", Omegainf - pi/2)` | yes |
| A7 | sympy | 76 | `expect_zero("linear coefficient - (4-pi)/(2pi)", expected_linear - (4 - pi)/(2*pi))` | no (tautological) |
| A8 | mathematica | 42 | `expectZero["Omega_max - Pi/2", omegaMax - Pi/2]` | yes |
| A9 | mathematica | 43 | `expectZero["A_I,max - Pi^2/4", aIMax - Pi^2/4]` | yes |
| A10 | mathematica | 58 | `expectZero["same total source strength", sigmaTotal - ell]` | yes |
| A11 | mathematica | 59 | `expectZero["Omega_alpha closed form", omegaAlpha - omegaAlphaSimple]` | yes |
| A12 | mathematica | 70 | `expectZero["uniform branch value", omega0 - 1]` | yes |
| A13 | mathematica | 71 | `expectZero["antinode concentration limit", omegaInf - Pi/2]` | yes |
| A14 | mathematica | 72 | `expectZero["linear coefficient - (4-Pi)/(2Pi)", linearCoeff - (4 - Pi)/(2 Pi)]` | no (tautological) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py:72-76`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl:63-72`

**What's wrong:**
In the SymPy script lines 72-76:
```
series_small = sp.series(Omega_alpha, alpha, 0, 3).removeO()
print("small-alpha series =", series_small)
expected_linear = sp.simplify(2 / pi - sp.Rational(1, 2))
print("linear coefficient =", expected_linear)
expect_zero("linear coefficient - (4-pi)/(2pi)", expected_linear - (4 - pi) / (2 * pi))
```
The variable `expected_linear` is defined as `2/pi - 1/2`. The subsequent assertion then checks `expected_linear - (4-pi)/(2*pi) == 0`, which is the pure algebraic identity `2/pi - 1/2 == (4-pi)/(2pi)`. This holds by elementary fraction arithmetic, completely independent of `Omega_alpha`, `chi0`, `sigma_alpha`, or the actual small-alpha series. The `series_small` object is computed and printed but never compared to anything in an assertion. Thus the "linear coefficient" claim — that the linear-in-alpha term of `Omega_alpha` equals `(4-pi)/(2pi)` — is not actually verified by an assertion.

The Mathematica script lines 63-72 mirror this exactly:
```
seriesSmall = FullSimplify[Normal[Series[omegaAlpha, {alpha, 0, 2}]], Assumptions -> ell > 0];
linearCoeff = FullSimplify[2/Pi - 1/2];
...
expectZero["linear coefficient - (4-Pi)/(2Pi)", linearCoeff - (4 - Pi)/(2 Pi)];
```
Here `linearCoeff` is a hardcoded literal `2/Pi - 1/2`, never extracted from `seriesSmall`. The assertion is again the trivial identity `2/Pi - 1/2 == (4-Pi)/(2Pi)`.

**Why this matters:**
If the actual coefficient of `alpha^1` in the series expansion of `Omega_alpha` were wrong (e.g., from an integration error or a typo in `Omega_alpha_simpler`), this assertion would still pass. The whole purpose of the small-alpha series check is to verify that the typed-in closed form has the correct linear behavior near the uniform branch — but as written, the assertion cannot detect such an error.

**Required change:**
Replace the hardcoded `expected_linear`/`linearCoeff` with the linear coefficient actually extracted from the symbolic series of `Omega_alpha`, then assert that this extracted coefficient equals `(4-pi)/(2pi)`.

SymPy (`scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py:72-76`):
- Compute `series_small = sp.series(Omega_alpha, alpha, 0, 3).removeO()` (already there).
- Replace `expected_linear = sp.simplify(2 / pi - sp.Rational(1, 2))` with `linear_coeff = sp.simplify(series_small.coeff(alpha, 1))`.
- Replace the trivial assert with `expect_zero("linear coefficient - (4-pi)/(2pi)", linear_coeff - (4 - pi) / (2 * pi))`.

Mathematica (`mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl:63-72`):
- Keep `seriesSmall = FullSimplify[Normal[Series[omegaAlpha, {alpha, 0, 2}]], Assumptions -> ell > 0];`.
- Replace `linearCoeff = FullSimplify[2/Pi - 1/2];` with `linearCoeff = FullSimplify[Coefficient[seriesSmall, alpha, 1], Assumptions -> ell > 0];`.
- The existing `expectZero["linear coefficient - (4-Pi)/(2Pi)", linearCoeff - (4 - Pi)/(2 Pi)]` then becomes substantive.

**Verification:**
After the fix, the assertion `linear coefficient - (4-pi)/(2pi) = 0` still holds (the actual linear coefficient is `2/pi - 1/2 = (4-pi)/(2pi)`), but it now genuinely depends on the integration that produced `Omega_alpha`. The verifier confirms the new code lines (`series_small.coeff(alpha, 1)` / `Coefficient[seriesSmall, alpha, 1]`) appear in place of the literal `2/pi - 1/2` / `2/Pi - 1/2`, and the script still exits 0.

## Independent-derivation check (Mathematica)

Both scripts share the same physical premises (the sine eigenfunction `chi0` and the exponential family `sigma_alpha`) and the same independently-typed closed form `Omega_alpha_simpler`/`omegaAlphaSimple`. The two engines independently perform the integrals `I_W = integral(chi0)` and `I_alpha = integral(sigma_alpha*chi0)` — SymPy via `sp.integrate`, Mathematica via `Integrate`. They also independently evaluate the alpha→0 and alpha→infinity limits via their respective limit primitives. This is genuine redundant computation, not transliteration; the heavy lifting (integrating an exponential-times-sine over a finite interval, and the limit computations) is done from scratch by each engine.

The variable choreography is parallel (same names camelCased) and the closed-form check is typed in identically, but this is the expected pattern for an independent re-verification of the same symbolic identity. Not a transliteration finding.

## Engine cross-check

Both engines produce numerically/symbolically identical final outputs:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `I_W` | `2*sqrt(2)*sqrt(L)/pi` | `(2*Sqrt[2]*Sqrt[ell])/Pi` |
| `Omega_max` | `pi/2` | `Pi/2` |
| `A_I,max` | `pi**2/4` | `Pi^2/4` |
| `Sigma_total - L` | `0` | `0` |
| `Omega_alpha closed form` residual | `0` | `0` |
| `Omega_alpha(alpha->0)` | `1` | `1` |
| `Omega_alpha(alpha->+infty)` | `pi/2` | `Pi/2` |
| small-alpha series linear term | `-1/2 + 2/pi` | `-1/2 + 2/Pi` |

All assertions pass in both engines (exit 0). No engine disagreement.

Output freshness:
- SymPy script mtime: 2026-04-01 12:39; SymPy output mtime: 2026-05-11 12:43 (fresh).
- Mathematica script mtime: 2026-05-11 11:56; Mathematica output mtime: 2026-05-11 12:52 (fresh).

## Verdict justification

The substantive claims (Dirac-delta bound `Omega_max = pi/2`, exponential-family integral closed form, alpha->0 and alpha->infinity limits, total mass `L`) are each backed by non-tautological assertions in both engines, and both engines independently confirm them. The only defect is the "small-alpha linear coefficient" check, which in both engines compares a hardcoded `2/pi - 1/2` against the algebraically equivalent `(4-pi)/(2pi)` — never touching the actual series expansion of `Omega_alpha`. The fix is mechanical (replace the literal with `series_small.coeff(alpha, 1)` / `Coefficient[seriesSmall, alpha, 1]`). No stop-cold conditions apply: the fix does not change any quoted downstream numerical result (the linear coefficient really is `(4-pi)/(2pi)` once correctly extracted).

## Self-test notes

I checked: (1) variable independence — the proposed `series_small.coeff(alpha, 1)` and `Coefficient[seriesSmall, alpha, 1]` operate on the existing `series_small`/`seriesSmall` which already depend on `alpha` via `Omega_alpha`, so the extraction is well-defined and non-trivial; (2) trivial-case pre-check — substituting the printed `series_small = alpha**2*(...) + alpha*(-1/2 + 2/pi) + 1` into the new `coeff(alpha, 1)` yields `-1/2 + 2/pi`, and `-1/2 + 2/pi - (4-pi)/(2pi) = (-pi + 4 - 4 + pi)/(2pi) = 0`, so the rewritten assertion still passes on the current correct integration; (3) path specifications — both files exist at the named paths and only require in-place edits, no new file creation. Parity / integral-symmetry traps do not apply (no unbounded integrals; the symmetric-domain pitfall is absent here).
