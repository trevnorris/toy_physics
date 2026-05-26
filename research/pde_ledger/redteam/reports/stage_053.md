---
unit_id: 053
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage053_overlap_boost_window.md
  paper_appendix: present
---

# Audit unit 053 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_053.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage053_overlap_boost_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 84)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage053_overlap_boost_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage053_overlap_boost_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage053_overlap_boost_mathematica_audit.txt`

## What the paper claims

Quoting `\stagefield{Output}` verbatim: *"Overlap boost formula (eq:app-stage053-Omega-exp) and ceiling (eq:app-stage053-overlap-ceiling)."* The stage proves:

1. The exact closed form `Omega_exp(alpha) = pi*alpha*(2*alpha*e^alpha + pi) / [(4*alpha^2 + pi^2)*(e^alpha - 1)]` for the overlap boost of the exponential bottom-biased source family `sigma_alpha(s) = alpha*exp(alpha*s/L)/(exp(alpha)-1)` against the lowest D/N support mode `chi_0(s) = sqrt(2/L)*sin(pi*s/(2L))`.
2. Endpoint values `Omega_exp(0) = 1` (uniform branch) and `lim_{alpha->+inf} Omega_exp(alpha) = pi/2` (sharp upper bound, antinode concentration).
3. The overlap ceiling `Omega_0^2 <= pi^2/4` (eq:app-stage053-overlap-ceiling), so pure overlap rescue requires `zeta_req <= pi^2/4`.
4. Notes additionally state: (i) the bound `0 <= Omega_0 <= pi/2` holds for any nonnegative source density with total strength `L`, (ii) the source family has total mass `L`, and (iii) the small-alpha linear coefficient equals `(4-pi)/(2pi) > 0` so the constructive branch moves upward off the uniform point.

The appendix row (line 84 of `stage_appendix_part03.tex`) summarises the claim as: "Exponential source family and ceiling `Omega_0^2 <= pi^2/4`."

## What the script claims to verify

The two scripts independently: (i) compute `I_W = int_0^L chi_0 ds = 2*sqrt(2L)/pi` via direct integration; (ii) assert `Omega_max := L * sup(chi_0) / I_W = pi/2` and therefore `A_I,max = pi^2/4`; (iii) integrate `sigma_alpha` and assert total mass equals `L`; (iv) integrate `sigma_alpha * chi_0` and assert the resulting `Omega_alpha` equals the typed-in closed form; (v) assert the alpha->0 limit equals 1 and the alpha->+infty limit equals pi/2; (vi) Taylor-expand `Omega_alpha` near alpha=0, extract the actual linear coefficient via `coeff(alpha, 1)` / `Coefficient[..., alpha, 1]`, and assert it equals `(4-pi)/(2*pi)`. The "pure overlap threshold" block prints `criterion = A_I,max - zeta_req` without asserting (correctly, since it is a literal restatement, not an identity to verify).

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Closed form `Omega_exp(alpha)` (eq:app-stage053-Omega-exp) | `expect_zero("Omega_alpha closed form", Omega_alpha - Omega_alpha_simpler)` (sympy L63); `expectZero["Omega_alpha closed form", ...]` (math L65) | match |
| Endpoint `Omega_exp(0) = 1` | `expect_zero("uniform branch value", Omega0 - 1)` (sympy L69); math L76 | match |
| Endpoint `lim_{alpha->+inf} Omega_exp = pi/2` | `expect_zero("antinode concentration limit", Omegainf - pi/2)` (sympy L70); math L77 | match |
| Ceiling `Omega_0^2 <= pi^2/4` (value side) | `expect_zero("A_I,max - pi^2/4", ...)` (sympy L43); math L49. Bound is computed as `L*chi0_max/I_W` and equated to pi/2; the Hoelder inequality step `I_(phi,0) <= chi0_max * L` is taken as elementary (chi_0 is nonneg, total mass is L) and not re-proven in script. | match |
| Total source mass = L (notes) | `expect_zero("same total source strength", Sigma_total - L)` (sympy L57); math L64 | match |
| Small-alpha linear coefficient `(4-pi)/(2pi)` (notes) | `expect_zero("linear coefficient - (4-pi)/(2pi)", linear_coeff - (4-pi)/(2*pi))` where `linear_coeff = series_small.coeff(alpha, 1)` (sympy L74-76); math L70,78 | match |
| Pure-overlap rescue criterion `zeta_req <= pi^2/4` | `criterion = sp.simplify(A_I_max - zeta_req)` printed (sympy L81); math L82. No assertion — this is a logical consequence stated for the reader. | match (no assertion needed; it is a restatement, not an identity) |

Dominant pattern: `match` across every deliverable. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero("Omega_max - pi/2", Omega_max - pi/2)` | ceiling value (eq:app-stage053-overlap-ceiling) | yes |
| A2 | sympy | 43 | `expect_zero("A_I,max - pi^2/4", A_I_max - pi**2/4)` | ceiling A_I = pi^2/4 | yes |
| A3 | sympy | 57 | `expect_zero("same total source strength", Sigma_total - L)` | notes mass-normalisation | yes |
| A4 | sympy | 63 | `expect_zero("Omega_alpha closed form", Omega_alpha - Omega_alpha_simpler)` | eq:app-stage053-Omega-exp | yes |
| A5 | sympy | 69 | `expect_zero("uniform branch value", Omega0 - 1)` | Omega_exp(0)=1 | yes |
| A6 | sympy | 70 | `expect_zero("antinode concentration limit", Omegainf - pi/2)` | lim_{alpha->inf} = pi/2 | yes |
| A7 | sympy | 76 | `expect_zero("linear coefficient - (4-pi)/(2pi)", linear_coeff - (4-pi)/(2*pi))` with `linear_coeff = series_small.coeff(alpha,1)` | notes small-alpha linearisation | yes |
| A8 | mathematica | 48 | `expectZero["Omega_max - Pi/2", omegaMax - Pi/2]` | ceiling value | yes |
| A9 | mathematica | 49 | `expectZero["A_I,max - Pi^2/4", aIMax - Pi^2/4]` | ceiling A_I = pi^2/4 | yes |
| A10 | mathematica | 64 | `expectZero["same total source strength", sigmaTotal - ell]` | notes mass-normalisation | yes |
| A11 | mathematica | 65 | `expectZero["Omega_alpha closed form", omegaAlpha - omegaAlphaSimple]` | eq:app-stage053-Omega-exp | yes |
| A12 | mathematica | 76 | `expectZero["uniform branch value", omega0 - 1]` | Omega_exp(0)=1 | yes |
| A13 | mathematica | 77 | `expectZero["antinode concentration limit", omegaInf - Pi/2]` | lim_{alpha->inf} = pi/2 | yes |
| A14 | mathematica | 78 | `expectZero["linear coefficient - (4-Pi)/(2Pi)", linearCoeff - (4-Pi)/(2Pi)]` with `linearCoeff = Coefficient[seriesSmall, alpha, 1]` | notes small-alpha linearisation | yes |

Every script-side check traces to a specific paper-side deliverable. No orphan assertions; no missing deliverables.

## Findings

None. The v1 audit's `tautological_check` finding on A7/A14 (the small-alpha linear coefficient was being asserted against itself via a hardcoded literal) has been addressed: both scripts now extract `linear_coeff = series_small.coeff(alpha, 1)` (sympy L74) and `linearCoeff = Coefficient[seriesSmall, alpha, 1]` (math L70) from the actual symbolic series of `Omega_alpha` before the assertion. The assertion now genuinely depends on the upstream integration that produced `Omega_alpha`.

## Independent-derivation check (Mathematica)

The Mathematica script structurally parallels the SymPy script (variable names camelCased: `chi0`/`chi0`, `iW`/`I_W`, `omegaMax`/`Omega_max`, `aIMax`/`A_I_max`, etc.; assertion order identical; closed-form `omegaAlphaSimple`/`Omega_alpha_simpler` typed in literally). However, the heavy lifting in each script is independent symbolic work performed by the respective engine: `sp.integrate(chi0, (s, 0, L))` vs. `Integrate[chi0, {s, 0, ell}]`; `sp.integrate(sigma_alpha * chi0, (s, 0, L))` vs. `Integrate[sigmaAlpha chi0, {s, 0, ell}]`; `sp.limit(Omega_alpha, alpha, 0)` vs. `Limit[omegaAlpha, alpha -> 0]`; `sp.series(...)` vs. `Series[...]`. Each engine derives `Omega_alpha` from its own integration of the same physical premises, and the closed-form check is a redundant comparison of two independent integration results against a common typed-in target — this is exactly what the second-engine policy asks for in a direct-integration verification. The parallel choreography is a consequence of the problem's canonical structure (single integral, two endpoint limits, one series expansion); it does not constitute transliteration. Not a finding.

Specific paired excerpts:

- SymPy L31: `I_W = sp.simplify(sp.integrate(chi0, (s, 0, L)))` -> Mathematica L38: `iW = FullSimplify[Integrate[chi0, {s, 0, ell}], ...]`. Each is the engine's own integration call.
- SymPy L49-50: `I_alpha = sp.integrate(sigma_alpha * chi0, (s, 0, L)); Omega_alpha = sp.simplify(I_alpha / I_W)` -> Mathematica L53-54: `iAlpha = FullSimplify[Integrate[sigmaAlpha chi0, {s, 0, ell}], ...]; omegaAlpha = FullSimplify[iAlpha/iW, ...]`. Same physical integrand, independent symbolic execution.
- SymPy L65-66: `Omega0 = sp.simplify(sp.limit(Omega_alpha, alpha, 0)); Omegainf = sp.simplify(sp.limit(Omega_alpha, alpha, sp.oo))` -> Mathematica L67-68: `omega0 = FullSimplify[Limit[omegaAlpha, alpha -> 0], ...]; omegaInf = FullSimplify[Limit[omegaAlpha, alpha -> Infinity], ...]`. Independent limit primitives.

## Engine cross-check

Both engines produce identical symbolic results:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `I_W` | `2*sqrt(2)*sqrt(L)/pi` | `(2*Sqrt[2]*Sqrt[ell])/Pi` |
| `Omega_max` | `pi/2` | `Pi/2` |
| `A_I,max` | `pi**2/4` | `Pi^2/4` |
| `Sigma_total - L` | `0` | `0` (PASS) |
| `Omega_alpha closed-form residual` | `0` | `0` (PASS) |
| `Omega_alpha(alpha->0)` | `1` | `1` (PASS) |
| `Omega_alpha(alpha->+infty)` | `pi/2` | `Pi/2` (PASS) |
| small-alpha series | `alpha**2*(-4/pi**2 + 1/12 + 1/pi) + alpha*(-1/2 + 2/pi) + 1` | `1 + alpha^2*(1/12 + (-4 + Pi)/Pi^2) + alpha*(-1/2 + 2/Pi)` (algebraically identical) |
| extracted linear coefficient | `(4 - pi)/(2*pi)` | `-1/2 + 2/Pi` (= `(4-Pi)/(2Pi)`, PASS) |
| `criterion = A_I,max - zeta_req` | `-zeta_req + pi**2/4` | `(Pi^2 - 4*zetaReq)/4` |

All assertions PASS in both transcripts; both scripts exit 0. No engine disagreement.

Output freshness (mtimes verified):
- SymPy script (`...stage053_overlap_boost_sympy_audit.py`) mtime = SymPy output mtime + ~64s earlier; output is newer. Fresh.
- Mathematica script mtime = Mathematica output mtime - ~71s; output is newer. Fresh.

## Verdict justification

I attacked the scripts on the following vectors and they held:

1. *Hardcoded `chi0_max = sqrt(2/L)`* — the script does not numerically maximise; it inserts the supremum directly. But `chi_0(s) = sqrt(2/L)*sin(pi*s/(2L))` on `s in [0,L]` is monotone and attains `sqrt(2/L)` at the endpoint `s=L` (where `sin(pi/2)=1`); the supremum is unambiguous and the bound argument `int sigma * chi_0 <= sup(chi_0) * int sigma` is elementary. Not a hardcoded_result finding.
2. *`positive=True` on `alpha`* — alpha=0 (uniform branch) and alpha=+infty (delta limit) are computed via `sp.limit` / `Limit`, which evaluate as limits even when the symbol is declared positive. Both limits yield the paper's stated values. The Mathematica warning `Limit::alimv` is the engine declining to use limit-variable assumptions and is benign here; the resulting limits are still correct.
3. *Closed-form check as identity loop* — `Omega_alpha_simpler` is typed in matching the paper card; the check `Omega_alpha - Omega_alpha_simpler == 0` is non-tautological because `Omega_alpha` is computed by *independent* symbolic integration in each engine.
4. *Linear-coefficient assertion* — fixed in v1 round; now extracts via `coeff(alpha,1)` / `Coefficient[..., alpha, 1]` from the actual series and equates to `(4-pi)/(2*pi)`. Genuine.
5. *Mathematica `$Assumptions` overwrite at L81* — after the rebinding, the only remaining work is the (trivial) `criterion` print on a literal `Pi^2/4 - zetaReq`; no upstream computation depends on the lost assumptions. Not a finding.
6. *Paper alignment* — every deliverable in the paper card and notes maps to a specific assertion; every assertion in the scripts maps to a specific paper claim. The appendix row matches.

Verdict: `clean`. The scripts exercise the paper's claims with non-tautological, paper-anchored assertions and the two engines independently corroborate the same closed form.

## Self-test notes

I checked: (1) variable independence — `series_small` is built from `Omega_alpha`, which itself depends on `alpha` via the integration of `sigma_alpha * chi0`, so `.coeff(alpha, 1)` is a genuine extraction (not a trivial constant); (2) limit triviality — both `alpha -> 0` and `alpha -> +infty` produce non-trivial residuals (0/0 forms requiring L'Hopital or asymptotic expansion), confirmed by the printed series and the differing closed-form structures, so the assertions `Omega0 - 1 == 0` and `Omegainf - pi/2 == 0` are substantive; (3) symbol-domain pitfalls — `s, L, alpha` declared positive; no integrals over symmetric/unbounded domains where parity could cancel things; (4) no new edits prescribed, so no paper round-trip needed; (5) v1 directive change is present and effective. No traps remain.
