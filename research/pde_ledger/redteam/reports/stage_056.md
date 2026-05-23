---
unit_id: 056
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 056 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.txt`

## What the script claims to verify

The pair of scripts claim that (1) the ansatz `sigma = exp(v_sigma·s/D_sigma)` makes the steady-state drift-diffusion flux `J = -D_sigma·sigma' + v_sigma·sigma` vanish; (2) the normalized exponential source `sigma_Pe = Pe·exp(Pe·s/L)/(exp(Pe)-1)` integrates to `L` on `[0,L]`; (3) the lowest D/N mode overlap `Omega_Pe = I_Pe/I_W` with `chi0 = sqrt(2/L)·sin(pi·s/(2L))` reduces in closed form to `pi·Pe·(2·Pe·exp(Pe)+pi)/((4·Pe^2+pi^2)·(exp(Pe)-1))`; (4) `Omega_Pe(0) = 1` and `lim_{Pe→∞} Omega_Pe = pi/2`; (5) `dOmega_Pe/dPe = Cov(chi0, s)/I_W` under the normalized density `p_Pe = sigma_Pe/L`; and (6) the small-Pe linear coefficient is `(4-pi)/(2*pi)` while the large-Pe expansion gives `Omega_Pe ~ pi/2 - pi^3/(8*Pe^2)`.

## Assertion inventory

| #   | Script      | Line | Form                                                                                | Anchored to claim? |
| --- | ----------- | ---- | ----------------------------------------------------------------------------------- | ------------------ |
| A1  | sympy       | 40   | `expect_zero("zero-flux transport residual", J)`                                    | yes                |
| A2  | sympy       | 45   | `expect_zero("normalization int sigma_Pe ds - L", norm - L)`                        | yes                |
| A3  | sympy       | 61   | `expect_zero("Omega_Pe - expected formula", Omega_Pe - Omega_expected)`             | yes                |
| A4  | sympy       | 68-69| `if Omega0 != 1: raise`                                                             | yes                |
| A5  | sympy       | 70-71| `if simplify(OmegaInf - pi/2) != 0: raise`                                          | yes                |
| A6  | sympy       | 79   | `expect_zero("dOmega/dPe - Cov/I_W", diff(Omega_Pe, Pe) - Cov/I_W)`                 | yes                |
| A7  | sympy       | 85-88| `expect_zero("small-Pe linear coefficient", ...)`                                   | yes                |
| A8  | sympy       | 92-95| `expect_zero("large-Pe asymptotic through O(Pe^-2)", ...)`                          | yes                |
| M1  | mathematica | 34   | `expectZero["zero-flux transport residual", jFlux]`                                 | yes                |
| M2  | mathematica | 50   | `expectZero["normalization int sigma_Pe ds - ell", norm - ell]`                     | yes                |
| M3  | mathematica | 51   | `expectZero["Omega_Pe - expected formula", omegaPe - omegaExpected]`                | yes                |
| M4  | mathematica | 57   | `expectZero["twin baseline", omega0 - 1]`                                           | yes                |
| M5  | mathematica | 58   | `expectZero["upper finite-throat overlap limit", omegaInf - Pi/2]`                  | yes                |
| M6  | mathematica | 65   | `expectZero["dOmega/dPe - Cov/I_W", D[omegaPe,pe] - cov/iW]`                        | yes                |
| M7  | mathematica | 75-78| `expectZero["small-Pe linear coefficient", ...]`                                    | yes                |
| M8  | mathematica | 79-82| `expectZero["large-Pe second-order coefficient + Pi^3/8", largeCoeff + Pi^3/8]`     | yes                |

## Findings

No findings.

## Independent-derivation check (Mathematica)

The Mathematica script uses the same variable choreography as the SymPy script (`sigma = Exp[vSigma s/dSigma]`, `sigmaPe = pe Exp[pe s/ell]/(Exp[pe]-1)`, identical `chi0`, identical hardcoded `omegaExpected`, same intermediate quantities `iW, iPe, pPe, eChi, eS, eChiS, cov`, same assertion order). This is structurally a port. However, each engine independently performs the symbolic integrations (`Integrate` vs `sp.integrate`) and the algebraic simplification used to bring `iPe/iW` into the expected form is delegated to each engine's own simplifier (`FullSimplify` vs `simplify(expand(...))`). The two engines also diverge on the large-Pe verification: SymPy uses `sp.series(Omega_Pe, Pe, oo, 3).removeO()` and matches against `pi/2 - pi**3/(8*Pe**2)`, whereas Mathematica uses `Limit[pe^2 (omegaPe - Pi/2), pe -> Infinity]` and matches against `-Pi^3/8`. The closed-form `omegaExpected` is hardcoded in both engines but is verified against the integrated `iPe/iW`, so neither engine confirms the result merely by quoting it. On balance, the residual independent step (each engine doing its own definite integral of `sigma_Pe · chi0` and its own simplification, plus the different mechanism for the large-Pe coefficient) is enough to avoid a `mathematica_transliteration` finding, though the surface structural similarity is high.

## Engine cross-check

Both engines agree:
- `I_W = 2*sqrt(2)*sqrt(L)/pi` (sympy) / `(2*Sqrt[2]*Sqrt[ell])/Pi` (mathematica).
- `I_Pe = 2*sqrt(2)*sqrt(L)*Pe*(2*Pe*exp(Pe)+pi)/((4*Pe^2+pi^2)*(exp(Pe)-1))` (sympy) / `(2*Sqrt[2]*Sqrt[ell]*pe*(2*E^pe*pe+Pi))/((-1+E^pe)*(4*pe^2+Pi^2))` (mathematica).
- `Omega_Pe = pi*Pe*(2*Pe*exp(Pe)+pi)/((4*Pe^2+pi^2)*(exp(Pe)-1))` (sympy) / `(pe*Pi*(2*E^pe*pe+Pi))/((-1+E^pe)*(4*pe^2+Pi^2))` (mathematica).
- `Omega_Pe(0) = 1`, `Omega_Pe(∞) = pi/2` in both engines.
- Small-Pe series: sympy `-4*Pe**2/pi**2 + Pe**2/12 + Pe**2/pi - Pe/2 + 2*Pe/pi + 1`; mathematica `1 - pe/2 + pe^2/12 - (4*pe^2)/Pi^2 + (2*pe)/Pi + pe^2/Pi`. Term-by-term identical.
- Large-Pe leading correction: both confirm coefficient `-pi^3/8` on `Pe^-2` (sympy via direct series match against `pi/2 - pi^3/(8*Pe^2)`; mathematica via `Limit[pe^2(Omega-pi/2)]`).
- `dOmega/dPe - Cov/I_W = 0` in both.

All cross-checks consistent. Output files are fresh: sympy output mtime 2026-05-11 12:43 is newer than script mtime 2026-04-01 12:39; mathematica output mtime 2026-05-11 12:53 is newer than script mtime 2026-05-11 11:56.

## Verdict justification

The audit verdict is `clean`. Attacks tried and failed: (i) the zero-flux residual check looks tautological at first because `sigma = exp(v·s/D)` is the exact solution, but it would catch a typo in the exponent or sign, so it is a meaningful sanity check rather than a no-op; (ii) the closed-form `Omega_expected` is hardcoded but is verified by both engines against `Integrate[sigma_Pe · chi0]/I_W`, not against itself; (iii) the small-Pe assertion only extracts the linear coefficient, but the expected linear coefficient `(4-pi)/(2*pi) = 2/pi - 1/2` is matched non-trivially by the engine-produced series (`-Pe/2 + 2*Pe/pi`); (iv) the large-Pe assertion uses two distinct mechanisms across engines (series-removeO vs `Limit[pe^2 (Omega-pi/2)]`) and both arrive at `-pi^3/8`, eliminating the risk of a coincidental simplifier artifact; (v) the covariance identity `dOmega/dPe = Cov(chi0, s)/I_W` is derivable from the exponential-family score and is a non-trivial structural identity, not a definitional rewriting; (vi) symbol domains are uniformly positive-real and consistent with the integration intervals; (vii) the limits `Omega(0)=1` and `Omega(∞)=pi/2` both depend on the full closed-form structure and are independently attainable from the explicit `Omega_Pe`. The structural similarity between the two engines is real but does not cross into transliteration: the verification mechanisms (integration, series vs limit, simplifier paths) are engine-native.

## Self-test notes

Variable independence: every `sp.diff` / `D[...]` argument depends on the differentiation variable (`sigma` depends on `s`, `Omega_Pe` depends on `Pe`). Parity: all integrations are over `[0, L]` (asymmetric), so no even/odd cancellation trap applies. Trivial-case pre-check: substituting `Pe = 0` into `Omega_Pe` requires the limit (0/0), but the explicit limit and the series both give the constant term `1`, and the large-Pe coefficient `-pi^3/8` follows from `Limit[pe^2 (Omega - pi/2)]` independently of the series.
