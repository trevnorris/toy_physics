---
unit_id: 056
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage056_transport_source_asymmetry.md
  paper_appendix: present
---

# Audit unit 056 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_056.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage056_transport_source_asymmetry.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 90 and `\input{stages/stage_056}` at line 230)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_056.tex` lines 7-40) states Stage 056 gives the exponential source profile a physical transport origin via a one-dimensional drift-diffusion conservation law on the finite throat. The `\stagefield{Output}` (line 39) is verbatim: "Transport interpretation \eqref{eq:app-stage056-Pe} and physical overlap factor \eqref{eq:app-stage056-Omega-Pe}", which concretely sets:

- `Pe = v_sigma L / D_sigma` (boxed eq `app-stage056-Pe`)
- `Omega_Pe = pi*Pe*(2*Pe*exp(Pe) + pi) / [(4*Pe^2 + pi^2)*(exp(Pe) - 1)]` (boxed eq `app-stage056-Omega-Pe`)

The card also states that on the zero-flux stationary branch `J = 0`, `sigma(s) ∝ exp(v_sigma s / D_sigma)`, and that the normalized stationary family equals Stage 053's exponential source with `alpha = Pe`. The notes (`moving_throat_pde_stage056_transport_source_asymmetry.md`, sections 2-5) augment the paper card with: normalization to total source strength L; the explicit `sigma_Pe = Pe exp(Pe s/L)/(exp(Pe)-1)`; the lowest D/N mode `chi_0 = sqrt(2/L) sin(pi s/(2L))`; endpoint values `Omega_Pe(0) = 1` and `lim_{Pe → ∞} Omega_Pe = pi/2`; the covariance/monotonicity identity `dOmega_Pe/dPe = Cov_Pe(chi_0,s)/I_W`; and the asymptotics `Omega_Pe = 1 + ((4-pi)/(2 pi)) Pe + O(Pe^2)` and `Omega_Pe = pi/2 - pi^3/(8 Pe^2) + O(Pe^-3)`. The part-03 appendix row 90 summarizes: "Drift-diffusion zero-flux branch and Peclet-number identification."

## What the script claims to verify

Both scripts verify, in lockstep: (1) the zero-flux residual `J = -D_sigma diff(sigma,s) + v_sigma sigma` vanishes on the ansatz `sigma = exp(v_sigma s/D_sigma)`; (2) the normalized branch `sigma_Pe = Pe exp(Pe s/L)/(exp(Pe)-1)` integrates to L over [0,L]; (3) the overlap with chi_0 yields exactly the paper's closed-form `Omega_Pe`; (4) `Omega_Pe(0) = 1` and `lim_{Pe → ∞} Omega_Pe = pi/2`; (5) the covariance identity `dOmega_Pe/dPe = Cov(chi_0,s)/I_W` holds; (6) the small-Pe linear coefficient is `(4-pi)/(2 pi)`; (7) the large-Pe correction is `-pi^3/(8 Pe^2)`. The Mathematica script uses `Limit[pe^2 (omegaPe - Pi/2), pe -> Infinity]` for the large-Pe coefficient instead of SymPy's `series(...).removeO()` form.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Zero-flux ⇒ sigma ∝ exp(v_sigma s / D_sigma) (card eq `app-stage056-transport`) | sympy:38-40, wl:32-34 | match (ansatz verified on the ODE) |
| Normalization to L (notes §2) | sympy:43-45, wl:36-37,50 | match |
| `Omega_Pe = pi Pe (2 Pe e^Pe + pi)/[(4 Pe^2 + pi^2)(e^Pe-1)]` (boxed in card) | sympy:52-61, wl:40-51 | match |
| `Omega_Pe(0) = 1`, `lim → pi/2` (notes §5) | sympy:64-71, wl:53-58 | match |
| Covariance/monotonicity identity (notes §4) | sympy:74-79, wl:60-65 | match |
| Small-Pe linear coeff `(4-pi)/(2 pi)` (notes §5) | sympy:82-88, wl:67-68,75-78 | match |
| Large-Pe `−pi^3/(8 Pe^2)` (notes §5) | sympy:90-95, wl:69-70,79-82 | match |
| `Pe = v_sigma L / D_sigma` (definitional, boxed in card) | implicit in the substitution `v_sigma s / D_sigma → Pe s / L` between the J=0 branch and `sigma_Pe` | match (definitional, no numerical identity to check) |

Every paper/notes deliverable has a non-trivial script-side check. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `expect_zero("zero-flux transport residual", J)` | zero-flux ⇒ exponential | partial (ansatz verification, not derivation) |
| A2 | sympy | 45 | `expect_zero("normalization int sigma_Pe ds - L", norm - L)` | normalization | yes |
| A3 | sympy | 61 | `expect_zero("Omega_Pe - expected formula", Omega_Pe - Omega_expected)` | Omega_Pe closed form | yes |
| A4 | sympy | 68-69 | `if Omega0 != 1: raise` | endpoint Omega(0)=1 | yes |
| A5 | sympy | 70-71 | `if simplify(OmegaInf - pi/2) != 0: raise` | endpoint Omega(∞)=pi/2 | yes |
| A6 | sympy | 79 | `expect_zero("dOmega/dPe - Cov/I_W", ...)` | covariance identity | yes |
| A7 | sympy | 85-88 | `expect_zero("small-Pe linear coefficient", ...)` | small-Pe asymptote | yes |
| A8 | sympy | 92-95 | `expect_zero("large-Pe asymptotic through O(Pe^-2)", ...)` | large-Pe asymptote | yes |
| M1 | mathematica | 34 | `expectZero["zero-flux transport residual", jFlux]` | zero-flux ⇒ exponential | partial (ansatz verification) |
| M2 | mathematica | 50 | `expectZero["normalization int sigma_Pe ds - ell", norm - ell]` | normalization | yes |
| M3 | mathematica | 51 | `expectZero["Omega_Pe - expected formula", omegaPe - omegaExpected]` | Omega_Pe closed form | yes |
| M4 | mathematica | 57 | `expectZero["twin baseline", omega0 - 1]` | endpoint Omega(0)=1 | yes |
| M5 | mathematica | 58 | `expectZero["upper finite-throat overlap limit", omegaInf - Pi/2]` | endpoint Omega(∞)=pi/2 | yes |
| M6 | mathematica | 65 | `expectZero["dOmega/dPe - Cov/I_W", D[omegaPe,pe] - cov/iW]` | covariance identity | yes |
| M7 | mathematica | 75-78 | `expectZero["small-Pe linear coefficient", ...]` | small-Pe asymptote | yes |
| M8 | mathematica | 79-82 | `expectZero["large-Pe second-order coefficient + Pi^3/8", largeCoeff + Pi^3/8]` | large-Pe asymptote | yes |

## Findings

None.

Items considered and rejected:

- **A1/M1 zero-flux residual as `tautological_check`**: Picking `sigma = exp(v_sigma s/D_sigma)` and verifying `J = 0` is verifying an ansatz solves the ODE rather than deriving the solution from `J = 0`. In isolation it would be a weak `tautological_check`, but the immediately following normalization check (A2/M2) pins the normalizer `Pe/(exp(Pe)-1)` uniquely. The pair exercises the full claim. Not a finding.

- **Banner/docstring legacy labels "Stage 39"**: The SymPy script docstring says "Moving-Throat PDE — Stage 39 SymPy audit" (line 2) and the banner says "STAGE 39 — TRANSPORT ORIGIN OF THE SOURCE-SHAPE ASYMMETRY" (line 31); the final print says "Stage 39 audit passed." (line 97). The Mathematica banner says "STAGE 039 — TRANSPORT ORIGIN OF THE SOURCE-SHAPE ASYMMETRY" (line 26) though the final print says "Stage 056 Mathematica audit passed." (line 85). The notes themselves use legacy "Stage 38/39" terminology in headings, and the paper card uses "Stage 056". This is cosmetic banner/comment drift from a renumbering, not a math defect; substantive identities verified by the script match the paper's `\stagefield{Output}`. Not a `paper_misalignment` (no claim mismatch). Future cleanup, but not a v2 finding.

- **Symbol assumption `Pe > 0`**: The notes mention both Pe ≥ 0 (constructive) and Pe < 0 (destructive) branches. The script restricts to `Pe > 0`. The paper card's `\stagefield{Output}` makes no Pe-sign claim beyond what the analytic formula expresses, and downstream Stage 057 carries Pe ≥ 0 (per `\stagefield{Downstream use}` line 40). Restricting to `Pe > 0` matches the paper's actual downstream use. Not a finding.

- **`Omega_small` series carrying spurious quadratic terms**: SymPy output shows `-4*Pe**2/pi**2 + Pe**2/12 + Pe**2/pi - Pe/2 + 2*Pe/pi + 1` (sympy_output line 22) but the assertion only checks the *linear* coefficient, which is what notes §5 asserts (`O(Pe^2)` is the leading correction order, not a specific coefficient). Script does not over-claim against the notes. No finding.

- **`omegaExpected` hardcoded in both engines**: Both scripts hardcode `omegaExpected` and then check `omegaPe - omegaExpected == 0`. But `omegaPe` is computed independently via `Integrate[sigmaPe chi0]/iW`, so the assertion is non-trivial (it could fail if the integration didn't match the hardcoded formula). Not a `hardcoded_result` finding.

## Independent-derivation check (Mathematica)

The `.wl` mirrors the `.py` in variable choreography (`sigma`, `sigmaPe`, `chi0`, `iW`, `iPe`, `pPe`, `eChi`, `eS`, `eChiS`, `cov`; same assertion order). However the substantive computations are engine-native:

- Definite integrals (`Integrate` vs `sp.integrate`) are performed by each engine independently from the same physical premises (`sigma_Pe`, `chi_0`). Both arrive at matching closed forms.
- For the large-Pe asymptotic coefficient, the engines diverge in *technique*: `.py` uses `sp.series(Omega_Pe, Pe, sp.oo, 3).removeO()` and matches against `pi/2 - pi**3/(8 Pe**2)`; `.wl` uses `Limit[pe^2 (omegaPe - Pi/2), pe -> Infinity]` and matches against `-Pi^3/8`. Genuinely different mechanisms.
- For the small-Pe coefficient, both use Series-then-extract-coefficient, which is unavoidable since both engines support this idiom natively as the cleanest path.
- The covariance identity uses the definitional split `Cov = E[chi_0 s] - E[chi_0] E[s]` in both, which is the *mathematical* definition, not a transliteration of algebra.

The two engines reach matching closed forms via at-least-partially-distinct simplification paths, and the large-Pe verification uses different idioms. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical closed forms:

- `I_W = 2*sqrt(2)*sqrt(L)/pi` (sympy txt line 15) ↔ `(2*Sqrt[2]*Sqrt[ell])/Pi` (wl txt line 15). Agree.
- `I_Pe = 2*sqrt(2)*sqrt(L)*Pe*(2*Pe*exp(Pe) + pi)/((4*Pe**2 + pi**2)*(exp(Pe) - 1))` (sympy txt 16) ↔ `(2*Sqrt[2]*Sqrt[ell]*pe*(2*E^pe*pe + Pi))/((-1 + E^pe)*(4*pe^2 + Pi^2))` (wl txt 16). Agree.
- `Omega_Pe(0) = 1`, `lim Omega_Pe = pi/2` in both (sympy txt 19-20; wl txt 26-27).
- Small-Pe series (sympy txt 22) ↔ wl txt 38: term-by-term identical.
- Large-Pe coefficient `-pi^3/8` (sympy txt 25 verifies the full `O(Pe^-2)` form; wl txt 40 verifies the coefficient via Limit).
- `dOmega/dPe - Cov/I_W = 0` in both (sympy txt 21; wl txt 33).

All assertions PASS in both transcripts. `engines_agree: true`.

Output freshness:
- SymPy script mtime 2026-04-01 12:39; SymPy output mtime 2026-05-11 12:43 — fresh.
- Mathematica script mtime 2026-05-11 11:56; Mathematica output mtime 2026-05-11 12:53 — fresh.

`outputs_fresh: true`.

## Verdict justification

The scripts encode the paper's two boxed equations (Pe definition and Omega_Pe closed form) and all the additional deliverables enumerated in the notes (normalization, endpoints, covariance/monotonicity, small/large-Pe asymptotics). Each assertion is non-tautological except for the constructive ansatz check on J=0, which the immediately following normalization assertion makes substantive (a wrong normalizer would fail A2/M2). Both engines independently integrate `sigma_Pe * chi_0` and produce identical closed forms; the Mathematica script uses a Limit-based technique for the large-Pe coefficient that differs from SymPy's Series-removeO, providing some independent-derivation hardening. Both transcripts are fresh relative to script mtimes. I attacked: (a) the zero-flux residual as tautology — defused by A2/M2 normalizer pin; (b) the `Pe > 0` symbol restriction — defused by paper's constructive-branch-only downstream use; (c) the small-Pe series including spurious higher-order terms — defused by the assertion checking only the linear coefficient that the notes assert; (d) potential Mathematica transliteration — defused by engine-native Integrate and the divergent large-Pe Limit technique; (e) `omegaExpected` hardcoded — defused by independent integration on the LHS. Nothing breaks. Verdict: `clean`. I confirm I read the paper stage card, the notes file, and the appendix row, and the script's verified claim matches the paper's `\stagefield{Output}`.

## Self-test notes

(1) Variable independence: every `sp.diff(Omega_Pe, Pe)` / `D[omegaPe, pe]` is on an expression that explicitly depends on Pe via `Pe exp(Pe)/(exp(Pe)-1)` and `Pe/(4 Pe^2 + pi^2)` — no silent zero-derivative trap. (2) Parity: all integrals are over [0,L] (asymmetric); no even/odd cancellation trap applies. (3) Trivial-case substitution: at Pe → 0 the formula gives `lim pi Pe (2 Pe e^Pe + pi)/[(4 Pe^2 + pi^2)(e^Pe - 1)] = pi * pi / (pi^2 * 1) = 1`, matching the twin baseline; at Pe → ∞, dominant term `pi Pe * 2 Pe e^Pe / (4 Pe^2 e^Pe) = pi/2`, matching the upper limit. (4) Path specifications: N/A (no missing-script directive). (5) Paper round-trip: N/A (no fix prescribed; verdict is clean).
