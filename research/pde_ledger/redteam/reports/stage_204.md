---
unit_id: 204
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor.md]
  paper_appendix: present
---

# Audit unit 204 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_204.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows for stage 204: line 39 summary table; §"Free-quintuple log rays" lines 594-653; §"First-order and second-order root predictors" lines 655-674)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` states: "Graph-invariant log-ray family, primitive free-direction table, scalar root equation \(\Phi_{\mathbf s}(\tau)=1\), monotone-ray uniqueness theorem, and first local root predictors." The notes (§Purpose, items 1-7) enumerate seven deliverables: (1) the exact free-quintuple log ray \(\mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s}\); (2) the exact finite graph lift carrying the dependent triple via four explicit exponents \(\sigma_\delta,\sigma_T,\sigma_{K_\eta},\sigma_\mu\); (3) the exact theorem that the graph-lifted ray stays on the target orbit \(\mathcal O_*\) for all \(\tau\) (finite monomial invariance, equivalently \(M_*\dot{\Delta\mathbf x}=0\)); (4) the primitive direction table for the five coordinate directions; (5) the scalarized closure function \(\Phi_{\mathbf s}(\tau)\); (6) the monotone-ray uniqueness theorem; (7) the affine/log-linear root predictors \(\tau_{\rm aff}=(1-\Phi_0)/\Phi_1\), \(\tau_{\log}=-\ln\Phi_0/L_0\), with first-order agreement \(\tau_{\log}-\tau_{\rm aff}=-\varepsilon^2/(2L_0)+O(\varepsilon^3)\). The appendix (lines 594-674) gives the same exponent formulas verbatim: \(\mathfrak a_*=(1+\delta_{U,*})/(1+\chi_{0,*})\), \(\sigma_\delta=-\mathfrak a_*(s_\gamma+s_c-s_U)\), \(\sigma_T=s_U+\sigma_\delta\), \(\sigma_{K_\eta}=2s_c-s_U\), \(\sigma_\mu=2s_c-s_U+2s_W-2s_\lambda-E_*(2s_\gamma+2s_\lambda-s_U-s_W)+F_*\sigma_\delta\). Items 6 (monotone-ray uniqueness, an IVT/monotonicity existence theorem) is not a computational identity and is not expected to be a script assertion.

## What the script claims to verify

The SymPy script defines the constant log-slopes (Sec. I), constructs each dependent graph quantity (`delta_graph`, `T_graph`, `Keta_graph`, `mu_graph`) two ways — by direct substitution of \(\mathbf y_{\mathbf s}(\tau)\) into the closed graph formula, and as `const_0 * exp(sigma * tau)` — and asserts the two agree (Sec. II), thereby deriving each dependent exponent. It then substitutes the graph lift into the three target monomials \(\mathfrak C_{\rm tr},\mathfrak C_{\rm nt},\epsilon_\eta\) and asserts each equals its target constant for all \(\tau\) (Sec. III), checks the hardcoded quotient matrix \(M_*\) annihilates the ray tangent (Sec. IV), reproduces the five-row primitive direction table against expected literals (Sec. V), verifies the log-ray reproduces a general path to first order (Sec. VI), and series-expands \(\tau_{\log}-\tau_{\rm aff}\) confirming it starts at \(-\varepsilon^2/(2L_0)\) (Sec. VII). All checks route through `expect_zero`, which runs `simplify(expand(...))` and raises if nonzero. The script exits 0 with all checks passing.

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) log ray + constant slopes | Sec. I, lines 71-91 (`d/dtau ln X - s_X == 0`) | match |
| (2) finite graph lift, four exponents | Sec. II, lines 94-138 (direct vs exp form for δ,T,K_η,μ) | match |
| (3) monomial invariance / \(M_*\dot{\Delta x}=0\) | Sec. III lines 141-156, Sec. IV lines 159-171 | match |
| (4) primitive direction table | Sec. V, lines 174-197 (all 5 rows × 4 exponents) | match |
| (5) scalar closure \(\Phi_{\mathbf s}(\tau)\) | implicit (monomial invariance is the operative content; \(\widehat\chi_Q\) not symbolically reconstructed) | partial |
| (6) monotone-ray uniqueness theorem | none (analytic existence theorem, not a symbolic identity) | n/a (not script-expected) |
| (7) affine/log-linear predictors + first-order agreement | Sec. VII, lines 207-223 | match |
| Mathematica second engine | (absent) | missing |

`paper_alignment: aligned` — every computational deliverable maps to a faithful, non-tautological SymPy check, with the only gaps being deliverables (5)/(6) that are correctly not symbolic-identity targets (item 5's operative content is carried by the Sec. III monomial-invariance checks; item 6 is an IVT/monotonicity existence theorem with no closed-form residual).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1-A5 | sympy | 87-91 | `simplify(d/dtau ln X - s_X) == 0` | claim 1 | yes |
| A6 | sympy | 109 | `delta_direct - delta_exp == 0` | claim 2 (σ_δ) | yes |
| A7 | sympy | 116 | `T_direct - T_exp == 0` | claim 2 (σ_T) | yes |
| A8 | sympy | 123 | `Keta_direct - Keta_exp == 0` | claim 2 (σ_Kη) | yes |
| A9 | sympy | 138 | `mu_direct - mu_exp == 0` | claim 2 (σ_μ) | yes |
| A10-A12 | sympy | 154-156 | `monomial(tau) - target == 0` | claim 3 (invariance) | yes |
| A13 | sympy | 171 | `M_* · dx_ray == 0` (matrix) | claim 3 (kernel form) | yes |
| A14-A33 | sympy | 197 | 5 rows × 4 exponents `actual - expected == 0` | claim 4 (table) | yes |
| A34 | sympy | 204 | `series_match - (Y0 + Y1·eps) == 0` | claim 1 (local completeness) | yes |
| A35 | sympy | 220-223 | `series[τ_log-τ_aff] + eps²/(2L0) - 2eps³/(3L0) == 0` | claim 7 | yes |

Every row is non-tautological: the LHS and RHS in each `expect_zero` are produced from genuinely different constructions (direct closed-form substitution vs. exponential reconstruction; matrix product vs. zero; series coefficient vs. expected literal). No row defines `x=expr` then asserts `x==expr`.

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- target: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 204 is non-status-only (`is_status_only_candidate: False`) and a non-checkpoint, and it computes a large body of independently verifiable symbolic math: four closed-form dependent exponents, a matrix-kernel identity, a five-row primitive-exponent table, and two series expansions. The dual-engine project contract requires a Mathematica `.wl` wherever Mathematica *can* independently verify the stage. Mathematica can verify every one of these claims natively (D[], Series + Coefficient/SeriesCoefficient, matrix multiplication, Simplify/FullSimplify). No `.wl` exists. The card itself states "Mathematica audit: none yet." (`stage_204.tex:11`), confirming the gap is known.

**Why this matters:**
A single-engine stage has no cross-check against a SymPy transcription or simplification artifact. The whole point of the second engine is to re-derive the same physical results by a different route; without it, an error masked by SymPy's `simplify` (e.g., a branch/assumption artifact in the `(...)**(1/(1+chi0))` powers in Sec. II) would go undetected.

**Required change:**
Add the Mathematica audit script per the directive's claim manifest (M1-M7). It must derive the results independently (different decomposition than the .py), not transliterate the SymPy algebra.

**Verification:**
`redteam exec-mathematica 204` produces a `.wl` at the target path that exits 0 with substantive `expectZero`-style checks for M1-M7, and the engine cross-check (Sec. "Engine cross-check" below) then becomes applicable.

## Independent-derivation check (Mathematica)

No `.wl` exists, so no transliteration assessment applies yet. The directive's anti-transliteration guard requires the new script to use a different decomposition than the SymPy file (see directive F1).

## Engine cross-check

Not applicable — only one engine present. This is itself finding F1.

## Verdict justification

The SymPy script holds up under attack against the paper. I tried to break it on: (a) tautology — every `expect_zero` compares two independently-built expressions, not `x` against itself; (b) the τ-exponent algebra for all four dependent quantities — δ, T, K_η, μ all reduce to the paper's σ formulas (verified by hand: the μ exponent collects to `2s_c-s_U+2s_W-2s_λ-E_*(2s_γ+2s_λ-s_U-s_W)+F_*σ_δ`, matching appendix line 620-624); (c) the predictor series — `τ_log-τ_aff` leading term is `-eps²/(2L0)`, matching notes §8.1 and appendix; (d) the primitive table — `e_λ`,`e_W` give `σ_δ=σ_T=σ_Kη=0` and only move `σ_μ`, matching the notes' "two immediate consequences"; (e) symbol domains — `chi0,deltaU,E,F` positive and the base point positive, which legitimizes the real-power simplifications in Sec. II. Paper alignment is exact for every computational deliverable. The only finding is the missing Mathematica engine (F1), mandated by the dual-engine contract because Mathematica genuinely can verify this stage. Verdict: `findings`, no stop-cold.

## Self-test notes

Trap 1 (variable independence): the Sec. I `diff(log(X), tau)` derivatives all act on `X = X0*exp(s*tau)`, which genuinely depends on `tau`, so each derivative is a real nonzero slope — no identically-zero-derivative trap. Trap 2 (parity): no unbounded-domain integrals in this stage, so no parity hazard. Trap 3 (trivial-case): hand-substituted the primitive directions (`e_λ → σ_μ=-2-2E_*`, `e_W → σ_μ=2+E_*`) and confirmed they reduce to the expected literals; confirmed `series[τ_log-τ_aff]` reduces to `-eps²/(2L0)+2eps³/(3L0)`. Trap 4 (paths): the directive names the full `.wl` target under `mathematica/`. Trap 5 (paper round-trip): the proposed fix adds an engine only and prescribes the same paper-stated constants (`a_* = (1+deltaU)/(1+chi0)` etc.), introducing no new paper_misalignment.
