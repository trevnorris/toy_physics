---
unit_id: 152
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage152_family1_actual_correction.md]
  paper_appendix: present
---

# Audit unit 152 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_152.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage152_family1_actual_correction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_152}` at line 1338; no separate row content)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.txt`

## What the paper claims

The stage computes the *actual* finite-mouth-profile first-order correction to the Family-1 canonical point, using the full compensated residual profile `R_*(x)` instead of the linearized tangent. The notes box four classes of deliverables: (1) the weighted residual covariances `Cov_*(c,R_*) ≈ 0.0648069687666328` and `Cov_*(K_q,R_*) ≈ 0.0388718368650403`, giving the moment shifts `δg_act ≈ -0.0648…`, `δS_act ≈ -0.0388…`; (2) the retuned point `δΠ_act ≈ 0.907084414842908`, `δT̂_{m,act} ≈ 0.271653979462338`, hence `Π_corr ≈ 2.41591392833607` and `T̂_{m,corr} ≈ 1.17313803363654`; (3) the one-step nonlinear Picard iterate moments `g_1 ≈ 0.684423574065325`, `S_1 ≈ 0.616333130570251` and the resulting `Π_1 ≈ 2.53914847609768`, `T̂_{m,1} ≈ 1.21036942084359`; (4) the effective interpolation parameters `λ_eff^(Π) ≈ 0.380487632771110`, `λ_eff^(T) ≈ 0.378939241176339`, i.e. `λ_eff ≈ 0.38` (broadening fraction `1-λ_eff ≈ 0.62`). The card's quoted Output is "Full mouth profile broadens the source and shifts the point to `(Π_corr, T̂_{m,corr})`." Status is ExactClosure/Numerical.

## What the script claims to verify

Both scripts numerically construct `R_*(x) = Σ_m^*·(4 T_s(x) − T_q(x)) − Π_* x`, build the normalized exponential weight `Σ_*(x)`, and compute the mean-centered weighted covariances `Cov_*(c,R_*)`, `Cov_*(K_q,R_*)` by quadrature. From these they form `δg, δS, δΠ = Cov_cR/g'_*, δT = A_T δg + B_T δS`, the corrected point `(Π_corr, T_corr)`, the one-step Picard moments `g_1, S_1` and `(Π_1, T_1)`, and the effective interpolation `λ_eff^(Π)`, `λ_eff^(T)`. The load-bearing assertions (`expect_close`/`expectApprox`) pin `δΠ_act`, `δT_act`, both `λ_eff` against the boxed targets (tol 1e-6), and a linearized-covariance consistency identity (tol 1e-9) cross-checking that the hardcoded `g_star` equals the quadrature `∫Σ_* c`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `Cov_*(c,R_*) ≈ 0.0648…` / `δg_act` | py L75/78, wl L81/84 printed; exercised through δΠ assertion | match |
| `Cov_*(K_q,R_*) ≈ 0.0389…` / `δS_act` | py L76/79, wl L82/85 printed; exercised through δT assertion | match |
| `δΠ_act ≈ 0.907084414842908` | py L117 / wl L126 `expect ... tol 1e-6` | match |
| `δT̂_{m,act} ≈ 0.271653979462338` | py L118 / wl L127 | match |
| `Π_corr ≈ 2.41591392833607` | py L89 / wl L95 print = Π_*+δΠ (asserted via δΠ) | match |
| `T̂_{m,corr} ≈ 1.17313803363654` | py L90 / wl L96 print = T_*+δT (asserted via δT) | match |
| `g_1 ≈ 0.684…`, `S_1 ≈ 0.616…` | py L100-101 / wl L100-101 print | match (print-only, but Π_1/T_1 derive from them) |
| `Π_1 ≈ 2.53914847609768`, `T̂_{m,1} ≈ 1.21036942084359` | py L108-109 / wl L107-108 print | match (print-only) |
| `λ_eff^(Π) ≈ 0.380487632771110` | py L119 / wl L128 assert | match |
| `λ_eff^(T) ≈ 0.378939241176339` | py L120 / wl L129 assert | match |

Dominant pattern: aligned. The one-step Picard outputs (`g_1, S_1, Π_1, T_1`) are print-only rather than asserted against their boxed targets, but they are not independent literals — they flow from the same `R_*` and `Σ_*` machinery whose load-bearing endpoints (`δΠ`, `λ_eff`) ARE asserted, and both engines reproduce them to ~13 digits (cross-engine agreement substitutes for a missing literal assert). Not a finding: the stage is Numerical/ledger and the core retuning claim is asserted; the Picard check is explicitly framed as "a useful nonlinear check," not the primary output.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 117 | `expect_close(deltaPi, 0.907084414842908, 1e-6)` | claim 2 (δΠ_act → Π_corr) | yes |
| A2 | sympy | 118 | `expect_close(deltaT, 0.271653979462338, 1e-6)` | claim 2 (δT → T_corr) | yes |
| A3 | sympy | 119 | `expect_close(lam_Pi, 0.380487632771110, 1e-6)` | claim 4 (λ_eff^Π) | yes |
| A4 | sympy | 120 | `expect_close(lam_T, 0.378939241176339, 1e-6)` | claim 4 (λ_eff^T) | yes |
| A5 | sympy | 123-125 | linearized-cov identity == δg, 1e-9 | claim 1 consistency + anchors hardcoded g_star | yes |
| A6 | mma | 126 | `expectApprox[deltaPi, 0.907084414842908, 1e-6]` | claim 2 | yes |
| A7 | mma | 127 | `expectApprox[deltaT, 0.271653979462338, 1e-6]` | claim 2 | yes |
| A8 | mma | 128 | `expectApprox[lambdaEffPi, 0.380487632771110, 1e-6]` | claim 4 | yes |
| A9 | mma | 129 | `expectApprox[lambdaEffT, 0.378939241176339, 1e-6]` | claim 4 | yes |
| A10 | mma | 125 | `expectApprox[linearizedDeltaG, deltaG, 1e-9]` | claim 1 consistency | yes |

Every assertion traces to a paper deliverable. No orphaned/extra checks.

## Findings

None. The stage holds up against all attacks attempted below.

## Independent-derivation check (Mathematica)

This is the key scrutiny for the 105–175 orchestrator-direct band. **The `.wl` is NOT a transliteration of the `.py`; it is materially more independent.** Two structural divergences prove genuine re-derivation rather than a line-by-line port:

1. **Canonical-branch constants.** The SymPy script *hardcodes* the upstream constants as `mp.mpf` string literals: `Pi_star`, `Sigma_m_star`, `g_star`, `S_star`, `gprime_star`, `AT`, `BT`, `Tstar` (py L27-34). The Mathematica script *re-derives all of them from first principles*: it defines closed-form `gFormula[p]`, `sFormula[p]`, the Family-1 radius `rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1]` and `gMinus = rF1 - Sqrt[1+rF1^2]/2`, then solves `FindRoot[gFormula[p] == gMinus]` for `pStar` (wl L31-36) and computes `gPrimeStar = D[gFormula,p]`, `aT`, `bT`, `tStar`, `sigmaMStar` analytically (wl L37-47). The two engines therefore start from different inputs (literals vs derivation) and agree — this is the opposite of the band's transliteration defect.

2. **λ_eff interpolation.** SymPy uses opaque hardcoded Stage-148 endpoints: `lam_Pi = (1.69941496131430 - deltaPi)/2.08240814741023` (py L112-113). Mathematica instead re-derives the two interpolation endpoints from closed forms `gU=2/π`, `sU=2 tanh(π/2)/π`, `gD=π/4`, `sD=(1+sinh(π/2))/(2cosh(π/2))`, forms `deltaPiU/deltaPiD`, and computes `lambdaEffPi = (deltaPi - deltaPiU)/(deltaPiD - deltaPiU)` (wl L110-119) — a structurally different (and even differently-oriented) formula that reproduces the same `0.380487632771111`. This cross-validates SymPy's hardcoded endpoints against an analytic route.

The residual `R_*`, weight `Σ_*`, and covariance machinery are necessarily the same physical objects in both engines, but that is the shared *physics*, not shared *algebra*. Verdict: genuine second engine, no `mathematica_transliteration`.

Canonical-radius watch (per prompt): `rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1]`. Algebra: `12·(37/20)² = 12·1369/400 = 4107/100`, so `rF1 = Sqrt[(4107 − 100π²)/(100π²)] = √(4107 − 100π²)/(10π)` — exactly the canonical Family-1 radius. No stale `168π²`/`168%` anywhere; the `≈0.38`/`≈0.62` percentages are the correct `λ_eff`/broadening values, unrelated to the `100π²` watch item.

## Engine cross-check

Both engines agree to ~13 significant digits on every deliverable:

| Quantity | SymPy (.txt) | Mathematica (.txt) |
|---|---|---|
| Cov_*(c,R_*) | 0.0648069687666328 | 0.06480696876663268… |
| δΠ_act | 0.907084414842908 | 0.90708441484290542… |
| Π_corr | 2.41591392833607 | 2.41591392833606100… |
| T_corr | 1.17313803363654 | 1.17313803363654159… |
| g_1 | 0.684423574065325 | 0.68442357406532438… |
| Π_1 | 2.53914847609768 | 2.53914847609767349… |
| λ_eff^(Π) | 0.380487632771111 | 0.38048763277111113… |
| λ_eff^(T) | 0.378939241176339 | 0.37893924117633769… |

All five SymPy `expect_close` and all five Mathematica `expectApprox` checks PASS in the committed transcripts (residuals ≤ ~2.6e-15). No `engine_disagreement`.

## Verdict justification

`clean`. I read the paper card, the notes, and the appendix include first, then attacked the scripts. Attacks tried and failed: (a) **transliteration** — refuted; the `.wl` re-derives the canonical constants and λ_eff endpoints by genuinely different routes than the `.py`'s hardcoded literals. (b) **tautology on the linearized-cov check** — it reduces to `−Cov(c,R) = δg`, but because `g_star` enters as an independently-hardcoded literal that the integral `∫Σ_*c` must reproduce, the check actually anchors the hardcoded `g_star` and would fail on a transcription error; not vacuous. (c) **hardcoded targets** — the assertion targets are the stage's own boxed deliverables, and the asserted quantities (`δΠ`, `δT`, `λ_eff`) are computed by quadrature from `R_*`, not re-stated literals; the cross-engine independent derivation removes the "number against itself" risk. (d) **measure normalization** — `∫₀¹Σ_* dx = 1` exactly, so the covariance weight is a proper probability measure and the mean-centering (`f − ⟨f⟩`) in `cov` correctly centers deformations (satisfies the card's "centered before covariance formulas" check). (e) **symbol/domain** — no symbolic `simplify` under risky assumptions; everything is `mpmath`/`NIntegrate` numerics with `Pi_star>0`, integrand smooth on [0,1]; the wl `kap²−pStar²` denominators in `tq`/`sFormula` are nonzero at `pStar≈1.509` (kap=π/2≈1.571). (f) **radius/percentage drift** — `rF1` is the correct `√(4107−100π²)/(10π)`; no `168π²` contamination. The paper claim and the script's verified claim match exactly. Outputs are fresh (sympy .txt 2026-05-28 11:30 > .py 10:01; mma .txt 2026-05-28 11:31 > .wl 2026-05-11).

## Self-test notes

Checked: (1) Variable independence — `D[gFormula[p],p]` and `D[sFormula[p],p]` (wl L39-40) depend on `p`, derivative nonzero at `pStar`; no zero-derivative trap. (2) Symmetry/parity — integrals are over the finite [0,1] domain, not symmetric/unbounded, so parity-cancellation traps don't apply; `Σ_*` normalizes to 1 (verified analytically). (3) Trivial-case — `R_*(0)=0` confirmed in both transcripts; mean-centering makes `Cov` a genuine non-trivial integral (≈0.0648, nonzero). No directive needed (zero findings).

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value emitted by the scripts was located in the stage notes (`.md`); the terse `.tex` card legitimately omits the numerics (it points to the block narrative), so the `.md` is the natural carrier and all matches there are MATCH per the guards.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Cov_*(c,R_*) = 0.0648069687666328 | py L75/out L9, wl L81/out L9 | notes:25 | MATCH |
| Cov_*(K_q,R_*) = 0.0388718368650403 | py L76/out L10, wl L82/out L10 | notes:30 | MATCH |
| δg_act = -0.0648069687666328 | py L78/out L11, wl L84/out L11 | notes:37 | MATCH |
| δS_act = -0.0388718368650403 | py L79/out L12, wl L85/out L12 | notes:39 | MATCH |
| δΠ_act = 0.907084414842908 | py L80/out L13 (assert L117), wl L86/out L13 (assert L126) | notes:62 | MATCH |
| δT̂_{m,act} = 0.271653979462338 | py L81/out L14 (assert L118), wl L87/out L14 (assert L127) | notes:68 | MATCH |
| Π_corr = 2.41591392833607 | py L89/out L15, wl L95/out L15 | notes:78 | MATCH |
| T̂_{m,corr} = 1.17313803363654 | py L90/out L16, wl L96/out L16 | notes:86 | MATCH |
| g_1 = 0.684423574065325 | py L100/out L17, wl L100/out L17 | notes:112 | MATCH |
| S_1 = 0.616333130570251 | py L101/out L18, wl L101/out L18 | notes:114 | MATCH |
| Π_1 = 2.53914847609768 | py L108/out L19, wl L107/out L19 | notes:121 | MATCH |
| T̂_{m,1} = 1.21036942084359 | py L109/out L20, wl L108/out L20 | notes:123 | MATCH |
| λ_eff^(Π) = 0.380487632771110 | py L112/out L21 (assert L119), wl L118/out L21 (assert L128) | notes:137 | MATCH |
| λ_eff^(T) = 0.378939241176339 | py L113/out L22 (assert L120), wl L119/out L22 (assert L129) | notes:139 | MATCH |
| Π_* = 1.50882951349316 | py L27 (wl re-derives via FindRoot) | notes:14 | MATCH |
| Σ_m^* = 0.451485277739090 | py L28 (wl L41 derives) | notes:16 | MATCH |
| g_* = 0.758035078944663 | py L29 (wl L38 derives) | notes:18 | MATCH |
| S_* = 0.658075937605429 | py L30 (wl L37 derives) | notes:20 | MATCH |
| g_*' = 0.0714453558083195 | py L31 (wl L39 derives) | notes:52 | MATCH |
| A_T = -4.27263956256927 | py L32 (wl L43 derives) | notes:54 | MATCH |
| B_T = 0.134875005736706 | py L33 (wl L47 derives) | notes:56 | MATCH |

INTERNAL (genuine scaffolding / imported upstream constants not expected as a stage-152 prose numeric): `Tstar = 0.901484054174205` (upstream `T̂_{m,*}`, referenced symbolically in notes:85, not a stage-152 deliverable literal); the Stage-148 interpolation endpoint literals `1.69941496131430`, `2.08240814741023`, `0.508756302215085`, `0.625700104366894` (upstream affine-law constants, correctly not printed as stage-152 results — the `.wl` re-derives the equivalent endpoints analytically); `R(0)=0`, `R(1)=-0.9048…`, `R'(0)` forward-diffs (residual sanity prints); `Z1` (Picard normalizer); pass/fail diffs and tolerances.

reconciliation: complete; 21 deliverable values checked, 0 misaligned
