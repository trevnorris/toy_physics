---
unit_id: 002
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage002_breathing_reduction.md]
  paper_appendix: present
---

# Audit unit 002 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_002.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage002_breathing_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 26; intro paragraph line 9; `\input` line 83)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage002_breathing_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.txt`

## What the paper claims

Stage 002 is a checkpoint `\StatusExactClosure{}` stage proving that the distributed wall theory of Stage 001 reduces, in the lowest axisymmetric sector, back to the old two-coordinate `(a,L)` collective-geometry closure, and that the grouped real `P_2` lane is degenerate at the same conservative level. `\stagefield{Output}{The stage outputs the \((\delta a,\delta L)\) matrix system \eqref{eq:app-stage002-el}, the effective matrices \eqref{eq:app-stage002-mass-matrix}--\eqref{eq:app-stage002-stiffness-matrix}, and the uncoupled grouped \(P_2\) mass/stiffness pair \eqref{eq:app-stage002-m2}--\eqref{eq:app-stage002-k2}.}` Distinct deliverables: (1) the `Y_00` normalization bridge `Y_00=1/(2√π)`, `∫Y_00²dΩ=1`, `(1/4π)∫Y_00 dΩ=1/(2√π)` underlying `q_00=2√π δa`; (2) the `4π` overlap prefactor from the `2√π` ansatz normalization; (3) boxed `M_AB=4π∫dw μ_η α_A α_B`; (4) boxed `K_AB=4π∫dw[T_w α_A'α_B'+K_0 α_A α_B]` with `K_0=K_η`; (5) the conservative EL matrix system `M_AB Q̈^B+K_AB Q^B=0`; (6) boxed `M_2=∫dw μ_η β_2²` (NO 4π); (7) boxed `K_2=∫dw[T_w β_2'²+(K_η+6T_Ω)β_2²]`; (8) the single-component degenerate EOM `M_2 q̈_{2m}+K_2 q_{2m}=0` for every real `P_2`; plus the explicit check that the `l=2` angular stiffness shift is `l(l+1)T_Ω=6T_Ω`. The notes are fully consistent with the card and add the `q_00(0,t)=2√π δa(t)` mouth-average identification and the densitized `dw` (no extra `√γ_0`) measure convention.

## What the script claims to verify

Both scripts verify, in three sections, exactly these deliverables. Section I integrates the genuine real harmonic `Y_00` to confirm the average/normalization facts, the `q_00`-to-`q_00/(2√π)` mouth average, and that `(2√π Y_00)²` integrates to `4π`. Section II inserts the two-mode ansatz `η_0=2√π(α_a δa+α_L δL)Y_00` into the Stage-001 angular wall action density, performs the `dΩ` integral, and (SymPy) checks the reduced density equals the hand-built `4π` boxed `M`/`K` quadratic form, then checks the formal `dw`-integrated Lagrangian equals the boxed matrix form, then derives the EL equations and checks they equal `M Q̈+K Q`. The Mathematica version is stronger here: it EXTRACTS `M`,`K` from the reduced density via `Coefficient` and compares to the boxed forms, rather than asserting a pre-built form. Section III builds the five explicit real `l=2` harmonics, verifies orthonormality (5×5 identity), the angular gradient matrix `=6·Id`, the eigenvalue `-Δ_{S²}Y=6Y` for each, the per-component reduced density collapses to the common `(M_2,K_2)` form, and the single EOM `M_2 q̈+K_2 q=0`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `Y_00=1/(2√π)`, `∫Y_00²=1`, `(1/4π)∫Y_00=1/(2√π)`, `q_00=2√π δa` | py 76-89 / wl 83-99 | match |
| (2) `4π` overlap prefactor from `2√π` norm | py 91-97 / wl 101-105 | match |
| (3) boxed `M_AB=4π∫μ_η α_A α_B` | py 137-142,157,173 / wl 134-151,170 | match |
| (4) boxed `K_AB=4π∫[T_w α'α'+K_0 αα]`, `K_0=K_η` | py 143-154,157,173 / wl 146-151,170 | match |
| (5) EL matrix system `M Q̈+K Q=0` | py 194-204 / wl 189-200 | match |
| (6) boxed `M_2=∫μ_η β_2²` (no 4π) | py 286,288 / wl 303 | match |
| (7) boxed `K_2=∫[T_w β_2'²+(K_η+6T_Ω)β_2²]` | py 287 / wl 304 | match |
| (8) degenerate EOM `M_2 q̈_{2m}+K_2 q_{2m}=0`, all 5 comps | py 247-290 / wl 248-307 | match |
| Check: `l(l+1)T_Ω=6T_Ω` angular shift | py 248,250 / wl 249,256 + per-comp 6 | match |

Dominant pattern: every paper-side deliverable has a faithful, non-tautological script-side check in BOTH engines. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 83 | `simplify(avg_y00 - 1/(2√π))==0` | claim 1 | yes |
| A2 | sympy | 84 | `simplify(norm_y00 - 1)==0` | claim 1 | yes |
| A3 | sympy | 88 | `q00·avg - q00/(2√π)==0` | claim 1 | yes |
| A4 | sympy | 97 | `∫(2√π Y00)² - 4π==0` | claim 2 | yes |
| A5 | sympy | 157 | reduced density − boxed `4π` quad form `==0` | claims 3,4 | yes |
| A6 | sympy | 173 | `Integral(lw) − boxed matrix form==0` | claims 3,4 | yes |
| A7 | sympy | 197-204 | EL eqs `+ (MQ̈+KQ)==0` | claim 5 | yes |
| A8 | sympy | 247 | `norm_matrix − I₅==0` | claim 8 (orthonorm) | yes |
| A9 | sympy | 248 | `grad_matrix − 6 I₅==0` | claim 7,8,check | yes |
| A10 | sympy | 250 | `−Δ_{S²}Y − 6Y==0` (×5) | check `6=l(l+1)` | yes |
| A11 | sympy | 282 | P2 density − degenerate form `==0` | claims 6,7,8 | yes |
| A12 | sympy | 290 | single P2 EL `+ (M₂q̈+K₂q)==0` | claim 8 | yes |
| B1-B4 | wl | 94,95,98,105 | mirror of A1-A4 (genuine `SphericalHarmonicY`) | claims 1,2 | yes |
| B5 | wl | 150-151 | `Coefficient`-extracted M,K − boxed `==0` | claims 3,4 | yes (stronger) |
| B6 | wl | 170 | `Integral(lw) − boxed matrix==0` | claims 3,4 | yes |
| B7 | wl | 193-200 | EL eqs − `(MQ̈+KQ)==0` | claim 5 | yes |
| B8 | wl | 216-217 | phase-shift relations among real P2 | claim 8 (extra) | yes |
| B9 | wl | 228,229,248,249 | per-comp + 5×5 norm/grad matrices | claims 7,8,check | yes |
| B10 | wl | 252-256 | `−Δ_{S²}Y−6Y==0` (×5) | check `6=l(l+1)` | yes |
| B11 | wl | 292 | per-comp density − degenerate form (×5) | claims 6,7,8 | yes |
| B12 | wl | 307 | single P2 EL − `(M₂q̈+K₂q)==0` | claim 8 | yes |

## Findings

None. Attacks tried and the reasons they failed are in the verdict justification.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. The two engines diverge in their derivation strategy at the load-bearing step:

- Section I: SymPy hardcodes `y00 = 1/(2√π)` (py:74); Mathematica derives it from the library `y00 = FullSimplify[SphericalHarmonicY[0,0,theta,phi]]` (wl:83). This is an independent source for the same normalization.
- Section II (the core reduction): SymPy hand-constructs `M_integrand`/`K_integrand` as `4π·[...]` matrices (py:137-154) and asserts the angular-reduced density equals that quadratic form (py:157). Mathematica instead EXTRACTS the matrices from the reduced density via `Coefficient[lw, dadt, 2]` etc. (wl:134-141) and only then compares to the boxed `4π` form (wl:150-151). Extraction-then-compare is a genuinely different (and stronger) route than assert-against-hand-built-form.
- Section III: SymPy uses `−Δ_{S²}` and explicit real `Y_2m` to build the matrices; Mathematica additionally checks the real-harmonic phase-shift construction (`Y21s = Y21c|_{φ→φ−π/2}`, wl:216-217) which SymPy does not — an extra independent structural check, not an echo.

This satisfies the second-engine policy.

## Engine cross-check

Both saved outputs show all checks resolving to literal `0` / `PASS`, with no FAIL and a clean `Exit[0]` path implied (Mathematica) and no `AssertionError` (SymPy). The recovered overlap matrices printed by each engine agree symbolically: SymPy prints `4π∫α_a²μ_η dw` etc.; Mathematica prints `Integrate[4*Pi*alphaAF[w]^2*muEtaF[w], ...]` etc. The off-diagonal `K` entry in the Mathematica print appears as `Integrate[2*Pi*(2*alphaAF*alphaLF*K0F + 2*TwF*α_a'*α_L'), ...]` = `4π∫[...]` after pulling the factor — algebraically identical to the SymPy `4π∫[K0 α_a α_L + T_w α_a'α_L']` entry. No sign, factor-of-2, or `4π` discrepancy between engines. `engines_agree: true`.

## Verdict justification

Verdict `clean`. Attacks attempted, all failed:
(1) **`4π` vs no-`4π` mismatch trap.** The most likely defect for this stage is a stray (or missing) `4π` on `M_2`/`K_2`. The paper deliberately boxes `M_AB,K_AB` WITH `4π` (from `(2√π)²∫Y_00²=4π`) but `M_2,K_2` WITHOUT `4π` (bare `∫Y_2m²=1`). Both scripts honor this exactly: `M2 = Integral(mu_eta*beta2**2,...)` carries no `4π` (py:286, wl:303), and the angular integral in the P2 density (py:265-273, wl:272-284) uses normalized `Y_2m` whose `∫Y²=1` produces coefficient 1, not `4π`. Internally consistent and paper-aligned.
(2) **EL sign convention.** Both assert the EL output `+ (MQ̈+KQ)==0` (SymPy `euler_equations`, lhs `=Eq(...,0)`) / `− (MQ̈+KQ)==0` (Mathematica `EulerEquations` returning `lhs==rhs`, then `elA = rhs−lhs`). The opposite sign in the comparison is exactly the convention difference between the two libraries' EL outputs, and each yields residual 0 — non-tautological and correct.
(3) **Tautology in Section II.** SymPy A5 could be suspected tautological (assert hand-built `M_integrand` equals what you inserted), but `lw` is independently produced by performing the `dΩ` integral of the action density, so A5 genuinely tests that the angular reduction yields the `4π` form. Mathematica removes even this suspicion by extracting via `Coefficient`. Not tautological.
(4) **Symbol-domain / assumptions.** `0<θ<π`, `wL<wR`, reals — all physically justified (polar angle range, ordered integration limits, real collective coordinates). No `positive` assumption is used to force a `simplify` outcome; the `4π` and `6` are derived from genuine integrals.
(5) **Degeneracy breadth.** Claim 8 ("every real P2 component") is exercised over all five components in both engines (5×5 matrices + per-component density loops), not just one — no `missing_branch`/`insufficient_verification`.
Outputs are fresh (both `.txt` mtimes newer than their scripts). Paper card, notes, and appendix row were all read; the script's verified claim matches the paper's stated claim in every deliverable.

## Value Reconciliation (pass-2 augmentation)

Stage 002 emits no benchmark/figure-of-merit numbers; its deliverables are symbolic closed forms plus two structural constants (`4π`, `6`). Reconciliation is therefore over the symbolic boxed results and those two constants. Saved outputs are fresh and confirm every emitted form.

| value | source (py / wl + output line) | .tex / .md location | status |
|---|---|---|---|
| `Y_00 = 1/(2√π)` | py:74 / wl:83 (out py:5, wl:9-10) | tex:18 eq; md:18 | MATCH |
| `∫Y_00²dΩ = 1` | py:84 / wl:95 (out py:6) | tex:19; md:19 | MATCH |
| `(1/4π)∫Y_00 dΩ = 1/(2√π)` | py:83 / wl:94 (out py:5) | tex:20-22; md:20-21 | MATCH |
| `q_00 = 2√π δa` (mouth bridge) | py:89 print / wl:99 print (out py:8) | tex:26-27 eq:two-mode-ansatz; md:24-26 | MATCH |
| `4π` overlap prefactor | py:97 / wl:105 (out py:9, wl:16-17) | tex:48,54 boxed; md:83,87 | MATCH |
| `M_AB = 4π∫μ_η α_A α_B` | py:137-142 / wl:142-145,150 (out py:17-29, wl:31) | tex:47-48 (boxed eq:mass-matrix); md:82-83 | MATCH |
| `K_AB = 4π∫[T_w α'α'+K_0 αα]` | py:143-154 / wl:146-149,151 (out py:30-64, wl:32) | tex:53-56 (boxed eq:stiffness-matrix); md:86-92 | MATCH |
| `K_0 = K_η` | py:109 (`K_0` symbol) / wl:114 | tex:58; md:59-60 | MATCH |
| EL system `M_AB Q̈^B+K_AB Q^B=0` | py:197-204 / wl:193-200 (out py:65-66, wl:33-36) | tex:64-65 (boxed eq:el); md:104 | MATCH |
| `M_2 = ∫μ_η β_2²` (no 4π) | py:286 / wl:303 | tex:89-90 (boxed eq:m2); md:138-139 | MATCH |
| `K_2 = ∫[T_w β_2'²+(K_η+6T_Ω)β_2²]` | py:287 / wl:304 | tex:95-98 (boxed eq:k2); md:142-148 | MATCH |
| `−Δ_{S²}Y_2m = 6Y_2m` (l(l+1)=6) | py:250 / wl:252-256 (out py:91-95, wl:71-80) | tex:79,97,113; md:126,147 | MATCH |
| P2 EOM `M_2 q̈_{2m}+K_2 q_{2m}=0` (×5, degenerate) | py:282,290 / wl:292,307 (out py:96-97, wl:81-92) | tex:104-105 (boxed eq:p2-eom); md:152 | MATCH |
| real P2 5×5 orthonormality = `I₅` | py:247 / wl:248 (out py:71-80, wl:65-67) | tex:114 (degeneracy via label-independence); md:154 | MATCH |

INTERNAL (scaffolding, no finding): per-entry residuals all `0`; `norm_matrix`/`grad_matrix` intermediate forms; `lw`/`lw_target`, `lred`/`lred_target` comparison residuals; `Coefficient`-extracted `MaaExt…KaLExt` intermediates; phase-shift checks `Y21s−Y21c|_{φ→φ−π/2}`, `Y22s−Y22c|_{φ→φ−π/4}` (structural, mirrors the real-harmonic construction — present in `.wl` only, supports degeneracy but is not a separately boxed deliverable); pass/PASS flags; banners.

reconciliation: complete; 14 deliverable values checked, 0 misaligned.

## Self-test notes

Checked: (1) variable-independence — every `diff`/`D` in the EL and density checks is w.r.t. a variable the integrand actually depends on (`t` for fields `q_a(t),q_L(t),q(t)`; `w` for profiles `α,β_2`), so no spuriously-zero derivative drives a trivial pass. (2) Parity/symmetry of the angular integrals — `Y_2m²` and `gradS2Inner` integrands are even under the `dΩ` measure and give the nonzero `1`/`6`; off-diagonal `Y_iY_j` (i≠j) products are orthogonal and correctly give `0` (the 5×5 minus-identity residual is the all-zero matrix in both outputs). (3) Trivial-case: the `4π` arises only from `(2√π)²·∫Y_00²dΩ=4π·1`, and `M_2`'s coefficient is `1·∫Y_2m²=1` — confirming the deliberate presence/absence of `4π` between the two sectors is correct, not an oversight. No traps triggered; no directive written.
</content>
</invoke>
