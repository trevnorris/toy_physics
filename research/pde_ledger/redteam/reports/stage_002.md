---
unit_id: 002
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00-06:00
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage002_breathing_reduction.md
  paper_appendix: present
---

# Audit unit 002 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_002.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage002_breathing_reduction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 26 covers stage 002)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage002_breathing_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.txt`

## What the paper claims

Stage 002 verifies that the distributed wall lift, restricted to the lowest axisymmetric truncation, reproduces the old `(δa, δL)` matrix system, and that the grouped real P_2 sector becomes a one-mode degenerate family of the same wall PDE. The paper's `\stagefield{Output}` states: "The stage outputs the `(δa, δL)` matrix system (eq:app-stage002-el), the effective matrices (eq:app-stage002-mass-matrix)-(eq:app-stage002-stiffness-matrix), and the uncoupled grouped P_2 mass/stiffness pair (eq:app-stage002-m2)-(eq:app-stage002-k2)." Concretely the deliverables are: (D1) the Y_00 normalization bridge `q_00 = 2 sqrt(pi) δa`; (D2) the boxed `M_AB = 4π ∫dw μ_η α_A α_B` and `K_AB = 4π ∫dw [T_w α_A' α_B' + K_0 α_A α_B]`; (D3) the conservative Euler-Lagrange matrix equation `M_AB Q̈^B + K_AB Q^B = 0`; (D4) the grouped real P_2 ansatz with `-Δ_{S²}Y_{2m}^real = 6 Y_{2m}^real`; (D5) the boxed P_2 pair `M_2 = ∫dw μ_η β_2²` and `K_2 = ∫dw [T_w β_2'² + (K_η + 6 T_Ω) β_2²]` (no 4π prefactor on the P_2 side, since the real Y_{2m}^real basis is L² normalized); (D6) the P_2 Euler-Lagrange equation `M_2 q̈_{2m} + K_2 q_{2m} = 0`; (D7) the five-component degeneracy of grouped real P_2. The `\stagefield{Checks}` list also requires explicit verification of the 4π prefactor's provenance from the 2√π normalization, of the `l(l+1) = 6` angular eigenvalue, and of the no-dependence on real component label of `M_2/K_2`.

## What the script claims to verify

The SymPy and Mathematica scripts share a common three-section structure: (Section I) the Y_00 monopole normalization bridge — they verify `∫_{S²} Y_00 dΩ / (4π) = 1/(2 sqrt(π))`, `∫_{S²} Y_00² dΩ = 1`, that mouth-averaging `q_00 Y_00` yields `q_00/(2 sqrt(π))`, and that `(2 sqrt(π) Y_00)² ∫dΩ = 4π`. (Section II) the two-mode axisymmetric reduction — they form the reduced 1D Lagrangian density by integrating the Stage 001 quadratic wall action over Ω with the two-mode ansatz, then check that it equals the bilinear form built from the boxed `4π` overlap matrices, and additionally derive the Euler-Lagrange equations for `q_a, q_L`. The Mathematica script adds an independent `Coefficient[...]`-based extraction of the M and K integrands and compares them directly to the boxed forms. (Section III) the grouped real P_2 sector — they construct the explicit five real Y_{2m} basis, verify pairwise orthonormality (norm matrix = I_5), verify the angular stiffness matrix equals 6 I_5, verify `-Δ_{S²} basis_i = 6 basis_i` per component, and assemble the reduced density to check the boxed `M_2/K_2` form and the per-component degeneracy. Both engines then derive the single-component P_2 Euler-Lagrange equation.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 — Y_00 bridge `q_00 = 2√π δa` | normalization_bridge_audit (sympy 69-97); subbanner I (math 77-105) | match |
| D2 — boxed `M_AB`, `K_AB` with 4π | `lw - lw_target` (sympy 156-157); `lred - lred_target` (sympy 171-173); `MintegrandExtracted - MintegrandBoxed` and `LredTarget` checks (math 134-170) | match |
| D3 — `M_AB Q̈^B + K_AB Q^B = 0` | `euler_equations(lred_time, qa(t), [t])` checks (sympy 194-204); `EulerEquations` (math 189-200) | match |
| D4 — `-Δ_{S²} Y_{2m}^real = 6 Y_{2m}^real` | `lap_s2` per-basis check (sympy 250); `lapS2` per-basis check (math 251-257) | match |
| D5 — boxed `M_2 = ∫dw μ_η β_2²`, `K_2 = ∫dw [T_w β_2'² + (K_η + 6 T_Ω) β_2²]` | `lw_p2 - lw_p2_target` (sympy 264-282); per-i `lwP2i - lwP2Target` (math 268-297); `M2`, `K2` Integrals (sympy 286-287; math 303-304) | match |
| D6 — `M_2 q̈ + K_2 q = 0` | `single-component P2 Euler-Lagrange equation` (sympy 289-290; math 306-307) | match |
| D7 — five-component degeneracy | `norm_matrix - I_5`, `grad_matrix - 6 I_5` (sympy 247-248); `normMatrix5 - IdentityMatrix[5]`, `gradMatrix5 - 6 IdentityMatrix[5]` (math 234-249); per-component reduced density (math 268-297) | match |
| Check: 4π provenance from 2√π normalization | `angular prefactor from (2 sqrt(pi) Y00)^2 - 4 pi` (sympy 91-97; math 101-105) | match |
| Check: l(l+1) = 6 angular stiffness | `grad_matrix - 6 I_5` and per-basis `-Δ_{S²} y - 6y` (sympy 247-250; math 219-257) | match |
| Check: M_2/K_2 independent of real component label | per-i `lwP2i - lwP2Target` over all 5 components (math 268-297); collective form via norm/grad-matrix-as-I_5/6I_5 (sympy 247-282) | match |

Dominant pattern: `match` across all deliverables. Setting `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 83 | `Y00 mouth average - 1/(2 sqrt(pi)) = 0` | D1 (normalization bridge) | yes |
| A2 | sympy | 84 | `norm(Y00) - 1 = 0` | D1 | yes |
| A3 | sympy | 88 | `mouth average from q00 Y00 - q00/(2 sqrt(pi)) = 0` | D1 | yes |
| A4 | sympy | 97 | `angular prefactor (2 sqrt(pi) Y00)^2 - 4 pi = 0` | D2 prefactor check | yes |
| A5 | sympy | 157 | `lw - lw_target = 0` (reduced Lagrangian density vs boxed overlap form) | D2 | yes |
| A6 | sympy | 173 | `lred - lred_target = 0` (formal reduced Lagrangian vs boxed matrix form) | D2 | yes |
| A7 | sympy | 197-200 | Euler-Lagrange for q_a | D3 | yes |
| A8 | sympy | 201-204 | Euler-Lagrange for q_L | D3 | yes |
| A9 | sympy | 247 | `norm_matrix - eye(5) = 0` | D7 (and prerequisite for D5/D6) | yes |
| A10 | sympy | 248 | `grad_matrix - 6 eye(5) = 0` | D7 / l(l+1)=6 check | yes |
| A11 | sympy | 250 | `-Δ_S2 basis[i] - 6 basis[i] = 0` (5 checks) | D4 | yes |
| A12 | sympy | 282 | `lw_p2 - lw_p2_target = 0` (grouped reduced density) | D5 + D7 | yes |
| A13 | sympy | 290 | `single-component P2 Euler-Lagrange = 0` | D6 | yes |
| B1 | math | 94 | `Y00 mouth average - 1/(2 sqrt(pi))` | D1 | yes |
| B2 | math | 95 | `norm(Y00) - 1` | D1 | yes |
| B3 | math | 98 | `mouth avg q00 Y00 - q00/(2 sqrt(pi))` | D1 | yes |
| B4 | math | 105 | `angular prefactor (2 sqrt(pi) Y00)^2 - 4 pi` | D2 prefactor | yes |
| B5 | math | 150 | `extracted M - boxed M (4π overlap)` | D2 (M structure) | yes |
| B6 | math | 151 | `extracted K - boxed K (4π overlap)` | D2 (K structure) | yes |
| B7 | math | 170 | `Lred - LredTarget = 0` (formal reduced Lagrangian) | D2 | yes |
| B8 | math | 193-196 | Euler-Lagrange for q_a | D3 | yes |
| B9 | math | 197-200 | Euler-Lagrange for q_L | D3 | yes |
| B10 | math | 216 | phase shift Y21s vs Y21c | D7 (independent basis check) | yes |
| B11 | math | 217 | phase shift Y22s vs Y22c | D7 (independent basis check) | yes |
| B12 | math | 228 | per-i `norm(Y_i) - 1` (5 checks) | D7 | yes |
| B13 | math | 229 | per-i `angular energy(Y_i) - 6` (5 checks) | D4 / l(l+1)=6 | yes |
| B14 | math | 248 | `normMatrix5 - I_5` | D7 | yes |
| B15 | math | 249 | `gradMatrix5 - 6 I_5` | D7 | yes |
| B16 | math | 252-255 | per-i `-Δ_S2 Y_i - 6 Y_i` (5 checks) | D4 | yes |
| B17 | math | 292-295 | per-i `lwP2i - lwP2Target` (5 checks) | D5 + D7 | yes |
| B18 | math | 307 | `single-component P2 Euler-Lagrange` | D6 | yes |

Every script-side assertion traces to a specific paper-side deliverable. No orphaned scaffolding.

## Findings

(none)

## Independent-derivation check (Mathematica)

The Mathematica script is not a transliteration. Key independent moves:

- **Y_00 source**: Mathematica uses the builtin `SphericalHarmonicY[0, 0, theta, phi]` (line 83) while SymPy constructs `1/(2 sqrt(pi))` by hand (line 74). Both agree on the result, but the derivations come from independent definitional sources.
- **Section II M/K extraction**: Mathematica uses `Coefficient[lw, dadt, 2]`, `Coefficient[Coefficient[lw, dadt], dLdt]` etc. (lines 134-141) to extract the M and K integrands directly from the reduced Lagrangian density, then compares against the boxed form (lines 150-151). SymPy instead constructs the boxed `M_integrand`, `K_integrand` matrices explicitly and checks that the integrated reduced Lagrangian equals the bilinear form built from them (sympy line 156-157, 171-173). These are genuinely different algebraic angles: SymPy verifies "boxed M gives back lw"; Mathematica verifies "lw gives back boxed M". Either direction would catch a mismatched coefficient.
- **Section III per-component checks**: Mathematica iterates over `basis` with a `Do` loop and verifies per-component norm, angular energy, reduced density, and Laplacian eigenvalue (lines 219-297). SymPy uses matrix-level constructs (`norm_matrix`, `grad_matrix` as 5×5 sympy matrices). Mathematica also includes phase-shift identities (lines 216-217) that are absent from SymPy and provide an additional sanity check on the explicit basis encoding.
- **Euler-Lagrange**: Mathematica uses `EulerEquations` from `VariationalMethods` (loaded at line 3); SymPy uses `sympy.calculus.euler.euler_equations`. Different library implementations.

These differences support the claim of independent re-derivation rather than line-by-line port.

## Engine cross-check

Both engines run to completion and all assertions pass:

- SymPy output (`scripts/output/moving_throat_pde_stage002_breathing_reduction_sympy_audit.txt`): every `expect_zero` reports `= 0`. Final section summarises the four claim families.
- Mathematica output (`mathematica/output/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.txt`): every `expectZero` reports `PASS: ...`. Final summary identical in substance.

Concrete symbolic side-by-side checks where both engines compute the same quantity:

| Quantity | SymPy result | Mathematica result | Agree? |
|---|---|---|---|
| `∫ Y_00 dΩ / (4π) - 1/(2√π)` | 0 | 0 | yes |
| `∫ Y_00² dΩ - 1` | 0 | 0 | yes |
| `∫ (2√π Y_00)² dΩ - 4π` | 0 | 0 | yes |
| reduced Lagrangian = boxed bilinear form | 0 | 0 (also verified by Coefficient extraction) | yes |
| Euler-Lagrange for q_a | 0 | 0 | yes |
| Euler-Lagrange for q_L | 0 | 0 | yes |
| `norm_matrix - I_5` | 0 (matrix) | 0 (matrix) | yes |
| `grad_matrix - 6 I_5` | 0 (matrix) | 0 (matrix) | yes |
| `-Δ_S² Y_2m - 6 Y_2m` (×5) | 0 each | 0 each | yes |
| Grouped P_2 reduced density vs boxed form | 0 | 0 (per component i) | yes |
| Single-component P_2 EL | 0 | 0 | yes |

No engine disagreement.

Output freshness: sympy script mtime `2025-04-21 17:04`, sympy output mtime `2026-05-21 11:25`; mathematica script mtime `2026-05-21 00:39`, mathematica output mtime `2026-05-21 11:50`. Outputs are newer than the corresponding scripts. No `stale_output`.

## Verdict justification

The audit attempted several attacks and all failed:

1. **Tautology attempt**: Could the section-II `lw - lw_target` check pass for any wrong M_AB / K_AB? No — `lw` is computed from the action density via independent ansatz substitution and angular integration; `lw_target` is built from the separately-declared `M_integrand` and `K_integrand`. The match is a structural verification of the reduction, not a self-comparison.
2. **Hidden component degeneracy**: Could the per-i `lwP2i - lwP2Target` checks in Mathematica be passing because of common factoring? No — each iteration substitutes the specific Y_{2m}^real basis function explicitly into the action integrand and integrates over dΩ. The five components have visibly different angular forms (Y20 ∝ 3cos²θ − 1; Y22c ∝ sin²θ cos2φ; etc.), so a structural error in the angular orthonormality would surface in at least one component.
3. **4π prefactor cancellation**: Could the boxed 4π factor be a free parameter the script forgets to constrain? No — `angular prefactor (2 sqrt(pi) Y00)^2 - 4π = 0` directly verifies that the prefactor follows from the 2√π normalization in the ansatz, and the reduced-Lagrangian check then propagates that 4π into the boxed M_AB/K_AB.
4. **Symbol-assumption errors**: SymPy declares all relevant symbols `real=True`; Mathematica declares the corresponding `Element[..., Reals]` plus `0 < theta < Pi`, `wL < wR`. The assumptions match the physical setup (real perturbations, finite axial window, polar angle in (0, π)).
5. **Convention disagreement (paper vs script)**: The paper writes `η_0 = 2√π [α_a δa + α_L δL]` and the script writes `eta = 2√π · axisym · Y_00`. These differ only in whether the Y_00 factor is implicit (paper) or explicit (script); after the angular integration both produce the same 4π prefactor on M_AB/K_AB. The script's explicit treatment is mathematically equivalent (since Y_00 = 1/(2√π) makes `2√π · axisym · Y_00 = axisym`, and `∫dΩ axisym² = 4π · axisym²`).
6. **M_2/K_2 prefactor**: The paper's P_2 boxed forms have NO 4π prefactor (since `∫ (Y_{2m}^real)² dΩ = 1` by the script-verified normalization). The scripts correctly omit the 4π in M_2 and K_2 (sympy lines 286-287; mathematica lines 303-304). This is the asymmetry that distinguishes the monopole branch from the P_2 branch and the scripts honor it.
7. **Mathematica transliteration**: Mathematica uses `SphericalHarmonicY` builtin, `Coefficient`-based M/K extraction, and per-i Do-loop checks (phase shifts included) that have no sympy counterpart. Not a transliteration.

Paper-side and script-side claims align item-for-item across all seven deliverables and all four explicit `\stagefield{Checks}` items. Both engines pass. Outputs are fresh. No findings.

## Self-test notes

- Variable independence: every derivative used in the EL checks operates on functions that genuinely depend on `t` (qa(t), qL(t), q(t)), and every `D[..., w]` operates on functions of `w` (alphaA, alphaL, beta2, mu_eta, T_w, K_eta, T_Omega). No identically-zero derivatives lurking in `assert_zero` blocks.
- Symmetry/parity: the angular integrals over the explicit Y_{2m}^real basis use orthogonal combinations of sin/cos of φ and θ; cross-products (e.g., `∫ Y_21c Y_22s dΩ`) involve odd-in-φ integrands over [0, 2π] which vanish, consistent with the off-diagonal zeros in the verified `norm_matrix`.
- Trivial-case pre-check: substituting `alpha_a = 1`, `alpha_L = 0`, mu_eta=1, T_w=0, K_0=0 reduces the reduced Lagrangian density to `(1/2) · 4π · δȧ²`, matching boxed `M_aa = 4π ∫dw`, so the structure works out by hand. The check is non-trivial.
- Paper round-trip: re-reading the paper card and notes against the scripts shows no `paper_misalignment`. No new finding to introduce.
