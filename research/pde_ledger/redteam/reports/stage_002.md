---
unit_id: 002
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: false
  outputs_fresh: true
---

# Audit unit 002 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage002_breathing_reduction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.txt`

(Both output files have mtimes later than their script mtimes: scripts dated 2026-04-21 17:04, outputs dated 2026-05-11 12:38 / 12:42. Outputs are fresh.)

## What the script claims to verify

The two scripts together claim to verify, for the breathing reduction stage of the moving-throat PDE: (i) a monopole normalization bridge identifying `q_00 = 2 sqrt(pi) delta_a` via the `Y_00` mouth average and the resulting `4 pi` angular prefactor; (ii) that inserting the axisymmetric two-mode ansatz `eta = 2 sqrt(pi)(alpha_a delta_a + alpha_L delta_L) Y_00` into the Stage 001 wall action reproduces the boxed overlap-integral form of the kinetic matrix `M_AB` and stiffness matrix `K_AB`, with the resulting two-variable conservative Euler-Lagrange system `M_AB Q_ddot^B + K_AB Q^B = 0`; (iii) orthonormality of the grouped real P2 multiplet, a common Laplace-Beltrami eigenvalue 6 across the multiplet, and the resulting isotropic five-component degeneracy with single-component reduced equation `M2 q_ddot + K2 q = 0`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 83 | `simplify(avg_y00 - 1/(2 sqrt(pi))) == 0` | yes |
| A2 | sympy | 84 | `simplify(norm_y00 - 1) == 0` | yes |
| A3 | sympy | 88 | `simplify(q00 * avg_y00 - q00/(2 sqrt(pi))) == 0` | partial (algebraic rescaling of A1) |
| A4 | sympy | 97 | `simplify(angular_prefactor - 4 pi) == 0` | yes |
| A5 | sympy | 157 | `simplify(lw - lw_target) == 0` (angular integration of two-mode ansatz vs. hand-built overlap matrices with 4 pi factor) | yes |
| A6 | sympy | 173 | `simplify(lred - lred_target) == 0` (integration commutes with matrix bilinear-form expansion) | partial (follows from A5 by linearity of integral) |
| A7 | sympy | 197-200 | Euler-Lagrange for q_a vanishes when matched against `Maa qa'' + MaL qL'' + Kaa qa + KaL qL` | yes |
| A8 | sympy | 201-204 | Euler-Lagrange for q_L vanishes (mirror of A7) | yes |
| A9 | sympy | 247 | `(norm_matrix - I_5) == 0` (real-P2 orthonormality across all five components) | yes |
| A10 | sympy | 248 | `(grad_matrix - 6 I_5) == 0` (common stiffness eigenvalue 6 across all five components) | yes |
| A11 | sympy | 250 | `-Delta_S2 basis[i] - 6 basis[i] == 0` for i = 0..4 (Laplace-Beltrami eigen-check) | yes |
| A12 | sympy | 282 | grouped-five-component reduced density - diagonal target form | yes |
| A13 | sympy | 290 | single-component P2 Euler-Lagrange equation matches `M2 q'' + K2 q` | yes |
| B1 | mathematica | 93 | `Y00 mouth average - 1/(2 sqrt(pi))` | yes |
| B2 | mathematica | 94 | `norm(Y00) - 1` | yes |
| B3 | mathematica | 97 | `mouth average from q00 Y00 - q00/(2 sqrt(pi))` | partial (algebraic rescaling of B1) |
| B4 | mathematica | 104 | `angular prefactor - 4 Pi` | yes |
| B5 | mathematica | 146 | `lw - lwTarget == 0` (same construction as A5) | yes |
| B6 | mathematica | 158 | `Lred - LredTarget == 0` | partial (follows from B5) |
| B7 | mathematica | 180-183 | Euler-Lagrange q_a check | yes |
| B8 | mathematica | 184-187 | Euler-Lagrange q_L check | yes |
| B9 | mathematica | 203 | `Y21s - Y21c(phi - Pi/2) == 0` (substitute for testing y21s) | partial (phase identity, not orthonormality) |
| B10 | mathematica | 204 | `Y22s - Y22c(phi - Pi/4) == 0` (substitute for testing y22s) | partial (phase identity, not orthonormality) |
| B11 | mathematica | 215 | `norm(basis[i]) - 1` for i in {Y20, Y21c, Y22c} only | partial (only 3/5 of multiplet) |
| B12 | mathematica | 216 | `angular energy(basis[i]) - 6` for i in {Y20, Y21c, Y22c} only | partial (only 3/5) |
| B13 | mathematica | 222-225 | `-Delta_S2 basis[i] - 6 basis[i]` for i in {Y20, Y21c, Y22c} only | partial (only 3/5) |
| B14 | mathematica | 262-265 | single-component reduced density per basis member, three components | partial (only 3/5; SymPy does sum-over-5 in one shot, MA does 3 separately) |
| B15 | mathematica | 277 | single-component P2 Euler-Lagrange | yes |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl:74-187`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py:69-204`

**What's wrong:**
The Mathematica script is a line-by-line port of the SymPy script's algebra rather than an independent second-engine derivation. Direct evidence:

1. **Identical variable choreography in Section I.** Compare:
   - SymPy lines 74-78: `y00 = sp.Rational(1,2)/sp.sqrt(sp.pi)`; `avg_y00 = ... integrate(y00 * dOmega ...) / (4 pi)`.
   - WL lines 82-87: `y00 = 1/(2 Sqrt[Pi])`; `avgY00 = Integrate[y00 dOmega ...] / (4 Pi)`.
   The intermediate names (`y00`, `avg_y00`/`avgY00`, `mouthAvg`, `angularPrefactor`) and their order are identical.

2. **Identical Section II construction.** Both scripts define the two-mode ansatz with the same `2 sqrt(pi)` prefactor wrapping `Y_00`, the same `M_integrand` / `Mintegrand` with `4 pi` pulled out front by hand, and the same target check `lw - lw_target == 0`. Compare:
   - SymPy lines 137-154 build `M_integrand` and `K_integrand` as explicit `4*pi*[[...]]` 2x2 matrices.
   - WL lines 133-140 build `Mintegrand` and `Kintegrand` as the same `4 Pi {{...}}` 2x2 matrices, in the same order, with the same off-diagonal symmetrization.

3. **Identical EL choreography.** SymPy lines 185-204 build `lred_time` then call `euler_equations(...)`; WL lines 172-187 build `lredTime` then call the helper `eulerLagrange1D`. Both then compare against a hand-constructed `Maa qa'' + MaL qL'' + Kaa qa + KaL qL` expression with identical signs.

4. **Identical name strings for every assertion.** Every `expectZero[...]` in the WL has a matching `expect_zero(...)` in the SymPy with byte-identical name strings (e.g. `"Y00 mouth average - 1/(2 sqrt(pi))"`, `"reduced Lagrangian density from the action - target overlap form"`, `"Euler-Lagrange equation for q_a"`). Engines that derived the result independently would not converge to identical naming.

An independent Mathematica derivation would, for example, use `SphericalHarmonicY[0,0,theta,phi]` instead of hardcoding `1/(2 Sqrt[Pi])`, use `EulerEquations` from `VariationalMethods` instead of reproducing the helper `D[D[L, D[field,t]],t] - D[L,field]` formula, and arrive at the matrices via `Coefficient[]` extraction from the integrated Lagrangian rather than postulating them and checking they match. The current WL does none of this.

**Why this matters:**
The second-engine policy exists to catch errors that propagate consistently through one engine's algebraic conventions. A line-by-line transliteration cannot catch them; both engines execute the same algebraic plan, so they will both fail or both pass for the same reason. The PASS in both outputs is then evidence of internal consistency only, not of independent verification.

**Required change:**
Refactor the Mathematica script so its Section I and Section II derivations do not mirror the SymPy structure. Concretely:
- Use `SphericalHarmonicY[0, 0, theta, phi]` to obtain `Y00` (Mathematica returns `1/(2 Sqrt[Pi])`, so the numeric value is verified by the engine, not hardcoded).
- For Section II, do not predeclare `Mintegrand` and `Kintegrand` with `4 Pi` factored out. Instead, perform the full angular integral of the two-mode ansatz Lagrangian density and then `Coefficient[lw, dadt^2]`, `Coefficient[lw, dadt dLdt]`, `Coefficient[lw, dLdt^2]` to extract the components of `M`, and analogously for `K`. Then the check becomes whether the extracted matrix equals the boxed form, not whether two hand-built bilinear expansions agree.
- For the EL step, call `VariationalMethods` `EulerEquations[lredTime, {qa, qLfun}, t]` rather than reproducing the formula by hand.

These three changes break the line-by-line correspondence with the SymPy script without changing the physics tested. Lines affected: WL 82-104 (Section I), 110-187 (Section II).

**Verification:**
After the refactor, diff-style inspection should show no shared intermediate variable names with the SymPy script beyond the user-facing physical symbols (`mu_eta`, `T_w`, `K_0`, `alpha_a`, `alpha_L`, `delta_a`, `delta_L`, `q00`). The final assertions (the matrix-element values of M, K and the EL equation form) should still pass.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl:200-267`

**What's wrong:**
The Mathematica script's Section III claims to verify "real P2 orthonormality, common angular stiffness 6, and the resulting isotropic grouped-real P2 degeneracy before coupling" (script summary, WL lines 285-287), but actually exercises only 3 of the 5 P2 basis elements. Specifically:

- Line 200: `basis = {y20, y21c, y22c};` — `y21s` and `y22s` are defined (lines 197, 199) but omitted from the test basis.
- Lines 206-227: `Do[ ... , {i, 1, Length[basis]}]` only loops over the three-element basis; orthonormality, angular-energy = 6, and `-Delta_S2 y - 6 y = 0` checks are not performed on `y21s` and `y22s`.
- Lines 238-267: the single-component reduced-density verification likewise loops over only those three components.
- Lines 203-204 do verify `y21s = y21c(theta, phi - Pi/2)` and `y22s = y22c(theta, phi - Pi/4)` as substitutes. Phase-shift identities prove that one Y_2m is a rigid rotation of another in phi, which is fine for arguing the multiplet inherits the eigenvalue and norm — but the inheritance argument is not what the script claims; the script claims it verified the multiplet directly.

The SymPy script (lines 219-250) does test all five components in the orthonormality matrix (5x5) and stiffness matrix (5x5) and runs the eigenvalue check on each. So the two engines disagree on the breadth of the verification, which is also relevant under engine cross-check: the WL output's "isotropic 5-component degeneracy" PASS reflects 3 components, not 5.

**Why this matters:**
The off-diagonal entries of the 5x5 orthonormality and stiffness matrices (e.g. ⟨Y21c, Y22s⟩, ⟨Y20, Y21s⟩, the four entries pairing a `c` with an `s` member) are never integrated in Mathematica. The claim being asserted is that the entire 5x5 norm matrix equals `I_5` and the entire 5x5 stiffness matrix equals `6 I_5`. Using a 3x3 sub-block plus two phase-shift identities is a weaker test: it shows the diagonal of the 5x5 is correct and that the `s` rows are rotations of the `c` rows, but it does not directly verify that, e.g., `⟨Y20, Y21s⟩_{S^2}` vanishes via the angular integral the script purports to use.

**Required change:**
Extend the Mathematica `basis` and the loops to cover all five P2 components. Edit:
- WL line 200: change to `basis = {y20, y21c, y21s, y22c, y22s};` and `basisNames = {"Y20", "Y21c", "Y21s", "Y22c", "Y22s"};`.
- WL lines 206-219: keep the `Do[ ... ]` loop structure but also add an outer-product 5x5 norm matrix and a 5x5 angular-stiffness matrix, asserting they equal `IdentityMatrix[5]` and `6 IdentityMatrix[5]` respectively (parallel to SymPy lines 221-247 but written with `Mathematica`-idiomatic `Table[...]` and `IdentityMatrix[5]`, not as a literal port of the Python).
- WL lines 235-267: extend `qvec` to five symbols (`{q20, q21c, q21s, q22c, q22s}`), extend `qdotvec` similarly, and adjust the per-component loop to iterate over all five.
- Keep the phase-shift identities on lines 203-204 as supplementary checks; they are not the problem.

The two existing phase-shift assertions can stay; they just shouldn't be a substitute for the missing direct checks.

**Verification:**
The refreshed Mathematica output should show `norm(...)` and `angular energy(...)` PASS lines for `Y21s` and `Y22s` as well, and a 5x5 matrix identity assertion alongside (or in place of) the per-component checks. Total assertion count in Section III should increase from 14 to at least 16 (two new norm + two new angular energy) and ideally include explicit matrix-form `5x5 norm - I_5 == 0` and `5x5 stiffness - 6 I_5 == 0` lines.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a transliteration of the SymPy script's algebra. The evidence is enumerated in F1; the short version is that every meaningful step is reproduced in the same order with the same intermediate names and identical assertion strings, and no `SphericalHarmonicY`, `EulerEquations`, or `Coefficient`-based extraction is used despite all three being available in Mathematica's standard libraries. The one structural divergence (Section III uses only 3 components in Mathematica vs. 5 in SymPy) is itself a defect rather than evidence of independence.

## Engine cross-check

Both engines emit PASS for the assertions they perform. However, they do not test the same set of statements in Section III: SymPy verifies orthonormality and angular stiffness across the full 5x5 P2 matrix, while Mathematica only tests 3x3 (Y20, Y21c, Y22c). For the assertions that do overlap (Sections I and II, plus the diagonal of Section III's basis), the engines agree. The overall claim of "five-component isotropic degeneracy" is verified only by the SymPy engine; the Mathematica engine verifies a weaker three-component statement.

| Statement | SymPy | Mathematica |
|---|---|---|
| `q_00 = 2 sqrt(pi) delta_a` bridge | PASS | PASS |
| 4 pi angular prefactor | PASS | PASS |
| Two-mode `lw = lw_target` | PASS | PASS |
| Two-mode `lred = lred_target` | PASS | PASS |
| EL `(q_a, q_L)` matrix form | PASS | PASS |
| P2 orthonormality (full 5x5) | PASS | NOT TESTED (3x3 only) |
| P2 angular stiffness (full 5x5 = 6 I) | PASS | NOT TESTED (3x3 only) |
| `-Delta_S2 Y_2m = 6 Y_2m` for m in {0,±1,±2} | PASS (all 5) | PASS (3 of 5) |
| Five-component reduced density target | PASS (sum-of-5) | PASS (per-component, 3 of 5) |
| Single-component EL `M2 q'' + K2 q = 0` | PASS | PASS |

`engines_agree: false` in the front-matter because the engines diverge on whether the 5-component multiplet has been verified.

## Verdict justification

`findings` (not `clean`, not `stop_cold`). The algebra in both scripts holds up under direct inspection: the bilinear forms expand correctly, the `4 pi` angular factor is consistent with the `2 sqrt(pi) Y_00` ansatz, the Euler-Lagrange residuals are zero with the standard SymPy `∂L/∂q − d/dt(∂L/∂q̇)` sign convention, the P2 spherical harmonics in SymPy do satisfy `−Δ_{S^2} Y = 6 Y` and form an orthonormal set under direct integration, and the saved outputs are newer than the scripts (no `stale_output`). I attempted to attack the EL sign convention, the angular normalization, the use of `simplify` under `real=True` assumptions (no positivity claims hidden), and the spherical harmonic eigenvalue formulae; none broke. The two findings are at a higher structural level: the Mathematica script is a port rather than an independent re-derivation (F1), and that port quietly reduces the multiplet coverage from 5 to 3 (F2). Neither defect renders the unit's math wrong, but each violates the verification policy this stage is supposed to satisfy.
