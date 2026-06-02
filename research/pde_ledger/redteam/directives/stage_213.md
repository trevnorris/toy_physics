---
unit_id: 213
batch: VI.1
created_at: 2026-06-02T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T16:19:31Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 213

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the new script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl`

(`.wl` files live in `mathematica/`. The filename mirrors the SymPy sibling with `_sympy_audit.py` → `_mathematica_audit.wl`, matching the existing `...stage218..._mathematica_audit.wl` convention.)

**Issue:** Stage 213 is non-status-only and non-checkpoint but has only a SymPy engine; the paper card itself notes "Mathematica audit: none yet." Every load-bearing result of this stage is a closed-form symbolic identity that Mathematica can independently verify, so the dual-engine contract requires a `.wl`. This directive states *what must hold* (the claim manifest below) and the acceptance criteria; it does NOT prescribe the implementation. Codex designs and writes the script.

**Required change:** Author a new, self-contained Mathematica audit script at the Target path that independently verifies the claim manifest M1–M9 below. The script must print a labeled residual (or boolean) for each claim and must `Exit[1]` on any nonzero residual / False (so a clean run exits 0). Use a helper such as `expectZero[name_, expr_]` that `FullSimplify`s the expression, prints it, and `Exit[1]`s if it is not identically zero; preemptively strip `ConditionalExpression[0, ...]` from any `Solve`/`Reduce`/`Maximize` output before the zero test.

**Anti-transliteration guard (mandatory):** The `.wl` must derive each result independently using native Mathematica primitives — `Solve`/`Reduce`/`Maximize`, `Series`+`Coefficient`, `D[]`, `Integrate`, `FindRoot`, matrix ops — via a DIFFERENT decomposition than the `.py`. A line-by-line port of the SymPy algebra (same variable choreography, same intermediate `subs`/`expand` steps rewritten in Wolfram syntax) is rejected as transliteration. In particular: derive the gradient optimum and the cross-leverage bound by constrained optimization (`Maximize`/Lagrange) rather than re-asserting the closed forms the `.py` posits; obtain the discriminant numerator by `Expand` + `Coefficient` extraction from `(1+r^2+s^2+t^2) (kSimplex^2 - 2 H0 kappaSimplex)` rather than re-typing the ten coefficients A…J by hand.

**Claim manifest** (each must be independently verified; symbols: `a_i,a_j,a_k,a_l ≥ 0`, `k_i,k_j,k_k,k_l > 0`, `H0 > 0`, `r,s,t,u,v ≥ 0`, symmetric `4×4` Hessian entries `u_ii … u_ll` real):

- **M1 — combinatorial ledger.** With axes `{λ,c,γ,U,W}`: the number of 4-element subsets is `Binomial[5,4] = 5`; every 3-element subset is contained in exactly 2 of those quadruples; every single axis appears in exactly 4. (Verify by enumeration over `Subsets`, asserting the three counts.)

- **M2 — face slope reductions.** For the simplex slope `k_simplex(a) = a_i k_i + a_j k_j + a_k k_k + a_l k_l` evaluated on each unit-norm face parametrization
  - `a_ijk(r,s) = (1,r,s,0)/√(1+r²+s²)` → `(k_i + r k_j + s k_k)/√(1+r²+s²)`,
  - `a_ijl(r,t) = (1,r,0,t)/√(1+r²+t²)` → `(k_i + r k_j + t k_l)/√(1+r²+t²)`,
  - `a_ikl(s,t) = (1,0,s,t)/√(1+s²+t²)` → `(k_i + s k_k + t k_l)/√(1+s²+t²)`,
  - `a_jkl(u,v) = (0,1,u,v)/√(1+u²+v²)` → `(k_j + u k_k + v k_l)/√(1+u²+v²)`,
  each with the face vector unit-normalized. (Residual = evaluated slope minus stated RHS = 0; and `aᵀa − 1 = 0` per face.)

- **M3 — gradient-optimal ray (by optimization, not by assertion).** Maximizing `k_simplex(a)` subject to `a_i²+a_j²+a_k²+a_l² = 1` yields the maximizer `a_grad = (k_i,k_j,k_k,k_l)/√(k_i²+k_j²+k_k²+k_l²)` with maximum value `√(k_i²+k_j²+k_k²+k_l²)`. Verify the maximum value equals `‖k‖` via `Maximize`/Lagrange (independent route), and that the optimal interior ratios are `a_j/a_i = k_j/k_i`, `a_k/a_i = k_k/k_i`, `a_l/a_i = k_l/k_i`.

- **M4 — synergy gaps.** `‖k‖² − ‖k_ijk‖² = k_l²`, `‖k‖² − ‖k_ijl‖² = k_k²`, `‖k‖² − ‖k_ikl‖² = k_j²`, `‖k‖² − ‖k_jkl‖² = k_i²`, where `‖k_ijk‖² = k_i²+k_j²+k_k²` etc. (Four zero residuals.)

- **M5 — cross-leverage identity and bound.** `w_Σ(a) := 2(a_i a_j + a_i a_k + a_i a_l + a_j a_k + a_j a_l + a_k a_l)` equals `(a_i+a_j+a_k+a_l)² − (a_i²+a_j²+a_k²+a_l²)` identically; and the Cauchy slack identity `4(a_i²+a_j²+a_k²+a_l²) − (Σa)² = Σ_{p<q}(a_p−a_q)²` holds. Additionally, the constrained maximum of `w_Σ` over `a_i²+…+a_l²=1, a≥0` is `3`, attained at `a=(1/2,1/2,1/2,1/2)` (verify the bound `≤ 3` by `Maximize`, not by substitution alone), and `w_Σ = 2` at the unit-norm triple-equal-mix point `(1,1,1,0)/√3`, `w_Σ = 1` at the pair-equal-edge `(1,1,0,0)/√2`.

- **M6 — curvature quadratic form and diagonal-neutral reduction.** `κ(a) := aᵀ H a` for the symmetric `4×4` block `H` with entries `u_••`; when all six off-diagonal entries vanish, `κ` reduces to `u_ii a_i² + u_jj a_j² + u_kk a_k² + u_ll a_l²`. (Zero residual after `Expand`.)

- **M7 — discriminant numerator.** On the ratio patch `a(r,s,t) = (1,r,s,t)/√(1+r²+s²+t²)`, with the ten coefficients
  `A=k_i²−2H₀u_ii, B=2k_ik_j−4H₀u_ij, C=2k_ik_k−4H₀u_ik, D=2k_ik_l−4H₀u_il, E=k_j²−2H₀u_jj, F=2k_jk_k−4H₀u_jk, G=2k_jk_l−4H₀u_jl, Ĥ=k_k²−2H₀u_kk, I=2k_kk_l−4H₀u_kl, J=k_l²−2H₀u_ll`,
  the identity `(1+r²+s²+t²)(k_simplex² − 2H₀κ) = A + Br + Cs + Dt + Er² + Frs + Grt + Ĥs² + Ist + Jt²` holds. Codex must obtain the RHS polynomial by `Expand`+`Coefficient` from the LHS (independent extraction), then confirm it matches the boxed `Δ♯` coefficient set above — NOT by typing `Δ♯` and subtracting.

- **M8 — certified τ bracket and ratio form.** `τ(a) := 2H₀/(k_simplex + √(k_simplex² − 2H₀κ))` evaluated on `a(r,s,t)` equals `2H₀√(1+r²+s²+t²)/(k_i + r k_j + s k_k + t k_l + √Δ♯)`. (Zero residual.)

- **M9 — face collapses of the bracket.** The ratio-form τ reduces correctly on each face:
  - `t=0` → `2H₀√(1+r²+s²)/(k_i+r k_j+s k_k + √(A+Br+Cs+Er²+Frs+Ĥs²))`,
  - `s=0` → `2H₀√(1+r²+t²)/(k_i+r k_j+t k_l + √(A+Br+Dt+Er²+Grt+Jt²))`,
  - `r=0` → `2H₀√(1+s²+t²)/(k_i+s k_k+t k_l + √(A+Cs+Dt+Ĥs²+Ist+Jt²))`,
  - and the jkl face `a_jkl(u,v)` gives `2H₀√(1+u²+v²)/(k_j+u k_k+v k_l + √(E+Fu+Gv+Ĥu²+Iuv+Jv²))`.
  (Four zero residuals.)

Optionally (mirrors `.py` Section VI, but the interval-gate checks are trivial ordering relations and need not be reproduced; M1–M9 are the load-bearing symbolic content).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 213` and confirm the new checks M1–M9 appear AND the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl`
- summary: Added the missing Mathematica audit script verifying M1-M9 with independent combinatorial, optimization, matrix, coefficient-extraction, and tau face-collapse checks.
- deviation: none
