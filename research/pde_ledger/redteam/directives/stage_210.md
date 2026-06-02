---
unit_id: 210
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T10:14:55-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 210

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose documents — the red-team only modifies scripts.

After creating the new script, RUN it (`math -script <path>`, wrapped in `timeout 600`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F1 — missing_verification_script

**Subtype:** missing_mathematica

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl`

**Issue:** Stage 210 is non-status-only and non-checkpoint and performs fully verifiable symbolic algebra (boundary reductions, a Lagrange/Cauchy-optimal point structure, a matrix quadratic-form curvature scalar, a square-root root map, and a polynomial discriminant that clears the square root), yet it has no Mathematica `.wl`. Mathematica can independently verify all of it, so the project dual-engine contract requires a second-engine script. The existing passing SymPy script does not discharge this requirement. The paper card itself records "Mathematica audit: none yet."

**Required change:** Create the `.wl` at the Target path. It must independently verify the claim manifest M1–M9 below. State each result symbolically; the script must establish each one and abort with a nonzero exit (e.g. `Exit[1]`) if any residual is nonzero. Symbol domains: `k_i,k_j,k_k > 0`; `a_i,a_j,a_k ≥ 0`; `H_0 > 0`; `r,s,nu ≥ 0`; the six Hessian-envelope entries `u_ii,u_ij,u_ik,u_jj,u_jk,u_kk` real. Use `Assuming`/`$Assumptions` to carry these, exactly as the physical setup demands.

**Anti-transliteration guard:** This must be an INDEPENDENT re-derivation, not a port of the `.py`. Derive each result using native Mathematica primitives — `Solve`/`Reduce`, `Series` + `Coefficient`, `D[]`, `Integrate`, `FindRoot`, `Eigenvalues`/matrix ops, `CoefficientList` — via a DIFFERENT decomposition than the SymPy script. A line-by-line port of the SymPy `subs`/`simplify` choreography (same intermediate variables, same substitution order) is rejected as `mathematica_transliteration`. Concretely, the SymPy script verifies the discriminant identity (M7) by substituting the ratio parametrization, clearing `(1+r²+s²)`, and `simplify`-ing to zero; the `.wl` should instead obtain `Δ^♯` by `CoefficientList`/`Series`+`Coefficient` extraction of the cleared numerator in the monomials `{1, r, s, r², rs, s²}` and check the extracted coefficients equal `A,B,C,D,E,F` term-by-term. Similarly, verify M3/M4 (gradient-optimal point) via Lagrange stationarity `Solve` on `maximize a·k subject to a·a==1`, not by asserting the closed form; and obtain the edge reductions (M2, M8) via `Limit`/substitution along a DIFFERENT edge parametrization than the `.py` uses where practical.

**Claim manifest** (each must be independently verified):

- **M1** — Edge normalization: with `a_ij(r)=(1,r,0)/√(1+r²)`, `a_ik(s)=(1,0,s)/√(1+s²)`, `a_jk(ν)=(0,1,ν)/√(1+ν²)`, each satisfies `a·a = 1`.
- **M2** — Boundary slope reduction: `(a·k)|_{a=a_ij(r)} = (k_i + r k_j)/√(1+r²)`, and the analogous identities for `a_ik(s)` → `(k_i + s k_k)/√(1+s²)` and `a_jk(ν)` → `(k_j + ν k_k)/√(1+ν²)`. (Each boundary edge reduces to a pairwise cone.)
- **M3** — Gradient-optimal point: the maximizer of `a·k` on the unit sphere `a·a=1` (a_i,a_j,a_k ≥ 0) is `a^grad = (k_i,k_j,k_k)/√(k_i²+k_j²+k_k²)`; verify it satisfies `a^grad·a^grad = 1` AND arises as the Lagrange stationary point (do NOT just substitute the closed form).
- **M4** — Maximum slope value: `a^grad·k = √(k_i²+k_j²+k_k²)`; and the interior optimal ratios are `a^grad_j/a^grad_i = k_j/k_i`, `a^grad_k/a^grad_i = k_k/k_i`.
- **M5** — Strict interior gain (Pythagorean decomposition): `(k_i²+k_j²+k_k²) − (k_i²+k_j²) − k_k² = 0`, and the two cyclic analogues `… −(k_i²+k_k²) − k_j² = 0`, `… −(k_j²+k_k²) − k_i² = 0` (these encode that the triple max strictly exceeds each pairwise edge max since each `k>0`).
- **M6** — Cross-leverage identity and bound: `w_Σ := 2(a_i a_j + a_i a_k + a_j a_k) = (a_i+a_j+a_k)² − (a_i²+a_j²+a_k²)`; the Cauchy slack identity `3(a_i²+a_j²+a_k²) − (a_i+a_j+a_k)² = (a_i−a_j)² + (a_i−a_k)² + (a_j−a_k)²`; and the screen values `w_Σ((1,1,1)/√3) = 2` (equal-mix barycenter, the unique maximizer) and `w_Σ((1,1,0)/√2) = 1` (pairwise equal edge).
- **M7** — Curvature law and discriminant numerator (load-bearing): with `H` the symmetric 3×3 envelope matrix `[[u_ii,u_ij,u_ik],[u_ij,u_jj,u_jk],[u_ik,u_jk,u_kk]]` and the ratio parametrization `(a_i,a_j,a_k)=(1,r,s)/√(1+r²+s²)`, the curvature scalar `κ=aᵀHa` and slope `k_simplex=a·k` satisfy `(1+r²+s²)·(k_simplex² − 2 H_0 κ) = A + B r + C s + D r² + E r s + F s²`, where `A=k_i²−2H_0 u_ii`, `B=2k_i k_j−4H_0 u_ij`, `C=2k_i k_k−4H_0 u_ik`, `D=k_j²−2H_0 u_jj`, `E=2k_j k_k−4H_0 u_jk`, `F=k_k²−2H_0 u_kk`. Verify by extracting the six monomial coefficients of the left side (in `{1,r,s,r²,rs,s²}`) and matching `A..F` individually.
- **M8** — Curvature edge reduction and diagonal-neutral case: `κ|_{a=a_ij(r)} = (u_ii + 2 u_ij r + u_jj r²)/(1+r²)` and the two cyclic analogues for `a_ik(s)`, `a_jk(ν)`; and `κ|_{u_ij=u_ik=u_jk=0} = u_ii a_i² + u_jj a_j² + u_kk a_k²`.
- **M9** — Fixed-simplex root map and its boundary reductions: with `τ(a) = 2H_0/(k_simplex + √(k_simplex² − 2H_0 κ))`, the ratio-coordinate form is `τ = 2H_0 √(1+r²+s²)/(k_i + r k_j + s k_k + √(Δ^♯))`; and `τ` restricted to `s=0` reduces to the pairwise-ij form `2H_0 √(1+r²)/(k_i + r k_j + √(k_i²−2H_0 u_ii + (2k_i k_j−4H_0 u_ij) r + (k_j²−2H_0 u_jj) r²))`, with the analogous reductions on `r=0` (pairwise ik) and along the `a_jk(ν)` edge (pairwise jk).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 210` (or `math -script` on the Target path) and confirms the new checks M1–M9 appear AND the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage210_three_coordinate_mixed_simplex_audit_mathematica_audit.wl`
- summary: Created the Mathematica audit script verifying manifest checks M1-M9 with symbolic residual checks and nonzero failure exits.
- deviation: none
