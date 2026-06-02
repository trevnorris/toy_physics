---
unit_id: 208
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T10:04:30-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 208

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named here.

After creating the script, RUN it (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose document. The red-team only modifies scripts.

## F1 — missing_verification_script

**Subtype:** missing_mathematica

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_mathematica_audit.wl`
(`.wl` files live in `mathematica/`. Create this file; do not place it in `scripts/`.)

**Issue:**
Stage 208 is non-status-only and non-checkpoint, so the dual-engine contract requires both a SymPy and a Mathematica audit; only the SymPy script exists. All of this stage's results are elementary symbolic linear algebra and single-variable calculus, which Mathematica can verify independently. Add a Mathematica audit that re-derives and verifies the claim manifest below.

**Required change:**
Create the `.wl` at the Target path. It must independently establish each manifest claim M1–M9 and assert each is zero (or the stated identity holds) using an `expectZero`-style helper that strips `ConditionalExpression[0, ...]` and `Exit[1]`s on any nonzero residual, printing a clear PASS banner at the end. Do not hardcode the answers — derive them in-script from the same physical premises the SymPy script uses (the mixed-ray direction, the gradient vector, the symmetric Hessian, the entrywise envelopes, and the carried root map), but via a DIFFERENT decomposition (see anti-transliteration guard).

Symbol setup (state, do not let it drift from the paper): real symbols `ki, kj` with `ki>0, kj>0`; real `r` with `r>=0`; real `hii, hij, hjj` (off-diagonal Hessian entry symmetric); real envelope symbols `hiiLo, hijLo, hjjLo, hiiHi, hijHi, hjjHi`; positive real `H0`. Mixed-ray direction `s(r) = {1, r}/Sqrt[1+r^2]`; gradient `g = {-ki, -kj}`; pair Hessian `Hpair = {{hii, hij},{hij, hjj}}`.

**Claim manifest** (each M must be independently verified by the new `.wl`):

- **M1 — oriented slope law.** `K_ij(r) := g . s(r)` simplifies to `-(ki + r*kj)/Sqrt[1+r^2]`, and the positive slope magnitude `k_ij(r) := -K_ij(r) = (ki + r*kj)/Sqrt[1+r^2]`. (paper: notes §2; appendix eq. `app-part06-pairwise-slope`.)
- **M2 — slope derivative law.** `D[k_ij(r), r]` simplifies to `(kj - ki*r)/(1+r^2)^(3/2)`. (notes §2.)
- **M3 — gradient-optimal ratio and value.** Solving `D[k_ij(r), r] == 0` for `r` over `r>=0` yields the unique root `r_grad = kj/ki`, and `k_ij(r_grad)` simplifies to `Sqrt[ki^2 + kj^2]` (equivalently `k_ij(r_grad)^2 - ki^2 - kj^2 == 0`). Derive `r_grad` via `Solve`/`Reduce` on the stationarity equation, NOT by substituting a pre-baked `kj/ki`. (notes §2, §5.1; appendix eq. `app-part06-pairwise-gradient-ratio`.)
- **M4 — mixed curvature decomposition and diagonal neutrality.** `H1(r) := s(r) . Hpair . s(r)` simplifies to `(hii + 2*r*hij + r^2*hjj)/(1+r^2)`; and setting `hij -> 0` gives `(hii + r^2*hjj)/(1+r^2)`. Obtain the cross weight as the `hij`-coefficient of `Expand[H1(r)*(1+r^2)]` divided by `(1+r^2)` (i.e., via `Coefficient[...]`), not by re-stating `2r/(1+r^2)`. (notes §3; appendix eq. `app-part06-pairwise-curvature`.)
- **M5 — cross-weight extremum.** `w_x(r) := 2*r/(1+r^2)` (as recovered in M4) satisfies `D[w_x(r), r] == 2*(1-r^2)/(1+r^2)^2`, `D[w_x(r), r] /. r->1` is `0`, and `w_x(1) == 1`. (notes §3; appendix eq. `app-part06-pairwise-cross-weight`.)
- **M6 — canonical screen rays.** The gradient-optimal direction `s(r_grad)` equals `{ki, kj}/Sqrt[ki^2+kj^2]`; its slope `g . ({ki,kj}/Sqrt[ki^2+kj^2])` equals `-Sqrt[ki^2+kj^2]`; its curvature (Rayleigh quotient) equals `(ki^2*hii + 2*ki*kj*hij + kj^2*hjj)/(ki^2+kj^2)`. The equal-mix direction `s(1) = {1,1}/Sqrt[2]`; its slope equals `-(ki+kj)/Sqrt[2]`; its curvature equals `(hii + 2*hij + hjj)/2`. (notes §5.1–5.2; appendix eq. `app-part06-pairwise-canonical-screens`.)
- **M7 — coincidence condition.** The two canonical rays coincide iff `ki == kj`; verify `r_grad - 1 == (kj - ki)/ki` and that `s(r_grad) - s(1)` simplifies to the zero vector under the substitution `kj -> ki`. (notes §5.3.)
- **M8 — mixed curvature envelopes.** `kappa_lo(r) := (hiiLo + 2*r*hijLo + r^2*hjjLo)/(1+r^2)` and `kappa_hi(r) := (hiiHi + 2*r*hijHi + r^2*hjjHi)/(1+r^2)` each equal the weight-form `w_i*h..Lo + w_x*h..ij + w_j*h..jj` (with `w_i = 1/(1+r^2)`, `w_j = r^2/(1+r^2)`). (notes §4.)
- **M9 — certified bracket root relation.** With the carried root map `T(H0,k,c) := 2*H0/(k + Sqrt[k^2 - 2*c*H0])`, set `tau_lo := T(H0, k_ij(r), kappa_lo(r))` and `tau_hi := T(H0, k_ij(r), kappa_hi(r))`. Verify each satisfies the closure quadratic `H0 - k_ij(r)*tau + (1/2)*kappa_*(r)*tau^2 == 0`. Confirm this by `Simplify`-ing the residual to 0 (or by `Reduce`-ing the quadratic and confirming `T` is the smaller root), NOT by algebraically re-deriving `T` from the same quadratic and pattern-matching. (notes §4; appendix eqs. `app-part06-pairwise-root-function`.)

**Anti-transliteration guard:**
The `.wl` must derive each claim independently using native Mathematica primitives — `Solve`/`Reduce`, `Series`+`Coefficient`, `D[]`, `Integrate`, `FindRoot`, matrix/dot operations, `Simplify`/`FullSimplify` — via a DIFFERENT decomposition than the SymPy script. In particular: recover the cross weight from `Coefficient[...]` of the expanded quadratic form (M4) rather than re-stating `2r/(1+r^2)`; obtain `r_grad` by `Solve`-ing the stationarity equation (M3) rather than substituting `kj/ki`; verify the bracket by reducing the closure quadratic (M9) rather than re-substituting a closed-form `T`. A line-by-line port of the `.py` algebra (same variable choreography, same intermediate steps rewritten in Mathematica syntax) is rejected as transliteration.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 208` and confirm the new `.wl` appears at the Target path, contains substantive checks for M1–M9, and exits 0 with all checks passing.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_mathematica_audit.wl`
- summary: Created the Mathematica audit for claims M1-M9 using native matrix, Solve/Reduce, coefficient extraction, and closure-residual checks.
- deviation: none
