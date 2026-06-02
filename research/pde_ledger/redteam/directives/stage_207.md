---
unit_id: 207
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-02T10:05:52-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 207

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond authoring the one new file named below.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts. You may NOT edit the existing SymPy script (it is sound and aligned); this directive only adds the missing Mathematica engine.

## F1 — missing_verification_script

**Subtype:** missing_mathematica

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_mathematica_audit.wl`
(`.wl` files live in `mathematica/`. Match the existing naming convention `moving_throat_pde_stage<NNN>_<slug>_mathematica_audit.wl`; the slug is `primitive_ray_hessian_envelopes_and_certified_ray_table`.)

**Issue:** Stage 207 is non-status-only and non-checkpoint and computes multiple independently-verifiable closed-form identities, but only a SymPy engine exists; the paper card itself records "Mathematica audit: none yet." Per the dual-engine contract a second engine is required wherever Mathematica *can* independently verify, which it can here. Author a new `.wl` that re-derives and verifies the claim manifest M1–M6 below.

**Required change:**
Create the target `.wl`. It must contain real in-file checks (each manifest item must `Print` a residual that is `0`/`True` and trigger `Exit[1]` on failure — e.g. an `expectZero`-style helper, with the [[feedback-mathematica-script-idioms]] guard: strip `ConditionalExpression[0, …]` from any `Solve`/`Reduce`/`Sign` output before the zero test, and check poles via `1/expr == 0` rather than `=!= Infinity`). It must verify every manifest item M1–M6. Use symbols with the same domains the physics requires: `H0, k, k_i > 0`; `\(\kappa_{\rm lo},\kappa_{\rm hi}\)` real (turning-row curvatures negative); `chi0_star, deltaU_star, E_star, F_star > 0`; the Hessian entries real.

**Anti-transliteration guard:** The `.wl` must derive each result *independently* using native Mathematica primitives (`Solve`/`Reduce`, `Series`+`Coefficient`, `D[]`, `Integrate`, `FindRoot`, matrix ops) via a DIFFERENT decomposition than the `.py`. In particular, for the certified brackets (M4) the SymPy script *constructs* the closed-form root `2H0/(k+Sqrt[k^2-2 c H0])` and substitutes it back into the comparison quadratic; the `.wl` must instead `Solve`/`Reduce` the comparison quadratic `H0 - k*tau + (1/2)*c*tau^2 == 0` for `tau`, SELECT the physically-correct (smaller, `tau>0`) root, and confirm it equals the paper's closed form — not echo the construct-then-substitute choreography. A line-by-line port of the SymPy algebra (same variable choreography, same intermediate steps rewritten in Wolfram syntax) is rejected as transliteration.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1 — diagonal Hessian reduction.** For a symbolic real symmetric 5x5 Hessian `H` over axes `\(\{\lambda,c,\gamma,U,W\}\)` with unit basis vectors `\(\mathbf e_i\)`: for every `i`, `\(\mathbf e_i^{\mathsf T} H\,\mathbf e_i = H_{ii}\)` and `\((-\mathbf e_i)^{\mathsf T} H\,(-\mathbf e_i) = H_{ii}\)`. (Verify all five diagonal entries `\(H_{\lambda\lambda},H_{cc},H_{\gamma\gamma},H_{UU},H_{WW}\)`; the residual `\(\mathbf e_i^{\mathsf T} H\,\mathbf e_i - H_{ii}\)` is identically `0`.)

- **M2 — mixed-ray quadratic form.** For a genuine two-coordinate ray `\(\mathbf s=a\,\mathbf e_i+b\,\mathbf e_j\)` (`\(i\neq j\)`): `\(\mathbf s^{\mathsf T} H\,\mathbf s = a^2 H_{ii} + 2ab\,H_{ij} + b^2 H_{jj}\)`. Extract the cross coefficient via `Coefficient[..., a b]` and confirm it equals `\(2 H_{ij}\)` (factor of 2 is load-bearing — this is where off-diagonal data first enters).

- **M3 — canonical orientation law.** With `\(\varepsilon:=-\operatorname{Sign}(\Gamma)\)` and oriented slope `\(K:=\varepsilon\Gamma\)` for `\(\Gamma\neq0\)`: `\(K + |\Gamma| = 0\)`, i.e. the oriented primitive slope is `\(-|\Gamma|<0\)`. (Use `Sign`/`Abs`; strip any `ConditionalExpression` before the zero test.)

- **M4 — certified monotone bracket root map.** The Stage 206/240 comparison quadratic on an oriented primitive row is `\(H_0 - k\,\tau + \tfrac12 c\,\tau^2 = 0\)` with `\(k=|\Gamma_i|>0\)`. Solve/Reduce for `\(\tau\)`, select the smaller positive root, and confirm it equals `\(\mathcal T(H_0,k;c)=\dfrac{2H_0}{k+\sqrt{k^2-2cH_0}}\)`. Verify for BOTH envelope labels: the lower `\(c=\kappa_{\rm lo}\)` giving `\(\tau_{\rm lo}\)` and the upper `\(c=\kappa_{\rm hi}\)` giving `\(\tau_{\rm hi}\)`. Equivalently/additionally confirm each closed form satisfies the quadratic (residual `0`). Convention check: the slope enters the quadratic as `\(-k\)` (negative, matching `\(K_i=-|\Gamma_i|<0\)`), and the `\(+k\)` form of `\(\mathcal T\)` matches the notes (lines 240, 261-263); do not flip the sign.

- **M5 — certified turning bracket.** For a turning row (`\(\Gamma_i=0\)`, curvature `\(\kappa<0\)`), the turning comparison is `\(H_0 - \tfrac12\,\kappa\,\tau^2 = 0\)` (here `\(\kappa\)` is `a_turn`/`b_turn` in the SymPy script, kept positive there as `\(-\)curvature`). Solve for `\(\tau>0\)` and confirm `\(\tau^{\rm(tp)}=\sqrt{2H_0/\kappa_{\rm pos}}\)` (equivalently `\(\sqrt{-2H_0/\kappa}\)` with `\(\kappa<0\)`), for both `\(\tau_{\rm lo}^{\rm(tp)}\)` and `\(\tau_{\rm hi}^{\rm(tp)}\)`.

- **M6 — sign-adapted primitive drift table (paper Section 4 / appendix Section 4).** Define the generic dependent-exponent functions of the free log-ray direction `\(\mathbf s=(s_\lambda,s_c,s_\gamma,s_U,s_W)\)`, with `\(\mathfrak a_*:=\dfrac{1+\delta_{U,*}}{1+\chi_{0,*}}\)`:
  - `\(\sigma_\delta = -\mathfrak a_*(s_\gamma+s_c-s_U)\)`
  - `\(\sigma_T = s_U + \sigma_\delta\)`
  - `\(\sigma_{K_\eta} = 2s_c - s_U\)`
  - `\(\sigma_\mu = 2s_c - s_U + 2s_W - 2s_\lambda - E_*(2s_\gamma+2s_\lambda-s_U-s_W) + F_*\sigma_\delta\)`

  Evaluate `\((\sigma_\delta,\sigma_T,\sigma_{K_\eta},\sigma_\mu)\)` on each oriented primitive vector `\(\mathbf s=\varepsilon_i\,\mathbf e_i\)` and confirm it equals `\(\varepsilon_i\)` times the tabulated row:
  - `\(\lambda\)`: `\((0,\,0,\,0,\,\varepsilon_\lambda(-2-2E_*))\)`
  - `\(c\)`: `\((-\varepsilon_c\mathfrak a_*,\,-\varepsilon_c\mathfrak a_*,\,2\varepsilon_c,\,\varepsilon_c(2-F_*\mathfrak a_*))\)`
  - `\(\gamma\)`: `\((-\varepsilon_\gamma\mathfrak a_*,\,-\varepsilon_\gamma\mathfrak a_*,\,0,\,\varepsilon_\gamma(-2E_*-F_*\mathfrak a_*))\)`
  - `\(U\)`: `\((+\varepsilon_U\mathfrak a_*,\,\varepsilon_U(1+\mathfrak a_*),\,-\varepsilon_U,\,\varepsilon_U(-1+E_*+F_*\mathfrak a_*))\)`
  - `\(W\)`: `\((0,\,0,\,0,\,\varepsilon_W(2+E_*))\)`

  (Each residual `row(\(\varepsilon_i\mathbf e_i\)) - \(\varepsilon_i\)·expected` is identically `0`. Use a `FullSimplify` over `\(\chi_{0,*}+1\neq0\)` so the `\(\mathfrak a_*\)` denominator clears.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 207` and confirm the new checks (M1–M6) appear AND the script exits 0; the Mathematica residuals must agree with the SymPy output for the shared identities (M1, M2, M3, M4, M5, M6).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_mathematica_audit.wl`
- summary: Added the missing Mathematica audit verifying M1-M6 with independent matrix, coefficient, sign, Solve-selected root, turning-root, and primitive drift-table checks.
- deviation: none
