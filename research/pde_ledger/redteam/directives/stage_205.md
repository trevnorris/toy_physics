---
unit_id: 205
batch: VI.1
created_at: 2026-06-01T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-02T11:09:56-06:00
findings_applied: 5
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 205

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

There are no `paper_misalignment` findings in this unit; no user resolution is pending.

---

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl` (NEW FILE — `.wl` files live in `mathematica/`).

**Issue:** Stage 205 is non-status-only and non-checkpoint but has only a SymPy engine. Its results are entirely CAS-verifiable. A transliteration-free Mathematica audit is required as the independent second engine.

**Requirement and acceptance criteria (Codex designs and writes the script):**
Create a new Mathematica audit that independently verifies the claim manifest below. Codex chooses the implementation route; do NOT port the SymPy algebra. The script must end in a clear PASS/`Exit[0]` on success and `Exit[1]` on any failed check, and each manifest item must have its own pass/fail check (use a strict zero-residual / boolean-equivalence pattern, e.g. an `expectZero`-style helper that strips `ConditionalExpression[0, …]` and treats a nonzero remainder as failure).

**Anti-transliteration guard (binding):** The `.wl` must derive each result independently using native Mathematica primitives — `Solve`/`Reduce`, `Series`+`Coefficient`, `D[]`, `Limit`, `Refine`/`Simplify` with `Assumptions`, and matrix operations (`Transpose`, `Dot`) — via a DIFFERENT decomposition than the `.py`. A line-by-line port of the SymPy variable choreography (same `Hchi/Phi0 - g g^T/Phi0^2` assembly order, same intermediate names) is rejected as transliteration. Concretely: build the log-Hessian from the actual second-derivative definition of `ln χ̂_Q` where practical, drive the predictors by `Solve`-ing the quadratic rather than substituting a pre-built closed form, and obtain the agreement coefficient by `Series` + `Coefficient` rather than dividing by `eps^3`.

**Claim manifest** (each must be an independent check):

- **M1 — log-curvature identity.** For a symmetric `5×5` matrix `H`, gradient vector `g`, direction `s`, and `Φ₀>0`, with `Φ₁ = sᵀg`, `Φ₂ = sᵀHs`, `L₀ = Φ₁/Φ₀`, and the log-Hessian `Hlog = H/Φ₀ − (g gᵀ)/Φ₀²`, `L₁ = sᵀ Hlog s`: verify
  `L₁ − (Φ₂/Φ₀ − Φ₁²/Φ₀²) = 0` (identically in all symbols).
- **M2 — ordinary/log bridge.** With the same definitions: `Φ₂ − Φ₀(L₁ + L₀²) = 0` (identically).
- **M3 — quadratic affine predictor residual, both slope signs.** With `Δ_aff = Φ₁² − 2Φ₂(Φ₀−1)` and `τ_quad = 2(1−Φ₀)/(Φ₁ + sgn(Φ₁)√Δ_aff)`: verify `(Φ₀−1) + Φ₁ τ_quad + ½ Φ₂ τ_quad² = 0` for `Φ₁>0` AND for `Φ₁<0` (treat `Φ₁` as nonzero real; do not assume positive). Obtain `τ_quad` by `Solve`-ing the quadratic and selecting the branch that reduces to `τ_aff` as `Φ₂→0`, not by hard-substituting the closed form.
- **M4 — affine zero-curvature limit.** `Limit[τ_quad, Φ₂ → 0] = (1−Φ₀)/Φ₁ = τ_aff` (on the selected branch, both slope signs).
- **M5 — quadratic log predictor residual, both slope signs.** With `Δ_log = L₀² − 2L₁ lnΦ₀` and `τ_log,2 = −2 lnΦ₀/(L₀ + sgn(L₀)√Δ_log)`: verify `lnΦ₀ + L₀ τ_log,2 + ½ L₁ τ_log,2² = 0` for `L₀>0` AND `L₀<0`.
- **M6 — log zero-curvature limit.** `Limit[τ_log,2, L₁ → 0] = −lnΦ₀/L₀ = τ_log` (both slope signs).
- **M7 — turning-point reality criterion (the headline theorem).** For real `Φ₀,Φ₂` with `Φ₁=0,Φ₀≠1`: verify the radicand `2(1−Φ₀)/Φ₂ ≥ 0` is logically equivalent (via `Reduce`) to `(1−Φ₀)Φ₂ ≥ 0`, i.e. real `τ_± = ±√(2(1−Φ₀)/Φ₂)` exist **iff** `(1−Φ₀)Φ₂ > 0`, and there is no real root when `(1−Φ₀)Φ₂ < 0`. Also verify each `τ_±` solves `(Φ₀−1) + ½Φ₂τ² = 0` only on the region where the criterion holds (the residual identity alone is by-construction and is NOT sufficient — the reality equivalence is the load-bearing check).
- **M8 — tangency model.** Substituting `Φ₀→1, Φ₁→0` into `(Φ₀−1) + Φ₁τ + ½Φ₂τ²` yields exactly `½Φ₂τ²`: verify `[(Φ₀−1)+Φ₁τ+½Φ₂τ²]|_{Φ₀=1,Φ₁=0} − ½Φ₂τ² = 0`.
- **M9 — ordinary curvature correction.** With `Φ₀=1+ε`, `Φ₁=Φ₀L₀`, `Φ₂=Φ₀(L₁+L₀²)`: the `Series` of `τ_quad − τ_aff + (Φ₂/(2Φ₁))τ_aff²` in `ε` about `0` to `O(ε²)` is `0` (i.e. `τ_quad = τ_aff − (Φ₂/2Φ₁)τ_aff² + O(τ_aff³)` holds).
- **M10 — logarithmic curvature correction.** Same setup: `Series` of `τ_log,2 − τ_log + (L₁/(2L₀))τ_log²` in `ε` to `O(ε²)` is `0`.
- **M11 — quadratic predictor agreement (cubic coefficient).** With `Φ₀=1+ε`: the `ε³` coefficient of `Series[τ_log,2 − τ_quad, {ε,0,3}]` equals `(L₀²+3L₁)/(6L₀³)`, and the `ε⁰,ε¹,ε²` coefficients are all `0` (so `τ_log,2 − τ_quad = O(ε³)`). Extract via `Coefficient[Normal[Series[…]], ε, k]`, not by dividing by `ε³`.

**Verification command:** verifier runs `redteam exec-mathematica 205`; the new `.wl` must contain checks M1–M11, exit 0, and pass the independent-derivation (non-transliteration) review.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl`
- summary: Created a Mathematica audit covering M1-M11 with derivative-built log Hessian checks, Solve-selected predictor branches, Reduce-based reality checks, tangency derivation, and series coefficient checks.
- deviation: none

---

## F2 — missing_branch

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py` — Section II (lines 84,88,97) and Section III (lines 104,107,116).

**Issue:** The affine and log predictors are tested only on the positive-slope branch because `Phi1a` (line 84) and `L0l` (line 104) are declared `positive=True`, so `sgn(·)` collapses to `+1`. The paper's predictors (appendix lines 710, 714; notes §5/§6) are general in the sign of the slope.

**Required change (Codex chooses exact code; the requirement is what must hold):**
Add negative-slope-branch coverage for both predictors so that the paper's `sgn`-general form is exercised. Acceptance criteria:
- A check that, with `Φ₁ < 0` (nonzero real, not positive) and the `sgn`-correct denominator `Φ₁ + sgn(Φ₁)√Δ_aff = Φ₁ − √Δ_aff`, the affine quadratic residual `(Φ₀−1) + Φ₁τ + ½Φ₂τ²` is zero, AND the `Φ₂→0` limit recovers `τ_aff = (1−Φ₀)/Φ₁` on that branch.
- The analogous check for the log predictor with `L₀ < 0`: denominator `L₀ − √Δ_log`, residual `lnΦ₀ + L₀τ + ½L₁τ²` zero, and `L₁→0` limit recovers `τ_log = −lnΦ₀/L₀`.
Either add a second block with negative-slope symbols, or generalize the existing symbols to nonzero-real and carry an explicit `sgn` factor; the positive-branch checks already present must remain.

**Verification command:** `redteam exec-sympy 205`; new negative-branch `expect_zero` lines appear and pass, script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py`
- summary: Added explicit negative-slope affine and logarithmic predictor residual and zero-curvature limit checks while preserving the existing positive-branch checks.
- deviation: none

---

## F3 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py:123–130`.

**Issue:** `tau_tp = sqrt(2*(1-Phi0t)/Phi2t)` then `res_tp_plus/minus = (Phi0t-1) + ½*Phi2t*tau_tp**2`. Since `tau_tp**2 = 2(1-Phi0t)/Phi2t` by construction, the residual is `0` identically regardless of physics or reality of the root. The actual theorem (paper §7.1 / appendix lines 717–721) is the reality criterion `(1−Φ₀)Φ₂>0`, which is untested.

**Required change (Codex chooses exact code; requirement stated):**
Test the turning-point **reality criterion**, not just the by-construction root identity. Acceptance criteria with `Phi0t,Phi2t` real:
- Establish (e.g. via `sp.simplify`/sign reasoning or comparing the radicand) that the radicand `2(1−Φ₀t)/Φ₂t ≥ 0` holds exactly on `(1−Φ₀t)Φ₂t ≥ 0`, i.e. the symmetric real predictors `τ_± = ±√(2(1−Φ₀t)/Φ₂t)` exist **iff** `(1−Φ₀t)Φ₂t > 0`.
- Confirm that under `(1−Φ₀t)Φ₂t < 0` the radicand is negative (no real closure point).
- A single non-tautological confirmation that `τ_±` solves `(Φ₀−1)+½Φ₂τ²=0` may remain, but it must not be the only/load-bearing check.
Remove or downgrade the current by-construction residuals at lines 129–130 so they are no longer the substantive verification of the theorem.

**Verification command:** `redteam exec-sympy 205`; a new assertion exercising the `(1−Φ₀)Φ₂>0 ⟺ real root` equivalence appears and passes.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py`
- summary: Replaced the load-bearing turning-point residual checks with sign-bridge and positive/negative product radicand checks, leaving one residual check as secondary confirmation on the real-root branch.
- deviation: none

---

## F4 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py:132–135`.

**Issue:** The tangency model (paper §7.2: `Δ_s(τ)=½Φ₂τ²+O(τ³)` at `Φ₀=1,Φ₁=0`) is only printed (`Delta_tangent = ½*Phi2g*tau**2`); the hard-written form is not derived from the quadratic model, so nothing is asserted.

**Required change:**
Add an `expect_zero` that derives the tangency form: substitute `Φ₀→1, Φ₁→0` into the ordinary quadratic model `(Φ₀−1) + Φ₁τ + ½Φ₂τ²` and assert the result minus `½Φ₂τ²` is zero. Use the curvature symbol consistently (the substituted model's `Φ₂` should match the printed `Phi2g`, or unify the symbol). Keep the existing print if desired, but it must be backed by this assertion.

**Verification command:** `redteam exec-sympy 205`; new `expect_zero("tangency model …", …)` appears and passes.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py`
- summary: Added an `expect_zero` deriving the tangency form by substituting `Phi0=1` and `Phi1=0` into the ordinary quadratic closure model.
- deviation: none

---

## F5 — banner mislabel (low; mechanical, apply alongside)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py:35,174`.

**Issue:** The banners read `STAGE 188 — …` and `STAGE 188 SYMPY AUDIT PASSED`, but this is stage 205. The wrong stage label in the printed transcript is misleading. (Not a math finding; mechanical correction only — change the printed banner text to reference stage 205. Do not alter any computation.)

**Required change:** Update the two banner string literals to read `STAGE 205 …` / `STAGE 205 SYMPY AUDIT PASSED`. If the new `.wl` (F1) prints a banner, label it `STAGE 205` as well.

**Verification command:** `redteam exec-sympy 205`; transcript header reads STAGE 205.

## Applied: F5

- files_changed:
  - `scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl`
- summary: Updated the SymPy banner strings to stage 205 and labeled the new Mathematica audit banner as stage 205.
- deviation: none
