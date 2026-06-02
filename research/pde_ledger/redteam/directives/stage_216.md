---
unit_id: 216
batch: VI.1
created_at: 2026-06-01T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T11:23:26-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 216

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

The shared claim manifest M1-M6 below is referenced by BOTH F1 and F2. It is reconstructed from the paper card (`paper/stages/stage_216.tex`), the source notes (`notes/stages/moving_throat_pde_stage216_*.md`), and the appendix narrative (`paper/appendices/stage_appendix_part06.tex`, lines 1093-1154), cross-checked against the existing SymPy computations.

## Claim manifest (shared by F1 and F2)

Let `\(k_\lambda,k_c,k_\gamma,k_U,k_W>0\)` be the five positive primitive slope magnitudes and let `\(S:=k_\lambda^2+k_c^2+k_\gamma^2+k_U^2+k_W^2\)`. Let `\(a_1,\dots,a_5\)` be free real simplex weights.

- **M1 (gradient-optimal unit norm).** With `\(\mathbf a_5^{\rm grad}=(k_\lambda,k_c,k_\gamma,k_U,k_W)/\sqrt{S}\)`, the squared norm satisfies `\(\|\mathbf a_5^{\rm grad}\|^2 = (k_\lambda^2+k_c^2+k_\gamma^2+k_U^2+k_W^2)/S = 1\)`. Must verify the simplified residual `\(\|\mathbf a_5^{\rm grad}\|^2 - 1 = 0\)`.
- **M2 (max slope magnitude).** The slope at the gradient-optimal point `\(k_5(\mathbf a_5^{\rm grad})=\sum_p (k_p/\sqrt S)\,k_p = S/\sqrt S = \sqrt{S}\)`. Must verify `\(k_5(\mathbf a_5^{\rm grad}) - \sqrt{S} = 0\)` (under `\(k_p>0\)` so the positive square root is selected). Equivalently verify that `\(\sqrt S\)` is the constrained maximum of `\(\sum_p a_p k_p\)` subject to `\(\sum_p a_p^2 = 1\)` (Cauchy-Schwarz / Lagrange).
- **M3 (per-face first-order gap).** For each primitive axis `\(p\in\{\lambda,c,\gamma,U,W\}\)`, the gap between the five-coordinate squared max slope and the corresponding quadruple-face squared max slope equals exactly `\(k_p^2\)`: `\(S - (S - k_p^2) = k_p^2 > 0\)`. Must verify all five residuals `\((S) - (\text{face\_max\_sq}_p) - k_p^2 = 0`, and that each `\(k_p^2 > 0\)` under `\(k_p>0\)` (strict improvement over every face).
- **M4 (cross-leverage identity + Cauchy slack).** Define `\(w_\Sigma := 2\sum_{p<q} a_p a_q\)`. Two identities must hold for ALL real `\(a_1,\dots,a_5\)` (no simplex constraint needed for the identity itself): (a) `\(w_\Sigma - \big[(\sum_p a_p)^2 - \sum_p a_p^2\big] = 0\)`; (b) `\(\big[5\sum_p a_p^2 - (\sum_p a_p)^2\big] - \sum_{p<q}(a_p - a_q)^2 = 0\)`. Both must reduce to the zero polynomial under `Expand`/full simplification.
- **M5 (barycenter leverage = 4).** At `\(\mathbf a_5^{\rm eq}=(1,1,1,1,1)/\sqrt5\)` (which lies on the simplex, `\(\sum a_p^2 = 1\)`), `\(w_\Sigma(\mathbf a_5^{\rm eq}) = 2\cdot\binom{5}{2}\cdot\tfrac15 = 2\cdot10\cdot\tfrac15 = 4\)`. Must verify `\(w_\Sigma(\mathbf a_5^{\rm eq}) - 4 = 0\)`. (Optional but in-scope: verify the bound `\(w_\Sigma(\mathbf a)\le 4\)` on the simplex follows from M4(b)'s nonnegative slack, i.e. that the maximizer is the barycenter.)
- **M6 (fixed-simplex certified bracket form).** With `\(H_0>0\)`, `\(k>0\)`, `\(\kappa>0\)`, the certified root `\(\tau(H_0,k,\kappa)=\dfrac{2H_0}{k+\sqrt{k^2-2H_0\kappa}}\)` is the smaller positive root of the quadratic `\(\tfrac12 \kappa\,\tau^2 - k\,\tau + H_0 = 0\)`. Must verify that substituting `\(\tau\)` into `\(\tfrac12 \kappa \tau^2 - k\tau + H_0\)` simplifies to `0` (i.e. the closed form actually solves the bracketing quadratic), the genuine, non-tautological check on the rationalized-denominator form.

Each `M*` value above is the paper's stated value; do NOT substitute a different constant.

## F1 — missing_verification_script (subtype: script_doesnt_cover_claim)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.py`

**Issue:** The script computes the load-bearing quantities (lines 27-29 `a_grad`/`grad_norm_sq`/`k_grad`; line 43 `diffs`; lines 69/72 `identity_check`/`slack_check`; line 79 `w_eq`; line 95 `tau`) but contains zero assertions — only `out.append(...)` prints. Its "PASS" is unconditional. Add real failing checks so the stage is actually defended.

**Required change (state-what-must-hold, not a prescribed implementation):**
Add explicit assertions (use `assert` or a small `expect_zero`/`require` helper that raises / `sys.exit(1)` on failure) that enforce, at minimum, the manifest claims:
- M1: assert `grad_norm_sq` simplifies to `1` (the existing `grad_norm_sq` on line 28).
- M2: assert `k_grad - sqrt(S)` simplifies to `0` under the `positive=True` assumptions already declared on `k_lam,...,k_W` (the existing `k_grad` on line 29; `S = k_norm_sq`).
- M3: assert each of the five `diffs` (line 43) simplifies to the corresponding `k_p**2`, and that each `k_p**2` is strictly positive given the positive symbols.
- M4: assert `identity_check == 0` (line 69) and `slack_check == 0` (line 72) — these are already computed; just assert them.
- M5: assert `w_eq - 4` simplifies to `0` (the existing `w_eq` on line 79).
- M6: assert that substituting the existing `tau` (line 95) into `Rational(1,2)*kappa*tau**2 - k*tau + H0` simplifies to `0`. (This is the non-tautological check that the closed form solves the bracketing quadratic; do NOT merely re-print `tau`.)
Keep all existing prints; only ADD checks. Do not change any computed expression's form, and do not introduce any new physical constant — assert against the paper values already in the manifest. Do not weaken any assertion with assumptions broader than the already-declared `positive=True` on the slopes / `positive=True` on `H0,kappa,k`.

**Verification command:**
The verifier will run `redteam exec-sympy 216` and confirm the new assertions (covering M1-M6) appear in the file AND the script exits 0 with all checks passing.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.py`
- summary: Added failing SymPy checks for manifest claims M1-M6 while preserving the existing printed transcript.
- deviation: none

## F2 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl`

**Issue:** No Mathematica `.wl` exists for stage 216. The stage is non-status-only and non-checkpoint, and its math (gradient optimum, max-slope identity, per-face gap, cross-leverage + Cauchy-slack identities, barycenter leverage, bracket closed form) is fully and independently verifiable in native Mathematica. The dual-engine contract requires a second engine wherever Mathematica CAN verify; it can here. Create the `.wl`.

**Required change (acceptance criteria + claim manifest only — Codex designs the route):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl` that independently verifies claims M1-M6 (above). Each claim must be checked by a failing test — e.g. a project `expectZero[...]` helper or an `If[FullSimplify[...] =!= 0, Print[...]; Exit[1]]` guard — so the script exits nonzero if any identity breaks. The script must print a short labeled transcript and `Exit[0]` only when all checks pass.

**Anti-transliteration guard (MANDATORY):** The `.wl` must derive each result independently using native Mathematica primitives — `Solve`/`Reduce` (e.g. obtain `\(\mathbf a_5^{\rm grad}\)` and the max slope `\(\sqrt S\)` as the Lagrange-multiplier / constrained-maximum solution of `Maximize`/`Reduce` over `\(\sum a_p^2=1\)`, rather than asserting the closed form by hand; obtain the bracket `\(\tau\)` for M6 from `Solve` of the quadratic `\(\tfrac12\kappa\tau^2 - k\tau + H_0 = 0\)` and confirm it equals the stated closed form), `Series`+`Coefficient`, `D[]`, `Integrate`, `FindRoot`, matrix/quadratic-form ops (e.g. realize `\(w_\Sigma\)` and the Cauchy slack via the quadratic form `\(\mathbf a^T(\,\cdot\,)\mathbf a\)` / `Sum` and `Expand`) — via a DIFFERENT decomposition than the `.py`. A line-by-line port of the SymPy `out.append` choreography (same variable order, same intermediate `diffs`/`identity_check`/`slack_check` scaffolding mirrored into Mathematica syntax) is REJECTED as transliteration. In particular, M2 and M6 should be *derived* (via constrained maximization and via solving the quadratic, respectively), not merely re-stated and `FullSimplify`-checked the way the SymPy script forms them by hand.

Declare the slopes `\(k_\lambda,\dots,k_W\)` and `\(H_0,\kappa,k\)` as positive reals (matching the SymPy `positive=True`) where needed for the square-root branch in M2/M6. Use only the paper values in the manifest; introduce no new constant.

**Claim manifest:** M1-M6 as listed in the shared "Claim manifest" block above. All six must each have a corresponding independent, failing check in the `.wl`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 216` and confirm the new `.wl` exists, contains independent checks for M1-M6, and `math -script <path>` exits 0 with all checks passing.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_mathematica_audit.wl`
- summary: Created a native Mathematica audit that independently verifies manifest claims M1-M6 with failing checks and explicit pass/fail exits.
- deviation: none
