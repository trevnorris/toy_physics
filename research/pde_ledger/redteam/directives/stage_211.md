---
unit_id: 211
batch: VI.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-02T10:45:38-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 211

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl` (new file — `.wl` lives under `mathematica/`, NOT `scripts/`)

**Issue:** Stage 211 is non-status, non-checkpoint, and currently has only a SymPy engine. The project dual-engine contract requires an independent Mathematica verification because every claim here is a closed-form symbolic identity that Mathematica can derive natively. The new script must verify the claim manifest M1–M6 below and exit 0.

**Required change:**
Author a new Mathematica script at the target path that independently verifies M1–M6. Each claim must terminate with a hard guard (e.g. `If[Simplify[lhs - rhs] =!= 0, Print["FAIL: <name>"]; Exit[1]]`, or an `expectZero`-style helper), so the script fails loudly if any identity breaks. Print a clear per-section banner and a final "All Stage 211 identities verified." line. Use `Stage 211` (not 194) in the banner.

**Anti-transliteration guard (REQUIRED):** The `.wl` must derive these results independently using native Mathematica primitives — `D[]`, `Simplify`/`FullSimplify`/`Together`/`Numerator`, `Expand`, `Exponent`/`MonomialList`/`CoefficientRules` for total degree, and `/.` substitution — via a DIFFERENT decomposition than the `.py`. Specifically:
- For M1, do NOT re-type the SymPy expressions for `M_r, M_s, L_r, L_s, N_r, N_s`. Instead define `Phi` and `tau` directly from the physical formulas (paper eqs. for `Phi` and `tau`), compute `D[Phi, r]` and `D[Phi, s]` symbolically, put each over a common denominator with `Together`, extract the numerator with `Numerator`, and check that this numerator equals the paper's `N_*` form `2 M_* Sqrt[Delta] + L_*` where `M_*, L_*` are defined from the paper's boxed formulas (notes section 2). The check is "the directly-differentiated numerator equals the paper's claimed stationary numerator," derived in Mathematica's own algebra.
- For the degree checks (M2/M3), obtain total degrees with `Max[Total /@ CoefficientRules[poly, {r, s}][[All, 1]]]` or `Exponent`-based logic, not a line-for-line port of `sp.Poly(...).total_degree()`.

A line-by-line port of the SymPy choreography is rejected as `mathematica_transliteration`.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1 — Exact stationary numerator law.** With
  `Delta(r,s) = A + B r + C s + D r^2 + E r s + F s^2`,
  `Phi(r,s) = (k_i + k_j r + k_k s + Sqrt[Delta]) / Sqrt[1 + r^2 + s^2]`,
  `M_r = (1+r^2+s^2) k_j - r(k_i + k_j r + k_k s)`,
  `M_s = (1+r^2+s^2) k_k - s(k_i + k_j r + k_k s)`,
  `L_r = (1+r^2+s^2) D_r[Delta] - 2 r Delta`,  `L_s = (1+r^2+s^2) D_s[Delta] - 2 s Delta`,
  `N_r = 2 M_r Sqrt[Delta] + L_r`,  `N_s = 2 M_s Sqrt[Delta] + L_s`,
  prove
  `D[Phi, r] == N_r / (2 (1+r^2+s^2)^(3/2) Sqrt[Delta])` and
  `D[Phi, s] == N_s / (2 (1+r^2+s^2)^(3/2) Sqrt[Delta])`
  (i.e. `Simplify[D[Phi,r] - N_r/(2(1+r^2+s^2)^(3/2) Sqrt[Delta])] == 0`, same for `s`).

- **M2 — Quartic cross-consistency.** With `C_cross = M_s L_r - M_r L_s`, prove the square-root-free identity
  `M_s N_r - M_r N_s - C_cross == 0`, AND that `C_cross` has total degree exactly 4 in `(r,s)`.

- **M3 — Sextic square eliminants.** With `S_r = L_r^2 - 4 M_r^2 Delta` and `S_s = L_s^2 - 4 M_s^2 Delta`, prove
  `N_r (N_r - 4 M_r Sqrt[Delta]) - S_r == 0`,  `N_s (N_s - 4 M_s Sqrt[Delta]) - S_s == 0`,
  AND that both `S_r` and `S_s` have total degree exactly 6 in `(r,s)`.

- **M4 — Bezout bound.** `(total degree of C_cross) * (total degree of S_r) == 24`. Assert the product equals 24 (this must be computed from the degrees obtained in M2/M3, not hardcoded as `4*6`).

- **M5 — Diagonal-isotropic reduction.** With the substitution
  `A -> k_i^2 - 2 H_0 u`, `B -> 2 k_i k_j`, `C -> 2 k_i k_k`, `D -> k_j^2 - 2 H_0 u`, `E -> 2 k_j k_k`, `F -> k_k^2 - 2 H_0 u`
  (i.e. `u_{ij}=u_{ik}=u_{jk}=0`, `u_{ii}=u_{jj}=u_{kk}=u`), prove
  `Delta_iso == (k_i + k_j r + k_k s)^2 - 2 H_0 u (1 + r^2 + s^2)`, AND, with `k_rs = (k_i + k_j r + k_k s)/Sqrt[1+r^2+s^2]`,
  `tau_iso == 2 H_0 / (k_rs + Sqrt[k_rs^2 - 2 H_0 u])`
  (where `tau = 2 H_0 / Phi` and `tau_iso` is `tau` under the isotropic substitution).

- **M6 — Full-symmetry equal-mix stationarity.** With `k_i = k_j = k_k = k` and the permutation-symmetric envelope
  `A -> k^2 - 2 H_0 u_d`, `B -> 2 k^2 - 4 H_0 u_x`, `C -> 2 k^2 - 4 H_0 u_x`, `D -> k^2 - 2 H_0 u_d`, `E -> 2 k^2 - 4 H_0 u_x`, `F -> k^2 - 2 H_0 u_d`,
  prove `N_r(1,1) == 0` and `N_s(1,1) == 0` (evaluate the symmetric `N_r, N_s` at `r=1, s=1`).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 211` and confirms the new `.wl` exists, contains a hard check for each of M1–M6, and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_mathematica_audit.wl`
- summary: Added the Stage 211 Mathematica audit that independently verifies M1-M6 with hard-failing symbolic checks.
- deviation: none

## F2 — stale_output (cosmetic stage-number banner)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py:35`

**Issue:** The SymPy banner prints `STAGE 194 — ...` (line 35), a leftover from an earlier renumbering, while the rest of the script and output correctly say "Stage 211". Cosmetic only; no assertion or math changes.

**Required change:**
At line 35, replace
`banner("STAGE 194 — FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION")`
with
`banner("STAGE 211 — FULL INTERIOR SIMPLEX OPTIMIZER AND FINITE ALGEBRAIC CANDIDATE REDUCTION")`.
Change nothing else. Do NOT touch the notes prose (out of red-team scope).

**Verification command:**
After Codex re-runs `python3 .../moving_throat_pde_stage211_..._sympy_audit.py`, the first banner line of the output reads `STAGE 211 — ...` and the script exits 0 with all checks passing.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage211_full_interior_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py`
- summary: Corrected the SymPy audit banner from Stage 194 to Stage 211.
- deviation: none
