---
unit_id: 240
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-03T08:35:57-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 240

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.py:113-115`

**Issue:** The notes §2.1 claim that the pole location `Omega_Q` does not enter the static loading-ratio extraction. The script "verifies" this with `assert_zero(sp.diff(c0_expr, Omega_Q))` and `assert_zero(sp.diff(c1_expr, Omega_Q))`. But `c0_expr = sp.simplify(1/rho_alpha) = alpha_mix/alpha_req` (line 88) and `c1_expr = (alpha_req - alpha_mix)/alpha_req` (line 89) — neither contains the symbol `Omega_Q`. The derivatives are therefore identically zero *by construction*, so both `assert_zero` calls are no-ops that cannot fail and provide no evidence for the §2.1 claim.

**Required change:**
Make the `Omega_Q`-independence check operate on weights *extracted from the `Omega_Q`-bearing precursor*, not on expressions already written without `Omega_Q`. Specifically:

1. Start from the precursor object that actually contains `Omega_Q` — `Y_support` (defined at line 83 as `alpha_mix/alpha_req + (alpha_req - alpha_mix)/alpha_req * pole`, with `pole = 1/(1 - omega**2/Omega_Q**2)`).
2. Extract the static contact weight and the pole weight *from* `Y_support` by a route that itself sees `Omega_Q` (e.g. the `omega -> 0` value gives the full static sum, the `pole`-coefficient gives the pole weight; or a partial-fraction split of `Y_support` in `omega`). Call these `c0_static`, `c1_static`. (You design the extraction; do not just re-type `alpha_mix/alpha_req`.)
3. Assert:
   - `assert_zero(sp.diff(c0_static, Omega_Q), "Omega_Q independence of extracted c0")`
   - `assert_zero(sp.diff(c1_static, Omega_Q), "Omega_Q independence of extracted c1")`
   where `c0_static`/`c1_static` are the extracted objects (so the derivative is a *non-trivial* zero — the extraction itself touched `Omega_Q`).
   - `assert_zero(sp.simplify(c0_static - c0_expr), "extracted c0 matches compiler c0")`
   - `assert_zero(sp.simplify(c1_static - c1_expr), "extracted c1 matches compiler c1")`

Self-test trap to honor: the whole point is that the differentiated expression must syntactically depend on `Omega_Q` *before* the static extraction collapses it — if your extraction yields something that never contained `Omega_Q`, you have reproduced the bug. Confirm `c0_static` is derived from `Y_support` (which contains `Omega_Q`) and that, prior to simplification, the extraction path referenced `Omega_Q`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 240` and confirms the new labels (`Omega_Q independence of extracted c0/c1`, `extracted c0/c1 matches compiler`) appear, the extracted weights are derived from an `Omega_Q`-bearing object, and the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.py`
- summary: Replaced the tautological Omega_Q derivatives with c0/c1 extracted from the Omega_Q-bearing Y_support precursor via static and pole limits.
- deviation: none

## F2 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_mathematica_audit.wl`

(Note: filename mandated by policy — must be the `scripts/` stem with `_mathematica_audit.wl` and live in `mathematica/`.)

**Issue:** Unit 240 has only a SymPy engine. Every claim is a closed-form algebraic / rational-function identity that Mathematica can independently verify. The dual-engine policy requires a second engine wherever it is *possible*. Write a NEW Mathematica audit script.

**Independence constraint:** This must be an INDEPENDENT re-derivation via native Mathematica primitives using a DIFFERENT decomposition than the `.py` — not a transliteration. Do not mirror the SymPy variable choreography line-by-line. Use Mathematica-native routes, e.g.: build the contact/pole split with `Apart`/`Together` rather than re-typing the rational forms; obtain the static weights via `Limit`/`Series` in `omega` (around `omega -> 0` and the pole); verify the regime interval with `Reduce`/`Resolve` or exact numeric comparison. Use an `expectZero[...]` helper (with `FullSimplify` and `ConditionalExpression`-stripping per project idiom) and exit nonzero on any failure.

**Claim manifest** (each must have its own in-file check):

- **M1** — Product-ratio identity: with `Pi_tr = NQ*alpha_req/beta0`, `C_mix = NQ*alpha_mix/beta0`, show `Pi_tr/C_mix == alpha_req/alpha_mix =: rho_alpha`. Also verify the spectral substitution `NQ = mhat_-^2 * beta0 * s_-/lambda_-` leaves the ratio unchanged.
- **M2** — Contact-plus-pole compiler: with `Y = alpha_mix/alpha_req + ((alpha_req-alpha_mix)/alpha_req)/(1 - omega^2/Omega_Q^2)`, extract (independently, via `Apart`/`Limit`, NOT by re-typing) `c0 = 1/rho_alpha`, `c1 = (rho_alpha-1)/rho_alpha`; verify `c0 + c1 == 1`, the inverses `rho_alpha == 1/c0 == 1/(1-c1)`, and `zeta_req == c1/c0 == rho_alpha - 1`. Verify the extracted `c0,c1` are free of `Omega_Q` (derivative of the *extracted* weight w.r.t. `Omega_Q` is zero — same substantive point as SymPy F1, but via the Mathematica extraction).
- **M3** — Minimal isotropic specialization: `c0 = 3/4, c1 = 1/4` give `rho_alpha == 4/3` (from both `1/c0` and `1/(1-c1)`) and `zeta_req == 1/3`.
- **M4** — Selected demand product: `Pi_tr == (4/3) C_mix` and `S_req == Pi_tr/C_mix == 4/3` on the minimal-isotropic branch. Derive `Pi_tr` from `rho_alpha * C_mix` with `rho_alpha` taken from the M3 extraction (do NOT define `Pi_tr := (4/3) C_mix` and then "confirm" it — that reproduces the SymPy A9 redundancy).
- **M5** — Support-selector reduction: with `C_mix = 8 Lambda (1 - eps_star)/Pi^2` and `Pi_tr = (4/3) C_mix`, show `varrho == Pi^2 Pi_tr/(16 Lambda) == 2(1 - eps_star)/3`.
- **M6** — Regime classification: confirm `1 < 4/3 < 2`, equivalently `C_mix < Pi_tr < 2 C_mix` strictly, placing the branch in the symmetric-lowest-twin sector.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 240` and confirms the new `.wl` exists at the path above, each of M1-M6 has a corresponding `expectZero`/`If[...Exit[1]]` check, the derivation route differs from the `.py` decomposition, and the script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_mathematica_audit.wl`
- summary: Added a Mathematica audit with independent M1-M6 algebraic checks using Together/Apart, limit-based coefficient extraction, and strict regime validation.
- deviation: none
