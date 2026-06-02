---
unit_id: 123
batch: retro
created_at: 2026-06-01T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-01T21:53:51Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 123 (retro-sweep: dual-engine .wl)

> This is a RETRO-SWEEP directive. Stage 123 was audited + verified in batch IV.3 (SymPy-only;
> its IV.3 notes-side `228 → 160` typo fix and banner relabel are already done). Under the
> dual-engine rule (a Mathematica `.wl` is REQUIRED on every stage Mathematica CAN independently
> verify), it is missing its second engine. The ONLY change is to ADD the `.wl`.
> (The prior IV.3 directive content for this unit is preserved in git history. The IV.3
> paper_misalignment is resolved — do NOT reopen it, do NOT touch notes/paper.)

The SymPy audit script for this stage is correct and is the REFERENCE engine. Do NOT modify it. Do NOT touch `paper.tex` or `notes/`. The only required change is the dual-engine gap below.

After creating the script, RUN it (`timeout 600 math -script <path>`) and iterate until it exits 0 with all checks passing. A timeout (exit 124) is a FAILURE — reformulate the math, never raise the cap. The orchestrator independently re-runs afterward.

## F1 — missing_mathematica

**Issue:** Stage 123 ("parent-normalized branch values") is dual-engine-capable — every load-bearing claim is a closed-form algebraic identity (two single-variable `Solve` inversions, radical simplification, and exact-numeric evaluation at explicit branch points); no integral, no BVP, no transcendental root — but it has no Mathematica `.wl`. Under the dual-engine rule an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the SymPy script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_123.tex` and the stage notes file are the source of truth for the math. (The card / notes use older "Stage 140 / 221" numbering — anchor to the live `scripts/` owners listed below. Do NOT edit the card or notes.)
- Use Mathematica-NATIVE primitives (`Solve`/`Reduce` with reality/branch filtering, `FullSimplify`/`RootReduce`/`PowerExpand` under `$Assumptions`, exact `N[...]`) via a DIFFERENT derivation route than the SymPy script (which uses `sp.solve(...)[0]`) — NOT a line-by-line port. Reference an existing verified `.wl` ONLY for house idioms (the `expectZero` helper that `Exit[1]`s on failure, `$Assumptions`, `stripCE`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion the `.wl` derives must match the SymPy result.

**Upstream-owned anchors (import as given; cite, do NOT re-derive here):**
- `K_s = 3π a²ℏ²/(5 m ρ ℓ)`, `K_q = (Z_q/μ₀)π²c_s²/(4 L²)`, `J_s = 4π a²ℓ/3` — owned by stage 118.
- `λ = −(8 Sqrt(2)/3) q v_w0 a² ℓ Sqrt(L)` — sign owned by stage 118 (λ < 0).
- `r_F1 = Sqrt(12 R²/π² − 1)`, `R = 37/20` — owned by stage 121.
- `g_± = (2 Sqrt(4107−100π²) ± 37 Sqrt(3))/(20π)` — owned by stage 122.
- healing lock `c_s → ℏ/(2 m ℓ)` (applied only inside the `Ξ_T` inversion) — owned by stage 118.

**Claim manifest** (match the SymPy PASS-label so the verifier can pair them):
- **M1** — `Xi_v law`: inverting `𝔯 = λ/Sqrt(K_s K_q)` for `v_w0`, substituting into the `Ξ_v` definition, and simplifying yields `Ξ_v(𝔯) = −(3 Sqrt(30) π^(3/2)/160)·𝔯`. **The leading minus sign carries the un-squared λ-sign — it must survive.**
- **M2** — `Xi_T law`: inverting `𝔤 = Sqrt(2 Z_q K_s)/(𝒯_m J_s c_s Sqrt(μ₀ L))` for `𝒯_m`, applying the healing lock `c_s → ℏ/(2 m ℓ)`, substituting into the `Ξ_T` definition, and simplifying yields `Ξ_T(𝔤) = (3 Sqrt(30)/(10 Sqrt(π)))·(1/𝔤)`.
- **M3** — `Xi_v(F1)` numeric: evaluating M1 at `𝔯 = r_F1` gives `Ξ_v^{F1} ≈ −1.01675633282526`. (Printed in the `.py`; the `.wl` must ASSERT it against the derived law AND reproduce the ≥15-digit decimal via `RootReduce` then exact `N[..., 20]`, making it a real PASS check.)
- **M4** — `Xi_T(nat/-/+)` numerics: `Ξ_T^{nat} = 3 Sqrt(30)/(10 Sqrt(π)) ≈ 0.927058084855655`, `Ξ_T^{(−)} = Ξ_T^{nat}/g_− ≈ 1.22297517701464`, `Ξ_T^{(+)} = Ξ_T^{nat}/g_+ ≈ 0.331334521644609`. (Same: assert against derived forms + reproduce decimals.)

**Route / correctness requirements (acceptance criteria — CRITICAL, you choose how to satisfy them):**
- **λ-sign (load-bearing).** `𝔯 = λ/Sqrt(K_s K_q)` is LINEAR in λ (un-squared), so the negative λ-sign propagates into M1's leading `−`. Use the SAME negative-λ convention as stage 118. Do NOT use λ² and do NOT drop the minus — either would silently flip `Ξ_v^{F1}` to `+1.0168` and mis-pass. (Contrast: `r_c = λ²/(K_s K_q)` upstream IS squared — do not confuse the two.)
- **`v_w0` is REAL, not positive** (λ < 0 forces a sign relationship). Mirror the SymPy declaration: `v_w0` real-only; all other symbols (`a, L, ℓ, μ₀, Z_q, m, ρ, q, ℏ, c_s, 𝒯_m, 𝔯, 𝔤`) positive. Do NOT assume `v_w0 > 0` or the inversion branch / radical simplification can pick the wrong sign.
- **Branch selection.** Both laws come from a single picked `Solve` root in the `.py`; the Mathematica route must pick the corresponding real branch deterministically (reality/positivity filtering), and anchor the expected closed form so a wrong-branch pick FAILS the identity.
- The healing lock `c_s → ℏ/(2 m ℓ)` is applied ONLY inside the `Ξ_T` inversion (it gives `Ξ_T` its `Sqrt(π)` instead of `c_s`) — cite the stage-118 relation, do not improvise a different `c_s` substitution.

**Anti-transliteration:** the `.wl` must NOT re-type `Ξ_v = −(3 Sqrt(30) π^(3/2)/160)·𝔯`, `Ξ_T = (3 Sqrt(30)/(10 Sqrt(π)))·(1/𝔤)`, or the four numeric targets and check them against themselves. The laws must EMERGE from the `v_w0` / `𝒯_m` inversions + simplification, and the numerics must be COMPUTED from the derived laws at the upstream-anchored `r_F1`/`g_±`. Importing the upstream building blocks (`K_s, K_q, J_s, λ, r_F1, g_±`) as givens is expected (stages 118/121/122 own + verify them); the independent work is the two inversions + exact evaluation via a route distinct from `sp.solve(...)[0]`.

**Comment hygiene:** avoid any `*)` substring inside Mathematica comments (premature comment close).

**Verification command:** the verifier runs `redteam exec-mathematica 123`, confirms exit 0 with all PASS lines (≥ 6 substantive checks: M1, M2, M3, M4×3), and reviews that the `.wl` is a genuinely independent route (native primitives, different decomposition; M1 retains the leading minus from un-squared λ) whose conclusions agree with the SymPy engine.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl`
- summary: Added a Mathematica dual-engine audit that derives the Xi_v and Xi_T laws via Reduce-based inversions and checks the F1, natural, minus, and plus branch values exactly with high-precision prints.
- deviation: none
