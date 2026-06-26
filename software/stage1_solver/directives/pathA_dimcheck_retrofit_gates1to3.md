# Directive — Dimensional-check retrofit for moving-throat Gates 1–3 (`pathA_30/31/32`)

**Status:** READY TO RUN (cold-start safe). Standalone; **independent of Gate 5** — does not touch any Gate-1–4 *physics* verdict.
**Tracking:** task **#98**. **Owner:** a fresh session (per user, 2026-06-26).
**One-line goal:** replace the **vacuous-but-green** dimensional checks in the three moving-throat ladder gate scripts with the
**genuine, free-carrier-independent, able-to-fail, dual-engine** check that Gate 4 (`pathA_33`) already uses — so a green dim check
*means* something.

---

## 0. Read first (cold-start context)

- **The regression that prompted this** (full statement): memory `feedback-dimensional-consistency-check` — the
  "⚠️ REGRESSION FOUND (2026-06-26)" + "⭐ Gate-4 RESOLUTION + the SHARPER trap" paragraphs. Also memory
  `feedback-negative-verdict-short-circuit` (the free-carrier-back-solve tautology + constant-`self_ablation` patterns).
- **The bar (the template to copy):** the Gate-4 dim check.
  - SymPy: `software/stage1_solver/tools/pathA_33_quadrupole_normalization_sympy.py` — `dim_of()` recursive walker (≈L313),
    `P0_physical = (c_s/a)**2 * (N0/D0)`, `dimensional_ok = (dim_of(P0_physical, …) == ZERO_DIM)` (≈L401–414), the
    corrupt-`[N0]` able-to-fail probe (≈L423, L500–503), the μ̂₀ kept ONLY as a labeled `"non-able-to-fail (mu_hat0 free
    carrier)"` diagnostic (≈L477), `base_verdict` mapping `("dimensional_ok", "QUAD_FAIL_DIMENSIONAL")` (≈L682), and a real
    two-verdict `self_ablation` via `ablation(...)` (≈L760).
  - Mathematica: `software/stage1_solver/tools/pathA_33_quadrupole_normalization.wl` — `dimOf[expr_, dims_]` recursive walker
    (≈L104), `zeroDim`, `dimAdd`, `dimScale`.
- **The real harness (reuse, don't reinvent):** `software/stage1_solver/src/stage1_solver/dimensional_check.py` (built in
  `pathA_18`). The Gate scripts must NOT hand-roll a fresh vacuous ledger again; reuse this harness OR the Gate-4 `dim_of`/`dimOf`
  pattern verbatim.
- **Why now / scope honesty:** the *physics* OOM-risk is LOW — the reference algebra in these gates is dimensionless-by-construction
  on a constant-coefficient reference, so we do not expect a hidden order-of-magnitude error. **This is a check-QUALITY / honesty
  retrofit:** the current green checks false-advertise. Do not over-claim a physics finding; the deliverable is checks that *can fail*.

---

## 1. What is wrong right now (precise, per script — verified 2026-06-26)

All three reinvented a typed-`{M,L,T}`-tuple ledger **structurally disconnected from the actual computed `sympy` expressions** (the
`dsolve` solutions, the `M_AB/K_AB` integrals, the `D_{A,n}` responses are never fed in). All three are **single-engine** (the `.wl`
does zero dimensional analysis). Specifics:

- **Gate 1 — `pathA_30_dn_unit_test_sympy.py`** (the least-bad): real tuple relations via `dimension_add`/`dimension_scale`, but
  `dim_K` is **reverse-engineered from the EOS answer** → tautological; **no able-to-fail row**; `.wl` does no dim analysis.
- **Gate 2 — `pathA_31_scalar_breathing_sympy.py`**: typed `Dim` tuples with TWO perturbation rows
  (`perturbed_Keta_dimension_fails`, `perturbed_alpha_L_dimension_fails`) — partially able-to-fail, but still **typed constants, not
  the assembled `M_AB/K_AB` expressions**; `.wl` does no dim analysis.
- **Gate 3 — `pathA_32_grouped_p2_isotropy_sympy.py`** (the worst): `build_dimensional_table()` just emits a typed `dimension_MLT`
  table per quantity. **No able-to-fail row; no `a`/`c_s`/angular-Jacobian symbol anywhere** — in the gate that most loudly demanded
  the `a⁵`/angular-Jacobian check; `.wl` does no dim analysis.

---

## 2. Requirement (what "fixed" means — acceptance criteria)

For **each** of `pathA_30`, `pathA_31`, `pathA_32`, replace the dim check so it satisfies ALL of:

1. **Operates on the ACTUAL computed expressions.** Give the model's symbols real dimensions
   (`a, c_s, ħ, m, K, ρ0, L0, μ_η, T_w, K_η, T_Ω, β₂, …` — sourced from the action / EOS `P=Kρ⁵`, NOT back-solved), then run the
   recursive `dim_of`/`dimOf` walker **on the assembled headline `sympy`/`wl` objects** the gate actually computes (see §3 for the
   per-gate headline quantity), not on a hand-typed tuple ledger.
2. **Verdict-bearing gate is FREE-CARRIER-INDEPENDENT.** The boolean that gates the verdict must check a quantity whose dimension is
   **fully fixed by sourced inputs** (homogeneity / a required dimensionless argument / a required `[ω²]=T⁻²` ratio — see §3). It must
   NOT be made to pass by solving a free/unconstrained parameter for an unknown dimension. (This is the Gate-4 lesson: μ̂₀-style
   back-solve = dead code. Keep any such diagnostic ONLY with the explicit label `"non-able-to-fail"`.)
3. **Carries an able-to-fail PROBE that corrupts a SOURCED input.** Perturb one *sourced* symbol's dimension (e.g. `[T_w]→[T_w]·L`)
   and assert the gate flips to FAIL; with the perturbation removed, it must pass. (Detection test for the free-carrier trap: corrupt
   a *sourced input dimension* — not a carrier — and the gate must FAIL. If corrupting a sourced input passes silently, the teeth are
   dead.) Prefer a real `self_ablation` (two recomputed verdicts `with_mutation` / `without_mutation`), exactly like Gate-4's
   `ablation(...)` — not a constant field.
4. **DUAL-ENGINE.** The `.wl` must independently run the `dimOf` walker on the headline quantity and reach the same dimensional
   verdict (copy the Gate-4 `.wl` `dimOf/zeroDim/dimAdd/dimScale` pattern). No more single-engine dim checks.
5. **The verdict string actually depends on it.** Wire `dimensional_ok` into the gate's `base_verdict`/`compute_verdict` so a FALSE
   value yields a real `FAIL_DIMENSIONAL`-class verdict (prove load-bearing by killing the residual verdict path and confirming the
   verdict still flips, per the standing adversarial-ablation rule).

**Out of scope / do NOT change:** the gate *physics* verdicts stand (`DN_UNITTEST_BC_DEPENDENT`, `BREATHING_CALIBRATED`,
`ISOTROPY_CALIBRATED`) and all their other gates/probes. If a genuine dimensional inhomogeneity is discovered (not expected), STOP and
escalate — that would be a real physics finding, handled separately, not silently folded.

---

## 3. The load-bearing quantity to check, per gate (requirement-level; Codex designs the route)

State the target the genuine gate must check; **do not pre-design the script** (Codex codes, Claude reviews —
`feedback-claude-reviews-codex-codes`). Each target's dimension is fixed by sourced inputs, so each is genuinely able-to-fail.

- **Gate 1 (`pathA_30`):** the speed law `c_s² = 5Kρ0⁴/m` must carry `[c_s²]=L²T⁻²` **from independently-sourced `[K],[ρ0],[m]`**
  (source `[K]` from `[P]=[energy density]` via `P=Kρ⁵` — do NOT reverse-engineer `[K]` from the c_s answer), AND the DtN
  `Z₀₀(ω)=−(ω/c_s)tan(ωL0/c_s)` must have a **dimensionless `tan` argument** `ωL0/c_s`. Able-to-fail: corrupt `[K]`, `[L0]`, or
  `[ρ0]` → the `tan` argument is no longer dimensionless / `[c_s²]≠L²T⁻²` → FAIL.
- **Gate 2 (`pathA_31`):** the assembled mass/stiffness matrices `M_AB`, `K_AB` (from operator projection) must each be
  dimensionally homogeneous across their entries, AND the eigen-ratio `K_AB/M_AB` must carry `[ω²]=T⁻²`. Able-to-fail: corrupt
  `[μ_η]` or `[T_w]` → `M_AB`/`K_AB` inhomogeneous or `K/M≠T⁻²` → FAIL. (Run `dim_of` on the REAL assembled `M_aa, M_aL, M_LL,
  K_aa, …` objects, not typed tuples.)
- **Gate 3 (`pathA_32`):** the angular-sector integrals `M₂=∫μ_η β₂² dV` and `K₂=∫[T_w β₂'² + (K_η+6 T_Ω) β₂²] dV` must be
  **term-by-term homogeneous** (every term in `K₂` shares one dimension; the `l(l+1)=6` angular-stiffness term `T_Ω` matches the
  `T_w β₂'²` and `K_η β₂²` terms), the **angular/radial measure (the Jacobian, incl. any `a`-power)** must carry its dimension
  explicitly, and `K₂/M₂` must be `[ω²]=T⁻²`. Able-to-fail: corrupt `[T_Ω]` or `[T_w]`, or drop the Jacobian dimension → `K₂`
  inhomogeneous → FAIL. **This is the gate the original check most failed** (no `a`/`c_s`/Jacobian symbol at all) — make the
  Jacobian/measure explicit here.

---

## 4. Process (the standing gauntlet — unchanged)

1. **Directive is this file.** A fresh session reads §0–§3, then drafts/edits the three scripts via Codex.
2. **Codex codes + runs + iterates to exit 0**, dual-engine, `--sandbox danger-full-access` for any `.wl` that must RUN; `codex exec
   … -c model_reasoning_effort=xhigh`, backgrounded, never wrapped in shell `timeout`; per-script `timeout 600` on the scripts.
3. **Claude reviews** (review only; Codex applies fixes). Review ordering: iterate Codex to green → one GLM pass (user-run) → fold →
   Codex to green.
4. **Tri-review per script** (clean agents): orchestrator arbiter re-run (both engines) + transliteration-fidelity audit +
   **adversarial-with-ablation**. **NEW explicit tri-review item (do not wave through):** *"the dim check operates on real assembled
   expressions, is free-carrier-independent, has an able-to-fail probe that fires on a corrupted SOURCED input, and is dual-engine."*
   The adversarial leg must itself corrupt a sourced input and confirm the verdict flips (the arbiter + fidelity legs systematically
   miss pass-by-construction — `feedback-negative-verdict-short-circuit`).
5. **User gate**, then commit (commit only when the user asks; stage explicit paths). One commit for the retrofit is fine
   (`pathA_30/31/32 dim-check retrofit: real-expression + able-to-fail + dual-engine`), or squash per-script — orchestrator's call
   with the user.
6. **On completion:** flip task #98 → done; update `STATUS.md` (the "⚠️ Tracked debt" bullet in the ⭐⭐ LATEST block and the
   parallel line in `decisions/13` §0) from "VACUOUS / retrofit debt" → "retrofitted (real-expression, able-to-fail, dual-engine)";
   update the Gate-1 "Known NITs" line (2) in the completion ladder (`research/pde_ledger/notes/stages/
   moving_throat_pde_completion_ladder.md`) which currently blesses the hand-coded tuples.

---

## 5. Related (do NOT do here unless asked) — the χ_Q reconciliation

STATUS.md flags a *separate* tracked doc item alongside this retrofit: `pathA_22b` Gate 3 computed `χ_Q ≈ 0.712` (numeric,
minimal-combination) while moving-throat Gate 4 (`pathA_33`) derived `χ_Q = 1` (outgoing-DtN Hankel). These are different
computations bearing the same name. Reconciling/relabeling them is a documentation task — **out of scope for this retrofit**; noted
here only so the next session knows it is tracked and distinct.

---

## 6. Cross-refs

- Memories: `feedback-dimensional-consistency-check` (the rules + the free-carrier trap), `feedback-negative-verdict-short-circuit`
  (pass-by-construction patterns + ablation backstop), `feedback-claude-reviews-codex-codes`, `feedback-dual-engine-required`.
- Template: `software/stage1_solver/tools/pathA_33_quadrupole_normalization_{sympy.py,.wl}` (§2.5 of its directive
  `directives/pathA_33_quadrupole_normalization.md`).
- Harness: `software/stage1_solver/src/stage1_solver/dimensional_check.py`.
- Targets: `software/stage1_solver/tools/pathA_30_dn_unit_test_*`, `pathA_31_scalar_breathing_*`,
  `pathA_32_grouped_p2_isotropy_*`.
- Front door: `STATUS.md` (⭐⭐ LATEST — the "Tracked debt" bullet); `decisions/13` §0.
