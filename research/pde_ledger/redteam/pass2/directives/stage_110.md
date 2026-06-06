---
unit_id: 110
batch: IV.2
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T01:53:25-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 110

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

You DESIGN the route. This directive states the requirement and acceptance criteria only; choose the actual independent construction yourself. Do NOT introduce unrelated features, refactors, or stylistic changes.

After editing, RUN the script (`timeout 600 math -script <path>`) and iterate until it exits 0 with all in-file checks passing. A timeout (exit 124) is a FAILURE — reformulate, never raise the cap. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, prose documents, or any numbering label. Do NOT edit the SymPy `.py` (it is the independent reference engine and stays unchanged).

## Context

The user has AUTHORIZED re-authoring this `.wl` to a genuinely independent route (dual-engine rule; IV.1-100 / V.3 / VI.1 precedent). The math is correct and the emitted values are right — the only defect is that the Mathematica engine re-types the SymPy computation.

## F1 — mathematica_transliteration

**Target:** `mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl:31-53`

**Issue:** The `.wl` is a line-by-line port of `scripts/...stage110..._sympy_audit.py:8-30`: identical `lambdaOut`, `lambdaR = lambdaOut + rho`, `yR = (-3 + rho)/lambdaR`, the SAME `Series[yR, {z, 0, 5}]` expansion on the same ratio (`sp.series(Y_R, z, 0, 6)`), the same `c2/c4/c5` (with `c5 = .../I`), `chiR = c5/(1/27)`, `chiRLinear = Series[chiR, {rho, 0, 2}]`, and the same five `expectZero` asserts in identical order with identical RHS. No structurally distinct route to `chi_Q^R`.

**Requirement:** Re-author the `.wl` so the Robin-outlet coefficients (`c2`, `c4`, `c5`) and `chi_Q^R = 3/(3 - rho)` (with its `rho`-linearization `1 + rho/3 + rho^2/9`) are reached by a construction genuinely independent of the `.py`. The `.wl` must NOT obtain `Y_R` via `Series[(-3+rho)/lambdaR, {z,0,5}]` on that ratio — that is the shared black-box. Use a structurally different native primitive (Codex's choice — e.g. an explicit geometric/polynomial inversion of `lambdaR * Y_R = (-3 + rho)` order-by-order via `Solve`/`CoefficientList`; or derive `chi_Q^R = 3/(3 - rho)` directly from the pole structure of `Y_R` / a residue or `Together`-then-coefficient route — note `chi_Q^R` is the z^5 odd coefficient normalized by `1/27`). The `rho != 3` non-degeneracy assumption stays.

**Acceptance criteria:**
1. The `.wl` no longer contains `Series[yR, {z, 0, 5}]` (or any equivalent series-of-the-`(-3+rho)/lambdaR`-ratio) as the route to the coefficients.
2. `c2`, `c4`, `c5`, `chi_Q^R = 3/(3 - rho)`, and the linearization `1 + rho/3 + rho^2/9` are all still asserted and still pass, with the SAME emitted values.
3. The route uses native Mathematica primitives in a way that is NOT a transliteration of the `.py` (the verifier checks independence: distinct primitive, distinct intermediate structure).
4. The committed Mathematica output remains byte-identical in its deliverable lines (method-only change); the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl`
- summary: Replaced the direct ratio series with a Mathematica polynomial-jet solve for the Robin outlet coefficients and a separate rho-jet solve for the linearization.
- deviation: none
