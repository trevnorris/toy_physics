---
unit_id: 107
batch: IV.2
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T07:53:37Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 107

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

You DESIGN the route. This directive states the requirement and acceptance criteria only; choose the actual independent construction yourself. Do NOT introduce unrelated features, refactors, or stylistic changes.

After editing, RUN the script (`timeout 600 math -script <path>`) and iterate until it exits 0 with all in-file checks passing. A timeout (exit 124) is a FAILURE — reformulate, never raise the cap. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, prose documents, or any numbering label. Do NOT edit the SymPy `.py` (it is the independent reference engine and stays unchanged).

## Context

The user has AUTHORIZED re-authoring this `.wl` to a genuinely independent route (dual-engine rule; IV.1-100 / V.3 / VI.1 precedent). The math is correct and the emitted values are right — the only defect is that the Mathematica engine re-types the SymPy computation, so it is not a true second engine.

## F1 — mathematica_transliteration

**Target:** `mathematica/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.wl:33-76`

**Issue:** The `.wl` is a line-by-line port of `scripts/...stage107..._sympy_audit.py:27-74`: identical `lambdaDef = Expand[sNorm*(lambdaOut /. z -> beta*z) + ...]`, identical `l0/l2/l4/l5` coefficient extraction, the SAME `Series[l0/lambdaDef, {z, 0, 5}]` normalized-ratio expansion the `.py` performs (`sp.series(L0/Lambda_def, z, 0, 6)`), the identically-typed `yFormula`, and the same `Solve[{m2==1/9, m4==4/81}, {sigma2, sigma4}]` and assertion order. There is no structurally distinct route to `chi_Q`.

**Requirement:** Re-author the `.wl` so the normalized deformed branch `Y_def` and the general `chi_Q = 3(sNorm*beta^5 + 9 sigma5)/(3 sNorm - sigma0)` are reached by a construction genuinely independent of the `.py`. In particular the `.wl` must NOT obtain `Y_def` via `Series[l0/lambdaDef, ...]` on the normalized ratio — that is the shared black-box. Use a structurally different native primitive (Codex's choice — e.g. an explicit polynomial inversion of the operator identity `Lambda_def * Y = L0` order-by-order via `Solve`/`CoefficientList` for the branch coefficients, the same independence pattern accepted at stage 105's deformed branch and stage 109; or a residue/partial-fraction route to the same coefficients). The even-matching `Solve` for `(sigma2, sigma4)` may remain, but the branch-coefficient / `chi_Q` derivation feeding it must be the independent route.

**Acceptance criteria:**
1. The `.wl` no longer contains `Series[l0/lambdaDef, ...]` (or any equivalent series-of-the-normalized-ratio) as the route to `Y_def`/the branch coefficients.
2. `Y_def`, `m2`, `m4`, `chi_Q`, the `Sigma2`/`Sigma4` exact formulas, and the general `chi_Q` factorization are all still asserted and still pass, with the SAME emitted values.
3. The route uses native Mathematica primitives in a way that is NOT a transliteration of the `.py` (the verifier will check independence: distinct primitive, distinct intermediate structure, not just renamed symbols).
4. The committed Mathematica output remains byte-identical in its deliverable lines (method-only change); the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.wl`
- summary: Replaced the normalized-ratio series with an order-by-order Mathematica coefficient solve for the deformed branch while preserving the existing assertions and emitted output.
- deviation: none
