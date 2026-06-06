---
unit_id: 114
batch: IV.2
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T07:45:21Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 114

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

You DESIGN the route. This directive states the requirement and acceptance criteria only; choose the actual independent construction yourself. Do NOT introduce unrelated features, refactors, or stylistic changes.

After editing, RUN the script (`timeout 600 math -script <path>`) and iterate until it exits 0 with all in-file checks passing. A timeout (exit 124) is a FAILURE — reformulate, never raise the cap. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, prose documents, or any numbering label. Do NOT edit the SymPy `.py` (it is the independent reference engine and stays unchanged).

## Context

The user has AUTHORIZED re-authoring this `.wl` to a genuinely independent route (dual-engine rule; IV.1-100 / V.3 / VI.1 precedent). The math is correct and the emitted values are right — the only defect is that the Mathematica engine re-types the SymPy computation. (Note: the downstream stage 137 mirrors this stage; 137's independence comes from using a DIFFERENT primitive than its own `.py` — apply the same idea here.)

## F1 — mathematica_transliteration

**Target:** `mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl:33-52`

**Issue:** The `.wl` reaches the Schur complement `delta_Lambda(D)` via `deltaD = Apart[c.Inverse[m].c, dSym]` — the same built-in-matrix-inverse primitive on the same matrix `m = {{kS,lam},{lam,-kQ*dSym}}` as `scripts/...stage114..._sympy_audit.py:30` `delta_D = apart((c.T * M.inv() * c)[0], D)`. The downstream `rho_c`/`r_c`/`sigma_tilde`/`sigma_c`/`kappa_c`/`gamma_c` definitions, the `target_D` Schur form, and the two `expectZero` asserts are all in identical order. The only "independence" is two CAS each calling their built-in matrix inverse — not a structurally distinct derivation.

**Requirement:** Re-author the `.wl` so `delta_Lambda(D)` (the core-level outlet response `c^T M^{-1} c`) is derived WITHOUT the built-in matrix-inverse primitive. Use an explicit elimination instead (Codex's choice — e.g. introduce the core unknowns `(s, q)`, write the 2×2 linear system `M.(s,q) = (g_s, g_q)` explicitly, `Solve` it for `(s, q)` by elimination, and form the response `delta = g_s*s + g_q*q`; or eliminate `q` by hand and back-substitute to build the Schur complement directly). The deliverable identifications (`rho_c`, `sigma_c`, `kappa_c`, `gamma_c`, `r_c`) and both identity asserts (`Schur form identity`, `low-frequency normalized outlet identity`) must be preserved.

**Acceptance criteria:**
1. The `.wl` no longer uses `Inverse[m]` (or `LinearSolve`/`PseudoInverse` as a thin wrapper that is just the matrix inverse) as the route to `delta_Lambda(D)`; it reaches the Schur complement by an explicit elimination/`Solve` of the core linear system.
2. The `Schur form identity` (`deltaD - targetD == 0`) and the `low-frequency normalized outlet identity` (`deltaZ - targetZ == 0`) both still pass, with the SAME `rho_c`/`sigma_c`/`kappa_c`/`gamma_c`/`r_c` identifications emitted.
3. The route uses native Mathematica primitives in a way that is NOT a transliteration of the `.py` (the verifier checks independence: distinct primitive, distinct intermediate structure).
4. The committed Mathematica output remains byte-identical in its deliverable lines (method-only change); the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl`
- summary: Replaced the matrix-inverse Schur derivation with an explicit scalar elimination of the two core equations before forming `delta_Lambda(D)`.
- deviation: none
