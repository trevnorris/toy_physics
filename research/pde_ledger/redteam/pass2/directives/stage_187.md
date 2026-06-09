---
unit_id: 187
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-09T00:13:48-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 187

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").
If a change is genuinely ambiguous or unsafe, append `## Blocked: F<n>` with a specific
question and continue with the rest.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. After editing, RUN the
affected scripts (`python3 <py>`, `math -script <wl>`, under `timeout 600`) and iterate until
they exit 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate,
never raise the cap. The orchestrator independently re-runs afterward.

## F1 — insufficient_verification (SymPy engine parity)

**Target:** `scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py:106`

**Issue:** The selected-minor determinant (the load-bearing non-singularity premise D2 from
notes §4: "the selected minor has determinant `1+χ₀,*>0`, this solve is exact and unique") is
only `print`ed in SymPy, never asserted. The Mathematica engine asserts it (wl:79). A
mis-transcription of `minor` in the `.py` would not be caught.

**Required change:** Immediately after the existing `print(... minor.det() ...)` at py:106, add
the matching assertion `expect_zero("selected minor determinant", minor.det() - (1 + chi))`,
where `chi` is the already-declared `chi0_star` (py:43). Change nothing else for F1.

**Acceptance:** a new `selected minor determinant = 0` line appears in the SymPy transcript and
the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage187_orbit_quotient_closure_sympy_audit.py`
- summary: Added the missing selected-minor determinant assertion immediately after the determinant print.
- deviation: none

## F2 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl`. NOT a paper_misalignment.

**Issue:** The `.wl` is a line-by-line port of the `.py`: the hand-coded finite log-ratio rows
`rowTr`/`rowNt`/`rowEta` (wl:34-36 ↔ py:46-48), the hand-coded monomial ratios (wl:48-50 ↔
py:70-74), the same `m`/`Dx` matrix, and the same `Solve[{rowTr==0,rowNt==0,rowEta==0},{…}]`
(wl:81 ↔ py:109). Both engines solve the same hand-coded system; a transcription error in the
shared rows/matrix would pass both.

**Requirement (you design the route):** Re-author the `.wl` so the central object (the finite
`M_*` rows / monomial log-ratio equalities) is DERIVED from the physical monomial definitions
rather than hand-coded identically to the `.py`, so a wrong exponent fails. You choose the
route. Do NOT change the `.py`. Preserve EXACTLY:
- D1: the three finite log-ratio equalities = `M_* Δx = 0` (from the `C_tr`/`C_nt`/`ε_η`
  monomial ratios);
- D2: the selected-minor determinant `= 1+χ₀,*`;
- D3: the exact finite fibre solve for `(Δ_η, Δ_T, Δ_μ)` reproducing the Stage-186 orbit laws;
- D4: substitution-back annihilation of the three invariant-fibre equations.

**Acceptance:** the verifier confirms (a) the finite rows / matrix are no longer a hand-coded
mirror of the `.py` solved by the same `Solve` — an independent derivation of the rows from the
monomial ratios is visibly present; (b) all deliverables still verify (banner `STAGE 187`);
(c) `math -script` exits 0; (d) `.py` unchanged.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage187_orbit_quotient_closure_mathematica_audit.wl`
- summary: Reauthored the Mathematica audit to derive finite rows and the matrix from physical monomial ratios, then solve the fibre equations by triangular elimination.
- deviation: none
