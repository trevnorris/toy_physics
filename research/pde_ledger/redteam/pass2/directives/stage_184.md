---
unit_id: 184
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-08T22:45:50-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 184

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").
If a change is genuinely ambiguous or unsafe, append `## Blocked: F<n>` with a specific
question and continue with the rest.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. After editing, RUN the
affected scripts (`python3 <py>`, `math -script <wl>`, under `timeout 600`) and iterate until
they exit 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate,
never raise the cap. The orchestrator independently re-runs afterward.

## F1 — stale_output

**Target:** the two committed `.txt` (scripts/output/ and mathematica/output/).

**Issue:** Both predate the scripts and line 3 prints the stale banner `STAGE 167`; the
current scripts print `STAGE 184`.

**Required change:** None in source. The F2/F3 re-runs regenerate fresh transcripts; the
orchestrator refreshes both committed `.txt`. Confirm line 3 reads
`STAGE 184 — EXACT BRANCH-INVARIANT COORDINATES`.

## Applied: F1

- files_changed:
  - `scripts/output/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.txt`
- summary: Regenerated both transcripts and confirmed line 3 prints the Stage 184 banner.
- deviation: none

## F2 — tautological_check (both engines)

**Target:**
- `scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py:59-62`
- `mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl:47-48`

**Issue:** `Rtarget` is defined as `Lam0·(1-eps_eta_var)/T2` (py:59 / wl:47), then the very
next "product identity" check subtracts `Lam0·(1-eps_eta_var)` from `Rtarget·T2` (py:62 /
wl:48). By construction `Rtarget·T2 ≡ Lam0·(1-eps_eta_var)`, so the residual is an
identical-term cancellation that can never be nonzero — it confirms a definition against
itself, not the physics.

**Required change (both engines):** Replace the self-referential product check with the
falsifiable *drift* relation the notes assert (notes §1.3): `δln(R_target·T²) = δln(1-ε_η)`.
To make it falsifiable, `R_target` must enter as an INDEPENDENT perturbed object (e.g. its own
`Rtarget0·exp(small·R1)` with a free `R1`), NOT as `Lam0·(1-eps_eta_var)/T2` — so the residual
is a genuine function of `R1`, `Ξ₁`, `Σ_η`, `ε_η` that can fail. Preserve the downstream
selected-branch complement check (`δln[(R_target·T²)/Λ0] = -(ε_η/(1-ε_η))·Σ_η`) by routing it
through a *separately-named* closed-form identity object `(1-ε_η,var)` (the notes' closed form
for `(R_target·T²)/Λ0`), not through the now-independent `R_target`. The deliverable values are
unchanged; only the check is made non-tautological.

**Acceptance:** the printed residual for the new product-drift check is a symbolic expression
in `R1`/`Ξ₁`/`Σ_η` before simplification (not a literal `0` from identical-term cancellation);
the complement law still verifies; both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage184_branch_invariant_coordinates_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl`
- summary: Replaced the self-referential product identity with an independent `R_target` drift residual and routed the complement check through the closed-form `1 - eps_eta` branch object.
- deviation: none

## F3 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl`. NOT a paper_misalignment.

**Issue:** Every logarithmic drift is extracted with `SeriesCoefficient[Log[ratio],{small,0,1}]`
— the exact transcription of SymPy's `series(log(ratio),small,0,2).coeff(small,1)` (py:67/74/79/85)
— on the same exponential-ansatz objects (wl:44-47 ≡ py:56-59), with 1:1 assertion order/names.
The second engine is a port.

**Requirement (you design the route):** Re-author the `.wl` so it verifies the SAME
deliverables by a derivation path the `.py` does NOT use — genuinely independent, can-fail,
not a re-typing of the `Series`-coefficient choreography with a different API. You choose the
route. Preserve EXACTLY:
- D1: the branch product identity (now the F2 falsifiable drift relation, mirrored here);
- D2: `δln T_* = Σ_tr` with `T_* = R_tr^(-C_*)`;
- D3: `δln N_* = Σ_nt` with `N_* = T²·R_tr^(B_*)`;
- D4: `δln ε_η = Σ_η` and the selected-branch complement law.

**Acceptance:** the verifier confirms (a) `SeriesCoefficient[Log[…],{small,0,1}]` is no longer
the drift extractor mirroring the `.py` — an independent construction is visibly present;
(b) all deliverables still verify (residuals → 0, banner `STAGE 184`); (c) `math -script`
exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage184_branch_invariant_coordinates_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit to compute deliverable drifts by first variations on branch variables instead of `SeriesCoefficient[Log[...]]`.
- deviation: none
