---
unit_id: 180
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-08T22:28:01-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 180

Apply each finding below. After applying, append an `## Applied: F<n>` block under that
finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").
If a finding's required change is genuinely ambiguous or unsafe, append `## Blocked: F<n>`
with a specific question instead — skip that finding, continue with the rest.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. Do NOT change the
SymPy `.py`. After editing, RUN `math -script <wl>` (under `timeout 600`) and iterate until
it exits 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate,
never raise the cap. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl` (the re-author-vs-accept
call was escalated and resolved in favour of re-author). This is NOT a paper_misalignment —
the math is correct and aligned; the defect is that the second engine is not independent.

**Issue:** The `.wl` is a line-by-line port of the `.py`. It mirrors, in the same order:
- `teff2 = t1^2 Exp[2 eps lam tau1] + t2^2 Exp[2 eps lam tau2]` then
  `xiEff = (D[Log[teff2],eps]/.eps->0)/lam` — identical object and identical slope operator
  to py:41-42 (`sp.diff(sp.log(Teff2),eps).subs(eps,0)/lam`);
- `t2Direct = beta0/k0` with the same hand-form `beta0`, `k0` (wl:43-45 ↔ py:61-65);
- the perturbed `t2Pert` ansatz with `(D[Log[t2Pert],e]/.e->0)/lam` (wl:71-74 ↔ py:100-105).
A transcription error in the shared hand-forms would pass both engines identically.

**Requirement (you design the route):** Re-author the Mathematica audit so it verifies the
SAME four deliverables by a derivation path the SymPy script does NOT use — a genuinely
independent route that could FAIL if the shared closed forms were mis-transcribed. It must
not be a re-typing of the `.py`'s construction with a different simplifier/extraction API.
You (Codex) choose the independent route. The deliverables to preserve EXACTLY:
- D1: multi-port effective-shape collapse `Ξ_eff = 2(ρ₁τ₁ + ρ₂τ₂)` with `ρ_r = T_r²/ΣT²`;
- D2: one-port continuum `T² = Z_W(1+ρ)²/[Ω_W²(1-ε_W)²]` (and the `β₀/K₀ → μ_W/K_W` form);
- D3: selected-branch `T² = (27π²G c_s⁵/(20 a⁵ c⁵))(1-ε_η)/R_target`;
- D4: the two weak-axisymmetric slope laws (direct and selected-branch).

**Acceptance:** the verifier confirms (a) the slope/collapse extractions are no longer the
`D[Log[…]]/.→0` (or `Series`-coefficient) mirror of the `.py` on the same hand-coded objects
— i.e. an independent construction is visibly present; (b) all four deliverables still verify
(residuals → 0, banner `STAGE 180`); (c) `math -script` exits 0; (d) the committed `.py` is
unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit to verify D1-D4 through port-load first variations, algebraic branch equation solving, and scaled differential slope checks rather than the SymPy perturb-and-log construction.
- deviation: none

## F2 — stale_output

**Target:** the two committed `.txt` (scripts/output/ and mathematica/output/).

**Issue:** Both predate the scripts and line 3 still prints the pre-renumber banner
`STAGE 163`. The current scripts print `STAGE 180`.

**Required change:** None in source. Your required re-run of the re-authored `.wl` (F1)
regenerates the Mathematica transcript with the `STAGE 180` banner; the orchestrator
refreshes both committed `.txt` from the fresh runs. Confirm after re-run that line 3 reads
`STAGE 180 — EFFECTIVE TRANSFER-SHAPE COLLAPSE`.

## Applied: F2

- files_changed:
  - `mathematica/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.txt`
  - `scripts/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.txt`
- summary: Regenerated both committed transcripts from fresh successful runs and confirmed line 3 reads `STAGE 180 — EFFECTIVE TRANSFER-SHAPE COLLAPSE`.
- deviation: none
