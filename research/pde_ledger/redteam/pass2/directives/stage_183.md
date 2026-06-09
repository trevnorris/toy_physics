---
unit_id: 183
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-08T22:36:47-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 183

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").
If a change is genuinely ambiguous or unsafe, append `## Blocked: F<n>` with a specific
question and continue with the rest.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. Do NOT change the
SymPy `.py`. After editing, RUN `math -script <wl>` (under `timeout 600`) and iterate until it
exits 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate, never
raise the cap. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl`. NOT a paper_misalignment.

**Issue:** The `.wl` is a line-by-line port of the `.py`: identical `eps` definition
(wl:38 ↔ py:46); identical hand-coded Stage-182 `theta1`/`xi1`/`r1` (wl:39-50 ↔ py:57-71)
with `sigmaTr` carried as an opaque symbol in both; identical `sigmaNTDef` via the same
term-deletion (wl:53-58 ↔ py:74-81); the SAME nine assertions in the SAME order with the SAME
label strings (wl:63-92 ↔ py:88-122). The `.wl` cannot catch a transcription error in the
shared hand-coded Stage-182 inputs.

**Requirement (you design the route):** Re-author the Mathematica audit so it verifies the
SAME deliverables by a route the SymPy script does NOT use — genuinely independent, can-fail.
You choose the route (a purely-numeric spot-check on a single substitution is NOT sufficient
on its own — the independent route must exercise the symbolic structure). Do NOT change the
`.py`. Preserve EXACTLY:
- D1: the branch-adapted `Σ_nt` definition;
- D2: `Ξ₁ = A_tr·Σ_tr + Σ_nt` and `R₁ + Ξ₁ + (ε_η/(1-ε_η))·Σ_η = 0` (triangular ledger);
- D3: the three exact inverse reconstructions (`Σ_tr`, `Σ_nt`, `Σ_eta` from the observables);
- D4: `A_tr/C_tr = 2(1+χ₀+δ_U)/δ_U`;
- D5: triple-rigidity — each diagonal prefactor (`C_tr`, `A_tr`, `ε_η/(1-ε_η)`) nonzero on
  the branch `χ₀,δ_U>0, 0<ε_η<1`.
Note the inverse reconstructions in the `.py` are round-trips by construction (`Σ_tr` →
`Θ₁` → `Σ_tr`); the independent `.wl` route should instead establish the map and its inverse
in a way that genuinely tests the prefactors (so a wrong `A_tr`/`C_tr`/`Σ_nt` coefficient
fails).

**Acceptance:** the verifier confirms (a) the `.wl` no longer mirrors the `.py`'s hand-coded
`theta1`/`xi1`/`r1` + nine-assert sequence as its verification mechanism — an independent
derivation is visibly present; (b) all deliverables still verify; (c) `math -script` exits 0;
(d) `.py` unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit around raw-slope row extraction, matrix factorization, symbolic inverse solving, and branch zero-set checks instead of the previous hand-coded observable sequence.
- deviation: none

## F2 — stale_output

**Target:** the two committed `.txt` (scripts/output/ and mathematica/output/).

**Issue:** Both predate the scripts and line 3 prints the stale banner `STAGE 166`; the
current scripts print `STAGE 183`.

**Required change:** None in source. Your required re-run of the re-authored `.wl` (F1)
regenerates the Mathematica transcript; the orchestrator refreshes both committed `.txt`.
Confirm line 3 reads `STAGE 183 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT`.

## Applied: F2

- files_changed:
  - `mathematica/output/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.txt`
  - `scripts/output/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.txt`
- summary: Regenerated both saved transcripts and confirmed line 3 carries the Stage 183 triangular-normal-form banner.
- deviation: none
