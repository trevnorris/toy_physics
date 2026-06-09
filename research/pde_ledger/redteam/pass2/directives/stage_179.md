---
unit_id: 179
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-08T22:19:56-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 179

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").
If a change is genuinely ambiguous or unsafe, append `## Blocked: F<n>` with a specific
question and continue with the rest.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. After editing, RUN the
affected scripts (`python3 <py>`, `math -script <wl>`, under `timeout 600`) and iterate until
they exit 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate,
never raise the cap. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl`. NOT a paper_misalignment.

**Issue:** The `.wl` mirrors the `.py`: the same hand-coded perturbation `n0A` (wl:51-60 ↔
py:61-70) with the slope `nuDirect = (D[Log[n0A],eps]/.eps->0)/lam` transcribing
`sp.series(sp.log(N0A),eps,0,2).coeff(eps,1)/lam` (py:73); the same hand-coded `T`/`tau`/
slippage-form constructions; the same checks in the same order.

**Requirement (you design the route):** Re-author the `.wl` so it verifies the SAME
deliverables by a derivation path the `.py` does NOT use — genuinely independent, can-fail.
You choose the route. Do NOT change the `.py`. Preserve EXACTLY:
- D1: the factorization `N0/K = T²` with `T = (Ĝ_W + R̂ Ĝ_U)/(1-R̂²)`;
- D2: the weak-axisymmetric slope identity `ν_r = κ₁ + 2τ_r`;
- D3: the equivalence to the slippage form `τ_r = m_r + I_r/(1+I_r) i_r + H_r/(1-H_r) h_r`;
- D4: the weighted defect identity `Ξ₁ = 2 Σ_r ρ_r τ_r` (strengthened per F2).

**Acceptance:** the verifier confirms (a) the `.wl` slope/identity extractions are no longer
the `D[Log[…]]/.→0` (or `Series`-coeff) mirror of the `.py` on the same hand-coded objects —
an independent construction is visibly present; (b) all deliverables still verify (banner
`STAGE 179`); (c) `math -script` exits 0; (d) `.py` unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit around wall-normalized shape variables and a logarithmic directional-slope operator instead of the primitive epsilon expansion.
- deviation: none

## F2 — insufficient_verification (both engines)

**Target:**
- `scripts/moving_throat_pde_stage179_transfer_shape_sympy_audit.py:105-112`
- `mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl:87-93`

**Issue:** The "Weighted defect identity" block builds BOTH `Xi` and `Xi_expected` from the
same freshly-declared free symbols `tau1,tau2,tau3` via `kappa1 + 2*tau_i`, so
`Xi - Xi_expected == 0` reduces to the `Σρ=1` algebra tautology. It does not connect the
per-port `ν_r` to the closed-form `τ` validated by the slope-identity check (py:84-91 /
wl:70-75). The substantive collapse `Ξ₁ = Σ_r ρ_r(ν_r-κ₁) = 2Σρτ` with `ν_r` the ACTUAL port
slope is not exercised.

**Required change (both engines):** Make the weighted-defect check load-bearing — feed the
per-port `ν_r` from the validated relation `ν_r = κ₁ + 2τ_r` (the closed-form `τ` established
by the slope-identity check), so the collapse rides on that validated identity rather than on
two copies of the same template. If the only honest strengthening requires reusing the
closed-form `τ` under independent per-port slope substitutions and that proves impractical,
append `## Blocked: F2` with that assessment rather than applying a still-tautological change.
(In the `.wl` this is folded into the F1 re-author.)

**Acceptance:** the `Ξ₁ - 2 weighted τ` check still reports `= 0` / PASS, and the weighted block
references per-port `ν_i = κ₁ + 2τ_i` rather than constructing both sides from the same template;
both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage179_transfer_shape_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl`
- summary: Strengthened the weighted-defect checks to build per-port `nu_i` and `tau_i` from validated port identities before collapsing to `Xi_1 = 2 weighted tau`.
- deviation: none

## F3 — stale_output

**Target:** the two committed `.txt` (scripts/output/ and mathematica/output/).

**Issue:** Both predate the scripts and line 3 prints the stale banner `STAGE 162` (and an
out-of-date `nu_direct`). The current scripts print `STAGE 179`.

**Required change:** None in source. The F1/F2 re-runs regenerate fresh transcripts; the
orchestrator refreshes both committed `.txt`. Confirm line 3 reads
`STAGE 179 — WALL-NORMALIZED TRANSFER-SHAPE THEOREM`.

## Applied: F3

- files_changed:
  - `scripts/output/moving_throat_pde_stage179_transfer_shape_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage179_transfer_shape_mathematica_audit.txt`
- summary: Regenerated both transcripts from the updated scripts and confirmed line 3 is the Stage 179 banner.
- deviation: none
