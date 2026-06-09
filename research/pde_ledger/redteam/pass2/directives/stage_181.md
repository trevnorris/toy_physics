---
unit_id: 181
batch: V.2
created_at: 2026-06-08T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-08T22:36:53-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 181

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under
that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").
If a change is genuinely ambiguous or unsafe, append `## Blocked: F<n>` with a specific
question and continue with the rest.

Do NOT touch paper.tex, notes/, or any prose document — scripts only. After editing, RUN the
affected scripts (`python3 <py>`, `math -script <wl>`, under `timeout 600`) and iterate until
they exit 0 with every in-file check PASS. A timeout (exit 124) is a FAILURE — reformulate,
never raise the cap. The orchestrator independently re-runs afterward.

## F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Target:** `mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`

**Status:** The user has AUTHORIZED a full re-author of this `.wl`. NOT a paper_misalignment —
the math is correct; the second engine is just not independent.

**Issue:** The `.wl` is a line-for-line port of the `.py`: identical symbol names and object
order; the identical hand-transcribed closed-form targets `eps1Expected` (wl:75-78) and
`theta1Expected` (wl:114-118) compared against `D[…]·coeff` constructions that mirror the
SymPy ones (py:76-80, py:128-132); the identical `s`-perturbation slope route
(`D[Log[…],s]/.s->0` ↔ `sp.diff(sp.log(…),s).subs(s,0)`); the identical `bad`-spoiler; the
same assertion order/labels. A transcription error baked identically into both targets would
pass both engines.

**Requirement (you design the route):** Re-author the Mathematica audit so it verifies the
SAME deliverables by a derivation path the SymPy script does NOT use — genuinely independent,
can-fail, not a re-typing with a different API. You choose the route. Do NOT change the `.py`.
Preserve EXACTLY:
- D1: the direct↔selected transfer-shape identity;
- D2: support-blindness — `d/dζ ln T² = 0` and `d/dζ ln R_target = 0` on the support-loaded
  route, AND the spoiled-packet counter-check (a deliberately wrong support packet must break
  blindness);
- D3: the split-blocking drift `ε₁`, the defect law `Ξ₁`, the selected-branch identity
  `Ξ₁ + η₁/(1-ε_η) + R₁ = 0`, and the tracking-factor drift `Θ₁`.
The load-bearing closed forms (`ε₁`, `Θ₁`) must be DERIVED by the independent route, not
re-typed as `eps1Expected`/`theta1Expected` literals.

**Acceptance:** the verifier confirms (a) the `eps1Expected`/`theta1Expected` hand-literals
and the `D[Log[…],s]/.s->0` mirror are no longer the verification mechanism — an independent
derivation is visibly present; (b) all deliverables still verify (residuals → 0); (c) the
spoiled-packet check still fails-on-purpose; (d) `math -script` exits 0; (e) `.py` unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit around directional differentials and product/factor ledgers so epsilon_1 and Theta_1 are derived rather than hand-literal expected forms.
- deviation: none

## F2 — tautological_check (both engines)

**Target:**
- `scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py:57`
- `mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl:49`

**Issue:** The "support-loaded branch product law" check
`R_target_loaded·M_tr − product_loaded == 0` (py:57) /
`rTargetLoaded·loadMassFromSupport − productLoaded == 0` (wl:49) is a round-trip of the
immediately-preceding definition `R_target_loaded := product_loaded/M_tr` (py:51) /
`rTargetLoaded := productLoaded/loadMassFromSupport` (wl:43). It reduces to
`product_loaded − product_loaded ≡ 0` and cannot fail.

**Required change:** Remove this single redundant check in BOTH engines (py:57 and wl:49).
The genuine support-blindness content (the `d/dζ` checks, the reconstruction identities, and
the spoiled-packet counter-check) is the load-bearing evidence and stays. (In the `.wl` this
removal is folded into the F1 re-author.)

**Acceptance:** the "support-loaded branch product law" line no longer appears in either
transcript; both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage181_coherent_tracking_defect_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage181_coherent_tracking_defect_mathematica_audit.wl`
- summary: Removed the redundant support-loaded branch product-law assertion from both engines while leaving the reconstruction, support-blindness, and spoiled-packet checks intact.
- deviation: none
