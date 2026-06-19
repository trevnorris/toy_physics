# Directive pathA_20b — Emergent constants, Step 2b: derive the gauge cone `c_γ` vs `c_s` (bulk + brane)

**Status:** Design-reviewed (Codex `gpt-5.5` → SOUND-WITH-FIXES; all fixes applied 2026-06-19; pending a confirm-pass
before execution — gated by the user). Resolves/sharpens the `C_GAMMA_RATIO_UNDERDETERMINED` residual that `pathA_20`
carried. Consumes the `pathA_19` base + `pathA_20` outputs: `c_s²=5Kρ0⁴/m_GNLS`, `[c_s]=[c_γ]=L T⁻¹`, the carried
`λ_γ=c_γ/c_s` and tail `(c/c_s)³=λ_γ³`.
**Date:** 2026-06-19
**Owner:** Codex (DERIVES + codes; iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** user fork choice B (decision-13 §0/§9, 2026-06-19). Chain: `pathA_19` → `pathA_20` → **this (`c_γ/c_s`)** →
`pathA_21` (`G` + mass-bridge) → `pathA_22` (scale-map → `m̂0²·S_port` → re-run B2c).

## Why this step (context)
`pathA_20` established (adversarially verified) that `Z(w)` is an OVERALL prefactor on `F_{MN}F^{MN}` — it cancels from
the gauge principal symbol and renormalizes only the coupling (`μ0_eff=μ0/Z_int`, `q_eff`), NOT the cone. So the gauge
field's characteristic speed is set by the **time-vs-space coefficient ratio** of the Maxwell principal symbol (the
parent metric normalization), while `c_s` is the **emergent acoustic cone** (`c_s²=5Kρ0⁴/m_GNLS`). The parent action
never relates them; `em_fields.tex` only ASSERTS `c=c_s` weak-field by reusing the acoustic d'Alembertian by fiat.
Whether `c_γ=c_s` is the linchpin of the "particle = standing photon wave" ontology (ceiling `c=c_γ`; if `c_γ≠c_s` the
photon is NOT sound). Derive it from the coupled action — separating the BULK principal cone (tractable here) from the
OBSERVED BRANE cone (may need the profile).

## Scope & stance
DERIVATION via linearization on the **HOMOGENEOUS background** (uniform `ρ0`, `v_b=0`, asymptotic). **Scope of what the
homogeneous calc can settle:** the BULK Maxwell principal cone, and — ONLY if a controlled far-field zero-mode reduction
is explicitly established — the asymptotic brane photon cone. If the observed brane photon requires mixed-sector
localization, projection data, or the solved throat profile, the final `c_γ/c_s` handoff remains
`C_GAMMA_RATIO_STILL_UNDERDETERMINED` with a named brane-cone residual. No model-formula changes, no freeze touch, no
`m̂0²·S_port`, no `M`-collapse, no `α_J`/`G`. Extend `dimensional_check.py` side-by-side (new group; leave pathA_18/19/20
groups + `dimensional_dictionary()`/D1-D3 intact). Same infra constraints as `pathA_20` (read-only sim dir; never touch
`physical_export_permitted`; no `$RT exec`; `timeout 600`; ≤2 MMA seats; YAML/md). Mass-symbol discipline: `m_GNLS`.
**Dual-engine** required for the dispersion/principal-symbol ALGEBRA only; metric-inheritance, zero-mode-reduction
validity, and "which branch is the physical photon" are non-algebraic judgments → human-readable residuals unless reduced
to an explicit quadratic operator. **Do NOT assume `c_γ=c_s`. Shared dimensions, natural-unit pins, `em_fields.tex`, and
reused acoustic-d'Alembertian notation are explicit FAILs as proof of equality.**

## Work items

### L1 — background validity + the coupled linearization setup
- **Background-validity check FIRST:** verify the homogeneous background actually SOLVES the coupled equations — specify
  `ψ0`, `A_{M0}`, the matter current `J_ψ0`, any external `J_ext0`, and the neutrality / current-cancellation condition
  (the linearized Maxwell source includes `δJ⁰=q_★ δρ` and a `ρ0 δA^i` term; `pde.tex:370-374, 912-925`). If a legal
  uniform charged background cannot be established from the cited sources, STOP with `HOMOGENEOUS_CHARGE_NEUTRALITY_UNSPECIFIED`.
- Linearize the COUPLED GNLS + gauge (Maxwell-sector) action (`part01_parent_geometry.tex:225-247`; `pde.tex:357-416`)
  about that background. Perturbations `(δρ, δθ)` (phonon) and `δA_M` (gauge).
- **`c_bulk` definition:** restore the time-derivative (`C_E`) and spatial-gradient (`C_B`) coefficients of the Maxwell
  PRINCIPAL symbol, including the parent metric temporal/spatial normalization; define `c_bulk² = C_B/C_E`. **Overall
  `Z/μ0` factors do NOT count as speed-setting data.** If the parent action only gives `F_{MN}F^{MN}` in natural units
  (`η_MN`, implicit `c=1`), then `c_bulk` carries no source-derived value relative to `c_s` → residual
  `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED` / `PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING`.
- The photon-cone claim must be made from GAUGE-INVARIANT field strengths / transverse physical modes; gauge-fixing
  terms may be consistency-checked but cannot be the source of the speed verdict.

### L1b — the coupled principal symbol (do NOT derive the two sectors in isolation)
- Form the coupled principal symbol for `(δρ, δθ, δA_M)`. EITHER show the off-diagonal principal GNLS↔gauge couplings
  vanish on the declared background, OR compute the coupled characteristic determinant and identify its branches.
- Separate PRINCIPAL-symbol terms (set the characteristic cone) from LOWER-ORDER current/London/plasma/mass terms
  (`pde.tex:843-852, 903-929`): the latter can gap/mix modes WITHOUT changing the high-frequency cone — record them
  separately and do NOT use them as the characteristic cone.

### L2 — the two dispersion relations
- PHONON dispersion (GNLS sector) → confirm the acoustic cone `c_s=√(5Kρ0⁴/m_GNLS)`; machine-check `[c_s]=L T⁻¹`.
- GAUGE dispersion → `c_γ` = the characteristic speed of the MASSLESS/TRANSVERSE gauge branch; machine-check
  `[c_γ]=L T⁻¹`. If no massless transverse branch is established, record `MASSLESS_GAUGE_BRANCH_NOT_ESTABLISHED`.
  Label any plasma/gapped/longitudinal branches separately.

### L3 — the TWO-LAYER verdict
- **`bulk_verdict ∈ { C_GAMMA_BULK_EQUALS_C_S | C_GAMMA_BULK_NE_C_S (closed-form `c_bulk/c_s` + ρ-dependence) |
  C_GAMMA_BULK_UNDERDETERMINED (named residual) }`** — from the homogeneous principal symbol. `EQUALS` requires a
  SOURCE-DERIVED equation identifying the gauge kinetic metric with the acoustic metric (`c_bulk²=5Kρ0⁴/m_GNLS`), NOT a
  shared symbol / natural-unit `1=1`.
- **`brane_verdict ∈ { C_GAMMA_EQUALS_C_S | C_GAMMA_NE_C_S (closed-form `λ_γ` + ρ-dependence) |
  C_GAMMA_RATIO_STILL_UNDERDETERMINED (named residual) }`** — the OBSERVED photon cone after the (explicitly declared)
  far-field zero-mode reduction; if that reduction is not established or needs the profile, this is
  `C_GAMMA_RATIO_STILL_UNDERDETERMINED` with `BRANE_ZERO_MODE_REDUCTION_UNDERIVED` or `BRANE_PHOTON_CONE_REQUIRES_PROFILE`.
- **pathA_21 consumes ONLY `brane_verdict`.** Also report whether `λ_γ` is a pure number or drifts with `ρ0` (since
  `c_s∝ρ0²`, an independent `c_bulk` ⟹ `ρ0`-dependent ratio).

### L4 — implications recorded (NOT derived here)
- Consequence for the standing-wave ontology (ceiling `c=c_γ` vs `c_s`) and the `(c/c_s)³` tail factor (`R_tail=λ_γ³`).
- The brane-localization sub-question explicitly belongs here as a possible BLOCKING brane-cone residual (not an
  afterthought), per the scope statement.

## Acceptance criteria (PASS/FAIL with NAMED RESIDUALS; script exit-0 is NECESSARY, NOT SUFFICIENT)
1. Background validity established (or `HOMOGENEOUS_CHARGE_NEUTRALITY_UNSPECIFIED`); coupled linearization set up with the
   gauge PRINCIPAL coefficients (`C_E`, `C_B`) restored as-written; photon cone from gauge-invariant/transverse modes.
2. The coupled principal symbol is formed (off-diagonal vanishing shown, or coupled determinant computed); lower-order
   plasma/London/mass terms separated from the characteristic cone.
3. Both dispersions derived + machine-checked (`[c_s]=[c_γ]=L T⁻¹`); the phonon branch reproduces `c_s²=5Kρ0⁴/m_GNLS`.
4. A TWO-LAYER verdict: `bulk_verdict` AND `brane_verdict` (each from its taxonomy), with `λ_γ` ρ-dependence stated.
   `pathA_21` consumes only `brane_verdict`.
5. **Negative control:** a harness fixture keeps `c_bulk` and `c_s0` as INDEPENDENT symbols; any `EQUALS` verdict FAILS
   unless a source-derived equation sets `c_bulk²=5Kρ0⁴/m_GNLS`. Machine checks of `[c_γ]=[c_s]=L/T` are explicitly
   NON-evidentiary for equality.
6. Dual-engine `.wl` agreement for the algebra; new harness group passes; scripts exit 0 within `timeout 600`;
   pathA_18/19/20 groups untouched.

**Allowed sharpened residuals (each needs a source anchor + downstream consequence):**
`BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`, `PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING`,
`HOMOGENEOUS_CHARGE_NEUTRALITY_UNSPECIFIED`, `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`, `BRANE_PHOTON_CONE_REQUIRES_PROFILE`,
`MASSLESS_GAUGE_BRANCH_NOT_ESTABLISHED`. A `STILL_UNDERDETERMINED` verdict MUST name one of these (sharper than
pathA_20's generic wall), not re-punt.

## Out of scope
The throat/brane-localization renormalization IF it needs the solved profile (carry as a blocking brane-cone residual);
`G` + the mass-bridge + `α_J` (`pathA_21`); the scale-map → `m̂0²·S_port` (`pathA_22`); B2c rerun; freeze amendment. Do
NOT collapse `M` or change the base set.

## Review (orchestrator, after Codex)
Transliteration-fidelity (one clean agent per new script); independent re-derivation of the phonon + gauge principal
symbols, `c_bulk`, and `λ_γ`; adversarial pass (distrust-all-clean AND distrust-premature-underdetermined) — did Codex
genuinely form the COUPLED principal symbol and derive both cones, or punt / derive sectors in isolation? Is any `EQUALS`
verdict backed by a source-derived `c_bulk²=5Kρ0⁴/m_GNLS` (not a shared symbol / em_fields fiat)? Does the negative
control actually fail a forced equality? Claude reads only residuals. Then gate to `pathA_21`.
