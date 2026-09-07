# S11c-c2 N6 covariance (R_cov) directive — decision-review adjudication (2026-09-06)

Orchestrator-written physics-bearing pre-builder directive → 2 decision legs (Codex-sol xhigh + Grok, identical
prompt). Both EXIT=0. Reports `scratchpad/{codex,grok}_N6cov_decision.log`. **Both: FOLD REQUIRED.** 3 verified folds
(F1 caught by BOTH, independently grounded); folded once → re-review (physics-bearing).

## F1 (BOTH legs) — Φ must be prolonged to the jet order present in `μ_E`, not the energy-level 0+1 dict
VERIFIED. `μ_E = EL(E) = ∂E/∂θ − Σ_i D_i(∂E/∂(∂_iθ))` (diagnostic `constitutive` :345-347); the `D_i` of a first-jet
flux writes **second** jets (`theta_second`/`theta_didi`, `e_W_didi`) into `μ_E` (`DERIVATIVE_MAP` b:733; Codex
confirmed the imported composite actually contains `theta_didi`/`e_W_didi` at `S11c_b_exports.py:5629`). The frozen
`Φ` dict (`frozen_relations` reconcile :53-57 / `material_pullback` b:1946-1968) is the **energy-level** map (fields +
first jets only) — complete for the scalar density `E`, INCOMPLETE for `μ`. Reusing it in `μ_E.subs(Φ)` leaves the
second jets unshifted, omits `D_iD_j a_ρ` / `D_iD_j h_α` (`O(σ_W)`, LIVE at retained order) ⇒ would **mislabel a
covariant construction as a defect**. **Fold:** `μ_E.subs(Φ)` = the differential prolongation of the field map through
EVERY θ/e_W jet in imported `μ_E` (currently rank-2: `θ_I ↦ D_I(θ+a_ρ)`, `e_{W,I} ↦ D_I(e_W+h_α)`), generated via the
live `b.total_derivative(...,background_depth=3)` + `DERIVATIVE_MAP` chain, applied `simultaneous=True`, with a
**pre-substitution domain-coverage census** proving every imported θ/e_W jet has a map entry. ⛔ Do NOT reuse the
energy-level 0+1 dict.

## F2a (Grok) — pin `R_COV_INCREMENT`'s `B` to `closed_response`, ⛔ not `build_increment`/`I`
VERIFIED (the same affine trap the reconcile forbade). `I(C_M,R_cov) = −C_M·p + B(C_M,R_cov)`; if astra implements the
guard as `build_increment(C_M,R_cov)`, then at `R_cov=0` it emits `−C_M·p` (certified NONZERO) — reporting a
discrepancy that is not there. **Fold:** `R_COV_INCREMENT` = reconcile `closed_response` on `(m_coeff, R_cov)`
(signatures 6/9/12 only), same keying as `SOURCE_CHANNEL`. ⛔ NOT `build_increment`/`I`.

## F2b (Codex) — the θ-independent-junk knife must be concretely executable
VERIFIED. `wave_terms` requires every nonzero term to carry exactly one recognized wave (diagnostic :363); a literal
θ-independent constant errors, or an unfortunate term is projected away — leaving that route-independence control
vacuous. **Fold:** prescribe a concrete actual-only perturbation after `μ_M` and before `source_terms`, e.g.
`κ_j·J_μ·e_W` with fresh nonzero-sampled `J_μ`, `DIMENSION_SCHEMA[J_μ]=(-1,-2,1)`, baseline `κ_j=0`, corrupted
`κ_j=1`, `μ_E∘Φ` unchanged; emit the source-level control delta. (Wave-linear so `wave_terms` accepts it.)

## Held sound (both legs)
Non-circularity (imported `μ_E` + supplied `Φ`, forbidding `material_pullback`/`ms`/`ΔS`; no hidden re-import); the
`a_ρ`/`h_α` formulas (route-2 §2); `Φ`-as-substitution; the object `R_cov = ms − ms_pred` at `source_terms`
(carrier-class, Z stripped); velocities isolate μ (V_E≡V_M SHA); the one-sided Φ-coefficient knife; no expected value /
no residual-zero exit; conditional PIT semantics; strict-vs-covariance withheld. ⚠ Note: the vet `.log` files are in
/tmp (Codex read-only couldn't reach them) — the committed adjudication record carried the needed content; for the
re-review, reference the committed adjudication, not /tmp logs.
