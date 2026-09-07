# S11c-c2 N6 (representation invariance) — RESOLVED, per-engine SymPy (2026-09-06)

**Verdict (user-adopted): per-engine N6 PASSES under Reading B — COMPLETE OPERATOR COVARIANCE** — at retained order
`(η^{≤1}, σ_W^{≤1})`, under the declared-map premise, with the PIT qualification. ⛔ NOT strict equality of the raw
operands. This adopts the independent physics judgment of **Codex-astra** (`gpt-6-astra`, on explicit user OK for one
non-coding physics question; prompt `directives/_legs/S11c_c2_N6_interpretation_astra.md`, log
`scratchpad/astra_N6_interpretation.log`), which converged with Grok and reconciled Codex-sol's strict objection.

## What the object is (why covariance, not strict invariance)
The c2 self-energy increment is a **correction to an interface response operator** — a coupling block mapping
perturbation amplitudes → interface equation/force contributions (`§5c` object def) — NOT a scalar energy observable.
Replacing an Eulerian density perturbation by a material one adds the background advection `a_ρ=u·∇ρ/ρ`, which changes
how a given physical disturbance is DECOMPOSED into (density, thickness `e_W`, displacement `θ`) amplitudes. Its
coefficients therefore legitimately change with the field-variable frame; the same physical interface response
remains. Calling it "self-energy" does not make its operator components invariant numbers ⇒ it is expected to
transform covariantly (like a vector's components under a change of basis).

## The evidence chain (per-engine SymPy; instruments astra-built, all leg-cleared)
1. **Geometry reconciles.** Carrier `C_E = C_M` (no-nonzero, all 4 cases; live control). No representation dependence
   in the mechanical/geometric part. (`_measurements/S11c_c2_N6_reconcile_adjudication.md`.)
2. **The residual is purely constitutive.** `R_N6 = I_E − I_{M→E}` nonzero (~18 cols, 3 of 4 cases) localizes ENTIRELY
   to the source channel `B(C_M, ΔS)` (SPLIT_CHECK=0, carrier/cross=0) — i.e. to `μ`, where the advection enters.
3. **The frame-change is implemented CORRECTLY (the decisive test).** `R_cov = ms − source_terms(μ_E.subs(Φ), V_E)`
   (the actual material source vs the prediction from transforming the Eulerian source by the DECLARED, independently
   prolonged map Φ — built without `material_pullback`) shows **no nonzero found in all 4 cases** (conditional
   δ≈2.6e-22), knives bite (84, 4). ⇒ strong evidence the material construction faithfully implements Φ; the nonzero
   `R_N6` is (to that bound) the Φ-image of the
   source. (`_measurements/S11c_c2_N6_covariance_build_clearance.md`.)

⭐ **The reconciliation of the two readings (astra):** a comparison pushed all the way to genuinely COMMON variables
must still vanish. `R_N6 ≠ 0` is acceptable ONLY because `R_cov` shows no nonzero found (conditional δ≈2.6e-22) —
strong evidence that difference is the transformation content itself, not an additional physical response. Reading A (Codex-sol) is wrong as a demand for identical
coefficients across DIFFERENT field definitions, but its kernel — equality in common variables — is mandatory and IS
met. Reading B is the correct notion for this object.

## ⚠ CARRY-FORWARD caveats (astra) — Reading B does NOT close these; they are PREMISE checks, not defects
`R_cov` no-nonzero (conditional δ≈2.6e-22) is strong evidence the material builder IMPLEMENTS the declared Φ; it cannot exclude an error SHARED by the declared
premise AND both routes. Still owed (route to the WL engine / comparator / step record, ⛔ not pre-cleared):
1. **Is Φ itself physically correct?** Derive Φ (signs, thickness anchoring `h_α`) from the actual material motion +
   density/measure transformation — not merely confirm the builder reproduces the declared Φ.
2. **Does the face velocity `V` transform correctly?** `V_E ≡ V_M` (SHA-equal) establishes the two BUILDERS agree, not
   that both are the correct physical velocity (the prediction uses `V_E`, not `Φ(V_E)`; null here only because
   `V_E≡V_M`). Derive the velocity transformation.
3. **Extracted-block leakage.** The increment is an EXTRACTED coupling block; confirm the transformation does not
   import contributions from OMITTED blocks at retained order.

## Downstream (astra's guidance)
- **Blind Wolfram N6** must INDEPENDENTLY reproduce BOTH the carrier and the source-covariance checks, with
  able-to-fail controls (⛔ not import the SymPy result).
- **c2 comparator/reconcile:** demand cross-engine agreement in MATCHED representations; reconcile via the
  independently-specified transformation. ⛔ Neither engine may subtract its own observed discrepancy and call it
  covariance.
- **c2 step record:** preserve BOTH findings (raw `R_N6` nonzero in 3 cases; transformed `R_cov` no-nonzero), carry
  the 3 caveats, and explicitly resolve the **misleading `I_{M→E}` "mapped-operand" terminology** (the label is not
  the fully-mapped-to-common-variables operator). **Cross-engine agreement remains OWED.**

⇒ Per-engine (SymPy) N6 CLOSED as covariance-satisfied. NEXT = blind Wolfram engine N6.

## Owed / durable-evidence
- The ~499 MB reviewed `.out` (and the covariance/reconcile per-case `.out`, ~13-17 MB each) are ephemeral in `/tmp`
  — reproducible from the committed instruments; preserve or regenerate if needed for the step record.
- F/G numeric re-grounding PAUSED INDEFINITELY (user 2026-09-06: effort too large, not blocking physics).
