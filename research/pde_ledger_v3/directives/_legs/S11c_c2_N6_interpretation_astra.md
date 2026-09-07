# Physics-interpretation question (reasoning only, ⛔ no code) — what does per-engine N6 assert, and does it pass?

You are `gpt-6-astra`, asked here for an **independent physics judgment**, NOT a coding task. ⛔ Do not write or run
code; ⛔ do not modify the tree. You MAY read the cited files to ground yourself. Working dir `/var/projects/toy_physics`;
paths under `research/pde_ledger_v3/`. Reason from the physics; the two candidate readings below are presented
neutrally — ⛔ do not assume either is "the house answer."

## The physical system (conceptual)
S11c-c is the curved two-face **interface** (the throat / brane walls) sitting in a superfluid medium. "c2" is the
interface's **self-energy**: a perturbation of the interface radiates a disturbance into the bulk, it scatters and
returns, and it dresses the interface's own response (a self-energy in the field-theory sense). The perturbation is
carried by fields — a phase/displacement `θ`, a thickness `e_W` — over a background density `ρ` and flow `u`.

The same physics can be written in **Eulerian** coordinates (fields at fixed lab positions) or **material**
coordinates (co-moving with the medium — the Lagrangian view). The map relating them is the field redefinition `Φ`,
and the distinguishing piece is the **advection** `a_ρ = u·∇ρ/ρ` (plus a thickness shift `h_α = u·∇W_bg/W_bg` at one
anchoring) — the Eulerian-vs-Lagrangian bookkeeping difference. `N6` (representation invariance) asks whether the
interface self-energy computed the two ways is consistent.

## The established computed facts (per-engine SymPy, retained order (η^{≤1}, σ_W^{≤1}); finite-field PIT, δ≈2.6e-22)
Grounding: `_measurements/S11c_c2_N6_reconcile_adjudication.md`, `_measurements/S11c_c2_N6_sufficient_test_vet_adjudication.md`;
spec §5c `directives/S11c_c2_SHARED_PHYSICS.md`; parent N4 `directives/S11c_a_SHARED_PHYSICS.md` §5a +
`directives/S11c_decisions.md`.
1. **The geometry reconciles exactly.** The mechanical/geometric "carrier" of the response is identical in the two
   frames (`C_E = C_M`, a build-verified live control). No representation dependence in the geometry.
2. **The constitutive part differs by exactly the frame-change.** The residual `R_N6 = I_E − I_{M→E}` is nonzero
   (~18 columns) and localizes ENTIRELY to the constitutive **source** channel (`R_N6 = B(C_M, ΔS)`), i.e. to `μ`
   (the medium's energy response), which is where the advection enters.
3. **The frame-change is implemented CORRECTLY (the decisive new result).** The source-naturality residual
   `R_cov = ms − source_terms(μ_E.subs(Φ), V_E)` — the actual material source vs the prediction obtained by
   transforming the Eulerian source by the *declared* map `Φ` (built independently of the material builder, prolonged
   through every jet) — is **zero in all four cases** (no nonzero found at δ≈2.6e-22), while two able-to-fail knives
   (a wrong `Φ` coefficient `2·a_ρ`; θ-independent junk in `μ`) both move it (84 and 4 columns). ⇒ the material
   construction faithfully implements the Eulerian↔material transformation `Φ`; `R_N6` is exactly the `Φ`-image of the
   source, not a defect or a genuine non-covariance.

## The interpretation question (the split — decide it on the physics)
- **Reading A — strict invariance.** N6 means the self-energy increment is the SAME object (same value) in both
  frames: `R_N6 = 0`. Here `R_N6 ≠ 0`, so the increment is NOT strictly frame-invariant; it changes by the (correct)
  covariant `Φ`-channel. Under this reading per-engine N6 is "not strictly satisfied," though `R_cov=0` certifies the
  change is the correct/sanctioned one.
- **Reading B — covariance.** N6 means the two frames agree AFTER the declared field redefinition: `I_M = I_E ∘ Φ`.
  `R_cov = 0` closes this ⇒ per-engine N6 SATISFIED; the nonzero `R_N6` is the sanctioned `Φ`-transformation content.
  (Spec support: `S11c_decisions.md` N4 says the routes "must agree after that field redefinition"; §5c says "the same
  operator in two representations.")

## What I actually need from you
1. **Physically, what KIND of object is this self-energy increment** — a scalar invariant that ought to be numerically
   identical in every field-variable frame, or a component/density that is EXPECTED to transform covariantly under a
   change of field variables? Ground the answer in what the increment IS (a correction to an interface operator built
   from `μ` and the face geometry), not in the spec wording alone.
2. **Given (1) and the facts, which reading (A or B) is the physically correct notion of N6 for THIS object, and does
   per-engine N6 pass?** If neither is quite right, state the correct notion.
3. **Is there a physical hazard in Reading B** — e.g. could "covariant under the declared `Φ`" hide a real problem
   (a wrong choice of `Φ` itself, an object that SHOULD be invariant being allowed to transform, a missing channel
   such as the velocity `V` not being `Φ`-transformed while `V_E ≡ V_M` was only established, not derived)? Say what
   would still need checking, if anything.
4. **Downstream:** does the choice change what the blind **Wolfram** cross-engine N6, the c2 comparator/reconcile, or
   the c2 step record must assert or test? (Reason about it; ⛔ do not build anything.)

## Output
A direct conceptual answer to 1-4, with your reasoning grounded in the physics of the object. End with your **verdict**:
Reading A, Reading B, or a corrected notion — and whether per-engine N6 passes. Brief, physics-first.
