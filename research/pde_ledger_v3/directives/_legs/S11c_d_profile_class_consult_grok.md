# Grok (grok-4.6 high) — S11c-d profile-class design consult (recommendation)

Advisory design consult (NOT a review leg). Verbatim.

---

I'll read the S11c decision and physics context first, then give an independent recommendation on the profile class rather than echoing the docs.I have the family split and power counting. Next I’ll pull the S11c-a/S11b coupling structure and the slit-edge falsification notes so the recommendation is grounded in the physics, not the decision list.**Recommendation: localized profile, weak-contrast Born kernel, with `σ_W` kept live. The leaning is the right object and the wrong name.** Call it “weak-gradient Born” and you will quietly expand away the slit.

---

## 1. What S11c-d should fix

**Class:** localized — `∇W₀` (equivalently `∇μ_R`) of compact or rapidly decaying support. That includes a defect bump *and* a weak edge; both source mixing only where the gradient lives.

**Regime:** Born in the contrast `η`, **not** in the gradient. Leakage at retained order is the on-shell matrix element of `CLOSED_COUPLING_KERNEL` between the sector modes/resolvent of `CLOSED_SLAB_OPERATOR`. Coupling stays `O(εη)`, rate `O(ε²η²)`.

**Output object:** the linear mixing — a profile-conditioned Born scattering-and-leakage kernel — plus the two N13 channels as separate objects: continuum conversion, and capture into a bound thickness/breathing pole. Not a recipe (“do perturbation theory”), and not `ω(k)` for generic `W₀(x)`.

**Why this, not the other two named classes**

- **WKB is the wrong geometry for the endgame.** A slit edge is anti-adiabatic. WKB conversion is Landau–Zener: exponentially small in `1/σ_W` when gapped, or order-unity at a crossing. That is the FORM of a taper, not of an edge. Using it would compare the lab knife-edge to an adiabatic waveguide.
- **Bloch is a sibling evaluation, not the first object.** A grating is periodic, and “gratings work” is a periodic *witness*. At small `η` a Bloch gap is the same kernel’s Fourier coefficient at Bragg kinematics. You get that later by feeding a periodic `w₁` to the same Born object. Choosing Bloch first throws away N13’s isolated bound breathing mode (bands, not a well). In 1D a weak attractive well always binds; that photon-kill channel is real at small `η` and is not a grating bandgap.

Localized is the only class that matches both the defect-conversion object (a particle/defect on the brane) and the slit-edge observable, without inventing a global dispersion.

**Why Born in `η`, not order-unity now**

The frozen exports are first shape-order truncations. The background ansatz is `W_bg = W̄₀[1+η w₁(ξ)]` with independent first-jet scale `σ_W = η W̄₀/L_W`, truncated at first order in each of `ε, η, σ_W`. Those operators are the first variation of the uniform theory in the background. Setting `η → O(1)` inside them is not “going non-perturbative”; it is using an incomplete operator at the order where the missing `O(η²)` pieces compete. That is the bookkeeping constraint, and it is physics.

Both exports are used, as-is, with no new upstream derivation:

- diagonal `CLOSED_SLAB_OPERATOR` → sector modes, Green function, thickness-sector poles (N13);
- off-diagonal `CLOSED_COUPLING_KERNEL` → the mixing vertex.

At retained leakage order `O(η²)`, distorted-wave corrections (diagonal inhomogeneity dressing `G₀`) are higher order. Leading conversion is the uniform-mode matrix element of the closed coupling kernel. That is exactly the linear mixing the family exists to compute.

---

## 2. The Born vs slit-edge tension — address it, don’t paper over it

A lab slit edge is **localized + order-unity contrast + sharp on the optical scale**. The exports can honestly do the first of those, and can *shape*-see the third only through a live `σ_W` on a smooth profile. They cannot do the second.

Do **not** evaluate the Born coefficient at `η = 1` and compare that number to the lab. That is the invalid extrapolation N7 already named.

The bridge that does **not** invalidate the lab bound:

1. **S11c-d** emits the kernel with `η` a bookkeeping parameter, never set to 1 in the build. The profile class is localized; the shape `w₁` stays general inside that class (or is a named representative, not a frozen “the slit”).
2. **S11c-e** emits the flux-normalized dimensionless **FORM** of that kernel: `∝ k·a`, supported where `∇μ_R ≠ 0`, with the `σ_W` form factor (Fourier content of the localized gradient). Magnitude stays blocked on throat-interior physics.
3. **The order-unity slit is a separate downstream reduction, withheld from the builder.** Two inferences are legitimate; a third is not.
   - **Legitimate, and non-perturbative:** uniform regions cannot mix (S11b, exact). Mixing is supported only on `∇μ_R ≠ 0`. A slit edge is where conversion *can* occur. That is FORM, not a number.
   - **Legitimate as a withheld OOM reductio:** if the unknown interpolating function `F(η)` with `F(η) ∼ η²` as `η → 0` were `O(1)` conversion per edge at `η ∼ 1`, gratings would fail. The lab bounds `F(1)`, not the Born prefactor evaluated at 1.
   - **Illegitimate:** “S11c-d computed the slit” or “the lab bounds our small-`η` coefficient.” Those identify `F'(0)` with `F(1)`.

What weak-`η` Born actually drops, and should be named so they are not silently lost:

| Dropped | Why it can wait |
|---|---|
| Unitarity / saturation | Born can exceed 1; exact `T` cannot. If Born is already `≪ 1`, this is usually conservative. It can *miss* a finite-`η` resonance that converts `O(1)` in a narrow band — that is a later risk, not a reason to fake `η = O(1)` with a truncated operator. |
| Coefficient renormalization | Changes magnitude, not the `k·a` / `∇μ_R` FORM. Magnitude is R1-blocked anyway. |
| True step / all jets | A knife-edge is not first-jet. The exports do not contain `∇²W, ∇³W, …`. A smooth localized `w₁` with live `σ_W` is the sharpest object they can see. |
| Domain change | A hole in the brane (`W → 0`, terminus) is not a thickness perturbation on a fixed slab. That is new geometry. |

What it does **not** drop, despite the tension:

- The FORM. Decoupling at `η = 0` is exact. Support on `∇μ_R` is exact. Higher orders in `η` can renormalize the prefactor and open weaker channels; they do not make a uniform patch convert.
- The bound-thickness channel at *small* `η` (N13). You do not need order-unity contrast to see it.
- The edge form factor, **if `σ_W` is kept live.** The mixing density is `∝ ∇μ_R`. The gradient of a kink *is* a bump, so a weak edge and a defect enter the vertex the same way; the slit shape is a choice of `w₁`, not a different construction.

---

## 3. Staged / hybrid — yes, with one naming constraint

The right split is:

- **S11c-d:** localized, weak-`η` Born scattering-and-leakage kernel (`σ_W` live). Object, not recipe.
- **S11c-e:** flux-normalized conversion FORM of that kernel. Magnitude out of scope.
- **Order-unity slit:** orchestrator-side reduction / withheld OOM, not a builder target and not a second profile class inside S11c-d.
- **Periodic:** later evaluation of the *same* kernel on a periodic `w₁` (grating kinematics). Do not make Bloch the S11c-d class.

That is the user’s hybrid, and it is the right architecture. It is not a dodge if the lab comparison is explicitly `F(1)` vs Born `F(η)∼η²`, never a number-to-number from `η=1` in the truncated kernel.

---

## 4. The leaning is mis-framed; one omitted class is real

**The regime axis as posed collapses two independent grades that S11c-a already split.**

- `η` = contrast (Born / scatterer strength).
- `σ_W = η W̄₀/L_W` = first-jet / sharpness (WKB / adiabaticity).
- `k L_W` is a third, kinematic, scale (whether the edge is sharp on the *wave*).

“Weak-gradient (perturbative in `η`)” identifies gradient with contrast. They are not the same. A slit is large contrast *and* large gradient. WKB is small gradient, possibly large contrast. Born is small contrast, gradient optional.

If S11c-d is specified as **weak-gradient Born**, a builder will expand in `σ_W` or send `L_W → ∞`. That kills the Fourier form factor of a localized edge and slides the object toward WKB — the one class that does *not* match a slit. That is the quiet throw-away, not the deferral of `η = O(1)`.

**Say instead:** localized class, Born in `η`, `σ_W` an independent live grade. Do not expand the profile in derivatives. Do not freeze a global `ω(k)`.

**Omitted class worth naming (and rejecting for S11c-d):** piecewise-uniform jump — two *uniform* S11b slabs of different `W` joined by matching. Each side is `η`-exact (uniform). That is the honest order-unity localized edge. It is not S11c-d: matching at finite `ΔW` is a new construction, the first-jet kernel integrated across a true discontinuity is again a thin-layer/small-integrated-strength approximation, and a real slit may be a terminus (`W→0`), which is a change of domain. Park it as the later non-perturbative reduction, not as a way to “keep `(ε,η)` honest” while setting `η = O(1)` in the truncated operator.

A specific shape (Gaussian, `sech²`, `tanh`) is a representative, not the class. The mixing vertex depends on `∇w₁`; a kink and a bump both give a localized vertex. If a representative is needed for tractability, a `tanh` kink is closer to an edge than a bump; do not freeze it as “the slit.”

---

## Direct answer to the leaning

The object is right: localized Born kernel now, order-unity slit later, `(ε,η)` polynomial intact, both frozen exports consumed as the vertex plus the sector resolvent.

The call is wrong if “weak-gradient” is taken literally. That expands `σ_W` and throws away the edge form factor the slit comparison needs. Defer **order-unity contrast**. Do not defer **sharpness** by Taylor-expanding it. Do not commit to order-unity `η` now: the exports are first shape order, and there is no honest `(ε,η)` bookkeeping in which `η = O(1)` in those operators.

The thing order-unity would capture that Born cannot — unitarity, finite-`η` resonances, true-step matching, a hole rather than a bump — is real, and is not secretly sitting inside the truncated kernel. Naming those as later objects is honest. Pretending the Born coefficient at `η = 1` is the slit is the error that would invalidate the lab bound.
