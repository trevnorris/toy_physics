# S11b-B — pre-registered predictions

**Written 2026-08-03, BEFORE any S11b-B script exists.** Committed for timestamp priority, then moved out
of the tree.

⚠ These come from the orchestrator's hand derivation during moves 1–2 of the walk. ⭐ **Some of that
derivation is solid and some is a sketch, and the difference is marked** — a sketch that comes back
confirmed proves much less than a derivation that comes back refuted.

---

## Derived during the walk — I expect these and would be surprised to be wrong

**P1 · The constraint.** Linearised slab-mass conservation gives, per unit `x`-3-volume,
`e_W + e_ρ = e − g·u` with `e = −∇·u`, `g = ∇ρ_br⁰/ρ_br⁰`. ⭐ **Exactly ONE internal degree of freedom
survives**, because the relation is algebraic and fixes the densification once thickness and in-plane
strain are known. **Confidence: high.**

**P2 · The effective modulus.** With `β ≡ B_ρ⁽³⁾/W₀²`:

```
B_eff = B_ρ⁽³⁾ (−μ_W ω² + k_W − iωZ) / (−μ_W ω² + β + k_W − iωZ)
```

**Confidence: high on the structure, medium on the exact placement of `W₀` factors.**

**P3 · Limits.** `ω→∞ ⇒ B_eff → B_ρ⁽³⁾`. `ω→0 ⇒ 1/B_eff = 1/B_ρ⁽³⁾ + 1/(k_W W₀²)`, i.e. compliances add.
⭐ **A modulus measured with the thickness NOT a degree of freedom is the `ω→∞` (stiff) limit.**
**Confidence: high on the limits, medium-high on the identification.**

**P4 · Effective inertia.** `μ_W + ρ_m/α` in the evanescent regime — ⭐ **not** `μ_W + 2ρ_m/α`, because
each face moves by half the thickness perturbation and the factors cancel. **Confidence: medium-high.**

**P5 · The dispersion is not a cone.** `ρ_br ω² = B_eff(ω,k)k²` is transcendental — `B_eff` depends on `ω`
directly and through `α(ω,k)` inside `Z`. **Prediction: no closed form for `ω(k)` in general.**
**Confidence: medium-high.**

**P6 · The root goes complex when `Re Z ≠ 0`**, and `Im ω` is the leakage rate. **Confidence: high.**

## ⭐⭐ Explicit NON-predictions — where a match must NOT be read as confirmation

**P7 · The transverse coupling's scaling.** ⛔ **I do NOT predict what power of the background gradient
appears.** My walk sketch said the coupling goes as `g` and the loss as `|g|²`, ⚠ **but the back-reaction
path was sketched, not derived**, and it could change the power. **Confidence: none on the exponent.**
⭐ I predict only that **the channel exists for a transverse mode when `g` has a component along the
polarization**, and that it **vanishes identically when `g = 0`**. **Confidence: high on existence,
high on the uniform limit.**

**P8 · Whether the transverse dissipation is large enough to matter.** ⛔ No prediction. It requires a
magnitude I do not have.

**P9 · The validity failure region.** ⛔ No prediction. §3 trap 3 says the discarded convective term
exceeds first order somewhere in the evanescent regime; ⛔ I have not worked out where the boundary sits.

## Predictions from the preceding sub-step's results

**P10 · Transverse dissipation, if it exists, inherits the `ωτ` structure** — suppressed as `ωτ → ∞`
because the conversion channel freezes out. **Confidence: medium.** ⚠ This is the physically interesting
one: it would mean **light near matter converts less at higher frequency.**

**P11 · Impermeable control** (`Λ_p⁰ = Λ_V⁰ = 0`): the **evanescent** imaginary part vanishes, the
**propagating** one survives as radiation resistance. **Confidence: high** — this is the previous
sub-step's result restated, so a failure here is a wiring bug, not physics.

**P12 · Uniform control:** transverse coupling vanishes **identically**, and `B_eff` is unchanged.
**Confidence: high.**

---

## ⚠ What I most expect to be wrong

**P4 and P2's `W₀` bookkeeping.** Factor placement between `B_ρ`, `B_ρ⁽³⁾`, `β` and `k_W W₀²` is exactly
the kind of thing I get wrong by hand and the engines get right. ⭐ A disagreement there is a **correction**,
not a finding about the physics.

⛔ **And the one I would bet against hardest is P7's absent exponent being absent for the right reason.**
I withheld it because the derivation was a sketch — ⚠ but I *do* privately expect `|g|²`, and if the
engines return exactly that I must not treat it as confirmation of a derivation I never performed.
