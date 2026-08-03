# S11b-A — pre-registered predictions (addendum)

**Written 2026-08-03, BEFORE any S11b-A script exists.** The original S11b pre-registration
(`S11b_PREREGISTERED_PREDICTION.md`, blob `e6d40a29`, commit `604890b6`) still stands and covers
`A1 ≡ T1`, `A2 ≡ T2`, `A4 ≡ T4`, `A5 ≡ T5`, `A6 ≡ T8`. ⭐ This addendum covers **only what S11b-A adds**:
the parity split, the algebraic closure, the relative-flux identification, and the split control.

---

## P-A · Impedance by parity

**Prediction: the per-face response `Z` is the SAME for the `δW` and `h` combinations**, because each face
sees an identical acoustic half-space and the half-space problem does not know which combination it
belongs to. ⇒ the parity difference should appear **only** in how the two face pressures combine into a
net force on the slab — `h` driving opposite-sign pressures above and below, `δW` driving same-sign ones.
**Confidence: medium-high.**

⚠ If `Z` itself comes out parity-dependent, my picture of the two half-spaces as independent is wrong,
and that is more interesting than the prediction.

## P-B · ⭐⭐ Which coefficient carries the dissipation

With `J_± = Λ_p δp|_face + Λ_ζ ∂_t ζ_±`:

- **`Λ_p` carries the dissipation.** It relates a flux to a *force*, which is the signature of an
  irreversible kinetic coefficient. **Prediction: `Λ_p` contributes a real, in-phase-with-velocity part to
  `Z` in ALL THREE regimes — including the evanescent one**, i.e. **loss without far-field radiation**.
- **`Λ_ζ` does not.** It relates a flux to a *velocity*, which reads as a kinematic slip that renormalises
  the boundary condition without adding loss.

**Confidence: medium on `Λ_p`, low-medium on `Λ_ζ`.** ⚠ This is the sharpened form of the original
pre-registration's explicit **non**-prediction `T5`, and it remains the single most important output of the
S11b programme. ⛔ If `Λ_ζ` also produces a real part, my reading of it as purely kinematic is wrong.

## P-C · Relative flux

**Prediction: the SYMMETRIC combination of `J₊` and `J₋` is net accretion by the slab, and the
ANTISYMMETRIC combination is through-flow.** ⛔ I do **not** predict the overall sign; my `§3` definition
carries `∓` and `(±1)` factors and I have not tracked them by hand.
**Confidence: high on the symmetric/antisymmetric assignment, none on the sign.**

## P-D · ⭐ The split control — I expect BOTH halves to move

`A7-A` (asymmetric window, symmetric interval) and `A7-B` (even window, asymmetric interval).

**Prediction: BOTH move `A2`'s result.** ⇒ the parity result is **not** a property of the window alone; it
requires the interval to be symmetric as well. **Confidence: medium-high**, and it is why the control was
split — a single combined control would have let a domain artifact pass as a selection rule.

⭐ **If only the window control moves it, the selection rule is stronger than I think** and `A2`'s
"exact" claim holds on any interval. That would be the better outcome and I am not predicting it.

## P-E · Grazing case

**Prediction: `Z` and the inertial loading are singular as the bulk normal wavenumber → 0** (the
denominator vanishes), and the singularity is **physical, not an artifact** — it is the resonance where
the disturbance runs exactly along the interface at the bulk sound speed.
**Confidence: medium on the singularity, low on the interpretation.**

---

## ⚠ What I most expect to be wrong here

**P-B's clean split.** "Force-coupling is dissipative, velocity-coupling is kinematic" is a tidy general
principle of exactly the kind that has produced three corrected errors in this programme already. ⛔ It is
the least earned item on this page and the one I would bet against first.
