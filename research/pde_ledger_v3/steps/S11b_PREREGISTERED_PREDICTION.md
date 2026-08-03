# S11b — pre-registered predictions

**Written 2026-08-03, BEFORE any S11b script exists.** Orchestrator's expectations, committed so the
timestamp proves priority, then moved out of the tree for the build's duration.

⚠ **Confidence is stated per item and it varies a lot.** Several entries are explicit *non*-predictions.
⛔ A non-prediction is not hedging — it marks where I genuinely do not know, so that a matching result
there cannot later be read as confirmation of anything.

---

## T1 · The projection identity

Integrating 4D continuity against a window `W(w)` and integrating by parts:

```
finite window:    S_leak = −[W j^w]_{w1}^{w2} + ∫_{w1}^{w2} W′(w) j^w(w) dw
infinite window:  S_leak = ∫_{−∞}^{∞} W′(w) j^w(w) dw          (fast decay kills the boundary term)
```

**Confidence: high.** This is integration by parts; I expect exact agreement.

## T2 · The parity selection rule

For a window **even** in `w` (`W(w) = W(−w)`, so `W′` is **odd**):

- `j^w` **even** in `w` — uniform normal translation, i.e. the **branon** — ⇒ `S_leak = 0` **exactly**.
- `j^w` **odd** in `w` — the **breathing / width** mode — ⇒ `S_leak ≠ 0` generically.

**Confidence: high on the mechanism, medium on exactness.** I predict the annihilation is *exact* for any
even window, not merely leading-order.

⚠ **I predict this DEPENDS on the window being even.** An asymmetric window should break it and let the
branon leak. That is the FORM control in T9.

## T3 · The dynamical window

Redoing T1 with `W = W(w; x, t)` should produce **additional** terms absent from the static derivation.

⛔ **I do NOT predict the complete enumeration.** I predict only: **at least one new term, carrying a time
derivative of the window**, and that these terms are what couple brane compression to the bulk.
**Confidence: medium on existence, none on the list.**

## T4 · The bulk radiation impedance

```
κ² = ω²/c_s0² − k²          Z = ρ_m ω / κ           ρ_m = m·ρ0   [M L⁻⁴]
```

| branch | `Z` | consequence |
|---|---|---|
| `κ² > 0`, `κ` real, outgoing | `ρ_m ω/κ`, **real** | resistive ⇒ radiates |
| `κ² < 0`, `κ = i\|κ\|`, decaying | `−i ρ_m ω/\|κ\|`, **imaginary** | reactive ⇒ added mass `m_add = ρ_m/\|κ\|` |

**Confidence: high.** Checks: at `k → 0`, `κ → ω/c_s0` and `Z → ρ_m c_s0` (plane-wave impedance). At
`k ≫ ω/c_s0`, `|κ| → k` and `m_add → ρ_m/k`.

⚠ **This is the IMPERMEABLE-face result.** See T5.

## T5 · ⛔ THE LEAK CORRECTION — an explicit NON-PREDICTION

How a permeable face changes `Z` is **the single most important output of this step and I do not know it.**

⭐ My only stated expectation, and it is a **guess at low confidence**: because material crossing the face
is an irreversible phase conversion, I would guess it contributes a **resistive** part to `Z` **even in the
evanescent branch** — i.e. loss without far-field radiation. ⛔ I am **not** predicting the magnitude, the
functional form, or that it is nonzero at linear order.

⚠⚠ **If that guess is right it is dangerous, not reassuring** — see T7.

## T6 · The coupled system

Assembled from the brane longitudinal equation, the width mode (symbolic inertia `μ_W`, symbolic wall
stiffness), and the bulk loading:

- **P6a — the feedback is self-limiting.** Added mass softens the width mode ⇒ `B_wall` drops ⇒ `B_comp`
  drops (series) ⇒ `c_L` drops ⇒ `κ²` more negative ⇒ `|κ|` larger ⇒ `m_add = ρ_m/|κ|` **smaller**.
  ⇒ negative feedback with a **unique fixed point**. **Confidence: medium.**
- **P6b — the bare inequality `c_L ⋛ c_s0` is the WRONG question**; the fixed point's location is the
  right one. **Confidence: medium-high.**
- **P6c — at long wavelength the bulk DOMINATES the width mode's inertia**, `m_add/μ_face ~ 1/(kW₀)` for
  comparable densities. **Confidence: medium-high.**

## T7 · ⭐⭐ The transverse check — the falsifier

The coupling coefficient of the **transverse** brane mode to `δW` must be **computed**, ⛔ not argued from
`∇·u = 0`.

**Prediction: it is ZERO at linear order in a homogeneous brane.** **Confidence: medium-high on the
result, but this MUST be able to come out nonzero** or the check is worthless.

⛔⛔ **Why this is the falsifier and not a formality.** If T5's guess holds (face motion ⇒ irreversible
conversion ⇒ loss) **and** the transverse mode reaches the face at linear order, then **photons acquire a
loss rate** — and photons from distant galaxies arrive. That would be a **no-go against observation**, and
it must be reported as one, ⛔ never softened.

## T8 · Dimensions

```
[Z] = M L⁻³T⁻¹        [m_add] = M L⁻³       [μ_W] = M L⁻³
[σ_wall] = M L⁻¹T⁻²   [W₀] = L              [ρ_m] = M L⁻⁴
[S_leak] = M L⁻³T⁻¹  (mass convention, projected)   or  L⁻⁴T⁻¹ (number convention, 4D)
```

## T9 · Controls

- ⭐ **FORM control A — give the bulk a shear modulus.** The transverse mode should then couple to the
  bulk **directly**, bypassing the width channel entirely. This tests the *shape* claim "the width is the
  only linear channel." Predicted: transverse acquires a bulk-dependent term it does not otherwise have.
- ⭐ **FORM control B — break the window's evenness.** T2's annihilation of the branon should **fail**.
- ⛔ **A coefficient control (rescaling `σ_wall`, `ρ_m`) cannot test either claim** — scaling never leaves
  the family, so it cannot test a channel *count* or a parity *selection rule*.

## T10 · Registry

**Predicted new:** `σ_wall` and `W₀`, both **postulated with a retirement condition at S5–S7** (the same
pattern `B_comp` got at S11). ⇒ ambient `12 → 14`.

⭐ **And a predicted REDUCTION, at medium confidence:** move 4's springs-in-series should become explicit,
expressing `B_comp` as a series combination of a density channel and a wall channel. If the density
channel is itself fixed by the medium block (`K`, `n_eos`, `ρ0`), then **`B_comp` retires from `parameter`
to `intermediate`** and its S11 retirement condition is partially discharged.

⇒ **Residue prediction: `7 → 9` if `B_comp` survives as a knob, `7 → 8` if it retires.** I predict `8`,
at medium confidence.

---

## ⚠ What I most expect to be wrong

**The `k^{3/2}` exponent.** Bulk-loaded, `ω_W² ~ σ_wall k²/(ρ_m/k)` gives `ω ~ √(σ_wall/ρ_m)·k^{3/2}`,
which would **correct move 5's `k²`**. ⛔ That is also the textbook capillary-wave/ripplon exponent, so I
may simply be pattern-matching to a formula I already know rather than deriving it. **Confidence: low, and
deliberately so.** ⭐ If the engines return `k²`, or an exponent that depends on something I have not
considered, that is the more informative outcome.
