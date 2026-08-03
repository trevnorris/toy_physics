# S11b-A — pre-registration addendum 2: the relaxation time

**Written 2026-08-03, BEFORE any S11b-A script exists.** The blind build was launched, then **stopped
before completion** when the user objected that removing time dynamics has repeatedly broken derivations
in this programme. The closure was generalised from memoryless to one relaxation time `τ`; ⭐ nothing was
read from the aborted run, and no artifact from it survives.

⚠ Prior pre-registrations stand: `e6d40a29` (S11b) and `0acb69cf` (S11b-A).

---

## P-F · The `τ → 0` limit

**Prediction: it reproduces the memoryless closure exactly.** ⛔ This is true by construction and is
therefore **not evidence of anything** — it is a wiring check, and I record it only so that a *failure*
here would be recognised as a bug rather than a result. **Confidence: certain, and worthless.**

## P-G · ⭐⭐ High-frequency freeze-out

**Prediction: as `ωτ → ∞` the conversion channel shuts off** — the interface cannot convert fast enough to
follow the wave — **and `Z` approaches the impermeable result of A4.**
**Confidence: medium-high.**

⭐ **Why this matters well beyond A5.** If it holds, the permeable channel is a **low-frequency**
phenomenon: slow disturbances see the leak, fast ones do not. ⚠ Light is fast. ⇒ a frequency-dependent
leak would be a mechanism by which the interface can be lossy for the background transfer and still
lossless for photons — **which the memoryless closure could never have expressed**, because it makes the
loss frequency-flat.

⛔ **I am not predicting that this rescues photon stability.** S11b-A cannot reach that question; the
transverse channel is S11b-B's. I am recording only that the generalisation *makes the question askable*,
where `τ = 0` foreclosed it.

## P-H · Where the dissipation sits in `ωτ`

`Re[Λ(ω)] = Λ⁰/(1 + ω²τ²)` and `Im[Λ(ω)] = Λ⁰ωτ/(1 + ω²τ²)`, so the in-phase part of the **flux response**
falls monotonically with `ωτ` while the out-of-phase part peaks at `ωτ = 1`.

⛔ **I do NOT predict that `Z`'s dissipative part follows either shape.** `Z` is a different combination —
the closure feeds through the bulk solution — and I have not worked it out. **Confidence: none.** ⭐ This
is deliberately an open prediction: the shape of dissipation versus `ωτ` is the most informative thing A5
can return, and I want no anchor on it.

## P-I · Whether one `τ` suffices

**Prediction: a single shared `τ` will prove sufficient for A5's tasks**, i.e. no reported quantity will
require the two coefficients to relax at different rates. **Confidence: low-medium.** ⛔ If a result
genuinely needs two, that is a finding about the interface, not a defect in the spec.

---

## ⚠ What I most expect to be wrong here

**P-G's clean freeze-out.** "The interface can't keep up, so the channel closes" is the same species of
tidy general argument that has produced every corrected error in this programme so far. ⭐ The honest
position is that `τ` belongs in the closure because the conversion has a finite rate — ⛔ **not** because I
know what the rate does.
