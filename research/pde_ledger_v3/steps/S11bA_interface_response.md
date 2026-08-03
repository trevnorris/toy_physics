# S11b-A · What the bulk does when the faces move

**Sector 1 (light), step 4a.** Walked side by side with the user, 2026-08-03.

---

## What it is

S11 closed with three open questions that reduce to one unbuilt object: the linear brane–bulk interface
coupling law. **S11b** is that object. ⭐ **S11b-A is its first half** — the bulk's response to moving
faces, and the identity that projects a four-dimensional balance law onto the slab. ⛔ The assembly —
stress balance, the brane's elastic sector, the coupled dispersion, the transverse channel, and any
radiate/bound/confined/observable verdict — is **S11b-B's**.

## ⚠ Why it is a half, and not by choice

Two attempts to specify the whole interface in one pass were **rejected by review before any build ran**.
The first specified a system with **no linear coupling at all** — both engines would have returned a
block-diagonal result and concluded *"the longitudinal does not couple,"* an artifact of an omitted input,
identical in both and invisible to their agreement. The second introduced new defects **caused by its own
fixes**.

⭐ **Two failed designs means the shape is wrong.** The split was the fix, and the evidence that it was the
right cut is that the next review round produced findings that were **all local** — a sign, a missing
radiation condition, an undefined quantity — where the previous two had been structural. ⚠ Dropping the
assembly also removed two defects **outright rather than by patching**: `ρ_br` was defined twice
incompatibly and vanished with the brane sector, and `B_comp⁰` collided with a registry symbol only one
engine could read, which vanished with the registry.

---

## The computed result — two blind engines, written independently

```
q² = ω²/c_s0² − k²                        bulk normal wavenumber
Z  = ω ρ_m / q_out                        pressure at a face ÷ its OUTWARD normal velocity
```

| regime | `q_out` | `Z` | meaning |
|---|---|---|---|
| `q² > 0` propagating | `q` | real, `ωρ_m/q` | **radiation resistance** |
| `q² < 0` evanescent | `iα` | imaginary, `−iωρ_m/α` | **inertial loading**, `m_add = ρ_m/α` per face |
| `q² = 0` grazing | `0` | **singular** | `Z → ∞`, `m_add → ∞`; driven face has no solution, undriven has a free amplitude |

**Permeable faces**, with a memory-carrying closure `J_± = Λ_p(ω)δp + Λ_V(ω)V_±`, `Λ(ω) = Λ⁰/(1 − iωτ)`:

```
Z_perm = (ρ_m r + Λ_V⁰) / (y r − Λ_p⁰) ,     r = 1 − iωτ ,     y = q_out/ω
```

**Projection identity**, from integrating 4D continuity against a window `Ω(w)`:

```
∂_t ∫Ω ρ dw + ∇_x·∫Ω j_x dw = Ω(w₁)j^w(w₁) − Ω(w₂)j^w(w₂) + ∫ j^w Ω′ dw
```

**Parity**, `Ω` even, on the **symmetric** interval `[−L, L]`: `j^w` even ⇒ source `= 0` **exactly**;
`j^w` odd ⇒ `−2Ω(L)j(L) + 2∫₀ᴸ j Ω′ du`, not identically zero.

**Flux channels:** net accretion `= −(J₊ + J₋)` (symmetric) · through-flow `= (J₊ − J₋)/2` (antisymmetric),
with `J_± = ρ_m(v_w − ∂_tζ_±)(±1)`.

**Two-face structure:** two amplitudes survive the radiation condition, two face conditions fix them, and
⭐ **the per-face `Z` is the same for the thickness and centre-shift combinations** — the half-spaces
decouple, so each face sees only its own.

---

## ⭐⭐⭐ The headline — the leak is FREQUENCY-DEPENDENT

**1. Evanescent dissipation is created entirely by the interface closure.** At `Λ_p⁰ = Λ_V⁰ = 0` the
propagating `Re Z` survives as `ρ_m/a_q` — that is radiation resistance, present with **impermeable**
faces — while the evanescent `Re Z` vanishes identically. ⇒ ⭐ *bound-vs-radiating* is really *"what is
`Re Z` in the evanescent regime"*, and that is now a formula rather than a question.

**2. ⭐⭐ The velocity coupling dissipates if and only if `ωτ ≠ 0`.**

```
Λ_p⁰ = 0 , ωτ ≠ 0  ⇒  Re Z = Λ_V⁰ ωτ / (b_q(1 + (ωτ)²))     NONZERO
Λ_p⁰ = 0 , ωτ = 0  ⇒  Re Z = 0
```

⇒ A velocity-coupled flux does no net work when it responds **instantly**, because it stays exactly in
quadrature. Give the conversion a **finite rate** and it lags — and the lag does the work.

**3. The channel freezes out at high frequency.** As `ωτ → ∞`, `Z_perm` reduces to the **impermeable**
result in both regimes. ⇒ ⭐ **slow disturbances see the leak; fast ones do not.**
⛔ This does **not** establish anything about photons — the transverse channel is S11b-B's. ⭐ It
establishes that the question is *askable*, which a frequency-flat loss could not express.

---

## ⛔⛔ CORRECTION — the discarded convective order, fixed HERE

Both engines reported the rest-frame linearisation's discarded correction as `O(v₀/c_s0)`. ⛔ **That is
right only on the sound cone.** From `(∂_t + v₀∂_w)²φ = c_s0²∇₄²φ` with `φ ∝ exp[i(k·x + q_n w − ωt)]`:

```
(ω − v₀q_n)² = c_s0²(k² + q_n²)   ⇒   relative correction = −2v₀q_n/ω + (v₀q_n/ω)²
```

⇒ **the correct general statement is `O(v₀|q_n|/ω)`**, reducing to `O(v₀/c_s0)` only for `|q_n| ~ ω/c_s0`.
In the **evanescent** regime `q_n = iα` it is `2(v₀/c_s0)(αc_s0/ω)`, which **exceeds first order whenever
`k c_s0 ≫ ω`** — ⚠ **precisely the regime carrying `m_add` and the evanescent dissipation.**

⭐⭐ **This was a defect in the SPECIFICATION, not in either engine.** §4/A8 presupposed the answer was a
power of `v₀/c_s0`; both engines implemented that faithfully and agreed. ⇒ **the one failure mode
dual-engine structurally cannot catch**, and three review rounds missed it because each checked whether
the script matched the spec. ⇒ **S11b-B's directive review must ask: does the specification presuppose the
form of any answer?**

⚠ **Consequence for S11b-B:** it works in the evanescent regime, so it ⛔ must not inherit the sound-cone
answer. The background flow `v₀` is the user's dark-energy leak, and its neglect is **least** justified
exactly where the interesting physics sits.

---

## Pre-registration scorecard

Predictions committed **before any script existed** — `e6d40a29`, `0acb69cf`, `d5732c46`.

| | prediction | outcome |
|---|---|---|
| **P-A** | per-face `Z` same for both combinations | ✓ confirmed |
| **P-B** | `Λ_p` dissipates, `Λ_V` does not | ⛔ **REFUTED** — see headline 2 |
| **P-C** | symmetric combination = net accretion | ✓ confirmed (sign explicitly not predicted) |
| **P-D** | **both** controls move the parity result | ✓ confirmed — parity needs an even window **and** a symmetric interval |
| **P-E** | `Z`, `m_add` singular at grazing | ✓ confirmed |
| **P-F** | `τ→0` reproduces memoryless | ✓ confirmed — ⛔ true by construction, worthless as declared |
| **P-G** | `ωτ→∞` freezes the channel out | ✓ confirmed |
| **P-H** | ⭐ *deliberately none* | the shape returned rich, and unanchored |

⭐⭐ **The refutation is worth more than the six confirmations, and it nearly did not happen.** The closure
was originally written **memoryless** — under which `Λ_V`'s dissipative term is identically zero and `P-B`
would have scored as confirmed, banking a false mechanism into S11b-B where the dissipation question is
decided. ⚠ **The user objected** that removing time dynamics has repeatedly broken derivations in this
programme; the closure was generalised to carry a relaxation time, and that is the only reason the result
exists. ⇒ `feedback_static_or_instantaneous_check`.

⚠ **And `P-A` was nearly a false confirmation of a different kind:** the first build produced it from a
ratio returning its own input. Two review legs caught it, the two-face problem was actually posed, and only
then did the agreement mean anything. ⇒ *a matching number is not evidence.*

---

## Verification

**Two engines, written independently, agreeing on every value.** The blind Mathematica audit was built
**first**; the SymPy audit was built while the `.wl` was **out of the tree**, restored **byte-identical**
to its committed blob (`cc62f0a6`), and the build log carries **zero** barred-path accesses. ⭐ For this
sub-step the SymPy engine was **also barred from the registry**, so both engines faced identical
information.

**Cross-engine:** ~43 tags each; every physics value agrees. `Z_perm` differs by an overall factor of `i`
in numerator and denominator, verified algebraically identical.

**Review legs: six.** Two on the directives (which rejected two whole-interface revisions before any build
ran), two on the `.wl`, two on the `.py`. ⚠ A seventh died mid-run at exit 144 from **Mathematica seat
contention** — two legs ablating the same script exceed a 2-seat licence. ⇒ such legs are now serialized.

**Independent re-derivation:** three review legs re-derived the values **by hand before reading the
scripts** — 25, 43 and ~30 values respectively, **zero disagreements**.

### ⚠ Known limits — recorded, ⛔ not fixed

- ⛔ **19 of 23 gated tags in the `.wl` cannot fail under upstream corruption.** Its `VERDICT: PASS`
  certifies less than it appears to. ⭐ The four that hold compare against a **physics-independent
  constant**. ⚠ A repair pass attempted this and made it worse — it added 183 lines whose checks were of
  the form `assoc["key"] === variable` where the association was built as `"key" -> variable`, i.e.
  `x === x`. ⇒ **the directive mandated a tautology while banning independent re-derivation.**
  ⭐ **Not evidence against the values** — five independent derivations back them. ⇒ red-team item.
- ⭐⭐ **The reusable error: trying to make ONE engine self-certifying.** The architecture's answer to
  *"is this value right?"* is the **blind second engine that can disagree**, ⛔ not a longer check list.
- **The convective-order tag text still reads `O(v₀/c_s0)`** in both scripts. Corrected above; the script
  text is a red-team item.
- **One engine renamed a tag** (`RELATIVE_FLUX_CHANNELS` → `RELATIVE_FLUX_COMBINATIONS`); values agree.
- **The amplitude and condition counts are not printed**, though the solve is genuine.

---

## What S11b-A hands S11b-B

| object | value |
|---|---|
| bulk response | `Z = ωρ_m/q_out`, three regimes, grazing singular |
| inertial loading | `m_add = ρ_m/α` per face |
| permeable response | `Z_perm(Λ_p⁰, Λ_V⁰, ωτ)` |
| projection source | `Ω j^w\|_{w₁} − Ω j^w\|_{w₂} + ∫ j^w Ω′` |
| parity rule | requires an even window **and** a symmetric interval |

⛔⛔ **Three traps for B, all measured here:**
1. **`m_add`'s minus sign is a CONVENTION artifact, ⛔ not negative inertia.** Against the *outward*
   acceleration both faces give `+ρ_m/α`. ⛔ Do not consume `(ρ_m/α, −ρ_m/α)` as an inertia.
2. ⛔ **Do not inherit `O(v₀/c_s0)`.** Use `O(v₀|q_n|/ω)`; B lives in the evanescent regime.
3. ⛔ **Do not read propagating `Re Z` as evidence of leakage** — it is radiation resistance and exists
   with impermeable faces. Only the **evanescent** `Re Z` is created by the closure.

**Registry:** ⭐ **no rows earned, by design.** This sub-step derives no brane quantity; `Λ_p⁰`, `Λ_V⁰` and
`τ` are closure coefficients of an interface law that S11b-B will assemble. ⇒ registry work belongs there.
