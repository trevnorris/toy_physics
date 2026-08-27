# One-shot CAS consult — S11c-a projection object: two representational choices

You will not remember this conversation after you answer. Give a **complete, self-contained** answer in
one shot. Do not ask follow-up questions.

## What this is

Two independent computer-algebra audits ("engines") re-derived the **same projection object** from the
supplied spec quoted below. They made two **different representational choices** while doing so. Your job:
determine **by computation** whether those two choices produce the **same** projection object or
**different** ones. There are two sub-questions. For **each**, you must produce a runnable CAS script and
its **literal stdout**; a prose-only derivation will be **discarded** (a hand derivation that agrees with
itself is worthless here). Do **not** assume the answer — compute the difference and report whatever it is.

## Source of truth — the supplied spec (verbatim)

### §1b · Bulk acoustics and the projection law
```
v_bulk = ∇₄φ ,   δp = −ρ_m ∂_tφ ,   ∂_t²φ = c_s0² ∇₄²φ .
j = ρ_4D v_bulk ,   ∂_tρ_4D + ∇₄·j = 0 .
```
`Ω` is a smooth slab window, approximately one inside and zero outside. The projection object is the
result of integrating this conservation law against `Ω` and **integrating by parts in `w`** (the normal
coordinate). It uses the dynamic, anchored window supplied in §3; it may not be replaced by a static
window before differentiation.

### §3a · Face maps and slab geometry (the window's arguments)
```
R_s^L(X,t) = ( x, s W_bg(x)/2 + ζ_s ) ,   R_s^M(X,t) = ( x, s W_bg(X)/2 + ζ_s ) ,   ζ_s = ζ_c + s δW/2 .
h_s^α(x,t) = s W_bg/2 + ζ_s ,   F_s^α(x,w,t) ≡ w − h_s^α(x,t) .
W^α ≡ h_+^α − h_−^α ,   ρ_4D^α ≡ ρ_4D,bg⁰ (1+θ) .
```
`w` is the normal coordinate; the slab spans a finite thickness `W_bg` about the face.

### §3c · The dynamic window and the current background
Let `G_s^α ≡ s F_s^α`. Supply one fixed smooth two-argument window `𝒪` and define
```
Ω^α(x,w,t) ≡ 𝒪( G_+^α(x,w,t), G_-^α(x,w,t) ) ,
```
with `𝒪` ≈ 1 when both arguments are negative and → 0 outside. **Its time dependence is retained until
after the projection identity and its shape derivative have been computed.** The rest-frame background
current `ρ_4D⁰ v_bulk⁰` vanishes and the supplied density background depends on the in-plane anchor, not
on `w`. `V_s⁰ = J_s⁰ = 0` (§2d).

The perturbation bulk current `δj = (δj_1, δj_2, δj_3, δj_w)` is the first-order (in the wave/perturbation
amplitude) part of `j = ρ_4D v_bulk`; it has zero background but is otherwise a generic bulk field.

## The two representational choices (stated neutrally — assume nothing about which is correct)

**On the perturbation current `δj`:**
- Choice ①: each component is represented as a **constant in `w`** — its value at the background face
  `h⁰` — carried together with its in-plane derivatives `∂_i δj_i` (i = 1,2,3). No `w`-dependence.
- Choice ②: each component is represented as a **full bulk field `δj_a(x, w, t)` that depends on the
  normal coordinate `w`**, kept **inside** a projection integral `∫ … dw` that is left formally
  un-evaluated (the engine differentiates the integrand but never integrates by parts and never
  substitutes the conservation law).

**On the dynamic window `Ω(G₊, G₋)`:**
- Choice Ⓐ: the shape derivative `d/dε` is taken **analytically**, expanded via the chain rule through
  `G₊, G₋` (producing explicit in-plane-profile-gradient × window-second-derivative terms).
- Choice Ⓑ: `d/dε` is applied to the **un-evaluated integrand**, and the `w`-integral is kept formal.

## Compute — two sub-questions, each needs a script + its literal stdout

**Q1 — does the current's `w`-dependence change the projection?**
With a generic smooth window `𝒪(g₊, g₋)` that → 0 as its arguments → +∞ (decaying/compact support in `w`),
and a generic perturbation current with genuine `w`-dependence, compute the **first-shape-order** projection
of `∂_tδρ + ∇₄·δj` against the window, integrating by parts in `w` per §1b. Then recompute with every
current component replaced by its value at the background face `h⁰` (constant in `w`). **Report the
difference (choice ② minus choice ①) as a CAS expression.** State whether it is identically zero; if not,
name exactly which term survives and what it depends on. State explicitly any window boundary property you
relied on (e.g. `𝒪 → 0` at the ends) — do not silently drop boundary terms.

**Q2 — does the analytic window shape-derivative equal the differentiate-under-the-integral form?**
For the same generic window, compute `d/dε [ ∫ 𝒪(G₊(ε), G₋(ε)) dw ]` two ways: (a) the analytic chain-rule
expansion through `G₊, G₋`, then integrate; (b) differentiate under the integral sign, keeping the integral
formal. **Report the difference.** State whether `d/dε` and `∫dw` commute here and whether the chain-rule
expansion is exact.

## Output
For each sub-question: the script path, its literal stdout, and a one-line **verdict**: `SAME object` or
`DIFFERENT` (with the surviving term). Do not assume boundary/limit behaviour without stating the property
used. Do not tailor the computation toward any particular result.
