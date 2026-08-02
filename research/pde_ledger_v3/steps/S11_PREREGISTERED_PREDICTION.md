# S11 — PRE-REGISTERED PREDICTIONS

⛔ **Written by the orchestrator BEFORE any S11 script exists.** Committed for timestamp priority, then
moved out of the tree for the duration of the build so no builder or review leg can read it.

⚠ Everything below is a prediction from the side-by-side walk of 2026-08-02, moves 1–6. Where I was
unsure I say so — an uncertainty named here and resolved by the engines is worth more than a lucky hit.

---

## The Lagrangian I expect the engines to work from

```
L = ½ ρ_br (∂_t u)²  −  ½ μ_R (curl u)²  −  ½ B_comp (∇·u)²
```

on a `D`-dimensional brane, `u` an in-plane `D`-vector, with S10's convention
`(curl u)² ≡ ½ Σ_{i,j} (∂_i u_j − ∂_j u_i)²`.

The symmetric-traceless (Cauchy) invariant is set to **zero** — `μ_br = 0`, the user's provisional call
of 2026-08-02, reopenable if a derivation forces it.

---

## P1 · The invariant count

A quadratic isotropic form in `∂_i u_j` on a `D`-brane has **exactly three** independent coefficients
(trace, symmetric-traceless, antisymmetric), for every `D ≥ 2`. S10 kept one of the three.

## P2 · The dynamical matrix and the spectrum

```
ρ_br ω² a  =  μ_R [ k²a − k(k·a) ]  +  B_comp k (k·a)
```

| root | nullity | eigenvectors |
|---|---|---|
| `ω² = (μ_R/ρ_br) k²` | `D − 1` | ⊥ `k` |
| `ω² = (B_comp/ρ_br) k²` | `1` | ∥ `k` |

Sum `= D`, for `D = 2, 3, 4, 5`. At `D = 3`: nullities **2** and **1**.

⭐ **P2a — the load-bearing one.** The transverse root is **exactly** S9/S10's, unchanged:
`c_γ² = μ_R/ρ_br`. The `B_comp` term contributes **identically zero** to the transverse sector, because
it enters only through `(k·a)`.

⭐ **P2b.** The two roots coincide **iff `B_comp = μ_R`**. Generic `B_comp` gives two distinct speeds and
no accidental degeneracy.

## P3 · Dimensions

```
[B_comp] = (2 − D, −2, 1)          at D = 3:  (−1, −2, 1)   — identical to [μ_R]
[c_L]    = (1, −1, 0)
```

All three stiffness coefficients share `[μ_R]`'s dimension, because each multiplies a **dimensionless**
object (`u` is a length, `∂` an inverse length).

⚠ **Uncertainty I am flagging, not predicting away:** whether S10's factor-of-½ convention on
`(curl u)²` forces a compensating factor on the trace term. I predict **it does not** — that
`c_L² = B_comp/ρ_br` comes out clean with the Lagrangian exactly as written above. If a stray 2 appears,
I was wrong about the convention, ⛔ not about the physics.

## P4 · The bulk matching

Brane mode `ω = c_L k`; bulk sound `ω² = c_s0²(k² + k_w²)`; `ω` and the in-plane `k` shared by symmetry.

```
k_w²  =  ω²/c_s0² − k²  =  k² ( c_L²/c_s0² − 1 )
```

- `c_L > c_s0` ⇒ `k_w` **real** ⇒ propagating bulk wave ⇒ the brane mode **radiates** and is a leaky
  resonance, not a bound mode.
- `c_L < c_s0` ⇒ `k_w` **imaginary** ⇒ evanescent `e^{−|k_w||w|}` ⇒ **bound**.

⭐ **P4a.** In registry quantities: **bound ⟺ `B_comp < ρ_br c_s0²`**. Dimensionally
`[ρ_br c_s0²] = (−1,−2,1) = [μ_R] = [B_comp]`, so the inequality is between commensurate objects.

⭐ **P4b.** The same matching applied to the **transverse** mode has no solution to attempt — the bulk
has no shear dispersion relation. Light's confinement is **unconditional**, not an inequality that
happens to hold.

## P5 · Controls

**FORM control** — replace the trace invariant with the symmetric-traceless one
(`μ_br S̃:S̃` in place of `½B_comp(∇·u)²`).
⭐ **Prediction: the TRANSVERSE speed CHANGES.** For `k·a = 0`, `S̃:S̃ = ½k²|a|² ≠ 0`, so the
symmetric-traceless invariant *does* charge transverse waves. ⇒ The trace invariant is the **unique**
one that lifts the longitudinal while leaving light untouched. That is S11's shape claim, and this is
what tests it.

**COEFFICIENT control** — rescale `B_comp → 2 B_comp`.
Prediction: `c_L²` doubles; `c_γ²` **unchanged**; both nullities **unchanged**. ⛔ Tests arithmetic
only — it cannot test the shape claim, because scaling never leaves the family.

## P6 · Registry

New rows: `Q.brane.B_comp` (`kind: parameter`, `counting_axis: continuous-model`, **postulated**) and
`Q.brane.c_L` (`kind: intermediate`). New relation `R5`: `c_L − √(B_comp/ρ_br) = 0`, designated output
`Q.brane.c_L`, denominator guard on `ρ_br` — structurally a twin of `R4`.

**Prediction: ambient continuous 10 → 12, residue 6 → 7**, the residue gaining `B_comp` and nothing
else: `{hbar, mass, K, rho0, rho_br, mu_R, B_comp}`. Discrete payload `{n_eos = 5, D_brane = 3}`
**unchanged**.

⚠ **Uncertainty:** the acceptance line reports three numbers (`10→6, 10→6, 10→5`) and I have not read
which measure each is. I predict the first two go to `12→7`; I am **not** predicting the third.

⭐ **P6a.** I predict the schema **cannot** express P4a's inequality as a relation — relations are
residuals equal to zero. It will have to live in the step record and, at most, as an assumption
predicate. ⛔ If a builder manufactures an equality to make it fit, that is a defect, not a success.

## P7 · Naming — a collision I expect to have to defend

`B_comp` is ⛔ **not** `K_br` (bulk modulus of the **rejected** Cauchy branch), ⛔ **not** `B_eff`
(`= ρ_B0²/χ_c`, the committed Lagrangian's bulk piece), and ⛔ **not** the bulk medium's modulus. Same
hazard shape as `μ_br` vs `μ_R`. I predict a review leg raises this; the answer is that it is a distinct
object and the name was chosen to avoid all three.

---

## What I expect to be WRONG about

1. Factors of 2 anywhere. I have derived none of the normalisations.
2. The third acceptance number.
3. Possibly the sign convention on `k_w²` if the bulk dispersion is set up with a different metric
   signature — the **condition** should survive any such choice, the intermediate expression may not.

## What would make me think the walk was wrong

- The transverse root moving at all when `B_comp` is switched on (contradicts P2a) — that would mean
  compression touches light, and moves 3–6 would all need redoing.
- The longitudinal nullity coming out `≠ 1` in any `D`.
- `[B_comp] ≠ [μ_R]` — that would mean the trace term is not the invariant I think it is.
