# S10 · Two transverse photons — and the count comes from the brane's own dimension

**Sector 1 (light), step 2.** Walked side by side with the user, 2026-08-02.

---

## What it is

S9 gave light a cone. S10 asks how many independent waves the brane's transverse sector actually
carries, and — the part that matters — **why that number**.

## What it does

Establishes that the mode count is `D_brane − 1`, computed rather than inferred, and registers the
brane's dimension as the input it has silently been all along.

## The walk, move by move

**Move 1 — decompose the amplitude about `k`.** For a plane wave `u = a·exp(i(k·x − ωt))`, the constant
amplitude `a` splits into a piece along `k̂` (one number) and a piece perpendicular to it. In `D`
dimensions the perpendicular subspace has `D − 1` directions:

```
D = 1 + (D − 1)
```

⚠ **The identification made here, flagged before it was used:** `u` is the **in-plane** displacement, a
`D`-vector on the brane — ⛔ **not** a `D+1`-vector including motion out of the brane into `±w`. The
out-of-plane motion is the **`h`-branon**, a separate object carrying charge (`h ≠ u_L`, stated in three
places in the corpus). ⭐ **User confirmed the picture.** ⛔ If that were wrong — if out-of-plane motion
belonged to the same elastic sector — the transverse count at `D=3` would be **3**, not 2.

**Move 2 — the curl sees only the perpendicular part.** `∇×u → i k × a`, so a longitudinal amplitude has
zero curl. ⭐ The clean algebraic form of this (contributed by a review leg, and it is the whole reason):

```
(1/2) Σ_{i,j} (k_i a_j − k_j a_i)²  =  |k|²|a|² − (k·a)²
```

⇒ the stiffness quadratic form **penalises only `a_⊥`**.

**Move 3 — the equation of motion, per piece.** `ρ_br ω² a = μ_R[k²a − k(k·a)]`, giving `ω² = 0` for the
longitudinal and `ω² = (μ_R/ρ_br)k²` for the perpendicular directions.

⭐ **The longitudinal mode is not forbidden — it is FREE.** Zero curl means zero energy cost, so no
restoring force. ⭐ That is precisely what *curl-only* means: the medium resists twisting and is
indifferent to compression along the wave. An ordinary elastic solid resists both, and carries a
longitudinal wave.

⭐ **The two transverse modes are degenerate by SYMMETRY**, not coincidence: nothing in the equation
refers to *which* perpendicular direction `a_⊥` points, so both carry the same `ω`.

## ⭐⭐ The computed result — and the claim it settles

⛔ **The count is not read off a multiplicity in a factorised determinant** (which is how the
orchestrator had it by hand). Both engines compute the **nullity of the dynamical matrix at each root**.

| `D` | `ω² = 0` | `ω² = (μ_R/ρ_br)k²` | sum |
|---|---|---|---|
| 2 | nullity 1, ∥ `k` | nullity 1, ⊥ `k` | 2 |
| **3** | nullity 1, ∥ `k` | **nullity 2, ⊥ `k`** | 3 |
| 4 | nullity 1, ∥ `k` | nullity 3, ⊥ `k` | 4 |
| 5 | nullity 1, ∥ `k` | nullity 4, ⊥ `k` | 5 |

⛔⛔ **THE COUNT IS NOT THE CODIMENSION — and the sharp form of that is stronger than "wrong number".**
⭐ **The bulk never enters this computation at all.** There is no codimension anywhere in the
calculation; the nullity depends only on the brane's own dimension. ⇒ Codimension is not merely the wrong
answer (`4−3 = 1`), it is an **absent quantity**.

⇒ ⭐ **Light having exactly two polarisations is a statement that our space is three-dimensional.**

## The equations and dimensions

```
L = ½ ρ_br (∂_t u)² − ½ μ_R (curl u)²          an energy density on a D-dimensional brane
(curl u)² ≡ ½ Σ_{i,j} (∂_i u_j − ∂_j u_i)²      the antisymmetric derivative, defined in every D
```

⭐ **The reduction to the ordinary `|∇×u|²` at `D=3` is COMPUTED** (residual `0`), ⛔ not assumed.

**Dimensions as a closed function of `D`** — new here, and it does retroactive work:

```
[ρ_br] = (−D,   0, 1)          at D=3:  (−3,  0, 1)  ✓ registry
[μ_R]  = (2−D, −2, 1)          at D=3:  (−1, −2, 1)  ✓ registry
```

⭐⭐ **This EXPLAINS an S9 loose end.** An S9 ablation showed that changing the assumed brane dimension
corrupts both constituent dimensions while the speed's dimension check stays green — filed then as a
blind spot in the audit. It is not a blind spot but an **identity**:

```
[μ_R] − [ρ_br] = (2−D) − (−D) = 2   in the length slot, for EVERY D
```

⇒ The speed's dimension **cannot** see the brane dimension. A weak-looking check is a property of the
physics.

## What's new — the introduction inventory

| item | class | why |
|---|---|---|
| `D_brane = 3` | **postulated** | what light *requires* of the medium; ⛔ not derived from anything |
| the mode count | **derived** | nullity of the dynamical matrix; `D_brane − 1` transverse, 1 longitudinal |

⛔ **`n_pol` is deliberately NOT a registry row.** Established by test, ⛔ not preference:
`constraint_dimension` **raises** on a relation whose residual involves discrete symbols, and the schema
(`additionalProperties: false`) has no vocabulary for *"this discrete value derives from that one"*. ⇒ A
bare discrete row would make a **derived** count look like a **free choice** — worse than omitting it.

⭐ **`D_brane` closes a real hole.** S9's dimensions depend on it, and it was **declared nowhere**.

## Inputs, by locus

| input | from |
|---|---|
| `L`, `ρ_br`, `μ_R`, the curl-only form | **S9** (`steps/S9_light_requires_shear.md`) |
| `u` is the in-plane displacement | this step's flagged identification, user-confirmed |
| ⛔ nothing from the bulk | it does not enter |

## Regime

Linear response about an unstrained brane; small oscillations — half one's scope.

## Departure

⛔ **None at this step.** ⚠ But note what S10 does *not* deliver: `ω² = 0` means the longitudinal is
**non-propagating**, ⛔ **not absent**. Maxwell has no longitudinal photon at all. The slot still exists
here, and once the medium is allowed to be **compressible** that zero lifts. → **S11**.

⛔⛔ **And per the user (2026-08-02), that is NOT a defect to be fixed.** Exact Maxwell would be the
**failure**: Maxwell puts charge in by hand, so a model matching it exactly would have **no way to
physically anchor charge**. The extra mode is the anchor — it is what made the drum-head picture click.
→ `CHARTER.md#falsification-standard`, `V3_STEP_PLAN.md` § S11.

## ⭐⭐ VERIFICATION — rebuilt 2026-08-05 under the script-rebuild programme

⚠ **The engines were renamed off their result.** They are now
`mathematica/S10_brane_mode_spectrum_mathematica_audit.wl` (2983 tags) and
`scripts/S10_brane_mode_spectrum_sympy_audit.py` (4279 tags), rebuilt from a written specification
(`directives/S10_SHARED_PHYSICS.md`) under the standing rule that **a script may print computed objects
and may not state conclusions**. ⭐ Their stdout is committed under `mathematica/out/` and `scripts/out/`.

### ⭐ The rebuild confirms this step's recorded result

| claim in this record | rebuild |
|---|---|
| the mode-count table, all 8 cells (`D = 2,3,4,5` × 2 roots) | ✅ **agrees, both engines, every cell** |
| `[ρ_br] = (−D, 0, 1)`, `[μ_R] = (2−D, −2, 1)` | ✅ **agrees**, solved symbolically in `D` |
| the `D=3` reduction to `\|∇×u\|²`, residual `0` | ✅ **agrees, both engines** — ⚠ see the Q7 caveat |
| full-gradient stiffness collapses the roots to **one**, nullity `D` | ✅ **agrees, both engines** |
| rescaling `μ_R` moves the root and leaves the nullities | ✅ **agrees, both engines** |

⛔ **Nothing in the recorded result disagreed.**

### ⭐⭐ What the rebuild ADDED — ⛔ and the first item changes what the numbers mean

⛔ **The original classified each null vector `∥ k` or `⊥ k` by inspecting it.** The rebuild **computes**
the transverse count: `N3` is the rank of the dynamical matrix **stacked on the row `kᵀ`** — basis
independent — and `N7` counts the same thing by a **second, independent** routine.

⚠ At **every** `D`, the `ω² = 0` root has **corank 1 but transverse nullity 0.**
⇒ ⛔⛔ **a naive corank test counts that null vector as a mode, and would report 3 transverse modes at
`D=3` instead of 2.** ⭐ `N7`'s independent count agrees with `N3` in every package and every `D`
(residual `0` throughout).

⭐⭐ **The form control is sharper than recorded — it MIRRORS the step.** `XFORM_DIVONLY` (stiffness
penalising **divergence only**) returns **the same two roots**, `{0, μ_R k²/ρ_br}`, with the mode
assignment **exchanged**: the `ω²=0` root now carries corank 2 / **transverse 2** and the propagating root
carries corank 1 / **transverse 0**.
⇒ ⭐ **curl-only ⇒ longitudinal free, transverse on the cone; divergence-only ⇒ transverse free,
longitudinal on the cone.** ⛔ That is not an arithmetic check — it is the claim of this step, exhibited by
swapping the one structure it rests on.

⭐⭐ **And the degeneracy is confirmed to be a SYMMETRY.** `XFORM_ANISO` returns **three** distinct roots
— the degenerate transverse pair **splits** — and one branch stops being purely transverse
(corank 1, transverse nullity 0), because anisotropic stiffness no longer aligns the eigenvectors with the
`∥/⊥` decomposition. ⇒ break the isotropy, and the two-fold degeneracy lifts exactly as move 4 predicts.

⭐ `XFORM_SIGNFLIP` and `XCOEF_SCALE` reproduce `MAIN` cell for cell — ⚠ **a coefficient control cannot
test this step**, which is the S9 lesson applied.

### ⛔⛔ SUPPLIED RATHER THAN DERIVED — ⛔ do not read these as verified

- ⛔⛔ **`Q7`'s VERDICT was supplied by the specification.** §Q7 stated that for non-curl packages the
  residual is *expected to be nonzero* and that this is *the control working*, and **both engines read
  it**. ⭐ The residual is a computed polynomial and the engines agree on it; ⛔ the judgement of what it
  ought to be was given. ⇒ `directives/S10_SPEC_CONSISTENCY_FINDINGS.md`.
- ⛔ **The dimensional-homogeneity booleans are structurally vacuous.** `[u]` enters every coefficient with
  the same weight and cancels in every ratio, so ⛔ no cleverer version of the check exists. ⚠ Measured:
  perturbing `[u]` moves hundreds of dimension payloads and **not one** homogeneity boolean.
  ⇒ they are a **classification**, ⛔ not a test.
- ⚠ **The brane dimension is an INPUT**, swept. The step establishes `D_brane − 1`; ⛔ it does **not**
  establish `D_brane = 3`.

### ⛔ The two-engine comparison — what it caught, and what it covers

⭐⭐ **The cross-engine layer found a real defect that eight review legs had not.** `Q7`'s form control was
**dead in one engine**: its stiffness payload was **byte-identical across all six packages**, including
the two whose only purpose is to change the stiffness form, so its residual was `0` **by construction**.
⚠ **The cause was a defect in the shared specification** — it named one object in its instruction list and
a different one in an appended amendment, and each engine followed a different sentence.
✅ **Repaired**; all 18 `Q7` objects at `D3` now match across engines.

⛔⛔ **Coverage is a SUBSET, and the disagree count is NOT a quality metric.** Of ~2983 and ~4279 emitted
tags, **562 names exist in both engines** and **677 pairs** are configured for comparison — §8 pinned the
tag **prefix** and not the **quantity name**, so both engines matched all 13 `(package, D)` prefixes
exactly and then diverged on nearly everything to the right.
⚠ **The count did not move across the `Q7` repair** — `446/100` before, `446/100` after — because those
rows shifted from *different objects* to *the same object under a different symbol spelling*.
⇒ ⭐ **the buckets are the result.**
