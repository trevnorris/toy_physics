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

## Verification

**Two engines, written independently, agreeing on every value.** The blind Mathematica audit
(`mathematica/S10_two_transverse_photons_mathematica_audit.wl`) reads nothing and was built first; the
SymPy audit (`scripts/S10_two_transverse_photons_sympy_audit.py`) imports the shared registry and was
built while the `.wl` was **quarantined out of the tree**, so it could not transcribe from it. The
quarantine held — the restored `.wl` is byte-identical to its committed blob.

⭐ **Every value also matches predictions committed BEFORE either script ran**
(`steps/S10_PREREGISTERED_PREDICTION.md`, `bc276485`), including both flagged uncertainties: the curl
reduction holds at `D=3`, and `D=2` is **not** degenerate in kind.

**Form control** — replacing the antisymmetric-derivative stiffness with full gradient-squared collapses
the two roots into **one**, nullity `D`, and the `ω²=0` root disappears: the longitudinal propagates too.
⭐ **A coefficient control cannot test this step** — rescaling `μ_R` moves the roots and leaves the
nullities at 1 and 2. ⚠ That is the S9 review's lesson applied: **a coefficient control tests the
arithmetic; only a FORM control tests the physics.**

Registry gates, all green with the payload **unchanged**: acceptance `MATCH` (10→6, 10→6, 10→5),
dimensional homogeneity `PASS`, able-to-fail `PASS`, 11 tests.
