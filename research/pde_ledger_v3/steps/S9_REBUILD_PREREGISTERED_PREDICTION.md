# S9 REBUILD — pre-registered predictions

**Written 2026-08-04 by the orchestrator, BEFORE running any script**, old or new, and before reading
either review leg's report. Committed so the timestamp establishes priority.

⭐ Under the three clauses a builder cannot type a conclusion, so this file no longer exists to blind the
builder. **It exists to test me.** If the rebuilt engine disagrees with anything below, the disagreement
is a finding and I record it — I do not reconcile it.

Derived by hand from `L = ½ ρ_br (∂_t u)·(∂_t u) − ½ μ_R (∇×u)·(∇×u)` with `u` a 3-component field.
Notation: `T = k² I − k kᵀ`, `k² = k·k`.

---

## The main derivation

| # | object | prediction |
|---|---|---|
| **A** | equation of motion | `ρ_br ∂_t²u = μ_R (∇²u − ∇(∇·u))`, i.e. `−μ_R ∇×(∇×u)` |
| **B** | dynamical matrix | `M = ρ_br ω² I − μ_R T` |
| **C** | `det M`, factored | `ρ_br ω² · (ρ_br ω² − μ_R k²)²` at `D=3` — in general `ρ_br ω² (ρ_br ω² − μ_R k²)^{D−1}` |
| **D** | complete root set for `ω²` | exactly two distinct roots: `ω² = 0` and `ω² = μ_R k² / ρ_br` |
| **E** | transverse test `M·T = 0` | passes **only** at `ω² = μ_R k²/ρ_br`; at `ω² = 0` it gives `M·T = −μ_R T² = −μ_R k² T ≠ 0` |
| **F** | transverse subset | a single root, `μ_R k²/ρ_br` |
| **G** | speed squared | `μ_R / ρ_br` |
| **H** | cone residual `ω² − c²k²` | `0` |
| **I** | `D[speedSquared, k]` | `0` (no residual `k` dependence) |

⭐ **E is the load-bearing prediction of the transverse-selection machinery.** I derived it as:
`T² = k⁴I − k²kkᵀ = k² T`, and at the propagating root `M = μ_R k kᵀ`, whose product with `T` vanishes
because `kᵀT = k²kᵀ − k²kᵀ = 0`.

## Dimensions, symbolic in `D`

Convention `[L, T, M]`. Energy density on a `D`-sheet `= [2, −2, 1] − D·[1,0,0] = [2−D, −2, 1]`.

| # | object | prediction |
|---|---|---|
| **J** | `[ρ_br]` | `(−D, 0, 1)` |
| **K** | `[μ_R]` | `(2−D, −2, 1)` |
| **L** | `[μ_R] − [ρ_br]` | `(2, −2, 0)` — ⭐ **independent of `D`** |
| **M** | dimension implied for `speed²` | `(2, −2, 0)` |
| **N** | dimension of a squared velocity | `(2, −2, 0)` |
| **O** | difference M − N | `(0, 0, 0)` |
| **P** | J and K at `D = 3` | `(−3, 0, 1)` and `(−1, −2, 1)` |

⭐ **L is the prediction I care most about.** If it holds, S9's recorded "blind to the assumed brane
dimension" limitation is not a blind spot in the audit but an **identity of the physics** — the speed's
dimension *cannot* see `D`. That reclassifies a filed defect as a property.

## The controls

| id | prediction for the emitted objects |
|---|---|
| **X1** `ρ_br → λ_ρ ρ_br` | speed² `= μ_R/(λ_ρ ρ_br)`; difference from main `= (μ_R/ρ_br)(1/λ_ρ − 1) ≠ 0` |
| **X2** `μ_R → λ_μ μ_R` | speed² `= λ_μ μ_R/ρ_br`; difference `= (μ_R/ρ_br)(λ_μ − 1) ≠ 0` |
| **X3** divergence-only | `M = ρ_br ω² I − μ_R k kᵀ`; `det = (ρ_br ω²)^{D−1}(ρ_br ω² − μ_R k²)`. Roots `{0, μ_R k²/ρ_br}` again — ⚠ **but the transverse test now selects `ω² = 0`**, not the propagating root. Difference from main `= −μ_R k²/ρ_br`. ⭐ The transverse sector loses its restoring force. |
| **X4** isotropic gradient | `M = (ρ_br ω² − μ_R k²) I`; `det = (ρ_br ω² − μ_R k²)^D`. ⭐ **ONE root only**, `μ_R k²/ρ_br`, with nullity `D` — the `ω²=0` root **disappears** and longitudinal and transverse degenerate. |
| **X5** flipped stiffness sign | `M = ρ_br ω² I + μ_R T`; transverse root `= −μ_R k²/ρ_br` **< 0** ⇒ `ω` imaginary ⇒ instability. Difference from main `= −2μ_R k²/ρ_br`. |
| **X6** anisotropic inertia | degeneracy broken for generic `k`; ⚠ **the transverse subset comes out EMPTY** — with `diag(ρ_br,ρ_br,ρ_z)` no single `ω²` annihilates the whole transverse subspace. |

---

## ⛔ A DEFECT I FOUND IN MY OWN DIRECTIVE WHILE DERIVING THIS

§5 asks, for **X4**, for *"the difference between its transverse root and its **non**-transverse root."*

⛔ **That is ill-posed.** Under X4 the root set **collapses to a single root** (prediction above) — which
is precisely the physics X4 exists to expose — so there is no non-transverse root to difference against.
A builder told to emit it would have to **invent** something, and a check whose only honest outcome is an
invented value is worse than no check.

⭐ **The well-posed replacement**, which works in every control including the collapsed and empty cases:
emit the **complete root multiset with multiplicities**, and the **nullity of `M` at each root**, and
difference the *multisets* against the main derivation's.

⚠ I am holding this repair until both review legs report, so I can see whether either finds it
independently. That is the honest test of the review, and folding my own fix in first would destroy it.

## What this build does NOT test

Recorded so a passing build is not read as more than it is.

- **P2** — that a scalar superfluid carries no transverse mode. Cited, never computed, and this build
  does not compute it either.
- **P4** — bulk shear-freeness. The bulk is absent from the action, so nothing here can test it.
- The wall width, background flow and dissipation are all sent to zero (directive §6).
- ⛔ S9 has **one** engine, so there is no cross-engine disagreement available at this step. The only
  external checks are the two review legs.
