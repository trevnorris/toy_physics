# S10 — my predictions, written BEFORE the scripts landed

⚠ Pre-registered so the comparison means something. These are the orchestrator's hand-derived
expectations from the walk with the user (moves 1–3). ⛔ None is verified. If a script disagrees,
**the script wins and the disagreement is the finding.**

## Main computation

| D | root | predicted nullity | predicted orientation |
|---|---|---|---|
| 2 | `ω² = 0` | 1 | parallel to k |
| 2 | `ω² = (μ_R/ρ_br)k²` | 1 | perpendicular to k |
| 3 | `ω² = 0` | 1 | parallel to k |
| 3 | `ω² = (μ_R/ρ_br)k²` | **2** | perpendicular to k |
| 4 | `ω² = 0` | 1 | parallel |
| 4 | `ω² = (μ_R/ρ_br)k²` | 3 | perpendicular |
| 5 | `ω² = 0` | 1 | parallel |
| 5 | `ω² = (μ_R/ρ_br)k²` | 4 | perpendicular |

⇒ predicted pattern: **longitudinal nullity 1 in every D; transverse nullity D−1**; sum = D.

⭐ **The claim that matters:** the count is `D_brane − 1`, so at D=3 it is 2. ⛔ It is NOT the
codimension (4−3=1). If the table shows transverse nullity tracking D−1, the "not codimension" point is
a computed result rather than my prose.

## Form control — gradient-squared stiffness `Σ(∂_i u_j)²`

Predicted: the dynamical matrix becomes isotropic (`∝ μ_R k² I`), so **one** root `ω² = (μ_R/ρ_br)k²`
with nullity **D** — the longitudinal mode propagates too, and the `ω²=0` root disappears entirely.
⇒ The medium resists every deformation, not only twisting.

## Coefficient control

Predicted: roots rescale; **nullities and orientations unchanged**. ⚠ This is exactly why a coefficient
control cannot test this step — the count is what S10 claims, and a coefficient never moves it.

## What I am least sure of

- whether the antisymmetric-derivative stiffness reduces to `|∇×u|²` at D=3 with the factor I assumed
  (the `1/2` in `(1/2)Σ_{i,j}(∂_i u_j − ∂_j u_i)²`);
- whether D=2 behaves in kind with the rest — "curl" in 2D is a scalar, and the antisymmetric derivative
  has only one independent component there, so that row may be degenerate in a way the others are not.
