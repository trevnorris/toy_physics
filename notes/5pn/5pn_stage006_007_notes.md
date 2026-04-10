# 5PN continuation notes — stages 6 and 7

These two stages take the Stage-5 primitive-deformation problem and turn it into a genuine mechanism sieve.

## Stage 6: which primitive sectors are dead, and which corridor survives?

The Stage-5 continuation point was to test

- `K1 = D21 + D01/9`,
- `Xi_load = N01/N0 - D01/D0`,
- `H_even = D41 - (2/3) D21 - D01/27`.

The results are:

1. **Wall-only weak-axisymmetric anisotropy is dead.**
   With only `(dK,dM)` active,
   `D01 = dK`, `D21 = -dM`, `D41 = 0`, `N01 = 0`,
   and the full solve of `(K1, Xi_load, H_even) = 0` gives only the trivial branch
   `dK = dM = 0`.

2. **Pure BdG weak-axisymmetric anisotropy is also dead.**
   In logarithmic form with
   `x_c = δ ln c`, `x_varpi = δ ln varpi`,
   the full solve of `(K1, Xi_load, H_even) = 0` again gives only the trivial branch.

3. **BdG self-similarity kills only the load defect, not the full 5PN triplet.**
   On a wall-fixed pure-BdG branch, Stage-157 self-similarity reduces to
   `x_c = x_varpi`.
   Then `Xi_load = 0`, but generically `K1 != 0` and `H_even != 0` unless the branch is trivial.

So after the primitive sieve, the nontrivial surviving corridor is **not** wall-only or BdG-only.

## Stage 6: exact self-similarity and outgoing-load theorems

The exact Stage-157 decomposition is

`Xi_load = (delta_N - delta_K) + omega_B (delta_B - delta_K) + omega_Z (delta_Z - delta_K)`.

Equivalently, in wall-referenced defect fields,

`Xi_load = Sigma^(N) + omega_B Sigma^(B) + omega_Z Sigma^(Z)`

with the understood normalized sums over modes/ports.

So if the weak-axisymmetric branch is **statically self-similar** relative to the wall baseline,

- `Sigma^(B) = 0`,
- `Sigma^(Z) = 0`,
- `Sigma^(N) = 0`,

then automatically

`Xi_load = 0`.

Stage 158 sharpens that further. On conservative-shape-preserving branches,

`Xi_load = 2 sum_r rho_r^(N) δ ln Λ_r - δK`,

where `Λ_r = P_r / Δ_r` is the outgoing load factor.

A key no-go follows immediately:

- if all wall-normalized shapes are frozen, including `δ ln Λ_r = 0`, then
  `Xi_load = -δK`,
  so naive common self-similarity fails on any nontrivial wall-loading branch.

Therefore the outgoing sector must actively load with the wall baseline.

## Stage 6: exact surviving outgoing corridor

Stage 159 factors the outgoing load as

`Λ_r^2 / K = M_r^2 (1 + I_r)^2 / (1 - H_r)^2`

with

- `M_r = G_W / (Ω_W^2 sqrt(K))`,
- `I_r = R G_U / (Ω_U^2 G_W)`,
- `H_r = R^2 / (Ω_U^2 Ω_W^2)`.

The outgoing defect splits exactly into

`Sigma_r^(N) = 2 δ ln M_r + 2 I_r/(1+I_r) δ ln I_r + 2 H_r/(1-H_r) δ ln H_r`.

So if the interference and hybridization ratios are rigid,

- `δ ln I_r = 0`,
- `δ ln H_r = 0`,

then the whole defect collapses to the raw mixed leg:

`Sigma_r^(N) = 2 δ ln M_r`.

This yields the exact **square-root mixed-leg law**

`G_W / Ω_W^2 ∝ sqrt(K)`

as the surviving nontrivial first-order cancellation condition.

## Stage 7: one scalar amplitude controls the remaining weak-axisymmetric defect

Stage 160 then shows that on the weak-axisymmetric grouped branch,

- every microscopic outgoing slippage inherits the grouped signature
  `(1, 1/2, -1)`,
- each outgoing port collapses to one scalar amplitude

`σ_r = 2 m_r + 2 I_r/(1+I_r) i_r + 2 H_r/(1-H_r) h_r`,

with

- `m_r = g_W - o_W - κ1/2`,
- `i_r = r + g_U - o_U - g_W`,
- `h_r = 2 r - o_U - o_W`.

The full remaining grouped defect is therefore one weighted scalar

`Xi_1 = sum_r rho_r^(N) σ_r`.

And the grouped pattern is fixed exactly:

- `Xi^(20) = ε Xi_1`,
- `Xi^(21) = ε Xi_1 / 2`,
- `Xi^(22) = - ε Xi_1`.

So its grouped anisotropy automatically obeys

`b = 3 a`.

Most importantly, the same scalar is the physical outgoing-prefactor slope:

`Xi_1 = P1 / P0`.

So after Stage 7 the remaining weak-axisymmetric grouped problem is no longer “compute all microscopic drifts.”
It is:

> compute the single microscopic amplitude `Xi_1 = P1/P0` on the actual moving-throat branch.

That is the direct continuation point.
