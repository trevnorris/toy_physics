# Moving-Throat PDE — Stage 14: Selected-Branch Reachability, Monotonicity, and the Stable-Side Normalization Window

## Purpose

Stage 13 translated the selected lower quadrupole mode into the same normalized-response language used by the earlier grouped-real `P2` work and reduced the remaining theorem gap to one exact spectral ratio,

`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-`.

The next question is now unavoidable:

> on the stable branch, is this selected prefactor a wildly oscillatory function of the loading, or is the normalization target a clean one-parameter crossing problem?

The answer is much better than it had any right to be.
For the exact `2 x 2` loaded wall problem, the selected prefactor is **strictly monotone increasing** on the stable side.
It starts at the flat-branch value, diverges at the softening threshold, and therefore hits any target above its starting value exactly once before instability.

So the normalization question is no longer an arbitrary tuning problem.
It is a stable-side spectral crossing theorem.

---

## 1. Exact selected overlap derivative

Keep the Stage-13 notation

`A = K_0 - Xi_0`,

`B = K_1 - Xi_0 = A + Delta K_ax`,

`s_- = (v.e_-)^2`,

`lambda_- = ( A + B - alpha_0 sigma - R ) / 2`,

`R = sqrt( (Delta K_ax + alpha_0 delta_kappa)^2 + 4 alpha_0^2 KappaProd )`.

The first new exact identity is

`d s_- / d alpha_0 = 2 Delta K_ax^2 KappaProd / R^3`.

Because

`Delta K_ax > 0`,

`KappaProd = kappa_0^2 kappa_1^2 > 0`,

and `R > 0`,

this derivative is strictly positive on the whole stable branch.

So the selected overlap with the loading vector does not meander.
It grows monotonically as the directional loading is increased.

---

## 2. Exact monotonicity of the selected static prefactor

The selected prefactor is

`P_{0,-} = beta_0 s_- / lambda_-`.

Using the Hellmann–Feynman identity

`d lambda_- / d alpha_0 = - s_-`,

its derivative is exactly

`d P_{0,-} / d alpha_0`
`= beta_0 [ (d s_- / d alpha_0) lambda_- + s_-^2 ] / lambda_-^2`.

On the stable branch,

- `beta_0 > 0`,
- `lambda_- > 0`,
- `d s_- / d alpha_0 > 0`,
- `s_-^2 > 0`.

Therefore

`d P_{0,-} / d alpha_0 > 0`.

So the selected static prefactor is strictly monotone increasing all the way up to the softening threshold.

This is the core theorem of the stage.
It turns the selected-branch normalization problem from an opaque reduced-parameter search into a one-dimensional ordered crossing problem.

---

## 3. Starting value at zero loading

At `alpha_0 = 0`, the lower wall mode is just the flat `u_0` branch, so

`lambda_-(0) = A = K_0 - Xi_0`,

`s_-(0) = kappa_0^2`.

Therefore the selected prefactor starts at

`P_{0,-}(0) = beta_0 kappa_0^2 / (K_0 - Xi_0)`.

This is the exact stable-side starting value against which the universal target must be compared.

---

## 4. Exact softening threshold and divergence of `P_{0,-}`

The exact determinant identity is

`lambda_- lambda_+ = A B - alpha_0 ( B kappa_0^2 + A kappa_1^2 )`.

So the refined softening threshold is

`alpha_crit = A B / ( B kappa_0^2 + A kappa_1^2 )`.

At that point,

`lambda_-(alpha_crit) = 0`,

while the selected overlap `s_-` stays finite and positive.

Because

`P_{0,-} = beta_0 s_- / lambda_-`,

it follows that

`P_{0,-} -> +infinity`

as `alpha_0 -> alpha_crit^-` from the stable side.

So the stable branch spans the full interval

`P_{0,-}(0) <= P_{0,-} < +infinity`.

---

## 5. Unique stable-side crossing theorem

Combine the last two sections.

- `P_{0,-}(alpha_0)` is continuous on `0 <= alpha_0 < alpha_crit`,
- strictly increasing on that interval,
- starts at `P_{0,-}(0)`,
- and diverges at `alpha_crit^-`.

Therefore:

> for every target value `P_target > P_{0,-}(0)`, there exists a **unique** stable-side loading `alpha_*` with `0 < alpha_* < alpha_crit` such that `P_{0,-}(alpha_*) = P_target`.

Applied to the 2.5PN normalization target,

`P_target = N_Q^(target) / mhat_-^2`

with

`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`,

this means the selected-branch quadrupole condition is not a discrete impossibility.
It is a unique crossing problem on the stable branch.

The only remaining question is whether the physical moving-throat branch puts the system in the right region of parameter space so that the crossing occurs on the natural passive/outgoing branch before other approximations fail.

---

## 6. Best current theorem gate after Stage 14

After Stage 14, the normalization bottleneck can be stated in its sharpest current form.

The selected-branch 2.5PN problem is no longer:

- “derive some outgoing coefficient,”
- or “guess the right axial profile,”
- or “hope the branch lands near the right value.”

It is now the following exact test:

1. compute the physical stable-branch data `(Xi_0, beta_0, alpha_0)` from the coupled moving-throat operator,
2. evaluate the selected prefactor
   `P_{0,-}(alpha_0) = beta_0 (v.e_-)^2 / lambda_-`,
3. compare it to the target `N_Q^(target) / mhat_-^2`,
4. and check whether the resulting stable-side crossing sits on the natural passive/outgoing branch with `alpha_0 < alpha_crit`.

So the theorem gap is now a controlled spectral-placement problem, not an algebraic unknown.
