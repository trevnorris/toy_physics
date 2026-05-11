# Moving-Throat PDE — Stage 030: Selected-Mode Normalized Response, Exact Static Prefactor, and the Selected-Branch Quadrupole Target

## Purpose

Stage 12 reached the first exact operator-selected odd quadrupole coefficient,

`delta D_-^(odd)(omega) = - i beta_5 (v.e_-)^2 omega^5 + O(omega^7)`.

That was the right operator-level result, but it still stopped one step short of the language used by the earlier grouped-real `P2` normalization bridge.

The next missing step is therefore very specific:

> translate the selected lower wall mode into the same normalized-response convention used in Stages 5 and 8, identify the exact selected-branch static prefactor, and write the surviving `2G/(5 c^5)` test as a single spectral condition on the selected conservative eigenvalue.

This stage does exactly that.

The main result is simple and important.
If the conservative selected lower wall eigenvalue is written as `lambda_-(omega)`, then the normalized selected-mode response has the low-frequency form

`Y_-(omega) = 1 + u_{2,-} omega^2 + u_{4,-} omega^4 + i Gamma_{5,-} omega^5 + ...`,

with

`Gamma_{5,-} = Gamma_2^(port) * P_{0,-}`,

where the exact selected static prefactor is

`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-(0)`

and, equivalently,

`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.

So the selected-branch normalization problem is no longer “some outgoing coefficient somewhere.”
It is one exact ratio of a conservative overlap factor to the conservative selected wall stiffness.

---

## 1. Selected conservative eigenvalue and notation

Keep the Stage-12 conservative wall matrix

`K_eff^(0) = [[A, 0], [0, B]] - alpha_0 v v^T`,

with

`A = K_0 - Xi_0`,

`B = K_1 - Xi_0 = A + Delta K_ax`,

`v = (kappa_0, kappa_1)^T`,

`sigma = kappa_0^2 + kappa_1^2`,

`delta_kappa = kappa_0^2 - kappa_1^2`,

`KappaProd = kappa_0^2 kappa_1^2`.

The exact selected lower conservative eigenvalue is

`lambda_-`
`= ( A + B - alpha_0 sigma - R ) / 2`,

where

`R = sqrt( (Delta K_ax + alpha_0 delta_kappa)^2 + 4 alpha_0^2 KappaProd )`.

The corresponding upper branch is

`lambda_+ = ( A + B - alpha_0 sigma + R ) / 2`.

Throughout this stage,

`D_{-0} = lambda_-`

is the selected conservative wall stiffness at zero frequency.

---

## 2. Exact selected overlap from Hellmann–Feynman

Stage 12 already established the exact Hellmann–Feynman identity

`(v.e_-)^2 = - d lambda_- / d alpha_0`.

Writing

`s_- = (v.e_-)^2`,

the explicit closed form is

`s_-`
`= (1/2) [ sigma`
`          + ( (Delta K_ax + alpha_0 delta_kappa) delta_kappa + 4 alpha_0 KappaProd ) / R ]`.

This formula immediately passes the expected checks:

- weak loading `alpha_0 -> 0`:
  `s_- -> kappa_0^2`,
- strong loading `alpha_0 -> +infinity`:
  `s_- -> sigma = |v|^2`.

So the selected overlap interpolates smoothly from the flat branch to the Stage-10 max-coupling branch.

---

## 3. From the selected operator to the normalized selected response

Write the selected wall operator in the same low-frequency form used in Stage 022,

`D_-(omega)`
`= D_{-0} + D_{-2} omega^2 + D_{-4} omega^4 - i C_{5,-} omega^5 + O(omega^6)`.

Here the Stage-12 odd coefficient is

`C_{5,-} = beta_5 s_-`.

The normalized selected-mode response is therefore

`Y_-(omega) = D_{-0} / D_-(omega)`
`          = 1 + u_{2,-} omega^2 + u_{4,-} omega^4 + i Gamma_{5,-} omega^5 + O(omega^6)`.

Exactly as in Stage 022, the even coefficients are

`u_{2,-} = - D_{-2} / D_{-0}`,

`u_{4,-} = ( D_{-2}^2 - D_{-0} D_{-4} ) / D_{-0}^2`.

The selected odd coefficient is

`Gamma_{5,-} = C_{5,-} / D_{-0}`
`           = beta_5 s_- / lambda_-`.

So the Stage-12 operator-level result is now fully translated into the Stage-022 normalized-response language.

---

## 4. Exact selected static prefactor `P_{0,-}`

Using the compact outgoing fingerprint

`beta_5 = Gamma_2^(port) beta_0`,

`Gamma_2^(port) = a^5 / (27 c_s^5)`,

we obtain

`Gamma_{5,-} = Gamma_2^(port) P_{0,-}`

with the exact selected static prefactor

`P_{0,-} = beta_0 s_- / lambda_-`.

This is the selected-mode analogue of the isotropic Stage-022 prefactor `P_0 = N_0 / D_0`.

Using the Hellmann–Feynman identity, the same prefactor can be rewritten as

`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.

This is a very useful compression of the theorem gap.
It says the passive/outgoing quadrupole normalization is controlled by how sensitively the conservative selected wall stiffness responds to the directional loading parameter.

---

## 5. The selected-branch `2G/(5 c^5)` target

The invariant 2.5PN quadrupole condition is still

`mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5)`.

Substituting the selected prefactor gives

`mhat_-^2 * Gamma_2^(port) * P_{0,-} = 2 G / (5 c^5)`.

Using

`Gamma_2^(port) = a^5 / (27 c_s^5)`

this is equivalent to

`mhat_-^2 P_{0,-} = 54 G c_s^5 / (5 a^5 c^5)`.

So the selected branch is required to hit exactly the same normalization stack as the isotropic Stage-022/8 branch, but now with the selected-mode prefactor in place of the old isotropic lane prefactor.

Equivalently, the target becomes a direct conservative spectral condition:

`lambda_- = mhat_-^2 beta_0 s_- / N_Q^(target)`,

where

`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`.

On the natural source-map branch,

`mhat_- = 1 + O(a^2/r^2)`,

so at leading point-particle order,

`lambda_- = beta_0 s_- * 5 a^5 c^5 / (54 G c_s^5)`.

This is the selected-mode generalization of the Stage-9 stiffness test.
The difference is that the “stiffness” is no longer a hand-picked wall constant.
It is the exact selected conservative eigenvalue of the loaded wall/BdG/Maxwell/mixed operator.

---

## 6. Exact determinant identity and a useful spectral rewrite

Because the selected mode comes from a rank-1 update of a diagonal `2 x 2` wall block, its determinant factor is exact:

`lambda_- lambda_+ = A B - alpha_0 ( B kappa_0^2 + A kappa_1^2 )`.

So the exact softening threshold is

`alpha_crit = A B / ( B kappa_0^2 + A kappa_1^2 )`,

which is just the Stage-12 refined threshold written in the compact `A,B` notation.

This means the selected-branch quadrupole target can also be read as a condition on **how close** the physical branch sits to the softening surface.
If `lambda_-` is too large, the selected outgoing coefficient is too small.
If `lambda_-` is too small, the branch is close to instability.

So the remaining PDE question is now sharply spectral:

> what conservative selected eigenvalue `lambda_-` and selected overlap `s_-` does the physical moving-throat branch produce, and does their ratio land on the universal target before the wall mode softens?

---

## 7. Best current summary after Stage 13

The theorem gap has narrowed again.

We no longer need to talk vaguely about an operator-level odd coefficient.
The selected-branch 2.5PN problem is now the single exact quantity

`P_{0,-} = beta_0 (v.e_-)^2 / lambda_-`

or, equivalently,

`P_{0,-} = - beta_0 d(ln lambda_-)/d alpha_0`.

And the full selected-branch normalization theorem is simply

`mhat_-^2 P_{0,-} = 54 G c_s^5 / (5 a^5 c^5)`.

That is the cleanest selected-mode formulation we have reached so far.
The next honest step is no longer to invent another reduced parameter.
It is to determine whether the physical stable branch can actually reach this target in a controlled way.
