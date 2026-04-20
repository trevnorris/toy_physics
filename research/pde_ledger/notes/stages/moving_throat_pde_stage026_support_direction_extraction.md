# Moving-Throat PDE — Stage 26: Continuum Extraction of the Actual Support/BdG Loading Direction

## Purpose

Stage 24 and Stage 25 reduced the next bottleneck to one sharp question:

> does the physical support/BdG kernel follow the deformed mixed direction `z`, or does it stay tied to the original source direction `v`?

The minimal Stage-20 continuum operator implicitly answered that only in the most degenerate case: because the support mode `phi` coupled directly to the wall through the same overlap vector `v`, the support direction was exactly source-tied.

But that was still not the first genuinely nontrivial support operator.

The first symmetry-allowed extension is to turn on the bilinear `U/phi` coupling. Once that is done, the support direction can be extracted exactly from the continuum kernel rather than guessed.

The main result is:

1. the effective support/BdG loading remains rank-1,
2. its actual wall-basis direction is an exact rotated vector `y`,
3. the rotation is controlled by one new interference ratio `sigma_0`,
4. the source-tied closure is the exact minimal-kernel limit `sigma_0 = 0`,
5. the tracking closure is an exact codimension-one interference-match condition,
6. and the generic first extended continuum kernel lands on an **intermediate support direction**, neither source-tied nor tracking.

So the Stage-24 binary fork is now resolved more sharply:

- the **minimal** continuum kernel is source-tied,
- the **first extended** continuum kernel is generically intermediate,
- and exact tracking requires one special interference match.

---

## 1. First symmetry-allowed support extension of the continuum kernel

Keep the Stage-22 split-`U` continuum operator, but add the first local bilinear support coupling

`L_(Uphi) = + c_(Uphi) int_0^L U phi ds.`

After the same N/N and D/N mode projection used in Stages 20–22, the reduced static couplings are

`L_(eta U)   = - g_U Q.U,`

`L_(eta phi) = - g_B (v.Q) Phi,`

`L_(Uphi)    = + g_S (v.U) Phi,`

with wall basis `Q = (Q_0,Q_1)^T`, D/N overlap vector

`v = (kappa_0, kappa_1)^T,`

and the same split internal stiffnesses

`A_(U0) = K_U,`

`A_(U1) = K_U (1 + delta_U).`

So the support sector now sees the wall both

- directly through `g_B v`, and
- indirectly through `U` via `g_U g_S`.

That is the first honest continuum mechanism that can rotate the support direction away from the original source vector.

---

## 2. Exact effective support vector

Eliminate the split `U` doublet first.

The exact effective wall-to-support coupling becomes

`y = g_B v + g_U g_S D_U v,`

with

`D_U = diag( 1/K_U, 1/[K_U(1+delta_U)] ).`

So the wall-basis components are

`y_0 = kappa_0 [ g_B + g_U g_S / K_U ],`

`y_1 = kappa_1 [ g_B + g_U g_S / ( K_U (1+delta_U) ) ].`

Define the exact support-interference ratio

`sigma_0 := g_U g_S / (K_U g_B)`

or equivalently, in continuum couplings,

`sigma_0 = c_(etaU) c_(Uphi) / ( K_U c_(etaphi) ).`

Then the support vector takes the exact factorized form

`y_0 = kappa_0 g_B (1 + sigma_0),`

`y_1 = kappa_1 g_B [ 1 + sigma_0/(1+delta_U) ].`

So the actual support/BdG direction ratio is

`y_1 / y_0 = (kappa_1 / kappa_0) R_phi,`

with the exact support-direction factor

`R_phi := [ 1 + sigma_0/(1+delta_U) ] / (1 + sigma_0).`

This is the exact support analogue of the mixed-direction factor from Stage 22,

`R_U = [ 1 + rho_0/(1+delta_U) ] / (1 + rho_0).`

---

## 3. Exact support direction-splitting theorem

The support direction is collinear with the original source vector `v` iff

`R_phi = 1.`

The exact support direction-splitting invariant is

`D_phi := kappa_0 y_1 - kappa_1 y_0`

`      = - kappa_0 kappa_1 g_B sigma_0 delta_U / (1 + delta_U).`

Therefore

`D_phi = 0  <=>  sigma_0 = 0  or  delta_U = 0.`

So the exact theorem is:

> **The support/BdG loading remains source-tied only in the minimal kernel limit `sigma_0 = 0` (or in the unsplit limit `delta_U = 0`).**

That identifies the precise structural role of the minimal Stage-20 operator: it was not just simple; it was exactly the source-tied special case.

---

## 4. Exact support pole shift and actual support baseline

Eliminating `U` also renormalizes the support pole.

The exact effective support denominator is

`A_phi^(eff) = K_phi^(eff) - g_S^2 v.D_U.v,`

with

`K_phi^(eff) := K_phi + pi^2 T_phi / (4 L^2).`

Using

`sigma = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2),`

and `kappa_1^2/sigma = 2/11`, the overlap contraction is

`v.D_U.v = (sigma / K_U) [ 1 - (2/11) delta_U/(1+delta_U) ].`

So if we define the flat support-blocking ratio

`eps_phi := c_(Uphi)^2 sigma / ( K_U K_phi^(eff) ),`

then the exact split support-blocking ratio is

`eps_phi^(split) = eps_phi [ 1 - (2/11) delta_U/(1+delta_U) ].`

The physical support loading baseline on the selected wall problem is therefore

`M_supp`
`= 8 Z_phi (1 + sigma_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_phi^(split)) ],`

with

`Z_phi := c_(etaphi)^2 / ( K_eta^(eff) K_phi^(eff) ).`

So the support sector is now no longer just a direction question. Its actual baseline strength is an exact continuum output, parallel to the mixed baseline `M_mix` of Stage 22.

---

## 5. Exact tracking theorem relative to the mixed vector

Stage 22 already gave the mixed loading vector

`z_0 = kappa_0 g_W (1 + rho_0),`

`z_1 = kappa_1 g_W [ 1 + rho_0/(1+delta_U) ].`

So the support vector `y` tracks the mixed vector `z` iff `R_phi = R_U`.

The exact mixed-support collinearity invariant is

`D_(phi z) := y_0 z_1 - y_1 z_0`

`          = - delta_U kappa_0 kappa_1 g_U ( g_B g_R - g_W g_S ) / [ K_U (1+delta_U) ].`

Therefore

`y || z  <=>  g_B g_R = g_W g_S`

or equivalently

`sigma_0 = rho_0.`

So exact tracking is not generic. It is a codimension-one interference match between

- the `eta/W`–`U/W` mixed lane, and
- the `eta/phi`–`U/phi` support lane.

---

## 6. Exact mismatch formula and branch interpretation

The support and mixed direction factors differ by

`R_phi - R_U`
`= delta_U (rho_0 - sigma_0)`
`  / [ (1+delta_U)(1+rho_0)(1+sigma_0) ].`

So the sign of the direction mismatch is set by the sign of `rho_0 - sigma_0`.

This yields the clean continuum interpretation:

- `sigma_0 = 0`  ->  exact source-tied closure,
- `sigma_0 = rho_0`  ->  exact tracking closure,
- generic `sigma_0`  ->  intermediate support direction.

So the physical kernel does **not** generically choose between the two Stage-24 extremes. It generically interpolates between them.

---

## 7. Best current theorem statement after Stage 26

The support-direction bottleneck is no longer open in principle.

### What is now exact

- the first symmetry-allowed continuum extension of the support sector,
- the actual effective support vector `y`,
- the exact support-direction factor `R_phi`,
- the exact support direction-splitting invariant `D_phi`,
- the exact support-blocking ratio `eps_phi^(split)`,
- the exact physical support baseline `M_supp`,
- and the exact tracking condition `g_B g_R = g_W g_S`.

### What this means physically

- The **minimal** continuum kernel lands on the **source-tied** closure.
- The **first extended** continuum kernel is generically **intermediate**.
- Exact **tracking** requires a special interference match.

So the next theorem step is no longer “guess the support direction.”
It is:

> insert the exact continuum-selected quantities `(M_mix, M_supp, R_U, R_phi)` into the Stage-24/25 selected-branch formulas and determine the physical selected branch.
