# Moving-Throat PDE — Stage 10: First Nonconstant Finite-Throat Wall/Brane Family, Exact Overlap Law, and the Profile-Selection Theorem Gate

## Purpose

Stage 9 showed that on the simplest finite-throat branch — constant N/N wall and brane-like zero modes coupled to the lowest D/N half-wave — the whole radial/axial overlap problem collapses to one exact geometric constant,

`kappa_0 = 2 sqrt(2) / pi`.

That was useful, but it still left one obvious concern:

> was the Stage-9 normalization equation an artifact of choosing the completely flat axial zero mode for the wall and brane-like internal coordinate?

So the next honest step is to replace the constant profile by the **first nontrivial finite-throat profile family** and see what survives.

The cleanest way to do that is to stay in the exact N/N and D/N finite-throat bases and turn on the first N/N axial excitation. The main result is that the Stage-9 branch does survive, but it deforms in a very specific way:

- the whole branch is controlled by a single profile-dependent overlap
  `kappa(theta)`,
- the wall stiffness picks up the expected axial-gradient penalty,
- there is an exact **blind angle** where the D/N support/mixed route is shut off,
- and there is also an exact **max-coupling angle** where the wall/support overlap is stronger than on the constant branch.

So after this stage the next theorem gate is no longer “does any nonconstant profile destroy the branch?”
It is:

> which axial profile is selected by the real moving-throat eigenproblem, and does it lie near the Stage-9 branch, the max-coupling branch, or the blind-angle no-go branch?

---

## 1. First nonconstant finite-throat family

Keep the finite throat interval

`s in [0,L]`.

Use the exact N/N basis for the wall and brane-like internal coordinate,

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`.

These satisfy

`-u_0'' = 0`,

`-u_1'' = (pi^2 / L^2) u_1`,

with

`u_0'(0)=u_0'(L)=0`,

`u_1'(0)=u_1'(L)=0`,

and are orthonormal on `[0,L]`.

Keep the same exact D/N half-wave branch for the trapped support and mixed channels,

`f_0(s) = sqrt(2/L) sin(pi s / (2L))`.

The minimal coherent nonconstant family is then

`chi_theta(s) = cos(theta) u_0(s) + sin(theta) u_1(s)`.

In this stage we use that same profile family for both

- the wall quadrupole axial shape, and
- the brane-like internal gauge/support coordinate,

because it is the smallest nonconstant replacement of the old “constant wall/brane zero-mode” branch that still preserves one common coherent axial deformation parameter.

So the branch choice is

`chi_eta = chi_theta`,

`u = chi_theta`,

`phi = f_0`,

`w = f_0`.

This has one immediate exact benefit:

`int_0^L ds chi_eta u = 1`.

So the direct wall/brane overlap does not degrade just because we turned on the first N/N profile correction.

---

## 2. Exact overlap law with the D/N half-wave

Define the two basic constants

`kappa_0 = int_0^L ds u_0 f_0 = 2 sqrt(2) / pi`,

`kappa_1 = int_0^L ds u_1 f_0 = -4 / (3 pi)`.

Then the coherent profile family has the exact overlap

`kappa(theta) = int_0^L ds chi_theta f_0`
`             = kappa_0 cos(theta) + kappa_1 sin(theta)`
`             = 2 ( -2 sin(theta) + 3 sqrt(2) cos(theta) ) / (3 pi)`.

So the whole nonconstant branch still collapses to one profile function, but now it is angle dependent rather than constant.

A useful amplitude form is

`kappa(theta) = rho cos(theta - theta_max)`,

with

`rho = sqrt(kappa_0^2 + kappa_1^2) = 2 sqrt(22) / (3 pi)`.

That means the family can reach a maximal overlap magnitude

`|kappa|_max = 2 sqrt(22) / (3 pi) ~= 0.994031`.

This is larger than the constant-branch value

`kappa_0 = 2 sqrt(2) / pi ~= 0.900316`.

So the flat Stage-9 branch is not the strongest-coupled member of the first nonconstant family.

### 2.1 Exact blind angle

The support/mixed D/N route shuts off when `kappa(theta)=0`, i.e.

`tan(theta_blind) = 3 sqrt(2) / 2`.

So

`theta_blind = arctan( 3 sqrt(2) / 2 )`.

At this angle the coherent N/N profile is orthogonal to the lowest D/N half-wave.

That is the first exact profile-selection no-go theorem of the moving-throat PDE program:

> if the actual wall/brane profile is driven onto the blind angle, the Stage-9 support/mixed quadrupole branch cannot realize the outgoing quadrupole normalization target.

### 2.2 Exact max-coupling angle

The positive overlap is maximized at

`tan(theta_max) = kappa_1 / kappa_0 = - sqrt(2) / 3`,

so

`theta_max = - arctan( sqrt(2) / 3 )`.

At that angle,

`sin^2(theta_max) = 2 / 11`,

and

`kappa(theta_max) = 2 sqrt(22) / (3 pi)`.

So the strongest-coupled member of the first N/N family is obtained by a small **negative** admixture of the first excited N/N mode, and the price paid is a fixed axial-gradient occupancy of `2/11`.

---

## 3. Exact wall stiffness on the nonconstant family

The linear wall operator carried from the earlier stages is

`G_eta = -T_w d_s^2 + K_eta + 6 T_Omega`

in the `l=2` quadrupole sector.

Evaluating it on the coherent family gives

`K_geo(theta) = int_0^L ds chi_theta G_eta chi_theta`
`             = K_eta + 6 T_Omega + T_w (pi^2 / L^2) sin^2(theta)`.

So the first nonconstant family changes the Stage-9 wall branch in exactly the way one would expect physically:

- the constant branch is recovered at `theta=0`,
- the first excited N/N admixture adds an axial-gradient penalty,
- and the gradient cost is quadratic in the profile angle.

At the max-coupling point,

`K_geo(theta_max) = K_eta + 6 T_Omega + 2 T_w pi^2 / (11 L^2)`.

So the first nonconstant family makes the trade-off completely explicit:

- larger D/N overlap,
- but larger axial wall stiffness.

---

## 4. Exact substitution into the Stage-8/9 minimal isotropic module

Using the same reduced couplings as before,

`C   = lambda_B int chi_eta phi`,

`G_U = lambda_U int chi_eta u`,

`G_W = lambda_W int chi_eta w`,

`R   = lambda_R int u w`,

the coherent nonconstant family gives

`C   = lambda_B kappa(theta)`,

`G_U = lambda_U`,

`G_W = lambda_W kappa(theta)`,

`R   = lambda_R kappa(theta)`.

So the Stage-8/9 quantities become

`Delta(theta) = Omega_U^2 Omega_W^2 - lambda_R^2 kappa(theta)^2`,

`Q(theta) = lambda_U^2 Omega_W^2`
`         + 2 lambda_U lambda_W lambda_R kappa(theta)^2`
`         + lambda_W^2 Omega_U^2 kappa(theta)^2`,

`P(theta) = kappa(theta) ( Omega_U^2 lambda_W + lambda_R lambda_U )`.

Therefore

`B_0(theta) = lambda_B^2 kappa(theta)^2 / varpi^2`,

`Z_0(theta) = Q(theta) / Delta(theta)`,

`N_0(theta) = kappa(theta)^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta(theta)^2`,

`D_0(theta) = K_geo(theta) - B_0(theta) - Z_0(theta)`.

This is the cleanest possible outcome: the whole Stage-9 branch survives with the single replacement

`kappa_0  ->  kappa(theta)`

and the single geometry lift

`K_eta + 6 T_Omega  ->  K_eta + 6 T_Omega + T_w pi^2 sin^2(theta) / L^2`.

So the Stage-9 branch was not an algebraic accident of the constant profile. It is the `theta=0` member of a whole exact finite-throat family.

---

## 5. Exact branch-level normalization condition

Substituting the nonconstant family into the Stage-8 target gives

`mhat_rad^2 kappa(theta)^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`/ [ Delta(theta) ( K_geo(theta) Delta(theta)`
`                - Delta(theta) lambda_B^2 kappa(theta)^2 / varpi^2`
`                - Q(theta) ) ]`
`= 54 G c_s^5 / (5 a^5 c^5)`.

Solving for the required geometry-side stiffness gives

`K_req(theta)`
` = lambda_B^2 kappa(theta)^2 / varpi^2`
` + Q(theta) / Delta(theta)`
` + mhat_rad^2 kappa(theta)^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`   / [ (54 G c_s^5 / (5 a^5 c^5)) Delta(theta)^2 ]`.

So the full theorem gate on this family is now

`K_eta + 6 T_Omega + T_w pi^2 sin^2(theta) / L^2 = K_req(theta)`.

This is sharper than the Stage-9 gate, not weaker.

Stage 9 asked whether the constant wall stiffness matched one algebraic target.
Stage 10 asks whether the actual wall eigenprofile chosen by the PDE lands on an angle `theta` for which the geometry-side stiffness exactly equals the profile-dressed support/mixed normalization load.

---

## 6. Exact consequences of the first nonconstant family

### 6.1 Stage 9 is recovered exactly

At `theta=0`,

`chi_theta -> u_0`,

`kappa(theta) -> kappa_0 = 2 sqrt(2) / pi`,

`K_geo(theta) -> K_eta + 6 T_Omega`.

So the entire Stage-9 branch is recovered exactly.

### 6.2 Blind-angle no-go theorem

At the blind angle,

`kappa(theta_blind)=0`.

Then

`C = G_W = R = 0`,

`P(theta_blind)=0`,

`N_0(theta_blind)=0`.

Therefore the left-hand side of the normalization equation vanishes, while the GR target on the right-hand side remains strictly positive.

So the blind-angle branch is an exact no-go for the outgoing quadrupole normalization bridge.

This is important: it means not every nonconstant wall profile is acceptable. The PDE must choose a profile that keeps nonzero projection onto the D/N support/mixed half-wave.

### 6.3 Max-coupling branch is allowed but not free

At `theta=theta_max`, the D/N overlap is maximized,

`kappa = 2 sqrt(22) / (3 pi)`.

But the wall stiffness also increases to

`K_geo = K_eta + 6 T_Omega + 2 T_w pi^2 / (11 L^2)`.

So the first nonconstant family exposes an exact trade-off:

- stronger coupling to the D/N support/mixed branch,
- but a larger geometry-side axial-gradient cost.

That is the first place where profile selection becomes a real dynamical competition rather than just an overlap calculation.

---

## 7. Best current summary after Stage 10

The first nonconstant finite-throat family does **not** destroy the Stage-9 route.
Instead it does something better: it turns Stage 9 into one member of a whole exact finite-throat profile family and isolates the first genuine profile-selection theorem gate.

The remaining question is no longer

> do nonconstant wall profiles matter?

They do.

The sharper question is now

> what axial profile does the real moving-throat eigenproblem actually select, and does that profile keep enough overlap with the D/N support/mixed half-wave to satisfy the outgoing quadrupole normalization equation?

The next honest derivation step is therefore to stop choosing the profile family by hand and solve the first actual axial wall eigenproblem with the matter/gauge loading included, so that the profile angle `theta` is no longer free but is an output of the reduced moving-throat operator itself.
