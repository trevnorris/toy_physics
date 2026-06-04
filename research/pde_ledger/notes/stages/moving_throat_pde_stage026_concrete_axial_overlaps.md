# Moving-Throat PDE — Stage 026: Concrete Finite-Throat Axial Modes, Exact Overlap Constants, and the First Branch-Level Normalization Test

## Purpose

Stage 025 reduced the remaining theorem gap to one scalar expression in the radial/axial overlap data,

`mhat_rad^2 P^2 / [ Delta ( K Delta - Delta C^2 / varpi^2 - Q ) ] = 54 G c_s^5 / (5 a^5 c^5)`.

But the quantities

- `I_(eta,phi)`,
- `I_(eta,u)`,
- `I_(eta,w)`,
- `I_(u,w)`

were still formal integrals.

The next honest step is therefore to pick a **concrete finite-throat axial mode problem**, compute those overlaps exactly, and substitute them into the Stage-025 formula.

The main result of this stage is that on the first natural mixed axial branch the whole overlap problem collapses to a single geometric constant

`kappa = 2 sqrt(2) / pi`.

That turns the Stage-025 normalization target into one explicit algebraic relation among

- the constant wall quadrupole stiffness,
- the BdG support coupling,
- the brane-like internal gauge frequency,
- the mixed-channel frequency,
- and the three reduced couplings `(lambda_B, lambda_U, lambda_W, lambda_R)`.

So after this stage the remaining gap is not “unknown overlaps in general.”
It is the value of a small number of branch parameters on one concrete finite-throat mode family.

---

## 1. Concrete finite-throat axial mode problem

Take the finite throat interval

`s in [0,L]`.

Use the flat axial measure

`mu_s(s) = 1`.

### 1.1 Constant N/N zero mode

For collective wall shape and the brane-like internal gauge coordinate, take the lowest Neumann/Neumann mode

`u_0(s) = 1 / sqrt(L)`.

It obeys

`-u_0'' = 0`,

`u_0'(0) = u_0'(L) = 0`,

and is normalized by

`int_0^L ds u_0^2 = 1`.

This is the simplest axial profile for a throat deformation or brane-like zero mode that remains nonzero at the mouth.

### 1.2 D/N half-wave ladder

For the trapped support channel and the mixed `A_w / F_(mu w) / J^w` channel, use the exact finite-throat D/N ladder already isolated in the frozen-wall benchmark:

`f_n'' + k_n^2 f_n = 0`,

`f_n(0) = 0`,

`f_n'(L) = 0`,

with normalized solutions

`f_n(s) = sqrt(2/L) sin( (n + 1/2) pi s / L )`,

`k_n = (n + 1/2) pi / L`.

So the lowest trapped axial support/mixed profile is

`f_0(s) = sqrt(2/L) sin( pi s / (2L) )`.

The associated support and mixed frequencies are taken to be

`varpi_n^2 = M_phi^2 + c_phi^2 k_n^2`,

`Omega_(W,n)^2 = M_W^2 + c_W^2 k_n^2`.

On the minimal branch we keep only `n=0`.

### 1.3 Wall and brane-like profile choice on the first concrete branch

The first concrete isotropic branch is therefore

- wall quadrupole axial profile: `chi_eta = u_0`,
- BdG support axial profile: `phi = f_0`,
- brane-like internal gauge profile: `u = u_0`,
- mixed-channel axial profile: `w = f_0`.

On that branch we write

`varpi = varpi_0`,

`Omega_W = Omega_(W,0)`,

and keep `Omega_U` as the reduced internal restoring frequency of the brane-like zero mode `u_0` (so its frequency need not come from an axial gradient term).

This is the simplest branch that

- keeps the wall deformation nonzero at the mouth,
- reuses the exact D/N finite-throat support ladder,
- and keeps the brane-like gauge coordinate on the natural zero mode.

---

## 2. Exact axial overlap constants

The exact overlap of the constant zero mode with the D/N ladder is

`kappa_n = int_0^L ds u_0(s) f_n(s) = sqrt(2) / ( (n + 1/2) pi )`.

So the lowest branch value is

`kappa = kappa_0 = 2 sqrt(2) / pi`.

Numerically,

`kappa ~= 0.900316316`.

The needed minimal-branch overlap integrals are therefore

`I_(eta,phi) = int_0^L ds u_0 f_0 = kappa`,

`I_(eta,u)   = int_0^L ds u_0 u_0 = 1`,

`I_(eta,w)   = int_0^L ds u_0 f_0 = kappa`,

`I_(u,w)     = int_0^L ds u_0 f_0 = kappa`.

So on this concrete branch all radial/axial overlap information collapses to the single geometric constant `kappa`.

There is also an exact axial selection rule lurking here:

- the constant wall/brane mode couples to the `n`th D/N throat wave with strength `kappa_n ~ 1/(n+1/2)`,
- so the lowest trapped half-wave `n=0` is automatically the dominant branch,
- and higher support/mixed D/N levels are axially suppressed.

That is the first real finite-throat hierarchy statement for the overlap problem.

---

## 3. Explicit minimal-branch coefficients

Define the reduced couplings exactly as in Stage 025,

`C   = lambda_B I_(eta,phi)`,

`G_U = lambda_U I_(eta,u)`,

`G_W = lambda_W I_(eta,w)`,

`R   = lambda_R I_(u,w)`.

On the concrete branch this becomes

`C   = kappa lambda_B`,

`G_U = lambda_U`,

`G_W = kappa lambda_W`,

`R   = kappa lambda_R`.

So the minimal-branch Stage-025 quantities are

`Delta = Omega_U^2 Omega_W^2 - kappa^2 lambda_R^2`,

`Q = lambda_U^2 Omega_W^2 + 2 kappa^2 lambda_U lambda_W lambda_R + kappa^2 lambda_W^2 Omega_U^2`,

`P = kappa ( Omega_U^2 lambda_W + lambda_R lambda_U )`.

The BdG softening becomes

`X = C^2 / varpi^2 = kappa^2 lambda_B^2 / varpi^2`.

Therefore

`B_0 = kappa^2 lambda_B^2 / varpi^2`,

`Z_0 = Q / Delta`,

`N_0 = kappa^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta^2`.

The conservative wall coefficient is

`D_0 = K - kappa^2 lambda_B^2 / varpi^2 - Q / Delta`.

---

## 4. First branch-level normalization test

Substituting the concrete overlaps into the Stage-025 target gives

`mhat_rad^2 kappa^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`/ [ Delta ( K Delta - Delta kappa^2 lambda_B^2 / varpi^2 - Q ) ]`
`= 54 G c_s^5 / (5 a^5 c^5)`.

This is the first fully explicit branch-level normalization equation derived from an actual finite-throat mode problem.

It can be solved exactly for the required wall stiffness:

`K_req = kappa^2 lambda_B^2 / varpi^2 + Q / Delta`
`      + mhat_rad^2 kappa^2 ( Omega_U^2 lambda_W + lambda_R lambda_U )^2`
`        / [ (54 G c_s^5 / (5 a^5 c^5)) Delta^2 ]`.

So on this branch the GR quadrupole target is equivalent to one concrete statement:

> the constant wall quadrupole stiffness must equal the support softening plus the conservative Maxwell/mixed self-load plus one explicit outgoing-normalization load.

That is a much sharper theorem target than we had before.

---

## 5. Interpreting the wall stiffness on this branch

On the constant wall profile `chi_eta = u_0`, the axial gradient term vanishes. In the same modal normalization used in the earlier wall reduction, the bare quadrupole wall stiffness is therefore

`K = K_eta + 6 T_Omega`.

This is not a new wall-side derivation in Stage 026. It is the `l=2`
specialization of the Stage-1 modal operator, equivalently the one-mode grouped
`P2` reduction written explicitly in Stage 2.

So the normalization test can be rewritten directly as

`K_eta + 6 T_Omega = K_req`.

This is important for the roadmap because the conservative 3PN program already isolated a separate geometry lane, and the Stage-025 gap still involved the unsolved wall stiffness. On the present branch, the higher-order quadrupole normalization target directly constrains that same geometry-side combination.

So the moving-throat PDE job has tightened again.
It is not to produce an arbitrary number called `K`.
It is to determine whether the actual quadrupolar wall stiffness of the throat satisfies the explicit algebraic relation above after support and internal gauge loading are included.

---

## 6. What is exact and what is still conditional here

### Exact in this stage

The following are exact on the chosen finite-throat branch:

- the N/N zero mode `u_0`,
- the D/N ladder `f_n`,
- the overlap law `kappa_n = sqrt(2) / ( (n+1/2) pi )`,
- the minimal-branch value `kappa = 2 sqrt(2) / pi`,
- and the fully substituted Stage-025 normalization equation.

### Still conditional

What is **not** yet fixed from first principles are the branch parameters

- `lambda_B`,
- `lambda_U`,
- `lambda_W`,
- `lambda_R`,
- `Omega_U`,
- `Omega_W`,
- `varpi`,
- and `mhat_rad`.

So this is not yet the final normalization theorem.
But it is no longer an abstract overlap problem either.
It is one explicit finite-throat branch equation in a short list of reduced parameters.

---

## 7. Best current summary after Stage 026

The overlap problem is now materially tighter.

Before this stage, the remaining gap was

> compute the radial/axial overlap integrals on the true branch.

After this stage, the first concrete finite-throat branch shows that

- the overlap problem can collapse to one exact geometric constant,
- the lowest D/N half-wave is automatically the dominant support/mixed branch,
- and the quadrupole normalization target can be rewritten as an explicit condition on the wall quadrupole stiffness.

So the next honest step is no longer to “invent more overlap algebra.”
It is to decide whether the real moving-throat branch is better approximated by this N/N–D/N axial family or whether a different wall-profile family is forced by the full PDE.

Either way, the theorem target is now much sharper than it was at the end of Stage 025.
