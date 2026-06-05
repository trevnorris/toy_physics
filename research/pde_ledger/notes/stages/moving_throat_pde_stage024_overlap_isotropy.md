# Moving-Throat PDE — Stage 024: Explicit Overlap Integrals, the `O(3)` Isotropy Theorem, and the First Axisymmetric Splitting Law

## Purpose

Stage 023 reduced the entire grouped real `20/21/22` bundle to one exact question:

`mhat_0^2 N_0 / (K - B_0 - Z_0) = 54 G c_s^5 / (5 a^5 c^5)`.

But Stage 023 still treated the microscopic inputs

- `B_(A,n)`,
- `Z_(A,n)`,
- `N_(A,n)`

as formal coefficients.

The next honest step is therefore **not** another abstract ratio manipulation. It is to write the first actual overlap-integral model for the grouped real `P2` bundle and ask what the angular geometry already forces before the unsolved radial/axial PDE is touched.

This stage does exactly that.

The main outputs are:

1. an exact normalized real-STF harmonic basis on the throat sphere,
2. an exact angular source-map identity for the grouped `20/21/22` quadrupole ports,
3. an `O(3)` isotropy theorem showing that every isotropic reduced kernel collapses the three grouped lanes to one common value,
4. explicit common radial/axial overlap formulas for `B_n`, `Z_n`, and `N_n`,
5. the first exact symmetry-breaking law for a weak **axisymmetric quadrupole** perturbation,
6. and the resulting first-order transport law for the normalization ratio `P_A = N_A / D_A`.

So Stage 024 is the first point where the words “actual overlap integrals on the moving-throat branch” become mathematically concrete.

---

## 1. Normalized real STF harmonics on the throat sphere

Let `n^i` be the unit direction on the throat sphere and let

`E_(20), E_(21c), E_(21s), E_(22c), E_(22s)`

be the canonical real STF basis already used in the 2.5PN audit. Define the real angular functions

`Y_A(n) = sqrt(15/(8 pi)) E_A^{ij} n_i n_j`,

with `A in {20,21c,21s,22c,22s}`.

The exact fourth-moment identity on the unit sphere gives

`int dOmega n_i n_j n_k n_l = (4 pi / 15) (delta_ij delta_kl + delta_ik delta_jl + delta_il delta_jk)`.

Since the STF basis satisfies `Tr(E_A E_B) = delta_AB`, the normalized real harmonics obey

`int dOmega Y_A Y_B = delta_AB`.

So this basis is orthonormal without any further angular normalization choices.

That one fact already matters for the theorem gate, because it means the grouped real quadrupole bundle has a **canonical angular basis** on the isotropic throat.

---

## 2. Exact angular source-map identity

Write the orbital/worldtube STF quadrupole source on the throat sphere as

`S(n) = sum_A S_A Y_A(n)`.

Projecting onto the same port basis gives

`S_A^(port) = int dOmega Y_A(n) S(n) = sum_B S_B int dOmega Y_A Y_B = S_A`.

So on the natural isotropic source map the **angular** matching matrix is exactly the identity.

This means that the still-open source normalization `mhat_0` factors as

`mhat_0 = mhat_rad * mhat_ang`,

with

`mhat_ang = 1`

exactly on the canonical isotropic branch.

So the remaining normalization gap is not an angular mismatch between the orbital STF source and the grouped real throat ports. It is a **radial/axial and port-amplitude** issue.

That is a real narrowing of the theorem problem.

---

## 3. Explicit overlap-integral factorization on the isotropic branch

Take the grouped real wall deformation to be expanded as

`eta_A(s,Omega,t) = q_A(t) chi_eta(s) Y_A(Omega)`.

Take the stable BdG support sector in the same `l=2` bundle as

`X_(alpha,A)(s,Omega,t) = X_(alpha,A)(t) phi_alpha(s) Y_A(Omega)`.

Take the conservative brane-like gauge and mixed `A_w/F_(mu w)/J^w` sectors as

`U_(r,A)(s,Omega,t) = U_(r,A)(t) u_r(s) Y_A(Omega)`,

`W_(r,A)(s,Omega,t) = W_(r,A)(t) w_r(s) Y_A(Omega)`.

Assume the reference throat and the reduced kernels are `O(3)` invariant, so that the angular dependence enters only through scalar contractions.

Then all angular overlaps collapse by orthonormality, and every microscopic coupling becomes lane-independent:

`c_(A,alpha) = C_alpha = lambda_(B,alpha) I_(eta,alpha)`,

`g_(U,A,r) = G_(U,r) = lambda_(U,r) I_(eta,u,r)`,

`g_(W,A,r) = G_(W,r) = lambda_(W,r) I_(eta,w,r)`,

`R_(A,r) = R_r = lambda_(R,r) I_(u,w,r)`.

The radial/axial overlap integrals are

`I_(eta,alpha) = int ds mu_s(s) chi_eta(s) phi_alpha(s)`,

`I_(eta,u,r)  = int ds mu_s(s) chi_eta(s) u_r(s)`,

`I_(eta,w,r)  = int ds mu_s(s) chi_eta(s) w_r(s)`,

`I_(u,w,r)    = int ds mu_s(s) u_r(s) w_r(s)`.

These are exactly the Stage-023 radial/axial overlap objects specialized to an
`O(3)`-invariant kernel. This stage adds the angular closure and the symmetry
consequences; it does not introduce a different radial/axial model.

So on the isotropic branch the full Stage-023 coefficients become true scalar lane-independent objects:

`B_(A,0) = B_0 = sum_alpha C_alpha^2 / varpi_alpha^2`,

`B_(A,2) = B_2 = sum_alpha C_alpha^2 / varpi_alpha^4`,

`B_(A,4) = B_4 = sum_alpha C_alpha^2 / varpi_alpha^6`.

For each conservative Maxwell/mixed pair `r`, define

`Delta_r = Omega_(U,r)^2 Omega_(W,r)^2 - R_r^2`,

`S_r = Omega_(U,r)^2 + Omega_(W,r)^2`,

`Q_r = G_(U,r)^2 Omega_(W,r)^2 + 2 G_(U,r) G_(W,r) R_r + G_(W,r)^2 Omega_(U,r)^2`,

`H_r = G_(U,r)^2 + G_(W,r)^2`,

where `H_r` is just the Stage-023 combined gauge/mixed coupling strength written
with a new letter to avoid collisions with Newton's `G`.

`P_r = Omega_(U,r)^2 G_(W,r) + R_r G_(U,r)`.

Then

`Z_0 = sum_r Q_r / Delta_r`,

`Z_2 = sum_r [ Q_r S_r - H_r Delta_r ] / Delta_r^2`,

`Z_4 = sum_r [ Q_r (S_r^2 - Delta_r) - S_r H_r Delta_r ] / Delta_r^3`,

and

`N_0 = sum_r P_r^2 / Delta_r^2`,

`N_2 = sum_r 2 P_r (P_r S_r - Delta_r G_(W,r)) / Delta_r^3`,

`N_4 = sum_r [ Delta_r^2 G_(W,r)^2 - 2 Delta_r P_r^2 - 4 Delta_r P_r S_r G_(W,r) + 3 P_r^2 S_r^2 ] / Delta_r^4`.

So the conservative wall operator is exactly

`D_A(omega) = D_0 + D_2 omega^2 + D_4 omega^4 + O(omega^6)`

with

`D_0 = K - B_0 - Z_0`,

`D_2 = -(M + B_2 + Z_2)`,

`D_4 = -(B_4 + Z_4)`.

The grouped-lane consequence is immediate:

`D_(20,n) = D_(21,n) = D_(22,n) = D_n`,

`N_(20,n) = N_(21,n) = N_(22,n) = N_n`,

and therefore

`a_(D,n) = b_(D,n) = 0`,

`a_(N,n) = b_(N,n) = 0`.

So inside any truly `O(3)`-invariant reduced kernel, the grouped real `20/21/22` bundle is forced to be isotropic.

This is the first honest PDE-side isotropy theorem reached in the program, even though it is still a reduced-sector theorem rather than a full nonlinear moving-throat theorem.

---

## 4. The first allowed symmetry-breaking channel: a weak axisymmetric quadrupole background

Once isotropy is understood, the next question is the first way it can fail.

The leading symmetry-breaking correction that talks directly to the grouped real `P2` bundle is a weak quadrupolar anisotropy. In the axisymmetric frame, write the perturbation as

`delta K = eps kappa(s) Y_20(Omega)`.

The angular matrix that perturbs any `l=2` overlap is then

`M_AB^(20) = int dOmega Y_A Y_20 Y_B`.

Using the exact sixth moment of the unit sphere,

`int dOmega n_i n_j n_k n_l n_m n_n = (4 pi / 105) sum_pairings delta delta delta`,

one finds the exact five-mode result

`M^(20) = kappa_* diag(1, 1/2, 1/2, -1, -1)`,

with

`kappa_* = sqrt(5) / (7 sqrt(pi))`.

So the weak axisymmetric quadrupole background does **not** produce an arbitrary lane splitting. It produces one universal angular fingerprint:

- `20` lane shift proportional to `+1`,
- `21c,21s` shifts proportional to `+1/2`,
- `22c,22s` shifts proportional to `-1`.

After regrouping the `c/s` pairs into the three grouped lanes, the pattern is

`lambda_(20) = 1`,

`lambda_(21) = 1/2`,

`lambda_(22) = -1`.

Any first-order microscopic coefficient on that branch therefore has the form

`x_(20) = x^(0) + eps x^(1)`,

`x_(21) = x^(0) + (eps/2) x^(1)`,

`x_(22) = x^(0) - eps x^(1)`.

The weighted grouped trace/anomaly variables are then

`xbar = x^(0)`,

`a_x = (eps/4) x^(1)`,

`b_x = (3 eps / 4) x^(1)`.

So the first axisymmetric symmetry-breaking law is

`b_x = 3 a_x`.

This is a strong and very usable result. It means that if a future PDE computation shows a weak grouped-lane anisotropy produced by an axisymmetric `l=2` distortion of the isotropic branch, the defects are **not** free. They must sit on the one-dimensional line

`b = 3 a`.

If that relation fails, the symmetry breaking is not a pure axisymmetric quadrupole perturbation. It must involve either

- non-axisymmetric `m != 0` structure,
- higher-rank angular content,
- or a more complicated non-separable reduced kernel.

---

## 5. First-order normalization transport on the weak axisymmetric branch

Now apply the same axisymmetric law to the Stage-023 normalization ratio.

Suppose

`D_A = D_0 + eps lambda_A D_1 + O(eps^2)`,

`N_A = N_0 + eps lambda_A N_1 + O(eps^2)`,

with

`lambda_(20)=1`,

`lambda_(21)=1/2`,

`lambda_(22)=-1`.

Then the grouped-lane prefactor becomes

`P_A = N_A / D_A = P_0 + eps lambda_A P_1 + O(eps^2)`

with

`P_0 = N_0 / D_0`,

`P_1 = (N_1 D_0 - N_0 D_1) / D_0^2`.

So the first normalization anisotropies obey the same exact defect law:

`abar_P = 0`,

`a_P = (eps/4) P_1`,

`b_P = (3 eps / 4) P_1`,

`b_P = 3 a_P`.

This is the first actual transport law for the grouped normalization test once symmetry is weakly broken.

It says that the isotropic normalization target from Stage 023 is stable in the expected way: weak axisymmetric anisotropy does not create an arbitrary three-parameter deformation of the normalization stack. It creates one universal first-order lane pattern.

---

## 6. What Stage 024 changes in the theorem problem

Stage 024 narrows the remaining gap in three important ways.

### 6.1 The angular part of the source map is no longer open

On the natural isotropic grouped real basis,

`mhat_ang = 1`

exactly.

So the remaining source normalization issue is radial/axial and dynamical, not angular.

### 6.2 The isotropy theorem is now explicit on the reduced PDE side

If the reference throat and reduced kernels are truly `O(3)` invariant, then the grouped real `20/21/22` bundle is forced to collapse exactly:

`u_2^(20) = u_2^(21) = u_2^(22)`,

`u_4^(20) = u_4^(21) = u_4^(22)`,

`P_0^(20) = P_0^(21) = P_0^(22)`.

So on that branch the remaining theorem gap is not isotropy itself. It is the actual scalar value of the common radial/axial overlap data.

### 6.3 The first symmetry-breaking pattern is now diagnostic

A weak axisymmetric `l=2` deformation predicts one exact grouped signature:

`(20,21,22) ~ (1, 1/2, -1)`,

or equivalently

`b = 3 a`.

So future PDE data can be classified immediately:

- if weak anisotropy obeys `b = 3 a`, it is consistent with a pure axisymmetric quadrupole perturbation of the isotropic branch;
- if not, the symmetry breaking must be more complicated.

---

## 7. Best current summary after Stage 024

The moving-throat PDE problem is now narrower than it was at the end of Stage 023.

The remaining higher-order bridge is no longer:

> somehow compute all grouped-lane coefficients.

It is now:

1. compute the **common radial/axial isotropic overlap integrals** on the natural branch,
2. insert them into
   `B_n`, `Z_n`, `N_n`,
3. test the single isotropic ratio
   `mhat_rad^2 N_0 / (K - B_0 - Z_0)`,
4. and then, only after that, study symmetry breaking around that branch using the exact Stage-024 angular fingerprints.

That is a real tightening of the theorem target.
