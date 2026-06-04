# Moving-Throat PDE — Stage 032: Explicit Finite-Throat Mode Integrals, Kernel-Level Couplings, and Elimination of the Abstract Selected-Branch Source Map

## Purpose

After Stage 031 the normalization bottleneck had been reduced to the selected-branch quantity

`mhat_-^2 P_{0,-}`,

but the source-map factor `mhat_-` was still being carried as an abstract branch datum.
That was honest, but it also meant one important part of the selected-branch normalization problem had not yet been tied back to an explicit finite-throat mode model.

The next natural question is therefore:

> in the first explicit finite-throat axial source model, does the same D/N overlap structure that generates the mixed-sector loading also determine the selected source map, so that `mhat_-` stops being an independent parameter?

This stage answers yes.
In the first local isotropic kernel model built on the exact N/N wall basis and the exact D/N half-wave, every relevant wall–internal and source–wall coupling is determined by the same overlap vector

`v = (kappa_0, kappa_1)^T`.

As a result, on the natural D/N source branch the selected-mode source map is

`mhat_- = (v.e_-)/kappa_0`,

so

`mhat_-^2 = s_- / kappa_0^2`,

with `s_- = (v.e_-)^2`.

This removes the abstract source-map factor from the selected-branch quadrupole target.
The full selected-branch normalization product becomes

`mhat_-^2 P_{0,-} = beta_0 s_-^2 / (kappa_0^2 lambda_-)`,

which is now completely determined by the conservative selected-mode spectral data and the outgoing transfer factor.

There is also a useful side result:
under this natural axial source model the source-map amplification is modest and exact,

`1 <= mhat_-^2 < 11/9`.

So the normalization burden is not being carried by a large hidden source renormalization.
It is being carried mainly by the selected stiffness `lambda_-` and the mixed-sector transfer factor `beta_0`.

---

## 1. Exact finite-throat axial basis and the overlap vector `v`

Keep the same finite-throat axial interval `s in [0,L]` and the same wall basis used from Stage 027 onward.
The first two exact N/N wall modes are

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`.

The natural compact outgoing/internal half-wave is the exact D/N mode

`f_0(s) = sqrt(2/L) sin(pi s / (2L))`.

The wall-to-D/N overlap vector is therefore

`v_i = int_0^L u_i(s) f_0(s) ds`.

The exact values are

`kappa_0 = 2 sqrt(2) / pi`,

`kappa_1 = - 4 / (3 pi)`,

`v = (kappa_0, kappa_1)^T`.

So the squared norm is

`sigma = v.v = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2)`.

It is also useful to record the flat-branch normalization constant

`kappa_0^2 = 8 / pi^2`,

and therefore the exact ratio

`sigma / kappa_0^2 = 11 / 9`.

That number will become the maximal source-map enhancement on the selected branch.

---

## 2. Local isotropic kernel model and exact reduced couplings

Now make the first explicit finite-throat bilinear-kernel choice that is still compatible with the project ontology:

- a wall field `eta(s)` expanded on `{u_0,u_1}`,
- a brane-like internal doublet `U(s)` expanded on the same N/N basis,
- one BdG support half-wave `phi(s) = phi f_0(s)`,
- one mixed `A_w/F_(mu w)/J^w` half-wave `W(s) = W f_0(s)`.

Take the simplest local isotropic couplings,

`L_(eta U) = g_U int_0^L eta(s) U(s) ds`,

`L_(eta phi) = g_B int_0^L eta(s) phi(s) ds`,

`L_(eta W) = g_W int_0^L eta(s) W(s) ds`,

`L_(U W) = - g_R int_0^L U(s) W(s) ds`.

Expanding in modes gives the exact reduced couplings

`L_(eta U) = g_U q.u`,

`L_(eta phi) = g_B (v.q) phi`,

`L_(eta W) = g_W (v.q) W`,

`L_(U W) = - g_R (v.u) W`.

So in the first explicit finite-throat kernel model,

`lambda_U = g_U`,

`lambda_B = g_B`,

`lambda_W = g_W`,

`lambda_R = g_R`,

and the entire directional structure is carried by the same overlap vector `v`.

This is exactly the pattern that had been assumed abstractly in Stage 029; now it is derived from the explicit finite-throat basis and local isotropic kernels.

---

## 3. Exact Schur-complement reduction with explicit kernel-level couplings

Let the reduced frequency-domain kernels be

`A_phi(omega) = varpi^2 - omega^2`,

`A_U(omega)   = Omega_U^2 - omega^2`,

`A_W(omega)   = Omega_W^2 - omega^2 - Pi_out(omega)`.

With the local couplings above, elimination of the internal block `(u, phi, W)` gives the same exact wall self-energy structure found in Stage 029,

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T`,

but now with explicit kernel-level meaning:

`Xi(omega) = g_U^2 / A_U(omega)`,

`alpha(omega) = g_B^2 / A_phi(omega)`
`             + ( A_U(omega) g_W + g_R g_U )^2`
`               / [ A_U(omega) Delta_UW(omega) ]`,

`Delta_UW(omega) = A_U(omega) A_W(omega) - g_R^2 sigma`.

So the static selected-branch data are no longer abstract at this level.
They are

`Xi_0 = g_U^2 / Omega_U^2`,

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`,

`alpha_0 = g_B^2 / varpi^2 + ( Omega_U^2 g_W + g_R g_U )^2 / ( Omega_U^2 Delta_0 )`,

`beta_0 = ( Omega_U^2 g_W + g_R g_U )^2 / Delta_0^2`.

This is the first point where the Stage-029/031 selected-branch spectral data have been written directly in terms of explicit finite-throat mode integrals and local bilinear kernel strengths.

---

## 4. Natural D/N source branch and the exact selected source map

Now attach the orbital/worldtube STF quadrupole source through the same D/N mouth load,

`L_src = g_Q Q_STF int_0^L eta(s) f_0(s) ds`.

Expanding `eta(s) = q_0 u_0(s) + q_1 u_1(s)` gives

`L_src = g_Q Q_STF (v.q)`.

So the external source vector in wall-coordinate space is simply

`J_src = g_Q Q_STF v`.

Projecting onto the selected lower eigenvector `e_-` therefore gives the selected source amplitude

`J_- = g_Q Q_STF (v.e_-) = g_Q Q_STF sqrt(s_-)`,

where

`s_- = (v.e_-)^2`.

At zero loading the lower mode is just the flat branch, so

`J_-(0) = g_Q Q_STF kappa_0`.

This defines the natural selected-branch source map,

`mhat_- = J_-(alpha_0) / J_-(0) = (v.e_-) / kappa_0`,

and hence the exact squared source-map factor is

`mhat_-^2 = s_- / kappa_0^2`.

So on the natural D/N source branch the abstract factor `mhat_-` is gone.
It is completely fixed by the same selected overlap `s_-` that already controls the wall loading.

---

## 5. Exact bound on the source-map factor

Because the selected overlap grows monotonically from the flat branch to the max-coupling branch,

`kappa_0^2 <= s_- < sigma`,

and therefore

`1 <= mhat_-^2 < sigma / kappa_0^2 = 11/9`.

So the natural source-map factor stays in the exact window

`1 <= mhat_-^2 < 11/9`,

or

`1 <= mhat_- < sqrt(11/9)`.

This is a useful structural simplification.
The selected-branch normalization problem cannot be hidden inside a huge undetermined source-map amplification on the natural D/N source branch.
The source factor is real, positive, monotone, and modest.

---

## 6. Elimination of the abstract source-map factor from the quadrupole target

Stage 030 wrote the selected-branch normalization quantity as

`mhat_-^2 P_{0,-}`

with

`P_{0,-} = beta_0 s_- / lambda_-`.

Substituting the source-map result from Section 4 gives the exact product

`mhat_-^2 P_{0,-} = beta_0 s_-^2 / ( kappa_0^2 lambda_- )`.

So the invariant 2.5PN target becomes

`beta_0 s_-^2 / ( kappa_0^2 lambda_- ) = 54 G c_s^5 / (5 a^5 c^5)`.

This is sharper than the Stage-030 formulation because there is no longer an independent `mhat_-` datum on the natural D/N source branch.
Everything is now carried by:

- the explicit mixed-sector transfer factor `beta_0`,
- the selected overlap `s_-`,
- the selected conservative eigenvalue `lambda_-`,
- and the known flat-branch overlap `kappa_0`.

---

## 7. Best current summary after Stage 032

The selected-branch theorem gap has narrowed again.

The first explicit finite-throat mode-integral model now does three things at once:

1. it derives the Stage-029 loading structure from explicit local isotropic kernels,
2. it writes `Xi_0`, `alpha_0`, and `beta_0` directly in terms of those kernel strengths,
3. and it removes the abstract selected-branch source-map factor by showing that on the natural D/N source branch
   `mhat_-^2 = s_- / kappa_0^2`.

So the remaining normalization problem is no longer

- “selected stiffness plus an unknown source map,”

but rather

- “selected stiffness plus selected overlap plus the explicit mixed-sector transfer factor.”

The next honest step is therefore to write the full selected-branch normalization equation directly in microscopic coupling language and see what exact stability and reachability constraints it imposes.
