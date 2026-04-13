# 5PN progress notes — explicit overlap model and weak-axisymmetric transport

## What was added

### Stage 3 — explicit isotropic finite-throat overlap model

Files:
- `5pn_stage3_isotropic_overlap_model.py`
- `5pn_stage3_isotropic_overlap_model_output.txt`

Model choice:
- finite interval `s in [0,L]`
- wall profile
  `chi_eta(s) = sqrt(2/L) sin(pi s / L)`
- D/N-like half-wave
  `phi_DN(s) = sqrt(2/L) sin(pi s / (2L))`
- simplest brane-like gauge profile `u = chi_eta`
- simplest mixed profile `w = phi_DN`

Exact overlaps:
- `I_eta_u = 1`
- `I_eta_phi = I_eta_w = I_uw = 8/(3 pi)`

Using one BdG mode and one conservative Maxwell/mixed pair, the script constructs
- `B0,B2,B4`
- `Z0,Z2,Z4`
- `N0,N2,N4`
- `D0,D2,D4`
- `u2,u4`
- `P0,P2,P4`

and verifies exact grouped isotropy when the three grouped lanes share the same radial/axial data.

### Most useful Stage-3 result

In this explicit prototype, two important requirements both fix the same static wall stiffness `K`:

1. the 5PN conservative one-pole identity
   `u4 = 4 u2^2`
2. the universal 2.5PN / 4PN normalization target for `P0`

So simultaneous success reduces to one explicit compatibility equation:

`N0 / P0_target = 3 (M + B2 + Z2)^2 / (B4 + Z4)`

This is the first concrete radial/axial compatibility surface I have for the 5PN program.
It is the explicit prototype version of “the PDE must make the conservative and outgoing calibrations agree on one branch.”

## Stage 4 — weak-axisymmetric transport and microscopic grouped obstructions

Files:
- `5pn_stage4_axisymmetric_transport.py`
- `5pn_stage4_axisymmetric_transport_output.txt`

This stage starts from the exact weak-axisymmetric grouped signature
- `lambda_(20) = 1`
- `lambda_(21) = 1/2`
- `lambda_(22) = -1`

and carries the grouped operator / transfer slopes through the full Stage-154/155/156 logic.

### Exact microscopic slope decomposition

Using
- wall slopes `(K1_wall, M1_wall)`
- BdG slopes `(B01, B21, B41)`
- conservative Maxwell/mixed slopes `(Z01, Z21, Z41)`
- outgoing transfer slope `N01`

it verifies

`K_1 = D21 + D01/9 = W1 - Bcal1 - Zcal1`

and

`G_1 = N01 - P0 D01 = -P0 K1_wall + P0 B01 + Nbundle1`

with
- `D01 = K1_wall - B01 - Z01`
- `D21 = -(M1_wall + B21 + Z21)`
- `D41 = -(B41 + Z41)`

### Physical weak-axisymmetric slopes

The script derives

`u2^(1) = -(D21 + u2 D01)/D0`

`P1/P0 = N01/N0 - D01/D0`

and on the canonical compensated branch verifies the hidden-even relation

`u4^(1) = (8/9) u2^(1)`

iff

`D41 = (2/3) D21 + D01/27`

### Most useful Stage-4 result

On the even-preserving branch `u2^(1)=0`, the conservative grouped response collapses to

- `D21 = -D01/9`
- `D41 = -D01/27`

and the remaining linear grouped `2.5`PN defect becomes one scalar loading mismatch

`Xi_load = N01/N0 - D01/D0 = P1/P0`

with fixed grouped-lane signature

- `(20,21,22) ~ (1, 1/2, -1)`
- equivalently `b = 3 a`

So after Stage 4, the next theorem gate is no longer “solve all grouped anisotropies.”
It is much narrower:

> compute `D01` and `N01` — or directly `Xi_load` — from a primitive weak-axisymmetric deformation of the explicit finite-throat overlap model.

## Run status

Both new scripts completed here without running out of time or memory.

## Stage 5 — primitive deformation and exact compensation surfaces

Files:
- `5pn_stage5_primitive_deformation_compensation.py`
- `5pn_stage5_primitive_deformation_compensation_output.txt`

This stage takes the explicit isotropic overlap model from Stage 3 and adds a primitive weak-axisymmetric microscopic deformation through the slopes

- `dK`, `dM`
- `d(lambda_B)`, `d(varpi)`
- `d(lambda_U)`, `d(lambda_W)`, `d(lambda_R)`
- `d(Omega_U)`, `d(Omega_W)`

Then it computes the induced grouped-lane slope data

- `D01`, `D21`, `D41`
- `N01`

and from them the three key combinations

- `K1 = D21 + D01/9`
- `G1 = N01 - P0 D01`
- `Xi_load = N01/N0 - D01/D0 = G1/N0`

### Most useful Stage-5 results

1. **Even-preserving compensation** `K1 = 0` is algebraic and fixes the inertia-side slope `dM` exactly.
2. **Odd/normalization-preserving compensation** `Xi_load = 0` is algebraic and fixes the static loading slope `dK` exactly.
3. The remaining explicit **5PN hidden-even residual** is
   `D41 - (2/3) D21 - D01/27`.

So after Stage 5, the next theorem gate is now extremely sharp:

> choose a concrete primitive anisotropy mechanism — wall-only, BdG-only, Maxwell/mixed-only, or a mixed combination — substitute its slopes into Stage 5, and test
> `K1 = 0`, `Xi_load = 0`, and `D41 - (2/3) D21 - D01/27 = 0`.

## Run status

All three new scripts completed here without running out of time or memory.
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

# 5PN continuation notes — stages 8 through 11

These stages continue the weak-axisymmetric grouped-`P2` program from the Stage-7 scalar amplitude
\[
\Xi_1=\frac{P_1}{P_0}.
\]

The overall effect is that the remaining first-order grouped normalization problem is no longer a generic “compute lots of microscopic drifts” problem. It has collapsed to a small sequence of exact equivalent formulations.

## Stage 8 — direct outgoing-port co-loading

The remaining grouped weak-axisymmetric defect can be written directly in terms of the actual outgoing-port slopes:
\[
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1)
=
\bar\nu_N-\kappa_1,
\]
where
\[
\nu_r = 2(\mathfrak p_r-\mathfrak d_r)
\]
is the logarithmic slope of
\[
N_{0}^{(r)}=\frac{P_r^2}{\Delta_r^2}.
\]

Equivalently, if
\[
N_{A,0}^{(r)} = K_A \mathcal T_{A,r}^2,
\qquad
\delta\ln \mathcal T_{A,r}=\epsilon\lambda_A\,\tau_r,
\]
then
\[
\nu_r=\kappa_1+2\tau_r,
\qquad
\Xi_1 = 2\sum_r \rho_r^{(N)}\tau_r.
\]

So the exact zero-defect condition is now:
\[
\sum_r \rho_r^{(N)}\tau_r=0.
\]

A stronger sufficient condition is per-port co-loading:
\[
\tau_r=0
\qquad
\text{for every active outgoing port }r.
\]

Under upper-leg and coupling rigidity, the transfer-shape slope collapses to the raw mixed-leg slope and recovers the old square-root mixed-leg law.

## Stage 9 — one effective transfer shape and the actual continuum branch

The many-port weighted sum collapses exactly to one effective transfer shape:
\[
\mathcal T_{\mathrm{eff},A}^2
=
\sum_r \mathcal T_{A,r}^2
=
\frac{N_{A,0}}{K_A},
\qquad
\Xi_1
=
\frac{\delta\ln\mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}.
\]

On the actual one-port continuum branch,
\[
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2},
\]
so
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\rho_1}{1+\rho}+\frac{2\varepsilon_W}{1-\epsilon_W}.
\]

In selected-branch form,
\[
\mathcal T_A^2
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\,
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}},
\]
which yields
\[
\Xi_1
=
-\frac{\eta_1}{1-\epsilon_\eta}
-\mathcal R_1.
\]

On the coherent local D/N branch,
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\epsilon
=
\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
and the defect becomes
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}
+\frac{2\epsilon_1}{1-\epsilon}.
\]

The support ratio drops out identically:
\[
\partial_\zeta \ln \mathcal T^2 = 0.
\]

So the coherent defect is support-blind, and exact tracking rigidity by itself is not enough to kill it.

## Stage 10 — microscopic coherent-kernel slippages and exact triangular normal form

The coherent branch depends only on the microscopic slippages
\[
\Sigma_Z,\quad
\Sigma_\chi,\quad
\Sigma_\epsilon,\quad
\Sigma_\delta,
\]
with one additional dressing slippage
\[
\Sigma_\eta
\]
entering the selected-branch form.

The exact microscopic grouped-defect law is
\[
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta
\right].
\]

The exact tracking combination is
\[
\Sigma_{\rm tr}
=
(1+\chi_0)\Sigma_\delta
+
(1+\delta_U)\Sigma_\chi,
\]
with
\[
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
\]

Defining the genuine nontracking slippage
\[
\Sigma_{\rm nt},
\]
the coherent problem takes the exact triangular normal form
\[
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
\]

So on the constructive coherent branch,
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.
\]

This is the first exact three-scalar normal form of the full coherent weak-axisymmetric problem.

## Stage 11 — direct microscopic monomials, similarity orbit, and quotient closure

The three branch-adapted coordinates can be written as logarithmic drifts of three exact microscopic monomials:
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=\Sigma_\eta,
\]
with
\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\left(\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}\right)
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
\]

The exact monomial-drift map is the rank-3 matrix
\[
M_*\,\delta\mathbf x
=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix},
\]
with
\[
\det M_*^{(\tau,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0,
\]
so
\[
\dim\ker M_* = 5.
\]

There is an exact five-parameter multiplicative similarity orbit \(\mathcal G_*\) preserving the three monomials exactly, and the scripts show
\[
\ker M_* = T_{\mathrm{id}}\mathcal G_*.
\]

More strongly, the exact finite invariant-fibre equations are
\[
M_*\,\Delta \mathbf x = 0,
\]
and solving them reproduces the exact orbit exponents. Therefore the finite level sets of
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\]
are precisely the similarity orbits \(\mathcal G_*\).

So the coherent weak-axisymmetric zero-defect theorem can now be written as
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x \in T_{\mathrm{id}}\mathcal G_*,
\]
and, at finite level,
\[
\mathcal M_+/\mathcal G_*
\cong
(\mathbb R_{>0})^3
\]
with quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
\]

## What this means for the 5PN program

At this point the weak-axisymmetric normalization problem is no longer an algebraic bottleneck.

There are now three equivalent ways to state the next theorem gate:

1. compute the direct grouped scalar
   \[
   \Xi_1=\frac{P_1}{P_0},
   \]
2. compute the branch-adapted defect coordinates
   \[
   (\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta),
   \]
3. or test whether the actual moving-throat weak-axisymmetric branch is tangent to the exact monomial-preserving similarity orbit \(\mathcal G_*\).

That last form is the sharpest current continuation point.
# 5PN continuation notes — Stage 12–13 normalized monomial bridge

These two stages connect the earlier explicit Stage-5 primitive prototype

a) directly to the Stage-10/11 microscopic similarity-orbit package, and
b) into formulas that can be used on the actual moving-throat branch.

## Stage 12 — normalized prototype already contains the Stage-11 quotient coordinates

Using the coherent-kernel dictionary

- `K = K_eta^(eff)`
- `M = mu_eta`
- `G_U = c_etaU / sqrt(mu_eta mu_U)`
- `G_W = lambda_W / sqrt(mu_eta mu_W)`
- `R = gamma lambda_W / sqrt(mu_U mu_W)`
- `Omega_U^2 = K_U / mu_U`
- `Omega_W^2 = K_W^(eff) / mu_W`
- `delta_U = pi^2 T_U / (L^2 K_U)`

one gets the exact normalized coherent ratios

\[
\chi_0 = \frac{R G_U}{\Omega_U^2 G_W},
\qquad
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2},
\qquad
\epsilon_W = \frac{R^2\sigma}{\Omega_U^2\Omega_W^2},
\qquad
Z_W = \frac{M G_W^2}{K\Omega_W^2}.
\]

The direct Stage-168/169 monomials become

\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{R G_U}{\Omega_U^2 G_W}\right)^{1+\delta_{U,*}}
\delta_U^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\frac{M G_W^2}{K\Omega_W^4}
\left(\frac{R^2\sigma}{\Omega_U^2\Omega_W^2}\right)^{E_*}
\delta_U^{-F_*},
\]

\[
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2}.
\]

So the Stage-5-style normalized prototype already contains the exact Stage-11 quotient coordinates once the split-`U` variable `delta_U` is admitted.

The Stage-10 slippages collapse to the normalized drift formulas

\[
\Sigma_Z = d\ln M + 2 d\ln G_W - d\ln K - 4 d\ln\Omega_W,
\]
\[
\Sigma_\chi = d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U,
\]
\[
\Sigma_\eta = d\ln M + 2 d\ln G_U - d\ln K - 2 d\ln\Omega_U,
\]
\[
\Sigma_\epsilon = 2(d\ln R - d\ln\Omega_U - d\ln\Omega_W),
\qquad
\Sigma_\delta = d\ln\delta_U.
\]

So the extra raw mass/stiffness bookkeeping mostly cancels out of the actual defect coordinates.

## Stage 13 — exact zero-defect kernel in normalized Stage-5 variables

In the normalized variables

\[
(d\ln G_W,
 d\ln G_U,
 d\ln R,
 d\ln K,
 d\ln M,
 d\ln\Omega_U,
 d\ln\Omega_W,
 d\ln\delta_U),
\]

the direct monomial-drift matrix is

\[
M_{\rm norm}=
\begin{pmatrix}
-(1+\delta_U) & 1+\delta_U & 1+\delta_U & 0 & 0 & -2(1+\delta_U) & 0 & 1+\chi_0 \\
2 & 0 & 2E_* & -1 & 1 & -2E_* & -(4+2E_*) & -F_* \\
0 & 2 & 0 & -1 & 1 & -2 & 0 & 0
\end{pmatrix}.
\]

Its rank is exactly `3`, so the normalized zero-defect tangent space has dimension `5`.

The exact zero-defect equations are

\[
(1+\delta_U)(d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U)
+
(1+\chi_0)d\ln\delta_U = 0,
\]

\[
d\ln M + 2 d\ln G_U - d\ln K - 2 d\ln\Omega_U = 0,
\]

\[
d\ln M + 2 d\ln G_W - d\ln K - 4 d\ln\Omega_W
+
2E_*(d\ln R - d\ln\Omega_U - d\ln\Omega_W)
-
F_* d\ln\delta_U = 0.
\]

These solve triangularly:

### Tracking
\[
d\ln\delta_U
=
-\frac{1+\delta_U}{1+\chi_0}
\bigl(d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U\bigr).
\]

### Dressing
\[
d\ln M
=
 d\ln K - 2 d\ln G_U + 2 d\ln\Omega_U.
\]

### Nontracking
\[
d\ln\Omega_W
=
\frac{d\ln G_W - d\ln G_U + (1-E_*)d\ln\Omega_U + E_* d\ln R - \tfrac{F_*}{2}d\ln\delta_U}{E_*+2}.
\]

So once a candidate moving-throat branch gives the drifts of

- `K`
- `G_U`
- `G_W`
- `R`
- `Omega_U`

it automatically fixes the drifts required in

- `delta_U`
- `M`
- `Omega_W`

to stay tangent to the exact similarity orbit.

## Practical Stage-5 absolute-slope form

Writing the Stage-5 primitive slopes as

- `dK`
- `dM`
- `d(lambda_U)`
- `d(lambda_W)`
- `d(lambda_R)`
- `d(Omega_U)`
- `d(Omega_W)`
- `d(delta_U)`

one gets the exact compatibility formulas

\[
d(\delta_U)
=
-\delta_U\frac{1+\delta_U}{1+\chi_0}
\left[
\frac{d\lambda_R}{\lambda_R}
+
\frac{d\lambda_U}{\lambda_U}
-
\frac{d\lambda_W}{\lambda_W}
-
2\frac{d\Omega_U}{\Omega_U}
\right],
\]

\[
dM
=
M\left[
\frac{dK}{K} - 2\frac{d\lambda_U}{\lambda_U} + 2\frac{d\Omega_U}{\Omega_U}
\right],
\]

\[
\frac{d\Omega_W}{\Omega_W}
=
\frac{1}{E_*+2}
\left[
\frac{d\lambda_W}{\lambda_W}
-
\frac{d\lambda_U}{\lambda_U}
+
(1-E_*)\frac{d\Omega_U}{\Omega_U}
+
E_*\frac{d\lambda_R}{\lambda_R}
-
\frac{F_*}{2}\frac{d(\delta_U)}{\delta_U}
\right].
\]

So the Stage-10/11 similarity-orbit theorem is now directly usable in the Stage-5 primitive deformation language.

## Immediate consequence for the 5PN program

The next honest theorem gate is now smaller than before:

1. extract the actual branch drifts of
   `K, lambda_U, lambda_W, lambda_R, Omega_U`
   from the moving-throat PDE,
2. use the formulas above to predict the required co-drifts of
   `delta_U, M, Omega_W`
   if the branch is tangent to the exact zero-defect similarity orbit,
3. then compare those predictions with the actual reduced PDE branch.

If they agree, the weak-axisymmetric first-order obstruction is pure similarity-gauge and the calibration survives. If they fail, the moving-throat branch leaves the exact monomial-preserving orbit and the calibration breaks for a concrete microscopic reason.

## Stage 14 — the BdG primitive drifts are exactly support-blind for the Stage-11 monomials

If the primitive drift space is extended back to the full Stage-5 list

- `lambda_B`
- `varpi`
- `lambda_U`
- `lambda_W`
- `lambda_R`
- `K`
- `M`
- `Omega_U`
- `Omega_W`
- `delta_U`

then the direct Stage-11 monomial-drift matrix acquires **two exact zero columns** in the `dln lambda_B` and `dln varpi` directions.

So the weak-axisymmetric direct monomials are exactly support-blind:

\[
\partial_{\ln \lambda_B}
(\delta\ln \mathfrak C_{{\rm tr},*},
 \delta\ln \mathfrak C_{{\rm nt},*},
 \delta\ln \epsilon_\eta)
=0,
\]
\[
\partial_{\ln \varpi}
(\delta\ln \mathfrak C_{{\rm tr},*},
 \delta\ln \mathfrak C_{{\rm nt},*},
 \delta\ln \epsilon_\eta)
=0.
\]

That means the extended primitive zero-defect tangent space has dimension

\[
2 + 5 = 7,
\]

namely

1. two BdG-support-blind directions,
2. plus the five normalized similarity directions from Stage 13.

This is important conceptually.

The Stage-10/11 similarity-orbit theorem constrains only the mixed/wall/U normalization problem. It does **not** constrain the explicit BdG support drifts. So monomial-rigidity alone cannot finish the full conservative 5PN problem.

Those BdG directions must still be controlled, if at all, by the separate conservative front-end conditions:

- `K1 = 0`,
- the hidden-even consistency slot,
- and the `O(omega^4)` single-pole / grouped-response test.

So the Stage-5/6 conservative extraction theorem and the Stage-10/11 similarity-orbit theorem are complementary rather than redundant.
# 5PN continuation notes — stages 15 through 17

These stages finally splice the Stage-5 primitive obstruction triplet into the later Stage-10/11 monomial-orbit package.

The main point is that the three objects are **not** on the same footing:

- `Xi_load = N01/N0 - D01/D0 = P1/P0` is exactly the weak-axisymmetric normalization defect governed by the monomial/similarity theorem;
- `K1 = D21 + D01/9` and
- `H_even = D41 - (2/3) D21 - D01/27`

are separate conservative even gates that survive after the normalization defect is killed.

## Stage 15 — the exact `Xi_load` bridge

The exact Stage-5 prefactor identity is
\[
P_0 = \frac{N_0}{D_0},
\qquad
P_1 = \frac{N_{01}D_0 - N_0 D_{01}}{D_0^2},
\qquad
\Xi_{\rm load}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}
=
\frac{P_1}{P_0}.
\]
So the Stage-5 load defect is **exactly** the same weak-axisymmetric prefactor slope already isolated later as
\[
\Xi_1 = \frac{P_1}{P_0}.
\]

Rewriting Stages 10–14 in one place, the direct monomial drifts are
\[
\Sigma_{\rm tr}=\delta\ln \mathfrak C_{{\rm tr},*},
\qquad
\Sigma_{\rm nt}=\delta\ln \mathfrak C_{{\rm nt},*},
\qquad
\Sigma_\eta=\delta\ln \epsilon_\eta,
\]
and the weak-axisymmetric normalization defect has the exact triangular form
\[
\Xi = A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt}.
\]
The Stage-13 normalized kernel annihilates all three monomial drifts, and the Stage-14 extension shows the explicit BdG directions `dln lambda_B` and `dln varpi` are zero columns of the monomial matrix. Therefore
\[
\Xi_{\rm load}=0
\]
on:

- the full normalized similarity kernel,
- its injected support-blind extension,
- the explicit BdG amplitude/frequency directions,
- the matched wall-only co-scaling direction.

So the monomial/similarity theorem really is a theorem about `Xi_load`.

## Stage 16 — why that still does **not** solve 5PN

The next check is to compare `Xi` with the conservative even gates.

On the matched wall-only co-scaling direction,
\[
d\ln K = d\ln M = 1,
\]
we still get
\[
K1_{\rm wall}=\frac{K}{9}-M,
\qquad
H_{{\rm even,wall}}=-\frac{K}{27}+\frac{2M}{3},
\]
which are generically nonzero even though `Xi = 0`.

On the explicit support-blind BdG directions,
\[
d\ln \lambda_B = 1,
\qquad d\ln \varpi = 0,
\]
or
\[
d\ln \lambda_B = 0,
\qquad d\ln \varpi = 1,
\]
we again have `Xi = 0`, but both `K1_B` and `H_even,B` are nonzero.

Even on the pure BdG self-similar branch
\[
d\ln \lambda_B = d\ln \varpi = \delta,
\]
which is the branch that kills the BdG load defect, one still gets
\[
K1_B = \frac{2B_0\delta}{\varpi^2},
\qquad
H_{{\rm even},B} = \frac{4B_0\delta(3-\varpi^2)}{3\varpi^4},
\]
so the even gates remain open unless `\delta = 0` (or an extra tuning is imposed).

That is the cleanest executable statement yet that:
\[
\text{similarity-orbit rigidity} \neq \text{full 5PN closure}.
\]
The orbit theorem kills `Xi`; it does **not** automatically close `K1` or `H_even`.

## Stage 17 — the first lower-bound intersection calculation

Once the Stage-13/14 monomial-rigid kernel is parameterized by
\[
(\alpha_K,\alpha_W,\alpha_U,\alpha_R,\alpha_{\Omega_U},\beta_B,\beta_{\varpi}),
\]
we can impose the **lower-bound** conservative even gates obtained by keeping only the explicit wall-only and pure-BdG pieces of `K1` and `H_even`.

This gives an exact `2 x 7` gate matrix of generic rank `2`, so the lower-bound even-gate intersection has dimension `5`.

A convenient exact solve is to determine
\[
\alpha_K,
\qquad
\beta_B,
\]
in terms of the remaining five free directions.

The corresponding null basis has five directions:

1. pure `\alpha_W`,
2. pure `\alpha_R`,
3. matched `\alpha_U/\alpha_{\Omega_U}`,
4. BdG-amplitude deformation `\beta_B` with compensating `\alpha_K,\alpha_U`,
5. BdG-frequency deformation `\beta_{\varpi}` with compensating `\alpha_K,\alpha_U`.

The important caution is that this is only a **lower-bound** solve, because the conservative mixed-sector `Z_2,Z_4` moments have not yet been reinstated. In particular, the fact that `\alpha_W` and `\alpha_R` survive untouched here is telling us exactly where the omitted mixed-sector moments still have to act.

## Net result after Stage 17

The 5PN continuation point is now much sharper:

1. `Xi_load = P1/P0` is fully absorbed into the similarity-orbit / monomial-rigidity theorem;
2. the real conservative 5PN bottleneck is the pair of even gates `K1` and `H_even`;
3. the explicit wall-only and BdG-only pieces show those even gates survive after `Xi` is killed;
4. and the next honest task is therefore:

> reinstate the conservative mixed-sector `Z_2,Z_4` moments and determine how they cut the remaining lower-bound tangent family.

That is the cleanest next theorem gate we have reached so far.
# 5PN continuation notes — stages 18 through 20

These stages do two things in parallel.

First, they sharpen the **isotropic full-bundle target surface** that the actual
moving-throat PDE must hit if the calibrated 5PN / 2.5PN / 4PN branch is real.
Second, they reinstate the omitted conservative Maxwell/mixed `Z_2,Z_4` sector in
the weak-axisymmetric even-gate problem and show exactly how it removes the fake
freedom left in the Stage-17 lower-bound picture.

## Stage 18 — exact isotropic full-bundle target surface

Work on the isotropic grouped-lane branch with

a) conservative operator moments

a) `D0 = K - B0 - Z0`,

b) `D2 = -(M + B2 + Z2)`,

c) `D4 = -(B4 + Z4)`.

Then the normalized grouped-response and outgoing-prefactor moments are

`u2 = -D2 / D0`,

`u4 = (D2^2 - D0 D4) / D0^2`,

`P0 = N0 / D0`,

`P2 = (D0 N2 - 2 D2 N0) / D0^2`,

`P4 = (D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0) / D0^3`.

The exact isotropic 5PN conservative one-pole defect is

`u4 - 4 u2^2 = [ D0 (B4 + Z4) - 3 (M + B2 + Z2)^2 ] / D0^2`.

So the isotropic one-pole surface is exactly

`D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`.

The constant-prefactor outgoing branch is also exact:

`P2 = 0  ->  N2 = 2 D2 N0 / D0`,

`P4 = 0  ->  N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

Finally, the universal 2.5PN / 4PN normalization target is

`mhat0^2 P0 = 54 G c_s^5 / (5 a^5 c^5)`

or equivalently

`mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`.

So the full isotropic moving-throat bundle must land on one exact combined target
surface:

1. `D0 = K - B0 - Z0`,
2. `D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`,
3. `mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`,
4. `N2 = 2 D2 N0 / D0`,
5. `N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

That is the sharpest isotropic full-bundle statement we have so far.

## Stage 19 — exact `Z`-sector bridge back into the even gates

The missing conservative Maxwell/mixed moments are

`Z0 = Q / Delta`,

`Z2 = (Q S2 - H Delta) / Delta^2`,

`Z4 = [ Q (S2^2 - Delta) - S2 H Delta ] / Delta^3`.

Their exact first variations are

`dZ0 = (Delta dQ - Q dDelta) / Delta^2`,

`dZ2 = [ Delta (-Delta dH - H dDelta + Q dS2 + S2 dQ) + 2 dDelta (Delta H - Q S2) ] / Delta^3`,

`dZ4 = -[ Delta^2 H dS2 + Delta^2 S2 dH + Delta^2 dQ - 2 Delta H S2 dDelta`
`         - 2 Delta Q S2 dS2 - 2 Delta Q dDelta - Delta S2^2 dQ + 3 Q S2^2 dDelta ] / Delta^4`.

Therefore the conservative `Z`-sector contributions to the two live 5PN even gates are

`K1_Z = -dZ2 - dZ0/9`,

`H_even,Z = -dZ4 + (2/3) dZ2 + dZ0/27`.

After inserting the exact Stage-13 normalized similarity kernel, both become honest
functions of the mixed-sector similarity directions `alpha_W` and `alpha_R`.

On the positive constructive slice

- `G_U = 5`, `G_W = 7`, `R = 2`,
- `Omega_U = 3`, `Omega_W = 4`,
- `chi_0 = 3/2`, `delta_U = 2/3`,
- `E_* = 1/4`, `F_* = 5/6`,

we get

`K1_Z = (78623501/25004700) alpha_OmegaU + (733046/6251175) alpha_R`
`       - (59010631/25004700) alpha_U - (32134513/50009400) alpha_W`,

`H_even,Z = -(28906377971/21003948000) alpha_OmegaU`
`           - (1174937411/21003948000) alpha_R`
`           + (11102468471/10501974000) alpha_U`
`           + (5617869293/21003948000) alpha_W`.

In particular, the pure directions give

`alpha_W:  K1_Z = -32134513/50009400,   H_even,Z = 5617869293/21003948000`,

`alpha_R:  K1_Z = 733046/6251175,       H_even,Z = -1174937411/21003948000`.

So the omitted `Z_2,Z_4` block does exactly what Stage 17 said it still had to do:
it activates the previously untouched mixed directions.

## Stage 20 — full even-gate solve on the constructive branch

Once wall-only, pure-BdG, and the reinstated `Z`-sector are combined, the full
constructive-slice even-gate matrix is

`Gate_full(slice) =`

`[ -25/9,  -32134513/50009400,   91017569/25004700,    733046/6251175,`
`  -71404699/25004700,  -8/9,  4/3 ]`

`[  52/27,  5617869293/21003948000, -30905427529/10501974000,`
`  -1174937411/21003948000, 55109414029/21003948000, 32/81, -16/27 ]`.

The matrix still has rank `2`, so the full even-gate intersection is still
5-dimensional. That part is unsurprising: there are still only two even equations.

The new structural fact is the mixed-sector minor:

`det Gate_(alpha_W, alpha_R) = 942737330573 / 205838690400000 != 0`.

So on this branch the full even system can solve directly for the previously
untouched mixed directions:

`alpha_W =`
`  (14503089433000/942737330573) alpha_K`
`+ (30450672110098/942737330573) alpha_OmegaU`
`- (29120459867142/942737330573) alpha_U`
`- (18876066395200/25453907925471) beta_B`
`+ (9438033197600/8484635975157) beta_varpi`,

`alpha_R =`
`  (101802968743000/942737330573) alpha_K`
`+ (189815725996721/942737330573) alpha_OmegaU`
`- (188832473718440/942737330573) alpha_U`
`+ (89510801038400/25453907925471) beta_B`
`- (44755400519200/8484635975157) beta_varpi`.

So there are no longer pure `alpha_W` or pure `alpha_R` null directions on the full
constructive branch. The mixed-sector freedom was never genuinely unconstrained; it
was just hidden in the omitted `Z_2,Z_4` block.

## Net result after Stage 20

The continuation point is sharper again.

1. The isotropic full-bundle target surface is now exact and explicit.
2. The omitted conservative `Z_2,Z_4` sector has been reinstated exactly.
3. On a clean constructive branch it does the precise job Stage 17 predicted:
   it removes the fake freedom in the mixed directions `alpha_W, alpha_R`.
4. So the remaining work is no longer “some mixed-sector freedom somewhere.”
   It is to compute the actual overlap data on the true moving-throat branch and see
   whether they land on the isotropic full-bundle target surface from Stage 18.
# 5PN continuation notes — Stages 21–23 and 27–29

This note fills the earlier numbering gap and records the next continuation past the Stage 26 placement map.

## Missing stages 21–23

### Stage 21 — Dimensionless continuum placement map
The continuum-selected quadrupole placement problem compresses to the exact dimensionless ratios

a) `eps_eta = c_(etaU)^2 / (K_U K_eta^(eff))`

b) `eps_W = c_(UW)^2 sigma / (K_U K_W^(eff))`

c) `rho = c_(UW) c_(etaU) / (K_U c_(etaW))`

d) `Z_W = c_(etaW)^2 / (K_eta^(eff) K_W^(eff))`

e) `delta_0 = pi^2 T_w / (L^2 K_eta^(eff))`

with radiative scale

`Lambda = 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W)`.

The exact placement coordinates are

`delta = delta_0 / (1 - eps_eta)`,

`M_mix = 8 Z_W (1 + rho)^2 / [pi^2 (1 - eps_eta)(1 - eps_W)]`,

`R_target = Lambda (1 - eps_eta)(1 - eps_W)^2 / [Z_W (1 + rho)^2]`.

The key exact product law is

`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`

`= 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W)`.

So Stage 21 factorizes the continuum map into a geometry lane, a mixed-stability/product lane, and a redistribution lane.

### Stage 22 — First non-flat U doublet
Turning on the first axial structure of the internal `U` sector preserves the scalar placement map but breaks the old one-direction geometry.

The new exact quantities are

`delta_split = [delta_0 + eps_eta delta_U/(1+delta_U)] / (1-eps_eta)`,

`eps_W_split = eps_W [1 - (2/11) delta_U/(1+delta_U)]`,

`R_U = [1 + rho_0/(1+delta_U)] / (1 + rho_0)`.

The exact direction-splitting invariant is

`D_dir = - kappa_0 kappa_1 g_W rho_0 delta_U / (1+delta_U)`.

So collinearity survives iff `delta_U = 0` or `rho_0 = 0`.

### Stage 23 — Generalized selected-branch normalization
With distinct loading vector `z` and source vector `s`, the selected branch survives exactly, but the old flat-`U` function pair `(F,G)` deforms to

`F_(q,eta)(xi,delta)` and `G_q(xi,delta)`.

For the split-`U` continuum, the whole deformation collapses to the single ratio `R_U`, giving the exact one-parameter family

`F_U(xi,delta;R_U)`, `G_U(xi,delta;R_U)`.

Setting `R_U = 1` recovers the Stage 18–19 flat branch exactly.

The exact first deformation about the flat branch is

`F_U/F_flat = 1 + eps H_F + O(eps^2)`,

`G_U/G_flat = 1 + eps H_G + O(eps^2)`,

with the closed-form `H_F` and `H_G` verified in the companion SymPy audit.

## Continuation past Stage 26: Stages 27–29

### Stage 27 — Continuum-selected rank-2 closure
Once the actual support-direction data from Stage 26 are inserted, the physical selected wall branch is pinned by an exact quadratic for the softening depth `xi`:

`xi^2 + B_cont xi + C_cont = 0`,

with

`B_cont = delta - M_mix (1+t^2 R_U^2) - M_supp (1+t^2 R_phi^2)`,

`C_cont = - delta (M_mix + M_supp) + t^2 M_mix M_supp (R_U - R_phi)^2`.

So the continuum-selected normalization theorem gate becomes simply

`R_target = F_cont(xi_phys)`.

The exact special surfaces are:

- minimal-kernel source-tied surface: `sigma_0 = 0`,
- interference-matched tracking surface: `sigma_0 = rho_0`.

### Stage 28 — Coherent local D/N support kernel
The first coherent local D/N kernel forces

`g_B g_R = g_W g_S`,

so it lands exactly on the tracking surface. Therefore

`R_phi = R_U = R_tr`,

and the full Stage 27 rank-2 closure collapses to the one-parameter tracking laws

`M_tr = G_tr(xi,delta;R_tr)`,

`R_target = F_tr(xi,delta;R_tr)`.

The coherent branch also has the exact constructive bound

`1/(1+delta_U) < R_tr < 1`.

### Stage 29 — Tracking-branch bounds and residual comparison
On the coherent tracking branch, at fixed `(xi,delta)`,

`dG_tr/dR < 0`,

`dF_tr/dR > 0`.

So for the constructive split-`U` branch, where `R_tr < 1`, one has the exact inequalities

`G_tr > G_flat`,

`F_tr < F_flat`.

The endpoint formulas are

`G_tr(R=0) = xi`,

`F_tr(R=0) = 1/(1-xi)`.

The exact residual comparison theorem is

`E_tr - E_flat = F_flat - F_tr > 0`.

So the first coherent local D/N kernel does **not** rescue the Stage 18–19 normalization target. It makes the target harder: more total loading is required, while the normalized selected-branch response is lower than on the flat branch.

## Bottom line after this batch

The selected-branch problem is now narrower than before:

1. the physical softening depth is set by an exact quadratic;
2. the first coherent local D/N kernel forces exact tracking rather than a generic intermediate rank-2 closure;
3. and the resulting split-`U` deformation worsens the normalization residual relative to the old flat branch.

So the next theorem gate is no longer “which closure should we use?” It is whether the later support-enhancement / compensation mechanisms in the notes are strong enough to overcome this exact tracking-branch deficit.

## Further continuation — Stages 30–31

### Stage 30 — Coherent-kernel dimensionless map
The first coherent local D/N kernel compresses to a very small exact parameter set:

- `eps_eta`
- `delta_U`
- `chi_0`
- `eps_W`
- `Z_W`
- `zeta`
- `Lambda`

The exact coherent tracking factor is

`R_tr = [1 + chi_0/(1+delta_U)] / (1 + chi_0)`.

The total baseline splits as

`M_mix = 8 Z_W (1 + chi_0)^2 / [pi^2 (1 - eps_eta)(1 - eps)]`,

`M_supp = 8 zeta Z_W (1 + chi_0)^2 / [pi^2 (1 - eps_eta)(1 - zeta eps)]`,

with

`eps = eps_W [1 - (2/11) delta_U/(1+delta_U)]`.

So the support lane enters only through the exact enhancement factor

`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps)`,

and

`M_tr = M_mix S(zeta;eps)`.

The target ratio remains independent of `zeta`:

`R_target = Lambda (1 - eps_eta)(1 - eps)^2 / [Z_W (1 + chi_0)^2]`.

### Stage 31 — Support-compensation theorem
On the coherent tracking branch,

`M_tr = G_tr(xi,delta;R_tr)`

with critical load

`M_crit(delta,R_tr) = G_tr(1,delta;R_tr)`.

The support enhancement factor is strictly increasing and exactly invertible:

`zeta_req = (S_req - 1) / [1 + eps (S_req - 2)]`.

For every finite target ratio `R_target > 1`, there exists a stable-side branch point `xi_req in (0,1)` and a corresponding required load `M_req < M_crit`.

Therefore, if the mixed-only branch is below that required load, there is a unique coherent support ratio `zeta_req < zeta_crit < 1/eps` that reaches the target **before** the branch softens out.

So there is no reduced-level support no-go. The remaining question is whether the actual PDE produces a physical `zeta` large enough to meet that exact threshold.

## Further continuation — Stages 32–33

### Stage 32 — Explicit D/N overlap extraction of `zeta`
The coherent support ratio is no longer phenomenological. For the first explicit finite-throat D/N family,

`zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) (I_n / I_0)^2`,

with

`I_n = ∫_0^L ds sigma(s) chi_n(s)`,

`chi_n(s) = sqrt(2/L) sin[(n+1/2) pi s / L]`.

For the first uniform local source density `sigma(s)=1`,

`I_n / I_0 = 1 / (2n+1)`.

So the physical coherent support ratio becomes

`zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2`.

On the same-operator twin family,

`zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ]`.

The exact lowest-twin value is

`zeta_0^(twin) = 1`.

### Stage 33 — Exact comparison of `zeta_phys` against `zeta_req`
Because Stage 31 reduced support feasibility to `zeta_phys >= zeta_req`, the explicit D/N tower gives exact selection rules.

The lowest symmetric twin lane has

`zeta_0^(twin) = 1`,

so its enhancement factor is exactly

`S_0 = 2`,

independent of `eps`.

Therefore the lowest symmetric twin lane succeeds exactly when

`zeta_req <= 1`,

equivalently

`S_req <= 2`.

For higher D/N harmonics,

`zeta_n^(twin) < 1/(2n+1)^2`,

so they are ruled out immediately if

`zeta_req > 1/(2n+1)^2`.

When they are not ruled out immediately, the exact softness threshold is

`x <= x_max(n; zeta_req)`

with

`x_max(n; zeta_req) = [1 / ((2n+1)^2 zeta_req) - 1] / [n(n+1)]`.

So the explicit coherent support tower is strongly biased toward the lowest twin lane.

## Bottom line after Stages 21–33

The program is now much sharper than when the numbering gap first appeared.

1. Stages 21–23 show that split-`U` physics preserves the scalar placement map but deforms the selected-branch normalization geometry through an exact direction-splitting parameter.
2. Stages 27–29 show that once the actual continuum support data are inserted, the physical wall softening depth is fixed by an exact quadratic, and the first coherent local D/N kernel lands on a tracking branch whose split-`U` deformation worsens the normalization target relative to the old flat branch.
3. Stages 30–31 show there is no reduced-level support no-go: coherent support enhancement can in principle compensate the tracking-branch deficit before softening.
4. Stages 32–33 then make that compensation test microscopic: the lowest symmetric twin D/N lane is an exact universal doubling branch, and every higher support harmonic is strongly overlap- and stiffness-suppressed.

So the next really important question is no longer which closure to use. It is whether the actual moving-throat PDE places the physical coherent support sector on the lowest twin branch with a large enough `zeta_phys` to satisfy the exact threshold.
# 5PN / Moving-Throat continuation — Stages 24–26

These stages continue the explicit moving-throat spectral-placement branch after the concrete axial / loaded-profile / selected-prefactor scripts.

## Stage 24 — source map, microscopic normalization, softening depth

Files:
- `5pn_stage24_source_map_microscopic_normalization_softening.py`
- `5pn_stage24_source_map_microscopic_normalization_softening_output.txt`

Main results:

1. The natural D/N source-map factor is no longer independent:
   
   \[
   \hat m_-^2 = \frac{s_-(x)}{\kappa_0^2}.
   \]

2. The selected-branch normalization product becomes
   
   \[
   N_-(x) = \frac{\beta_0\, s_-(x)^2}{\kappa_0^2 (A-x)}.
   \]

3. The selected branch is exactly parameterized by the softening depth
   \[
   x = A-\lambda_-.
   \]

4. The exact secular loading law is
   
   \[
   \alpha_0(x)=\frac{x(x+\Delta K_{ax})}{\kappa_0^2(x+\Delta K_{ax})+\kappa_1^2 x}.
   \]

5. The exact support-loading requirement is
   
   \[
   \frac{g_{B,req}^2}{\varpi^2}=\alpha_0(x)-\alpha_{mix}.
   \]

So the selected quadrupole bridge is no longer an eigenvector problem. It is a one-variable spectral-placement problem in `x`.

---

## Stage 25 — universal D/N branch geometry

Files:
- `5pn_stage25_dimensionless_normalization_and_support_frontier.py`
- `5pn_stage25_dimensionless_normalization_and_support_frontier_output.txt`

Main results:

Using the exact D/N constants
\[
\kappa_0^2=\frac{8}{\pi^2},\qquad
\kappa_1^2=\frac{16}{9\pi^2},\qquad
\eta=\frac{\kappa_1^2}{\kappa_0^2}=\frac29,
\]
and the dimensionless variables
\[
\xi=\frac{x}{A},\qquad \delta=\frac{\Delta K_{ax}}{A},
\]
we get two universal branch functions:

1. The normalization function
   
   \[
   F(\xi,\delta)=
   \frac{(9\delta+11\xi)^4}{81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2}.
   \]

2. The support-feasibility function
   
   \[
   G(\xi,\delta)=\frac{9\xi(\xi+\delta)}{9\delta+11\xi}.
   \]

Exact geometric facts verified in SymPy:

- `F` is strictly increasing on the stable branch,
- `F(0,delta)=1`,
- `F -> +infinity` as `xi -> 1^-`,
- `G` is strictly increasing,
- `G(0,delta)=0`,
- `G_max(delta)=9(1+delta)/(9 delta + 11)`.

So for fixed `delta`, the selected branch is controlled by a unique normalization locus `F = R_target` and a sharp support-feasibility frontier `M_mix <= G`.

---

## Stage 26 — continuum-kernel extraction and placement map

Files:
- `5pn_stage26_continuum_kernel_extraction_and_placement_map.py`
- `5pn_stage26_continuum_kernel_extraction_and_placement_map_output.txt`

Main results:

From the first explicit finite-throat continuum kernel, the reduced branch data are extracted exactly:

\[
A=\frac{K_U K_{\eta}^{eff}-c_{\eta U}^2}{\mu_\eta K_U},
\qquad
\Delta K_{ax}=\frac{\pi^2 T_w}{\mu_\eta L^2},
\]
\[
\alpha_{mix}=
\frac{(K_U c_{\eta W}+c_{UW} c_{\eta U})^2}
{\mu_\eta K_U (K_U K_W^{eff}-c_{UW}^2\sigma)},
\]
\[
\beta_0=
\frac{\mu_W}{\mu_\eta}
\frac{(K_U c_{\eta W}+c_{UW} c_{\eta U})^2}
{(K_U K_W^{eff}-c_{UW}^2\sigma)^2}.
\]

Then the full placement problem compresses to the exact dimensionless kernel ratios

\[
\epsilon_\eta = \frac{c_{\eta U}^2}{K_U K_\eta^{eff}},
\qquad
\epsilon_W = \frac{c_{UW}^2\sigma}{K_U K_W^{eff}},
\]
\[
\rho = \frac{c_{UW} c_{\eta U}}{K_U c_{\eta W}},
\qquad
Z_W = \frac{c_{\eta W}^2}{K_\eta^{eff} K_W^{eff}},
\qquad
\delta_0 = \frac{\pi^2 T_w}{L^2 K_\eta^{eff}},
\]
plus the radiative demand scale
\[
\Lambda = \frac{27 \pi^2 G c_s^5 K_W^{eff}}{20 a^5 c^5 \mu_W}.
\]

The exact placement formulas are

\[
\delta = \frac{\delta_0}{1-\epsilon_\eta},
\qquad
M_{mix} = \frac{8 Z_W (1+\rho)^2}{\pi^2 (1-\epsilon_\eta)(1-\epsilon_W)},
\]
\[
R_{target} = \frac{\Lambda (1-\epsilon_\eta)(1-\epsilon_W)^2}{Z_W (1+\rho)^2}.
\]

And the exact product law is

\[
R_{target} M_{mix} = \frac{8\Lambda(1-\epsilon_W)}{\pi^2}.
\]

This cleanly separates the placement problem into three lanes:

1. geometry lane: `delta = delta_0 / (1-eps_eta)`
2. mixed-stability/product lane: `R_target M_mix = 8 Lambda (1-eps_W)/pi^2`
3. redistribution lane: `(eps_eta, Z_W, rho)` move the branch along the fixed product curve.

---

## What these stages change in the 5PN program

The remaining bridge is no longer a vague “solve more PDE” task.
It is now a very narrow branch-placement problem:

1. compute the actual continuum ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)`,
2. map them to `(delta, M_mix, R_target)`,
3. find the unique `xi` solving `F(xi,delta)=R_target`,
4. and check the support-feasibility condition `M_mix <= G(xi,delta)`.

If that succeeds, the selected quadrupole branch is admissible on the natural finite-throat continuum placement map.
# 5PN stages 34–40 notes

This bundle continues the post-Stage-33 support/normalization program by turning the
`zeta_req` threshold into exact operator criteria on the moving-throat branch.

## Stage 34 — exact lowest-twin sufficiency criterion

Using the tracking-branch functions
\[
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]
\[
F_{\rm tr}(\xi,\delta;R)
=
\frac{\bigl[9\delta+(9+2R^2)\xi\bigr]^2\bigl[9\delta+(9+2R)\xi\bigr]^2}
{81(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2},
\]
the exact product is
\[
\Pi_{\rm tr}=F_{\rm tr}G_{\rm tr}
=
\frac{\xi(\xi+\delta)\bigl[9\delta+(9+2R)\xi\bigr]^2\bigl[9\delta+(9+2R^2)\xi\bigr]}
{9(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2}.
\]

With
\[
C_{\rm mix}:=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\qquad
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}},
\]
the symmetric lowest twin lane is sufficient iff
\[
\Pi_{\rm tr}(\xi_{\rm req},\delta;R_{\rm tr})
\le
\frac{16\Lambda(1-\epsilon)}{\pi^2}.
\]

Equivalent threshold scales:
\[
\Lambda_{\rm twin,req}=\frac{\pi^2}{16(1-\epsilon)}\Pi_{\rm tr},
\qquad
M_{\rm mix}^{(\rm twin,req)}=\frac{G_{\rm tr}}{2},
\]
\[
Z_W^{(\rm twin,req)}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)\,G_{\rm tr}}
{16(1+\chi_0)^2}.
\]

The exact twin-saturation depth at fixed mixed baseline is the unique root of
\[
G_{\rm tr}(\xi_{2\times},\delta;R)=2M_{\rm mix},
\]
namely
\[
\xi_{2\times}
=
\frac{2M_{\rm mix}(9+2R^2)-9\delta+\sqrt{(2M_{\rm mix}(9+2R^2)-9\delta)^2+648M_{\rm mix}\delta}}{18}.
\]

## Stage 35 — exact non-twin asymmetry requirement

Define
\[
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}
=
\frac{\Pi_{\rm tr}-C_{\rm mix}}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]
Then
\[
\frac{d\zeta_{\rm req}}{d\Pi_{\rm tr}}
=
\frac{C_{\rm mix}(1-\epsilon)}{\bigl[C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})\bigr]^2}>0,
\]
so the required coherent support ratio grows monotonically with the branch product.

Exact regime split:
\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\Rightarrow\quad
\text{mixed-only enough},
\]
\[
C_{\rm mix}<\Pi_{\rm tr}\le2C_{\rm mix}
\quad\Rightarrow\quad
\text{symmetric lowest twin enough},
\]
\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\Rightarrow\quad
\text{non-twin asymmetry required}.
\]

The exact excess beyond the symmetric twin branch is
\[
\Delta_\zeta:=\zeta_{\rm req}-1
=
\frac{(1-\epsilon)(\Pi_{\rm tr}-2C_{\rm mix})}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]

For a general lowest support lane,
\[
\zeta_0^{(\rm phys)}
=
\frac{K_W^{(\rm eff)}}{K_{\phi,0}^{(\rm eff)}}\Omega_0^2.
\]
So the two equivalent exact rescue thresholds are
\[
\Omega_0^2 \ge \zeta_{\rm req}\frac{K_{\phi,0}^{(\rm eff)}}{K_W^{(\rm eff)}},
\qquad
K_{\phi,0}^{(\rm eff)} \le K_W^{(\rm eff)}\frac{\Omega_0^2}{\zeta_{\rm req}}.
\]

## Stage 36 — exact overlap-boost window

For the D/N lowest support mode
\[
\chi_0(s)=\sqrt{\frac{2}{L}}\sin\!\frac{\pi s}{2L},
\qquad
I_W=\int_0^L\chi_0(s)\,ds=\frac{2\sqrt{2L}}{\pi},
\]
and the normalized exponential source family
\[
\sigma_\alpha(s)=\frac{\alpha e^{\alpha s/L}}{e^\alpha-1},
\qquad
\int_0^L \sigma_\alpha(s)\,ds=L,
\]
the overlap boost is
\[
\Omega_{\exp}(\alpha)
=
\frac{\int_0^L \sigma_\alpha(s)\chi_0(s)\,ds}{I_W}
=
\frac{\pi\alpha\bigl(2\alpha e^\alpha+\pi\bigr)}
{(4\alpha^2+\pi^2)(e^\alpha-1)}.
\]

Exact endpoint values:
\[
\Omega_{\exp}(0)=1,
\qquad
\lim_{\alpha\to+\infty}\Omega_{\exp}(\alpha)=\frac{\pi}{2}.
\]
Therefore
\[
0\le \Omega_0\le\frac{\pi}{2},
\qquad
A_I:=\Omega_0^2\le\frac{\pi^2}{4}.
\]

So pure overlap rescue alone is possible only if
\[
\zeta_{\rm req}\le\frac{\pi^2}{4}.
\]

## Stage 37 — Robin-compliance softening

Replacing the Dirichlet mouth by a Robin compliance gives the lowest-lane eigenvalue condition
\[
y\tan y=\eta,
\qquad
0<y<\frac{\pi}{2},
\]
with
\[
\eta:=hL.
\]

If
\[
x:=\frac{\pi^2 T_X}{L^2K_W^{(\rm eff)}},
\qquad
0<x<4,
\]
then the exact support-softening factor is
\[
A_K(\eta)
=
\frac{K_W^{(\rm eff)}}{K_{\phi,0}^{(\rm eff)}}
=
\frac{1}{1-x/4+xy^2/\pi^2}.
\]

Endpoint window:
\[
A_K=\;1\;\text{at}\;y=\frac{\pi}{2},
\qquad
A_{K,\max}=\frac{4}{4-x}\;\text{at}\;y\to0.
\]

So pure support softening alone can rescue the Stage-35 threshold only if
\[
\zeta_{\rm req}\le\frac{4}{4-x}.
\]

At fixed \(\zeta_{\rm req}\), the exact eigenvalue and Robin thresholds are
\[
y_{\rm req}^2=\frac{\pi^2}{x}\left(\frac{1}{\zeta_{\rm req}}-1+\frac{x}{4}\right),
\qquad
\eta_{\rm req}=y_{\rm req}\tan y_{\rm req}.
\]

## Stage 38 — explicit non-twin lowest-lane reachability

Combining the Stage-36 overlap family with Stage-37 Robin softening gives
\[
\zeta_0^{(\exp+R)}(\alpha,\eta)
=
\frac{\Omega_{\exp}(\alpha)^2}{1-x/4+xy(\eta)^2/\pi^2}.
\]

Its exact closure range is
\[
1\le \zeta_0^{(\exp+R)} \le \frac{\pi^2}{4-x}.
\]

So the explicit family reaches the Stage-35 threshold iff
\[
\zeta_{\rm req}\le\frac{\pi^2}{4-x}.
\]

This gives the exact three-regime split:
\[
\zeta_{\rm req}\le\frac{\pi^2}{4}
\quad\Rightarrow\quad
\text{overlap alone enough},
\]
\[
\frac{\pi^2}{4}<\zeta_{\rm req}\le\frac{\pi^2}{4-x}
\quad\Rightarrow\quad
\text{overlap + softening enough},
\]
\[
\zeta_{\rm req}>\frac{\pi^2}{4-x}
\quad\Rightarrow\quad
\text{even this explicit family fails}.
\]

## Stage 39 — transport origin of source asymmetry

The Stage-36 exponential family is exactly the stationary zero-flux branch of
\[
\partial_t \sigma + \partial_s J = 0,
\qquad
J=-D_\sigma \partial_s \sigma + v_\sigma \sigma.
\]

On the stationary recirculating branch \(J=0\),
\[
\sigma_{Pe}(s)
=
\frac{Pe\,e^{Pe s/L}}{e^{Pe}-1},
\qquad
Pe:=\frac{v_\sigma L}{D_\sigma}.
\]

The corresponding overlap boost is
\[
\Omega_{Pe}
=
\frac{\pi Pe\bigl(2Pe\,e^{Pe}+\pi\bigr)}
{(4Pe^2+\pi^2)(e^{Pe}-1)}.
\]

Exact identities:
\[
\Omega_{Pe}(0)=1,
\qquad
\lim_{Pe\to+\infty}\Omega_{Pe}=\frac{\pi}{2},
\]
and the score-function identity
\[
\partial_{Pe}p_{Pe}(x)=(x-\mathbb E_{Pe}[x])p_{Pe}(x)
\]
implies
\[
\frac{d}{dPe}\mathbb E_{Pe}[\chi_0]=\operatorname{Cov}_{Pe}(\chi_0,x),
\]
so the constructive branch is monotone increasing because \(\chi_0\) is increasing on \([0,1]\).

## Stage 40 — physical \((Pe,\kappa,\eta)\) placement map

Define the physical support ratios
\[
\kappa:=\frac{K_XL^2}{T_X},
\qquad
\eta:=hL=\frac{K_mL}{T_X}.
\]
Then
\[
x=\frac{\pi^2}{\kappa+\pi^2/4},
\]
and the Robin softening factor becomes
\[
A_K(\eta;\kappa)=\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2}.
\]

The explicit physical lowest-lane family is therefore
\[
\zeta_0^{(Pe+R)}(Pe,\eta;\kappa)
=
\Omega_{Pe}^2\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2}.
\]

This map is monotone:
- increasing in \(Pe\),
- decreasing in \(\eta\),
- decreasing in \(\kappa\).

Its exact constructive-branch ceiling is
\[
\zeta_{\max}(\kappa)=\frac{\pi^2}{4}\frac{\kappa+\pi^2/4}{\kappa}.
\]

So the Stage-35 demand is reachable on this first physical family iff
\[
\zeta_{\rm req}\le \zeta_{\max}(\kappa),
\]
equivalently, whenever \(\zeta_{\rm req}>\pi^2/4\),
\[
\kappa \le \kappa_{\max}(\zeta_{\rm req})
:= \frac{\pi^4}{4(4\zeta_{\rm req}-\pi^2)}.
\]

The exact physical threshold surfaces are
\[
\Omega_{\rm req}^2
=
\zeta_{\rm req}\frac{\kappa+y(\eta)^2}{\kappa+\pi^2/4},
\]
\[
y_{\rm req}^2
=
\frac{\Omega_{Pe}^2}{\zeta_{\rm req}}(\kappa+\pi^2/4)-\kappa,
\]
\[
\kappa_{\rm req}
=
\frac{\Omega_{Pe}^2\pi^2/4-\zeta_{\rm req}y(\eta)^2}{\zeta_{\rm req}-\Omega_{Pe}^2}.
\]

## Where the program stands after Stage 40

The support/normalization problem is no longer phrased in abstract deformation variables.
It has collapsed to three physical moving-throat operator ratios:

- axial source Peclet number `Pe`,
- mouth compliance `eta`,
- baseline support stiffness ratio `kappa`.

That makes the next clean move very sharp: derive the coupled support/source branch equation that selects `Pe`
from the same operator that already fixes `eta` and `kappa`.
# 5PN stages 41–51 notes

This bundle continues the post-Stage-40 support/normalization program by turning the
physical placement map into an operator-selected branch law and then projecting the
remaining theorem gap back into explicit parent-action overlap data.

## Stage 41 — coupled support/source operator and exact `Pe` branch equation

The first coupled axial operator is

a) source transport
\[
\partial_t\sigma + \partial_s J = 0,
\qquad
J = -D_\sigma \partial_s\sigma + v_\sigma \sigma,
\]

b) support field
\[
-T_X \partial_s^2 \phi + K_X \phi = \Lambda_\phi \sigma,
\]
with support boundary conditions
\[
T_X \phi_s(0)=K_m\phi(0),
\qquad
\phi_s(L)=0.
\]

After nondimensionalization,
\[
\kappa = K_X L^2/T_X,
\qquad
\eta = K_m L/T_X,
\qquad
Pe = v_\sigma L/D_\sigma,
\qquad
\Xi = \mu_\sigma \Lambda_\phi^2 L^2/(D_\sigma T_X).
\]

On the stationary zero-flux branch,
\[
\Sigma_{Pe}(x)=\frac{Pe\,e^{Pe x}}{e^{Pe}-1},
\qquad x=s/L\in[0,1].
\]
The exact support-drop kernel is
\[
K_{\kappa,\eta}(x)=\frac{\cosh(\alpha x)+(\eta/\alpha)\sinh(\alpha x)-\cosh(\alpha(1-x))}{\alpha\sinh\alpha+\eta\cosh\alpha},
\qquad \alpha=\sqrt\kappa,
\]
with derivative
\[
\frac{dK_{\kappa,\eta}}{dx}=
\frac{\alpha\sinh(\alpha x)+\eta\cosh(\alpha x)+\alpha\sinh(\alpha(1-x))}{\alpha\sinh\alpha+\eta\cosh\alpha}>0.
\]

So the exact dimensionless support drop is
\[
\Delta(Pe;\kappa,\eta)=\int_0^1 dx\;K_{\kappa,\eta}(x)\Sigma_{Pe}(x),
\]
with endpoint values
\[
\Delta_0(\kappa,\eta)=
\frac{\eta(\cosh\alpha-1)}{\alpha^2(\alpha\sinh\alpha+\eta\cosh\alpha)},
\]
\[
\Delta_\infty(\kappa,\eta)=
\frac{\cosh\alpha+(\eta/\alpha)\sinh\alpha-1}{\alpha\sinh\alpha+\eta\cosh\alpha}.
\]

The branch point is therefore selected by the exact fixed-point equation
\[
Pe = \Xi\,\Delta(Pe;\kappa,\eta),
\]
and every constructive branch root obeys
\[
\Xi\Delta_0(\kappa,\eta)
\le Pe_* \le
\Xi\Delta_\infty(\kappa,\eta).
\]
At weak coupling,
\[
Pe_* = \Xi \Delta_0(\kappa,\eta)+O(\Xi^2).
\]

## Stage 42 — exact residual bounds on the operator-selected branch

Evaluating the Stage-40 physical support ratio on the branch root gives
\[
\zeta_{\rm phys}(\Xi,\eta;\kappa)=
\Omega_{Pe_*}^2\,\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2},
\qquad y\tan y=\eta.
\]
Since \(\Omega_{Pe}\) is strictly increasing, the Stage-41 branch interval gives the exact support bracket
\[
\zeta_-(\Xi,\eta;\kappa)
\le \zeta_{\rm phys}(\Xi,\eta;\kappa)
\le \zeta_+(\Xi,\eta;\kappa),
\]
where
\[
\zeta_-=
\Omega_{\Xi\Delta_0}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2},
\qquad
\zeta_+=
\Omega_{\Xi\Delta_\infty}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2}.
\]

Inside the Stage-40 reachability window, define the unique constructive point \(Pe_{\rm req}\) by
\[
\Omega_{Pe_{\rm req}}^2=
\zeta_{\rm req}\,\frac{\kappa+y^2}{\kappa+\pi^2/4}.
\]
Then the exact coupling thresholds are
\[
\Xi_{\rm fail}=\frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)},
\qquad
\Xi_{\rm suff}=\frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\qquad
\Xi_{\rm fail}\le \Xi_{\rm suff}.
\]
So the branch has a sharp three-zone structure:

- \(\Xi\le\Xi_{\rm fail}\): impossible,
- \(\Xi\ge\Xi_{\rm suff}\): guaranteed,
- \(\Xi_{\rm fail}<\Xi<\Xi_{\rm suff}\): only here is the full root solve needed.

The exact residual envelope is
\[
R_- \le R_{\rm phys} \le R_+,
\qquad
R_{\rm phys}=\zeta_{\rm req}-\zeta_{\rm phys},
\]
with
\[
R_- = \zeta_{\rm req}-\zeta_+,
\qquad
R_+ = \zeta_{\rm req}-\zeta_-.
\]

Using
\[
\Omega_{Pe}^2 = 1 + \frac{4-\pi}{\pi}Pe + O(Pe^2),
\]
one gets the weak-coupling branch law
\[
\zeta_{\rm phys}=
A_K(\eta;\kappa)
\Bigl[1+\frac{4-\pi}{\pi}\,\Xi\Delta_0(\kappa,\eta)+O(\Xi^2)\Bigr].
\]

## Stage 43 — entropic source microclosure and microscopic gain

The first explicit source/support free energy is
\[
F[\sigma,\phi]=
\int_0^L ds\Bigl[
\Theta_\sigma\sigma(\log(\sigma/\sigma_*)-1)-\Lambda_\phi\sigma\phi
+\frac{T_X}{2}\phi_s^2+\frac{K_X}{2}\phi^2
\Bigr]+
\frac{K_m}{2}\phi(0)^2.
\]
Its exact variations are
\[
\mu_\sigma^{\rm chem}=
\frac{\delta F}{\delta\sigma}=
\Theta_\sigma\log(\sigma/\sigma_*)-\Lambda_\phi\phi,
\]
\[
-T_X\phi_{ss}+K_X\phi=\Lambda_\phi\sigma,
\qquad
T_X\phi_s(0)=K_m\phi(0),
\qquad
\phi_s(L)=0.
\]

The minimal positive-density Onsager current is
\[
J=-M_\sigma\sigma\partial_s\mu_\sigma^{\rm chem}
  =-D_\sigma\partial_s\sigma + M_\sigma\Lambda_\phi\sigma\partial_s\phi,
\]
with exact Einstein relation
\[
D_\sigma=M_\sigma\Theta_\sigma.
\]

Under the affine-drop reduction
\[
\phi(s)\approx \phi(0)+[\Delta\phi]s/L,
\]
the stationary zero-flux branch gives the exact exponential family
\[
\sigma(s)=C\exp\!igl[(\Lambda_\phi\Delta\phi)/(\Theta_\sigma L)\,s\bigr],
\]
so
\[
Pe=(\Lambda_\phi/\Theta_\sigma)\Delta\phi.
\]
Using
\[
\Delta\phi=(\Lambda_\phi L^2/T_X)\Delta(Pe;\kappa,\eta),
\]
one gets the exact microscopic coupling
\[
\Xi_{\rm micro}=
\frac{\Lambda_\phi^2 L^2}{\Theta_\sigma T_X}
=
\chi_\sigma\frac{\Lambda_\phi^2 L^2}{T_X},
\qquad
\chi_\sigma:=1/\Theta_\sigma.
\]

The closure is automatically passive:
\[
\frac{dF}{dt}=-\int_0^L ds\;\frac{J^2}{M_\sigma\sigma}\le 0
\]
under no-flux boundaries.

## Stage 44 — microscopic gain thresholds and exact phase diagram

Using \(\kappa=K_XL^2/T_X\), the Stage-43 coupling becomes
\[
\Xi_{\rm micro}=\kappa G_{\rm micro},
\qquad
G_{\rm micro}:=\chi_\sigma\Lambda_\phi^2/K_X.
\]
So the exact microscopic thresholds are
\[
G_{\rm fail}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_\infty(\kappa,\eta)},
\qquad
G_{\rm suff}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_0(\kappa,\eta)}.
\]
The exact phase diagram is:

- \(G_{\rm micro}\le G_{\rm fail}\): impossible,
- \(G_{\rm micro}\ge G_{\rm suff}\): guaranteed,
- only the bounded interval in between needs the full root solve.

Equivalent threshold surfaces are
\[
\chi_\sigma \le \frac{T_X Pe_{\rm req}}{\Lambda_\phi^2 L^2 \Delta_\infty}\quad\Rightarrow\quad\text{fail},
\]
\[
\chi_\sigma \ge \frac{T_X Pe_{\rm req}}{\Lambda_\phi^2 L^2 \Delta_0}\quad\Rightarrow\quad\text{succeed},
\]
and similarly for \(\Lambda_\phi^2\).

Soft-support limit:
\[
\Delta_0\to \frac12,
\qquad
\Delta_\infty\to 1,
\qquad
G_{\rm fail}\sim \frac{Pe_{\rm req}}{\kappa},
\qquad
G_{\rm suff}\sim \frac{2Pe_{\rm req}}{\kappa}.
\]
So very soft support is strongly disfavored.

Highly compliant-mouth limit:
\[
\Delta_0^{(\infty)}=
\frac{1-\operatorname{sech}(\sqrt\kappa)}{\kappa},
\qquad
\Delta_\infty^{(\infty)}=
\frac{\tanh(\sqrt\kappa)}{\sqrt\kappa},
\]
so
\[
G_{\rm fail}^{(\infty)}=
\frac{Pe_{\rm req}}{\sqrt\kappa\tanh(\sqrt\kappa)},
\qquad
G_{\rm suff}^{(\infty)}=
\frac{Pe_{\rm req}}{1-\operatorname{sech}(\sqrt\kappa)}.
\]
For \(\kappa\gg1\), these reduce to
\[
G_{\rm fail}^{(\infty)}\sim \frac{Pe_{\rm req}}{\sqrt\kappa},
\qquad
G_{\rm suff}^{(\infty)}\sim Pe_{\rm req}.
\]

## Stage 45 — parent-action projection of the microscopic gain

Starting from the parent matter energy
\[
H_\psi=
\int d^4X\left[\frac{\hbar^2}{2m}|D_i\psi|^2+V_{\rm conf}\rho+U(\rho)\right],
\]
with frozen EOS
\[
U(\rho)=K\rho^5/4,
\qquad
h(\rho)=\frac{5K}{4}\rho^4,
\qquad
h'(\rho)=5K\rho^3=\frac{m c_s^2(\rho)}{\rho},
\]
the local compressional quadratic energy is
\[
\delta H_{\rm comp}=
\frac12\int d^4X\;h'(\rho_*)(\delta\rho)^2.
\]

Project one source channel
\[
\delta\rho(s,y)=\sigma(s)\chi_\sigma(y)
\]
and one support channel entering the confinement as
\[
\delta V_{\rm conf}(s,y)=-g_\phi\chi_\phi(y)\phi(s).
\]
Then the exact reduced coefficients are
\[
\Theta_\sigma=h'(\rho_*)N_{\sigma\sigma},
\qquad
\Lambda_\phi=g_\phi O_{\sigma\phi},
\]
with parent overlap invariants
\[
N_{\sigma\sigma}=\int d^3y\,\chi_\sigma^2,
\qquad
N_{\phi\phi}=\int d^3y\,\chi_\phi^2,
\qquad
O_{\sigma\phi}=\int d^3y\,\chi_\sigma\chi_\phi.
\]

So the effective source susceptibility is
\[
\chi_\sigma^{\rm eff}=
\frac{1}{\Theta_\sigma}=
\frac{\rho_*}{m c_{s,*}^2 N_{\sigma\sigma}},
\]
and the microscopic gain becomes the explicit parent quantity
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 O_{\sigma\phi}^2}{m c_{s,*}^2 K_X N_{\sigma\sigma}}.
\]
Introducing the coherence factor
\[
C_{\sigma\phi}^2=
\frac{O_{\sigma\phi}^2}{N_{\sigma\sigma}N_{\phi\phi}},

aud via Cauchy–Schwarz, one gets the exact factorization
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
C_{\sigma\phi}^2,
\qquad 0\le C_{\sigma\phi}^2\le 1.
\]

## Stage 46 — parent-overlap threshold theorem

Combining the parent gain with the Stage-44 phase diagram gives exact parent thresholds
\[
g_{\phi,\rm fail}^2=
\frac{m c_{s,*}^2 K_X N_{\sigma\sigma} G_{\rm fail}}{\rho_* O_{\sigma\phi}^2},
\qquad
g_{\phi,\rm suff}^2=
\frac{m c_{s,*}^2 K_X N_{\sigma\sigma} G_{\rm suff}}{\rho_* O_{\sigma\phi}^2}.
\]
Equivalently, in coherence form,
\[
C_{\rm fail}^2=
\frac{m c_{s,*}^2 K_X G_{\rm fail}}{\rho_* g_\phi^2 N_{\phi\phi}},
\qquad
C_{\rm suff}^2=
\frac{m c_{s,*}^2 K_X G_{\rm suff}}{\rho_* g_\phi^2 N_{\phi\phi}}.
\]

There is an exact Cauchy no-go theorem:
if
\[
G_{\max}(g_\phi):=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
< G_{\rm fail}(\kappa,\eta),
\]
then no profile engineering of \(\chi_\sigma\) can rescue the branch.

Inserting
\[
G_{\rm fail}=
\frac{Pe_{\rm req}}{\kappa\Delta_\infty},
\qquad
G_{\rm suff}=
\frac{Pe_{\rm req}}{\kappa\Delta_0},
\qquad
\kappa=K_XL^2/T_X,
\]
one finds exact amplitude thresholds
\[
g_{\phi,\rm fail}^2=
\frac{m c_{s,*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2 \Delta_\infty},
\qquad
g_{\phi,\rm suff}^2=
\frac{m c_{s,*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2 \Delta_0}.
\]
So \(K_X\) cancels from the explicit prefactor and survives only through the geometry-shape functions \(\Delta_0,\Delta_\infty\).

## Stage 47 — parent equilibrium source/support alignment

The parent equilibrium law
\[
H(y)\,\delta\rho(s,y)+\delta V_{\rm conf}(s,y)=0,
\qquad H(y):=h'(\rho_*(y)),
\]
forces the aligned source profile
\[
\chi_\sigma(y)=g_\phi\chi_\phi(y)/H(y).
\]
So the overlap invariants become
\[
O_{\sigma\phi}=g_\phi I_1,
\qquad
N_{\sigma\sigma}=g_\phi^2 I_2,
\]
with
\[
I_1=\int d^3y\;\frac{\chi_\phi(y)^2}{H(y)},
\qquad
I_2=\int d^3y\;\frac{\chi_\phi(y)^2}{H(y)^2}.
\]
Therefore
\[
C_{\sigma\phi}^2=
\frac{I_1^2}{N_{\phi\phi} I_2}\le 1.
\]
In the thin active layer where \(H(y)\approx H_w\) is nearly constant,
\[
I_1=N_{\phi\phi}/H_w,
\qquad
I_2=N_{\phi\phi}/H_w^2,
\qquad
C_{\sigma\phi}^2=1.
\]
So the matched-layer branch is not arbitrary; it is the natural thin-layer limit of the parent equilibrium branch.

The exact eliminated-source support softening is
\[
\Delta K_X^{\rm (eq)}=g_\phi^2 I_1,
\qquad
G_{\rm eq}=\Delta K_X^{\rm(eq)}/K_X = g_\phi^2 I_1/K_X.
\]

## Stage 48 — explicit thin-wall confinement branch

For the explicit wall family
\[
V_{\rm conf}(r;a)=V_0 f\bigl((r-a)/\ell\bigr),
\]
with wall coordinate \(\xi=(r-a)/\ell\), a support displacement \(a\to a+\phi(s)\) gives
\[
\delta V_{\rm conf} = +\frac{V_0}{\ell} f'(\xi)\phi(s),
\]
so the parent loading amplitude is exactly
\[
g_\phi=V_0/\ell.
\]

The shell integral entering the equilibrium gain is
\[
I_1 = 4\pi\ell\bigl[a^2J_1+2a\ell J_2+\ell^2J_3\bigr],
\]
where
\[
J_n := \int d\xi\;\frac{\xi^n f'(\xi)^2}{H(\xi)}.
\]
For a centered symmetric wall layer, \(J_2=0\), so
\[
I_1 = 4\pi\ell\bigl[a^2J_1+\ell^2J_3\bigr].
\]
The exact equilibrium gain is
\[
G_{\rm eq}=4\pi V_0^2\Bigl[\frac{a^2J_1}{\ell}+2aJ_2+\ell J_3\Bigr]/K_X.
\]
In the thin-wall limit \(\ell\ll a\), the leading gain is
\[
G_{\rm eq}^{\rm(tw)}=\frac{4\pi a^2 V_0^2 J_1}{K_X\ell}.
\]

Comparing with the Stage-44 thresholds gives wall-amplitude surfaces
\[
V_{0,\rm fail}^2=\frac{K_X\ell G_{\rm fail}}{4\pi a^2J_1},
\qquad
V_{0,\rm suff}^2=\frac{K_X\ell G_{\rm suff}}{4\pi a^2J_1}.
\]
After inserting
\[
G_{\rm fail}=
\frac{Pe_{\rm req}}{\kappa\Delta_\infty},
\qquad
G_{\rm suff}=
\frac{Pe_{\rm req}}{\kappa\Delta_0},
\qquad
\kappa=K_XL^2/T_X,
\]
the explicit \(K_X\) factor cancels:
\[
V_{0,\rm fail}^2=
\frac{T_X\ell Pe_{\rm req}}{4\pi a^2L^2J_1\Delta_\infty},
\qquad
V_{0,\rm suff}^2=
\frac{T_X\ell Pe_{\rm req}}{4\pi a^2L^2J_1\Delta_0}.
\]

For an almost constant-compressibility active wall layer, \(H(\xi)\approx H_w\),
\[
J_1 = I_f/H_w,
qquad I_f:=\int d\xi\;f'(\xi)^2,
\]
so the thresholds reduce to
\[
V_{0,\rm fail}^2=
\frac{H_w T_X\ell Pe_{\rm req}}{4\pi a^2L^2I_f\Delta_\infty},
\qquad
V_{0,\rm suff}^2=
\frac{H_w T_X\ell Pe_{\rm req}}{4\pi a^2L^2I_f\Delta_0}.
\]

## Stage 49 — dimensionless wall figure of merit

For the same thin-wall matched branch, define
\[
W_{\rm wall}:=
\frac{4\pi a^2L^2J_1V_0^2}{T_X\ell}.
\]
Since \(\kappa=K_XL^2/T_X\), this is exactly
\[
W_{\rm wall}=\kappa G_{\rm eq}^{\rm(tw)}.
\]
The Stage-44 operator theorem therefore becomes
\[
W_{\rm fail}=\frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)},
\qquad
W_{\rm suff}=\frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\]
with exact wall-form theorem:

- \(W_{\rm wall}\le W_{\rm fail}\): fail,
- \(W_{\rm wall}\ge W_{\rm suff}\): succeed,
- only the narrow intermediate band still needs the full root solve.

If \(H(\xi)\approx H_w\), the wall control variable becomes
\[
W_H=
\frac{4\pi a^2L^2I_fV_0^2}{H_w T_X\ell}.
\]
So the explicit parent branch is now controlled by one dimensionless figure of merit rather than a diffuse set of amplitudes.

## Stage 50 — sech–Gaussian coherence resonance benchmark

For the explicit independent-profile benchmark
\[
\chi_\sigma(y)=\operatorname{sech}(y/w_f),
\qquad
\chi_\phi(y)=e^{-y^2/w_g^2},
\qquad
r:=w_g/w_f,
\]
the exact norms are
\[
N_{\sigma\sigma}=2w_f,
\qquad
N_{\phi\phi}=w_g\sqrt{\pi/2}.
\]
The overlap is
\[
O_{\sigma\phi}=w_f I(r),
\qquad
I(r):=\int_{-\infty}^{\infty} dx\;\operatorname{sech}(x)e^{-x^2/r^2}.
\]
So the coherence is
\[
C^2(r)=\frac{I(r)^2}{r\sqrt{2\pi}}.
\]

Parseval/Fourier-transform arguments give the exact duality
\[
I(r)=\frac{r}{\sqrt\pi}I(\pi/r),
\qquad
C^2(r)=C^2(\pi/r).
\]
Hence the self-dual stationary point is
\[
r_* = \sqrt\pi.
\]
Numerically,
\[
C_{\rm res}^2:=C^2(\sqrt\pi)=0.9944188364515293487\ldots,
\]
so the resonance penalty factor is
\[
P_{\rm res}:=1/C_{\rm res}^2=1.0056124877605762169\ldots.
\]
Thus the best independent sech–Gaussian mismatch branch misses the ideal matched-layer coherence only by about \(0.56\%\).

## Stage 51 — resonance-corrected thresholds

The general parent gain is
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
C_{\sigma\phi}^2.
\]
On the Stage-47 matched equilibrium branch, \(C_{\sigma\phi}^2=1\), so the matched gain is
\[
G_{\rm match}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}.
\]
Stage 49 repackaged this as the wall figure of merit
\[
W_{\rm wall}=
\frac{4\pi a^2L^2J_1V_0^2}{T_X\ell}.
\]
For the independent sech–Gaussian family,
\[
G_{\rm res}(r)=C^2(r)G_{\rm match},
\qquad
W_{\rm res}(r)=C^2(r)W_{\rm wall}.
\]
Therefore the exact profile-family thresholds are
\[
W_{\rm wall}\le \frac{Pe_{\rm req}}{C^2(r)\Delta_\infty}
\quad\Rightarrow\quad \text{fail},
\]
\[
W_{\rm wall}\ge \frac{Pe_{\rm req}}{C^2(r)\Delta_0}
\quad\Rightarrow\quad \text{succeed}.
\]
At resonance \(r=\sqrt\pi\), this becomes
\[
W_{\rm wall}\le \frac{P_{\rm res} Pe_{\rm req}}{\Delta_\infty}
\quad\Rightarrow\quad \text{fail on the resonance family},
\]
\[
W_{\rm wall}\ge \frac{P_{\rm res} Pe_{\rm req}}{\Delta_0}
\quad\Rightarrow\quad \text{succeed on the resonance family}.
\]
So the explicit independent-profile benchmark modifies the Stage-49 wall thresholds by only the tiny factor \(P_{\rm res}\approx1.00561249\).
# 5PN stages 52–72 notes

This bundle continues the moving-throat/support-source chain from the resonance-threshold point
through the explicit Family-1 branch closure and the first direct comparison with the minimal
isotropic passive/outgoing quadrupole module.

## Stage 52 — final reduced verdict for the support/source program

The matched-branch theorem remains

- `W_wall <= Pe_req / Delta_inf`  -> fail,
- `W_wall >= Pe_req / Delta_0`    -> succeed,

and the explicit sech–Gaussian resonance family only perturbs those thresholds by the tiny factor
`P_res = 1.005612487760576...`, so the genuinely profile-sensitive sub-bands are only about `0.56%`
wide.

## Stage 53 — explicit GNLS wall-shell reduction

Projecting the parent GNLS quadratic shell energy onto the wall-support mode `q(s) chi_phi(y)` gives

- `T_X = hbar^2 N_(phi phi) / (4 m rho_w)`,
- `K_X = H_w N_(phi phi) + hbar^2 G_(phi phi)/(4 m rho_w)`,

and on the matched thin-wall branch the support/source fixed-point coupling is exactly the wall
figure of merit:
`Xi = W_wall`.

## Stage 54 — canonical tanh-wall branch

For the canonical wall
`f(xi) = (1 + tanh xi)/2`, `chi_phi = f'(xi) = (1/2) sech^2 xi`,
the exact moments are

- `I_f = 1/3`,
- `I_g = 4/15`,
- `I_g / I_f = 4/5`.

With the natural local Robin closure `K_m = T_X / ell`, the explicit branch is

- `eta = L/ell`,
- `kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L/ell)^2`,
- `W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2)`.

## Stage 55 — explicit branch thresholds

Writing
`chi_s = m c_(s,w) L / hbar`,
`Lambda_ell = L/ell`,
`Upsilon_w = 4 rho_w^2 V0^2 /(hbar^2 c_(s,w)^2)`,
the branch map is

- `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`,
- `eta = Lambda_ell`,
- `W_wall = Upsilon_w Lambda_ell^2`.

So the explicit thresholds are

- `Upsilon_fail = Pe_req / [Lambda_ell^2 Delta_inf(kappa,eta)]`,
- `Upsilon_suff = Pe_req / [Lambda_ell^2 Delta_0(kappa,eta)]`.

The shell-gradient and compression-dominated asymptotics from the notes are verified directly.

## Stages 56–57 — Family-1 geometry map and healing lock

Using the carried reference values

- `L/a = 37/20`,
- `ell/a = 1/20`,

the Family-1 branch fixes

- `Lambda_ell = 37`,
- `eta = 37`.

With the healing-width closure `ell = hbar/(2 m c_(s,w))`, the same branch fixes

- `chi_s = 37/2`,
- `kappa = 12321/5`,
- `alpha = sqrt(kappa) = 111/sqrt(5) ≈ 49.6407091`.

So the only remaining explicit branch input is the wall-depth amplitude.

## Stage 58 — explicit Family-1 threshold window

At `(kappa,eta) = (12321/5, 37)`, the exact support/source scales are

- `Delta_0 ≈ 1.73302079021525e-4`,
- `Delta_inf ≈ 2.01447565540522e-2`.

Hence

- `Upsilon_fail ≈ 0.0362605617972939 Pe_req`,
- `Upsilon_suff ≈ 4.21495341569977 Pe_req`,
- `Theta_fail ≈ 3.62605617972939e-4 Pe_req`,
- `Theta_suff ≈ 4.21495341569977e-2 Pe_req`

after `Upsilon_w = 100 Theta_w`.

## Stages 59–60 — exact `n=5` wall-depth lock and shell-weighted extraction

For the frozen `n=5` EOS,
`h(rho) = m c_s(rho)^2 / 4`, so with
`mu_* = lambda_mu h_w`
and the healing-width lock one gets

- `Theta_w = lambda_mu^2 rho_w^2 /(16 ell^2)`.

On the normalized Family-1 wall this becomes
`Theta_w = 25 lambda_mu^2 rho_w^2`.

Using the explicit Family-1 wall profile and canonical support weight gives

- `<rho_r>_chi ≈ 0.192619005556493`,
- `<rho_r^2>_chi ≈ 0.162745294003265`,
- `Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2`,
- `Theta_w^(J)   ≈ 0.927552032539308 lambda_mu^2`.

## Stage 61 — explicit Family-1 wall-depth verdict

Comparing the extracted `Theta_w` values to the Stage-58 window gives

- `Pe_suff^(chi) ≈ 96.5285247264386 lambda_mu^2`,
- `Pe_fail^(chi) ≈ 11220.5441626259 lambda_mu^2`,
- `Pe_suff^(J)   ≈ 22.0062226330754 lambda_mu^2`,
- `Pe_fail^(J)   ≈ 2558.01892349205 lambda_mu^2`.

So the explicit Family-1 wall-depth supply is not naturally starved; the remaining question becomes the
quadrupole-side demand `Pe_req`.

## Stages 62–64 — Family-1 demand map in `Pe`, `zeta_req`, and `Pi_tr/C_mix`

For the Family-1 branch,

- `y_F1 tan y_F1 = 37`,
- `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2) ≈ 1.00005192880220`.

The exact constructive support map is
`zeta_F1(Pe) = A_F1 Omega_Pe^2`,
with hard ceiling
`zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601`.

Converting the Stage-61 `Pe_req` windows through this map gives, at `lambda_mu = 1`,

- `zeta_suff^(chi) ≈ 2.46622291347846`,
- `zeta_fail^(chi) ≈ 2.46752913273870`,
- `zeta_suff^(J)   ≈ 2.44257571477179`,
- `zeta_fail^(J)   ≈ 2.46752736855058`.

Using
`Pi_tr = C_mix Q(zeta;eps_blk)`,
`Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta]/[1 - eps_blk zeta]`,
the unblocked (`eps_blk = 0`) explicit product window becomes

- `Pi_tr/C_mix <= 3.46622291347846`  -> guaranteed success,
- `Pi_tr/C_mix >= 3.46752913273870`  -> guaranteed failure,
- `Pi_tr/C_mix <  3.46752922945601`  -> hard explicit ceiling.

## Stages 65–70 — master reduced quadrupole residual and pure loading-ratio collapse

The whole reduced moving-throat PDE is compressed to one scalar residual
`R_quad = zeta_req - zeta_phys(Pe_*)`,
with `Pe_*` selected by the support/source fixed-point law.

After the exact demand-cancellation step, the normalized support theorem depends only on the loading ratio

`rho_alpha = alpha_req / alpha_mix`,

with support demand
`zeta_req = (rho_alpha - 1) / [1 - eps_blk (2 - rho_alpha)]`.

So the explicit Family-1 theorem is finally

- `rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,
- `rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,
- `rho_alpha < rho_max^(F1)(eps_blk)`               -> hard ceiling,

and in the natural unblocked case

- `rho_alpha <= 3.46622291347846`  -> guaranteed success,
- `rho_alpha >= 3.46752913273870`  -> guaranteed failure,
- `rho_alpha <  3.46752922945601`  -> hard ceiling.

So by Stage 70 the explicit support/source side is completely reduced to one number:
the outgoing-branch loading ratio `rho_alpha`.

## Stages 71–72 — loading ratio from the minimal isotropic conservative quadrupole module

The natural contact-plus-pole reading of the minimal isotropic conservative precursor

`Y_Q^cons(omega) = c0 + c1/(1 - omega^2/Omega_Q^2)`

gives

- `c0 = 1/rho_alpha`,
- `c1 = (rho_alpha - 1)/rho_alpha`,
- `zeta_req = c1/c0`.

For the minimal isotropic module
`c0 = 3/4`, `c1 = 1/4`,
the exact loading data are

- `rho_alpha = 4/3`,
- `zeta_req = 1/3`,
- `Pi_tr = (4/3) C_mix`.

This lies in the symmetric-lowest-twin regime and, on the explicit Family-1 branch, it already satisfies

- `zeta_req = 1/3 < A_F1 ≈ 1.00005192880220`.

So the explicit Family-1 branch succeeds at **zero transport bias** on the natural minimal isotropic passive/outgoing quadrupole module.

## Bottom line from Stage 72

By the end of this bundle, the explicit Family-1 support/source side is no longer the active reduced bottleneck. On the natural minimal isotropic contact-plus-pole branch it is already comfortably satisfied. The remaining reduced theorem question becomes whether the actual grouped-`P2` / geometry moving-throat branch really realizes that minimal isotropic passive/outgoing quadrupole module, which is the right bridge back into the 5PN grouped-`P2` program.

# 5PN / Moving-Throat Continuation — Stages 73–90

This batch continues directly from the Stage 72 result that the explicit Family-1 support/source side is no longer the active bottleneck. The live question is whether the actual grouped-\(P_2\)/geometry branch realizes the minimal isotropic contact-plus-pole conservative quadrupole module, and then whether the passive/outgoing \(l=2\) branch carries the canonical outgoing normalization.

## Stage 73
Recasts the post-72 status in theorem language: the explicit Family-1 support/source branch already succeeds on the minimal isotropic branch with
\[
\rho_\alpha=\frac43,\qquad \zeta_{\rm req}=\frac13,\qquad Pe_{\rm req}=0,
\]
so the remaining reduced theorem gap is no longer on the support/source side.

## Stage 74
Derives the `3/4 + 1/4` conservative module directly from the 3PN conservative split.
If the isotropic grouped-\(P_2\) branch is carried by one effective pole and the geometry lane is static through \(O(\omega^4)\), then the minimal isotropic branch identity forces
\[
K_{\rm pole}=\frac{K_0}{4},\qquad K_{\rm geom}=\frac{3K_0}{4},
\]
hence
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
\]

## Stage 75
Allows the geometry lane to carry dynamic even moments and derives the exact obstruction formula
\[
c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\qquad
\epsilon_2=\frac{\Omega_Q^2 K_{(g,2)}}{K_{\rm pole}},
\qquad
\epsilon_4=\frac{\Omega_Q^4 K_{(g,4)}}{K_{\rm pole}}.
\]
So the `3/4 + 1/4` split is exact iff the geometry lane is static at the relevant orders.

## Stage 76
Freezes the updated reduced status: the only remaining reduced ambiguity is whether the actual moving-throat geometry lane is dynamically inert through \(O(\omega^4)\).

## Stage 77
Proves the isotropic geometry-decoupling theorem. In the exact isotropic quadratic wall theory, the \(l=0\) scalar/geometry lane is orthogonal to the grouped real \(l=2\) bundle, so no dynamic \(O(\omega^2)\) or \(O(\omega^4)\) geometry contamination can enter the grouped quadrupole carrier at linear order:
\[
K_{(g,2)}=K_{(g,4)}=0,\qquad \epsilon_2=\epsilon_4=0.
\]

## Stage 78
Shows that if weak anisotropy later induces an \(l=0\leftrightarrow l=2\) mixing source, the first nonzero geometry contamination appears only at second order in that mixing:
\[
\epsilon_2,\epsilon_4 = O(\chi^2).
\]
So the Stage-74 split is perturbatively stable.

## Stage 79
Records the actual reduced-branch verdict: on the natural isotropic branch,
\[
\epsilon_2=\epsilon_4=0,
\]
hence
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
\qquad
\rho_\alpha=\frac43,\qquad
\zeta_{\rm req}=\frac13.
\]

## Stage 80
Shows that once the actual isotropic grouped-\(P_2\) one-pole branch is accepted, the full reduced 2.5PN normalization problem collapses to one scalar defect
\[
N_Q:=\frac{\overline K_0}{\overline K_0^{\rm target}}.
\]
All low-frequency coefficients scale by the same factor on that branch:
\[
\overline K_2,\ \overline K_4,\ \overline\Gamma_5 \propto N_Q.
\]

## Stage 81
Proves the explicit Family-1 support/source theorem is automatic on the actual isotropic branch. Since
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk})=\frac{1}{3-2\epsilon_{\rm blk}},
\]
any explicit family with \(\zeta_{\max}>1\) already passes throughout the admissible blocked regime. Family-1 has
\[
\zeta_{\max}^{(F1)}\approx 2.46752922945601 > 1,
\]
so it is automatically safe.

## Stage 82
Freezes the reduced finish line: the grouped-\(P_2\) conservative branch is geometry-clean, the support/source side is automatic, and the only remaining reduced theorem gap is the single normalization defect \(N_Q-1\).

## Stage 83
Factors the last 2.5PN obstruction into a conservative and an outgoing piece. Introducing one outgoing-normalization factor \(\chi_Q\),
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q \sigma_Q^{\rm can}\omega^5}+O(\omega^6),
\]
one gets
\[
\frac{\overline\Gamma_5}{\overline\Gamma_5^{\rm target}} = \chi_Q N_Q,
\]
and the full odd normalization condition is
\[
\hat m_0^{\,2}\chi_Q N_Q = 1.
\]

## Stage 84
Uses the natural point-particle source map \(\hat m_0\to 1\) to show the remaining reduced obstruction is purely outgoing:
\[
N_Q=\frac{1}{\chi_Q}.
\]
So all remaining reduced uncertainty is now in the outgoing normalization factor \(\chi_Q\).

## Stage 85
Shows higher odd retarded data beginning at \(O(\omega^7)\) are irrelevant to the 2.5PN theorem. The only live retarded obstruction at 2.5PN is the leading \(\omega^5\) normalization factor \(\chi_Q\).

## Stage 86
States the conditional reduced closure:
if \(\chi_Q=1\), the reduced GR-like point-particle 2.5PN theorem is closed; if not, the whole remaining failure is measured by \(\Delta_Q=\chi_Q-1\).

## Stage 87
Computes the exact outgoing spherical \(l=2\) DtN fingerprint:
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}-\frac{2z^6}{27}-i\frac{z^7}{27}+O(z^8),
\]
and therefore
\[
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}-\frac{11z^6}{729}-i\frac{z^7}{243}+O(z^8).
\]

## Stage 88
Matches the Stage-87 DtN fingerprint to the retarded grouped-\(P_2\) one-pole-plus-contact module and fixes
\[
\chi_Q=1
\]
on the canonical compact passive/outgoing branch. A deformed DtN branch with
\[
\Lambda_2^{\rm def}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\xi_Q\frac{z^5}{9}+O(z^6)
\]
simply gives \(\chi_Q=\xi_Q\).

## Stage 89
Inserts \(\chi_Q=1\) into the reduced normalization stack and closes the reduced 2.5PN theorem on the canonical outgoing branch. In the strict point-particle limit,
\[
N_Q=1,
\]
and the canonical invariant coefficients are exactly
\[
\overline K_0=\frac{54Gc_s^5}{5a^5c^5},\qquad
\overline K_2=\frac{6Gc_s^3}{5a^3c^5},\qquad
\overline K_4=\frac{8Gc_s}{15ac^5},\qquad
\overline\Gamma_5=\frac{2G}{5c^5}.
\]

## Stage 90
Builds the first explicit isotropic DtN deformation algebra. For
\[
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6),
\]
the normalized outgoing factor is
\[
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0},
\]
with \(\Sigma_2,\Sigma_4\) fixed by the requirement that the canonical even moments remain unchanged. So the only isotropic branch data that can shift \(\chi_Q\) are \(\beta,\Sigma_0,\Sigma_5\).

## Net result of Stages 73–90

Inside the present reduced hierarchy, the conservative grouped-\(P_2\)/geometry branch now reaches the same minimal isotropic `3/4 + 1/4` module that the 2.5PN audit wanted, the explicit Family-1 support/source side is automatic, and the outgoing \(l=2\) DtN branch fixes the last reduced scalar to
\[
\chi_Q=1
\]
on the canonical compact branch. So the reduced nonspinning point-particle 2.5PN theorem is closed on that canonical branch; what remains genuinely PDE-facing is branch realization and explicit isotropic DtN deformation data.
# 5PN / Moving-Throat Continuation — Stages 91–100

This batch continues directly from Stage 90, where the isotropic outgoing normalization factor was reduced to
\[
\chi_Q=\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
\]
The goal of Stages 91–100 is to stop treating \((\beta,\Sigma_0,\Sigma_5)\) as abstract branch labels and push them through explicit outlet/core models.

## Stage 91
Classifies the robustness classes of \(\chi_Q\).

- Pure overall scaling is exactly invisible:
  \[
  \Lambda_2^{\rm def}=S\Lambda_2^{\rm out}
  \quad\Longrightarrow\quad
  \chi_Q=1.
  \]
- Pure scale+argument deformation preserves the already-fixed even moments only on the natural positive branch \(\beta=1\), so it is also harmless.
- A genuine additive isotropic throat-core channel can move \(\chi_Q\) while leaving the lower even moments canonical:
  \[
  \chi_Q=\frac{3(S+9\Sigma_5)}{3S-\Sigma_0}
  \qquad(\beta=1).
  \]
- The exact preservation submanifold is
  \[
  \Sigma_5=\frac{S(1-\beta^5)}{9}-\frac{\Sigma_0}{27}.
  \]

## Stage 92
Linearizes around the canonical outgoing branch:
\[
S=1+\varepsilon s,\qquad
\beta=1+\varepsilon b,\qquad
\Sigma_0=\varepsilon a_0,\qquad
\Sigma_5=\varepsilon a_5.
\]
Then
\[
\chi_Q
=
1+\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)+O(\varepsilon^2).
\]
So the minimal isotropic branch-selection data are the triple
\[
(b,a_0,a_5),
\]
and first-order preservation requires
\[
5b+\frac{a_0}{3}+9a_5=0.
\]

## Stage 93
Introduces the first explicit isotropic geometric outlet:
\[
\Lambda_2^{\rm R}(z)=\Lambda_2^{\rm out}(z)+\rho_R.
\]
The normalized response is
\[
\widehat Y_2^{\rm R}(z)
=
1+\frac{z^2}{9-3\rho_R}
+\frac{4-\rho_R}{9(3-\rho_R)^2}z^4
+i\frac{z^5}{27-9\rho_R}+O(z^6),
\]
so
\[
\chi_Q^{\rm R}=\frac{3}{3-\rho_R}.
\]
A pure Robin core therefore generically shifts both the even branch and the odd normalization.

## Stage 94
Adds the first explicit isotropic hidden mixed side-channel:
\[
\Lambda_2^{\rm mix}(z)
=
\Lambda_2^{\rm out}(z)
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}+O(z^6).
\]
The even-preserving conditions force
\[
\kappa_W=-\frac19,
\qquad
\sigma_W=0.
\]
So a standalone isotropic hidden pole cannot sit on the canonical even branch unless it is absent. Its normalization factor is
\[
\chi_Q^{\rm mix}
=
\frac{3(1-9\sigma_W\gamma_W)}{3+\sigma_W}.
\]

## Stage 95
Combines the Robin core and the mixed side-channel:
\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}+O(z^6).
\]
Imposing the canonical even branch yields exactly two solutions:
\[
(\rho_R,\kappa_W)=(\sigma_W,0)
\quad\text{or}\quad
(\rho_R,\kappa_W)=\left(4\sigma_W,\frac13\right).
\]
The second is the nontrivial compensated branch. On it,
\[
\chi_Q^{\rm hyb}=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
\]
and canonical odd normalization is preserved iff
\[
\gamma_W=\frac19.
\]
With that value the whole outlet collapses to a pure harmless scale deformation.

## Stage 96
Freezes the outlet-model classification:

1. pure Robin loading generically shifts \(\chi_Q\);
2. a standalone isotropic mixed pole is too rigid and cannot preserve the canonical even branch unless it vanishes;
3. the hybrid Robin–mixed outlet admits one nontrivial compensated canonical-even branch, with
   \[
   \rho_R=4\sigma_W,\qquad \kappa_W=\frac13,\qquad \gamma_W=\frac19
   \]
   on the fully preserved canonical branch.

## Stage 97
Replaces the reduced outlet coefficients by a concrete two-channel core model with internal variables

- \(s(\omega)\): static shell/compliance mode,
- \(q(\omega)\): mixed \(A_w/F_{\mu w}\)-type side-channel mode.

The linear core system
\[
\begin{pmatrix}
K_s & \lambda\\
\lambda & -K_q D_W^{\rm bare}(z)
\end{pmatrix}
\binom{s}{q}
=
u\binom{g_s}{g_q}
\]
gives the exact Schur-complement outlet
\[
\delta\Lambda_{\rm core}(z)
=
\frac{g_s^2}{K_s}
-
\frac{(K_s g_q-\lambda g_s)^2}
{K_s\big(K_sK_q D_W^{\rm bare}(z)+\lambda^2\big)}.
\]
Defining
\[
r_c=\frac{\lambda^2}{K_sK_q},
\]
the reduced Robin–mixed parameters are
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=\frac{(K_s g_q-\lambda g_s)^2}{K_s^2K_q(1+r_c)},
\qquad
\kappa_c=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_c=\frac{\gamma_0}{1+r_c}.
\]

## Stage 98
Determines exactly when the concrete core lands on the compensated canonical branch. The nontrivial compensation condition is
\[
\rho_c=4\sigma_c,\qquad \kappa_c=\frac13,\qquad \gamma_c=\frac19.
\]
This is equivalent to the exact coupling-balance law
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
\]
so
\[
g_q=
\frac{g_s}{2K_s}\left(2\lambda\pm\sqrt{K_sK_q+\lambda^2}\right).
\]
The bare mixed channel must itself be a scale-deformed copy of the canonical compact outgoing branch:
\[
\kappa_0=\frac{1+r_c}{3},\qquad
\gamma_0=\frac{1+r_c}{9}.
\]
On that surface,
\[
\delta\Lambda_{\rm core}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\qquad
\sigma_*=\frac{g_s^2}{4K_s},
\]
and the normalized outgoing fingerprint stays exactly canonical.

## Stage 99
Realizes the bare mixed side-channel geometrically as the first D/N half-wave on an auxiliary tube of length \(L_W\):
\[
k_W=\frac{\pi}{2L_W},
\qquad
\Omega_W=\frac{\pi c_s}{2L_W}.
\]
In outlet variables \(z=a\omega/c_s\),
\[
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
\]
The compensation condition \(\kappa_0=(1+r_c)/3\) fixes
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
\]
If the bare mixed outlet is a pure-scale deformation of the canonical compact outgoing \(l=2\) branch,
\[
D_W^{\rm bare}(z)=(1+r_c)\left(1-\frac{z^2}{3}-i\frac{z^5}{9}\right)+O(z^6),
\]
then the hybridization factor is removed exactly and
\[
\kappa_c=\frac13,\qquad \gamma_c=\frac19.
\]

## Stage 100
Freezes the concrete outlet-core status. The surviving PDE-facing question is no longer “is there some deformed outlet?” It is much sharper:

> Does the actual moving-throat core realize the concrete coupling-balance surface together with the auxiliary D/N-tube normalization?

On that surface the effective outlet preserves the canonical normalized outgoing fingerprint exactly.

## Net result of Stages 91–100

The isotropic outgoing-branch ambiguity is no longer open-ended. The first explicit outlet audit shows:

- pure scale deformations are harmless,
- pure Robin loading and a standalone isotropic mixed pole do not by themselves preserve the canonical branch,
- a specific compensated Robin–mixed outlet does preserve it,
- that compensated branch is realized by a concrete two-channel throat-core model,
- and the mixed side-channel can be given a direct finite D/N tube realization.

So the next honest step is Stage 101: extract the concrete core parameters \((K_s,K_q,\lambda,g_s,g_q)\) from a parent-action throat-core ansatz rather than leaving them as reduced variables.

# 5PN stages 101–120 notes

This bundle continues the chain immediately after the compensated outlet/core result at Stage 100.
The focus of this block is the next honest microscopic gate: replace the reduced outlet/core
coefficients by explicit parent-action overlaps, then carry that data all the way into the
mouth-layer fixed-point problem.

## Stage 101 — parent-action extraction of the core parameters

The reduced two-channel core variables from Stages 97–100 are replaced by explicit overlap
formulas from one concrete GNLS + localized-Maxwell throat-core ansatz.

The shell/compliance mode gives
\[
K_s
=
4\pi a^2\left(
\frac{H_w\ell}{3}
+
\frac{\hbar^2}{15m_\psi\rho_w\,\ell}
\right),
\]
and on the healing-locked shell branch
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\qquad
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
\]

The mixed D/N half-wave gives
\[
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi}{\sqrt2\,L_W^{3/2}}.
\]

The shell–mixed hybridization and shell mouth coupling become
\[
\lambda=-q_* v_{w0}\mathcal I_{sq},
\qquad
\mathcal I_{sq}=\frac{8\sqrt2}{3}a^2\ell\sqrt{L_W},
\qquad
g_s=\mathcal T_m \frac{4\pi a^2\ell}{3}.
\]

So the Stage-97 core matrix is no longer abstract; every entry now has an explicit parent
overlap meaning. This is exactly the next gate identified in the moving-throat notes. fileciteturn30file16

## Stages 102–103 — collapse to normalized parent ratios

The exact core-balance condition
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2
\]
collapses to the two dimensionless ratios
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\]
through
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.
\]

The D/N mixed-tube length becomes
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
\]

So the surviving outlet/core theorem gate is no longer “what are the reduced coefficients?”
It is only: which branch point \((\mathfrak r,\mathfrak g)\) the actual GNLS + localized-Maxwell
core selects. fileciteturn30file16turn30file5

## Stages 104–107 — explicit Family-1 bridge

To keep the executable chain sequential, I added the missing Family-1 bridge audits that the
later notes assume implicitly.

Using the carried Family-1 geometry \(L/a=37/20\) together with the D/N length law gives
\[
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
\approx 1.77799353547498.
\]

The two compensated coupling branches are
\[
\mathfrak g_\pm^{F1}
=
\mathfrak r_{F1}\pm\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
numerically
\[
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
\]

These bridge scripts also verify the useful ordering
\[
\frac{2}{\pi}<\mathfrak g_-^{F1}<\frac{\pi}{4}<1<\mathfrak g_+^{F1},
\]
which is exactly the window the later positive-source theorem needs.

## Stages 108–111 — positive mouth-source selection

For any positive normalized axial mouth source profile \(\sigma(z)\) on the first D/N interval,
the mouth-bias factor is
\[
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz,
\]
so positivity immediately forces
\[
0\le \mathfrak g[\sigma]\le 1.
\]

Since \(\mathfrak g_+^{F1}>1\), the upper compensated branch is impossible under any positive
localized mouth source, while \(\mathfrak g_-^{F1}\in(0,1)\), so the lower branch is the unique
physically admissible canonical candidate.

The first explicit positive families then show that the lower branch is easy to reach:

- self-matched derivative source:
  \[
  \mathfrak g_{\rm match}=\frac{\pi}{4},
  \]
  only \(3.61\%\) away in traction from exact lower-branch compensation;

- convex derivative/uniform family:
  \[
  \sigma_\xi=(1-\xi)k\cos(kz)+\xi/L,
  \]
  reaches the exact lower branch at
  \[
  \xi_*\approx 0.183918405511540;
  \]

- slab and truncated-exponential penetration laws reach the same branch at
  \[
  x_*^{\rm slab}\approx 0.797839360904564,
  \qquad
  x_*^{\exp}\approx 0.662765402623161.
  \]

So by Stage 111 the branch-selection ambiguity has collapsed: the lower compensated branch is
the unique positive-source branch and is reached with moderate penetration, not by exotic
sign-changing mouth forcing.

## Stages 112–115 — explicit GNLS + localized-Maxwell mouth boundary layer

The abstract positive-source family is replaced by the first explicit boundary-layer law.
With the mouth free energy
\[
F_{\rm mouth}[\sigma]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\!\left(\ln\frac{\sigma}{\sigma_*}-1\right)
+
V_m(z)\sigma
\Big],
\qquad
V_m(z)\approx V_1 z,
\]
and zero-flux Onsager current, the exact normalized source law is
\[
\sigma_\Pi(z)=\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})},
\qquad
\Pi=\frac{V_1L}{\Theta_\sigma}.
\]

Its exact mouth-bias factor is
\[
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}
{(4\Pi^2+\pi^2)(e^\Pi-1)},
\]
with range \(2/\pi\to1\) as \(\Pi:0^+\to\infty\).

Solving \(\mathfrak g_\Pi=\mathfrak g_-^{F1}\) gives the unique canonical Family-1 point
\[
\Pi_* \approx 1.50882951349316.
\]

So the parent threshold is now concrete:
\[
\Pi_m=\frac{L}{\Theta_\sigma}(T_m-q_*A_0')=\Pi_*,
\]
with local sensitivity
\[
\mathfrak g_*' \approx 0.0714453558083196.
\]

At this point the mouth problem is no longer branch sign or family choice. It is the
single microscopic bias question: does the real mouth layer select \(\Pi_m\approx1.51\)?

## Stages 116–120 — coupled mouth fixed point and explicit core-to-mouth gain map

The next honest step is to couple the shell/compliance and mixed Maxwell channels directly in
the mouth layer. The exact scalar D/N response kernel to the exponential source is
\[
\mathcal S(\Pi,\kappa)
=
\frac{
\Pi\left[\kappa\tanh\kappa+\Pi\left(e^{-\Pi}\operatorname{sech}\kappa-1\right)\right]
}{
(1-e^{-\Pi})(\kappa^2-\Pi^2)
},
\]
with \(\mathcal S(\Pi,0)=1\).

So the fully coupled mouth bias obeys
\[
\Pi = \sum_\alpha M_\alpha \mathcal S(\Pi,\kappa_\alpha).
\]

On the first explicit Family-1 branch with one static shell lane and one mixed D/N half-wave,
\[
\kappa_s=0,
\qquad
\kappa_q=\frac{\pi}{2},
\]
this reduces to
\[
\Pi = M_s + M_q \mathcal S_q(\Pi),
\qquad
\mathcal S_q(\Pi)=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
\]

At the canonical point,
\[
\mathcal S_q(\Pi_*)\approx 0.658075937605429,
\]
so the exact gain line is
\[
M_s=\Pi_* - M_q \mathcal S_q(\Pi_*).
\]

Imposing the outlet-consistent shell:mixed ratio \(4:-1\) collapses the mouth problem to
\[
\Pi = \Sigma_m[4-\mathcal S_q(\Pi)].
\]
At \(\Pi_*\),
\[
\Sigma_m^*\approx 0.451485277739088,
\qquad
M_s^*\approx 1.80594111095635,
\qquad
M_q^*\approx -0.451485277739088.
\]

Finally, Stage 120 removes the last abstract gain pair entirely. Using the exact Stage-97
Schur complement,
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=\frac{(K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)},
\]
the actual mouth-layer gains are
\[
M_s=\frac{L g_s^2}{K_s\Theta_\sigma},
\qquad
M_q=
-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
\]

So by the end of Stage 120 the mouth fixed-point ambiguity has collapsed from a free source
profile and a free gain pair to one explicit set of parent core quantities. The next clean step
is to normalize these gains and carry them into the self-consistent branch law beyond 120. fileciteturn30file11turn30file16
# Manifest — 5PN stages 121–140

Helper module:
- `fivepn_stage121_140_common.py`

Scripts:
- `5pn_stage121_normalized_mouth_gain_family.py`
- `5pn_stage122_family1_actual_mouth_gains.py`
- `5pn_stage123_selfmatched_mouth_susceptibility.py`
- `5pn_stage124_mouth_gain_status.py`
- `5pn_stage125_selfconsistent_mouth_branch.py`
- `5pn_stage126_equal_normalized_singular_limit.py`
- `5pn_stage127_unique_regular_canonical_branch.py`
- `5pn_stage128_mouth_branch_selection_status.py`
- `5pn_stage129_positive_deformation_expansion.py`
- `5pn_stage130_first_order_rigidity_kernel.py`
- `5pn_stage131_representative_positive_families.py`
- `5pn_stage132_mouth_rigidity_status.py`
- `5pn_stage133_full_profile_mouth_potential.py`
- `5pn_stage134_first_order_selected_correction.py`
- `5pn_stage135_family1_actual_correction.py`
- `5pn_stage136_full_mouth_correction_status.py`
- `5pn_stage137_coevolving_core_mouth_map.py`
- `5pn_stage138_frozen_traction_fixedpoint.py`
- `5pn_stage139_renormalized_canonical_branch.py`
- `5pn_stage140_core_mouth_coevolution_status.py`

Each script has a paired `_output.txt` file with its run output.
# 5PN continuation — Stages 121–140

This batch covers the next explicit mouth/core chain after Stage 120.

## What is now fixed

Stages 121–123 rewrite the mouth gains in normalized parent variables and collapse the explicit core-to-mouth map to

\[
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2},
\qquad
R_q:=-\frac{M_q}{M_s}=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
\]

On the exact compensation family
\[
\mathfrak g_c=\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2},
\]
that ratio collapses to
\[
R_q=\frac14.
\]
So the Stage-118 outlet-consistent closure is derived rather than assumed. On the self-matched mouth-susceptibility closure the overall shell gain becomes
\[
\Sigma_0=M_s=\frac{20}{9}\,\widehat T_m^2.
\]

Stage 122 evaluates this on the explicit Family-1 branch:

- natural equal-normalized branch:
  \[
  M_s\approx 1.66854,
  \qquad
  M_q\approx -0.24270;
  \]
- exact compensated branch:
  \[
  M_s\approx 1.80594,
  \qquad
  M_q\approx -0.45149.
  \]

The corresponding normalized traction amplitudes differ by only about 4.04%.

## Branch selection

Stages 125–128 close the self-consistent Family-1 mouth branch law,
\[
\Pi=\Sigma_0\bigl[1-R_q(\Pi)\,\mathcal S_q(\Pi)\bigr],
\qquad
R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

The positive exponential mouth family obeys
\[
0<\mathfrak g_\Pi<1
\quad\text{for every finite }\Pi>0,
\]
so the equal-normalized branch \(\mathfrak g_c=1\) is only a singular point-source limit. The upper compensated branch is impossible because \(\mathfrak g_+^{F1}>1\). The lower compensated branch is the unique regular finite branch, reached at
\[
\Pi_*\approx 1.50882951349316,
\qquad
\widehat T_{m,*}\approx 0.901484054174204.
\]

## Mouth rigidity under positive non-exponential deformations

Stages 129–132 derive the first-order rigidity kernel around the canonical point. For any positive normalized deformation, only two source moments matter,
\[
\bar g_\varsigma=\int_0^1\varsigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx,
\qquad
\bar S_\varsigma=\int_0^1\varsigma(x)K_q(x)dx,
\]
and the traction shift is
\[
\delta\widehat T_m
=
\epsilon\Bigl[
A_T(\bar g_\varsigma-\mathfrak g_*)
+
B_T(\bar S_\varsigma-\mathcal S_*)
\Bigr].
\]
The explicit coefficients are
\[
A_T\approx -4.27263956256927,
\qquad
B_T\approx 0.134875005736706,
\]
so overlap changes dominate mixed-kernel changes by a factor
\[
|A_T|/B_T\approx 31.68.
\]

Representative families:

- uniform broadening raises the canonical point,
  \[
  \frac{\delta\Pi_u}{\epsilon}\approx 1.69941,
  \qquad
  \frac{\delta\widehat T_{m,u}}{\epsilon}\approx 0.508756;
  \]
- self-matched derivative sharpening lowers it,
  \[
  \frac{\delta\Pi_d}{\epsilon}\approx -0.382993,
  \qquad
  \frac{\delta\widehat T_{m,d}}{\epsilon}\approx -0.116944.
  \]

So the mouth-side ambiguity is no longer branch choice; it is one explicit rigidity-kernel problem.

## Full-profile mouth correction

Stages 133–136 replace the tangent exponential potential by the full coupled shell + mixed D/N mouth profile. The exact residual
\[
R_*(x)=\Phi_*(x)-\Pi_*x
\]
is tangent matched,
\[
R_*(0)=0,
\qquad
R_*'(0)=0,
\]
and has negative curvature at the mouth,
\[
R_*''(0)=-3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0,
\]
so the actual full profile broadens the source relative to the tangent exponential branch.

Projecting that residual onto the Stage-130 rigidity kernel gives the actual first-order correction:
\[
\delta\mathfrak g_{\rm act}\approx -0.06480697,
\qquad
\delta\mathcal S_{\rm act}\approx -0.03887184,
\]
\[
\delta\Pi_{\rm act}\approx 0.907084,
\qquad
\delta\widehat T_{m,\rm act}\approx 0.271654.
\]
So the corrected canonical point is approximately
\[
\Pi_{\rm corr}\approx 2.41591,
\qquad
\widehat T_{m,\rm corr}\approx 1.17314.
\]
A one-step nonlinear iterate shifts slightly further but in the same direction.

## Full core–mouth co-evolution

Stages 137–140 promote the corrected mouth layer to a fully co-evolving fixed point,
\[
\Sigma(x)=\frac{e^{-\Phi_{\Sigma_0}[\Sigma](x)}}{\int_0^1 e^{-\Phi_{\Sigma_0}[\Sigma](y)}dy},
\qquad
\Phi_{\Sigma_0}[\Sigma](x)=\Sigma_0\bigl[\mathcal T_s[\Sigma](x)-\mathcal R[\Sigma]\mathcal T_q[\Sigma](x)\bigr],
\]
with
\[
\mathcal R[\Sigma]=\frac{(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

At the old canonical traction \(\Sigma_0^*\), the fixed point stays close in bias but drifts off exact compensation,
\[
\mathfrak g_{\rm fp}\approx 0.69336,
\qquad
\mathcal R_{\rm fp}\approx 0.28271,
\qquad
\Pi_{\rm fp}\approx 1.48857.
\]

Restoring exact lower-branch compensation requires a unique renormalized traction. With the reduced numerical fixed-point solver used in this batch, that root is
\[
\Sigma_0^{\rm can}\approx 4.65077,
\qquad
\widehat T_{m,\rm can}\approx 1.44667,
\qquad
\Pi_{\rm can}\approx 3.87134.
\]
These values are very close to the handoff’s quoted ones; the small differences come from the reduced collocation/iteration resolution used in the executable scripts.

## Bottom line after Stage 140

Inside the explicit Family-1 closure:

1. the upper compensated branch is impossible;
2. the equal-normalized branch is singular;
3. the lower compensated branch remains the unique regular canonical branch;
4. full self-consistency preserves that branch only after a finite traction/bias renormalization.

So the mouth/core side is no longer blocked by branch ambiguity. The remaining 2.5PN/5PN theorem gap is the projection of the residual deviation from this renormalized canonical Family-1 point into the outgoing quadrupole-normalization defect.

# 5PN / Moving-Throat continuation — Stages 141–150

This batch codifies the next handwritten moving-throat block after Stage 140: the linear defect-transport, hybrid-outlet projection, bare mixed-port slippage, D/N similarity decomposition, parent-family rigidity, microscopic off-family normal coordinate, explicit log-channel reduction, exact lower-branch drift laws, four-observable bundle inversion, and the tangent-compensation theorem.

## What is newly frozen in executable form

### Stage 141
The co-evolving Family-1 point is reduced to the linear mouth/core defect ledger
\[
\delta\mathcal R = -\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_*^2}},
\]
\[
\delta M_s = \delta\Sigma_0,
\qquad
\delta M_q = -\frac14\,\delta\Sigma_0 + \frac{\Sigma_{0,*}}{\sqrt{1+\mathfrak r_*^2}}\,\delta\mathfrak g,
\]
\[
\delta\Pi
=
\left(1-\frac14\mathcal S_*\right)\delta\Sigma_0
-\frac{\Sigma_{0,*}}{4}\,\delta\mathcal S
+\frac{\Sigma_{0,*}\mathcal S_*}{\sqrt{1+\mathfrak r_*^2}}\delta\mathfrak g.
\]

### Stage 142
Projecting the defect to the compensated Robin–mixed outlet gives the exact first-order outlet defects
\[
\delta E_2 = \frac{\delta\mathcal C - 9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\]
\[
\delta E_4 = \frac{5\delta\mathcal C - 72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]
\[
\Delta_Q = \frac{\delta\mathcal C - 27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)}.
\]
Imposing the canonical-even gate \(\delta E_2=\delta E_4=0\) yields
\[
\delta\mathcal C = 0,
\qquad
\delta\kappa_W = 0,
\qquad
\delta\mathfrak g = 0,
\]
and therefore
\[
\Delta_Q = -\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
\]

### Stage 143
The concrete core algebra collapses the remaining odd defect to the bare mixed-port slippage scalar
\[
\delta\mathfrak B_W := \delta\gamma_0 - \frac13\,\delta\kappa_0,
\qquad
\delta\gamma_W = \frac{\delta\mathfrak B_W}{1+r_{c,*}}.
\]
With the tangential susceptibility \(\delta\mathfrak B_W = \Upsilon_\Pi\,\delta\Pi_{\rm tan}\),
\[
\Delta_Q
=
-\frac{9\sigma_*\,\Upsilon_\Pi}{(1-\sigma_*)(1+r_{c,*})}\,\delta\Pi_{\rm tan}.
\]

### Stage 144
The black-box susceptibility is decomposed into the D/N similarity-slippage scalar
\[
\Xi_{\rm slip}:=\Xi_\gamma - 2\Xi_L,
\qquad
\Upsilon_\Pi = \frac{1+r_{c,*}}{9}\,\Xi_{\rm slip},
\]
so the reduced defect becomes
\[
\Delta_Q = -\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]
If the D/N similarity law
\[
\gamma_0 = \frac{4L_W^2}{3\pi^2 a^2}
\]
is preserved to first order, then \(\Xi_{\rm slip}=0\).

### Stage 145
On the exact parent compensation family
\[
1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2,
\qquad
\frac{L_W}{a} = \frac{\pi}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
\qquad
\gamma_0 = \frac{1+\mathfrak r^2}{9},
\]
automatic similarity preservation is exact:
\[
\delta\ln\gamma_0 - 2\,\delta\ln(L_W/a) = 0.
\]
On the lower branch, \(\delta\mathfrak g=0\) implies \(\delta\mathfrak r=0\), so every first-order similarity defect vanishes and
\[
\Delta_Q = 0.
\]

### Stage 146
The exact off-family normal coordinate is
\[
\delta_\perp := \delta\mathfrak g - \mathfrak g'_-(\mathfrak r_*)\,\delta\mathfrak r,
\]
with
\[
\delta\mathcal F = 4\sqrt{1+\mathfrak r_*^2}\,\delta_\perp,
\qquad
\delta R_q = -\frac{\delta_\perp}{\sqrt{1+\mathfrak r_*^2}}.
\]
Its explicit microscopic form is
\[
\delta_\perp
=
\mathfrak g_*\,
\delta\ln\!\left(\frac{g_q K_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}\,
\delta\ln\!\left(\frac{K_s K_q}{\lambda^2}\right).
\]

### Stage 147
Those two log channels are reduced to explicit wall/BdG/Maxwell/mixed drifts. Under the carried overlap and healing-lock closures:
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
\]
\[
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
\]
So \(\delta_\perp\) becomes an explicit linear combination of
\[
(\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W).
\]

### Stage 148
The exact lower compensated branch fixes
\[
\delta\ln L_W = \delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_s^2 c_{s,w}^3}{\rho_w a^5}\right),
\]
\[
\delta\ln\mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_{s,w}^3}{\rho_w c_s^2 a^3}\right).
\]
With the frozen \(n=5\) wall-EOS branch, the irreducible microscopic drift space collapses to
\[
(\delta\ln\mathcal Z_q,\ \delta\ln\rho_w,\ \delta\ln c_s,\ \delta\ln a).
\]

### Stage 149
Those four drifts are exactly inverted by the bundle observables
\[
(\Theta_w,\ K_s,\ K_q,\ P_0),
\qquad
P_0=\frac{N_0}{D_0}.
\]
The exact inversion is
\[
\delta\ln\rho_w = \frac12\,\delta\ln\Theta_w,
\]
\[
\delta\ln a = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln c_s = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w + \frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q = \delta\ln K_q - \frac25\,\delta\ln P_0.
\]

### Stage 150
Every remaining first-order mouth/background drift is an explicit algebraic image of \((\Theta_w,K_s,K_q,P_0)\):
\[
\delta\ln c_{s,w} = \delta\ln\Theta_w,
\qquad
\delta\ln\ell = -\delta\ln\Theta_w,
\qquad
\delta\ln L_W = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{13}{8}\,\delta\ln\Theta_w,
\]
\[
\delta\ln \mathcal T_m
=
-\frac54\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{15}{8}\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_q
=
-\frac34\,\delta\ln K_s
+\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln\lambda = \frac12(\delta\ln K_s+\delta\ln K_q).
\]
The tangent-compensation theorem then holds exactly:
\[
\delta\ln r_c = 0,
\qquad
\delta\ln\mathfrak r = 0,
\qquad
\delta\ln\mathfrak g = 0,
\qquad
\delta_\perp = 0.
\]

## What this means

The remaining first-order isotropic problem is now no longer “general branch drift.” It has collapsed to a bundle-observable calculation:

\[
(\Theta_w,\ K_s,\ K_q,\ P_0)
\longrightarrow
\text{all first-order isotropic mouth/core/background drifts}.
\]

And the executable result is stronger than expected: **arbitrary first-order isotropic bundle drift is tangent to the exact compensated Family-1 parent family.**

So the next live theorem gate after Stage 150 is not first-order isotropic transport anymore. It is the first correction that escapes this closed algebra, i.e. the first **off-bundle** slippage.
# Stage 151–160 continuation notes

This batch fills the missing bridge between the Stage 150 tangent-compensation theorem and the later outgoing-load / similarity-orbit chain.

The key logical contraction is:

1. **Stage 151**: first-order off-bundle motion is not a large microscopic vector; it collapses to three scalar slippages \((\varepsilon_L,\varepsilon_v,\varepsilon_T)\), and then to one weighted scalar \(\varepsilon_\perp\) with \(\delta_\perp=-\varepsilon_\perp\).

2. **Stage 152**: pure grouped real `P2` anisotropy cannot linearly feed those scalar slippages. Its scalar feed-down begins only at quadratic order through the grouped invariant
   \[
   \mathcal I[x,y]=\frac15\,\delta x^T G_{\rm grp}\,\delta y
   =4a_x a_y+\frac45 b_x b_y.
   \]
   On the weak-axisymmetric branch this becomes
   \[
   \mathcal I[x,y]=\frac{7}{10}\epsilon^2 x^{(1)}y^{(1)}.
   \]

3. **Stage 153**: the remaining **linear** grouped problem therefore lives only in the direct outlet coefficients. It collapses to
   \[
   \mathcal K_A=\delta D_{A,2}+\frac19\,\delta D_{A,0},
   \qquad
   \mathcal G_A=\delta N_{A,0}-P_0\,\delta D_{A,0},
   \]
   together with the hidden-even compatibility relation
   \[
   \delta D_{A,4}=\frac23\,\delta D_{A,2}+\frac1{27}\,\delta D_{A,0}.
   \]

4. **Stage 154**: those two grouped obstructions are not sourced by the whole microscopic bundle independently. They decompose into:
   \[
   \mathcal K_A=\mathcal W_A-\mathcal B_A-\mathcal Z_A,
   \qquad
   \mathcal G_A=-P_0\delta K_A+P_0\delta B_{A,0}+\mathcal N_A.
   \]
   So the linear grouped-even/odd problem is driven only by wall, BdG, conservative Maxwell/mixed, and outgoing-transfer anisotropies.

5. **Stage 155**: on the canonical compensated branch, the microscopic obstruction pair is just the physical slope pair
   \[
   \mathfrak K_1=-D_0 u_2^{(1)},
   \qquad
   \mathfrak G_1=D_0 P_1.
   \]
   The direct outlet amplitudes become
   \[
   \kappa_1=-\frac{3(1-\sigma_*)}{\sigma_*}u_2^{(1)},
   \qquad
   \gamma_1=-\frac{1-\sigma_*}{9\sigma_*}\frac{P_1}{P_0}.
   \]

6. **Stage 156**: expanding the actual grouped operator moments on the weak-axisymmetric branch yields
   \[
   u_2^{(1)}=-\frac{D_{21}+D_{01}/9}{D_0},
   \qquad
   \frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
   \]
   On the even-preserving branch,
   \[
   D_{21}=-\frac{D_{01}}{9},
   \qquad
   D_{41}=-\frac{D_{01}}{27},
   \]
   so the whole remaining linear grouped normalization defect is one static loading mismatch
   \[
   \Xi_{\rm load}:=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
   \]

7. **Stage 157**: that mismatch is the weighted failure of static self-similarity relative to the wall baseline:
   \[
   \Xi_{\rm load}
   =(\delta_N-\delta_K)+\omega_B(\delta_B-\delta_K)+\omega_Z(\delta_Z-\delta_K).
   \]

8. **Stage 158**: factoring by the wall baseline yields wall-normalized shape/load variables
   \[
   B_{0,\alpha}=K\chi_\alpha^2,
   \qquad
   Z_0^{(r)}=K\Upsilon_r,
   \qquad
   N_0^{(r)}=\Lambda_r^2,
   \]
   and the outgoing-load theorem
   \[
   2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r=\delta_K
   \]
   on conservative-shape-preserving branches.

9. **Stage 159**: the outgoing load factor splits exactly into
   \[
   \frac{\Lambda_r^2}{K}
   =
   \mathcal M_r^2\frac{(1+\mathcal I_r)^2}{(1-\mathcal H_r)^2},
   \]
   with
   \[
   \mathcal M_r=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
   \quad
   \mathcal I_r=\frac{R_rG_{U,r}}{\Omega_{U,r}^2G_{W,r}},
   \quad
   \mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
   \]
   Under interference/hybridization rigidity, zero defect is the square-root mixed-leg law
   \[
   \frac{G_{W,r}}{\Omega_{W,r}^2}\propto \sqrt K.
   \]

10. **Stage 160**: on the weak-axisymmetric branch the whole outgoing-slippage bundle collapses to one scalar amplitude
    \[
    \Xi_1
    =
    \sum_r \rho_r^{(N)}
    \left[
      2\mathfrak m_r
      +\frac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r
      +\frac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r
    \right],
    \]
    with the grouped pattern
    \[
    \Xi_{\rm load}^{(20)}=\epsilon\,\Xi_1,
    \qquad
    \Xi_{\rm load}^{(21)}=\frac{\epsilon}{2}\Xi_1,
    \qquad
    \Xi_{\rm load}^{(22)}=-\epsilon\,\Xi_1,
    \]
    and the exact physical identification
    \[
    \Xi_1=\frac{P_1}{P_0}.
    \]

So after Stage 160, the remaining weak-axisymmetric grouped `2.5`PN bottleneck is no longer a diffuse outlet-bundle problem. It is just the single microscopic amplitude \(\Xi_1 = P_1/P_0\) on the actual moving-throat branch.
# 5PN / moving-throat continuation — Stage 161–170 bundle

This bundle fills the previously missing numbered continuation after Stage 160.

## What is in this bundle

The scripts are split one stage at a time so the executable numbering now matches the note chain:

- Stage 161 — outgoing-port co-loading theorem
- Stage 162 — wall-normalized transfer-shape theorem
- Stage 163 — effective transfer-shape collapse
- Stage 164 — coherent tracking-branch defect law and support-blindness
- Stage 165 — microscopic coherent-kernel slippage decomposition
- Stage 166 — exact triangular normal form
- Stage 167 — branch-invariant coordinates
- Stage 168 — direct microscopic monomials and compatibility ledger
- Stage 169 — exact microscopic similarity orbit
- Stage 170 — exact orbit–quotient closure

Each script has a paired `_output.txt` file captured from a clean run in this environment.

## Structural summary

The chain now reads:

1. `Xi_1` is the mismatch between the outgoing-weighted static transfer slope and the wall-baseline slope.
2. That mismatch is exactly twice the weighted transfer-shape slope.
3. The whole many-port problem collapses to one effective transfer shape `T_eff^2 = N_0 / K`.
4. On the coherent branch, the support ratio drops out of the weak-axisymmetric defect exactly.
5. The remaining defect depends only on microscopic mixed/outgoing placement drifts.
6. Those drifts collapse to the three branch-adapted coordinates
   `Sigma_tr`, `Sigma_nt`, `Sigma_eta`.
7. Those coordinates are the logarithmic drifts of three exact direct microscopic monomials
   `C_tr,*`, `C_nt,*`, and `epsilon_eta`.
8. Their zero-defect set is the tangent space of an exact five-parameter multiplicative similarity orbit `G_*`.
9. At the finite level, the level sets of `(C_tr,*, C_nt,*, epsilon_eta)` are exactly the `G_*` orbits.
10. So the remaining microscopic question is purely branch-selective: whether the true moving-throat branch preserves those three quotient coordinates.

## Practical continuation point

The smallest honest next theorem gate after Stage 170 is:

- compute the actual branch drift of the three quotient coordinates
  `(C_tr,*, C_nt,*, epsilon_eta)`
  from the real moving-throat PDE branch;
- equivalently, determine whether the real branch stays on a single `G_*` similarity orbit.

If it does, the coherent first-order grouped weak-axisymmetric defect vanishes automatically.# 5PN Stages 171–174 — Adiabatic Wall, Elastic Branch Selection, and Orbit Locking

These stages implement the requested bridge from the `$g-2$` track back into the moving-throat / 5PN branch-selection problem.

## Stage 171 — Adiabatic wall bundle transport

Impose the adiabatic wall constraint
\[
\delta\ln\Theta_w=0.
\]
Using the exact Stage-149/150 inversion and bundle transport laws, the isotropic drift family collapses to
\[
\delta\ln\rho_w=\delta\ln c_{s,w}=\delta\ln\ell=0,
\]
\[
\delta\ln a=\delta\ln L_W=\frac12\,\delta\ln K_s,
\]
\[
\delta\ln c_s=\frac12\,\delta\ln K_s+\frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q=\delta\ln K_q-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln v_{w0}=-\frac34\,\delta\ln K_s+\frac12\,\delta\ln K_q,
\]
\[
\delta\ln\mathcal T_m=-\frac54\,\delta\ln K_s+\frac12\,\delta\ln K_q-\frac25\,\delta\ln P_0.
\]
The parent invariants remain tangent-compensated:
\[
\delta\ln r_c=\delta\ln\mathfrak r=\delta\ln\mathfrak g=0.
\]
So the adiabatic wall removes the wall-depth / thermal-fraying isotropic drift channel, but leaves a 3-parameter isotropic family labelled by
\[
(\delta\ln K_s,\ \delta\ln K_q,\ \delta\ln P_0).
\]

## Stage 172 — Adiabatic-elastic slippage collapse

The scalar off-bundle slippages are
\[
\varepsilon_L=\delta\ln L_W-\frac12\,\delta\ln K_s,
\]
\[
\varepsilon_v=\delta\ln v_{w0}+\frac34\,\delta\ln K_s-\frac12\,\delta\ln K_q,
\]
\[
\varepsilon_T=\delta\ln\mathcal T_m+\frac54\,\delta\ln K_s-\frac12\,\delta\ln K_q+\frac25\,\delta\ln P_0.
\]
On the adiabatic wall branch these need not vanish automatically, but if we impose the stronger elastic/no-fraying rule
\[
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
then the scalar normal defect vanishes identically:
\[
\delta_\perp=-\varepsilon_\perp=0.
\]
So the adiabatic-elastic boundary condition kills the first scalar off-bundle obstruction completely.

## Stage 173 — Why adiabatic wall alone is not enough for orbit locking

Stage-169/170 says the reduced weak-axisymmetric defect vanishes iff the microscopic grouped drift is tangent to the exact similarity orbit $\mathcal G_*$, equivalently iff the three quotient coordinates are preserved:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\qquad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]
These are encoded by the rank-3 monomial-drift map
\[
M_*\,\delta\mathbf x.
\]
Stage 173 shows explicitly that `\delta\ln\Theta_w=0` by itself does **not** imply
\[
M_*\delta\mathbf x=0.
\]
There are microscopic drift directions, such as a pure `\Delta\lambda_W` or pure `\Delta c_{\eta U}` perturbation, that leave the wall-depth condition untouched but still move the quotient coordinates. So the adiabatic-wall condition removes one failure channel, but does not by itself prove orbit locking.

## Stage 174 — Adiabatic-elastic orbit theorem

Combining the Stage-172 scalar result with the Stage-166 and Stage-169/170 quotient-closure theorem gives the clean unified statement:

If we impose
\[
\delta\ln\Theta_w=0,
\qquad
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
then the scalar off-bundle source vanishes, and the remaining first-order defect is zero **iff** the branch preserves the three quotient coordinates:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\qquad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]
Equivalently,
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
M_*\delta\mathbf x=0
\iff
\delta\mathbf x\in T_{\mathrm{id}}\mathcal G_*
\iff
\text{the branch stays on a single exact }\mathcal G_*\text{ orbit.}
\]

So the requested “unified loop” result is now explicit:

- the adiabatic wall condition freezes the thermal wall channel,
- the elastic/no-fraying condition removes the scalar off-bundle obstruction,
- and the remaining branch-selection test is **exactly** whether the moving-throat branch preserves
  \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\), i.e. stays on one \(\mathcal G_*\) orbit.

## Direct continuation point

The next clean theorem gate is no longer to manipulate isotropic bundle transport. That part is closed. The next gate is to compute the actual physical-branch drift vector and test whether its projection under `M_*` vanishes. In other words:
\[
M_*\delta\mathbf x\stackrel{?}=0.
\]
If yes, the adiabatic-elastic moving-throat branch is orbit-locked. If not, the failure is localized immediately into the tracking, nontracking-transfer, or dressing quotient directions.
# 5PN continuation — Stages 175–180

This block continues directly from the adiabatic-wall / adiabatic-elastic orbit theorem.
The earlier result had already reduced the live first-order branch-selection problem to the exact
Stage-170 quotient test
\[
M_*\,\delta\mathbf x = 0,
\]
with quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
\]

The missing step was to turn that theorem into an **exact branch-selection compiler**:
for any candidate microscopic drift, isolate the pure-similarity part, isolate the true quotient
failure, and write the unique microscopic correction that would restore single-orbit lock.

## Stage 175 — exact orbit/quotient projectors

Using the exact Stage-170 pivot block in the dependent coordinates
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
with
\[
P_{(T,K_\eta,\mu)} = M_*[:,(T,K_\eta,\mu)],
qquad
\det P = 1+\chi_{0,*}>0,
\]
one gets exact complementary projectors
\[
Q = E P^{-1} M_*,
\qquad
O = I-Q.
\]
They satisfy
\[
Q^2=Q,
\qquad
O^2=O,
\qquad
QO=OQ=0,
\qquad
M_*O=0,
\qquad
M_*Q=M_*.
\]
So every microscopic drift splits uniquely as
\[
\Delta\mathbf x = \Delta\mathbf x_{\rm orbit}+\Delta\mathbf x_{\rm fail},
\]
with
\[
\Delta\mathbf x_{\rm orbit}\in\ker M_*,
\qquad
\Delta\mathbf x_{\rm fail}=Q\Delta\mathbf x.
\]
A particularly sharp result is that the entire quotient-failure piece has support only in the
**dependent triple**
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
not in the five free similarity directions.

## Stage 176 — observable-to-microscopic correction compiler

The Stage-166/167 observable inversion gives
\[
\delta\ln\mathfrak C_{{\rm tr},*}
= -\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}}\,\Theta_1,
\]
\[
\delta\ln\mathfrak C_{{\rm nt},*}
= \Xi_1+\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}\,\Theta_1,
\]
\[
\delta\ln\epsilon_\eta
= -\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\,(\mathcal R_1+\Xi_1).
\]
Composing this with the Stage-175 projector gives the exact microscopic quotient correction
supported only on the dependent triple:
\[
\Delta_T^{(q)}=
\frac{\delta\ln\mathfrak C_{{\rm tr},*}}{1+\chi_{0,*}},
\]
\[
\Delta_{K_\eta}^{(q)}=-\delta\ln\epsilon_\eta,
\]
\[
\Delta_\mu^{(q)}
=
\delta\ln\mathfrak C_{{\rm nt},*}
-\delta\ln\epsilon_\eta
+\frac{F_*}{1+\chi_{0,*}}\,\delta\ln\mathfrak C_{{\rm tr},*}.
\]
So once the observable defect triple is known, the exact microscopic dependent-coordinate
correction needed to represent it is already fixed.

## Stage 177 — exact finite restoration operator

Because the Stage-170 fibre equations are linear in the finite log-ratios, the same projector
logic gives an exact **finite** restoration operator.
For any candidate finite drift \(\Delta\mathbf x\), define
\[
\mathbf q = M_*\Delta\mathbf x,
\qquad
\Delta\mathbf x_{\rm fail}=E P^{-1}\mathbf q,
\qquad
\Delta\mathbf x_{\rm restore}=-\Delta\mathbf x_{\rm fail}.
\]
Then
\[
M_*(\Delta\mathbf x+\Delta\mathbf x_{\rm restore})=0.
\]
So any candidate branch can be returned to a single exact \(\mathcal G_*\)-orbit by adjusting
only the dependent triple \((T_U,K_\eta,\mu_W)\).

## Stage 178 — adiabatic-elastic branch decomposition

Under the strengthened boundary rule
\[
\delta\ln\Theta_w=0,
\qquad
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
the thermal wall channel and the scalar off-bundle escape are both removed. The remaining
first-order branch-selection problem is therefore purely microscopic and eight-dimensional.
For any adiabatic-elastic candidate branch drift \(\Delta\mathbf x_{AE}\), the exact split is
\[
\Delta\mathbf x_{AE} = \Delta\mathbf x_{\rm orbit}+\Delta\mathbf x_{\rm fail},
\]
with
\[
M_*\Delta\mathbf x_{\rm orbit}=0,
\qquad
M_*\Delta\mathbf x_{\rm fail}=M_*\Delta\mathbf x_{AE}.
\]
And the weak-axisymmetric observables depend only on the quotient piece:
\[
\Theta_1
= -\frac{\chi_{0,*}\delta_{U,*}}{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
\,q_1,
\]
\[
\Xi_1
= \frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}\,q_1 + q_2,
\]
\[
\mathcal R_1+\Xi_1
= -\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_3,
\]
where \(\mathbf q=M_*\Delta\mathbf x_{AE}\).
So the adiabatic-elastic branch is orbit-locked iff
\[
\Delta\mathbf x_{\rm fail}=0
\iff
M_*\Delta\mathbf x_{AE}=0.
\]

## Stage 179 — exact dependent-coordinate mismatch formulas

Comparing an arbitrary candidate branch to the exact Stage-170 single-orbit law gives the three
microscopic mismatch coordinates
\[
m_T := \Delta_T-\Delta_T^{\rm orbit},
\qquad
m_{K_\eta}:=\Delta_{K_\eta}-\Delta_{K_\eta}^{\rm orbit},
\qquad
m_\mu:=\Delta_\mu-\Delta_\mu^{\rm orbit}.
\]
These are not arbitrary. They are exactly the quotient drifts:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=(1+\chi_{0,*})m_T,
\]
\[
\delta\ln\epsilon_\eta=-m_{K_\eta},
\]
\[
\delta\ln\mathfrak C_{{\rm nt},*}=m_\mu-F_*m_T-m_{K_\eta}.
\]
So the remaining dynamical theorem gap has been localized completely:
**the PDE only has to prove that the dependent microscopic coordinates follow the exact
Stage-170 orbit law generated by the five free similarity directions.**

## Stage 180 — consolidated adiabatic-elastic orbit-lock verdict

The strengthened boundary rule has now been compiled all the way down to a single exact
microscopic finish line.

- adiabatic wall freezes the thermal wall datum,
- elastic/no-fraying removes the scalar off-bundle escape,
- and the remaining first-order defect is nothing but the mismatch of the dependent triple
  \((\Delta_T,\Delta_{K_\eta},\Delta_\mu)\) from the exact single-orbit law.

Equivalently, the adiabatic-elastic moving-throat branch is locked to one exact
\(\mathcal G_*\)-orbit iff
\[
\Theta_1=\Xi_1=\mathcal R_1+\Xi_1=0,
\]
or microscopically iff
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)
\]
follow the Stage-170 orbit law.

## Bottom line after Stage 180

The theorem gap is now narrower than “solve the PDE.” It is:

> show that the completed moving-throat dynamics makes the dependent microscopic triple
> \((T_U,K_\eta,\mu_W)\) follow the exact Stage-170 orbit law on the adiabatic-elastic branch.

If that happens, the first-order weak-axisymmetric defect vanishes automatically and the branch
stays on a single exact \(\mathcal G_*\)-orbit.
# Stage 181–186 notes

This block continues the post-170 / post-180 branch-selection program in the most natural direction: it upgrades the **first-order** orbit-lock language to an **exact finite** law for the dependent microscopic triple.

## Main new result

The exact Stage-168 monomials

\[
\mathfrak C_{{\rm tr},*},\qquad
\mathfrak C_{{\rm nt},*},\qquad
\epsilon_\eta
\]

can be solved **exactly** for the dependent microscopic coordinates

\[
(T_U, K_\eta, \mu_W)
\]

once the five free microscopic coordinates

\[
(\lambda_W, c_{\eta U}, \gamma, K_U, K_W)
\]

and the invariant triple are fixed.

That gives the exact finite single-orbit law:

\[
K_\eta^{(\mathrm{orbit})} = \frac{c_{\eta U}^2}{K_U\,\epsilon_{\eta,*}},
\]

\[
T_U^{(\mathrm{orbit})}
= \frac{L^2 K_U}{\pi^2}
\left[
\frac{\mathfrak C_{{\rm tr},*}}
{(\gamma c_{\eta U}/K_U)^{1+\delta_{U,*}}}
\right]^{\!1/(1+\chi_{0,*})},
\]

\[
\mu_W^{(\mathrm{orbit})}
=
\frac{\mathfrak C_{{\rm nt},*} K_\eta^{(\mathrm{orbit})} K_W^2}{\lambda_W^2}
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W}\right)^{-E_*}
\left(\frac{\pi^2 T_U^{(\mathrm{orbit})}}{L^2 K_U}\right)^{F_*}.
\]

So the finite similarity orbit is no longer abstract: the dependent triple is an exact algebraic function of the free microscopic point and the invariant triple.

## Exact finite mismatch coordinates

For any candidate branch with the **same five free microscopic coordinates** as the orbit point, define

\[
T_U = m_T\,T_U^{(\mathrm{orbit})},\qquad
K_\eta = m_K\,K_\eta^{(\mathrm{orbit})},\qquad
\mu_W = m_\mu\,\mu_W^{(\mathrm{orbit})}.
\]

Then the exact invariant ratios are

\[
\frac{\mathfrak C_{{\rm tr},*}}{\mathfrak C_{{\rm tr},*}^{(\mathrm{orbit})}} = m_T^{1+\chi_{0,*}},
\qquad
\frac{\epsilon_\eta}{\epsilon_{\eta,*}} = \frac{1}{m_K},
\qquad
\frac{\mathfrak C_{{\rm nt},*}}{\mathfrak C_{{\rm nt},*}^{(\mathrm{orbit})}} = \frac{m_\mu}{m_K m_T^{F_*}}.
\]

So the finite branch-selection problem is exactly three-dimensional.

## Exact logarithmic chart

If we write

\[
\tau := \ln m_T,
\qquad
\kappa := \ln m_K,
\qquad
\mu := \ln m_\mu,
\]

then the quotient coordinates are **exactly**

\[
q_{\rm tr} = (1+\chi_{0,*})\tau,
\qquad
q_\eta = -\kappa,
\qquad
q_{\rm nt} = \mu - \kappa - F_*\tau.
\]

So the Stage-179 first-order formulas are not merely infinitesimal approximations; they are the exact logarithmic chart of the finite mismatch ratios.

## Exact restoration map

Given the finite quotient coordinates, the exact restoration to the same orbit is achieved by changing only the dependent triple:

\[
T_U^{(\mathrm{restore})} = T_U\,e^{-q_{\rm tr}/(1+\chi_{0,*})},
\]

\[
K_\eta^{(\mathrm{restore})} = K_\eta\,e^{q_\eta},
\]

\[
\mu_W^{(\mathrm{restore})}
= \mu_W\,e^{-q_{\rm nt}+q_\eta-F_* q_{\rm tr}/(1+\chi_{0,*})}.
\]

This returns the candidate branch to the exact orbit with the same free microscopic coordinates and the same invariant triple.

## Finite adiabatic-elastic orbit-lock criterion

Under the adiabatic-elastic boundary rule, the exact finite branch-selection problem is:

\[
\text{orbit lock}
\iff
m_T = m_K = m_\mu = 1
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

So after this block the remaining theorem gap is completely concrete: once the PDE gives the actual microscopic branch values, one can test orbit lock directly by comparing the candidate dependent triple to the exact orbit values above.
# 5PN Stages 187–192 — Finite Similarity-Orbit Action, Reference-Transport Laws, and Exact Residual Diagnostics

This block continues directly from the finite orbit interface at Stage 186.
The earlier stages had already shown two things:

1. the weak-axisymmetric zero-defect branch is exactly the finite similarity orbit \
   \(\mathcal G_*\), and
2. the full finite branch-selection problem can be written as a test on the dependent
   triple \((T_U,K_\eta,\mu_W)\).

What was still missing was the **finite transport law** itself: how the five free microscopic
coordinates move a branch along an exact orbit, how to compare an actual candidate branch to a
reference orbit point, and how to localize any failure into exact multiplicative residuals.

## Stage 187 — exact finite similarity-orbit action

The five-parameter multiplicative similarity orbit \(\mathcal G_*\) is now written explicitly as a
finite action on the full microscopic state. If
\[
(\lambda_W,c_{\eta U},\gamma,K_U,K_W)
\to
(e^{\Lambda}\lambda_W,e^C c_{\eta U},e^{\Gamma}\gamma,e^U K_U,e^W K_W),
\]
then the dependent triple coevolves by
\[
K_\eta' = e^{2C-U}K_\eta,
\]
\[
T_U' = e^{U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}T_U,
\]
\[
\mu_W'
=
\exp\!\Bigl[
2C-U+2W-2\Lambda
-E_*\bigl(2\Gamma+2\Lambda-U-W\bigr)
-F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)
\Bigr]\mu_W.
\]
The three quotient monomials
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\]
are preserved **exactly**, not only infinitesimally.

## Stage 188 — exact group law, inverse, and parameter recovery

The finite orbit action is an exact abelian five-parameter group:
composition adds the five generator exponents, and inversion negates them. If two states are on
one orbit, the orbit parameters are recovered exactly from the free-coordinate ratios,
\[
\Lambda=\ln\frac{\lambda_W'}{\lambda_W},
\qquad
C=\ln\frac{c_{\eta U}'}{c_{\eta U}},
\qquad
\Gamma=\ln\frac{\gamma'}{\gamma},
\qquad
U=\ln\frac{K_U'}{K_U},
\qquad
W=\ln\frac{K_W'}{K_W}.
\]
So once those five free-coordinate ratios are known, the dependent triple on the same orbit is
predicted uniquely.

## Stage 189 — exact reference-orbit transport laws

Relative to a reference orbit point, the exact dependent-coordinate transport laws are
\[
R_{K_\eta}^{(\mathrm{orbit})}=\frac{R_c^2}{R_U},
\]
\[
R_{T_U}^{(\mathrm{orbit})}=R_U\left(\frac{R_U}{R_\gamma R_c}\right)^{\frac{1+\delta_{U,*}}{1+\chi_{0,*}}},
\]
\[
R_{\mu_W}^{(\mathrm{orbit})}
=
\frac{R_{K_\eta}^{(\mathrm{orbit})}R_W^2}{R_\lambda^2}
\left(\frac{R_\gamma^2R_\lambda^2}{R_UR_W}\right)^{-E_*}
\left(\frac{R_{T_U}^{(\mathrm{orbit})}}{R_U}\right)^{F_*}.
\]
This is the exact finite coevolution law of the dependent triple along a fixed orbit.

## Stage 190 — exact dependent residual coordinates

A general candidate branch with the same five free-coordinate ratios can be factored as
\[
R_{T_U}^{(\mathrm{actual})}=R_{T_U}^{(\mathrm{orbit})}m_T,
\qquad
R_{K_\eta}^{(\mathrm{actual})}=R_{K_\eta}^{(\mathrm{orbit})}m_K,
\qquad
R_{\mu_W}^{(\mathrm{actual})}=R_{\mu_W}^{(\mathrm{orbit})}m_\mu,
\]
where \((m_T,m_K,m_\mu)\) is the **dependent residual mismatch triple**.
The invariant ratios then collapse exactly to
\[
\frac{\mathfrak C_{{\rm tr},*}^{\mathrm{actual}}}{\mathfrak C_{{\rm tr},*}^{\mathrm{ref}}}=m_T^{1+\chi_{0,*}},
\qquad
\frac{\epsilon_\eta^{\mathrm{actual}}}{\epsilon_{\eta}^{\mathrm{ref}}}=\frac{1}{m_K},
\qquad
\frac{\mathfrak C_{{\rm nt},*}^{\mathrm{actual}}}{\mathfrak C_{{\rm nt},*}^{\mathrm{ref}}}=\frac{m_\mu}{m_Km_T^{F_*}}.
\]
So the free-coordinate transport along an orbit drops out completely; the quotient coordinates are
nothing but the logarithmic chart of the dependent residual triple.

## Stage 191 — factorized actual-branch interface

The actual candidate branch now admits an exact factorization
\[
\text{actual branch}
=
(\text{reference orbit point})
\times
(\text{free-coordinate orbit transport})
\times
(\text{dependent residual mismatch}).
\]
Restoration to the same orbit at fixed free-coordinate ratios is achieved by dividing only by the
residual mismatch ratios:
\[
K_\eta^{(\mathrm{restore})}=\frac{K_\eta^{(\mathrm{actual})}}{m_K},
\qquad
T_U^{(\mathrm{restore})}=\frac{T_U^{(\mathrm{actual})}}{m_T},
\qquad
\mu_W^{(\mathrm{restore})}=\frac{\mu_W^{(\mathrm{actual})}}{m_\mu}.
\]
So orbit lock is exactly the statement
\[
m_T=m_K=m_\mu=1.
\]

## Stage 192 — diagnostic signatures of each failure channel

The three quotient coordinates now have a direct physical interpretation:
\[
q_{\rm tr}=(1+\chi_{0,*})\ln m_T,
\qquad
q_\eta=-\ln m_K,
\qquad
q_{\rm nt}=\ln m_\mu-\ln m_K-F_*\ln m_T.
\]
That gives three clean signatures:

- pure \(T_U\) residual mismatch turns on \(q_{\rm tr}\) and, via \(F_*\), also \(q_{\rm nt}\),
- pure \(K_\eta\) residual mismatch turns on \(q_\eta\) and also \(q_{\rm nt}\),
- pure \(\mu_W\) residual mismatch turns on \(q_{\rm nt}\) only.

So once the PDE delivers an actual candidate branch, the pattern of
\[
(q_{\rm tr},q_\eta,q_{\rm nt})
\]
identifies exactly which dependent coevolution law failed, if any.

## Bottom line after Stage 192

The branch-selection problem is no longer merely
“compute the actual branch and see if the quotient coordinates move.”
It is now a **factorized finite comparison problem**:

1. choose a reference point on the exact similarity orbit,
2. use the five free-coordinate ratios to predict the orbit-transported dependent triple,
3. compare the actual dependent triple to that prediction,
4. read off the residual mismatch ratios \((m_T,m_K,m_\mu)\),
5. and infer the quotient coordinates and restoration map immediately.

So the next PDE theorem gate is now even sharper:

> compute the actual branch values of the eight microscopic coordinates, form the free-coordinate
> ratios and the dependent residual mismatch triple, and test whether
> \(m_T=m_K=m_\mu=1\).

If yes, the branch stays on a single exact \(\mathcal G_*\) orbit. If not, the failure is localized
exactly into the \(T_U\), \(K_\eta\), and/or \(\mu_W\) transport laws.

# Stages 193–198 — Pairwise orbit-lock, cocycle laws, and minimal-data verdict

This block continues the finite orbit-lock chain after Stages 181–192 by making the
orbit criterion fully **reference-independent** and **PDE-ready**.

## Stage 193 — exact pairwise orbit criterion

Given two positive microscopic states `x` and `y`, with shared branch constants
`(chi0_*, deltaU_*, E_*, F_*)`, the five free-coordinate ratios
\[
R_\lambda,\ R_c,\ R_\gamma,\ R_U,\ R_W
\]
determine the exact orbit-predicted dependent-coordinate ratios
\[
R_{K_\eta}^{(\mathrm{orbit})}=\frac{R_c^2}{R_U},
\]
\[
R_{T_U}^{(\mathrm{orbit})}
=
R_U\left(\frac{R_U}{R_\gamma R_c}\right)^{\frac{1+\delta_{U,*}}{1+\chi_{0,*}}},
\]
\[
R_{\mu_W}^{(\mathrm{orbit})}
=
\frac{R_{K_\eta}^{(\mathrm{orbit})}R_W^2}{R_\lambda^2}
\left(\frac{R_\gamma^2R_\lambda^2}{R_UR_W}\right)^{-E_*}
\left(\frac{R_{T_U}^{(\mathrm{orbit})}}{R_U}\right)^{F_*}.
\]

Comparing with the actual dependent-coordinate ratios defines the pairwise residual triple
\[
m_T^{(x\to y)}=\frac{R_{T_U}^{(\mathrm{act})}}{R_{T_U}^{(\mathrm{orbit})}},
\qquad
m_K^{(x\to y)}=\frac{R_{K_\eta}^{(\mathrm{act})}}{R_{K_\eta}^{(\mathrm{orbit})}},
\qquad
m_\mu^{(x\to y)}=\frac{R_{\mu_W}^{(\mathrm{act})}}{R_{\mu_W}^{(\mathrm{orbit})}}.
\]

The invariant ratios between `y` and `x` are then exactly
\[
\frac{\mathfrak C_{{\rm tr},*}(y)}{\mathfrak C_{{\rm tr},*}(x)}
=
\bigl(m_T^{(x\to y)}\bigr)^{1+\chi_{0,*}},
\]
\[
\frac{\epsilon_\eta(y)}{\epsilon_\eta(x)}
=
\frac{1}{m_K^{(x\to y)}},
\]
\[
\frac{\mathfrak C_{{\rm nt},*}(y)}{\mathfrak C_{{\rm nt},*}(x)}
=
\frac{m_\mu^{(x\to y)}}{m_K^{(x\to y)}\bigl(m_T^{(x\to y)}\bigr)^{F_*}}.
\]

So `x` and `y` lie on the same exact similarity orbit iff
\[
m_T^{(x\to y)}=m_K^{(x\to y)}=m_\mu^{(x\to y)}=1,
\]
equivalently iff the three invariant ratios are all unity.

## Stage 194 — multiplicative cocycle and additive quotient law

For three states `x,y,z`, the residual ratios compose exactly:
\[
m_T^{(x\to z)}=m_T^{(x\to y)}m_T^{(y\to z)},
\qquad
m_K^{(x\to z)}=m_K^{(x\to y)}m_K^{(y\to z)},
\qquad
m_\mu^{(x\to z)}=m_\mu^{(x\to y)}m_\mu^{(y\to z)}.
\]

The logarithmic quotient coordinates therefore add:
\[
q_{\rm tr}^{(x\to z)}=q_{\rm tr}^{(x\to y)}+q_{\rm tr}^{(y\to z)},
\]
\[
q_{\rm nt}^{(x\to z)}=q_{\rm nt}^{(x\to y)}+q_{\rm nt}^{(y\to z)},
\]
\[
q_\eta^{(x\to z)}=q_\eta^{(x\to y)}+q_\eta^{(y\to z)},
\]
with inverse laws under reversal,
\[
q^{(y\to x)}=-q^{(x\to y)}.
\]

So a sequence of PDE branch snapshots can be tracked either multiplicatively in the
residual ratios or additively in the quotient coordinates.

## Stage 195 — pairwise quotient-to-observable compiler

From the residual ratios,
\[
q_{\rm tr}=(1+\chi_{0,*})\ln m_T,\qquad
q_\eta=-\ln m_K,\qquad
q_{\rm nt}=\ln m_\mu-\ln m_K-F_*\ln m_T.
\]

Composing with the Stage-170 linear observable map gives the small pairwise observable
signature
\[
\Theta_1^{(\mathrm{lin})}
=
-\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
\,q_{\rm tr},
\]
\[
\Xi_1^{(\mathrm{lin})}
=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}\,q_{\rm tr}+q_{\rm nt},
\]
\[
(\mathcal R_1+\Xi_1)^{(\mathrm{lin})}
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_\eta.
\]

The pure-channel signatures become:

- pure `T_U` mismatch:
  \[
  q_{\rm tr}=(1+\chi_{0,*})\ln z_T,\qquad
  q_{\rm nt}=-F_*\ln z_T,\qquad
  q_\eta=0;
  \]

- pure `K_\eta` mismatch:
  \[
  q_{\rm tr}=0,\qquad
  q_{\rm nt}=-\ln z_K,\qquad
  q_\eta=-\ln z_K;
  \]

- pure `\mu_W` mismatch:
  \[
  q_{\rm tr}=0,\qquad
  q_{\rm nt}=\ln z_M,\qquad
  q_\eta=0.
  \]

This is the finite pairwise uplift of the earlier Stage-192 channel diagnostics.

## Stage 196 — exact inversion from invariant ratios

The invariant-ratio packet
\[
(R_{\rm tr},R_{\rm nt},R_\eta)
:=
\left(
\frac{\mathfrak C_{{\rm tr},*}(y)}{\mathfrak C_{{\rm tr},*}(x)},
\frac{\mathfrak C_{{\rm nt},*}(y)}{\mathfrak C_{{\rm nt},*}(x)},
\frac{\epsilon_\eta(y)}{\epsilon_\eta(x)}
\right)
\]
already contains the full orbit-lock verdict.

The exact inversion is
\[
m_T = R_{\rm tr}^{1/(1+\chi_{0,*})},
\qquad
m_K = \frac{1}{R_\eta},
\qquad
m_\mu = R_{\rm nt}\,m_K\,m_T^{F_*}.
\]

So once the invariant ratios are known, the five free-coordinate ratios are not needed
for the verdict or for the dependent-coordinate restoration.

## Stage 197 — canonical orbit-distance quadratic form

Write the logarithmic residuals as
\[
t=\ln m_T,\qquad k=\ln m_K,\qquad \mu=\ln m_\mu.
\]

Then
\[
\begin{pmatrix}
q_{\rm tr}\\ q_{\rm nt}\\ q_\eta
\end{pmatrix}
=
A
\begin{pmatrix}
t\\k\\\mu
\end{pmatrix},
\qquad
A=
\begin{pmatrix}
1+\chi_{0,*} & 0 & 0\\
-F_* & -1 & 1\\
0 & -1 & 0
\end{pmatrix}.
\]

The canonical scalar orbit-distance is
\[
D^2=q_{\rm tr}^2+q_{\rm nt}^2+q_\eta^2
=
\begin{pmatrix}
t&k&\mu
\end{pmatrix}
Q
\begin{pmatrix}
t\\k\\\mu
\end{pmatrix},
\qquad
Q=A^TA,
\]
i.e.
\[
Q=
\begin{pmatrix}
(1+\chi_{0,*})^2+F_*^2 & F_* & -F_*\\
F_* & 2 & -1\\
-F_* & -1 & 1
\end{pmatrix}.
\]

Its principal minors are
\[
(1+\chi_{0,*})^2+F_*^2>0,
\]
\[
2(1+\chi_{0,*})^2+F_*^2>0,
\]
\[
(1+\chi_{0,*})^2>0,
\]
so `Q` is positive definite on the constructive branch.

Therefore
\[
D^2=0
\iff
m_T=m_K=m_\mu=1.
\]

So the entire pairwise orbit-lock failure can be summarized by one exact reference-free
positive scalar.

## Stage 198 — minimal-data orbit verdict

The full orbit-lock verdict can be reached from **any one** of three exact packets:

1. residual mismatch ratios:
   \[
   (m_T,m_K,m_\mu),
   \]

2. invariant ratios:
   \[
   (R_{\rm tr},R_{\rm nt},R_\eta),
   \]

3. quotient coordinates:
   \[
   (q_{\rm tr},q_{\rm nt},q_\eta).
   \]

They are exactly interconvertible:
\[
R_{\rm tr}=m_T^{1+\chi_{0,*}},\qquad
R_{\rm nt}=\frac{m_\mu}{m_K m_T^{F_*}},\qquad
R_\eta=\frac{1}{m_K},
\]
\[
q_{\rm tr}=\ln R_{\rm tr},\qquad
q_{\rm nt}=\ln R_{\rm nt},\qquad
q_\eta=\ln R_\eta,
\]
\[
m_T=e^{q_{\rm tr}/(1+\chi_{0,*})},\qquad
m_K=e^{-q_\eta},\qquad
m_\mu=e^{q_{\rm nt}-q_\eta+F_*q_{\rm tr}/(1+\chi_{0,*})}.
\]

So future PDE numerics only need to provide whichever packet is cleanest. From that
packet one can reconstruct the dependent-coordinate restoration map and the scalar
distance `D^2`.

# Stages 199–201 — Bringing the 5PN endgame home

This block stops the long exploratory chain and compresses the remaining theorem gap
into the smallest exact packets the completed moving-throat PDE must still supply.

## Stage 199 — exact final branch residual packet

Take the actual grouped-lane low-frequency bundle data
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad A\in\{20,21,22\},
\]
together with the source-map factor \(m_{\hat 0}\).

Compile the normalized grouped response moments
\[
u_2^{(A)}=-\frac{D_{A2}}{D_{A0}},
\qquad
u_4^{(A)}=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2},
\]
and the outgoing prefactor moments
\[
P_0^{(A)}=\frac{N_{A0}}{D_{A0}},
\]
\[
P_2^{(A)}=\frac{D_{A0}N_{A2}-2D_{A2}N_{A0}}{D_{A0}^2},
\]
\[
P_4^{(A)}=
\frac{D_{A0}^2N_{A4}-2D_{A0}(D_{A2}N_{A2}+D_{A4}N_{A0})+3D_{A2}^2N_{A0}}
{D_{A0}^3}.
\]

Then extract the grouped trace/anomaly data
\[
(\bar u_2,a_2,b_2),\qquad
(\bar u_4,a_4,b_4),\qquad
(\bar P_0,a_{P_0},b_{P_0}).
\]

The exact final branch residual packet is
\[
\Delta_{\rm branch}
=
\bigl(
a_2,\ b_2,\ a_4,\ b_4,\ a_{P_0},\ b_{P_0},\ \Delta_{\rm pole},\ \Delta_{\rm norm}
\bigr),
\]
with
\[
\Delta_{\rm pole}=\bar u_4-4\bar u_2^{\,2},
\]
\[
\Delta_{\rm norm}=m_{\hat 0}^{\,2}\bar P_0-\frac{54Gc_s^5}{5a^5c^5}.
\]

So the completed PDE no longer has to “show 5PN somehow.” It has to drive one exact
finite-dimensional residual packet to zero.

## Stage 200 — exact endgame compiler

The orbit side was already reduced in Stages 181–198 to any one of the equivalent packets
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

This stage combines that orbit packet with \(\Delta_{\rm branch}\) and shows that the
whole reduced 5PN / 2.5PN / 4PN closure problem depends only on

1. the grouped branch packet \(\Delta_{\rm branch}\),
2. the orbit-lock packet \(\Delta_{\rm orbit}\).

It also records one useful practical simplification: on the minimal isotropic conservative
module, the explicit Family‑1 support/source side is already above the required threshold,
because
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
A_{F1}\approx 1.00005192880220 > \frac13.
\]
So the support/source side is no longer the active bottleneck inside the current hierarchy.

The reduced closure problem is therefore exactly
\[
\Delta_{\rm branch}=0,
\qquad
\Delta_{\rm orbit}=0.
\]

## Stage 201 — home-stretch theorem and minimal PDE data packet

The final theorem gate is now stated as a minimal-data problem.

### Packet A — grouped bundle data
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad
A\in\{20,21,22\},
\]
plus \(m_{\hat 0}\).

### Packet B — orbit/invariant data
Any one of
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

Everything else in the reduced closure test is an exact compiler output of these packets.

So the remaining theorem gap is no longer diffuse. The completed moving-throat PDE only
has to supply the data needed to evaluate
\[
\Delta_{\rm branch}
\quad\text{and}\quad
\Delta_{\rm orbit}.
\]

If \(\Delta_{\rm branch}\neq0\), the branch fails the reduced GR test.
If \(\Delta_{\rm branch}=0\) but \(\Delta_{\rm orbit}\neq0\), the branch is isotropic
but off the exact similarity orbit.
Only if both vanish is the reduced closure complete inside the current hierarchy.
# 5PN progress notes — explicit overlap model and weak-axisymmetric transport

## What was added

### Stage 3 — explicit isotropic finite-throat overlap model

Files:
- `5pn_stage3_isotropic_overlap_model.py`
- `5pn_stage3_isotropic_overlap_model_output.txt`

Model choice:
- finite interval `s in [0,L]`
- wall profile
  `chi_eta(s) = sqrt(2/L) sin(pi s / L)`
- D/N-like half-wave
  `phi_DN(s) = sqrt(2/L) sin(pi s / (2L))`
- simplest brane-like gauge profile `u = chi_eta`
- simplest mixed profile `w = phi_DN`

Exact overlaps:
- `I_eta_u = 1`
- `I_eta_phi = I_eta_w = I_uw = 8/(3 pi)`

Using one BdG mode and one conservative Maxwell/mixed pair, the script constructs
- `B0,B2,B4`
- `Z0,Z2,Z4`
- `N0,N2,N4`
- `D0,D2,D4`
- `u2,u4`
- `P0,P2,P4`

and verifies exact grouped isotropy when the three grouped lanes share the same radial/axial data.

### Most useful Stage-3 result

In this explicit prototype, two important requirements both fix the same static wall stiffness `K`:

1. the 5PN conservative one-pole identity
   `u4 = 4 u2^2`
2. the universal 2.5PN / 4PN normalization target for `P0`

So simultaneous success reduces to one explicit compatibility equation:

`N0 / P0_target = 3 (M + B2 + Z2)^2 / (B4 + Z4)`

This is the first concrete radial/axial compatibility surface I have for the 5PN program.
It is the explicit prototype version of “the PDE must make the conservative and outgoing calibrations agree on one branch.”

## Stage 4 — weak-axisymmetric transport and microscopic grouped obstructions

Files:
- `5pn_stage4_axisymmetric_transport.py`
- `5pn_stage4_axisymmetric_transport_output.txt`

This stage starts from the exact weak-axisymmetric grouped signature
- `lambda_(20) = 1`
- `lambda_(21) = 1/2`
- `lambda_(22) = -1`

and carries the grouped operator / transfer slopes through the full Stage-154/155/156 logic.

### Exact microscopic slope decomposition

Using
- wall slopes `(K1_wall, M1_wall)`
- BdG slopes `(B01, B21, B41)`
- conservative Maxwell/mixed slopes `(Z01, Z21, Z41)`
- outgoing transfer slope `N01`

it verifies

`K_1 = D21 + D01/9 = W1 - Bcal1 - Zcal1`

and

`G_1 = N01 - P0 D01 = -P0 K1_wall + P0 B01 + Nbundle1`

with
- `D01 = K1_wall - B01 - Z01`
- `D21 = -(M1_wall + B21 + Z21)`
- `D41 = -(B41 + Z41)`

### Physical weak-axisymmetric slopes

The script derives

`u2^(1) = -(D21 + u2 D01)/D0`

`P1/P0 = N01/N0 - D01/D0`

and on the canonical compensated branch verifies the hidden-even relation

`u4^(1) = (8/9) u2^(1)`

iff

`D41 = (2/3) D21 + D01/27`

### Most useful Stage-4 result

On the even-preserving branch `u2^(1)=0`, the conservative grouped response collapses to

- `D21 = -D01/9`
- `D41 = -D01/27`

and the remaining linear grouped `2.5`PN defect becomes one scalar loading mismatch

`Xi_load = N01/N0 - D01/D0 = P1/P0`

with fixed grouped-lane signature

- `(20,21,22) ~ (1, 1/2, -1)`
- equivalently `b = 3 a`

So after Stage 4, the next theorem gate is no longer “solve all grouped anisotropies.”
It is much narrower:

> compute `D01` and `N01` — or directly `Xi_load` — from a primitive weak-axisymmetric deformation of the explicit finite-throat overlap model.

## Run status

Both new scripts completed here without running out of time or memory.

## Stage 5 — primitive deformation and exact compensation surfaces

Files:
- `5pn_stage5_primitive_deformation_compensation.py`
- `5pn_stage5_primitive_deformation_compensation_output.txt`

This stage takes the explicit isotropic overlap model from Stage 3 and adds a primitive weak-axisymmetric microscopic deformation through the slopes

- `dK`, `dM`
- `d(lambda_B)`, `d(varpi)`
- `d(lambda_U)`, `d(lambda_W)`, `d(lambda_R)`
- `d(Omega_U)`, `d(Omega_W)`

Then it computes the induced grouped-lane slope data

- `D01`, `D21`, `D41`
- `N01`

and from them the three key combinations

- `K1 = D21 + D01/9`
- `G1 = N01 - P0 D01`
- `Xi_load = N01/N0 - D01/D0 = G1/N0`

### Most useful Stage-5 results

1. **Even-preserving compensation** `K1 = 0` is algebraic and fixes the inertia-side slope `dM` exactly.
2. **Odd/normalization-preserving compensation** `Xi_load = 0` is algebraic and fixes the static loading slope `dK` exactly.
3. The remaining explicit **5PN hidden-even residual** is
   `D41 - (2/3) D21 - D01/27`.

So after Stage 5, the next theorem gate is now extremely sharp:

> choose a concrete primitive anisotropy mechanism — wall-only, BdG-only, Maxwell/mixed-only, or a mixed combination — substitute its slopes into Stage 5, and test
> `K1 = 0`, `Xi_load = 0`, and `D41 - (2/3) D21 - D01/27 = 0`.

## Run status

All three new scripts completed here without running out of time or memory.
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

# 5PN continuation notes — stages 8 through 11

These stages continue the weak-axisymmetric grouped-`P2` program from the Stage-7 scalar amplitude
\[
\Xi_1=\frac{P_1}{P_0}.
\]

The overall effect is that the remaining first-order grouped normalization problem is no longer a generic “compute lots of microscopic drifts” problem. It has collapsed to a small sequence of exact equivalent formulations.

## Stage 8 — direct outgoing-port co-loading

The remaining grouped weak-axisymmetric defect can be written directly in terms of the actual outgoing-port slopes:
\[
\Xi_1
=
\sum_r \rho_r^{(N)}(\nu_r-\kappa_1)
=
\bar\nu_N-\kappa_1,
\]
where
\[
\nu_r = 2(\mathfrak p_r-\mathfrak d_r)
\]
is the logarithmic slope of
\[
N_{0}^{(r)}=\frac{P_r^2}{\Delta_r^2}.
\]

Equivalently, if
\[
N_{A,0}^{(r)} = K_A \mathcal T_{A,r}^2,
\qquad
\delta\ln \mathcal T_{A,r}=\epsilon\lambda_A\,\tau_r,
\]
then
\[
\nu_r=\kappa_1+2\tau_r,
\qquad
\Xi_1 = 2\sum_r \rho_r^{(N)}\tau_r.
\]

So the exact zero-defect condition is now:
\[
\sum_r \rho_r^{(N)}\tau_r=0.
\]

A stronger sufficient condition is per-port co-loading:
\[
\tau_r=0
\qquad
\text{for every active outgoing port }r.
\]

Under upper-leg and coupling rigidity, the transfer-shape slope collapses to the raw mixed-leg slope and recovers the old square-root mixed-leg law.

## Stage 9 — one effective transfer shape and the actual continuum branch

The many-port weighted sum collapses exactly to one effective transfer shape:
\[
\mathcal T_{\mathrm{eff},A}^2
=
\sum_r \mathcal T_{A,r}^2
=
\frac{N_{A,0}}{K_A},
\qquad
\Xi_1
=
\frac{\delta\ln\mathcal T_{\mathrm{eff},A}^2}{\epsilon\lambda_A}.
\]

On the actual one-port continuum branch,
\[
\mathcal T_A^2
=
\frac{Z_{W,A}(1+\rho_A)^2}{\Omega_{W,A}^2(1-\epsilon_{W,A})^2},
\]
so
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\rho_1}{1+\rho}+\frac{2\varepsilon_W}{1-\epsilon_W}.
\]

In selected-branch form,
\[
\mathcal T_A^2
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\,
\frac{1-\epsilon_{\eta,A}}{R_{\mathrm{target},A}},
\]
which yields
\[
\Xi_1
=
-\frac{\eta_1}{1-\epsilon_\eta}
-\mathcal R_1.
\]

On the coherent local D/N branch,
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\epsilon
=
\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
and the defect becomes
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}
+\frac{2\epsilon_1}{1-\epsilon}.
\]

The support ratio drops out identically:
\[
\partial_\zeta \ln \mathcal T^2 = 0.
\]

So the coherent defect is support-blind, and exact tracking rigidity by itself is not enough to kill it.

## Stage 10 — microscopic coherent-kernel slippages and exact triangular normal form

The coherent branch depends only on the microscopic slippages
\[
\Sigma_Z,\quad
\Sigma_\chi,\quad
\Sigma_\epsilon,\quad
\Sigma_\delta,
\]
with one additional dressing slippage
\[
\Sigma_\eta
\]
entering the selected-branch form.

The exact microscopic grouped-defect law is
\[
\Xi_1
=
\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta
\right].
\]

The exact tracking combination is
\[
\Sigma_{\rm tr}
=
(1+\chi_0)\Sigma_\delta
+
(1+\delta_U)\Sigma_\chi,
\]
with
\[
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
\]

Defining the genuine nontracking slippage
\[
\Sigma_{\rm nt},
\]
the coherent problem takes the exact triangular normal form
\[
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
\]

So on the constructive coherent branch,
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.
\]

This is the first exact three-scalar normal form of the full coherent weak-axisymmetric problem.

## Stage 11 — direct microscopic monomials, similarity orbit, and quotient closure

The three branch-adapted coordinates can be written as logarithmic drifts of three exact microscopic monomials:
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=\Sigma_\eta,
\]
with
\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\left(\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}\right)
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
\]

The exact monomial-drift map is the rank-3 matrix
\[
M_*\,\delta\mathbf x
=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix},
\]
with
\[
\det M_*^{(\tau,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0,
\]
so
\[
\dim\ker M_* = 5.
\]

There is an exact five-parameter multiplicative similarity orbit \(\mathcal G_*\) preserving the three monomials exactly, and the scripts show
\[
\ker M_* = T_{\mathrm{id}}\mathcal G_*.
\]

More strongly, the exact finite invariant-fibre equations are
\[
M_*\,\Delta \mathbf x = 0,
\]
and solving them reproduces the exact orbit exponents. Therefore the finite level sets of
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\]
are precisely the similarity orbits \(\mathcal G_*\).

So the coherent weak-axisymmetric zero-defect theorem can now be written as
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x \in T_{\mathrm{id}}\mathcal G_*,
\]
and, at finite level,
\[
\mathcal M_+/\mathcal G_*
\cong
(\mathbb R_{>0})^3
\]
with quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
\]

## What this means for the 5PN program

At this point the weak-axisymmetric normalization problem is no longer an algebraic bottleneck.

There are now three equivalent ways to state the next theorem gate:

1. compute the direct grouped scalar
   \[
   \Xi_1=\frac{P_1}{P_0},
   \]
2. compute the branch-adapted defect coordinates
   \[
   (\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta),
   \]
3. or test whether the actual moving-throat weak-axisymmetric branch is tangent to the exact monomial-preserving similarity orbit \(\mathcal G_*\).

That last form is the sharpest current continuation point.
# 5PN continuation notes — Stage 12–13 normalized monomial bridge

These two stages connect the earlier explicit Stage-5 primitive prototype

a) directly to the Stage-10/11 microscopic similarity-orbit package, and
b) into formulas that can be used on the actual moving-throat branch.

## Stage 12 — normalized prototype already contains the Stage-11 quotient coordinates

Using the coherent-kernel dictionary

- `K = K_eta^(eff)`
- `M = mu_eta`
- `G_U = c_etaU / sqrt(mu_eta mu_U)`
- `G_W = lambda_W / sqrt(mu_eta mu_W)`
- `R = gamma lambda_W / sqrt(mu_U mu_W)`
- `Omega_U^2 = K_U / mu_U`
- `Omega_W^2 = K_W^(eff) / mu_W`
- `delta_U = pi^2 T_U / (L^2 K_U)`

one gets the exact normalized coherent ratios

\[
\chi_0 = \frac{R G_U}{\Omega_U^2 G_W},
\qquad
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2},
\qquad
\epsilon_W = \frac{R^2\sigma}{\Omega_U^2\Omega_W^2},
\qquad
Z_W = \frac{M G_W^2}{K\Omega_W^2}.
\]

The direct Stage-168/169 monomials become

\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{R G_U}{\Omega_U^2 G_W}\right)^{1+\delta_{U,*}}
\delta_U^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\frac{M G_W^2}{K\Omega_W^4}
\left(\frac{R^2\sigma}{\Omega_U^2\Omega_W^2}\right)^{E_*}
\delta_U^{-F_*},
\]

\[
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2}.
\]

So the Stage-5-style normalized prototype already contains the exact Stage-11 quotient coordinates once the split-`U` variable `delta_U` is admitted.

The Stage-10 slippages collapse to the normalized drift formulas

\[
\Sigma_Z = d\ln M + 2 d\ln G_W - d\ln K - 4 d\ln\Omega_W,
\]
\[
\Sigma_\chi = d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U,
\]
\[
\Sigma_\eta = d\ln M + 2 d\ln G_U - d\ln K - 2 d\ln\Omega_U,
\]
\[
\Sigma_\epsilon = 2(d\ln R - d\ln\Omega_U - d\ln\Omega_W),
\qquad
\Sigma_\delta = d\ln\delta_U.
\]

So the extra raw mass/stiffness bookkeeping mostly cancels out of the actual defect coordinates.

## Stage 13 — exact zero-defect kernel in normalized Stage-5 variables

In the normalized variables

\[
(d\ln G_W,
 d\ln G_U,
 d\ln R,
 d\ln K,
 d\ln M,
 d\ln\Omega_U,
 d\ln\Omega_W,
 d\ln\delta_U),
\]

the direct monomial-drift matrix is

\[
M_{\rm norm}=
\begin{pmatrix}
-(1+\delta_U) & 1+\delta_U & 1+\delta_U & 0 & 0 & -2(1+\delta_U) & 0 & 1+\chi_0 \\
2 & 0 & 2E_* & -1 & 1 & -2E_* & -(4+2E_*) & -F_* \\
0 & 2 & 0 & -1 & 1 & -2 & 0 & 0
\end{pmatrix}.
\]

Its rank is exactly `3`, so the normalized zero-defect tangent space has dimension `5`.

The exact zero-defect equations are

\[
(1+\delta_U)(d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U)
+
(1+\chi_0)d\ln\delta_U = 0,
\]

\[
d\ln M + 2 d\ln G_U - d\ln K - 2 d\ln\Omega_U = 0,
\]

\[
d\ln M + 2 d\ln G_W - d\ln K - 4 d\ln\Omega_W
+
2E_*(d\ln R - d\ln\Omega_U - d\ln\Omega_W)
-
F_* d\ln\delta_U = 0.
\]

These solve triangularly:

### Tracking
\[
d\ln\delta_U
=
-\frac{1+\delta_U}{1+\chi_0}
\bigl(d\ln R + d\ln G_U - d\ln G_W - 2 d\ln\Omega_U\bigr).
\]

### Dressing
\[
d\ln M
=
 d\ln K - 2 d\ln G_U + 2 d\ln\Omega_U.
\]

### Nontracking
\[
d\ln\Omega_W
=
\frac{d\ln G_W - d\ln G_U + (1-E_*)d\ln\Omega_U + E_* d\ln R - \tfrac{F_*}{2}d\ln\delta_U}{E_*+2}.
\]

So once a candidate moving-throat branch gives the drifts of

- `K`
- `G_U`
- `G_W`
- `R`
- `Omega_U`

it automatically fixes the drifts required in

- `delta_U`
- `M`
- `Omega_W`

to stay tangent to the exact similarity orbit.

## Practical Stage-5 absolute-slope form

Writing the Stage-5 primitive slopes as

- `dK`
- `dM`
- `d(lambda_U)`
- `d(lambda_W)`
- `d(lambda_R)`
- `d(Omega_U)`
- `d(Omega_W)`
- `d(delta_U)`

one gets the exact compatibility formulas

\[
d(\delta_U)
=
-\delta_U\frac{1+\delta_U}{1+\chi_0}
\left[
\frac{d\lambda_R}{\lambda_R}
+
\frac{d\lambda_U}{\lambda_U}
-
\frac{d\lambda_W}{\lambda_W}
-
2\frac{d\Omega_U}{\Omega_U}
\right],
\]

\[
dM
=
M\left[
\frac{dK}{K} - 2\frac{d\lambda_U}{\lambda_U} + 2\frac{d\Omega_U}{\Omega_U}
\right],
\]

\[
\frac{d\Omega_W}{\Omega_W}
=
\frac{1}{E_*+2}
\left[
\frac{d\lambda_W}{\lambda_W}
-
\frac{d\lambda_U}{\lambda_U}
+
(1-E_*)\frac{d\Omega_U}{\Omega_U}
+
E_*\frac{d\lambda_R}{\lambda_R}
-
\frac{F_*}{2}\frac{d(\delta_U)}{\delta_U}
\right].
\]

So the Stage-10/11 similarity-orbit theorem is now directly usable in the Stage-5 primitive deformation language.

## Immediate consequence for the 5PN program

The next honest theorem gate is now smaller than before:

1. extract the actual branch drifts of
   `K, lambda_U, lambda_W, lambda_R, Omega_U`
   from the moving-throat PDE,
2. use the formulas above to predict the required co-drifts of
   `delta_U, M, Omega_W`
   if the branch is tangent to the exact zero-defect similarity orbit,
3. then compare those predictions with the actual reduced PDE branch.

If they agree, the weak-axisymmetric first-order obstruction is pure similarity-gauge and the calibration survives. If they fail, the moving-throat branch leaves the exact monomial-preserving orbit and the calibration breaks for a concrete microscopic reason.

## Stage 14 — the BdG primitive drifts are exactly support-blind for the Stage-11 monomials

If the primitive drift space is extended back to the full Stage-5 list

- `lambda_B`
- `varpi`
- `lambda_U`
- `lambda_W`
- `lambda_R`
- `K`
- `M`
- `Omega_U`
- `Omega_W`
- `delta_U`

then the direct Stage-11 monomial-drift matrix acquires **two exact zero columns** in the `dln lambda_B` and `dln varpi` directions.

So the weak-axisymmetric direct monomials are exactly support-blind:

\[
\partial_{\ln \lambda_B}
(\delta\ln \mathfrak C_{{\rm tr},*},
 \delta\ln \mathfrak C_{{\rm nt},*},
 \delta\ln \epsilon_\eta)
=0,
\]
\[
\partial_{\ln \varpi}
(\delta\ln \mathfrak C_{{\rm tr},*},
 \delta\ln \mathfrak C_{{\rm nt},*},
 \delta\ln \epsilon_\eta)
=0.
\]

That means the extended primitive zero-defect tangent space has dimension

\[
2 + 5 = 7,
\]

namely

1. two BdG-support-blind directions,
2. plus the five normalized similarity directions from Stage 13.

This is important conceptually.

The Stage-10/11 similarity-orbit theorem constrains only the mixed/wall/U normalization problem. It does **not** constrain the explicit BdG support drifts. So monomial-rigidity alone cannot finish the full conservative 5PN problem.

Those BdG directions must still be controlled, if at all, by the separate conservative front-end conditions:

- `K1 = 0`,
- the hidden-even consistency slot,
- and the `O(omega^4)` single-pole / grouped-response test.

So the Stage-5/6 conservative extraction theorem and the Stage-10/11 similarity-orbit theorem are complementary rather than redundant.
# 5PN continuation notes — stages 15 through 17

These stages finally splice the Stage-5 primitive obstruction triplet into the later Stage-10/11 monomial-orbit package.

The main point is that the three objects are **not** on the same footing:

- `Xi_load = N01/N0 - D01/D0 = P1/P0` is exactly the weak-axisymmetric normalization defect governed by the monomial/similarity theorem;
- `K1 = D21 + D01/9` and
- `H_even = D41 - (2/3) D21 - D01/27`

are separate conservative even gates that survive after the normalization defect is killed.

## Stage 15 — the exact `Xi_load` bridge

The exact Stage-5 prefactor identity is
\[
P_0 = \frac{N_0}{D_0},
\qquad
P_1 = \frac{N_{01}D_0 - N_0 D_{01}}{D_0^2},
\qquad
\Xi_{\rm load}
=
\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}
=
\frac{P_1}{P_0}.
\]
So the Stage-5 load defect is **exactly** the same weak-axisymmetric prefactor slope already isolated later as
\[
\Xi_1 = \frac{P_1}{P_0}.
\]

Rewriting Stages 10–14 in one place, the direct monomial drifts are
\[
\Sigma_{\rm tr}=\delta\ln \mathfrak C_{{\rm tr},*},
\qquad
\Sigma_{\rm nt}=\delta\ln \mathfrak C_{{\rm nt},*},
\qquad
\Sigma_\eta=\delta\ln \epsilon_\eta,
\]
and the weak-axisymmetric normalization defect has the exact triangular form
\[
\Xi = A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt}.
\]
The Stage-13 normalized kernel annihilates all three monomial drifts, and the Stage-14 extension shows the explicit BdG directions `dln lambda_B` and `dln varpi` are zero columns of the monomial matrix. Therefore
\[
\Xi_{\rm load}=0
\]
on:

- the full normalized similarity kernel,
- its injected support-blind extension,
- the explicit BdG amplitude/frequency directions,
- the matched wall-only co-scaling direction.

So the monomial/similarity theorem really is a theorem about `Xi_load`.

## Stage 16 — why that still does **not** solve 5PN

The next check is to compare `Xi` with the conservative even gates.

On the matched wall-only co-scaling direction,
\[
d\ln K = d\ln M = 1,
\]
we still get
\[
K1_{\rm wall}=\frac{K}{9}-M,
\qquad
H_{{\rm even,wall}}=-\frac{K}{27}+\frac{2M}{3},
\]
which are generically nonzero even though `Xi = 0`.

On the explicit support-blind BdG directions,
\[
d\ln \lambda_B = 1,
\qquad d\ln \varpi = 0,
\]
or
\[
d\ln \lambda_B = 0,
\qquad d\ln \varpi = 1,
\]
we again have `Xi = 0`, but both `K1_B` and `H_even,B` are nonzero.

Even on the pure BdG self-similar branch
\[
d\ln \lambda_B = d\ln \varpi = \delta,
\]
which is the branch that kills the BdG load defect, one still gets
\[
K1_B = \frac{2B_0\delta}{\varpi^2},
\qquad
H_{{\rm even},B} = \frac{4B_0\delta(3-\varpi^2)}{3\varpi^4},
\]
so the even gates remain open unless `\delta = 0` (or an extra tuning is imposed).

That is the cleanest executable statement yet that:
\[
\text{similarity-orbit rigidity} \neq \text{full 5PN closure}.
\]
The orbit theorem kills `Xi`; it does **not** automatically close `K1` or `H_even`.

## Stage 17 — the first lower-bound intersection calculation

Once the Stage-13/14 monomial-rigid kernel is parameterized by
\[
(\alpha_K,\alpha_W,\alpha_U,\alpha_R,\alpha_{\Omega_U},\beta_B,\beta_{\varpi}),
\]
we can impose the **lower-bound** conservative even gates obtained by keeping only the explicit wall-only and pure-BdG pieces of `K1` and `H_even`.

This gives an exact `2 x 7` gate matrix of generic rank `2`, so the lower-bound even-gate intersection has dimension `5`.

A convenient exact solve is to determine
\[
\alpha_K,
\qquad
\beta_B,
\]
in terms of the remaining five free directions.

The corresponding null basis has five directions:

1. pure `\alpha_W`,
2. pure `\alpha_R`,
3. matched `\alpha_U/\alpha_{\Omega_U}`,
4. BdG-amplitude deformation `\beta_B` with compensating `\alpha_K,\alpha_U`,
5. BdG-frequency deformation `\beta_{\varpi}` with compensating `\alpha_K,\alpha_U`.

The important caution is that this is only a **lower-bound** solve, because the conservative mixed-sector `Z_2,Z_4` moments have not yet been reinstated. In particular, the fact that `\alpha_W` and `\alpha_R` survive untouched here is telling us exactly where the omitted mixed-sector moments still have to act.

## Net result after Stage 17

The 5PN continuation point is now much sharper:

1. `Xi_load = P1/P0` is fully absorbed into the similarity-orbit / monomial-rigidity theorem;
2. the real conservative 5PN bottleneck is the pair of even gates `K1` and `H_even`;
3. the explicit wall-only and BdG-only pieces show those even gates survive after `Xi` is killed;
4. and the next honest task is therefore:

> reinstate the conservative mixed-sector `Z_2,Z_4` moments and determine how they cut the remaining lower-bound tangent family.

That is the cleanest next theorem gate we have reached so far.
# 5PN continuation notes — stages 18 through 20

These stages do two things in parallel.

First, they sharpen the **isotropic full-bundle target surface** that the actual
moving-throat PDE must hit if the calibrated 5PN / 2.5PN / 4PN branch is real.
Second, they reinstate the omitted conservative Maxwell/mixed `Z_2,Z_4` sector in
the weak-axisymmetric even-gate problem and show exactly how it removes the fake
freedom left in the Stage-17 lower-bound picture.

## Stage 18 — exact isotropic full-bundle target surface

Work on the isotropic grouped-lane branch with

a) conservative operator moments

a) `D0 = K - B0 - Z0`,

b) `D2 = -(M + B2 + Z2)`,

c) `D4 = -(B4 + Z4)`.

Then the normalized grouped-response and outgoing-prefactor moments are

`u2 = -D2 / D0`,

`u4 = (D2^2 - D0 D4) / D0^2`,

`P0 = N0 / D0`,

`P2 = (D0 N2 - 2 D2 N0) / D0^2`,

`P4 = (D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0) / D0^3`.

The exact isotropic 5PN conservative one-pole defect is

`u4 - 4 u2^2 = [ D0 (B4 + Z4) - 3 (M + B2 + Z2)^2 ] / D0^2`.

So the isotropic one-pole surface is exactly

`D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`.

The constant-prefactor outgoing branch is also exact:

`P2 = 0  ->  N2 = 2 D2 N0 / D0`,

`P4 = 0  ->  N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

Finally, the universal 2.5PN / 4PN normalization target is

`mhat0^2 P0 = 54 G c_s^5 / (5 a^5 c^5)`

or equivalently

`mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`.

So the full isotropic moving-throat bundle must land on one exact combined target
surface:

1. `D0 = K - B0 - Z0`,
2. `D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`,
3. `mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`,
4. `N2 = 2 D2 N0 / D0`,
5. `N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

That is the sharpest isotropic full-bundle statement we have so far.

## Stage 19 — exact `Z`-sector bridge back into the even gates

The missing conservative Maxwell/mixed moments are

`Z0 = Q / Delta`,

`Z2 = (Q S2 - H Delta) / Delta^2`,

`Z4 = [ Q (S2^2 - Delta) - S2 H Delta ] / Delta^3`.

Their exact first variations are

`dZ0 = (Delta dQ - Q dDelta) / Delta^2`,

`dZ2 = [ Delta (-Delta dH - H dDelta + Q dS2 + S2 dQ) + 2 dDelta (Delta H - Q S2) ] / Delta^3`,

`dZ4 = -[ Delta^2 H dS2 + Delta^2 S2 dH + Delta^2 dQ - 2 Delta H S2 dDelta`
`         - 2 Delta Q S2 dS2 - 2 Delta Q dDelta - Delta S2^2 dQ + 3 Q S2^2 dDelta ] / Delta^4`.

Therefore the conservative `Z`-sector contributions to the two live 5PN even gates are

`K1_Z = -dZ2 - dZ0/9`,

`H_even,Z = -dZ4 + (2/3) dZ2 + dZ0/27`.

After inserting the exact Stage-13 normalized similarity kernel, both become honest
functions of the mixed-sector similarity directions `alpha_W` and `alpha_R`.

On the positive constructive slice

- `G_U = 5`, `G_W = 7`, `R = 2`,
- `Omega_U = 3`, `Omega_W = 4`,
- `chi_0 = 3/2`, `delta_U = 2/3`,
- `E_* = 1/4`, `F_* = 5/6`,

we get

`K1_Z = (78623501/25004700) alpha_OmegaU + (733046/6251175) alpha_R`
`       - (59010631/25004700) alpha_U - (32134513/50009400) alpha_W`,

`H_even,Z = -(28906377971/21003948000) alpha_OmegaU`
`           - (1174937411/21003948000) alpha_R`
`           + (11102468471/10501974000) alpha_U`
`           + (5617869293/21003948000) alpha_W`.

In particular, the pure directions give

`alpha_W:  K1_Z = -32134513/50009400,   H_even,Z = 5617869293/21003948000`,

`alpha_R:  K1_Z = 733046/6251175,       H_even,Z = -1174937411/21003948000`.

So the omitted `Z_2,Z_4` block does exactly what Stage 17 said it still had to do:
it activates the previously untouched mixed directions.

## Stage 20 — full even-gate solve on the constructive branch

Once wall-only, pure-BdG, and the reinstated `Z`-sector are combined, the full
constructive-slice even-gate matrix is

`Gate_full(slice) =`

`[ -25/9,  -32134513/50009400,   91017569/25004700,    733046/6251175,`
`  -71404699/25004700,  -8/9,  4/3 ]`

`[  52/27,  5617869293/21003948000, -30905427529/10501974000,`
`  -1174937411/21003948000, 55109414029/21003948000, 32/81, -16/27 ]`.

The matrix still has rank `2`, so the full even-gate intersection is still
5-dimensional. That part is unsurprising: there are still only two even equations.

The new structural fact is the mixed-sector minor:

`det Gate_(alpha_W, alpha_R) = 942737330573 / 205838690400000 != 0`.

So on this branch the full even system can solve directly for the previously
untouched mixed directions:

`alpha_W =`
`  (14503089433000/942737330573) alpha_K`
`+ (30450672110098/942737330573) alpha_OmegaU`
`- (29120459867142/942737330573) alpha_U`
`- (18876066395200/25453907925471) beta_B`
`+ (9438033197600/8484635975157) beta_varpi`,

`alpha_R =`
`  (101802968743000/942737330573) alpha_K`
`+ (189815725996721/942737330573) alpha_OmegaU`
`- (188832473718440/942737330573) alpha_U`
`+ (89510801038400/25453907925471) beta_B`
`- (44755400519200/8484635975157) beta_varpi`.

So there are no longer pure `alpha_W` or pure `alpha_R` null directions on the full
constructive branch. The mixed-sector freedom was never genuinely unconstrained; it
was just hidden in the omitted `Z_2,Z_4` block.

## Net result after Stage 20

The continuation point is sharper again.

1. The isotropic full-bundle target surface is now exact and explicit.
2. The omitted conservative `Z_2,Z_4` sector has been reinstated exactly.
3. On a clean constructive branch it does the precise job Stage 17 predicted:
   it removes the fake freedom in the mixed directions `alpha_W, alpha_R`.
4. So the remaining work is no longer “some mixed-sector freedom somewhere.”
   It is to compute the actual overlap data on the true moving-throat branch and see
   whether they land on the isotropic full-bundle target surface from Stage 18.
# 5PN continuation notes — Stages 21–23 and 27–29

This note fills the earlier numbering gap and records the next continuation past the Stage 26 placement map.

## Missing stages 21–23

### Stage 21 — Dimensionless continuum placement map
The continuum-selected quadrupole placement problem compresses to the exact dimensionless ratios

a) `eps_eta = c_(etaU)^2 / (K_U K_eta^(eff))`

b) `eps_W = c_(UW)^2 sigma / (K_U K_W^(eff))`

c) `rho = c_(UW) c_(etaU) / (K_U c_(etaW))`

d) `Z_W = c_(etaW)^2 / (K_eta^(eff) K_W^(eff))`

e) `delta_0 = pi^2 T_w / (L^2 K_eta^(eff))`

with radiative scale

`Lambda = 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W)`.

The exact placement coordinates are

`delta = delta_0 / (1 - eps_eta)`,

`M_mix = 8 Z_W (1 + rho)^2 / [pi^2 (1 - eps_eta)(1 - eps_W)]`,

`R_target = Lambda (1 - eps_eta)(1 - eps_W)^2 / [Z_W (1 + rho)^2]`.

The key exact product law is

`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`

`= 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W)`.

So Stage 21 factorizes the continuum map into a geometry lane, a mixed-stability/product lane, and a redistribution lane.

### Stage 22 — First non-flat U doublet
Turning on the first axial structure of the internal `U` sector preserves the scalar placement map but breaks the old one-direction geometry.

The new exact quantities are

`delta_split = [delta_0 + eps_eta delta_U/(1+delta_U)] / (1-eps_eta)`,

`eps_W_split = eps_W [1 - (2/11) delta_U/(1+delta_U)]`,

`R_U = [1 + rho_0/(1+delta_U)] / (1 + rho_0)`.

The exact direction-splitting invariant is

`D_dir = - kappa_0 kappa_1 g_W rho_0 delta_U / (1+delta_U)`.

So collinearity survives iff `delta_U = 0` or `rho_0 = 0`.

### Stage 23 — Generalized selected-branch normalization
With distinct loading vector `z` and source vector `s`, the selected branch survives exactly, but the old flat-`U` function pair `(F,G)` deforms to

`F_(q,eta)(xi,delta)` and `G_q(xi,delta)`.

For the split-`U` continuum, the whole deformation collapses to the single ratio `R_U`, giving the exact one-parameter family

`F_U(xi,delta;R_U)`, `G_U(xi,delta;R_U)`.

Setting `R_U = 1` recovers the Stage 18–19 flat branch exactly.

The exact first deformation about the flat branch is

`F_U/F_flat = 1 + eps H_F + O(eps^2)`,

`G_U/G_flat = 1 + eps H_G + O(eps^2)`,

with the closed-form `H_F` and `H_G` verified in the companion SymPy audit.

## Continuation past Stage 26: Stages 27–29

### Stage 27 — Continuum-selected rank-2 closure
Once the actual support-direction data from Stage 26 are inserted, the physical selected wall branch is pinned by an exact quadratic for the softening depth `xi`:

`xi^2 + B_cont xi + C_cont = 0`,

with

`B_cont = delta - M_mix (1+t^2 R_U^2) - M_supp (1+t^2 R_phi^2)`,

`C_cont = - delta (M_mix + M_supp) + t^2 M_mix M_supp (R_U - R_phi)^2`.

So the continuum-selected normalization theorem gate becomes simply

`R_target = F_cont(xi_phys)`.

The exact special surfaces are:

- minimal-kernel source-tied surface: `sigma_0 = 0`,
- interference-matched tracking surface: `sigma_0 = rho_0`.

### Stage 28 — Coherent local D/N support kernel
The first coherent local D/N kernel forces

`g_B g_R = g_W g_S`,

so it lands exactly on the tracking surface. Therefore

`R_phi = R_U = R_tr`,

and the full Stage 27 rank-2 closure collapses to the one-parameter tracking laws

`M_tr = G_tr(xi,delta;R_tr)`,

`R_target = F_tr(xi,delta;R_tr)`.

The coherent branch also has the exact constructive bound

`1/(1+delta_U) < R_tr < 1`.

### Stage 29 — Tracking-branch bounds and residual comparison
On the coherent tracking branch, at fixed `(xi,delta)`,

`dG_tr/dR < 0`,

`dF_tr/dR > 0`.

So for the constructive split-`U` branch, where `R_tr < 1`, one has the exact inequalities

`G_tr > G_flat`,

`F_tr < F_flat`.

The endpoint formulas are

`G_tr(R=0) = xi`,

`F_tr(R=0) = 1/(1-xi)`.

The exact residual comparison theorem is

`E_tr - E_flat = F_flat - F_tr > 0`.

So the first coherent local D/N kernel does **not** rescue the Stage 18–19 normalization target. It makes the target harder: more total loading is required, while the normalized selected-branch response is lower than on the flat branch.

## Bottom line after this batch

The selected-branch problem is now narrower than before:

1. the physical softening depth is set by an exact quadratic;
2. the first coherent local D/N kernel forces exact tracking rather than a generic intermediate rank-2 closure;
3. and the resulting split-`U` deformation worsens the normalization residual relative to the old flat branch.

So the next theorem gate is no longer “which closure should we use?” It is whether the later support-enhancement / compensation mechanisms in the notes are strong enough to overcome this exact tracking-branch deficit.

## Further continuation — Stages 30–31

### Stage 30 — Coherent-kernel dimensionless map
The first coherent local D/N kernel compresses to a very small exact parameter set:

- `eps_eta`
- `delta_U`
- `chi_0`
- `eps_W`
- `Z_W`
- `zeta`
- `Lambda`

The exact coherent tracking factor is

`R_tr = [1 + chi_0/(1+delta_U)] / (1 + chi_0)`.

The total baseline splits as

`M_mix = 8 Z_W (1 + chi_0)^2 / [pi^2 (1 - eps_eta)(1 - eps)]`,

`M_supp = 8 zeta Z_W (1 + chi_0)^2 / [pi^2 (1 - eps_eta)(1 - zeta eps)]`,

with

`eps = eps_W [1 - (2/11) delta_U/(1+delta_U)]`.

So the support lane enters only through the exact enhancement factor

`S(zeta;eps) = 1 + zeta(1-eps)/(1-zeta eps)`,

and

`M_tr = M_mix S(zeta;eps)`.

The target ratio remains independent of `zeta`:

`R_target = Lambda (1 - eps_eta)(1 - eps)^2 / [Z_W (1 + chi_0)^2]`.

### Stage 31 — Support-compensation theorem
On the coherent tracking branch,

`M_tr = G_tr(xi,delta;R_tr)`

with critical load

`M_crit(delta,R_tr) = G_tr(1,delta;R_tr)`.

The support enhancement factor is strictly increasing and exactly invertible:

`zeta_req = (S_req - 1) / [1 + eps (S_req - 2)]`.

For every finite target ratio `R_target > 1`, there exists a stable-side branch point `xi_req in (0,1)` and a corresponding required load `M_req < M_crit`.

Therefore, if the mixed-only branch is below that required load, there is a unique coherent support ratio `zeta_req < zeta_crit < 1/eps` that reaches the target **before** the branch softens out.

So there is no reduced-level support no-go. The remaining question is whether the actual PDE produces a physical `zeta` large enough to meet that exact threshold.

## Further continuation — Stages 32–33

### Stage 32 — Explicit D/N overlap extraction of `zeta`
The coherent support ratio is no longer phenomenological. For the first explicit finite-throat D/N family,

`zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) (I_n / I_0)^2`,

with

`I_n = ∫_0^L ds sigma(s) chi_n(s)`,

`chi_n(s) = sqrt(2/L) sin[(n+1/2) pi s / L]`.

For the first uniform local source density `sigma(s)=1`,

`I_n / I_0 = 1 / (2n+1)`.

So the physical coherent support ratio becomes

`zeta_n^(phys) = (K_W^(eff) / K_(phi,n)^(eff)) / (2n+1)^2`.

On the same-operator twin family,

`zeta_n^(twin) = 1 / [ (2n+1)^2 (1 + x n(n+1)) ]`.

The exact lowest-twin value is

`zeta_0^(twin) = 1`.

### Stage 33 — Exact comparison of `zeta_phys` against `zeta_req`
Because Stage 31 reduced support feasibility to `zeta_phys >= zeta_req`, the explicit D/N tower gives exact selection rules.

The lowest symmetric twin lane has

`zeta_0^(twin) = 1`,

so its enhancement factor is exactly

`S_0 = 2`,

independent of `eps`.

Therefore the lowest symmetric twin lane succeeds exactly when

`zeta_req <= 1`,

equivalently

`S_req <= 2`.

For higher D/N harmonics,

`zeta_n^(twin) < 1/(2n+1)^2`,

so they are ruled out immediately if

`zeta_req > 1/(2n+1)^2`.

When they are not ruled out immediately, the exact softness threshold is

`x <= x_max(n; zeta_req)`

with

`x_max(n; zeta_req) = [1 / ((2n+1)^2 zeta_req) - 1] / [n(n+1)]`.

So the explicit coherent support tower is strongly biased toward the lowest twin lane.

## Bottom line after Stages 21–33

The program is now much sharper than when the numbering gap first appeared.

1. Stages 21–23 show that split-`U` physics preserves the scalar placement map but deforms the selected-branch normalization geometry through an exact direction-splitting parameter.
2. Stages 27–29 show that once the actual continuum support data are inserted, the physical wall softening depth is fixed by an exact quadratic, and the first coherent local D/N kernel lands on a tracking branch whose split-`U` deformation worsens the normalization target relative to the old flat branch.
3. Stages 30–31 show there is no reduced-level support no-go: coherent support enhancement can in principle compensate the tracking-branch deficit before softening.
4. Stages 32–33 then make that compensation test microscopic: the lowest symmetric twin D/N lane is an exact universal doubling branch, and every higher support harmonic is strongly overlap- and stiffness-suppressed.

So the next really important question is no longer which closure to use. It is whether the actual moving-throat PDE places the physical coherent support sector on the lowest twin branch with a large enough `zeta_phys` to satisfy the exact threshold.
# 5PN / Moving-Throat continuation — Stages 24–26

These stages continue the explicit moving-throat spectral-placement branch after the concrete axial / loaded-profile / selected-prefactor scripts.

## Stage 24 — source map, microscopic normalization, softening depth

Files:
- `5pn_stage24_source_map_microscopic_normalization_softening.py`
- `5pn_stage24_source_map_microscopic_normalization_softening_output.txt`

Main results:

1. The natural D/N source-map factor is no longer independent:
   
   \[
   \hat m_-^2 = \frac{s_-(x)}{\kappa_0^2}.
   \]

2. The selected-branch normalization product becomes
   
   \[
   N_-(x) = \frac{\beta_0\, s_-(x)^2}{\kappa_0^2 (A-x)}.
   \]

3. The selected branch is exactly parameterized by the softening depth
   \[
   x = A-\lambda_-.
   \]

4. The exact secular loading law is
   
   \[
   \alpha_0(x)=\frac{x(x+\Delta K_{ax})}{\kappa_0^2(x+\Delta K_{ax})+\kappa_1^2 x}.
   \]

5. The exact support-loading requirement is
   
   \[
   \frac{g_{B,req}^2}{\varpi^2}=\alpha_0(x)-\alpha_{mix}.
   \]

So the selected quadrupole bridge is no longer an eigenvector problem. It is a one-variable spectral-placement problem in `x`.

---

## Stage 25 — universal D/N branch geometry

Files:
- `5pn_stage25_dimensionless_normalization_and_support_frontier.py`
- `5pn_stage25_dimensionless_normalization_and_support_frontier_output.txt`

Main results:

Using the exact D/N constants
\[
\kappa_0^2=\frac{8}{\pi^2},\qquad
\kappa_1^2=\frac{16}{9\pi^2},\qquad
\eta=\frac{\kappa_1^2}{\kappa_0^2}=\frac29,
\]
and the dimensionless variables
\[
\xi=\frac{x}{A},\qquad \delta=\frac{\Delta K_{ax}}{A},
\]
we get two universal branch functions:

1. The normalization function
   
   \[
   F(\xi,\delta)=
   \frac{(9\delta+11\xi)^4}{81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2}.
   \]

2. The support-feasibility function
   
   \[
   G(\xi,\delta)=\frac{9\xi(\xi+\delta)}{9\delta+11\xi}.
   \]

Exact geometric facts verified in SymPy:

- `F` is strictly increasing on the stable branch,
- `F(0,delta)=1`,
- `F -> +infinity` as `xi -> 1^-`,
- `G` is strictly increasing,
- `G(0,delta)=0`,
- `G_max(delta)=9(1+delta)/(9 delta + 11)`.

So for fixed `delta`, the selected branch is controlled by a unique normalization locus `F = R_target` and a sharp support-feasibility frontier `M_mix <= G`.

---

## Stage 26 — continuum-kernel extraction and placement map

Files:
- `5pn_stage26_continuum_kernel_extraction_and_placement_map.py`
- `5pn_stage26_continuum_kernel_extraction_and_placement_map_output.txt`

Main results:

From the first explicit finite-throat continuum kernel, the reduced branch data are extracted exactly:

\[
A=\frac{K_U K_{\eta}^{eff}-c_{\eta U}^2}{\mu_\eta K_U},
\qquad
\Delta K_{ax}=\frac{\pi^2 T_w}{\mu_\eta L^2},
\]
\[
\alpha_{mix}=
\frac{(K_U c_{\eta W}+c_{UW} c_{\eta U})^2}
{\mu_\eta K_U (K_U K_W^{eff}-c_{UW}^2\sigma)},
\]
\[
\beta_0=
\frac{\mu_W}{\mu_\eta}
\frac{(K_U c_{\eta W}+c_{UW} c_{\eta U})^2}
{(K_U K_W^{eff}-c_{UW}^2\sigma)^2}.
\]

Then the full placement problem compresses to the exact dimensionless kernel ratios

\[
\epsilon_\eta = \frac{c_{\eta U}^2}{K_U K_\eta^{eff}},
\qquad
\epsilon_W = \frac{c_{UW}^2\sigma}{K_U K_W^{eff}},
\]
\[
\rho = \frac{c_{UW} c_{\eta U}}{K_U c_{\eta W}},
\qquad
Z_W = \frac{c_{\eta W}^2}{K_\eta^{eff} K_W^{eff}},
\qquad
\delta_0 = \frac{\pi^2 T_w}{L^2 K_\eta^{eff}},
\]
plus the radiative demand scale
\[
\Lambda = \frac{27 \pi^2 G c_s^5 K_W^{eff}}{20 a^5 c^5 \mu_W}.
\]

The exact placement formulas are

\[
\delta = \frac{\delta_0}{1-\epsilon_\eta},
\qquad
M_{mix} = \frac{8 Z_W (1+\rho)^2}{\pi^2 (1-\epsilon_\eta)(1-\epsilon_W)},
\]
\[
R_{target} = \frac{\Lambda (1-\epsilon_\eta)(1-\epsilon_W)^2}{Z_W (1+\rho)^2}.
\]

And the exact product law is

\[
R_{target} M_{mix} = \frac{8\Lambda(1-\epsilon_W)}{\pi^2}.
\]

This cleanly separates the placement problem into three lanes:

1. geometry lane: `delta = delta_0 / (1-eps_eta)`
2. mixed-stability/product lane: `R_target M_mix = 8 Lambda (1-eps_W)/pi^2`
3. redistribution lane: `(eps_eta, Z_W, rho)` move the branch along the fixed product curve.

---

## What these stages change in the 5PN program

The remaining bridge is no longer a vague “solve more PDE” task.
It is now a very narrow branch-placement problem:

1. compute the actual continuum ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)`,
2. map them to `(delta, M_mix, R_target)`,
3. find the unique `xi` solving `F(xi,delta)=R_target`,
4. and check the support-feasibility condition `M_mix <= G(xi,delta)`.

If that succeeds, the selected quadrupole branch is admissible on the natural finite-throat continuum placement map.
# 5PN stages 34–40 notes

This bundle continues the post-Stage-33 support/normalization program by turning the
`zeta_req` threshold into exact operator criteria on the moving-throat branch.

## Stage 34 — exact lowest-twin sufficiency criterion

Using the tracking-branch functions
\[
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]
\[
F_{\rm tr}(\xi,\delta;R)
=
\frac{\bigl[9\delta+(9+2R^2)\xi\bigr]^2\bigl[9\delta+(9+2R)\xi\bigr]^2}
{81(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2},
\]
the exact product is
\[
\Pi_{\rm tr}=F_{\rm tr}G_{\rm tr}
=
\frac{\xi(\xi+\delta)\bigl[9\delta+(9+2R)\xi\bigr]^2\bigl[9\delta+(9+2R^2)\xi\bigr]}
{9(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2}.
\]

With
\[
C_{\rm mix}:=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\qquad
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}},
\]
the symmetric lowest twin lane is sufficient iff
\[
\Pi_{\rm tr}(\xi_{\rm req},\delta;R_{\rm tr})
\le
\frac{16\Lambda(1-\epsilon)}{\pi^2}.
\]

Equivalent threshold scales:
\[
\Lambda_{\rm twin,req}=\frac{\pi^2}{16(1-\epsilon)}\Pi_{\rm tr},
\qquad
M_{\rm mix}^{(\rm twin,req)}=\frac{G_{\rm tr}}{2},
\]
\[
Z_W^{(\rm twin,req)}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)\,G_{\rm tr}}
{16(1+\chi_0)^2}.
\]

The exact twin-saturation depth at fixed mixed baseline is the unique root of
\[
G_{\rm tr}(\xi_{2\times},\delta;R)=2M_{\rm mix},
\]
namely
\[
\xi_{2\times}
=
\frac{2M_{\rm mix}(9+2R^2)-9\delta+\sqrt{(2M_{\rm mix}(9+2R^2)-9\delta)^2+648M_{\rm mix}\delta}}{18}.
\]

## Stage 35 — exact non-twin asymmetry requirement

Define
\[
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}
=
\frac{\Pi_{\rm tr}-C_{\rm mix}}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]
Then
\[
\frac{d\zeta_{\rm req}}{d\Pi_{\rm tr}}
=
\frac{C_{\rm mix}(1-\epsilon)}{\bigl[C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})\bigr]^2}>0,
\]
so the required coherent support ratio grows monotonically with the branch product.

Exact regime split:
\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\Rightarrow\quad
\text{mixed-only enough},
\]
\[
C_{\rm mix}<\Pi_{\rm tr}\le2C_{\rm mix}
\quad\Rightarrow\quad
\text{symmetric lowest twin enough},
\]
\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\Rightarrow\quad
\text{non-twin asymmetry required}.
\]

The exact excess beyond the symmetric twin branch is
\[
\Delta_\zeta:=\zeta_{\rm req}-1
=
\frac{(1-\epsilon)(\Pi_{\rm tr}-2C_{\rm mix})}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]

For a general lowest support lane,
\[
\zeta_0^{(\rm phys)}
=
\frac{K_W^{(\rm eff)}}{K_{\phi,0}^{(\rm eff)}}\Omega_0^2.
\]
So the two equivalent exact rescue thresholds are
\[
\Omega_0^2 \ge \zeta_{\rm req}\frac{K_{\phi,0}^{(\rm eff)}}{K_W^{(\rm eff)}},
\qquad
K_{\phi,0}^{(\rm eff)} \le K_W^{(\rm eff)}\frac{\Omega_0^2}{\zeta_{\rm req}}.
\]

## Stage 36 — exact overlap-boost window

For the D/N lowest support mode
\[
\chi_0(s)=\sqrt{\frac{2}{L}}\sin\!\frac{\pi s}{2L},
\qquad
I_W=\int_0^L\chi_0(s)\,ds=\frac{2\sqrt{2L}}{\pi},
\]
and the normalized exponential source family
\[
\sigma_\alpha(s)=\frac{\alpha e^{\alpha s/L}}{e^\alpha-1},
\qquad
\int_0^L \sigma_\alpha(s)\,ds=L,
\]
the overlap boost is
\[
\Omega_{\exp}(\alpha)
=
\frac{\int_0^L \sigma_\alpha(s)\chi_0(s)\,ds}{I_W}
=
\frac{\pi\alpha\bigl(2\alpha e^\alpha+\pi\bigr)}
{(4\alpha^2+\pi^2)(e^\alpha-1)}.
\]

Exact endpoint values:
\[
\Omega_{\exp}(0)=1,
\qquad
\lim_{\alpha\to+\infty}\Omega_{\exp}(\alpha)=\frac{\pi}{2}.
\]
Therefore
\[
0\le \Omega_0\le\frac{\pi}{2},
\qquad
A_I:=\Omega_0^2\le\frac{\pi^2}{4}.
\]

So pure overlap rescue alone is possible only if
\[
\zeta_{\rm req}\le\frac{\pi^2}{4}.
\]

## Stage 37 — Robin-compliance softening

Replacing the Dirichlet mouth by a Robin compliance gives the lowest-lane eigenvalue condition
\[
y\tan y=\eta,
\qquad
0<y<\frac{\pi}{2},
\]
with
\[
\eta:=hL.
\]

If
\[
x:=\frac{\pi^2 T_X}{L^2K_W^{(\rm eff)}},
\qquad
0<x<4,
\]
then the exact support-softening factor is
\[
A_K(\eta)
=
\frac{K_W^{(\rm eff)}}{K_{\phi,0}^{(\rm eff)}}
=
\frac{1}{1-x/4+xy^2/\pi^2}.
\]

Endpoint window:
\[
A_K=\;1\;\text{at}\;y=\frac{\pi}{2},
\qquad
A_{K,\max}=\frac{4}{4-x}\;\text{at}\;y\to0.
\]

So pure support softening alone can rescue the Stage-35 threshold only if
\[
\zeta_{\rm req}\le\frac{4}{4-x}.
\]

At fixed \(\zeta_{\rm req}\), the exact eigenvalue and Robin thresholds are
\[
y_{\rm req}^2=\frac{\pi^2}{x}\left(\frac{1}{\zeta_{\rm req}}-1+\frac{x}{4}\right),
\qquad
\eta_{\rm req}=y_{\rm req}\tan y_{\rm req}.
\]

## Stage 38 — explicit non-twin lowest-lane reachability

Combining the Stage-36 overlap family with Stage-37 Robin softening gives
\[
\zeta_0^{(\exp+R)}(\alpha,\eta)
=
\frac{\Omega_{\exp}(\alpha)^2}{1-x/4+xy(\eta)^2/\pi^2}.
\]

Its exact closure range is
\[
1\le \zeta_0^{(\exp+R)} \le \frac{\pi^2}{4-x}.
\]

So the explicit family reaches the Stage-35 threshold iff
\[
\zeta_{\rm req}\le\frac{\pi^2}{4-x}.
\]

This gives the exact three-regime split:
\[
\zeta_{\rm req}\le\frac{\pi^2}{4}
\quad\Rightarrow\quad
\text{overlap alone enough},
\]
\[
\frac{\pi^2}{4}<\zeta_{\rm req}\le\frac{\pi^2}{4-x}
\quad\Rightarrow\quad
\text{overlap + softening enough},
\]
\[
\zeta_{\rm req}>\frac{\pi^2}{4-x}
\quad\Rightarrow\quad
\text{even this explicit family fails}.
\]

## Stage 39 — transport origin of source asymmetry

The Stage-36 exponential family is exactly the stationary zero-flux branch of
\[
\partial_t \sigma + \partial_s J = 0,
\qquad
J=-D_\sigma \partial_s \sigma + v_\sigma \sigma.
\]

On the stationary recirculating branch \(J=0\),
\[
\sigma_{Pe}(s)
=
\frac{Pe\,e^{Pe s/L}}{e^{Pe}-1},
\qquad
Pe:=\frac{v_\sigma L}{D_\sigma}.
\]

The corresponding overlap boost is
\[
\Omega_{Pe}
=
\frac{\pi Pe\bigl(2Pe\,e^{Pe}+\pi\bigr)}
{(4Pe^2+\pi^2)(e^{Pe}-1)}.
\]

Exact identities:
\[
\Omega_{Pe}(0)=1,
\qquad
\lim_{Pe\to+\infty}\Omega_{Pe}=\frac{\pi}{2},
\]
and the score-function identity
\[
\partial_{Pe}p_{Pe}(x)=(x-\mathbb E_{Pe}[x])p_{Pe}(x)
\]
implies
\[
\frac{d}{dPe}\mathbb E_{Pe}[\chi_0]=\operatorname{Cov}_{Pe}(\chi_0,x),
\]
so the constructive branch is monotone increasing because \(\chi_0\) is increasing on \([0,1]\).

## Stage 40 — physical \((Pe,\kappa,\eta)\) placement map

Define the physical support ratios
\[
\kappa:=\frac{K_XL^2}{T_X},
\qquad
\eta:=hL=\frac{K_mL}{T_X}.
\]
Then
\[
x=\frac{\pi^2}{\kappa+\pi^2/4},
\]
and the Robin softening factor becomes
\[
A_K(\eta;\kappa)=\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2}.
\]

The explicit physical lowest-lane family is therefore
\[
\zeta_0^{(Pe+R)}(Pe,\eta;\kappa)
=
\Omega_{Pe}^2\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2}.
\]

This map is monotone:
- increasing in \(Pe\),
- decreasing in \(\eta\),
- decreasing in \(\kappa\).

Its exact constructive-branch ceiling is
\[
\zeta_{\max}(\kappa)=\frac{\pi^2}{4}\frac{\kappa+\pi^2/4}{\kappa}.
\]

So the Stage-35 demand is reachable on this first physical family iff
\[
\zeta_{\rm req}\le \zeta_{\max}(\kappa),
\]
equivalently, whenever \(\zeta_{\rm req}>\pi^2/4\),
\[
\kappa \le \kappa_{\max}(\zeta_{\rm req})
:= \frac{\pi^4}{4(4\zeta_{\rm req}-\pi^2)}.
\]

The exact physical threshold surfaces are
\[
\Omega_{\rm req}^2
=
\zeta_{\rm req}\frac{\kappa+y(\eta)^2}{\kappa+\pi^2/4},
\]
\[
y_{\rm req}^2
=
\frac{\Omega_{Pe}^2}{\zeta_{\rm req}}(\kappa+\pi^2/4)-\kappa,
\]
\[
\kappa_{\rm req}
=
\frac{\Omega_{Pe}^2\pi^2/4-\zeta_{\rm req}y(\eta)^2}{\zeta_{\rm req}-\Omega_{Pe}^2}.
\]

## Where the program stands after Stage 40

The support/normalization problem is no longer phrased in abstract deformation variables.
It has collapsed to three physical moving-throat operator ratios:

- axial source Peclet number `Pe`,
- mouth compliance `eta`,
- baseline support stiffness ratio `kappa`.

That makes the next clean move very sharp: derive the coupled support/source branch equation that selects `Pe`
from the same operator that already fixes `eta` and `kappa`.
# 5PN stages 41–51 notes

This bundle continues the post-Stage-40 support/normalization program by turning the
physical placement map into an operator-selected branch law and then projecting the
remaining theorem gap back into explicit parent-action overlap data.

## Stage 41 — coupled support/source operator and exact `Pe` branch equation

The first coupled axial operator is

a) source transport
\[
\partial_t\sigma + \partial_s J = 0,
\qquad
J = -D_\sigma \partial_s\sigma + v_\sigma \sigma,
\]

b) support field
\[
-T_X \partial_s^2 \phi + K_X \phi = \Lambda_\phi \sigma,
\]
with support boundary conditions
\[
T_X \phi_s(0)=K_m\phi(0),
\qquad
\phi_s(L)=0.
\]

After nondimensionalization,
\[
\kappa = K_X L^2/T_X,
\qquad
\eta = K_m L/T_X,
\qquad
Pe = v_\sigma L/D_\sigma,
\qquad
\Xi = \mu_\sigma \Lambda_\phi^2 L^2/(D_\sigma T_X).
\]

On the stationary zero-flux branch,
\[
\Sigma_{Pe}(x)=\frac{Pe\,e^{Pe x}}{e^{Pe}-1},
\qquad x=s/L\in[0,1].
\]
The exact support-drop kernel is
\[
K_{\kappa,\eta}(x)=\frac{\cosh(\alpha x)+(\eta/\alpha)\sinh(\alpha x)-\cosh(\alpha(1-x))}{\alpha\sinh\alpha+\eta\cosh\alpha},
\qquad \alpha=\sqrt\kappa,
\]
with derivative
\[
\frac{dK_{\kappa,\eta}}{dx}=
\frac{\alpha\sinh(\alpha x)+\eta\cosh(\alpha x)+\alpha\sinh(\alpha(1-x))}{\alpha\sinh\alpha+\eta\cosh\alpha}>0.
\]

So the exact dimensionless support drop is
\[
\Delta(Pe;\kappa,\eta)=\int_0^1 dx\;K_{\kappa,\eta}(x)\Sigma_{Pe}(x),
\]
with endpoint values
\[
\Delta_0(\kappa,\eta)=
\frac{\eta(\cosh\alpha-1)}{\alpha^2(\alpha\sinh\alpha+\eta\cosh\alpha)},
\]
\[
\Delta_\infty(\kappa,\eta)=
\frac{\cosh\alpha+(\eta/\alpha)\sinh\alpha-1}{\alpha\sinh\alpha+\eta\cosh\alpha}.
\]

The branch point is therefore selected by the exact fixed-point equation
\[
Pe = \Xi\,\Delta(Pe;\kappa,\eta),
\]
and every constructive branch root obeys
\[
\Xi\Delta_0(\kappa,\eta)
\le Pe_* \le
\Xi\Delta_\infty(\kappa,\eta).
\]
At weak coupling,
\[
Pe_* = \Xi \Delta_0(\kappa,\eta)+O(\Xi^2).
\]

## Stage 42 — exact residual bounds on the operator-selected branch

Evaluating the Stage-40 physical support ratio on the branch root gives
\[
\zeta_{\rm phys}(\Xi,\eta;\kappa)=
\Omega_{Pe_*}^2\,\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2},
\qquad y\tan y=\eta.
\]
Since \(\Omega_{Pe}\) is strictly increasing, the Stage-41 branch interval gives the exact support bracket
\[
\zeta_-(\Xi,\eta;\kappa)
\le \zeta_{\rm phys}(\Xi,\eta;\kappa)
\le \zeta_+(\Xi,\eta;\kappa),
\]
where
\[
\zeta_-=
\Omega_{\Xi\Delta_0}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2},
\qquad
\zeta_+=
\Omega_{\Xi\Delta_\infty}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2}.
\]

Inside the Stage-40 reachability window, define the unique constructive point \(Pe_{\rm req}\) by
\[
\Omega_{Pe_{\rm req}}^2=
\zeta_{\rm req}\,\frac{\kappa+y^2}{\kappa+\pi^2/4}.
\]
Then the exact coupling thresholds are
\[
\Xi_{\rm fail}=\frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)},
\qquad
\Xi_{\rm suff}=\frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\qquad
\Xi_{\rm fail}\le \Xi_{\rm suff}.
\]
So the branch has a sharp three-zone structure:

- \(\Xi\le\Xi_{\rm fail}\): impossible,
- \(\Xi\ge\Xi_{\rm suff}\): guaranteed,
- \(\Xi_{\rm fail}<\Xi<\Xi_{\rm suff}\): only here is the full root solve needed.

The exact residual envelope is
\[
R_- \le R_{\rm phys} \le R_+,
\qquad
R_{\rm phys}=\zeta_{\rm req}-\zeta_{\rm phys},
\]
with
\[
R_- = \zeta_{\rm req}-\zeta_+,
\qquad
R_+ = \zeta_{\rm req}-\zeta_-.
\]

Using
\[
\Omega_{Pe}^2 = 1 + \frac{4-\pi}{\pi}Pe + O(Pe^2),
\]
one gets the weak-coupling branch law
\[
\zeta_{\rm phys}=
A_K(\eta;\kappa)
\Bigl[1+\frac{4-\pi}{\pi}\,\Xi\Delta_0(\kappa,\eta)+O(\Xi^2)\Bigr].
\]

## Stage 43 — entropic source microclosure and microscopic gain

The first explicit source/support free energy is
\[
F[\sigma,\phi]=
\int_0^L ds\Bigl[
\Theta_\sigma\sigma(\log(\sigma/\sigma_*)-1)-\Lambda_\phi\sigma\phi
+\frac{T_X}{2}\phi_s^2+\frac{K_X}{2}\phi^2
\Bigr]+
\frac{K_m}{2}\phi(0)^2.
\]
Its exact variations are
\[
\mu_\sigma^{\rm chem}=
\frac{\delta F}{\delta\sigma}=
\Theta_\sigma\log(\sigma/\sigma_*)-\Lambda_\phi\phi,
\]
\[
-T_X\phi_{ss}+K_X\phi=\Lambda_\phi\sigma,
\qquad
T_X\phi_s(0)=K_m\phi(0),
\qquad
\phi_s(L)=0.
\]

The minimal positive-density Onsager current is
\[
J=-M_\sigma\sigma\partial_s\mu_\sigma^{\rm chem}
  =-D_\sigma\partial_s\sigma + M_\sigma\Lambda_\phi\sigma\partial_s\phi,
\]
with exact Einstein relation
\[
D_\sigma=M_\sigma\Theta_\sigma.
\]

Under the affine-drop reduction
\[
\phi(s)\approx \phi(0)+[\Delta\phi]s/L,
\]
the stationary zero-flux branch gives the exact exponential family
\[
\sigma(s)=C\exp\!igl[(\Lambda_\phi\Delta\phi)/(\Theta_\sigma L)\,s\bigr],
\]
so
\[
Pe=(\Lambda_\phi/\Theta_\sigma)\Delta\phi.
\]
Using
\[
\Delta\phi=(\Lambda_\phi L^2/T_X)\Delta(Pe;\kappa,\eta),
\]
one gets the exact microscopic coupling
\[
\Xi_{\rm micro}=
\frac{\Lambda_\phi^2 L^2}{\Theta_\sigma T_X}
=
\chi_\sigma\frac{\Lambda_\phi^2 L^2}{T_X},
\qquad
\chi_\sigma:=1/\Theta_\sigma.
\]

The closure is automatically passive:
\[
\frac{dF}{dt}=-\int_0^L ds\;\frac{J^2}{M_\sigma\sigma}\le 0
\]
under no-flux boundaries.

## Stage 44 — microscopic gain thresholds and exact phase diagram

Using \(\kappa=K_XL^2/T_X\), the Stage-43 coupling becomes
\[
\Xi_{\rm micro}=\kappa G_{\rm micro},
\qquad
G_{\rm micro}:=\chi_\sigma\Lambda_\phi^2/K_X.
\]
So the exact microscopic thresholds are
\[
G_{\rm fail}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_\infty(\kappa,\eta)},
\qquad
G_{\rm suff}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_0(\kappa,\eta)}.
\]
The exact phase diagram is:

- \(G_{\rm micro}\le G_{\rm fail}\): impossible,
- \(G_{\rm micro}\ge G_{\rm suff}\): guaranteed,
- only the bounded interval in between needs the full root solve.

Equivalent threshold surfaces are
\[
\chi_\sigma \le \frac{T_X Pe_{\rm req}}{\Lambda_\phi^2 L^2 \Delta_\infty}\quad\Rightarrow\quad\text{fail},
\]
\[
\chi_\sigma \ge \frac{T_X Pe_{\rm req}}{\Lambda_\phi^2 L^2 \Delta_0}\quad\Rightarrow\quad\text{succeed},
\]
and similarly for \(\Lambda_\phi^2\).

Soft-support limit:
\[
\Delta_0\to \frac12,
\qquad
\Delta_\infty\to 1,
\qquad
G_{\rm fail}\sim \frac{Pe_{\rm req}}{\kappa},
\qquad
G_{\rm suff}\sim \frac{2Pe_{\rm req}}{\kappa}.
\]
So very soft support is strongly disfavored.

Highly compliant-mouth limit:
\[
\Delta_0^{(\infty)}=
\frac{1-\operatorname{sech}(\sqrt\kappa)}{\kappa},
\qquad
\Delta_\infty^{(\infty)}=
\frac{\tanh(\sqrt\kappa)}{\sqrt\kappa},
\]
so
\[
G_{\rm fail}^{(\infty)}=
\frac{Pe_{\rm req}}{\sqrt\kappa\tanh(\sqrt\kappa)},
\qquad
G_{\rm suff}^{(\infty)}=
\frac{Pe_{\rm req}}{1-\operatorname{sech}(\sqrt\kappa)}.
\]
For \(\kappa\gg1\), these reduce to
\[
G_{\rm fail}^{(\infty)}\sim \frac{Pe_{\rm req}}{\sqrt\kappa},
\qquad
G_{\rm suff}^{(\infty)}\sim Pe_{\rm req}.
\]

## Stage 45 — parent-action projection of the microscopic gain

Starting from the parent matter energy
\[
H_\psi=
\int d^4X\left[\frac{\hbar^2}{2m}|D_i\psi|^2+V_{\rm conf}\rho+U(\rho)\right],
\]
with frozen EOS
\[
U(\rho)=K\rho^5/4,
\qquad
h(\rho)=\frac{5K}{4}\rho^4,
\qquad
h'(\rho)=5K\rho^3=\frac{m c_s^2(\rho)}{\rho},
\]
the local compressional quadratic energy is
\[
\delta H_{\rm comp}=
\frac12\int d^4X\;h'(\rho_*)(\delta\rho)^2.
\]

Project one source channel
\[
\delta\rho(s,y)=\sigma(s)\chi_\sigma(y)
\]
and one support channel entering the confinement as
\[
\delta V_{\rm conf}(s,y)=-g_\phi\chi_\phi(y)\phi(s).
\]
Then the exact reduced coefficients are
\[
\Theta_\sigma=h'(\rho_*)N_{\sigma\sigma},
\qquad
\Lambda_\phi=g_\phi O_{\sigma\phi},
\]
with parent overlap invariants
\[
N_{\sigma\sigma}=\int d^3y\,\chi_\sigma^2,
\qquad
N_{\phi\phi}=\int d^3y\,\chi_\phi^2,
\qquad
O_{\sigma\phi}=\int d^3y\,\chi_\sigma\chi_\phi.
\]

So the effective source susceptibility is
\[
\chi_\sigma^{\rm eff}=
\frac{1}{\Theta_\sigma}=
\frac{\rho_*}{m c_{s,*}^2 N_{\sigma\sigma}},
\]
and the microscopic gain becomes the explicit parent quantity
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 O_{\sigma\phi}^2}{m c_{s,*}^2 K_X N_{\sigma\sigma}}.
\]
Introducing the coherence factor
\[
C_{\sigma\phi}^2=
\frac{O_{\sigma\phi}^2}{N_{\sigma\sigma}N_{\phi\phi}},

aud via Cauchy–Schwarz, one gets the exact factorization
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
C_{\sigma\phi}^2,
\qquad 0\le C_{\sigma\phi}^2\le 1.
\]

## Stage 46 — parent-overlap threshold theorem

Combining the parent gain with the Stage-44 phase diagram gives exact parent thresholds
\[
g_{\phi,\rm fail}^2=
\frac{m c_{s,*}^2 K_X N_{\sigma\sigma} G_{\rm fail}}{\rho_* O_{\sigma\phi}^2},
\qquad
g_{\phi,\rm suff}^2=
\frac{m c_{s,*}^2 K_X N_{\sigma\sigma} G_{\rm suff}}{\rho_* O_{\sigma\phi}^2}.
\]
Equivalently, in coherence form,
\[
C_{\rm fail}^2=
\frac{m c_{s,*}^2 K_X G_{\rm fail}}{\rho_* g_\phi^2 N_{\phi\phi}},
\qquad
C_{\rm suff}^2=
\frac{m c_{s,*}^2 K_X G_{\rm suff}}{\rho_* g_\phi^2 N_{\phi\phi}}.
\]

There is an exact Cauchy no-go theorem:
if
\[
G_{\max}(g_\phi):=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
< G_{\rm fail}(\kappa,\eta),
\]
then no profile engineering of \(\chi_\sigma\) can rescue the branch.

Inserting
\[
G_{\rm fail}=
\frac{Pe_{\rm req}}{\kappa\Delta_\infty},
\qquad
G_{\rm suff}=
\frac{Pe_{\rm req}}{\kappa\Delta_0},
\qquad
\kappa=K_XL^2/T_X,
\]
one finds exact amplitude thresholds
\[
g_{\phi,\rm fail}^2=
\frac{m c_{s,*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2 \Delta_\infty},
\qquad
g_{\phi,\rm suff}^2=
\frac{m c_{s,*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2 \Delta_0}.
\]
So \(K_X\) cancels from the explicit prefactor and survives only through the geometry-shape functions \(\Delta_0,\Delta_\infty\).

## Stage 47 — parent equilibrium source/support alignment

The parent equilibrium law
\[
H(y)\,\delta\rho(s,y)+\delta V_{\rm conf}(s,y)=0,
\qquad H(y):=h'(\rho_*(y)),
\]
forces the aligned source profile
\[
\chi_\sigma(y)=g_\phi\chi_\phi(y)/H(y).
\]
So the overlap invariants become
\[
O_{\sigma\phi}=g_\phi I_1,
\qquad
N_{\sigma\sigma}=g_\phi^2 I_2,
\]
with
\[
I_1=\int d^3y\;\frac{\chi_\phi(y)^2}{H(y)},
\qquad
I_2=\int d^3y\;\frac{\chi_\phi(y)^2}{H(y)^2}.
\]
Therefore
\[
C_{\sigma\phi}^2=
\frac{I_1^2}{N_{\phi\phi} I_2}\le 1.
\]
In the thin active layer where \(H(y)\approx H_w\) is nearly constant,
\[
I_1=N_{\phi\phi}/H_w,
\qquad
I_2=N_{\phi\phi}/H_w^2,
\qquad
C_{\sigma\phi}^2=1.
\]
So the matched-layer branch is not arbitrary; it is the natural thin-layer limit of the parent equilibrium branch.

The exact eliminated-source support softening is
\[
\Delta K_X^{\rm (eq)}=g_\phi^2 I_1,
\qquad
G_{\rm eq}=\Delta K_X^{\rm(eq)}/K_X = g_\phi^2 I_1/K_X.
\]

## Stage 48 — explicit thin-wall confinement branch

For the explicit wall family
\[
V_{\rm conf}(r;a)=V_0 f\bigl((r-a)/\ell\bigr),
\]
with wall coordinate \(\xi=(r-a)/\ell\), a support displacement \(a\to a+\phi(s)\) gives
\[
\delta V_{\rm conf} = +\frac{V_0}{\ell} f'(\xi)\phi(s),
\]
so the parent loading amplitude is exactly
\[
g_\phi=V_0/\ell.
\]

The shell integral entering the equilibrium gain is
\[
I_1 = 4\pi\ell\bigl[a^2J_1+2a\ell J_2+\ell^2J_3\bigr],
\]
where
\[
J_n := \int d\xi\;\frac{\xi^n f'(\xi)^2}{H(\xi)}.
\]
For a centered symmetric wall layer, \(J_2=0\), so
\[
I_1 = 4\pi\ell\bigl[a^2J_1+\ell^2J_3\bigr].
\]
The exact equilibrium gain is
\[
G_{\rm eq}=4\pi V_0^2\Bigl[\frac{a^2J_1}{\ell}+2aJ_2+\ell J_3\Bigr]/K_X.
\]
In the thin-wall limit \(\ell\ll a\), the leading gain is
\[
G_{\rm eq}^{\rm(tw)}=\frac{4\pi a^2 V_0^2 J_1}{K_X\ell}.
\]

Comparing with the Stage-44 thresholds gives wall-amplitude surfaces
\[
V_{0,\rm fail}^2=\frac{K_X\ell G_{\rm fail}}{4\pi a^2J_1},
\qquad
V_{0,\rm suff}^2=\frac{K_X\ell G_{\rm suff}}{4\pi a^2J_1}.
\]
After inserting
\[
G_{\rm fail}=
\frac{Pe_{\rm req}}{\kappa\Delta_\infty},
\qquad
G_{\rm suff}=
\frac{Pe_{\rm req}}{\kappa\Delta_0},
\qquad
\kappa=K_XL^2/T_X,
\]
the explicit \(K_X\) factor cancels:
\[
V_{0,\rm fail}^2=
\frac{T_X\ell Pe_{\rm req}}{4\pi a^2L^2J_1\Delta_\infty},
\qquad
V_{0,\rm suff}^2=
\frac{T_X\ell Pe_{\rm req}}{4\pi a^2L^2J_1\Delta_0}.
\]

For an almost constant-compressibility active wall layer, \(H(\xi)\approx H_w\),
\[
J_1 = I_f/H_w,
qquad I_f:=\int d\xi\;f'(\xi)^2,
\]
so the thresholds reduce to
\[
V_{0,\rm fail}^2=
\frac{H_w T_X\ell Pe_{\rm req}}{4\pi a^2L^2I_f\Delta_\infty},
\qquad
V_{0,\rm suff}^2=
\frac{H_w T_X\ell Pe_{\rm req}}{4\pi a^2L^2I_f\Delta_0}.
\]

## Stage 49 — dimensionless wall figure of merit

For the same thin-wall matched branch, define
\[
W_{\rm wall}:=
\frac{4\pi a^2L^2J_1V_0^2}{T_X\ell}.
\]
Since \(\kappa=K_XL^2/T_X\), this is exactly
\[
W_{\rm wall}=\kappa G_{\rm eq}^{\rm(tw)}.
\]
The Stage-44 operator theorem therefore becomes
\[
W_{\rm fail}=\frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)},
\qquad
W_{\rm suff}=\frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\]
with exact wall-form theorem:

- \(W_{\rm wall}\le W_{\rm fail}\): fail,
- \(W_{\rm wall}\ge W_{\rm suff}\): succeed,
- only the narrow intermediate band still needs the full root solve.

If \(H(\xi)\approx H_w\), the wall control variable becomes
\[
W_H=
\frac{4\pi a^2L^2I_fV_0^2}{H_w T_X\ell}.
\]
So the explicit parent branch is now controlled by one dimensionless figure of merit rather than a diffuse set of amplitudes.

## Stage 50 — sech–Gaussian coherence resonance benchmark

For the explicit independent-profile benchmark
\[
\chi_\sigma(y)=\operatorname{sech}(y/w_f),
\qquad
\chi_\phi(y)=e^{-y^2/w_g^2},
\qquad
r:=w_g/w_f,
\]
the exact norms are
\[
N_{\sigma\sigma}=2w_f,
\qquad
N_{\phi\phi}=w_g\sqrt{\pi/2}.
\]
The overlap is
\[
O_{\sigma\phi}=w_f I(r),
\qquad
I(r):=\int_{-\infty}^{\infty} dx\;\operatorname{sech}(x)e^{-x^2/r^2}.
\]
So the coherence is
\[
C^2(r)=\frac{I(r)^2}{r\sqrt{2\pi}}.
\]

Parseval/Fourier-transform arguments give the exact duality
\[
I(r)=\frac{r}{\sqrt\pi}I(\pi/r),
\qquad
C^2(r)=C^2(\pi/r).
\]
Hence the self-dual stationary point is
\[
r_* = \sqrt\pi.
\]
Numerically,
\[
C_{\rm res}^2:=C^2(\sqrt\pi)=0.9944188364515293487\ldots,
\]
so the resonance penalty factor is
\[
P_{\rm res}:=1/C_{\rm res}^2=1.0056124877605762169\ldots.
\]
Thus the best independent sech–Gaussian mismatch branch misses the ideal matched-layer coherence only by about \(0.56\%\).

## Stage 51 — resonance-corrected thresholds

The general parent gain is
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
C_{\sigma\phi}^2.
\]
On the Stage-47 matched equilibrium branch, \(C_{\sigma\phi}^2=1\), so the matched gain is
\[
G_{\rm match}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}.
\]
Stage 49 repackaged this as the wall figure of merit
\[
W_{\rm wall}=
\frac{4\pi a^2L^2J_1V_0^2}{T_X\ell}.
\]
For the independent sech–Gaussian family,
\[
G_{\rm res}(r)=C^2(r)G_{\rm match},
\qquad
W_{\rm res}(r)=C^2(r)W_{\rm wall}.
\]
Therefore the exact profile-family thresholds are
\[
W_{\rm wall}\le \frac{Pe_{\rm req}}{C^2(r)\Delta_\infty}
\quad\Rightarrow\quad \text{fail},
\]
\[
W_{\rm wall}\ge \frac{Pe_{\rm req}}{C^2(r)\Delta_0}
\quad\Rightarrow\quad \text{succeed}.
\]
At resonance \(r=\sqrt\pi\), this becomes
\[
W_{\rm wall}\le \frac{P_{\rm res} Pe_{\rm req}}{\Delta_\infty}
\quad\Rightarrow\quad \text{fail on the resonance family},
\]
\[
W_{\rm wall}\ge \frac{P_{\rm res} Pe_{\rm req}}{\Delta_0}
\quad\Rightarrow\quad \text{succeed on the resonance family}.
\]
So the explicit independent-profile benchmark modifies the Stage-49 wall thresholds by only the tiny factor \(P_{\rm res}\approx1.00561249\).
# 5PN stages 52–72 notes

This bundle continues the moving-throat/support-source chain from the resonance-threshold point
through the explicit Family-1 branch closure and the first direct comparison with the minimal
isotropic passive/outgoing quadrupole module.

## Stage 52 — final reduced verdict for the support/source program

The matched-branch theorem remains

- `W_wall <= Pe_req / Delta_inf`  -> fail,
- `W_wall >= Pe_req / Delta_0`    -> succeed,

and the explicit sech–Gaussian resonance family only perturbs those thresholds by the tiny factor
`P_res = 1.005612487760576...`, so the genuinely profile-sensitive sub-bands are only about `0.56%`
wide.

## Stage 53 — explicit GNLS wall-shell reduction

Projecting the parent GNLS quadratic shell energy onto the wall-support mode `q(s) chi_phi(y)` gives

- `T_X = hbar^2 N_(phi phi) / (4 m rho_w)`,
- `K_X = H_w N_(phi phi) + hbar^2 G_(phi phi)/(4 m rho_w)`,

and on the matched thin-wall branch the support/source fixed-point coupling is exactly the wall
figure of merit:
`Xi = W_wall`.

## Stage 54 — canonical tanh-wall branch

For the canonical wall
`f(xi) = (1 + tanh xi)/2`, `chi_phi = f'(xi) = (1/2) sech^2 xi`,
the exact moments are

- `I_f = 1/3`,
- `I_g = 4/15`,
- `I_g / I_f = 4/5`.

With the natural local Robin closure `K_m = T_X / ell`, the explicit branch is

- `eta = L/ell`,
- `kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L/ell)^2`,
- `W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2)`.

## Stage 55 — explicit branch thresholds

Writing
`chi_s = m c_(s,w) L / hbar`,
`Lambda_ell = L/ell`,
`Upsilon_w = 4 rho_w^2 V0^2 /(hbar^2 c_(s,w)^2)`,
the branch map is

- `kappa = 4 chi_s^2 + (4/5) Lambda_ell^2`,
- `eta = Lambda_ell`,
- `W_wall = Upsilon_w Lambda_ell^2`.

So the explicit thresholds are

- `Upsilon_fail = Pe_req / [Lambda_ell^2 Delta_inf(kappa,eta)]`,
- `Upsilon_suff = Pe_req / [Lambda_ell^2 Delta_0(kappa,eta)]`.

The shell-gradient and compression-dominated asymptotics from the notes are verified directly.

## Stages 56–57 — Family-1 geometry map and healing lock

Using the carried reference values

- `L/a = 37/20`,
- `ell/a = 1/20`,

the Family-1 branch fixes

- `Lambda_ell = 37`,
- `eta = 37`.

With the healing-width closure `ell = hbar/(2 m c_(s,w))`, the same branch fixes

- `chi_s = 37/2`,
- `kappa = 12321/5`,
- `alpha = sqrt(kappa) = 111/sqrt(5) ≈ 49.6407091`.

So the only remaining explicit branch input is the wall-depth amplitude.

## Stage 58 — explicit Family-1 threshold window

At `(kappa,eta) = (12321/5, 37)`, the exact support/source scales are

- `Delta_0 ≈ 1.73302079021525e-4`,
- `Delta_inf ≈ 2.01447565540522e-2`.

Hence

- `Upsilon_fail ≈ 0.0362605617972939 Pe_req`,
- `Upsilon_suff ≈ 4.21495341569977 Pe_req`,
- `Theta_fail ≈ 3.62605617972939e-4 Pe_req`,
- `Theta_suff ≈ 4.21495341569977e-2 Pe_req`

after `Upsilon_w = 100 Theta_w`.

## Stages 59–60 — exact `n=5` wall-depth lock and shell-weighted extraction

For the frozen `n=5` EOS,
`h(rho) = m c_s(rho)^2 / 4`, so with
`mu_* = lambda_mu h_w`
and the healing-width lock one gets

- `Theta_w = lambda_mu^2 rho_w^2 /(16 ell^2)`.

On the normalized Family-1 wall this becomes
`Theta_w = 25 lambda_mu^2 rho_w^2`.

Using the explicit Family-1 wall profile and canonical support weight gives

- `<rho_r>_chi ≈ 0.192619005556493`,
- `<rho_r^2>_chi ≈ 0.162745294003265`,
- `Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2`,
- `Theta_w^(J)   ≈ 0.927552032539308 lambda_mu^2`.

## Stage 61 — explicit Family-1 wall-depth verdict

Comparing the extracted `Theta_w` values to the Stage-58 window gives

- `Pe_suff^(chi) ≈ 96.5285247264386 lambda_mu^2`,
- `Pe_fail^(chi) ≈ 11220.5441626259 lambda_mu^2`,
- `Pe_suff^(J)   ≈ 22.0062226330754 lambda_mu^2`,
- `Pe_fail^(J)   ≈ 2558.01892349205 lambda_mu^2`.

So the explicit Family-1 wall-depth supply is not naturally starved; the remaining question becomes the
quadrupole-side demand `Pe_req`.

## Stages 62–64 — Family-1 demand map in `Pe`, `zeta_req`, and `Pi_tr/C_mix`

For the Family-1 branch,

- `y_F1 tan y_F1 = 37`,
- `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2) ≈ 1.00005192880220`.

The exact constructive support map is
`zeta_F1(Pe) = A_F1 Omega_Pe^2`,
with hard ceiling
`zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601`.

Converting the Stage-61 `Pe_req` windows through this map gives, at `lambda_mu = 1`,

- `zeta_suff^(chi) ≈ 2.46622291347846`,
- `zeta_fail^(chi) ≈ 2.46752913273870`,
- `zeta_suff^(J)   ≈ 2.44257571477179`,
- `zeta_fail^(J)   ≈ 2.46752736855058`.

Using
`Pi_tr = C_mix Q(zeta;eps_blk)`,
`Q(zeta;eps_blk) = [1 + (1 - 2 eps_blk) zeta]/[1 - eps_blk zeta]`,
the unblocked (`eps_blk = 0`) explicit product window becomes

- `Pi_tr/C_mix <= 3.46622291347846`  -> guaranteed success,
- `Pi_tr/C_mix >= 3.46752913273870`  -> guaranteed failure,
- `Pi_tr/C_mix <  3.46752922945601`  -> hard explicit ceiling.

## Stages 65–70 — master reduced quadrupole residual and pure loading-ratio collapse

The whole reduced moving-throat PDE is compressed to one scalar residual
`R_quad = zeta_req - zeta_phys(Pe_*)`,
with `Pe_*` selected by the support/source fixed-point law.

After the exact demand-cancellation step, the normalized support theorem depends only on the loading ratio

`rho_alpha = alpha_req / alpha_mix`,

with support demand
`zeta_req = (rho_alpha - 1) / [1 - eps_blk (2 - rho_alpha)]`.

So the explicit Family-1 theorem is finally

- `rho_alpha <= rho_suff^(chi)(lambda_mu;eps_blk)`  -> guaranteed success,
- `rho_alpha >= rho_fail^(chi)(lambda_mu;eps_blk)`  -> guaranteed failure,
- `rho_alpha < rho_max^(F1)(eps_blk)`               -> hard ceiling,

and in the natural unblocked case

- `rho_alpha <= 3.46622291347846`  -> guaranteed success,
- `rho_alpha >= 3.46752913273870`  -> guaranteed failure,
- `rho_alpha <  3.46752922945601`  -> hard ceiling.

So by Stage 70 the explicit support/source side is completely reduced to one number:
the outgoing-branch loading ratio `rho_alpha`.

## Stages 71–72 — loading ratio from the minimal isotropic conservative quadrupole module

The natural contact-plus-pole reading of the minimal isotropic conservative precursor

`Y_Q^cons(omega) = c0 + c1/(1 - omega^2/Omega_Q^2)`

gives

- `c0 = 1/rho_alpha`,
- `c1 = (rho_alpha - 1)/rho_alpha`,
- `zeta_req = c1/c0`.

For the minimal isotropic module
`c0 = 3/4`, `c1 = 1/4`,
the exact loading data are

- `rho_alpha = 4/3`,
- `zeta_req = 1/3`,
- `Pi_tr = (4/3) C_mix`.

This lies in the symmetric-lowest-twin regime and, on the explicit Family-1 branch, it already satisfies

- `zeta_req = 1/3 < A_F1 ≈ 1.00005192880220`.

So the explicit Family-1 branch succeeds at **zero transport bias** on the natural minimal isotropic passive/outgoing quadrupole module.

## Bottom line from Stage 72

By the end of this bundle, the explicit Family-1 support/source side is no longer the active reduced bottleneck. On the natural minimal isotropic contact-plus-pole branch it is already comfortably satisfied. The remaining reduced theorem question becomes whether the actual grouped-`P2` / geometry moving-throat branch really realizes that minimal isotropic passive/outgoing quadrupole module, which is the right bridge back into the 5PN grouped-`P2` program.

# 5PN / Moving-Throat Continuation — Stages 73–90

This batch continues directly from the Stage 72 result that the explicit Family-1 support/source side is no longer the active bottleneck. The live question is whether the actual grouped-\(P_2\)/geometry branch realizes the minimal isotropic contact-plus-pole conservative quadrupole module, and then whether the passive/outgoing \(l=2\) branch carries the canonical outgoing normalization.

## Stage 73
Recasts the post-72 status in theorem language: the explicit Family-1 support/source branch already succeeds on the minimal isotropic branch with
\[
\rho_\alpha=\frac43,\qquad \zeta_{\rm req}=\frac13,\qquad Pe_{\rm req}=0,
\]
so the remaining reduced theorem gap is no longer on the support/source side.

## Stage 74
Derives the `3/4 + 1/4` conservative module directly from the 3PN conservative split.
If the isotropic grouped-\(P_2\) branch is carried by one effective pole and the geometry lane is static through \(O(\omega^4)\), then the minimal isotropic branch identity forces
\[
K_{\rm pole}=\frac{K_0}{4},\qquad K_{\rm geom}=\frac{3K_0}{4},
\]
hence
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
\]

## Stage 75
Allows the geometry lane to carry dynamic even moments and derives the exact obstruction formula
\[
c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\qquad
\epsilon_2=\frac{\Omega_Q^2 K_{(g,2)}}{K_{\rm pole}},
\qquad
\epsilon_4=\frac{\Omega_Q^4 K_{(g,4)}}{K_{\rm pole}}.
\]
So the `3/4 + 1/4` split is exact iff the geometry lane is static at the relevant orders.

## Stage 76
Freezes the updated reduced status: the only remaining reduced ambiguity is whether the actual moving-throat geometry lane is dynamically inert through \(O(\omega^4)\).

## Stage 77
Proves the isotropic geometry-decoupling theorem. In the exact isotropic quadratic wall theory, the \(l=0\) scalar/geometry lane is orthogonal to the grouped real \(l=2\) bundle, so no dynamic \(O(\omega^2)\) or \(O(\omega^4)\) geometry contamination can enter the grouped quadrupole carrier at linear order:
\[
K_{(g,2)}=K_{(g,4)}=0,\qquad \epsilon_2=\epsilon_4=0.
\]

## Stage 78
Shows that if weak anisotropy later induces an \(l=0\leftrightarrow l=2\) mixing source, the first nonzero geometry contamination appears only at second order in that mixing:
\[
\epsilon_2,\epsilon_4 = O(\chi^2).
\]
So the Stage-74 split is perturbatively stable.

## Stage 79
Records the actual reduced-branch verdict: on the natural isotropic branch,
\[
\epsilon_2=\epsilon_4=0,
\]
hence
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
\qquad
\rho_\alpha=\frac43,\qquad
\zeta_{\rm req}=\frac13.
\]

## Stage 80
Shows that once the actual isotropic grouped-\(P_2\) one-pole branch is accepted, the full reduced 2.5PN normalization problem collapses to one scalar defect
\[
N_Q:=\frac{\overline K_0}{\overline K_0^{\rm target}}.
\]
All low-frequency coefficients scale by the same factor on that branch:
\[
\overline K_2,\ \overline K_4,\ \overline\Gamma_5 \propto N_Q.
\]

## Stage 81
Proves the explicit Family-1 support/source theorem is automatic on the actual isotropic branch. Since
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk})=\frac{1}{3-2\epsilon_{\rm blk}},
\]
any explicit family with \(\zeta_{\max}>1\) already passes throughout the admissible blocked regime. Family-1 has
\[
\zeta_{\max}^{(F1)}\approx 2.46752922945601 > 1,
\]
so it is automatically safe.

## Stage 82
Freezes the reduced finish line: the grouped-\(P_2\) conservative branch is geometry-clean, the support/source side is automatic, and the only remaining reduced theorem gap is the single normalization defect \(N_Q-1\).

## Stage 83
Factors the last 2.5PN obstruction into a conservative and an outgoing piece. Introducing one outgoing-normalization factor \(\chi_Q\),
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q \sigma_Q^{\rm can}\omega^5}+O(\omega^6),
\]
one gets
\[
\frac{\overline\Gamma_5}{\overline\Gamma_5^{\rm target}} = \chi_Q N_Q,
\]
and the full odd normalization condition is
\[
\hat m_0^{\,2}\chi_Q N_Q = 1.
\]

## Stage 84
Uses the natural point-particle source map \(\hat m_0\to 1\) to show the remaining reduced obstruction is purely outgoing:
\[
N_Q=\frac{1}{\chi_Q}.
\]
So all remaining reduced uncertainty is now in the outgoing normalization factor \(\chi_Q\).

## Stage 85
Shows higher odd retarded data beginning at \(O(\omega^7)\) are irrelevant to the 2.5PN theorem. The only live retarded obstruction at 2.5PN is the leading \(\omega^5\) normalization factor \(\chi_Q\).

## Stage 86
States the conditional reduced closure:
if \(\chi_Q=1\), the reduced GR-like point-particle 2.5PN theorem is closed; if not, the whole remaining failure is measured by \(\Delta_Q=\chi_Q-1\).

## Stage 87
Computes the exact outgoing spherical \(l=2\) DtN fingerprint:
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}-\frac{2z^6}{27}-i\frac{z^7}{27}+O(z^8),
\]
and therefore
\[
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}-\frac{11z^6}{729}-i\frac{z^7}{243}+O(z^8).
\]

## Stage 88
Matches the Stage-87 DtN fingerprint to the retarded grouped-\(P_2\) one-pole-plus-contact module and fixes
\[
\chi_Q=1
\]
on the canonical compact passive/outgoing branch. A deformed DtN branch with
\[
\Lambda_2^{\rm def}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\xi_Q\frac{z^5}{9}+O(z^6)
\]
simply gives \(\chi_Q=\xi_Q\).

## Stage 89
Inserts \(\chi_Q=1\) into the reduced normalization stack and closes the reduced 2.5PN theorem on the canonical outgoing branch. In the strict point-particle limit,
\[
N_Q=1,
\]
and the canonical invariant coefficients are exactly
\[
\overline K_0=\frac{54Gc_s^5}{5a^5c^5},\qquad
\overline K_2=\frac{6Gc_s^3}{5a^3c^5},\qquad
\overline K_4=\frac{8Gc_s}{15ac^5},\qquad
\overline\Gamma_5=\frac{2G}{5c^5}.
\]

## Stage 90
Builds the first explicit isotropic DtN deformation algebra. For
\[
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6),
\]
the normalized outgoing factor is
\[
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0},
\]
with \(\Sigma_2,\Sigma_4\) fixed by the requirement that the canonical even moments remain unchanged. So the only isotropic branch data that can shift \(\chi_Q\) are \(\beta,\Sigma_0,\Sigma_5\).

## Net result of Stages 73–90

Inside the present reduced hierarchy, the conservative grouped-\(P_2\)/geometry branch now reaches the same minimal isotropic `3/4 + 1/4` module that the 2.5PN audit wanted, the explicit Family-1 support/source side is automatic, and the outgoing \(l=2\) DtN branch fixes the last reduced scalar to
\[
\chi_Q=1
\]
on the canonical compact branch. So the reduced nonspinning point-particle 2.5PN theorem is closed on that canonical branch; what remains genuinely PDE-facing is branch realization and explicit isotropic DtN deformation data.
# 5PN / Moving-Throat Continuation — Stages 91–100

This batch continues directly from Stage 90, where the isotropic outgoing normalization factor was reduced to
\[
\chi_Q=\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
\]
The goal of Stages 91–100 is to stop treating \((\beta,\Sigma_0,\Sigma_5)\) as abstract branch labels and push them through explicit outlet/core models.

## Stage 91
Classifies the robustness classes of \(\chi_Q\).

- Pure overall scaling is exactly invisible:
  \[
  \Lambda_2^{\rm def}=S\Lambda_2^{\rm out}
  \quad\Longrightarrow\quad
  \chi_Q=1.
  \]
- Pure scale+argument deformation preserves the already-fixed even moments only on the natural positive branch \(\beta=1\), so it is also harmless.
- A genuine additive isotropic throat-core channel can move \(\chi_Q\) while leaving the lower even moments canonical:
  \[
  \chi_Q=\frac{3(S+9\Sigma_5)}{3S-\Sigma_0}
  \qquad(\beta=1).
  \]
- The exact preservation submanifold is
  \[
  \Sigma_5=\frac{S(1-\beta^5)}{9}-\frac{\Sigma_0}{27}.
  \]

## Stage 92
Linearizes around the canonical outgoing branch:
\[
S=1+\varepsilon s,\qquad
\beta=1+\varepsilon b,\qquad
\Sigma_0=\varepsilon a_0,\qquad
\Sigma_5=\varepsilon a_5.
\]
Then
\[
\chi_Q
=
1+\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)+O(\varepsilon^2).
\]
So the minimal isotropic branch-selection data are the triple
\[
(b,a_0,a_5),
\]
and first-order preservation requires
\[
5b+\frac{a_0}{3}+9a_5=0.
\]

## Stage 93
Introduces the first explicit isotropic geometric outlet:
\[
\Lambda_2^{\rm R}(z)=\Lambda_2^{\rm out}(z)+\rho_R.
\]
The normalized response is
\[
\widehat Y_2^{\rm R}(z)
=
1+\frac{z^2}{9-3\rho_R}
+\frac{4-\rho_R}{9(3-\rho_R)^2}z^4
+i\frac{z^5}{27-9\rho_R}+O(z^6),
\]
so
\[
\chi_Q^{\rm R}=\frac{3}{3-\rho_R}.
\]
A pure Robin core therefore generically shifts both the even branch and the odd normalization.

## Stage 94
Adds the first explicit isotropic hidden mixed side-channel:
\[
\Lambda_2^{\rm mix}(z)
=
\Lambda_2^{\rm out}(z)
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}+O(z^6).
\]
The even-preserving conditions force
\[
\kappa_W=-\frac19,
\qquad
\sigma_W=0.
\]
So a standalone isotropic hidden pole cannot sit on the canonical even branch unless it is absent. Its normalization factor is
\[
\chi_Q^{\rm mix}
=
\frac{3(1-9\sigma_W\gamma_W)}{3+\sigma_W}.
\]

## Stage 95
Combines the Robin core and the mixed side-channel:
\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}+O(z^6).
\]
Imposing the canonical even branch yields exactly two solutions:
\[
(\rho_R,\kappa_W)=(\sigma_W,0)
\quad\text{or}\quad
(\rho_R,\kappa_W)=\left(4\sigma_W,\frac13\right).
\]
The second is the nontrivial compensated branch. On it,
\[
\chi_Q^{\rm hyb}=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
\]
and canonical odd normalization is preserved iff
\[
\gamma_W=\frac19.
\]
With that value the whole outlet collapses to a pure harmless scale deformation.

## Stage 96
Freezes the outlet-model classification:

1. pure Robin loading generically shifts \(\chi_Q\);
2. a standalone isotropic mixed pole is too rigid and cannot preserve the canonical even branch unless it vanishes;
3. the hybrid Robin–mixed outlet admits one nontrivial compensated canonical-even branch, with
   \[
   \rho_R=4\sigma_W,\qquad \kappa_W=\frac13,\qquad \gamma_W=\frac19
   \]
   on the fully preserved canonical branch.

## Stage 97
Replaces the reduced outlet coefficients by a concrete two-channel core model with internal variables

- \(s(\omega)\): static shell/compliance mode,
- \(q(\omega)\): mixed \(A_w/F_{\mu w}\)-type side-channel mode.

The linear core system
\[
\begin{pmatrix}
K_s & \lambda\\
\lambda & -K_q D_W^{\rm bare}(z)
\end{pmatrix}
\binom{s}{q}
=
u\binom{g_s}{g_q}
\]
gives the exact Schur-complement outlet
\[
\delta\Lambda_{\rm core}(z)
=
\frac{g_s^2}{K_s}
-
\frac{(K_s g_q-\lambda g_s)^2}
{K_s\big(K_sK_q D_W^{\rm bare}(z)+\lambda^2\big)}.
\]
Defining
\[
r_c=\frac{\lambda^2}{K_sK_q},
\]
the reduced Robin–mixed parameters are
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=\frac{(K_s g_q-\lambda g_s)^2}{K_s^2K_q(1+r_c)},
\qquad
\kappa_c=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_c=\frac{\gamma_0}{1+r_c}.
\]

## Stage 98
Determines exactly when the concrete core lands on the compensated canonical branch. The nontrivial compensation condition is
\[
\rho_c=4\sigma_c,\qquad \kappa_c=\frac13,\qquad \gamma_c=\frac19.
\]
This is equivalent to the exact coupling-balance law
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
\]
so
\[
g_q=
\frac{g_s}{2K_s}\left(2\lambda\pm\sqrt{K_sK_q+\lambda^2}\right).
\]
The bare mixed channel must itself be a scale-deformed copy of the canonical compact outgoing branch:
\[
\kappa_0=\frac{1+r_c}{3},\qquad
\gamma_0=\frac{1+r_c}{9}.
\]
On that surface,
\[
\delta\Lambda_{\rm core}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\qquad
\sigma_*=\frac{g_s^2}{4K_s},
\]
and the normalized outgoing fingerprint stays exactly canonical.

## Stage 99
Realizes the bare mixed side-channel geometrically as the first D/N half-wave on an auxiliary tube of length \(L_W\):
\[
k_W=\frac{\pi}{2L_W},
\qquad
\Omega_W=\frac{\pi c_s}{2L_W}.
\]
In outlet variables \(z=a\omega/c_s\),
\[
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
\]
The compensation condition \(\kappa_0=(1+r_c)/3\) fixes
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
\]
If the bare mixed outlet is a pure-scale deformation of the canonical compact outgoing \(l=2\) branch,
\[
D_W^{\rm bare}(z)=(1+r_c)\left(1-\frac{z^2}{3}-i\frac{z^5}{9}\right)+O(z^6),
\]
then the hybridization factor is removed exactly and
\[
\kappa_c=\frac13,\qquad \gamma_c=\frac19.
\]

## Stage 100
Freezes the concrete outlet-core status. The surviving PDE-facing question is no longer “is there some deformed outlet?” It is much sharper:

> Does the actual moving-throat core realize the concrete coupling-balance surface together with the auxiliary D/N-tube normalization?

On that surface the effective outlet preserves the canonical normalized outgoing fingerprint exactly.

## Net result of Stages 91–100

The isotropic outgoing-branch ambiguity is no longer open-ended. The first explicit outlet audit shows:

- pure scale deformations are harmless,
- pure Robin loading and a standalone isotropic mixed pole do not by themselves preserve the canonical branch,
- a specific compensated Robin–mixed outlet does preserve it,
- that compensated branch is realized by a concrete two-channel throat-core model,
- and the mixed side-channel can be given a direct finite D/N tube realization.

So the next honest step is Stage 101: extract the concrete core parameters \((K_s,K_q,\lambda,g_s,g_q)\) from a parent-action throat-core ansatz rather than leaving them as reduced variables.

# 5PN stages 101–120 notes

This bundle continues the chain immediately after the compensated outlet/core result at Stage 100.
The focus of this block is the next honest microscopic gate: replace the reduced outlet/core
coefficients by explicit parent-action overlaps, then carry that data all the way into the
mouth-layer fixed-point problem.

## Stage 101 — parent-action extraction of the core parameters

The reduced two-channel core variables from Stages 97–100 are replaced by explicit overlap
formulas from one concrete GNLS + localized-Maxwell throat-core ansatz.

The shell/compliance mode gives
\[
K_s
=
4\pi a^2\left(
\frac{H_w\ell}{3}
+
\frac{\hbar^2}{15m_\psi\rho_w\,\ell}
\right),
\]
and on the healing-locked shell branch
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\qquad
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
\]

The mixed D/N half-wave gives
\[
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi}{\sqrt2\,L_W^{3/2}}.
\]

The shell–mixed hybridization and shell mouth coupling become
\[
\lambda=-q_* v_{w0}\mathcal I_{sq},
\qquad
\mathcal I_{sq}=\frac{8\sqrt2}{3}a^2\ell\sqrt{L_W},
\qquad
g_s=\mathcal T_m \frac{4\pi a^2\ell}{3}.
\]

So the Stage-97 core matrix is no longer abstract; every entry now has an explicit parent
overlap meaning. This is exactly the next gate identified in the moving-throat notes. fileciteturn30file16

## Stages 102–103 — collapse to normalized parent ratios

The exact core-balance condition
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2
\]
collapses to the two dimensionless ratios
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\]
through
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.
\]

The D/N mixed-tube length becomes
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
\]

So the surviving outlet/core theorem gate is no longer “what are the reduced coefficients?”
It is only: which branch point \((\mathfrak r,\mathfrak g)\) the actual GNLS + localized-Maxwell
core selects. fileciteturn30file16turn30file5

## Stages 104–107 — explicit Family-1 bridge

To keep the executable chain sequential, I added the missing Family-1 bridge audits that the
later notes assume implicitly.

Using the carried Family-1 geometry \(L/a=37/20\) together with the D/N length law gives
\[
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
\approx 1.77799353547498.
\]

The two compensated coupling branches are
\[
\mathfrak g_\pm^{F1}
=
\mathfrak r_{F1}\pm\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
numerically
\[
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
\]

These bridge scripts also verify the useful ordering
\[
\frac{2}{\pi}<\mathfrak g_-^{F1}<\frac{\pi}{4}<1<\mathfrak g_+^{F1},
\]
which is exactly the window the later positive-source theorem needs.

## Stages 108–111 — positive mouth-source selection

For any positive normalized axial mouth source profile \(\sigma(z)\) on the first D/N interval,
the mouth-bias factor is
\[
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz,
\]
so positivity immediately forces
\[
0\le \mathfrak g[\sigma]\le 1.
\]

Since \(\mathfrak g_+^{F1}>1\), the upper compensated branch is impossible under any positive
localized mouth source, while \(\mathfrak g_-^{F1}\in(0,1)\), so the lower branch is the unique
physically admissible canonical candidate.

The first explicit positive families then show that the lower branch is easy to reach:

- self-matched derivative source:
  \[
  \mathfrak g_{\rm match}=\frac{\pi}{4},
  \]
  only \(3.61\%\) away in traction from exact lower-branch compensation;

- convex derivative/uniform family:
  \[
  \sigma_\xi=(1-\xi)k\cos(kz)+\xi/L,
  \]
  reaches the exact lower branch at
  \[
  \xi_*\approx 0.183918405511540;
  \]

- slab and truncated-exponential penetration laws reach the same branch at
  \[
  x_*^{\rm slab}\approx 0.797839360904564,
  \qquad
  x_*^{\exp}\approx 0.662765402623161.
  \]

So by Stage 111 the branch-selection ambiguity has collapsed: the lower compensated branch is
the unique positive-source branch and is reached with moderate penetration, not by exotic
sign-changing mouth forcing.

## Stages 112–115 — explicit GNLS + localized-Maxwell mouth boundary layer

The abstract positive-source family is replaced by the first explicit boundary-layer law.
With the mouth free energy
\[
F_{\rm mouth}[\sigma]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\!\left(\ln\frac{\sigma}{\sigma_*}-1\right)
+
V_m(z)\sigma
\Big],
\qquad
V_m(z)\approx V_1 z,
\]
and zero-flux Onsager current, the exact normalized source law is
\[
\sigma_\Pi(z)=\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})},
\qquad
\Pi=\frac{V_1L}{\Theta_\sigma}.
\]

Its exact mouth-bias factor is
\[
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}
{(4\Pi^2+\pi^2)(e^\Pi-1)},
\]
with range \(2/\pi\to1\) as \(\Pi:0^+\to\infty\).

Solving \(\mathfrak g_\Pi=\mathfrak g_-^{F1}\) gives the unique canonical Family-1 point
\[
\Pi_* \approx 1.50882951349316.
\]

So the parent threshold is now concrete:
\[
\Pi_m=\frac{L}{\Theta_\sigma}(T_m-q_*A_0')=\Pi_*,
\]
with local sensitivity
\[
\mathfrak g_*' \approx 0.0714453558083196.
\]

At this point the mouth problem is no longer branch sign or family choice. It is the
single microscopic bias question: does the real mouth layer select \(\Pi_m\approx1.51\)?

## Stages 116–120 — coupled mouth fixed point and explicit core-to-mouth gain map

The next honest step is to couple the shell/compliance and mixed Maxwell channels directly in
the mouth layer. The exact scalar D/N response kernel to the exponential source is
\[
\mathcal S(\Pi,\kappa)
=
\frac{
\Pi\left[\kappa\tanh\kappa+\Pi\left(e^{-\Pi}\operatorname{sech}\kappa-1\right)\right]
}{
(1-e^{-\Pi})(\kappa^2-\Pi^2)
},
\]
with \(\mathcal S(\Pi,0)=1\).

So the fully coupled mouth bias obeys
\[
\Pi = \sum_\alpha M_\alpha \mathcal S(\Pi,\kappa_\alpha).
\]

On the first explicit Family-1 branch with one static shell lane and one mixed D/N half-wave,
\[
\kappa_s=0,
\qquad
\kappa_q=\frac{\pi}{2},
\]
this reduces to
\[
\Pi = M_s + M_q \mathcal S_q(\Pi),
\qquad
\mathcal S_q(\Pi)=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
\]

At the canonical point,
\[
\mathcal S_q(\Pi_*)\approx 0.658075937605429,
\]
so the exact gain line is
\[
M_s=\Pi_* - M_q \mathcal S_q(\Pi_*).
\]

Imposing the outlet-consistent shell:mixed ratio \(4:-1\) collapses the mouth problem to
\[
\Pi = \Sigma_m[4-\mathcal S_q(\Pi)].
\]
At \(\Pi_*\),
\[
\Sigma_m^*\approx 0.451485277739088,
\qquad
M_s^*\approx 1.80594111095635,
\qquad
M_q^*\approx -0.451485277739088.
\]

Finally, Stage 120 removes the last abstract gain pair entirely. Using the exact Stage-97
Schur complement,
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=\frac{(K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)},
\]
the actual mouth-layer gains are
\[
M_s=\frac{L g_s^2}{K_s\Theta_\sigma},
\qquad
M_q=
-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
\]

So by the end of Stage 120 the mouth fixed-point ambiguity has collapsed from a free source
profile and a free gain pair to one explicit set of parent core quantities. The next clean step
is to normalize these gains and carry them into the self-consistent branch law beyond 120. fileciteturn30file11turn30file16
# Manifest — 5PN stages 121–140

Helper module:
- `fivepn_stage121_140_common.py`

Scripts:
- `5pn_stage121_normalized_mouth_gain_family.py`
- `5pn_stage122_family1_actual_mouth_gains.py`
- `5pn_stage123_selfmatched_mouth_susceptibility.py`
- `5pn_stage124_mouth_gain_status.py`
- `5pn_stage125_selfconsistent_mouth_branch.py`
- `5pn_stage126_equal_normalized_singular_limit.py`
- `5pn_stage127_unique_regular_canonical_branch.py`
- `5pn_stage128_mouth_branch_selection_status.py`
- `5pn_stage129_positive_deformation_expansion.py`
- `5pn_stage130_first_order_rigidity_kernel.py`
- `5pn_stage131_representative_positive_families.py`
- `5pn_stage132_mouth_rigidity_status.py`
- `5pn_stage133_full_profile_mouth_potential.py`
- `5pn_stage134_first_order_selected_correction.py`
- `5pn_stage135_family1_actual_correction.py`
- `5pn_stage136_full_mouth_correction_status.py`
- `5pn_stage137_coevolving_core_mouth_map.py`
- `5pn_stage138_frozen_traction_fixedpoint.py`
- `5pn_stage139_renormalized_canonical_branch.py`
- `5pn_stage140_core_mouth_coevolution_status.py`

Each script has a paired `_output.txt` file with its run output.
# 5PN continuation — Stages 121–140

This batch covers the next explicit mouth/core chain after Stage 120.

## What is now fixed

Stages 121–123 rewrite the mouth gains in normalized parent variables and collapse the explicit core-to-mouth map to

\[
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2},
\qquad
R_q:=-\frac{M_q}{M_s}=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
\]

On the exact compensation family
\[
\mathfrak g_c=\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2},
\]
that ratio collapses to
\[
R_q=\frac14.
\]
So the Stage-118 outlet-consistent closure is derived rather than assumed. On the self-matched mouth-susceptibility closure the overall shell gain becomes
\[
\Sigma_0=M_s=\frac{20}{9}\,\widehat T_m^2.
\]

Stage 122 evaluates this on the explicit Family-1 branch:

- natural equal-normalized branch:
  \[
  M_s\approx 1.66854,
  \qquad
  M_q\approx -0.24270;
  \]
- exact compensated branch:
  \[
  M_s\approx 1.80594,
  \qquad
  M_q\approx -0.45149.
  \]

The corresponding normalized traction amplitudes differ by only about 4.04%.

## Branch selection

Stages 125–128 close the self-consistent Family-1 mouth branch law,
\[
\Pi=\Sigma_0\bigl[1-R_q(\Pi)\,\mathcal S_q(\Pi)\bigr],
\qquad
R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

The positive exponential mouth family obeys
\[
0<\mathfrak g_\Pi<1
\quad\text{for every finite }\Pi>0,
\]
so the equal-normalized branch \(\mathfrak g_c=1\) is only a singular point-source limit. The upper compensated branch is impossible because \(\mathfrak g_+^{F1}>1\). The lower compensated branch is the unique regular finite branch, reached at
\[
\Pi_*\approx 1.50882951349316,
\qquad
\widehat T_{m,*}\approx 0.901484054174204.
\]

## Mouth rigidity under positive non-exponential deformations

Stages 129–132 derive the first-order rigidity kernel around the canonical point. For any positive normalized deformation, only two source moments matter,
\[
\bar g_\varsigma=\int_0^1\varsigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx,
\qquad
\bar S_\varsigma=\int_0^1\varsigma(x)K_q(x)dx,
\]
and the traction shift is
\[
\delta\widehat T_m
=
\epsilon\Bigl[
A_T(\bar g_\varsigma-\mathfrak g_*)
+
B_T(\bar S_\varsigma-\mathcal S_*)
\Bigr].
\]
The explicit coefficients are
\[
A_T\approx -4.27263956256927,
\qquad
B_T\approx 0.134875005736706,
\]
so overlap changes dominate mixed-kernel changes by a factor
\[
|A_T|/B_T\approx 31.68.
\]

Representative families:

- uniform broadening raises the canonical point,
  \[
  \frac{\delta\Pi_u}{\epsilon}\approx 1.69941,
  \qquad
  \frac{\delta\widehat T_{m,u}}{\epsilon}\approx 0.508756;
  \]
- self-matched derivative sharpening lowers it,
  \[
  \frac{\delta\Pi_d}{\epsilon}\approx -0.382993,
  \qquad
  \frac{\delta\widehat T_{m,d}}{\epsilon}\approx -0.116944.
  \]

So the mouth-side ambiguity is no longer branch choice; it is one explicit rigidity-kernel problem.

## Full-profile mouth correction

Stages 133–136 replace the tangent exponential potential by the full coupled shell + mixed D/N mouth profile. The exact residual
\[
R_*(x)=\Phi_*(x)-\Pi_*x
\]
is tangent matched,
\[
R_*(0)=0,
\qquad
R_*'(0)=0,
\]
and has negative curvature at the mouth,
\[
R_*''(0)=-3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0,
\]
so the actual full profile broadens the source relative to the tangent exponential branch.

Projecting that residual onto the Stage-130 rigidity kernel gives the actual first-order correction:
\[
\delta\mathfrak g_{\rm act}\approx -0.06480697,
\qquad
\delta\mathcal S_{\rm act}\approx -0.03887184,
\]
\[
\delta\Pi_{\rm act}\approx 0.907084,
\qquad
\delta\widehat T_{m,\rm act}\approx 0.271654.
\]
So the corrected canonical point is approximately
\[
\Pi_{\rm corr}\approx 2.41591,
\qquad
\widehat T_{m,\rm corr}\approx 1.17314.
\]
A one-step nonlinear iterate shifts slightly further but in the same direction.

## Full core–mouth co-evolution

Stages 137–140 promote the corrected mouth layer to a fully co-evolving fixed point,
\[
\Sigma(x)=\frac{e^{-\Phi_{\Sigma_0}[\Sigma](x)}}{\int_0^1 e^{-\Phi_{\Sigma_0}[\Sigma](y)}dy},
\qquad
\Phi_{\Sigma_0}[\Sigma](x)=\Sigma_0\bigl[\mathcal T_s[\Sigma](x)-\mathcal R[\Sigma]\mathcal T_q[\Sigma](x)\bigr],
\]
with
\[
\mathcal R[\Sigma]=\frac{(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

At the old canonical traction \(\Sigma_0^*\), the fixed point stays close in bias but drifts off exact compensation,
\[
\mathfrak g_{\rm fp}\approx 0.69336,
\qquad
\mathcal R_{\rm fp}\approx 0.28271,
\qquad
\Pi_{\rm fp}\approx 1.48857.
\]

Restoring exact lower-branch compensation requires a unique renormalized traction. With the reduced numerical fixed-point solver used in this batch, that root is
\[
\Sigma_0^{\rm can}\approx 4.65077,
\qquad
\widehat T_{m,\rm can}\approx 1.44667,
\qquad
\Pi_{\rm can}\approx 3.87134.
\]
These values are very close to the handoff’s quoted ones; the small differences come from the reduced collocation/iteration resolution used in the executable scripts.

## Bottom line after Stage 140

Inside the explicit Family-1 closure:

1. the upper compensated branch is impossible;
2. the equal-normalized branch is singular;
3. the lower compensated branch remains the unique regular canonical branch;
4. full self-consistency preserves that branch only after a finite traction/bias renormalization.

So the mouth/core side is no longer blocked by branch ambiguity. The remaining 2.5PN/5PN theorem gap is the projection of the residual deviation from this renormalized canonical Family-1 point into the outgoing quadrupole-normalization defect.

# 5PN / Moving-Throat continuation — Stages 141–150

This batch codifies the next handwritten moving-throat block after Stage 140: the linear defect-transport, hybrid-outlet projection, bare mixed-port slippage, D/N similarity decomposition, parent-family rigidity, microscopic off-family normal coordinate, explicit log-channel reduction, exact lower-branch drift laws, four-observable bundle inversion, and the tangent-compensation theorem.

## What is newly frozen in executable form

### Stage 141
The co-evolving Family-1 point is reduced to the linear mouth/core defect ledger
\[
\delta\mathcal R = -\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_*^2}},
\]
\[
\delta M_s = \delta\Sigma_0,
\qquad
\delta M_q = -\frac14\,\delta\Sigma_0 + \frac{\Sigma_{0,*}}{\sqrt{1+\mathfrak r_*^2}}\,\delta\mathfrak g,
\]
\[
\delta\Pi
=
\left(1-\frac14\mathcal S_*\right)\delta\Sigma_0
-\frac{\Sigma_{0,*}}{4}\,\delta\mathcal S
+\frac{\Sigma_{0,*}\mathcal S_*}{\sqrt{1+\mathfrak r_*^2}}\delta\mathfrak g.
\]

### Stage 142
Projecting the defect to the compensated Robin–mixed outlet gives the exact first-order outlet defects
\[
\delta E_2 = \frac{\delta\mathcal C - 9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\]
\[
\delta E_4 = \frac{5\delta\mathcal C - 72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]
\[
\Delta_Q = \frac{\delta\mathcal C - 27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)}.
\]
Imposing the canonical-even gate \(\delta E_2=\delta E_4=0\) yields
\[
\delta\mathcal C = 0,
\qquad
\delta\kappa_W = 0,
\qquad
\delta\mathfrak g = 0,
\]
and therefore
\[
\Delta_Q = -\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
\]

### Stage 143
The concrete core algebra collapses the remaining odd defect to the bare mixed-port slippage scalar
\[
\delta\mathfrak B_W := \delta\gamma_0 - \frac13\,\delta\kappa_0,
\qquad
\delta\gamma_W = \frac{\delta\mathfrak B_W}{1+r_{c,*}}.
\]
With the tangential susceptibility \(\delta\mathfrak B_W = \Upsilon_\Pi\,\delta\Pi_{\rm tan}\),
\[
\Delta_Q
=
-\frac{9\sigma_*\,\Upsilon_\Pi}{(1-\sigma_*)(1+r_{c,*})}\,\delta\Pi_{\rm tan}.
\]

### Stage 144
The black-box susceptibility is decomposed into the D/N similarity-slippage scalar
\[
\Xi_{\rm slip}:=\Xi_\gamma - 2\Xi_L,
\qquad
\Upsilon_\Pi = \frac{1+r_{c,*}}{9}\,\Xi_{\rm slip},
\]
so the reduced defect becomes
\[
\Delta_Q = -\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]
If the D/N similarity law
\[
\gamma_0 = \frac{4L_W^2}{3\pi^2 a^2}
\]
is preserved to first order, then \(\Xi_{\rm slip}=0\).

### Stage 145
On the exact parent compensation family
\[
1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2,
\qquad
\frac{L_W}{a} = \frac{\pi}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
\qquad
\gamma_0 = \frac{1+\mathfrak r^2}{9},
\]
automatic similarity preservation is exact:
\[
\delta\ln\gamma_0 - 2\,\delta\ln(L_W/a) = 0.
\]
On the lower branch, \(\delta\mathfrak g=0\) implies \(\delta\mathfrak r=0\), so every first-order similarity defect vanishes and
\[
\Delta_Q = 0.
\]

### Stage 146
The exact off-family normal coordinate is
\[
\delta_\perp := \delta\mathfrak g - \mathfrak g'_-(\mathfrak r_*)\,\delta\mathfrak r,
\]
with
\[
\delta\mathcal F = 4\sqrt{1+\mathfrak r_*^2}\,\delta_\perp,
\qquad
\delta R_q = -\frac{\delta_\perp}{\sqrt{1+\mathfrak r_*^2}}.
\]
Its explicit microscopic form is
\[
\delta_\perp
=
\mathfrak g_*\,
\delta\ln\!\left(\frac{g_q K_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}\,
\delta\ln\!\left(\frac{K_s K_q}{\lambda^2}\right).
\]

### Stage 147
Those two log channels are reduced to explicit wall/BdG/Maxwell/mixed drifts. Under the carried overlap and healing-lock closures:
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
\]
\[
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
\]
So \(\delta_\perp\) becomes an explicit linear combination of
\[
(\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W).
\]

### Stage 148
The exact lower compensated branch fixes
\[
\delta\ln L_W = \delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_s^2 c_{s,w}^3}{\rho_w a^5}\right),
\]
\[
\delta\ln\mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_{s,w}^3}{\rho_w c_s^2 a^3}\right).
\]
With the frozen \(n=5\) wall-EOS branch, the irreducible microscopic drift space collapses to
\[
(\delta\ln\mathcal Z_q,\ \delta\ln\rho_w,\ \delta\ln c_s,\ \delta\ln a).
\]

### Stage 149
Those four drifts are exactly inverted by the bundle observables
\[
(\Theta_w,\ K_s,\ K_q,\ P_0),
\qquad
P_0=\frac{N_0}{D_0}.
\]
The exact inversion is
\[
\delta\ln\rho_w = \frac12\,\delta\ln\Theta_w,
\]
\[
\delta\ln a = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln c_s = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w + \frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q = \delta\ln K_q - \frac25\,\delta\ln P_0.
\]

### Stage 150
Every remaining first-order mouth/background drift is an explicit algebraic image of \((\Theta_w,K_s,K_q,P_0)\):
\[
\delta\ln c_{s,w} = \delta\ln\Theta_w,
\qquad
\delta\ln\ell = -\delta\ln\Theta_w,
\qquad
\delta\ln L_W = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{13}{8}\,\delta\ln\Theta_w,
\]
\[
\delta\ln \mathcal T_m
=
-\frac54\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{15}{8}\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_q
=
-\frac34\,\delta\ln K_s
+\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln\lambda = \frac12(\delta\ln K_s+\delta\ln K_q).
\]
The tangent-compensation theorem then holds exactly:
\[
\delta\ln r_c = 0,
\qquad
\delta\ln\mathfrak r = 0,
\qquad
\delta\ln\mathfrak g = 0,
\qquad
\delta_\perp = 0.
\]

## What this means

The remaining first-order isotropic problem is now no longer “general branch drift.” It has collapsed to a bundle-observable calculation:

\[
(\Theta_w,\ K_s,\ K_q,\ P_0)
\longrightarrow
\text{all first-order isotropic mouth/core/background drifts}.
\]

And the executable result is stronger than expected: **arbitrary first-order isotropic bundle drift is tangent to the exact compensated Family-1 parent family.**

So the next live theorem gate after Stage 150 is not first-order isotropic transport anymore. It is the first correction that escapes this closed algebra, i.e. the first **off-bundle** slippage.
# Stage 151–160 continuation notes

This batch fills the missing bridge between the Stage 150 tangent-compensation theorem and the later outgoing-load / similarity-orbit chain.

The key logical contraction is:

1. **Stage 151**: first-order off-bundle motion is not a large microscopic vector; it collapses to three scalar slippages \((\varepsilon_L,\varepsilon_v,\varepsilon_T)\), and then to one weighted scalar \(\varepsilon_\perp\) with \(\delta_\perp=-\varepsilon_\perp\).

2. **Stage 152**: pure grouped real `P2` anisotropy cannot linearly feed those scalar slippages. Its scalar feed-down begins only at quadratic order through the grouped invariant
   \[
   \mathcal I[x,y]=\frac15\,\delta x^T G_{\rm grp}\,\delta y
   =4a_x a_y+\frac45 b_x b_y.
   \]
   On the weak-axisymmetric branch this becomes
   \[
   \mathcal I[x,y]=\frac{7}{10}\epsilon^2 x^{(1)}y^{(1)}.
   \]

3. **Stage 153**: the remaining **linear** grouped problem therefore lives only in the direct outlet coefficients. It collapses to
   \[
   \mathcal K_A=\delta D_{A,2}+\frac19\,\delta D_{A,0},
   \qquad
   \mathcal G_A=\delta N_{A,0}-P_0\,\delta D_{A,0},
   \]
   together with the hidden-even compatibility relation
   \[
   \delta D_{A,4}=\frac23\,\delta D_{A,2}+\frac1{27}\,\delta D_{A,0}.
   \]

4. **Stage 154**: those two grouped obstructions are not sourced by the whole microscopic bundle independently. They decompose into:
   \[
   \mathcal K_A=\mathcal W_A-\mathcal B_A-\mathcal Z_A,
   \qquad
   \mathcal G_A=-P_0\delta K_A+P_0\delta B_{A,0}+\mathcal N_A.
   \]
   So the linear grouped-even/odd problem is driven only by wall, BdG, conservative Maxwell/mixed, and outgoing-transfer anisotropies.

5. **Stage 155**: on the canonical compensated branch, the microscopic obstruction pair is just the physical slope pair
   \[
   \mathfrak K_1=-D_0 u_2^{(1)},
   \qquad
   \mathfrak G_1=D_0 P_1.
   \]
   The direct outlet amplitudes become
   \[
   \kappa_1=-\frac{3(1-\sigma_*)}{\sigma_*}u_2^{(1)},
   \qquad
   \gamma_1=-\frac{1-\sigma_*}{9\sigma_*}\frac{P_1}{P_0}.
   \]

6. **Stage 156**: expanding the actual grouped operator moments on the weak-axisymmetric branch yields
   \[
   u_2^{(1)}=-\frac{D_{21}+D_{01}/9}{D_0},
   \qquad
   \frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
   \]
   On the even-preserving branch,
   \[
   D_{21}=-\frac{D_{01}}{9},
   \qquad
   D_{41}=-\frac{D_{01}}{27},
   \]
   so the whole remaining linear grouped normalization defect is one static loading mismatch
   \[
   \Xi_{\rm load}:=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
   \]

7. **Stage 157**: that mismatch is the weighted failure of static self-similarity relative to the wall baseline:
   \[
   \Xi_{\rm load}
   =(\delta_N-\delta_K)+\omega_B(\delta_B-\delta_K)+\omega_Z(\delta_Z-\delta_K).
   \]

8. **Stage 158**: factoring by the wall baseline yields wall-normalized shape/load variables
   \[
   B_{0,\alpha}=K\chi_\alpha^2,
   \qquad
   Z_0^{(r)}=K\Upsilon_r,
   \qquad
   N_0^{(r)}=\Lambda_r^2,
   \]
   and the outgoing-load theorem
   \[
   2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r=\delta_K
   \]
   on conservative-shape-preserving branches.

9. **Stage 159**: the outgoing load factor splits exactly into
   \[
   \frac{\Lambda_r^2}{K}
   =
   \mathcal M_r^2\frac{(1+\mathcal I_r)^2}{(1-\mathcal H_r)^2},
   \]
   with
   \[
   \mathcal M_r=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
   \quad
   \mathcal I_r=\frac{R_rG_{U,r}}{\Omega_{U,r}^2G_{W,r}},
   \quad
   \mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
   \]
   Under interference/hybridization rigidity, zero defect is the square-root mixed-leg law
   \[
   \frac{G_{W,r}}{\Omega_{W,r}^2}\propto \sqrt K.
   \]

10. **Stage 160**: on the weak-axisymmetric branch the whole outgoing-slippage bundle collapses to one scalar amplitude
    \[
    \Xi_1
    =
    \sum_r \rho_r^{(N)}
    \left[
      2\mathfrak m_r
      +\frac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r
      +\frac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r
    \right],
    \]
    with the grouped pattern
    \[
    \Xi_{\rm load}^{(20)}=\epsilon\,\Xi_1,
    \qquad
    \Xi_{\rm load}^{(21)}=\frac{\epsilon}{2}\Xi_1,
    \qquad
    \Xi_{\rm load}^{(22)}=-\epsilon\,\Xi_1,
    \]
    and the exact physical identification
    \[
    \Xi_1=\frac{P_1}{P_0}.
    \]

So after Stage 160, the remaining weak-axisymmetric grouped `2.5`PN bottleneck is no longer a diffuse outlet-bundle problem. It is just the single microscopic amplitude \(\Xi_1 = P_1/P_0\) on the actual moving-throat branch.
# 5PN / moving-throat continuation — Stage 161–170 bundle

This bundle fills the previously missing numbered continuation after Stage 160.

## What is in this bundle

The scripts are split one stage at a time so the executable numbering now matches the note chain:

- Stage 161 — outgoing-port co-loading theorem
- Stage 162 — wall-normalized transfer-shape theorem
- Stage 163 — effective transfer-shape collapse
- Stage 164 — coherent tracking-branch defect law and support-blindness
- Stage 165 — microscopic coherent-kernel slippage decomposition
- Stage 166 — exact triangular normal form
- Stage 167 — branch-invariant coordinates
- Stage 168 — direct microscopic monomials and compatibility ledger
- Stage 169 — exact microscopic similarity orbit
- Stage 170 — exact orbit–quotient closure

Each script has a paired `_output.txt` file captured from a clean run in this environment.

## Structural summary

The chain now reads:

1. `Xi_1` is the mismatch between the outgoing-weighted static transfer slope and the wall-baseline slope.
2. That mismatch is exactly twice the weighted transfer-shape slope.
3. The whole many-port problem collapses to one effective transfer shape `T_eff^2 = N_0 / K`.
4. On the coherent branch, the support ratio drops out of the weak-axisymmetric defect exactly.
5. The remaining defect depends only on microscopic mixed/outgoing placement drifts.
6. Those drifts collapse to the three branch-adapted coordinates
   `Sigma_tr`, `Sigma_nt`, `Sigma_eta`.
7. Those coordinates are the logarithmic drifts of three exact direct microscopic monomials
   `C_tr,*`, `C_nt,*`, and `epsilon_eta`.
8. Their zero-defect set is the tangent space of an exact five-parameter multiplicative similarity orbit `G_*`.
9. At the finite level, the level sets of `(C_tr,*, C_nt,*, epsilon_eta)` are exactly the `G_*` orbits.
10. So the remaining microscopic question is purely branch-selective: whether the true moving-throat branch preserves those three quotient coordinates.

## Practical continuation point

The smallest honest next theorem gate after Stage 170 is:

- compute the actual branch drift of the three quotient coordinates
  `(C_tr,*, C_nt,*, epsilon_eta)`
  from the real moving-throat PDE branch;
- equivalently, determine whether the real branch stays on a single `G_*` similarity orbit.

If it does, the coherent first-order grouped weak-axisymmetric defect vanishes automatically.# 5PN Stages 171–174 — Adiabatic Wall, Elastic Branch Selection, and Orbit Locking

These stages implement the requested bridge from the `$g-2$` track back into the moving-throat / 5PN branch-selection problem.

## Stage 171 — Adiabatic wall bundle transport

Impose the adiabatic wall constraint
\[
\delta\ln\Theta_w=0.
\]
Using the exact Stage-149/150 inversion and bundle transport laws, the isotropic drift family collapses to
\[
\delta\ln\rho_w=\delta\ln c_{s,w}=\delta\ln\ell=0,
\]
\[
\delta\ln a=\delta\ln L_W=\frac12\,\delta\ln K_s,
\]
\[
\delta\ln c_s=\frac12\,\delta\ln K_s+\frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q=\delta\ln K_q-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln v_{w0}=-\frac34\,\delta\ln K_s+\frac12\,\delta\ln K_q,
\]
\[
\delta\ln\mathcal T_m=-\frac54\,\delta\ln K_s+\frac12\,\delta\ln K_q-\frac25\,\delta\ln P_0.
\]
The parent invariants remain tangent-compensated:
\[
\delta\ln r_c=\delta\ln\mathfrak r=\delta\ln\mathfrak g=0.
\]
So the adiabatic wall removes the wall-depth / thermal-fraying isotropic drift channel, but leaves a 3-parameter isotropic family labelled by
\[
(\delta\ln K_s,\ \delta\ln K_q,\ \delta\ln P_0).
\]

## Stage 172 — Adiabatic-elastic slippage collapse

The scalar off-bundle slippages are
\[
\varepsilon_L=\delta\ln L_W-\frac12\,\delta\ln K_s,
\]
\[
\varepsilon_v=\delta\ln v_{w0}+\frac34\,\delta\ln K_s-\frac12\,\delta\ln K_q,
\]
\[
\varepsilon_T=\delta\ln\mathcal T_m+\frac54\,\delta\ln K_s-\frac12\,\delta\ln K_q+\frac25\,\delta\ln P_0.
\]
On the adiabatic wall branch these need not vanish automatically, but if we impose the stronger elastic/no-fraying rule
\[
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
then the scalar normal defect vanishes identically:
\[
\delta_\perp=-\varepsilon_\perp=0.
\]
So the adiabatic-elastic boundary condition kills the first scalar off-bundle obstruction completely.

## Stage 173 — Why adiabatic wall alone is not enough for orbit locking

Stage-169/170 says the reduced weak-axisymmetric defect vanishes iff the microscopic grouped drift is tangent to the exact similarity orbit $\mathcal G_*$, equivalently iff the three quotient coordinates are preserved:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\qquad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]
These are encoded by the rank-3 monomial-drift map
\[
M_*\,\delta\mathbf x.
\]
Stage 173 shows explicitly that `\delta\ln\Theta_w=0` by itself does **not** imply
\[
M_*\delta\mathbf x=0.
\]
There are microscopic drift directions, such as a pure `\Delta\lambda_W` or pure `\Delta c_{\eta U}` perturbation, that leave the wall-depth condition untouched but still move the quotient coordinates. So the adiabatic-wall condition removes one failure channel, but does not by itself prove orbit locking.

## Stage 174 — Adiabatic-elastic orbit theorem

Combining the Stage-172 scalar result with the Stage-166 and Stage-169/170 quotient-closure theorem gives the clean unified statement:

If we impose
\[
\delta\ln\Theta_w=0,
\qquad
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
then the scalar off-bundle source vanishes, and the remaining first-order defect is zero **iff** the branch preserves the three quotient coordinates:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\qquad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]
Equivalently,
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
M_*\delta\mathbf x=0
\iff
\delta\mathbf x\in T_{\mathrm{id}}\mathcal G_*
\iff
\text{the branch stays on a single exact }\mathcal G_*\text{ orbit.}
\]

So the requested “unified loop” result is now explicit:

- the adiabatic wall condition freezes the thermal wall channel,
- the elastic/no-fraying condition removes the scalar off-bundle obstruction,
- and the remaining branch-selection test is **exactly** whether the moving-throat branch preserves
  \((\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)\), i.e. stays on one \(\mathcal G_*\) orbit.

## Direct continuation point

The next clean theorem gate is no longer to manipulate isotropic bundle transport. That part is closed. The next gate is to compute the actual physical-branch drift vector and test whether its projection under `M_*` vanishes. In other words:
\[
M_*\delta\mathbf x\stackrel{?}=0.
\]
If yes, the adiabatic-elastic moving-throat branch is orbit-locked. If not, the failure is localized immediately into the tracking, nontracking-transfer, or dressing quotient directions.
# 5PN continuation — Stages 175–180

This block continues directly from the adiabatic-wall / adiabatic-elastic orbit theorem.
The earlier result had already reduced the live first-order branch-selection problem to the exact
Stage-170 quotient test
\[
M_*\,\delta\mathbf x = 0,
\]
with quotient coordinates
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
\]

The missing step was to turn that theorem into an **exact branch-selection compiler**:
for any candidate microscopic drift, isolate the pure-similarity part, isolate the true quotient
failure, and write the unique microscopic correction that would restore single-orbit lock.

## Stage 175 — exact orbit/quotient projectors

Using the exact Stage-170 pivot block in the dependent coordinates
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
with
\[
P_{(T,K_\eta,\mu)} = M_*[:,(T,K_\eta,\mu)],
qquad
\det P = 1+\chi_{0,*}>0,
\]
one gets exact complementary projectors
\[
Q = E P^{-1} M_*,
\qquad
O = I-Q.
\]
They satisfy
\[
Q^2=Q,
\qquad
O^2=O,
\qquad
QO=OQ=0,
\qquad
M_*O=0,
\qquad
M_*Q=M_*.
\]
So every microscopic drift splits uniquely as
\[
\Delta\mathbf x = \Delta\mathbf x_{\rm orbit}+\Delta\mathbf x_{\rm fail},
\]
with
\[
\Delta\mathbf x_{\rm orbit}\in\ker M_*,
\qquad
\Delta\mathbf x_{\rm fail}=Q\Delta\mathbf x.
\]
A particularly sharp result is that the entire quotient-failure piece has support only in the
**dependent triple**
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
not in the five free similarity directions.

## Stage 176 — observable-to-microscopic correction compiler

The Stage-166/167 observable inversion gives
\[
\delta\ln\mathfrak C_{{\rm tr},*}
= -\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}{\chi_{0,*}\delta_{U,*}}\,\Theta_1,
\]
\[
\delta\ln\mathfrak C_{{\rm nt},*}
= \Xi_1+\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}\,\Theta_1,
\]
\[
\delta\ln\epsilon_\eta
= -\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\,(\mathcal R_1+\Xi_1).
\]
Composing this with the Stage-175 projector gives the exact microscopic quotient correction
supported only on the dependent triple:
\[
\Delta_T^{(q)}=
\frac{\delta\ln\mathfrak C_{{\rm tr},*}}{1+\chi_{0,*}},
\]
\[
\Delta_{K_\eta}^{(q)}=-\delta\ln\epsilon_\eta,
\]
\[
\Delta_\mu^{(q)}
=
\delta\ln\mathfrak C_{{\rm nt},*}
-\delta\ln\epsilon_\eta
+\frac{F_*}{1+\chi_{0,*}}\,\delta\ln\mathfrak C_{{\rm tr},*}.
\]
So once the observable defect triple is known, the exact microscopic dependent-coordinate
correction needed to represent it is already fixed.

## Stage 177 — exact finite restoration operator

Because the Stage-170 fibre equations are linear in the finite log-ratios, the same projector
logic gives an exact **finite** restoration operator.
For any candidate finite drift \(\Delta\mathbf x\), define
\[
\mathbf q = M_*\Delta\mathbf x,
\qquad
\Delta\mathbf x_{\rm fail}=E P^{-1}\mathbf q,
\qquad
\Delta\mathbf x_{\rm restore}=-\Delta\mathbf x_{\rm fail}.
\]
Then
\[
M_*(\Delta\mathbf x+\Delta\mathbf x_{\rm restore})=0.
\]
So any candidate branch can be returned to a single exact \(\mathcal G_*\)-orbit by adjusting
only the dependent triple \((T_U,K_\eta,\mu_W)\).

## Stage 178 — adiabatic-elastic branch decomposition

Under the strengthened boundary rule
\[
\delta\ln\Theta_w=0,
\qquad
\varepsilon_L=\varepsilon_v=\varepsilon_T=0,
\]
the thermal wall channel and the scalar off-bundle escape are both removed. The remaining
first-order branch-selection problem is therefore purely microscopic and eight-dimensional.
For any adiabatic-elastic candidate branch drift \(\Delta\mathbf x_{AE}\), the exact split is
\[
\Delta\mathbf x_{AE} = \Delta\mathbf x_{\rm orbit}+\Delta\mathbf x_{\rm fail},
\]
with
\[
M_*\Delta\mathbf x_{\rm orbit}=0,
\qquad
M_*\Delta\mathbf x_{\rm fail}=M_*\Delta\mathbf x_{AE}.
\]
And the weak-axisymmetric observables depend only on the quotient piece:
\[
\Theta_1
= -\frac{\chi_{0,*}\delta_{U,*}}{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
\,q_1,
\]
\[
\Xi_1
= \frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}\,q_1 + q_2,
\]
\[
\mathcal R_1+\Xi_1
= -\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_3,
\]
where \(\mathbf q=M_*\Delta\mathbf x_{AE}\).
So the adiabatic-elastic branch is orbit-locked iff
\[
\Delta\mathbf x_{\rm fail}=0
\iff
M_*\Delta\mathbf x_{AE}=0.
\]

## Stage 179 — exact dependent-coordinate mismatch formulas

Comparing an arbitrary candidate branch to the exact Stage-170 single-orbit law gives the three
microscopic mismatch coordinates
\[
m_T := \Delta_T-\Delta_T^{\rm orbit},
\qquad
m_{K_\eta}:=\Delta_{K_\eta}-\Delta_{K_\eta}^{\rm orbit},
\qquad
m_\mu:=\Delta_\mu-\Delta_\mu^{\rm orbit}.
\]
These are not arbitrary. They are exactly the quotient drifts:
\[
\delta\ln\mathfrak C_{{\rm tr},*}=(1+\chi_{0,*})m_T,
\]
\[
\delta\ln\epsilon_\eta=-m_{K_\eta},
\]
\[
\delta\ln\mathfrak C_{{\rm nt},*}=m_\mu-F_*m_T-m_{K_\eta}.
\]
So the remaining dynamical theorem gap has been localized completely:
**the PDE only has to prove that the dependent microscopic coordinates follow the exact
Stage-170 orbit law generated by the five free similarity directions.**

## Stage 180 — consolidated adiabatic-elastic orbit-lock verdict

The strengthened boundary rule has now been compiled all the way down to a single exact
microscopic finish line.

- adiabatic wall freezes the thermal wall datum,
- elastic/no-fraying removes the scalar off-bundle escape,
- and the remaining first-order defect is nothing but the mismatch of the dependent triple
  \((\Delta_T,\Delta_{K_\eta},\Delta_\mu)\) from the exact single-orbit law.

Equivalently, the adiabatic-elastic moving-throat branch is locked to one exact
\(\mathcal G_*\)-orbit iff
\[
\Theta_1=\Xi_1=\mathcal R_1+\Xi_1=0,
\]
or microscopically iff
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)
\]
follow the Stage-170 orbit law.

## Bottom line after Stage 180

The theorem gap is now narrower than “solve the PDE.” It is:

> show that the completed moving-throat dynamics makes the dependent microscopic triple
> \((T_U,K_\eta,\mu_W)\) follow the exact Stage-170 orbit law on the adiabatic-elastic branch.

If that happens, the first-order weak-axisymmetric defect vanishes automatically and the branch
stays on a single exact \(\mathcal G_*\)-orbit.
# Stage 181–186 notes

This block continues the post-170 / post-180 branch-selection program in the most natural direction: it upgrades the **first-order** orbit-lock language to an **exact finite** law for the dependent microscopic triple.

## Main new result

The exact Stage-168 monomials

\[
\mathfrak C_{{\rm tr},*},\qquad
\mathfrak C_{{\rm nt},*},\qquad
\epsilon_\eta
\]

can be solved **exactly** for the dependent microscopic coordinates

\[
(T_U, K_\eta, \mu_W)
\]

once the five free microscopic coordinates

\[
(\lambda_W, c_{\eta U}, \gamma, K_U, K_W)
\]

and the invariant triple are fixed.

That gives the exact finite single-orbit law:

\[
K_\eta^{(\mathrm{orbit})} = \frac{c_{\eta U}^2}{K_U\,\epsilon_{\eta,*}},
\]

\[
T_U^{(\mathrm{orbit})}
= \frac{L^2 K_U}{\pi^2}
\left[
\frac{\mathfrak C_{{\rm tr},*}}
{(\gamma c_{\eta U}/K_U)^{1+\delta_{U,*}}}
\right]^{\!1/(1+\chi_{0,*})},
\]

\[
\mu_W^{(\mathrm{orbit})}
=
\frac{\mathfrak C_{{\rm nt},*} K_\eta^{(\mathrm{orbit})} K_W^2}{\lambda_W^2}
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W}\right)^{-E_*}
\left(\frac{\pi^2 T_U^{(\mathrm{orbit})}}{L^2 K_U}\right)^{F_*}.
\]

So the finite similarity orbit is no longer abstract: the dependent triple is an exact algebraic function of the free microscopic point and the invariant triple.

## Exact finite mismatch coordinates

For any candidate branch with the **same five free microscopic coordinates** as the orbit point, define

\[
T_U = m_T\,T_U^{(\mathrm{orbit})},\qquad
K_\eta = m_K\,K_\eta^{(\mathrm{orbit})},\qquad
\mu_W = m_\mu\,\mu_W^{(\mathrm{orbit})}.
\]

Then the exact invariant ratios are

\[
\frac{\mathfrak C_{{\rm tr},*}}{\mathfrak C_{{\rm tr},*}^{(\mathrm{orbit})}} = m_T^{1+\chi_{0,*}},
\qquad
\frac{\epsilon_\eta}{\epsilon_{\eta,*}} = \frac{1}{m_K},
\qquad
\frac{\mathfrak C_{{\rm nt},*}}{\mathfrak C_{{\rm nt},*}^{(\mathrm{orbit})}} = \frac{m_\mu}{m_K m_T^{F_*}}.
\]

So the finite branch-selection problem is exactly three-dimensional.

## Exact logarithmic chart

If we write

\[
\tau := \ln m_T,
\qquad
\kappa := \ln m_K,
\qquad
\mu := \ln m_\mu,
\]

then the quotient coordinates are **exactly**

\[
q_{\rm tr} = (1+\chi_{0,*})\tau,
\qquad
q_\eta = -\kappa,
\qquad
q_{\rm nt} = \mu - \kappa - F_*\tau.
\]

So the Stage-179 first-order formulas are not merely infinitesimal approximations; they are the exact logarithmic chart of the finite mismatch ratios.

## Exact restoration map

Given the finite quotient coordinates, the exact restoration to the same orbit is achieved by changing only the dependent triple:

\[
T_U^{(\mathrm{restore})} = T_U\,e^{-q_{\rm tr}/(1+\chi_{0,*})},
\]

\[
K_\eta^{(\mathrm{restore})} = K_\eta\,e^{q_\eta},
\]

\[
\mu_W^{(\mathrm{restore})}
= \mu_W\,e^{-q_{\rm nt}+q_\eta-F_* q_{\rm tr}/(1+\chi_{0,*})}.
\]

This returns the candidate branch to the exact orbit with the same free microscopic coordinates and the same invariant triple.

## Finite adiabatic-elastic orbit-lock criterion

Under the adiabatic-elastic boundary rule, the exact finite branch-selection problem is:

\[
\text{orbit lock}
\iff
m_T = m_K = m_\mu = 1
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

So after this block the remaining theorem gap is completely concrete: once the PDE gives the actual microscopic branch values, one can test orbit lock directly by comparing the candidate dependent triple to the exact orbit values above.
# 5PN Stages 187–192 — Finite Similarity-Orbit Action, Reference-Transport Laws, and Exact Residual Diagnostics

This block continues directly from the finite orbit interface at Stage 186.
The earlier stages had already shown two things:

1. the weak-axisymmetric zero-defect branch is exactly the finite similarity orbit \
   \(\mathcal G_*\), and
2. the full finite branch-selection problem can be written as a test on the dependent
   triple \((T_U,K_\eta,\mu_W)\).

What was still missing was the **finite transport law** itself: how the five free microscopic
coordinates move a branch along an exact orbit, how to compare an actual candidate branch to a
reference orbit point, and how to localize any failure into exact multiplicative residuals.

## Stage 187 — exact finite similarity-orbit action

The five-parameter multiplicative similarity orbit \(\mathcal G_*\) is now written explicitly as a
finite action on the full microscopic state. If
\[
(\lambda_W,c_{\eta U},\gamma,K_U,K_W)
\to
(e^{\Lambda}\lambda_W,e^C c_{\eta U},e^{\Gamma}\gamma,e^U K_U,e^W K_W),
\]
then the dependent triple coevolves by
\[
K_\eta' = e^{2C-U}K_\eta,
\]
\[
T_U' = e^{U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}T_U,
\]
\[
\mu_W'
=
\exp\!\Bigl[
2C-U+2W-2\Lambda
-E_*\bigl(2\Gamma+2\Lambda-U-W\bigr)
-F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)
\Bigr]\mu_W.
\]
The three quotient monomials
\[
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta)
\]
are preserved **exactly**, not only infinitesimally.

## Stage 188 — exact group law, inverse, and parameter recovery

The finite orbit action is an exact abelian five-parameter group:
composition adds the five generator exponents, and inversion negates them. If two states are on
one orbit, the orbit parameters are recovered exactly from the free-coordinate ratios,
\[
\Lambda=\ln\frac{\lambda_W'}{\lambda_W},
\qquad
C=\ln\frac{c_{\eta U}'}{c_{\eta U}},
\qquad
\Gamma=\ln\frac{\gamma'}{\gamma},
\qquad
U=\ln\frac{K_U'}{K_U},
\qquad
W=\ln\frac{K_W'}{K_W}.
\]
So once those five free-coordinate ratios are known, the dependent triple on the same orbit is
predicted uniquely.

## Stage 189 — exact reference-orbit transport laws

Relative to a reference orbit point, the exact dependent-coordinate transport laws are
\[
R_{K_\eta}^{(\mathrm{orbit})}=\frac{R_c^2}{R_U},
\]
\[
R_{T_U}^{(\mathrm{orbit})}=R_U\left(\frac{R_U}{R_\gamma R_c}\right)^{\frac{1+\delta_{U,*}}{1+\chi_{0,*}}},
\]
\[
R_{\mu_W}^{(\mathrm{orbit})}
=
\frac{R_{K_\eta}^{(\mathrm{orbit})}R_W^2}{R_\lambda^2}
\left(\frac{R_\gamma^2R_\lambda^2}{R_UR_W}\right)^{-E_*}
\left(\frac{R_{T_U}^{(\mathrm{orbit})}}{R_U}\right)^{F_*}.
\]
This is the exact finite coevolution law of the dependent triple along a fixed orbit.

## Stage 190 — exact dependent residual coordinates

A general candidate branch with the same five free-coordinate ratios can be factored as
\[
R_{T_U}^{(\mathrm{actual})}=R_{T_U}^{(\mathrm{orbit})}m_T,
\qquad
R_{K_\eta}^{(\mathrm{actual})}=R_{K_\eta}^{(\mathrm{orbit})}m_K,
\qquad
R_{\mu_W}^{(\mathrm{actual})}=R_{\mu_W}^{(\mathrm{orbit})}m_\mu,
\]
where \((m_T,m_K,m_\mu)\) is the **dependent residual mismatch triple**.
The invariant ratios then collapse exactly to
\[
\frac{\mathfrak C_{{\rm tr},*}^{\mathrm{actual}}}{\mathfrak C_{{\rm tr},*}^{\mathrm{ref}}}=m_T^{1+\chi_{0,*}},
\qquad
\frac{\epsilon_\eta^{\mathrm{actual}}}{\epsilon_{\eta}^{\mathrm{ref}}}=\frac{1}{m_K},
\qquad
\frac{\mathfrak C_{{\rm nt},*}^{\mathrm{actual}}}{\mathfrak C_{{\rm nt},*}^{\mathrm{ref}}}=\frac{m_\mu}{m_Km_T^{F_*}}.
\]
So the free-coordinate transport along an orbit drops out completely; the quotient coordinates are
nothing but the logarithmic chart of the dependent residual triple.

## Stage 191 — factorized actual-branch interface

The actual candidate branch now admits an exact factorization
\[
\text{actual branch}
=
(\text{reference orbit point})
\times
(\text{free-coordinate orbit transport})
\times
(\text{dependent residual mismatch}).
\]
Restoration to the same orbit at fixed free-coordinate ratios is achieved by dividing only by the
residual mismatch ratios:
\[
K_\eta^{(\mathrm{restore})}=\frac{K_\eta^{(\mathrm{actual})}}{m_K},
\qquad
T_U^{(\mathrm{restore})}=\frac{T_U^{(\mathrm{actual})}}{m_T},
\qquad
\mu_W^{(\mathrm{restore})}=\frac{\mu_W^{(\mathrm{actual})}}{m_\mu}.
\]
So orbit lock is exactly the statement
\[
m_T=m_K=m_\mu=1.
\]

## Stage 192 — diagnostic signatures of each failure channel

The three quotient coordinates now have a direct physical interpretation:
\[
q_{\rm tr}=(1+\chi_{0,*})\ln m_T,
\qquad
q_\eta=-\ln m_K,
\qquad
q_{\rm nt}=\ln m_\mu-\ln m_K-F_*\ln m_T.
\]
That gives three clean signatures:

- pure \(T_U\) residual mismatch turns on \(q_{\rm tr}\) and, via \(F_*\), also \(q_{\rm nt}\),
- pure \(K_\eta\) residual mismatch turns on \(q_\eta\) and also \(q_{\rm nt}\),
- pure \(\mu_W\) residual mismatch turns on \(q_{\rm nt}\) only.

So once the PDE delivers an actual candidate branch, the pattern of
\[
(q_{\rm tr},q_\eta,q_{\rm nt})
\]
identifies exactly which dependent coevolution law failed, if any.

## Bottom line after Stage 192

The branch-selection problem is no longer merely
“compute the actual branch and see if the quotient coordinates move.”
It is now a **factorized finite comparison problem**:

1. choose a reference point on the exact similarity orbit,
2. use the five free-coordinate ratios to predict the orbit-transported dependent triple,
3. compare the actual dependent triple to that prediction,
4. read off the residual mismatch ratios \((m_T,m_K,m_\mu)\),
5. and infer the quotient coordinates and restoration map immediately.

So the next PDE theorem gate is now even sharper:

> compute the actual branch values of the eight microscopic coordinates, form the free-coordinate
> ratios and the dependent residual mismatch triple, and test whether
> \(m_T=m_K=m_\mu=1\).

If yes, the branch stays on a single exact \(\mathcal G_*\) orbit. If not, the failure is localized
exactly into the \(T_U\), \(K_\eta\), and/or \(\mu_W\) transport laws.

# Stages 193–198 — Pairwise orbit-lock, cocycle laws, and minimal-data verdict

This block continues the finite orbit-lock chain after Stages 181–192 by making the
orbit criterion fully **reference-independent** and **PDE-ready**.

## Stage 193 — exact pairwise orbit criterion

Given two positive microscopic states `x` and `y`, with shared branch constants
`(chi0_*, deltaU_*, E_*, F_*)`, the five free-coordinate ratios
\[
R_\lambda,\ R_c,\ R_\gamma,\ R_U,\ R_W
\]
determine the exact orbit-predicted dependent-coordinate ratios
\[
R_{K_\eta}^{(\mathrm{orbit})}=\frac{R_c^2}{R_U},
\]
\[
R_{T_U}^{(\mathrm{orbit})}
=
R_U\left(\frac{R_U}{R_\gamma R_c}\right)^{\frac{1+\delta_{U,*}}{1+\chi_{0,*}}},
\]
\[
R_{\mu_W}^{(\mathrm{orbit})}
=
\frac{R_{K_\eta}^{(\mathrm{orbit})}R_W^2}{R_\lambda^2}
\left(\frac{R_\gamma^2R_\lambda^2}{R_UR_W}\right)^{-E_*}
\left(\frac{R_{T_U}^{(\mathrm{orbit})}}{R_U}\right)^{F_*}.
\]

Comparing with the actual dependent-coordinate ratios defines the pairwise residual triple
\[
m_T^{(x\to y)}=\frac{R_{T_U}^{(\mathrm{act})}}{R_{T_U}^{(\mathrm{orbit})}},
\qquad
m_K^{(x\to y)}=\frac{R_{K_\eta}^{(\mathrm{act})}}{R_{K_\eta}^{(\mathrm{orbit})}},
\qquad
m_\mu^{(x\to y)}=\frac{R_{\mu_W}^{(\mathrm{act})}}{R_{\mu_W}^{(\mathrm{orbit})}}.
\]

The invariant ratios between `y` and `x` are then exactly
\[
\frac{\mathfrak C_{{\rm tr},*}(y)}{\mathfrak C_{{\rm tr},*}(x)}
=
\bigl(m_T^{(x\to y)}\bigr)^{1+\chi_{0,*}},
\]
\[
\frac{\epsilon_\eta(y)}{\epsilon_\eta(x)}
=
\frac{1}{m_K^{(x\to y)}},
\]
\[
\frac{\mathfrak C_{{\rm nt},*}(y)}{\mathfrak C_{{\rm nt},*}(x)}
=
\frac{m_\mu^{(x\to y)}}{m_K^{(x\to y)}\bigl(m_T^{(x\to y)}\bigr)^{F_*}}.
\]

So `x` and `y` lie on the same exact similarity orbit iff
\[
m_T^{(x\to y)}=m_K^{(x\to y)}=m_\mu^{(x\to y)}=1,
\]
equivalently iff the three invariant ratios are all unity.

## Stage 194 — multiplicative cocycle and additive quotient law

For three states `x,y,z`, the residual ratios compose exactly:
\[
m_T^{(x\to z)}=m_T^{(x\to y)}m_T^{(y\to z)},
\qquad
m_K^{(x\to z)}=m_K^{(x\to y)}m_K^{(y\to z)},
\qquad
m_\mu^{(x\to z)}=m_\mu^{(x\to y)}m_\mu^{(y\to z)}.
\]

The logarithmic quotient coordinates therefore add:
\[
q_{\rm tr}^{(x\to z)}=q_{\rm tr}^{(x\to y)}+q_{\rm tr}^{(y\to z)},
\]
\[
q_{\rm nt}^{(x\to z)}=q_{\rm nt}^{(x\to y)}+q_{\rm nt}^{(y\to z)},
\]
\[
q_\eta^{(x\to z)}=q_\eta^{(x\to y)}+q_\eta^{(y\to z)},
\]
with inverse laws under reversal,
\[
q^{(y\to x)}=-q^{(x\to y)}.
\]

So a sequence of PDE branch snapshots can be tracked either multiplicatively in the
residual ratios or additively in the quotient coordinates.

## Stage 195 — pairwise quotient-to-observable compiler

From the residual ratios,
\[
q_{\rm tr}=(1+\chi_{0,*})\ln m_T,\qquad
q_\eta=-\ln m_K,\qquad
q_{\rm nt}=\ln m_\mu-\ln m_K-F_*\ln m_T.
\]

Composing with the Stage-170 linear observable map gives the small pairwise observable
signature
\[
\Theta_1^{(\mathrm{lin})}
=
-\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
\,q_{\rm tr},
\]
\[
\Xi_1^{(\mathrm{lin})}
=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}\,q_{\rm tr}+q_{\rm nt},
\]
\[
(\mathcal R_1+\Xi_1)^{(\mathrm{lin})}
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,q_\eta.
\]

The pure-channel signatures become:

- pure `T_U` mismatch:
  \[
  q_{\rm tr}=(1+\chi_{0,*})\ln z_T,\qquad
  q_{\rm nt}=-F_*\ln z_T,\qquad
  q_\eta=0;
  \]

- pure `K_\eta` mismatch:
  \[
  q_{\rm tr}=0,\qquad
  q_{\rm nt}=-\ln z_K,\qquad
  q_\eta=-\ln z_K;
  \]

- pure `\mu_W` mismatch:
  \[
  q_{\rm tr}=0,\qquad
  q_{\rm nt}=\ln z_M,\qquad
  q_\eta=0.
  \]

This is the finite pairwise uplift of the earlier Stage-192 channel diagnostics.

## Stage 196 — exact inversion from invariant ratios

The invariant-ratio packet
\[
(R_{\rm tr},R_{\rm nt},R_\eta)
:=
\left(
\frac{\mathfrak C_{{\rm tr},*}(y)}{\mathfrak C_{{\rm tr},*}(x)},
\frac{\mathfrak C_{{\rm nt},*}(y)}{\mathfrak C_{{\rm nt},*}(x)},
\frac{\epsilon_\eta(y)}{\epsilon_\eta(x)}
\right)
\]
already contains the full orbit-lock verdict.

The exact inversion is
\[
m_T = R_{\rm tr}^{1/(1+\chi_{0,*})},
\qquad
m_K = \frac{1}{R_\eta},
\qquad
m_\mu = R_{\rm nt}\,m_K\,m_T^{F_*}.
\]

So once the invariant ratios are known, the five free-coordinate ratios are not needed
for the verdict or for the dependent-coordinate restoration.

## Stage 197 — canonical orbit-distance quadratic form

Write the logarithmic residuals as
\[
t=\ln m_T,\qquad k=\ln m_K,\qquad \mu=\ln m_\mu.
\]

Then
\[
\begin{pmatrix}
q_{\rm tr}\\ q_{\rm nt}\\ q_\eta
\end{pmatrix}
=
A
\begin{pmatrix}
t\\k\\\mu
\end{pmatrix},
\qquad
A=
\begin{pmatrix}
1+\chi_{0,*} & 0 & 0\\
-F_* & -1 & 1\\
0 & -1 & 0
\end{pmatrix}.
\]

The canonical scalar orbit-distance is
\[
D^2=q_{\rm tr}^2+q_{\rm nt}^2+q_\eta^2
=
\begin{pmatrix}
t&k&\mu
\end{pmatrix}
Q
\begin{pmatrix}
t\\k\\\mu
\end{pmatrix},
\qquad
Q=A^TA,
\]
i.e.
\[
Q=
\begin{pmatrix}
(1+\chi_{0,*})^2+F_*^2 & F_* & -F_*\\
F_* & 2 & -1\\
-F_* & -1 & 1
\end{pmatrix}.
\]

Its principal minors are
\[
(1+\chi_{0,*})^2+F_*^2>0,
\]
\[
2(1+\chi_{0,*})^2+F_*^2>0,
\]
\[
(1+\chi_{0,*})^2>0,
\]
so `Q` is positive definite on the constructive branch.

Therefore
\[
D^2=0
\iff
m_T=m_K=m_\mu=1.
\]

So the entire pairwise orbit-lock failure can be summarized by one exact reference-free
positive scalar.

## Stage 198 — minimal-data orbit verdict

The full orbit-lock verdict can be reached from **any one** of three exact packets:

1. residual mismatch ratios:
   \[
   (m_T,m_K,m_\mu),
   \]

2. invariant ratios:
   \[
   (R_{\rm tr},R_{\rm nt},R_\eta),
   \]

3. quotient coordinates:
   \[
   (q_{\rm tr},q_{\rm nt},q_\eta).
   \]

They are exactly interconvertible:
\[
R_{\rm tr}=m_T^{1+\chi_{0,*}},\qquad
R_{\rm nt}=\frac{m_\mu}{m_K m_T^{F_*}},\qquad
R_\eta=\frac{1}{m_K},
\]
\[
q_{\rm tr}=\ln R_{\rm tr},\qquad
q_{\rm nt}=\ln R_{\rm nt},\qquad
q_\eta=\ln R_\eta,
\]
\[
m_T=e^{q_{\rm tr}/(1+\chi_{0,*})},\qquad
m_K=e^{-q_\eta},\qquad
m_\mu=e^{q_{\rm nt}-q_\eta+F_*q_{\rm tr}/(1+\chi_{0,*})}.
\]

So future PDE numerics only need to provide whichever packet is cleanest. From that
packet one can reconstruct the dependent-coordinate restoration map and the scalar
distance `D^2`.

# Stages 199–201 — Bringing the 5PN endgame home

This block stops the long exploratory chain and compresses the remaining theorem gap
into the smallest exact packets the completed moving-throat PDE must still supply.

## Stage 199 — exact final branch residual packet

Take the actual grouped-lane low-frequency bundle data
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad A\in\{20,21,22\},
\]
together with the source-map factor \(m_{\hat 0}\).

Compile the normalized grouped response moments
\[
u_2^{(A)}=-\frac{D_{A2}}{D_{A0}},
\qquad
u_4^{(A)}=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2},
\]
and the outgoing prefactor moments
\[
P_0^{(A)}=\frac{N_{A0}}{D_{A0}},
\]
\[
P_2^{(A)}=\frac{D_{A0}N_{A2}-2D_{A2}N_{A0}}{D_{A0}^2},
\]
\[
P_4^{(A)}=
\frac{D_{A0}^2N_{A4}-2D_{A0}(D_{A2}N_{A2}+D_{A4}N_{A0})+3D_{A2}^2N_{A0}}
{D_{A0}^3}.
\]

Then extract the grouped trace/anomaly data
\[
(\bar u_2,a_2,b_2),\qquad
(\bar u_4,a_4,b_4),\qquad
(\bar P_0,a_{P_0},b_{P_0}).
\]

The exact final branch residual packet is
\[
\Delta_{\rm branch}
=
\bigl(
a_2,\ b_2,\ a_4,\ b_4,\ a_{P_0},\ b_{P_0},\ \Delta_{\rm pole},\ \Delta_{\rm norm}
\bigr),
\]
with
\[
\Delta_{\rm pole}=\bar u_4-4\bar u_2^{\,2},
\]
\[
\Delta_{\rm norm}=m_{\hat 0}^{\,2}\bar P_0-\frac{54Gc_s^5}{5a^5c^5}.
\]

So the completed PDE no longer has to “show 5PN somehow.” It has to drive one exact
finite-dimensional residual packet to zero.

## Stage 200 — exact endgame compiler

The orbit side was already reduced in Stages 181–198 to any one of the equivalent packets
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

This stage combines that orbit packet with \(\Delta_{\rm branch}\) and shows that the
whole reduced 5PN / 2.5PN / 4PN closure problem depends only on

1. the grouped branch packet \(\Delta_{\rm branch}\),
2. the orbit-lock packet \(\Delta_{\rm orbit}\).

It also records one useful practical simplification: on the minimal isotropic conservative
module, the explicit Family‑1 support/source side is already above the required threshold,
because
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
A_{F1}\approx 1.00005192880220 > \frac13.
\]
So the support/source side is no longer the active bottleneck inside the current hierarchy.

The reduced closure problem is therefore exactly
\[
\Delta_{\rm branch}=0,
\qquad
\Delta_{\rm orbit}=0.
\]

## Stage 201 — home-stretch theorem and minimal PDE data packet

The final theorem gate is now stated as a minimal-data problem.

### Packet A — grouped bundle data
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad
A\in\{20,21,22\},
\]
plus \(m_{\hat 0}\).

### Packet B — orbit/invariant data
Any one of
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

Everything else in the reduced closure test is an exact compiler output of these packets.

So the remaining theorem gap is no longer diffuse. The completed moving-throat PDE only
has to supply the data needed to evaluate
\[
\Delta_{\rm branch}
\quad\text{and}\quad
\Delta_{\rm orbit}.
\]

If \(\Delta_{\rm branch}\neq0\), the branch fails the reduced GR test.
If \(\Delta_{\rm branch}=0\) but \(\Delta_{\rm orbit}\neq0\), the branch is isotropic
but off the exact similarity orbit.
Only if both vanish is the reduced closure complete inside the current hierarchy.
# 5PN continuation — normalized Packet-A / Packet-B bridge from the explicit overlap model

## What this session did

The previous session turned the Stage 199–201 endgame into working code for the **branch side** of the 5PN problem:

- Stage 202: exact Packet-A compiler,
- Stage 203: minimal isotropic single-mode overlap bridge,
- Stage 204: linear grouped outlet obstructions.

The main gap after that was structural rather than algebraic:

> Packet A and Packet B were still living in different microscopic languages.

Packet A was phrased in the explicit overlap variables
\((K,M,C,\varpi,G_U,G_W,R,\Omega_U,\Omega_W)\),
while the later weak-axisymmetric / similarity-orbit notes were phrased in the normalized coherent variables
\((\chi_0,\epsilon_\eta,Z_W,\delta_U,\dots)\).

This session closes that gap.

I built two new exact bridges:

1. a **normalized coherent isotropic Packet-A bridge**,
2. a **normalized coherent Packet-B / orbit bridge**.

So the explicit finite-throat overlap model and the later normalized monomial/orbit theorem now talk to each other directly.

---

## 1. Stage 205 — normalized coherent isotropic Packet-A bridge

### 1.1 Normalized coherent dictionary

Starting from the minimal isotropic overlap model, define the normalized coherent variables

\[
\chi_0 = \frac{R G_U}{\Omega_U^2 G_W},
\qquad
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2},
\qquad
Z_W = \frac{M G_W^2}{K\Omega_W^2}.
\]

Then the microscopic couplings are reconstructed exactly as

\[
G_U = \Omega_U\sqrt{\frac{\epsilon_\eta K}{M}},
\qquad
G_W = \Omega_W\sqrt{\frac{Z_W K}{M}},
\qquad
R = \chi_0\,\Omega_U\Omega_W\sqrt{\frac{Z_W}{\epsilon_\eta}}.
\]

So the isotropic overlap packet can now be written in the same normalized coherent variables used later in the 5PN notes.

---

### 1.2 The exact mixed-sector control parameter

The Stage-205 script shows that the conservative/outgoing mixed block enters through the single exact combination

\[
\epsilon_{\rm mix}
:=
\frac{\chi_0^2 Z_W}{\epsilon_\eta}.
\]

Then the one-port isotropic denominator becomes

\[
\Delta = \Omega_U^2\Omega_W^2\bigl(1-\epsilon_{\rm mix}\bigr).
\]

This is the cleanest exact reduction of the mixed block I have so far.

---

### 1.3 Exact isotropic conservative/outgoing coefficients in normalized form

The script gives the first exact normalized formulas for the isotropic Packet-A front end.

The BdG support moments remain

\[
B_0=\frac{C^2}{\varpi^2},
\qquad
B_2=\frac{C^2}{\varpi^4},
\qquad
B_4=\frac{C^2}{\varpi^6}.
\]

The mixed conservative static slot becomes

\[
Z_0
=
\frac{K}{M}
\frac{\epsilon_\eta + Z_W(1+2\chi_0)}{1-\epsilon_{\rm mix}}.
\]

The outgoing-transfer static slot becomes

\[
N_0
=
\frac{K}{M\Omega_W^2}
\frac{Z_W(1+\chi_0)^2}{(1-\epsilon_{\rm mix})^2}.
\]

The higher conservative mixed slots are packaged as

\[
Z_2 = \frac{K}{M}\,\Sigma_2,
\qquad
Z_4 = \frac{K}{M}\,\Sigma_4,
\]

with exact closed expressions \(\Sigma_2,\Sigma_4\) produced in the script.

So the isotropic Packet-A bundle is now explicitly controlled by

- the **support pair** \((C,\varpi)\),
- the **wall pair** \((K,M)\),
- the **transport scales** \((\Omega_U,\Omega_W)\),
- the **normalized mixed ratios** \((\chi_0,\epsilon_\eta,Z_W)\).

---

### 1.4 Exact support-blind theorem for the outgoing-transfer bundle

One of the strongest exact results of this session is:

\[
\frac{\partial N_0}{\partial C}
=
\frac{\partial N_0}{\partial \varpi}
=
\frac{\partial N_2}{\partial C}
=
\frac{\partial N_2}{\partial \varpi}
=
\frac{\partial N_4}{\partial C}
=
\frac{\partial N_4}{\partial \varpi}
=0.
\]

So on the isotropic branch, the entire outgoing-transfer packet
\((N_0,N_2,N_4)\)
is **exactly support-blind** in the explicit BdG pair \((C,\varpi)\).

That does **not** mean the full Packet-A verdict is support-blind, because

\[
D_0 = K - B_0 - Z_0
\]

still depends on \(B_0=C^2/\varpi^2\), and similarly for \(D_2,D_4\).

So the correct separation is:

- the **outgoing transfer side** is support-blind,
- the **conservative wall operator** is not.

This is precisely the kind of structural split we needed to make the roadmap honest.

---

### 1.5 Exact normalized compatibility surface

Combining the 5PN one-pole gate with the 2.5PN/4PN normalization target still gives the same exact scalar compatibility law,

\[
\frac{N_0}{P_{0,\rm target}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4},
\qquad
P_{0,\rm target}=
\frac{54Gc_s^5}{5a^5c^5m_{\hat 0}^{\,2}}.
\]

But now it is written in normalized coherent language:

\[
\frac{K}{M\Omega_W^2}
\frac{Z_W(1+\chi_0)^2}{(1-\epsilon_{\rm mix})^2}
\frac{1}{P_{0,\rm target}}
=
\frac{3\left(M+\frac{C^2}{\varpi^4}+\frac{K}{M}\Sigma_2\right)^2}
{\frac{C^2}{\varpi^6}+\frac{K}{M}\Sigma_4}.
\]

That is the first exact normalized form of the isotropic Packet-A theorem gate.

So the remaining isotropic Packet-A problem is no longer “solve all overlap data somehow.”
It is:

> compute the scalar set
> \((K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,Z_W,m_{\hat0})\)
> on the physical branch.

Then the entire isotropic Packet-A verdict is immediate.

---

## 2. Stage 206 — normalized coherent Packet-B / orbit bridge

### 2.1 Exact normalized monomials

The later weak-axisymmetric/similarity-orbit notes were encoded in the direct monomials

\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta.
\]

This session rewrites them directly in the normalized coherent variables.

With

\[
\chi_0 = \frac{R G_U}{\Omega_U^2 G_W},
\qquad
\epsilon_W = \frac{R^2\sigma}{\Omega_U^2\Omega_W^2},
\qquad
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2},
\]

the exact monomials are

\[
\mathfrak C_{\rm tr}
=
\chi_0^{1+\delta_{U,*}}\,
\delta_U^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{\rm nt}
=
\frac{M G_W^2}{K\Omega_W^4}
\epsilon_W^{E_*}
\delta_U^{-F_*},
\]

\[
\mathfrak C_\eta = \epsilon_\eta.
\]

So Packet B is now coded directly in the same microscopic language as Packet A.

---

### 2.2 Exact invariant packet and quotient coordinates

Relative to a reference branch,

\[
R_{\rm tr} = \frac{\mathfrak C_{\rm tr}}{\mathfrak C_{{\rm tr},\rm ref}},
\qquad
R_{\rm nt} = \frac{\mathfrak C_{\rm nt}}{\mathfrak C_{{\rm nt},\rm ref}},
\qquad
R_\eta = \frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

Then the exact quotient-coordinate packet is just

\[
q_{\rm tr} = \ln R_{\rm tr},
\qquad
q_{\rm nt} = \ln R_{\rm nt},
\qquad
q_\eta = \ln R_\eta.
\]

The new script verifies that this packet roundtrips exactly through the existing Stage-200/201 common compiler.

So Packet B is now operational in normalized coherent variables rather than only in abstract quotient notation.

---

### 2.3 Exact normalized monomial-drift matrix

The script then re-derives the direct monomial-drift matrix from first principles in the normalized variable set

\[
(d\ln G_W,
 d\ln G_U,
 d\ln R,
 d\ln K,
 d\ln M,
 d\ln\Omega_U,
 d\ln\Omega_W,
 d\ln\delta_U).
\]

The exact result is

\[
M_{\rm norm}=
\begin{pmatrix}
-(1+\delta_{U,*}) & 1+\delta_{U,*} & 1+\delta_{U,*} & 0 & 0 & -2(1+\delta_{U,*}) & 0 & 1+\chi_{0,*} \\
2 & 0 & 2E_* & -1 & 1 & -2E_* & -(4+2E_*) & -F_* \\
0 & 2 & 0 & -1 & 1 & -2 & 0 & 0
\end{pmatrix}.
\]

The script verifies this exactly against the Stage-13 matrix.

So the Stage-13 normalized similarity-orbit theorem is now implemented as code rather than only as notes.

---

### 2.4 Exact triangular zero-defect solve

Setting the three monomial drifts to zero and solving the system gives the same triangular co-drift laws as the notes.

The script recovers exactly:

\[
d\ln\delta_U
=
-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
\bigl(d\ln R + d\ln G_U - d\ln G_W - 2d\ln\Omega_U\bigr),
\]

\[
d\ln M
=
d\ln K - 2d\ln G_U + 2d\ln\Omega_U,
\]

\[
d\ln\Omega_W
=
\frac{d\ln G_W - d\ln G_U + (1-E_*)d\ln\Omega_U + E_* d\ln R - \tfrac{F_*}{2}d\ln\delta_U}{E_*+2}.
\]

So the zero-defect tangency problem is now executable:

> once the actual branch gives drifts of
> \((G_W,G_U,R,K,\Omega_U)\),
> the required co-drifts of
> \((\delta_U,M,\Omega_W)\)
> are fixed exactly.

That is the cleanest current reduced theorem gate on the Packet-B side.

---

## 3. Combined reading of Stages 205 and 206

The main gain of this session is that Packet A and Packet B now share one microscopic vocabulary.

### Shared microscopic sector

Both packets depend on the same normalized coherent core:

\[
(K,M,G_U,G_W,R,\Omega_U,\Omega_W).
\]

Equivalently, after normalization,

\[
(K,M,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,Z_W).
\]

### Packet-A-only data

Packet A still needs the explicit support pair

\[
(C,\varpi),
\]

because the conservative wall operator depends on

\[
B_0,B_2,B_4.
\]

### Packet-B-only data

Packet B needs the tracking variable

\[
\delta_U,
\]

plus the reference-branch exponents and calibration constants

\[
(\chi_{0,*},\delta_{U,*},E_*,F_*).
\]

So the full reduced 5PN data request has now sharpened further.

The completed moving-throat branch does **not** need to hand us unrelated large symbolic objects.
It needs to hand us:

1. the shared normalized coherent core,
2. the support pair \((C,\varpi)\),
3. the tracking variable \(\delta_U\),
4. and the source factor \(m_{\hat0}\).

From that, the present code can now produce:

- the isotropic Packet-A residual,
- the orbit Packet-B residual,
- and the exact zero-defect tangency conditions.

---

## 4. Best continuation point after this session

The next clean move is now very explicit.

### Step 1 — extract isotropic branch scalars
From the moving-throat branch, compute

\[
K,
\ M,
\ C,
\ \varpi,
\ \Omega_U,
\ \Omega_W,
\ \chi_0,
\ \epsilon_\eta,
\ Z_W,
\ \delta_U,
\ m_{\hat0}.
\]

### Step 2 — run Packet A
Feed those into Stage 205 to get

- \(B_n,Z_n,N_n,D_n\),
- \(u_2,u_4,P_0\),
- \(\Delta_{\rm pole},\Delta_{\rm norm}\).

### Step 3 — run Packet B
Feed the same normalized core plus \(\delta_U\) into Stage 206 to get

- \((R_{\rm tr},R_{\rm nt},R_\eta)\),
- \((q_{\rm tr},q_{\rm nt},q_\eta)\),
- and the exact zero-defect tangency test.

### Step 4 — only then reopen anisotropy
Once the isotropic packet and orbit packet are both computed, the weak-axisymmetric Stage-204 obstruction test becomes the right next refinement.

So after this session, the continuation point is no longer vague.
It is:

> compute the actual isotropic branch scalars in the shared normalized coherent language.

That is now the shortest honest route back into the full 5PN endgame.
# 5PN stages 207–208 — explicit overlap extractor and isotropic packet compiler

## What this turn accomplished

This turn pushed the Stage-205/206 normalized bridge one level deeper into the actual moving-throat overlap model.

The main new result is that the explicit finite-throat prototype already compresses to one **finite isotropic branch state**, and the two surviving endgame packets do **not** require the same pieces of that state.

That is the first honest extractor-level answer to the question “what data does the moving-throat branch really have to provide?”

---

## Stage 207 — explicit finite-throat overlap extractor

### Prototype used

I fixed the explicit finite-throat radial/axial prototype on
\[
s\in[0,L],
\]
with
\[
\chi_\eta(s)=\sqrt{\frac{2}{L}}\sin\!\frac{\pi s}{L},
\qquad
\phi_{\rm DN}(s)=\sqrt{\frac{2}{L}}\sin\!\frac{\pi s}{2L},
\]
and chose the simplest support/mixed profiles
\[
u(s)=\chi_\eta(s),
\qquad
w(s)=\phi_{\rm DN}(s).
\]

The exact overlaps are
\[
I_{\eta u}=1,
\qquad
I_{\eta\phi}=I_{\eta w}=I_{uw}=\frac{8}{3\pi}.
\]

### Overlap-renormalized couplings

If the raw microscopic inputs are
\[
(K,M,\lambda_B^{\rm raw},\varpi,c_{\eta U}^{\rm raw},\lambda_W^{\rm raw},\gamma^{\rm raw},K_U,\mu_U,K_W,\mu_W,T_U,\sigma),
\]
then the overlap extractor gives
\[
C=\frac{8}{3\pi}\lambda_B^{\rm raw},
\qquad
c_{\eta U}^{\rm eff}=c_{\eta U}^{\rm raw},
\qquad
\lambda_W^{\rm eff}=\frac{8}{3\pi}\lambda_W^{\rm raw},
\qquad
\gamma^{\rm eff}=\gamma^{\rm raw}.
\]

So in this explicit prototype the geometry renormalizes only the wall/support and wall/mixed amplitudes that actually feel the nontrivial half-wave overlap.

### Extracted isotropic branch state

The prototype then compresses to the finite isotropic state
\[
(K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,\epsilon_W,Z_W,\delta_U),
\]
with
\[
\Omega_U=\sqrt{K_U/\mu_U},
\qquad
\Omega_W=\sqrt{K_W/\mu_W},
\]
\[
\chi_0=\frac{\gamma^{\rm raw} c_{\eta U}^{\rm raw}}{K_U},
\qquad
\epsilon_\eta=\frac{(c_{\eta U}^{\rm raw})^2}{K K_U},
\]
\[
Z_W=\frac{(\lambda_W^{\rm eff})^2}{K K_W},
\qquad
\epsilon_W=\frac{(\gamma^{\rm raw})^2(\lambda_W^{\rm eff})^2\sigma}{K_U K_W},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U}.
\]

This is exactly the coherent local-kernel variable set, now obtained from an explicit overlap model rather than inserted abstractly.

### Selected-branch transfer identity

Writing
\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
the direct transfer shape is
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
\]

The selected-branch demand ratio becomes
\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\Lambda=
\frac{27\pi^2Gc_s^5K_W}{20a^5c^5\mu_W},
\]
so the exact product identity is
\[
R_{\rm target}\,\mathcal T^2
=
\Lambda_0(1-\epsilon_\eta),
\qquad
\Lambda_0=
\frac{27\pi^2Gc_s^5}{20a^5c^5}.
\]

That is the first explicit overlap-level version of the selected-branch identity.

---

## Stage 208 — isotropic Packet-A / Packet-B compiler from the extracted state

Stage 208 takes the Stage-207 state and compiles the endgame packets directly.

### Packet A on the isotropic branch

The conservative isotropic operator moments are
\[
D_0=K-B_0-Z_0,
\qquad
D_2=-(M+B_2+Z_2),
\qquad
D_4=-(B_4+Z_4),
\]
with the usual normalized branch defects
\[
\Delta_{\rm pole}
=
\bar u_4-4\bar u_2^2
=
-\frac{3D_2^2+D_0D_4}{D_0^2},
\]
\[
P_0=\frac{N_0}{D_0},
\qquad
\Delta_{\rm norm}
=
m_{\hat 0}^2P_0-rac{54Gc_s^5}{5a^5c^5}.
\]

Because the grouped lanes are already collapsed on the isotropic branch, the full Stage-199 branch packet reduces exactly to the scalar pair
\[
(\Delta_{\rm pole},\Delta_{\rm norm}).
\]

### Packet B from the same extracted state

Using the reference-branch exponents \((\chi_{0,*},\delta_{U,*},E_*,F_*)\), the exact direct monomials are
\[
\mathfrak C_{\rm tr}=\chi_0^{1+\delta_{U,*}}\,\delta_U^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{\rm nt}=\frac{Z_W}{\Omega_W^2}\,\epsilon_W^{E_*}\,\delta_U^{-F_*},
\qquad
\mathfrak C_\eta=\epsilon_\eta.
\]

So the orbit packet is compiled by
\[
(\mathfrak C_{\rm tr},\mathfrak C_{\rm nt},\mathfrak C_\eta)
\to
(R_{\rm tr},R_{\rm nt},R_\eta)
\to
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

### The main new structural theorem

The two packets separate cleanly.

**Packet A depends only on**
\[
(K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,Z_W)
\quad\text{plus}\quad m_{\hat 0},
\]
and is blind to
\[
(\epsilon_W,\delta_U).
\]

**Packet B depends only on**
\[
(\chi_0,\delta_U,\epsilon_\eta,\epsilon_W,Z_W,\Omega_W),
\]
and is blind to the explicit support pair
\[
(C,\varpi),
\]
as well as the wall inertia pair
\[
(M,\Omega_U)
\]
once the extracted invariants are formed.

So the corrected combined isotropic endgame state is
\[
(K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,\epsilon_W,Z_W,\delta_U)
\]
plus the source factor
\[
m_{\hat 0}.
\]

That is the cleanest current statement of what the actual moving-throat branch still has to provide.

---

## What changed in the roadmap

Before this turn, the next step was loosely “extract isotropic branch data from the moving-throat overlap model.”

After this turn, the next theorem gate is much sharper:

1. extract the exact 11-scalar isotropic branch state
   \[
   (K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,\epsilon_W,Z_W,\delta_U),
   \]
2. extract or compute the source factor \(m_{\hat 0}\),
3. compile
   \[
   \Delta_{\rm pole},\qquad \Delta_{\rm norm},\qquad (q_{\rm tr},q_{\rm nt},q_\eta),
   \]
4. and test whether they all vanish on the actual moving-throat branch.

That is the direct continuation point.
# 5PN continuation notes — stages 209 and 210

This session pushed the explicit isotropic overlap prototype into the first runnable
weak-axisymmetric mechanism sieve.

The clean continuation point after Stages 207–208 was to stop speaking abstractly
about “the weak-axisymmetric defect” and instead compile it directly from the same
finite-throat overlap model that already feeds Packet A and Packet B on the isotropic
branch.

The two new stages do exactly that.

---

## Stage 209 — primitive weak-axisymmetric extractor from the explicit overlap model

Files:
- `5pn_stage209_primitive_weak_axisymmetric_extractor.py`
- `5pn_stage209_primitive_weak_axisymmetric_extractor_output.txt`

### What was fixed

Start from the explicit isotropic overlap prototype already used in Stages 207–208:

- `I_cross = 8/(3 pi)`
- `C   = lambda_B I_cross`
- `G_U = lambda_U`
- `G_W = lambda_W I_cross`
- `R   = lambda_R I_cross`

and add the primitive weak-axisymmetric absolute slopes

- `dK`, `dM`
- `d(lambda_B)`, `d(varpi)`
- `d(lambda_U)`, `d(lambda_W)`, `d(lambda_R)`
- `d(Omega_U)`, `d(Omega_W)`.

From the isotropic bundle data

- `Delta = Omega_U^2 Omega_W^2 - R^2`
- `P = Omega_U^2 G_W + R G_U`
- `Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2`
- `H = G_U^2 + G_W^2`
- `S2 = Omega_U^2 + Omega_W^2`

this stage derives the exact first-order primitive variations

- `Delta01`
- `P01`
- `Q01`
- `H01`
- `S201`

and then the full slope packet

- `B01,B21,B41`
- `Z01,Z21,Z41`
- `N01`
- `D01,D21,D41`.

### Exact bundle decomposition

The script proves the exact weak-axisymmetric bundle identities

```text
D01 = dK - B01 - Z01
D21 = -(dM + B21 + Z21)
D41 = -(B41 + Z41)
Xi_load = N01/N0 - D01/D0 = P1/P0
```

so the full primitive deformation only reaches the grouped first-order obstruction triplet
through the four numbers

```text
(D01, D21, D41, N01).
```

### Compact BdG slope formulas

The explicit BdG support slopes reduce to the clean logarithmic forms

```text
B01 = 2 B0 (d(lambda_B)/lambda_B - d(varpi)/varpi)
B21 = 2 B2 (d(lambda_B)/lambda_B - 2 d(varpi)/varpi)
B41 = 2 B4 (d(lambda_B)/lambda_B - 3 d(varpi)/varpi).
```

So the support part of the primitive deformation is already exactly readable.

### Exact compensation surfaces

On the canonical compensated isotropic branch, the first two first-order gates are
algebraic.

1. **Even-preserving compensation** `K1 = 0` fixes the inertia-side slope exactly:


```text
K1 = D21 + D01/9 = 0
=> dM = D01/9 - B21 - Z21 = (dK - B01 - Z01)/9 - B21 - Z21.
```

2. **Odd/normalization-preserving compensation** `Xi_load = 0` fixes the static
   wall-loading slope exactly:

```text
Xi_load = N01/N0 - D01/D0 = 0
=> dK = B01 + Z01 + D0 (N01/N0).
```

So after Stage 209 the next live first-order gate is no longer vague. It is the
remaining hidden-even residual

```text
H_even = D41 - (2/3) D21 - D01/27.
```

That is exactly the explicit first-order 5PN gate left after the two compensation
surfaces are imposed.

---

## Stage 210 — mechanism sieve and the surviving outgoing corridor

Files:
- `5pn_stage210_mechanism_sieve_and_outgoing_corridor.py`
- `5pn_stage210_mechanism_sieve_and_outgoing_corridor_output.txt`

This stage turns the Stage-209 obstruction triplet into an actual sieve.

The three first-order gates are

```text
K1      = D21 + D01/9
Xi_load = N01/N0 - D01/D0
H_even  = D41 - (2/3) D21 - D01/27.
```

### Wall-only weak-axisymmetric anisotropy is dead

With only `(dK,dM)` active,

```text
D01 = dK
D21 = -dM
D41 = 0
N01 = 0
```

so

```text
K1_wall  = dK/9 - dM
Xi_wall  = -dK/D0
H_wall   = 2 dM/3 - dK/27.
```

The exact solve of

```text
(K1_wall, Xi_wall, H_wall) = 0
```

gives only the trivial branch

```text
dK = 0,
dM = 0.
```

So wall-only weak-axisymmetric anisotropy cannot carry the first-order 5PN closure.

### Pure BdG weak-axisymmetric anisotropy is also dead

Using logarithmic support drifts

```text
x_c     = delta ln C
x_varpi = delta ln varpi
```

the exact support slopes are

```text
B01 = 2 B0 (x_c - x_varpi)
B21 = 2 B2 (x_c - 2 x_varpi)
B41 = 2 B4 (x_c - 3 x_varpi)
```

and with wall/gauge/mixed data fixed,

```text
D01 = -B01
D21 = -B21
D41 = -B41
N01 = 0.
```

The exact solve of the full first-order system

```text
(K1, Xi_load, H_even) = 0
```

again gives only the trivial branch

```text
x_c = 0,
x_varpi = 0.
```

So pure BdG weak-axisymmetric anisotropy is also ruled out.

### BdG self-similarity kills only the load defect

On the pure-BdG branch, impose the natural self-similar condition

```text
x_c = x_varpi.
```

Then the script verifies

```text
Xi_load = 0,
K1      = 2 B2 x_varpi,
H_even  = (4/3) x_varpi (-B2 + 3 B4).
```

So BdG self-similarity kills only the normalization defect. It does **not** kill the
full first-order 5PN triplet unless the branch is trivial.

### Exact one-port outgoing-load factorization

The surviving nontrivial corridor is the outgoing mixed-sector load factor.
For one port,

```text
N0 = P^2 / Delta^2
```

with

```text
P     = Omega_U^2 G_W + R G_U
Delta = Omega_U^2 Omega_W^2 - R^2.
```

The script proves the exact factorization

```text
N0 / K = M_leg^2 (1 + I)^2 / (1 - H)^2
```

where

```text
M_leg = G_W / (Omega_W^2 sqrt(K))
I     = R G_U / (Omega_U^2 G_W)
H     = R^2 / (Omega_U^2 Omega_W^2).
```

### Exact one-port outgoing defect

The corresponding outgoing defect is

```text
Sigma^(N) = 2 m + 2 I/(1+I) i + 2 H/(1-H) h
```

with

- `m` the raw mixed-leg drift,
- `i` the interference-ratio drift,
- `h` the hybridization-ratio drift.

Under rigid interference and hybridization ratios,

```text
i = 0,
h = 0,
```

this collapses to

```text
Sigma^(N) = 2 m.
```

Writing the raw mixed-leg drift as

```text
m = g_W - 2 o_W - kappa_1/2,
```

the zero-defect condition becomes

```text
2 g_W - 4 o_W - kappa_1 = 0.
```

Equivalently,

```text
G_W / Omega_W^2 ∝ sqrt(K).
```

This is the exact **square-root mixed-leg law**. It is the first nontrivial
surviving first-order cancellation condition after the wall-only and BdG-only
primitive sectors are killed.

---

## Best current status after Stages 209–210

The weak-axisymmetric continuation is now much sharper.

The live first-order problem is no longer “find some anisotropy that works.”
What survives is a very short list:

1. the full primitive deformation only enters through
   ```text
   (D01, D21, D41, N01),
   ```
2. `K1 = 0` and `Xi_load = 0` are exact algebraic compensation surfaces,
3. wall-only and pure-BdG primitive sectors are dead,
4. BdG self-similarity kills only `Xi_load`, not the full 5PN triplet,
5. the surviving nontrivial corridor is outgoing mixed-sector co-loading,
6. and under rigid interference/hybridization the remaining first-order condition is
   exactly the square-root mixed-leg law
   ```text
   G_W / Omega_W^2 ∝ sqrt(K).
   ```

So the next honest theorem gate is now:

> choose a concrete mixed-sector weak-axisymmetric mechanism, substitute it into the
> Stage-209 compiler, and test whether it realizes the outgoing co-loading law while
> also killing the hidden-even residual.

That is the clean continuation point.


---

## Stage 211 — normalized monomial bridge and exact similarity kernel

Files:
- `5pn_stage211_normalized_similarity_kernel.py`
- `5pn_stage211_normalized_similarity_kernel_output.txt`

After Stage 210, the weak-axisymmetric continuation still needed the exact bridge
back to the later monomial/similarity-orbit package. Stage 211 now makes that bridge
runnable.

### Exact normalized monomial-drift matrix

In the normalized variable set

- `dln G_W`
- `dln G_U`
- `dln R`
- `dln K`
- `dln M`
- `dln Omega_U`
- `dln Omega_W`
- `dln delta_U`

and with constructive-branch parameters `(chi_0, delta_U, E_*, F_*)`, the direct
monomial-drift matrix is exactly

```text
M_norm =
[ -(1+delta_U),  1+delta_U,  1+delta_U,  0, 0, -2(1+delta_U),        0,      1+chi_0 ]
[            2,          0,      2E_*, -1, 1,        -2E_*, -(4+2E_*),        -F_* ]
[            0,          2,          0, -1, 1,           -2,         0,           0 ]
```

The script verifies this matrix exactly and shows

```text
rank(M_norm) = 3.
```

So the normalized zero-defect tangent space has dimension `5`.

### Exact triangular zero-defect solution

The zero-defect equations solve triangularly.

1. **Tracking** fixes `dln(delta_U)`:

```text
dln(delta_U)
= - (1+delta_U)/(1+chi_0)
  [ dln R + dln G_U - dln G_W - 2 dln Omega_U ].
```

2. **Dressing** fixes `dln(M)`:

```text
dln(M) = dln(K) - 2 dln(G_U) + 2 dln(Omega_U).
```

3. **Nontracking** fixes `dln(Omega_W)`:

```text
dln(Omega_W)
= [ dln(G_W) - dln(G_U) + (1-E_*) dln(Omega_U) + E_* dln(R)
    - (F_*/2) dln(delta_U) ] / (E_* + 2).
```

So once a candidate weak-axisymmetric branch supplies the five drifts

```text
(dln K, dln G_U, dln G_W, dln R, dln Omega_U),
```

the exact similarity-orbit theorem predicts the three co-drifts required to stay on
zero defect:

```text
(dln delta_U, dln M, dln Omega_W).
```

### Support-blind extension back to the explicit Stage-5 primitive space

Extending the primitive space back to

- `dln lambda_B`
- `dln varpi`
- plus the eight normalized drifts above,

the script proves that the first two columns of the extended monomial matrix are
exact zeros.

So the direct monomial package is exactly support-blind in the two explicit BdG
primitive directions:

```text
partial_{ln lambda_B} (dln C_tr, dln C_nt, dln epsilon_eta) = 0,
partial_{ln varpi}    (dln C_tr, dln C_nt, dln epsilon_eta) = 0.
```

The extended primitive monomial matrix still has rank `3`, so its nullity is

```text
10 - 3 = 7.
```

This splits into

1. the five normalized similarity directions, and
2. the two support-blind directions `(dln lambda_B, dln varpi)`.

### Why Stage 211 matters

This is the exact splice between the two halves of the 5PN continuation:

- Stages 209–210 now isolate which primitive first-order sectors are dead and which
  mixed/outgoing corridor survives,
- while Stage 211 shows exactly how the surviving mixed/wall/U normalization problem
  is constrained by monomial rigidity.

So the current weak-axisymmetric picture is now structurally clean:

- `Xi_load = P1/P0` lives on the same mixed/wall/U normalization problem controlled by
  the similarity-orbit theorem,
- but the explicit BdG-support directions remain outside that theorem because they are
  exact zero columns of the monomial map.

That is the sharpest current continuation point.
# 5PN continuation notes — stages 212 through 214

This session pushed the weak-axisymmetric program one step past the Stage-209–211 sieve.

The open question after Stage 211 was no longer whether the normalized similarity-orbit
package controls `Xi_load = P1/P0`. That part was already clean. The real next gate was the one
flagged in the notes after Stage 17:

> reinstate the conservative Maxwell/mixed `Z_2,Z_4` sector and see whether it removes the fake
> freedom left in the lower-bound even-gate picture.

These three stages do exactly that.

---

## Stage 212 — exact isotropic full-bundle target surface

Files:
- `5pn_stage212_isotropic_target_surface.py`
- `5pn_stage212_isotropic_target_surface_output.txt`

### What was frozen

Work on the isotropic grouped-lane branch with

```text
D0 = K - B0 - Z0
D2 = -(M + B2 + Z2)
D4 = -(B4 + Z4)
```

and compile the normalized grouped-response and outgoing-prefactor moments

```text
u2 = -D2 / D0
u4 = (D2^2 - D0 D4) / D0^2

P0 = N0 / D0
P2 = (D0 N2 - 2 D2 N0) / D0^2
P4 = (D0^2 N4 - 2 D0(D2 N2 + D4 N0) + 3 D2^2 N0) / D0^3.
```

### Exact one-pole defect

The script proves

```text
u4 - 4 u2^2 = [ D0 (B4 + Z4) - 3 (M + B2 + Z2)^2 ] / D0^2.
```

So the isotropic conservative one-pole surface is exactly

```text
D0 (B4 + Z4) = 3 (M + B2 + Z2)^2.
```

### Exact constant-prefactor branch

Imposing

```text
P2 = 0
P4 = 0
```

gives the exact outgoing conditions

```text
N2 = 2 D2 N0 / D0

N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2.
```

### Exact 2.5PN / 4PN normalization hit

The universal normalization condition is

```text
mhat0^2 P0 = 54 G c_s^5 / (5 a^5 c^5),
```

i.e.

```text
mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5).
```

### Most useful Stage-212 result

The completed isotropic moving-throat bundle must land on one exact algebraic target surface:

1. `D0 = K - B0 - Z0`,
2. `D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`,
3. `mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`,
4. `N2 = 2 D2 N0 / D0`,
5. `N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

So the actual PDE does not need to “show 5PN” in any vague sense. It needs to hit that surface.

---

## Stage 213 — exact `Z`-sector bridge back into the live even gates

Files:
- `5pn_stage213_z_sector_even_gate_bridge.py`
- `5pn_stage213_z_sector_even_gate_bridge_output.txt`

### Conservative Maxwell/mixed moments reinstated exactly

The omitted conservative block is

```text
Z0 = Q / Delta
Z2 = (Q S2 - H Delta) / Delta^2
Z4 = [ Q (S2^2 - Delta) - S2 H Delta ] / Delta^3
```

with

```text
Delta = Omega_U^2 Omega_W^2 - R^2
Q     = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2
H     = G_U^2 + G_W^2
S2    = Omega_U^2 + Omega_W^2.
```

The script differentiates these exactly and verifies

```text
dZ0 = (Delta dQ - Q dDelta) / Delta^2

dZ2 = [ Delta(-Delta dH - H dDelta + Q dS2 + S2 dQ)
        + 2 dDelta (Delta H - Q S2) ] / Delta^3

dZ4 = -[ Delta^2 H dS2 + Delta^2 S2 dH + Delta^2 dQ
          - 2 Delta H S2 dDelta - 2 Delta Q S2 dS2
          - 2 Delta Q dDelta - Delta S2^2 dQ + 3 Q S2^2 dDelta ] / Delta^4.
```

### Exact even-gate contribution

The missing `Z`-sector contribution to the two live even gates is

```text
K1_Z     = -dZ2 - dZ0/9
H_even,Z = -dZ4 + (2/3) dZ2 + dZ0/27.
```

### Constructive-slice evaluation

On the same constructive slice recorded in the notes,

```text
G_U = 5,  G_W = 7,  R = 2,
Omega_U = 3,  Omega_W = 4,
chi_0 = 3/2,  delta_U = 2/3,
E_* = 1/4,  F_* = 5/6,
```

the script gives

```text
K1_Z = (78623501/25004700) alpha_OmegaU + (733046/6251175) alpha_R
       - (59010631/25004700) alpha_U - (32134513/50009400) alpha_W

H_even,Z = -(28906377971/21003948000) alpha_OmegaU
           - (1174937411/21003948000) alpha_R
           + (11102468471/10501974000) alpha_U
           + (5617869293/21003948000) alpha_W.
```

In particular, the pure mixed directions are now active:

```text
alpha_W :  K1_Z = -32134513/50009400,
           H_even,Z = 5617869293/21003948000

alpha_R :  K1_Z = 733046/6251175,
           H_even,Z = -1174937411/21003948000.
```

### Most useful Stage-213 result

This is the first exact executable proof that the omitted `Z_2,Z_4` block does the precise job the
Stage-17 lower-bound picture predicted: it activates the previously untouched mixed directions
`alpha_W` and `alpha_R` in the even-gate problem.

---

## Stage 214 — full constructive-slice even-gate solve

Files:
- `5pn_stage214_full_even_gate_constructive_slice.py`
- `5pn_stage214_full_even_gate_constructive_slice_output.txt`

### Lower-bound matrix before `Z_2,Z_4`

Using only the matched wall block and the explicit one-mode BdG block with

```text
K = 2,
M = 3,
B0 = 2,
varpi = 3,
```

the lower-bound gate matrix in the seven directions

```text
(alpha_K, alpha_W, alpha_U, alpha_R, alpha_OmegaU, beta_B, beta_varpi)
```

has rank `2`, but its `alpha_W` and `alpha_R` columns vanish exactly.
So the lower-bound picture leaves those mixed directions completely untouched.

### Full matrix after exact `Z`-sector reinstatement

Adding the Stage-213 `Z`-sector contributions gives the full constructive-slice even-gate matrix

```text
[ -25/9,  -32134513/50009400,   91017569/25004700,    733046/6251175,
  -71404699/25004700,  -8/9,  4/3 ]

[  52/27,  5617869293/21003948000, -30905427529/10501974000,
  -1174937411/21003948000, 55109414029/21003948000, 32/81, -16/27 ]
```

with rank still `2`, so the full even-gate intersection remains five-dimensional.

### The structural fact that matters

The mixed-sector minor is now nonzero:

```text
det Gate_(alpha_W, alpha_R) = 942737330573 / 205838690400000 != 0.
```

So the apparently free mixed directions of the lower-bound picture were never genuine null
modes. They were only hidden in the omitted `Z_2,Z_4` block.

### Exact solve for the mixed directions

The full constructive-slice even-gate system solves directly for
`alpha_W` and `alpha_R`:

```text
alpha_W =
  (14503089433000/942737330573) alpha_K
+ (30450672110098/942737330573) alpha_OmegaU
- (29120459867142/942737330573) alpha_U
- (18876066395200/25453907925471) beta_B
+ (9438033197600/8484635975157) beta_varpi

alpha_R =
  (101802968743000/942737330573) alpha_K
+ (189815725996721/942737330573) alpha_OmegaU
- (188832473718440/942737330573) alpha_U
+ (89510801038400/25453907925471) beta_B
- (44755400519200/8484635975157) beta_varpi.
```

The script verifies by exact back-substitution that these kill both live even gates.

### Most useful Stage-214 result

There are no longer pure `alpha_W` or pure `alpha_R` null directions on the full constructive
branch. Once the conservative `Z_2,Z_4` sector is restored, the mixed-sector freedom is no longer
an open qualitative question. It is algebraically slaved to the remaining five directions.

---

## Best current status after stages 212–214

The next theorem gate is sharper again.

1. The exact isotropic full-bundle target surface is now explicit.
2. The omitted conservative `Z_2,Z_4` sector has been reinstated exactly.
3. On a clean constructive branch, that block removes the fake mixed-sector freedom left in the
   Stage-17 lower-bound picture.
4. So the remaining work is no longer “some mixed-sector freedom somewhere.” It is:

> compute the actual overlap data on the true moving-throat branch and test whether they land on
> the isotropic full-bundle target surface from Stage 212.

That is the cleanest next continuation point.
# 5PN continuation notes — stages 215 through 217

This session tightened the weak-axisymmetric 5PN program in a way that changes the next honest roadmap.

The key issue was that the recent mixed-sector continuation had still been phrased in the surrogate even-gate variables

- `K1 = D21 + D01/9`,
- `H_even = D41 - (2/3) D21 - D01/27`,
- `Xi_load = N01/N0 - D01/D0 = P1/P0`,

whereas the true reduced endgame from Stages 199–201 is formulated in the grouped Packet-A / Packet-B data.

So the first question for this session was:

> when does a weak-axisymmetric microscopic deformation actually remain tangent to the *true* first-order Packet-A/B endgame, rather than only solving a lower-bound surrogate gate?

That question is now answered exactly.

---

## Stage 215 — exact first-order Packet-A tangency theorem

Files:
- `5pn_stage215_first_order_packet_tangency_theorem.py`
- `5pn_stage215_first_order_packet_tangency_theorem_output.txt`

### Setup

Take an isotropic grouped-lane baseline with

```text
D_A0 = D0,
D_A2 = D2,
D_A4 = D4,
N_A0 = N0,
```

and add a weak-axisymmetric grouped perturbation with the exact signature

```text
lambda_(20) = 1,
lambda_(21) = 1/2,
lambda_(22) = -1.
```

So the perturbed grouped-lane data are

```text
D_A0 = D0 + epsilon lambda_A D01
D_A2 = D2 + epsilon lambda_A D21
D_A4 = D4 + epsilon lambda_A D41
N_A0 = N0 + epsilon lambda_A N01.
```

### Exact first-order grouped anomalies

Compiling the grouped response moments

```text
u2^(A) = -D_A2 / D_A0
u4^(A) = (D_A2^2 - D_A0 D_A4) / D_A0^2
P0^(A) = N_A0 / D_A0,
```

the script proves two exact facts.

First, the weighted grouped trace is invisible at first order:

```text
d/d epsilon ubar(u2)|0 = 0

d/d epsilon ubar(u4)|0 = 0

d/d epsilon ubar(P0)|0 = 0.
```

So on a genuine weak-axisymmetric grouped branch, the scalar trace defects
`Delta_pole` and `Delta_norm` are automatically invisible at `O(epsilon)`.

Second, the live first-order grouped anisotropy slots are one-scalar each:

```text
a(u2) =  epsilon u2^(1) / 4,
 b(u2) = 3 epsilon u2^(1) / 4,

a(u4) =  epsilon u4^(1) / 4,
 b(u4) = 3 epsilon u4^(1) / 4,

a(P0) =  epsilon P1 / 4,
 b(P0) = 3 epsilon P1 / 4.
```

So the entire first-order Packet-A tangency problem collapses to the three scalars

```text
u2^(1),  u4^(1),  P1.
```

### Exact operator formulas

The script derives the exact operator-slope formulas

```text
u2^(1) = (-D0 D21 + D2 D01) / D0^2
       = -(D21 + u2 D01) / D0
```

```text
u4^(1)
= [ D0(-D0 D41 + 2 D2 D21 - D4 D01) + 2 D01(D0 D4 - D2^2) ] / D0^3
= -(D41 + 2 u2 D21 + (u2^2 + u4) D01) / D0
```

and

```text
P1/P0 = N01/N0 - D01/D0 = Xi_load.
```

So the **true** first-order packet tangency conditions are exactly

```text
u2^(1) = 0,
nu4^(1) = 0,
P1 = 0.
```

### One-pole specialization

On an isotropic one-pole baseline `u4 = 4 u2^2`,

```text
u2^(1) = 0  ->  D21 = -u2 D01.
```

If this holds, then the second condition reduces to

```text
u4^(1) = 0  ->  D41 = -3 u2^2 D01.
```

So the exact one-pole tangent conditions are

```text
D21 = -u2 D01,
D41 = -3 u2^2 D01,
N01/N0 = D01/D0.
```

### Canonical compensated branch and the surrogate gates

On the **canonical compensated branch** used in the earlier notes,

```text
u2 = 1/9,
u4 = 4/81,
```

so the exact true tangency conditions become

```text
D21 = -D01/9,
D41 = -D01/27,
N01/N0 = D01/D0.
```

The script then proves the exact translation

```text
u2^(1) = -K1 / D0
```

and

```text
nu4^(1) = -(H_even + (8/9) K1) / D0.
```

Therefore, on the canonical compensated branch,

```text
u2^(1) = 0,
nu4^(1) = 0,
P1 = 0
```

is exactly equivalent to

```text
K1 = 0,
H_even = 0,
Xi_load = 0.
```

### Most useful Stage-215 result

This stage cleanly separates two things that had been easy to conflate:

1. **true first-order Packet-A tangency**, and
2. the earlier **surrogate even-gate variables**.

They coincide exactly on the canonical compensated branch, but not on a generic isotropic prototype with a different baseline `u2`.

That turns out to matter immediately in the next stage.

---

## Stage 216 — the support-blind mixed sector and the exact packet-null line

Files:
- `5pn_stage216_support_blind_packet_null_condition.py`
- `5pn_stage216_support_blind_packet_null_condition_output.txt`

### Explicit isotropic prototype branch

I then returned to the explicit positive isotropic finite-throat overlap prototype built from

```text
I_cross = 8 / (3 pi)
lambda_B = 1
varpi    = 2
G_U      = 1
G_W      = 1
R        = 1/2
Omega_U  = 3/2
Omega_W  = 2
M        = 1,
```

with `K` fixed on the exact isotropic one-pole surface

```text
D0 (B4 + Z4) = 3 (M + B2 + Z2)^2.
```

The first critical fact is that this explicit prototype does **not** sit on the canonical compensated branch. Its actual baseline value is

```text
u2 = (8575 + 12717 pi^2) / [210 (490 + 1503 pi^2)] ≈ 0.0416671714,
```

so

```text
u2 - 1/9 != 0.
```

Therefore the correct first-order packet test on this prototype is not the surrogate pair `(K1,H_even)`. It is the true Stage-215 packet system

```text
u2^(1) = 0,
nu4^(1) = 0,
Xi_load = 0.
```

### Support-blind mixed sector

I then imposed the same support-blind mixed-sector restriction at the input level:

- free drifts:
  ```text
  (alpha_K, alpha_GW, alpha_R),
  ```
- frozen support / upper-leg directions:
  ```text
  alpha_GU = 0,
  alpha_OU = 0,
  beta_B = 0,
  beta_varpi = 0,
  ```
- dependent co-drifts from the normalized orbit relations:
  ```text
  alpha_deltaU = -(1+deltaU)/(1+chi0) (alpha_R - alpha_GW),
  alpha_OW     = [alpha_GW + E_* alpha_R - (F_*/2) alpha_deltaU] / (E_* + 2).
  ```

So the support-blind mixed family is still parameterized by the three free directions

```text
(alpha_K, alpha_GW, alpha_R),
```

but now the **true** packet-null problem is the homogeneous `3 x 3` linear system generated by

```text
u2^(1),
nu4^(1),
Xi_load.
```

### Exact determinant and the packet-null line

The script computes the exact determinant of that `3 x 3` matrix and proves that it factors to a single affine line in the orbit constants `(E_*,F_*)`:

```text
det(M_packet) = 0
iff
F_* = F_crit(E_*),
```

with exact affine law

```text
F_crit(E_*) = -(A_E E_* + A_0) / A_F,
```

where

```text
A_E = 263797293760000
    + 1757766806455275 pi^2
    + 3339557838723645 pi^4
    + 1551622258297188 pi^6,

A_F = 48655861632000
    + 389171318788980 pi^2
    + 930178880126748 pi^4
    + 694451446430976 pi^6,

A_0 = 102703468791960 pi^6
    - 155911749769062 pi^4
    - 147002028439770 pi^2
    - 25886544768000.
```

So the support-blind mixed sector does **not** contain a generic true first-order Packet-A/B null direction.
It contains one only on the affine codimension-1 line

```text
F_* = F_crit(E_*).
```

### Concrete rank checks

The script then evaluates two concrete sample points.

At the representative constructive point already used in the notes,

```text
(E_*, F_*) = (1/4, 5/6),
```

the exact packet-null matrix has

```text
rank = 3,
```

so the support-blind packet-null system is **trivial** there.

At the illustrative point

```text
(E_*, F_*) = (1, 1),
```

the matrix again has

```text
rank = 3.
```

So the tempting support-blind mixed corridor that survived the lower-bound surrogate even gates does **not** survive the true first-order packet-null test on the explicit isotropic prototype at these positive orbit values.

### Most useful Stage-216 result

This is the first exact executable proof that the support-blind mixed corridor is much narrower than the lower-bound even-gate picture suggested.

It is not generically a first-order Packet-A/B null direction.
It exists only on one affine line in `(E_*,F_*)`.

---

## Stage 217 — positive-quadrant obstruction and the updated roadmap

Files:
- `5pn_stage217_positive_quadrant_obstruction.py`
- `5pn_stage217_positive_quadrant_obstruction_output.txt`

### Exact sign theorem for the packet-null line

Stage 217 turns the affine line from Stage 216 into a real obstruction theorem.

The exact line is

```text
F_* = F_crit(E_*),
```

and the script proves that both its slope and intercept are negative:

```text
dF_crit/dE_* < 0,
F_crit(0) < 0.
```

Numerically,

```text
F_crit(E_*) ≈ -0.1076895540 - 2.4072204236 E_*.
```

In particular,

```text
F_crit(0)   ≈ -0.1076895540,
F_crit(1/4) ≈ -0.7094946599,
F_crit(1/2) ≈ -1.3112997658,
F_crit(1)   ≈ -2.5149099776.
```

So for every `E_* >= 0`,

```text
F_crit(E_*) < 0.
```

### Consequence

The support-blind packet-null line never enters the constructive positive quadrant

```text
E_* >= 0,
F_* > 0.
```

So the support-blind mixed sector does admit a **mathematical** first-order packet-null corridor,
but only outside the positive constructive-orbit branch.

### Explicit mathematical null vector outside the constructive quadrant

To show the obstruction is physical rather than algebraic, the script evaluates one point on the line, namely

```text
E_* = 1/4,
F_* = F_crit(1/4) ≈ -0.7094946599.
```

There the packet-null matrix has

```text
rank = 2,
```

so a genuine one-parameter null corridor exists. A normalized null vector is

```text
(alpha_K, alpha_GW, alpha_R)
≈ (0.2698838731, 0.6492432127, 1).
```

So the support-blind mixed corridor is not algebraically impossible.
It is simply pushed out of the positive constructive-orbit branch.

### Most useful Stage-217 result

This is the new roadmap verdict:

1. the earlier support-blind mixed corridor was only a lower-bound even-gate survivor,
2. once the **true** first-order Packet-A/B tangency conditions are imposed on the explicit isotropic prototype,
   that corridor survives only on a negative-`F_*` affine line,
3. so it does **not** survive on the positive constructive-orbit branch.

That means the next honest 5PN continuation is now much sharper:

> restore the support-carrying directions
> `(alpha_GU, alpha_OU, beta_B, beta_varpi)`
> and recompute the full first-order packet-null matrix,
> or else move to a different isotropic baseline branch before re-testing the endgame.

This is a real change in the continuation point.
The support-blind mixed sector is no longer the default favorite corridor once the actual first-order Packet-A/B endgame is enforced.
# 5PN continuation notes — stages 218 through 220

These stages answer the fork opened by the Stage-217 support-blind obstruction.

The question was:

> does the failure of the support-blind mixed corridor force an isotropic-baseline pivot,
> or does it only mean that the next honest move is to restore the support carriers and
> solve the true first-order Packet-A/B master matrix on the same branch?

The answer is now sharp.

It does **not** force a pivot yet.

The support-blind corridor was killed, but the isotropic baseline itself survives.
Once the support-carrying directions are restored, the true first-order Packet-A/B master
matrix immediately regains nontrivial nullspace on the same explicit isotropic prototype.

So the correct next roadmap is:

1. restore the support variables,
2. solve the full first-order packet master matrix,
3. then ask which of those algebraic null directions can actually be realized by the
   moving-throat branch.

A baseline pivot is now the **backup** path, not the default one.

---

## Stage 218 — support-restored master matrix theorem

Files:
- `5pn_stage218_support_restored_master_matrix.py`
- `5pn_stage218_support_restored_master_matrix_output.txt`

### 1. Full true first-order packet system

Keep the same explicit positive isotropic finite-throat prototype used in Stages 215–217.
The true first-order packet-tangency system is still

```text
u2^(1) = 0,
u4^(1) = 0,
Xi_load = 0.
```

Restore the full support-carrying direction set

```text
(alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi),
```

with the orbit-lock-dependent drifts imposed exactly:

```text
alpha_deltaU = -(1+deltaU)/(1+chi0) (alpha_R + alpha_GU - alpha_GW - 2 alpha_OU),
```

```text
alpha_OW = [alpha_GW - alpha_GU + (1-E_*) alpha_OU + E_* alpha_R - (F_*/2) alpha_deltaU] / (E_* + 2),
```

```text
alpha_M = alpha_K - 2 alpha_GU + 2 alpha_OU.
```

So the full first-order Packet-A/B master matrix is the `3 x 7` Jacobian of

```text
(u2^(1), u4^(1), Xi_load)
```

with respect to those seven free directions.

### 2. Support-only 3x3 minor is already nonzero in the positive constructive quadrant

The key new theorem is that the support-only `3 x 3` minor in

```text
(alpha_K, alpha_GU, alpha_OU)
```

has exact determinant

```text
det M_support
= - pi^2 (8575 + 12717 pi^2)^2 [A_E E_* + A_F F_* + A_0]
  / [3961650000 (490 + 1503 pi^2)^6 (E_* + 2)],
```

with

```text
A_E = 36399941765600000
    + 252110109946857750 pi^2
    + 578056501787708475 pi^4
    + 433070868615970095 pi^6,
```

```text
A_F = 9400136276160000
    + 48552477588959400 pi^2
    + 110238264175381020 pi^4
    + 85279016350877148 pi^6,
```

```text
A_0 = 22826428394560000
    + 216604585787832900 pi^2
    + 597245590420119270 pi^4
    + 566382238265309238 pi^6.
```

All three coefficients are strictly positive.
So for every

```text
E_* >= 0,
F_* >= 0,
```

and in particular on the positive constructive branch `F_* > 0`,

```text
det M_support != 0.
```

### 3. Rank/nullity verdict

Because that support-only minor is already nonzero, the exact verdict is:

```text
rank(M_support) = 3,
nullity(M_support) = 2,
```

and therefore the full support-restored master matrix also has

```text
rank(M_master) = 3,
nullity(M_master) = 4.
```

So the Stage-217 obstruction has a very specific meaning:

- it kills the **support-blind** corridor,
- but it does **not** kill the current isotropic baseline.

That is the main decision result of this session.

---

## Stage 219 — support-only constructive rescue

Files:
- `5pn_stage219_support_only_constructive_rescue.py`
- `5pn_stage219_support_only_constructive_rescue_output.txt`

To make the Stage-218 theorem concrete, I then froze the representative constructive point

```text
(E_*, F_*) = (1/4, 5/6)
```

and restricted to the support-only slice

```text
alpha_GW = 0,
alpha_R  = 0.
```

So the remaining free support carriers are

```text
(beta_B, beta_varpi).
```

The exact solve of

```text
u2^(1) = 0,
u4^(1) = 0,
Xi_load = 0
```

for

```text
(alpha_K, alpha_GU, alpha_OU)
```

gives a two-parameter packet-null family:

```text
alpha_K  =  0.0898376063372746 beta_B - 0.183439324158911 beta_varpi,
alpha_GU =  0.0736940106940628 beta_B - 0.154862196153272 beta_varpi,
alpha_OU =  0.0385380345483458 beta_B - 0.105749511712718 beta_varpi.
```

So the constructive branch already contains two independent exact support-only null directions.
A convenient basis is:

### `beta_B` basis

```text
(alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi)
≈
( 0.0898376063,
  0,
  0.0736940107,
  0,
  0.0385380345,
  1,
  0).
```

### `beta_varpi` basis

```text
(alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi)
≈
(-0.1834393242,
  0,
 -0.1548621962,
  0,
 -0.1057495117,
  0,
  1).
```

Both were checked directly in the script and satisfy

```text
u2^(1) = 0,
u4^(1) = 0,
Xi_load = 0
```

exactly.

So on the same explicit isotropic prototype that killed the support-blind corridor,
restoring just the BdG support carriers is already enough to recover a true first-order
Packet-A/B null family.

---

## Stage 220 — full constructive support-restored master solve

Files:
- `5pn_stage220_full_constructive_master_solver.py`
- `5pn_stage220_full_constructive_master_solver_output.txt`

Finally, I solved the **full** support-restored constructive master system at

```text
(E_*, F_*) = (1/4, 5/6)
```

without freezing the mixed carriers.

A convenient free-carrier choice is

```text
(alpha_GW, alpha_R, beta_B, beta_varpi),
```

with the wall/support drifts

```text
(alpha_K, alpha_GU, alpha_OU)
```

fixed linearly by the true packet-null conditions.

Numerically, the exact constructive solve is

```text
alpha_K
=  0.729738128983757 alpha_GW
 - 0.779773212960421 alpha_R
 + 0.0898376063372746 beta_B
 - 0.183439324158911 beta_varpi,
```

```text
alpha_GU
= -0.269396497721430 alpha_GW
 + 0.367875865936610 alpha_R
 + 0.0736940106940628 beta_B
 - 0.154862196153272 beta_varpi,
```

```text
alpha_OU
= -0.216658637559893 alpha_GW
 + 0.268624898608163 alpha_R
 + 0.0385380345483458 beta_B
 - 0.105749511712718 beta_varpi.
```

So the full constructive branch carries a `4`-parameter exact packet-null family.
A canonical numerical basis is:

### `alpha_GW` basis

```text
( 0.729738129,
  1,
 -0.269396498,
  0,
 -0.216658638,
  0,
  0).
```

### `alpha_R` basis

```text
(-0.779773213,
  0,
  0.367875866,
  1,
  0.268624899,
  0,
  0).
```

### `beta_B` basis

```text
( 0.0898376063,
  0,
  0.0736940107,
  0,
  0.0385380345,
  1,
  0).
```

### `beta_varpi` basis

```text
(-0.1834393242,
  0,
 -0.1548621962,
  0,
 -0.1057495117,
  0,
  1).
```

The Stage-219 support-only rescue is recovered exactly by setting

```text
alpha_GW = alpha_R = 0.
```

So the full master solve confirms the roadmap verdict from Stage 218.

---

## What changed in the roadmap

Before these stages, the honest continuation fork was:

1. restore the support carriers and solve the master matrix, or
2. pivot the isotropic baseline.

After these stages, the updated verdict is:

### Not yet time to pivot

The current isotropic baseline remains live.
The support-blind obstruction was real, but it was only a slice obstruction.

### The correct next move is now explicit

The next theorem gate is:

> determine which of the support-restored first-order packet-null directions can actually
> be realized by the moving-throat branch.

That means the next work should be about **realizability**, not about reopening the baseline choice.

Concretely, the next continuation should now attack one of two closely related tasks:

1. derive the actual branch-induced weak-axisymmetric support drifts
   `(beta_B, beta_varpi)` and mixed drifts `(alpha_GW, alpha_R)` from the moving-throat overlap data,
   then project them onto the Stage-220 null basis;
2. or derive the branch tangent map from the moving-throat reduced PDE into
   `(alpha_K, alpha_GU, alpha_GW, alpha_R, alpha_OU, beta_B, beta_varpi)`
   and test whether its image intersects the exact packet-null family.

That is a much sharper continuation point than where Stage 217 left us.

---

## Honest caveat

This is a real rescue of the current isotropic prototype, but it is still only a
**first-order reduced packet-null result**.

What is still not proved is that the actual moving-throat branch realizes any of these null directions.
So the serious remaining problem is now:

> algebraic packet-null existence is established,
> but branch realization is still open.

That is the right place for the next push.
# 5PN continuation notes — stages 221 through 223

## Purpose

These stages take the support-restored packet-null family from Stages 218–220 and attack the next honest problem:

> how does an actual moving-throat overlap branch map into that packet-null family, and does the simplest coherent profile branch land inside it?

This is the realizability step that had remained open after the support-restoration rescue. The answer is now much sharper.

There are three distinct results:

1. an exact **realizability compiler** from overlap-state tangent variables into the Stage-220 packet-null carrier space,
2. a **no-go theorem** for the pure coherent-profile branch,
3. and a unique **compensated realization family** showing how the branch can still survive once the upper-leg and frequency variables are allowed to co-move.

---

## Stage 221 — exact overlap realizability map

Files:
- `5pn_stage221_overlap_realizability_map.py`
- `5pn_stage221_overlap_realizability_map_output.txt`

### Input

Start from the Stage-220 support-restored packet-null carrier variables

- `alpha_K`
- `alpha_GW`
- `alpha_GU`
- `alpha_R`
- `alpha_OU`
- `beta_B`
- `beta_varpi`

and restrict them to the moving-throat overlap dictionary supplied by the coherent finite-throat branch:

- `C = lambda_B * kappa`
- `G_U = lambda_U`
- `G_W = lambda_W * kappa`
- `R = lambda_R * kappa`

with primitive overlap-state drifts

- `sigma_K      = d ln K_geo`
- `sigma_kappa  = d ln kappa`
- `ell_B        = d ln lambda_B`
- `ell_W        = d ln lambda_W`
- `ell_R        = d ln lambda_R`

and compensation variables

- `ell_U        = d ln lambda_U`
- `ell_OmegaU   = d ln Omega_U`
- `ell_varpi    = d ln varpi`.

### Exact carrier map

The overlap-state tangent maps into the Stage-220 carrier space as

- `alpha_K     = sigma_K`
- `alpha_GW    = ell_W + sigma_kappa`
- `alpha_GU    = ell_U`
- `alpha_R     = ell_R + sigma_kappa`
- `alpha_OU    = ell_OmegaU`
- `beta_B      = ell_B + sigma_kappa`
- `beta_varpi  = ell_varpi`.

So the realizability question becomes purely linear:

> does there exist a choice of `(ell_U, ell_OmegaU, ell_varpi)` such that the induced carrier vector lies in the exact Stage-220 packet-null family?

### Exact realization solve

At the constructive point `(E_*,F_*) = (1/4, 5/6)`, solving the full packet-null equations gives

\[
ell_U
\approx
0.84421482069524\,\sigma_K
+0.13857149963136\,\sigma_\kappa
-0.002148228031649\,\ell_B
+1.02617196909895\,\ell_R
-0.885452241435934\,\ell_W,
\]

\[
ell_{\Omega_U}
\approx
0.57648223573430\,\sigma_K
+0.067558848518307\,\sigma_\kappa
-0.013251749605984\,\ell_B
+0.718150303781305\,\ell_R
-0.637339705657013\,\ell_W,
\]

\[
ell_{\varpi}
\approx
-5.45139383054918\,\sigma_K
+0.216979224836926\,\sigma_\kappa
+0.489740172938324\,\ell_B
-4.25085088235995\,\ell_R
+3.97808993425855\,\ell_W.
\]

So the overlap-state tangent is not generically free. But it does admit an exact three-parameter compensation solve.

### Induced orbit/companion drifts

The same solve fixes the dependent drifts required to stay on the Stage-220 packet-null manifold:

\[
d\ln\delta_U
\approx
0.50522670126549\,\sigma_K
-0.0056516769732235\,\sigma_\kappa
-0.039854080113249\,\ell_B
-0.965244046150368\,\ell_R
+0.999446449290394\,\ell_W,
\]

\[
d\ln M
\approx
0.46453483007811\,\sigma_K
-0.14202530222611\,\sigma_\kappa
-0.022207043148670\,\ell_B
-0.616043330635281\,\ell_R
+0.496225071557841\,\ell_W,
\]

\[
d\ln\Omega_W
\approx
-0.27660634196525\,\sigma_K
+0.51753444540572\,\sigma_\kappa
+0.0039179033515622\,\ell_B
+0.0731670124294903\,\ell_R
+0.440449529624671\,\ell_W.
\]

### Meaning

This is the first exact realizability compiler from the moving-throat overlap variables into the support-restored 5PN packet-null family.

So the old vague question

> “does the branch somehow realize the packet-null directions?”

has now become

> “does the actual branch generate the exact compensation ratios above?”

---

## Stage 222 — pure coherent-profile branch no-go

Files:
- `5pn_stage222_pure_profile_branch_nogo.py`
- `5pn_stage222_pure_profile_branch_nogo_output.txt`

### Setup

The cleanest candidate realization is the pure coherent-profile branch coming from the nonconstant finite-throat family already isolated in the moving-throat notes:

- only `kappa(theta)` and `K_geo(theta)` move,
- all couplings `lambda_B, lambda_W, lambda_R, lambda_U` are frozen,
- all frequencies are frozen.

In carrier language this is the tangent

\[
(\alpha_K,\alpha_{GW},\alpha_{GU},\alpha_R,\alpha_{OU},\beta_B,\beta_{\varpi})
=
(\sigma_K,\sigma_\kappa,0,\sigma_\kappa,0,\sigma_\kappa,0).
\]

### Exact result

Solving the full packet-null equations on this restricted branch gives only the trivial solution

\[
\sigma_K = \sigma_\kappa = 0.
\]

So the pure profile branch is an exact no-go for nontrivial packet-null realizability.

### Actual coherent family inserted

The nonconstant coherent family has

\[
\kappa(\theta)=\kappa_0\cos\theta+\kappa_1\sin\theta,
\qquad
\kappa_0=\frac{2\sqrt2}{\pi},
\qquad
\kappa_1=-\frac{4}{3\pi},
\]

and

\[
K_{\rm geo}(\theta)=K_0+T_{\rm grad}\sin^2\theta.
\]

Therefore

\[
\sigma_\kappa(\theta)=\partial_\theta\ln\kappa(\theta)\,d\theta
=
 d\theta\,\frac{3\sqrt2\sin\theta+2\cos\theta}{2\sin\theta-3\sqrt2\cos\theta},
\]

\[
\sigma_K(\theta)=\partial_\theta\ln K_{\rm geo}(\theta)\,d\theta
=
\frac{2T_{\rm grad}\sin(2\theta)}{2K_0-T_{\rm grad}\cos(2\theta)+T_{\rm grad}}\,d\theta.
\]

Two concrete points are especially useful:

- At `theta = 0`:
  \[
  \sigma_\kappa/d\theta = -\frac{\sqrt2}{3},
  \qquad
  \sigma_K/d\theta = 0.
  \]

- At the max-coupling point `theta = theta_max = -arctan(\sqrt2/3)`:
  \[
  \sigma_\kappa/d\theta = 0,
  \qquad
  \sigma_K/d\theta = -\frac{12\sqrt2\,T_{\rm grad}}{22K_0+4T_{\rm grad}}.
  \]

So the actual nonconstant overlap family does produce nontrivial profile drifts, but those drifts do **not** lie inside the packet-null family by themselves.

### Meaning

This is a real obstruction, but not a dead end.

It says only that profile motion alone cannot carry the reduced 5PN closure.
The upper-leg and frequency sector must co-evolve with it.

---

## Stage 223 — unique compensated profile family

Files:
- `5pn_stage223_profile_compensation_family.py`
- `5pn_stage223_profile_compensation_family_output.txt`

### Setup

Now keep the same coherent profile branch, but allow the upper-leg and frequency variables to compensate:

- `d ln lambda_U`
- `d ln Omega_U`
- `d ln varpi`

while still freezing `lambda_B`, `lambda_W`, and `lambda_R`.

This is the smallest compensated realization family worth testing after the pure-profile no-go.

### Exact compensation law

The Stage-221 realization compiler then collapses to

\[
d\ln\lambda_U
\approx
0.84421482069524\,\sigma_K
+0.13857149963136\,\sigma_\kappa,
\]

\[
d\ln\Omega_U
\approx
0.57648223573430\,\sigma_K
+0.067558848518307\,\sigma_\kappa,
\]

\[
d\ln\varpi
\approx
-5.45139383054918\,\sigma_K
+0.216979224836926\,\sigma_\kappa.
\]

So once the coherent profile chooses a tangent `(sigma_K, sigma_kappa)`, there is a **unique** compensation direction in `(lambda_U, Omega_U, varpi)` that restores packet-nullness.

### Induced companion drifts

The required co-moving branch also forces

\[
d\ln\delta_U
\approx
0.50522670126549\,\sigma_K
-0.0056516769732235\,\sigma_\kappa,
\]

\[
d\ln M
\approx
0.46453483007811\,\sigma_K
-0.14202530222611\,\sigma_\kappa,
\]

\[
d\ln\Omega_W
\approx
-0.27660634196525\,\sigma_K
+0.51753444540572\,\sigma_\kappa.
\]

So the realization family is not just “adjust three knobs.” It is a full co-evolving microscopic branch.

### Concrete evaluations

#### At `theta = 0`

Since `sigma_K = 0` and `sigma_kappa = -(sqrt(2)/3) dtheta`, the compensation slopes are

\[
\frac{d\ln\lambda_U}{d\theta}
\approx -0.0653232313790,
\qquad
\frac{d\ln\Omega_U}{d\theta}
\approx -0.0318475466110,
\qquad
\frac{d\ln\varpi}{d\theta}
\approx -0.102284987506.
\]

The induced companion drifts are

\[
\frac{d\ln\delta_U}{d\theta}
\approx 0.00266422607523,
\qquad
\frac{d\ln M}{d\theta}
\approx 0.0669513695361,
\qquad
\frac{d\ln\Omega_W}{d\theta}
\approx -0.243968077229.
\]

#### At `theta = theta_max`

Since `sigma_kappa = 0`, the whole family is driven only by `sigma_K`, so all compensation slopes become proportional to

\[
\frac{T_{\rm grad}}{2K_0 + 0.363636\ldots T_{\rm grad}}.
\]

The main ones are

\[
\frac{d\ln\lambda_U}{d\theta}
\approx
-1.30243641707\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln\Omega_U}{d\theta}
\approx
-0.889384359537\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln\varpi}{d\theta}
\approx
+8.41029282436\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}}.
\]

The induced drifts are

\[
\frac{d\ln\delta_U}{d\theta}
\approx
-0.779452857821\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln M}{d\theta}
\approx
-0.716674316608\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}},
\]

\[
\frac{d\ln\Omega_W}{d\theta}
\approx
+0.42674229845\,
\frac{T_{\rm grad}}{2K_0+0.363636\ldots T_{\rm grad}}.
\]

### Meaning

This is the strongest positive result of the session.

The pure profile branch is dead, but the actual coherent overlap family is **not**.
It survives as a unique compensated realization family inside the exact Stage-220 packet-null manifold.

So the remaining question is no longer algebraic and no longer generic. It is now extremely concrete:

> does the true moving-throat PDE branch dynamically generate the compensation ratios required by the Stage-223 family?

---

## Best current summary after Stages 221–223

The realizability problem is now sharply split.

### What is ruled out

- A pure coherent-profile deformation with frozen upper-leg and frequency data cannot realize the 5PN packet-null condition.

### What survives

- The coherent finite-throat overlap family does admit a unique support/frequency compensation direction that restores exact packet-nullness.

### What the PDE must now decide

The completed moving-throat branch has to choose between two possibilities:

1. it does **not** co-evolve in the Stage-223 ratios, in which case the present reduced 5PN closure fails on that branch;
2. it **does** co-evolve in those ratios, in which case the reduced 5PN obstruction survives only at the deeper level of actual branch existence and stability.

So the program is now exactly where it should be:

- the algebraic bottleneck is gone,
- the realization bottleneck is explicit,
- and the next test is a genuine microscopic branch-selection test rather than another symbolic repackaging.
# 5PN continuation notes — stages 224 through 227

These stages continue directly from the Stage-221/223 overlap realizability work.
The target was to answer the next honest question after the compensated coherent-profile family existed algebraically:

> do the **actual moving-throat coherent branch relations** force that compensation family,
> kill it, or leave it underdetermined?

The answer is now sharp:

- the current moving-throat coherent branch laws do **not** kill the compensated packet-null family,
- they do **not** force it either,
- they leave it **underdetermined by exactly four scalar conditions**.

So this is **not** a dead end. It is a precise realizability bottleneck.

---

## Stage 224 — exact overlap-image parameterization

Files:
- `5pn_stage224_overlap_image_parameterization.py`
- `5pn_stage224_overlap_image_parameterization_output.txt`

This stage takes the constructive support-restored packet-null family and rewrites it completely in the overlap-observable drift variables

\[
(d\ln K,\ d\ln M,\ d\ln C,\ d\ln\varpi,\ d\ln\Omega_U,\ d\ln\Omega_W,
\ d\ln\chi_0,\ d\ln\epsilon_\eta,\ d\ln Z_W,\ d\ln\epsilon_W,\ d\ln\delta_U).
\]

### Main result

Although the compensated family begins from five primitive overlap-side drifts,

\[
(\sigma_K,\ \sigma_\kappa,\ \ell_B,\ \ell_W,\ \ell_R),
\]

its observable overlap image has exact rank `4`.

A convenient exact coordinate system on that image is

\[
(d\ln\chi_0,\ d\ln Z_W,\ d\ln C,\ d\ln\varpi).
\]

So the packet-null family already splits naturally into

- a **2D orbit/shape layer**
  \[
  (d\ln\chi_0,\ d\ln Z_W),
  \]
- a **2D support layer**
  \[
  (d\ln C,\ d\ln\varpi).
  \]

The exact constructive image also gives the simple universal formulas

\[
d\ln\delta_U = -\frac{18}{11} d\ln\chi_0,
\]
\[
d\ln\Omega_W = \frac{41}{44} d\ln\chi_0 + \frac58 d\ln Z_W,
\]
\[
d\ln\epsilon_W = 2 d\ln\chi_0 + d\ln Z_W,
\]
\[
d\ln\epsilon_\eta = 0.
\]

The nontrivial remaining image coordinates

\[
d\ln K,
\qquad
d\ln M,
\qquad
d\ln\Omega_U,
\]

are exact linear functions of
\((d\ln\chi_0,d\ln Z_W,d\ln C,d\ln\varpi)\).
The script keeps them symbolically exact and prints useful numerical approximations.

---

## Stage 225 — exact similarity-orbit bridge in overlap variables

Files:
- `5pn_stage225_similarity_orbit_bridge.py`
- `5pn_stage225_similarity_orbit_bridge_output.txt`

This stage takes the direct coherent moving-throat branch laws from the compact PDE program and rewrites them in overlap variables.
The three exact coherent zero-defect residuals are

\[
\Sigma_\eta = d\ln\epsilon_\eta,
\]
\[
\Sigma_{\rm tr} = (1+\chi_0)\,d\ln\delta_U + (1+\delta_U)\,d\ln\chi_0,
\]
\[
\Sigma_{\rm nt} = (d\ln Z_W - 2 d\ln\Omega_W) + E_* d\ln\epsilon_W - F_* d\ln\delta_U.
\]

### Main result

On the constructive overlap image,

\[
\Sigma_\eta = \Sigma_{\rm tr} = \Sigma_{\rm nt} = 0
\]

identically.

So the full compensated overlap image already lies tangent to the exact monomial-preserving similarity orbit of the coherent moving-throat branch.

Just as important, the orbit verdict is **support-blind**:

\[
\partial_{d\ln C}\Sigma_\eta
=
\partial_{d\ln C}\Sigma_{\rm tr}
=
\partial_{d\ln C}\Sigma_{\rm nt}
=0,
\]
\[
\partial_{d\ln\varpi}\Sigma_\eta
=
\partial_{d\ln\varpi}\Sigma_{\rm tr}
=
\partial_{d\ln\varpi}\Sigma_{\rm nt}
=0.
\]

So the current moving-throat coherent/invariant laws only see the orbit/shape side, not the explicit support pair.

---

## Stage 226 — branch-selection theorem

Files:
- `5pn_stage226_branch_selection_theorem.py`
- `5pn_stage226_branch_selection_theorem_output.txt`

This is the actual answer to the fork question.

The full constructive overlap image is encoded by seven exact residuals

\[
(R_K,\ R_M,\ R_{\Omega_U},\ R_{\Omega_W},\ R_{\epsilon_\eta},\ R_{\epsilon_W},\ R_{\delta_U}).
\]

The actual moving-throat coherent branch laws are only the three monomial-preservation residuals

\[
(B_\eta,\ B_{\rm tr},\ B_{\rm nt})
=
(\Sigma_\eta,\Sigma_{\rm tr},\Sigma_{\rm nt}).
\]

### Exact equivalence proved by the script

The full overlap-image residual system is exactly equivalent to

1. the three actual branch-law residuals
   \[
   (B_\eta,B_{\rm tr},B_{\rm nt})=0,
   \]
2. plus one exact **selector quartet**
   \[
   (R_K,\ R_M,\ R_{\Omega_U},\ R_{\epsilon_W})=0.
   \]

Once those four selectors vanish, the remaining two image residuals vanish automatically:

\[
R_{\delta_U}=0,
\qquad
R_{\Omega_W}=0.
\]

### Dimension count

The script computes exact ranks:

- full overlap-state space: dimension `11`,
- actual coherent branch-law system: rank `3`, so branch-law manifold dimension `8`,
- constructive overlap image: rank `7`, so image dimension `4`.

Therefore the codimension gap is

\[
7-3 = 4.
\]

### Final verdict

The current moving-throat coherent branch laws

- do **not** forbid the compensation family,
- do **not** force it,
- and leave it **underdetermined by exactly four scalar conditions**.

That is the cleanest status statement we have had so far.

---

## Stage 227 — overlap realizability tester

Files:
- `5pn_stage227_overlap_realizability_tester.py`
- `5pn_stage227_overlap_realizability_tester_output.txt`

This stage packages the full overlap-image condition as a practical finite-dimensional tester.

### Input
An observed first-order overlap tangent

\[
(d\ln K,\ d\ln M,\ d\ln C,\ d\ln\varpi,\ d\ln\Omega_U,\ d\ln\Omega_W,
\ d\ln\chi_0,\ d\ln\epsilon_\eta,\ d\ln Z_W,\ d\ln\epsilon_W,\ d\ln\delta_U).
\]

### Output
The exact seven residuals

\[
(R_K,\ R_M,\ R_{\Omega_U},\ R_{\Omega_W},\ R_{\epsilon_\eta},\ R_{\epsilon_W},\ R_{\delta_U}).
\]

and the quadratic realizability distance

\[
D_{\rm real}^2 = \sum_i R_i^2.
\]

So once we have an actual moving-throat overlap extraction from a more expensive PDE or branch computation, there is now a direct yes/no test against the current constructive packet-null family.

---

## What changed in the roadmap

Before Stage 224–227, the continuation question was still phrased loosely as

> “Does the actual moving-throat branch generate the compensation ratios?”

Now it is much sharper.

The actual moving-throat coherent branch already guarantees the three monomial-preserving relations.
That is **not** the remaining bottleneck.

The remaining bottleneck is the exact selector quartet

\[
(R_K,\ R_M,\ R_{\Omega_U},\ R_{\epsilon_W}).
\]

So the next honest theorem gate is:

> extract the actual overlap tangent of the moving-throat branch and test whether those four remaining selector residuals vanish.

That is smaller and cleaner than the old “derive 5PN somehow” problem.
# 5PN continuation notes — stages 228 through 230

These stages continue directly from the Stage-224/226 overlap-image branch-selection result.

The previous verdict was:

- the actual coherent moving-throat branch laws do **not** kill the compensated overlap image,
- they do **not** force it either,
- and at the raw overlap-state level they leave the image underdetermined by a four-equation selector quartet.

The point of this continuation was to sharpen that further by working **inside the coherent branch-law manifold itself** and solving out the latent support pair.

The outcome is a real reduction of the theorem gate:

> once the actual branch is already coherent-law clean, overlap-image realizability is no longer a four-selector problem. It is an exact **two-scalar** problem.

---

## Stage 228 — exact branch-law reduction of the overlap-image test

Files:
- `5pn_stage228_branch_law_two_scalar_reduction.py`

Start from the full Stage-226 overlap-image residual vector
\[
(R_K,
 R_M,
 R_{\Omega_U},
 R_{\Omega_W},
 R_{\epsilon_\eta},
 R_{\epsilon_W},
 R_{\delta_U}).
\]

Restrict to the actual coherent moving-throat branch-law manifold
\[
B_\eta=0,
\qquad
B_{\rm tr}=0,
\qquad
B_{\rm nt}=0.
\]

On the constructive branch with
\[
\chi_0=\frac29,
\qquad
\delta_U=1,
\qquad
E_*=\frac14,
\qquad
F_*=\frac56,
\]
these branch laws imply
\[
\delta\ln\epsilon_\eta=0,
\qquad
\delta\ln\delta_U=-\frac{18}{11}\,\delta\ln\chi_0,
\]

and
\[
\delta\ln\Omega_W
=
\frac{1}{2}
\Bigl(
\delta\ln Z_W+
E_*\,\delta\ln\epsilon_W-
F_*\,\delta\ln\delta_U
\Bigr).
\]

After substitution, the seven overlap-image residuals collapse exactly to:

- one support residual triple
  \[
  (R_K,
   R_M,
   R_{\Omega_U}),
  \]
- one orbit/shape selector
  \[
  S_{\rm shape}
  :=
  \delta\ln\epsilon_W-
  \bigl(2\,\delta\ln\chi_0+\delta\ln Z_W\bigr),
  \]
- and the exact identities
  \[
  R_{\epsilon_\eta}=0,
  \qquad
  R_{\delta_U}=0,
  \qquad
  R_{\Omega_W}=\frac18 S_{\rm shape},
  \qquad
  R_{\epsilon_W}=S_{\rm shape}.
  \]

So on the coherent branch-law manifold, overlap-image membership is already equivalent to
\[
R_K=R_M=R_{\Omega_U}=0,
\qquad
S_{\rm shape}=0.
\]

That is already sharper than the Stage-226 selector quartet.

---

## Stage 229 — exact support-plane compiler

Files:
- `5pn_stage229_support_plane_compiler.py`
- `5pn_stage229_support_plane_compiler_output.txt`

The remaining latent support coordinates are
\[
(\delta\ln C,
\delta\ln\varpi).
\]

Inside the coherent branch-law manifold, the support residual pair
\[
R_K=0,
\qquad
R_M=0
\]
solves **uniquely** for that latent support pair in terms of the observed support-side pair
\[
(\delta\ln K,
\delta\ln M)
\]
and the orbit coordinates
\[
(\delta\ln\chi_0,
\delta\ln Z_W).
\]

Numerically, on the constructive branch,
\[
\delta\ln C
\approx
-22.0706526877\,\delta\ln K
+287.773569763\,\delta\ln M
+49.3529290037\,\delta\ln Z_W
+364.327814331\,\delta\ln\chi_0,
\]
\[
\delta\ln\varpi
\approx
-22.5706526877\,\delta\ln K
+191.654233426\,\delta\ln M
+32.3559287390\,\delta\ln Z_W
+234.769642713\,\delta\ln\chi_0.
\]

After that solve, the whole support-side realizability question collapses to one exact codimension-1 plane:
\[
S_{\rm support}
:=
\delta\ln\Omega_U-
\delta\ln\Omega_U^{({\rm pred})}=0,
\]
with
\[
\delta\ln\Omega_U^{({\rm pred})}
\approx
1.95241347905\,\delta\ln K
-12.3996352800\,\delta\ln M
-1.99787596774\,\delta\ln Z_W
-14.3144730533\,\delta\ln\chi_0.
\]

So the support residual triple is no longer a three-condition problem once the latent support pair is allowed to adjust. It is one exact scalar plane.

---

## Stage 230 — reduced two-scalar realizability tester and restoration map

Files:
- `5pn_stage230_reduced_branch_realizability_tester.py`
- `5pn_stage230_reduced_branch_realizability_tester_output.txt`

After imposing the coherent branch laws and solving out the latent support pair from
\[
(\delta\ln K,\delta\ln M),
\]
the full seven-component residual vector collapses exactly to
\[
(0,
 0,
 S_{\rm support},
 S_{\rm shape}/8,
 0,
 S_{\rm shape},
 0).
\]

So the entire overlap-image membership problem on the coherent branch-law manifold is now:
\[
\boxed{
S_{\rm shape}=0,
\qquad
S_{\rm support}=0.
}
\]

Equivalently, the actual moving-throat tangent no longer needs to be tested in the full 11-component overlap state.
It only needs the reduced observable sextuple
\[
(\delta\ln K,
\delta\ln M,
\delta\ln\Omega_U,
\delta\ln\chi_0,
\delta\ln Z_W,
\delta\ln\epsilon_W),
\]
and then the two exact selectors above.

### Exact restoration map

The minimal restoration back to the compensated overlap image is now explicit:

1. fix the latent support pair to the Stage-229 target values
   \[
   (\delta\ln C,\delta\ln\varpi),
   \]
2. shift
   \[
   \delta\ln\Omega_U \to \delta\ln\Omega_U-S_{\rm support},
   \]
3. shift the coherent shape pair by
   \[
   \delta\ln\epsilon_W \to \delta\ln\epsilon_W-S_{\rm shape},
   \qquad
   \delta\ln\Omega_W \to \delta\ln\Omega_W-\frac18 S_{\rm shape},
   \]
which preserves the coherent branch-law manifold and kills the remaining overlap residuals.

So there is now an exact reduced restoration operator on the coherent branch-law manifold.

---

## What this changes in the 5PN roadmap

The next theorem gate is now much smaller than before.

Before these stages, the remaining question after Stage 226 was roughly:

> does the actual moving-throat branch land in the full compensated overlap image?

After Stages 228–230, that has sharpened to:

> on the coherent branch-law manifold, does the actual moving-throat branch satisfy the two exact selectors
> \[
> S_{\rm shape}=0,
> \qquad
> S_{\rm support}=0?
> \]

So the live reduced task is now:

1. extract the reduced observable sextuple
   \[
   (\delta\ln K,
   \delta\ln M,
   \delta\ln\Omega_U,
   \delta\ln\chi_0,
   \delta\ln Z_W,
   \delta\ln\epsilon_W)
   \]
   from the actual moving-throat branch,
2. evaluate the two selectors above,
3. and check whether both vanish.

That is a real tightening of the branch-selection bottleneck.
# 5PN continuation notes — stages 231 through 233

These stages continue directly from the Stage-228/230 two-selector reduction, but now tie it back to the **actual microscopic coherent moving-throat variables** rather than leaving it as an abstract overlap-image test.

The payoff is strong:

> on the minimal coherent continuum branch, the Stage-230 overlap-image selectors are no longer abstract. They collapse to **two explicit microscopic defect scalars**.

Those two scalars are:

- the **dressing slippage**
  \[
  \Sigma_\eta = 2c_1-\kappa_U-\kappa_\eta,
  \]
- the **support-plane scalar**
  \[
  \Sigma_{\rm sup}
  =
  \frac{\kappa_U}{2}
  -c_K\kappa_\eta
  +c_\chi\Sigma_\chi
  +c_Z\zeta_Z,
  \]
  with
  \[
  c_K=\frac{8575}{4392}\approx 1.95241347905,
  \]
  \[
  c_\chi\approx 14.3144730533,
  \qquad
  c_Z\approx 1.99787596774,
  \]
  and
  \[
  \Sigma_\chi = c_1+\gamma_1-\kappa_U,
  \qquad
  \zeta_Z = 2\lambda_1-\kappa_\eta-\kappa_W.
  \]

So the next theorem gate is now much sharper than “extract the whole 11-component tangent.”
It is:

> does the actual moving-throat branch dynamically generate
> \[
> \Sigma_\eta = 0,
> \qquad
> \Sigma_{\rm sup}=0
> \]
> on the coherent branch-law manifold?

---

## Stage 231 — exact coherent selector bridge

Files:
- `5pn_stage231_coherent_selector_bridge.py`
- `5pn_stage231_coherent_selector_bridge_output.txt`

Using the Stage-230 selector formulas and substituting the actual minimal coherent-continuum observables
\[
 d\ln K = \kappa_\eta,
 \qquad
 d\ln M = 0,
 \qquad
 d\ln \Omega_U = \kappa_U/2,
\]
\[
 d\ln \chi_0 = \Sigma_\chi = c_1+\gamma_1-\kappa_U,
\]
\[
 d\ln Z_W = \zeta_Z = 2\lambda_1-\kappa_\eta-\kappa_W,
\]
\[
 d\ln \epsilon_W = \Sigma_\epsilon = 2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
\]
one gets the exact collapse
\[
 S_{\rm shape} = -\Sigma_\eta.
\]

This is already a major simplification. The first reduced selector is not some hidden overlap residual. It is exactly the microscopic dressing slippage.

Even better, Stage-166 had already shown
\[
 \mathcal R_1+\Xi_1
 =
 -\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
\]
Combining the two gives the exact selector bridge
\[
 \mathcal R_1+\Xi_1
 =
 \frac{\epsilon_\eta}{1-\epsilon_\eta}
 S_{\rm shape}.
\]
So on the physical branch \(0<\epsilon_\eta<1\), the following are equivalent:
\[
 S_{\rm shape}=0,
 \qquad
 \Sigma_\eta=0,
 \qquad
 \mathcal R_1+\Xi_1=0.
\]

This means the Stage-230 shape selector is the overlap-image form of the same microscopic dressing mismatch that already controlled the selected-branch demand law.

---

## Stage 232 — exact support-plane selector on the coherent continuum branch

Files:
- `5pn_stage232_coherent_support_plane_test.py`
- `5pn_stage232_coherent_support_plane_test_output.txt`

The second Stage-230 selector is
\[
 S_{\rm support}
 =
 d\ln \Omega_U - d\ln \Omega_U^{({\rm pred})}(d\ln K,d\ln M,d\ln\chi_0,d\ln Z_W).
\]

On the same minimal coherent-continuum branch, this becomes one exact microscopic scalar:
\[
 S_{\rm support} = \Sigma_{\rm sup}
 =
 \frac{\kappa_U}{2}
 -c_K\kappa_\eta
 +c_\chi\Sigma_\chi
 +c_Z\zeta_Z.
\]

So the support selector depends only on

- wall-stiffness drift \(\kappa_\eta\),
- internal-\(U\) frequency drift \(\kappa_U/2\),
- coherent interference-ratio drift \(\Sigma_\chi\),
- wall-to-mixed overlap drift \(\zeta_Z\).

And it is **independent of** \(d\ln\epsilon_W\) at this reduced level.

That is the second real narrowing of the theorem gap: the support part of the overlap-image test is now a concrete microscopic plane condition instead of an abstract overlap residual.

---

## Stage 233 — microscopic two-selector theorem and selected-spectral corollary

Files:
- `5pn_stage233_microscopic_two_selector_theorem.py`
- `5pn_stage233_microscopic_two_selector_theorem_output.txt`

Combining Stages 231 and 232 with the exact Stage-230 residual vector gives the full reduced microscopic form
\[
 (0,0,S_{\rm support},S_{\rm shape}/8,0,S_{\rm shape},0)
 =
 (0,0,\Sigma_{\rm sup},-\Sigma_\eta/8,0,-\Sigma_\eta,0).
\]

So overlap-image membership on the coherent branch-law manifold is now equivalent to exactly
\[
 \boxed{\Sigma_\eta=0,\qquad \Sigma_{\rm sup}=0.}
\]

This is the sharpest current 5PN microscopic continuation point.

### Exact restoration map

Because the Stage-230 restoration map is
\[
 d\ln\Omega_U \to d\ln\Omega_U-S_{\rm support},
\]
\[
 d\ln\epsilon_W \to d\ln\epsilon_W-S_{\rm shape},
\]
\[
 d\ln\Omega_W \to d\ln\Omega_W-S_{\rm shape}/8,
\]
it becomes, microscopically,
\[
 d\ln\Omega_U \to d\ln\Omega_U-\Sigma_{\rm sup},
\]
\[
 d\ln\epsilon_W \to d\ln\epsilon_W+\Sigma_\eta,
\]
\[
 d\ln\Omega_W \to d\ln\Omega_W+\Sigma_\eta/8.
\]
So the reduced overlap-image restoration has a direct microscopic meaning.

### Raw selected-spectral branch corollary

Applying the same selector formulas to the raw selected-spectral branch from the earlier spectral extractor,
\[
 d\ln K = \sigma_K^{\rm raw},
 \qquad
 d\ln M = 0,
 \qquad
 d\ln\Omega_U = 0,
\]
\[
 d\ln\chi_0 = 0,
 \qquad
 d\ln Z_W = 2\sigma_\kappa^{\rm raw}-\sigma_K^{\rm raw},
 \qquad
 d\ln\epsilon_W = 2\sigma_\kappa^{\rm raw},
\]
one gets the exact selector miss
\[
 S_{\rm shape}^{({\rm raw})} = \sigma_K^{\rm raw},
\]
\[
 S_{\rm support}^{({\rm raw})}
 =
 2c_Z\sigma_\kappa^{\rm raw}-(c_K+c_Z)\sigma_K^{\rm raw}.
\]
Numerically,
\[
 S_{\rm support}^{({\rm raw})}
 \approx
 3.99575193547\,\sigma_\kappa^{\rm raw}
 -3.95028944679\,\sigma_K^{\rm raw}.
\]

So the raw selected branch misses the packet-null image in exactly two ways:

1. a **shape/dressing miss** carried by \(\sigma_K^{\rm raw}\),
2. a **support-plane miss** carried by a fixed linear combination of \(\sigma_K^{\rm raw}\) and \(\sigma_\kappa^{\rm raw}\).

That is the cleanest current diagnostic of what actual moving-throat co-evolution the 5PN branch still has to generate.

---

## Best current continuation point after stages 231–233

The continuation point is now very narrow.

It is no longer:
- extract the whole overlap tangent,
- or solve a generic packet-null image problem,
- or guess how the selected branch should be corrected.

It is now exactly:

> compute the physical coherent-branch slippages that control
> \(\Sigma_\eta\) and \(\Sigma_{\rm sup}\),
> and test whether the actual moving-throat branch drives both to zero.

That is the next honest theorem gate.
# 5PN continuation notes — stages 234 through 236

These stages push the coherent tracking branch into the first explicit support-compensation and harmonic-selection theorem.

## Stage 234 — exact support-compensation theorem on the coherent tracking branch

Using the tracking-branch functions

\[
G_{\rm tr}(\xi,\delta;R)=\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]

\[
F_{\rm tr}(\xi,\delta;R)=
\frac{\bigl[9\delta+(9+2R^2)\xi\bigr]^2\bigl[9\delta+(9+2R)\xi\bigr]^2}
{81(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2},
\]

and the coherent support-enhancement factor

\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]

the exact support threshold is

\[
\zeta_{\rm req}=\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)},
\qquad
S_{\rm req}=\frac{M_{\rm req}}{M_{\rm mix}},
\qquad
M_{\rm req}=G_{\rm tr}(\xi_{\rm req},\delta;R).
\]

The derivative

\[
\frac{dS}{d\zeta}=\frac{1-\epsilon}{(1-\zeta\epsilon)^2}>0
\]

shows strict monotonicity, so the support threshold is unique once the stable-side branch point \(\xi_{\rm req}\) is fixed.

Because

\[
\frac1\epsilon-\zeta_{\rm req}=\frac{1-\epsilon}{\epsilon\,[1+\epsilon(S_{\rm req}-2)]}>0,
\]

there is no reduced-level support no-go for any finite \(S_{\rm req}>1\) on \(0<\epsilon<1\).

## Stage 235 — explicit D/N overlap extraction of the coherent support ratio

For the exact finite-throat D/N family

\[
\chi_n(s)=\sqrt{\frac{2}{L}}\sin\!\frac{(n+1/2)\pi s}{L},
\]

with the first uniform local source density \(\sigma(s)=1\), the exact overlap law is

\[
I_n=\int_0^L \chi_n(s)\,ds
=\frac{2\sqrt{2L}}{\pi(2n+1)},
\qquad
\frac{I_n}{I_0}=\frac{1}{2n+1}.
\]

So the physical coherent support ratio becomes

\[
\zeta_n^{(\rm phys)}=
\frac{K_W^{(\rm eff)}}{K_{\phi,n}^{(\rm eff)}}\,\frac{1}{(2n+1)^2}.
\]

On the same-operator twin family,

\[
K_{\phi,n}^{(\rm eff)}=K_W^{(\rm eff)}\bigl(1+x n(n+1)\bigr),
\]

this collapses to

\[
\zeta_n^{(\rm twin)}=
\frac{1}{(2n+1)^2\bigl(1+x n(n+1)\bigr)}.
\]

The exact lowest-twin value is therefore

\[
\zeta_0^{(\rm twin)}=1.
\]

## Stage 236 — exact lowest-twin sufficiency and higher-harmonic selection rules

Because \(\zeta_0^{(\rm twin)}=1\), the lowest symmetric twin lane has exact enhancement

\[
S_0=S(1;\epsilon)=2,
\]

independent of \(\epsilon\).

Therefore the lowest symmetric twin lane succeeds exactly when

\[
\zeta_{\rm req}\le 1,
\qquad\text{equivalently}\qquad
S_{\rm req}\le 2.
\]

For higher D/N harmonics,

\[
\zeta_n^{(\rm twin)}\le \frac{1}{(2n+1)^2},
\]

so harmonic \(n\) is ruled out immediately if

\[
\zeta_{\rm req}>\frac{1}{(2n+1)^2}.
\]

When that immediate impossibility bound is not violated, the exact twin softness threshold is

\[
x\le x_{\max}(n;\zeta_{\rm req})
:=
\frac{1/\bigl((2n+1)^2\zeta_{\rm req}\bigr)-1}{n(n+1)}.
\]

The exact enhancement ceiling at fixed \(\epsilon\) is

\[
S_n^{(\max)}
=
1+\frac{1-\epsilon}{(2n+1)^2-\epsilon}.
\]

So the explicit D/N support tower is strongly biased toward the lowest symmetric twin lane.

## Net result after stage 236

The coherent support question is no longer vague.

1. The tracking-branch support threshold is one exact scalar \(\zeta_{\rm req}\).
2. The explicit D/N overlap tower makes that threshold microscopic.
3. The lowest symmetric twin lane is an exact universal doubling branch with \(S_0=2\).
4. Every higher support harmonic is rapidly ruled out by exact overlap and stiffness suppression.

So the next clean theorem gate is now extremely narrow:

> Does the actual moving-throat PDE place the physical coherent support sector on the lowest twin D/N branch with \(\zeta_{\rm req}\le 1\), or is stronger-than-twin asymmetry required?
# 5PN continuation notes — stages 237 through 241

These stages fold the Stage-236 twin-selection theorem into the explicit moving-throat support/source operator from the compact PDE program.

The main structural change is that the coherent support problem is no longer phrased in abstract overlap variables.
It is now organized in three nested layers:

1. the **twin theorem** from Stage 236,
2. the first explicit **physical non-twin lowest-lane family** in terms of transport/compliance/stiffness variables,
3. the **operator-selected** and then **parent-projected** microscopic gain thresholds.

So after Stage 241 the next theorem gate is not “derive more support formulas somehow.”
It is:

> compute the actual parent confinement-loading amplitude and the equilibrium source/support alignment on the real moving-throat branch, then compare the resulting microscopic gain directly against the exact fail/succeed thresholds.

## Stage 237 — physical lowest-lane placement map

Starting from the first explicit non-twin lowest support family from the moving-throat notes, the physical support ratio is now

\[
\zeta_0^{(Pe+R)}(Pe,y;\kappa)
=
\Omega_{Pe}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2},
\]

with

\[
\Omega_{Pe}
=
\frac{\pi Pe\bigl(2Pe\,e^{Pe}+\pi\bigr)}{(4Pe^2+\pi^2)(e^{Pe}-1)}.
\]

Here:

- `Pe` is the physical axial source-drift Peclet number,
- `y` is the Robin eigenvalue parameter for the lowest support lane,
- `kappa = K_X L^2/T_X` is the baseline support stiffness ratio.

Exact identities proved in the script:

\[
\Omega_{Pe}(0)=1,
\qquad
\lim_{Pe\to +\infty}\Omega_{Pe}=\frac{\pi}{2}.
\]

So the Stage-236 symmetric twin baseline is recovered exactly at

\[
Pe=0,
\qquad
y=\frac{\pi}{2},
\qquad
\zeta_0^{(Pe+R)}=1.
\]

The exact closure ceiling of the first explicit non-twin family is

\[
\zeta_{\max}(\kappa)
=
\frac{\pi^2}{4}\,\frac{\kappa+\pi^2/4}{\kappa}.
\]

Therefore the Stage-236 twin theorem expands to the first physical phase split:

\[
\zeta_{\rm req}\le 1
\quad\Rightarrow\quad
\text{symmetric lowest twin enough},
\]

\[
1<\zeta_{\rm req}\le \zeta_{\max}(\kappa)
\quad\Rightarrow\quad
\text{first explicit non-twin family can in principle rescue},
\]

\[
\zeta_{\rm req}>\zeta_{\max}(\kappa)
\quad\Rightarrow\quad
\text{first explicit non-twin family fails identically}.
\]

Equivalently, when `\zeta_req > \pi^2/4`, the exact stiffness ceiling is

\[
\kappa \le \kappa_{\max}(\zeta_{\rm req})
:=
\frac{\pi^4}{4(4\zeta_{\rm req}-\pi^2)}.
\]

## Stage 238 — operator-selected support branch and exact `Xi` thresholds

The Stage-237 placement map is then evaluated on the first explicit coupled support/source branch.
The transport bias is no longer free; it is the root of the exact fixed-point equation

\[
Pe_* = \Xi\,\Delta(Pe_*;\kappa,\eta),
\]

with support-drop kernel endpoints

\[
\Delta_0(\kappa,\eta)
=
\frac{\eta(\cosh\sqrt{\kappa}-1)}{\kappa\bigl(\eta\cosh\sqrt{\kappa}+\sqrt{\kappa}\sinh\sqrt{\kappa}\bigr)},
\]

\[
\Delta_{\infty}(\kappa,\eta)
=
\frac{\eta\sinh\sqrt{\kappa}+\sqrt{\kappa}(\cosh\sqrt{\kappa}-1)}{\sqrt{\kappa}\bigl(\eta\cosh\sqrt{\kappa}+\sqrt{\kappa}\sinh\sqrt{\kappa}\bigr)}.
\]

Every constructive branch point obeys the exact bracket

\[
\Xi\Delta_0\le Pe_*\le \Xi\Delta_{\infty}.
\]

This immediately gives the exact support thresholds

\[
\Xi_{\rm fail} = \frac{Pe_{\rm req}}{\Delta_{\infty}(\kappa,\eta)},
\qquad
\Xi_{\rm suff} = \frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\]

where `Pe_req` is the unique constructive transport bias that hits the demanded support ratio.
So the operator-selected branch has the exact three-zone structure:

- `Xi <= Xi_fail` : impossible in this operator family,
- `Xi >= Xi_suff` : guaranteed success,
- `Xi_fail < Xi < Xi_suff` : only then is the full root solve still needed.

Useful exact endpoint data verified in the script:

\[
\lim_{\kappa\to 0^+}\Delta_0 = \frac12,
\qquad
\lim_{\kappa\to 0^+}\Delta_{\infty}=1,
\]

and in the highly compliant mouth limit,

\[
\lim_{\eta\to +\infty}\Delta_0 = \frac{1-\operatorname{sech}(\sqrt{\kappa})}{\kappa},
\qquad
\lim_{\eta\to +\infty}\Delta_{\infty} = \frac{\tanh(\sqrt{\kappa})}{\sqrt{\kappa}}.
\]

## Stage 239 — parent-action microscopic gain

The phenomenological operator strength is then pushed back to the parent 4D action.
Projecting the `n=5` GNLS/compressional sector onto one source channel and one support channel gives the exact microscopic gain

\[
G_{\rm micro}
=
\frac{\rho_* g_\phi^2 O_{\sigma\phi}^2}{m c_{s*}^2 K_X N_{\sigma\sigma}}.
\]

Equivalently, introducing the source/support coherence factor

\[
C_{\sigma\phi}^2 := \frac{O_{\sigma\phi}^2}{N_{\sigma\sigma}N_{\phi\phi}},
\]

one gets the exact factorization

\[
G_{\rm micro}
=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s*}^2 K_X}
\,C_{\sigma\phi}^2,
\qquad
0\le C_{\sigma\phi}^2\le 1.
\]

So the best-case gain at fixed confinement loading is

\[
G_{\max}(g_\phi)
=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s*}^2 K_X},
\]

reached only at perfect source/support alignment.

Using `\kappa = K_X L^2/T_X`, the exact fixed-point strength becomes

\[
\Xi_{\rm micro}=\kappa G_{\rm micro}
=
\frac{\rho_* g_\phi^2 O_{\sigma\phi}^2 L^2}{m c_{s*}^2 T_X N_{\sigma\sigma}}.
\]

So the support/source theorem gap is no longer an abstract gain problem.
It is a parent-overlap problem.

## Stage 240 — exact parent threshold decision theorem

Combining Stages 237–239 gives the exact microscopic phase diagram

\[
G_{\rm fail}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_{\infty}(\kappa,\eta)},
\qquad
G_{\rm suff}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_0(\kappa,\eta)}.
\]

So inside the non-twin reachability window:

- `G_micro <= G_fail` : fail,
- `G_micro >= G_suff` : success,
- only the band `G_fail < G_micro < G_suff` requires the full fixed-point solve.

In parent variables, this becomes exact threshold surfaces on the confinement-loading amplitude:

\[
g_{\phi,\rm fail}^2
=
\frac{m c_{s*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2\Delta_{\infty}(\kappa,\eta)},
\]

\[
g_{\phi,\rm suff}^2
=
\frac{m c_{s*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2\Delta_0(\kappa,\eta)},
\]

and on the coherence factor:

\[
C_{\rm fail}^2
=
\frac{m c_{s*}^2 K_X G_{\rm fail}}{\rho_* g_\phi^2 N_{\phi\phi}},
\qquad
C_{\rm suff}^2
=
\frac{m c_{s*}^2 K_X G_{\rm suff}}{\rho_* g_\phi^2 N_{\phi\phi}}.
\]

This also yields the first exact parent-overlap no-go theorem:
if

\[
G_{\max}(g_\phi) < G_{\rm fail}(\kappa,\eta),
\]

then the branch fails for **every possible** source profile in the chosen support channel.

So the full support decision now has a strict order:

1. check the Stage-236 twin theorem,
2. if needed, check the Stage-237 non-twin reachability ceiling,
3. only then compare the actual parent branch against the exact microscopic gain thresholds.

## Stage 241 — parent equilibrium alignment removes the free coherence datum

The parent equilibrium equations sharpen things one more step.
On the local static branch,

\[
H(y)\,\delta\rho(s,y)+\delta V_{\rm conf}(s,y)=0,
\qquad
H(y):=h'(\rho_*(y)),
\]

so the source profile induced by a support displacement is not free:

\[
\chi_\sigma(y)=\frac{g_\phi\chi_\phi(y)}{H(y)}.
\]

This forces the overlap invariants to become

\[
O_{\sigma\phi}=g_\phi I_1,
\qquad
N_{\sigma\sigma}=g_\phi^2 I_2,
\]

with

\[
I_1=\int d^3y\,\frac{\chi_\phi(y)^2}{H(y)},
\qquad
I_2=\int d^3y\,\frac{\chi_\phi(y)^2}{H(y)^2}.
\]

Therefore the coherence factor is now derived, not free:

\[
C_{\sigma\phi}^2=
\frac{I_1^2}{N_{\phi\phi} I_2}\le 1.
\]

The exact equilibrium gain becomes

\[
G_{\rm eq} = \frac{g_\phi^2 I_1}{K_X}.
\]

In the thin matched layer where `H(y)` is nearly constant,

\[
C_{\sigma\phi}^2=1,
\qquad
G_{\rm eq} = \frac{g_\phi^2 N_{\phi\phi}}{K_X H_w}
= \frac{\rho_w g_\phi^2 N_{\phi\phi}}{m c_{s,w}^2 K_X}.
\]

So the best-alignment formulas used above are not ad hoc. They are the natural thin-layer limit of the parent equilibrium branch.

## Net result after Stage 241

The support theorem is now much narrower than it was at Stage 236.

1. The symmetric lowest twin lane still succeeds exactly iff `zeta_req <= 1`.
2. If `zeta_req > 1`, the first explicit non-twin lowest-lane family is controlled by the physical placement map `zeta_0^(Pe+R)(Pe,y;kappa)`.
3. The transport bias `Pe` is not free; it is selected by the exact support/source fixed-point equation `Pe = Xi Delta(Pe;kappa,eta)`.
4. The microscopic gain is no longer phenomenological; it is a parent-action overlap quantity.
5. On the equilibrium branch, the source/support coherence is no longer independent either.

So the remaining theorem gap is now almost completely localized:

> compute the actual parent confinement-loading amplitude `g_phi`, the active-layer compressional stiffness profile `H(y)`, and the support profile `chi_phi(y)` on the real moving-throat branch, form the resulting equilibrium gain `G_eq`, and compare it directly to the exact fail/succeed thresholds inside the first explicit non-twin support family.
# 5PN continuation notes — stages 242 through 246

These stages pick up exactly from the Stage-241 equilibrium-alignment result and carry the support/source theorem through the first explicit parent confinement branch.

The structural change is that the support theorem is no longer phrased in terms of the abstract gain
`G_eq = g_phi^2 I_1 / K_X`
with free parent data `(g_phi, I_1, K_X)`.
It is now organized in five increasingly concrete layers:

1. the first explicit thin-wall confinement family,
2. the one-number wall figure of merit,
3. the explicit GNLS wall-shell reduction of `T_X`, `K_X`, and `kappa`,
4. the canonical `tanh` wall / local mouth closure branch,
5. the final explicit branch-placement problem in the three parent dimensionless controls
   `(chi_s, Lambda_ell, Upsilon_w)`.

So after Stage 246 the next theorem gate is no longer “compute `g_phi` somehow.”
It is:

> compute the actual moving-throat branch placement in `(chi_s, Lambda_ell, Upsilon_w)` and compare it directly to the exact support-threshold surfaces.

## Stage 242 — explicit thin-wall confinement branch

Take the first concrete parent wall family
\[
V_{\rm conf}(r;a) = V_0\, f\!\left(\frac{r-a}{\ell}\right),
\qquad
\xi := \frac{r-a}{\ell}.
\]

A support displacement `a -> a + phi(s)` gives
\[
\delta V_{\rm conf}
=
-\partial_a V_{\rm conf}\,\phi(s)
=
+\frac{V_0}{\ell} f'(\xi)\,\phi(s),
\]
so the parent support-loading amplitude is fixed exactly:
\[
g_\phi = \frac{V_0}{\ell}.
\]

Using the shell measure
\[
d^3y = 4\pi r^2 dr = 4\pi \ell (a+\ell \xi)^2 d\xi,
\]
the exact shell integral entering the equilibrium gain becomes
\[
I_1
=
4\pi \ell
\left[
a^2 J_1 + 2 a \ell J_2 + \ell^2 J_3
\right],
\]
where
\[
J_1 = \int d\xi \frac{f'(\xi)^2}{H(\xi)},
\qquad
J_2 = \int d\xi \frac{\xi f'(\xi)^2}{H(\xi)},
\qquad
J_3 = \int d\xi \frac{\xi^2 f'(\xi)^2}{H(\xi)}.
\]

For a centered symmetric wall layer, `J_2 = 0`, so
\[
I_1 = 4\pi \ell\left[a^2 J_1 + \ell^2 J_3\right].
\]

Therefore the exact equilibrium-aligned gain is
\[
G_{\rm eq}
=
\frac{g_\phi^2 I_1}{K_X}
=
\frac{4\pi V_0^2}{K_X}
\left[
\frac{a^2 J_1}{\ell}
+ 2 a J_2
+ \ell J_3
\right].
\]

On the centered branch, the thin-wall leading term is
\[
G_{\rm eq}^{(\rm tw)}
=
\frac{4\pi a^2 V_0^2 J_1}{K_X \ell}.
\]

So the parent support/source theorem is now a direct wall-amplitude theorem on the first explicit wall family.

## Stage 243 — wall figure of merit and direct fail/succeed thresholds

Introduce the support geometry parameter
\[
\kappa = \frac{K_X L^2}{T_X}.
\]
Then the natural dimensionless wall control variable is
\[
W_{\rm wall}
:=
\kappa G_{\rm eq}^{(\rm tw)}
=
\frac{4\pi a^2 L^2 J_1 V_0^2}{T_X \ell}.
\]

Using the operator-selected threshold pair
\[
G_{\rm fail} = \frac{Pe_{\rm req}}{\kappa \Delta_\infty(\kappa,\eta)},
\qquad
G_{\rm suff} = \frac{Pe_{\rm req}}{\kappa \Delta_0(\kappa,\eta)},
\]
the whole branch collapses to the exact support theorem
\[
W_{\rm wall} \le \frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)}
\quad\Rightarrow\quad
\text{fail},
\]
\[
W_{\rm wall} \ge \frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)}
\quad\Rightarrow\quad
\text{succeed}.
\]

The physical wall-amplitude thresholds are therefore
\[
V_{0,\rm fail}^2
=
\frac{T_X \ell\, Pe_{\rm req}}{4\pi a^2 L^2 J_1 \Delta_\infty},
\qquad
V_{0,\rm suff}^2
=
\frac{T_X \ell\, Pe_{\rm req}}{4\pi a^2 L^2 J_1 \Delta_0}.
\]

So the support/source gate is no longer an abstract gain comparison. It is one dimensionless wall figure-of-merit comparison.

## Stage 244 — explicit GNLS wall-shell reduction and exact `W_wall = Xi`

Now derive the support coefficients directly from the parent GNLS energy on a thin active shell around the throat wall.

With the shell support mode
\[
\delta\rho(s,y)=q(s)\chi_\phi(y),
\]
the quadratic reduced support energy is
\[
E^{(2)}[q]
=
\frac12 \int_0^L ds \left[T_X q'(s)^2 + K_X q(s)^2\right],
\]
where
\[
T_X = \frac{\hbar^2}{4m\rho_w} N_{\phi\phi},
\qquad
K_X = H_w N_{\phi\phi} + \frac{\hbar^2}{4m\rho_w} G_{\phi\phi},
\]
and
\[
H_w = \frac{m c_{s,w}^2}{\rho_w}.
\]

For the explicit thin shell with `chi_phi = f'(\xi)`,
\[
N_{\phi\phi} = 4\pi a^2 \ell I_f,
\qquad
G_{\phi\phi} = \frac{4\pi a^2 I_g}{\ell},
\]
so
\[
T_X = \frac{\pi a^2 \ell I_f \hbar^2}{m\rho_w},
\]
\[
K_X = 4\pi a^2 \ell I_f H_w + \frac{\pi a^2 I_g \hbar^2}{m\rho_w \ell}.
\]

Therefore
\[
\kappa = \frac{K_X L^2}{T_X}
=
4\left(\frac{m c_{s,w} L}{\hbar}\right)^2
+
\frac{I_g}{I_f}\left(\frac{L}{\ell}\right)^2.
\]

The thin-wall matched-layer wall figure collapses exactly to
\[
W_{\rm wall}
=
\frac{4\rho_w^2 V_0^2 L^2}{\hbar^2 c_{s,w}^2 \ell^2}.
\]

And the Stage-41/42 fixed-point coupling is
\[
\Xi
=
\frac{g_\phi^2 I_1 L^2}{T_X}
=
W_{\rm wall}
\]
exactly on this explicit matched thin-wall branch.

So the wall figure of merit is not merely analogous to the support/source fixed-point coupling. It is the same object.

## Stage 245 — canonical `tanh` wall branch and natural local mouth closure

Choose the canonical smooth wall
\[
f(\xi) = \frac{1+\tanh \xi}{2},
\qquad
f'(\xi)=\frac12 \operatorname{sech}^2\xi,
\qquad
f''(\xi)=-\operatorname{sech}^2\xi \tanh\xi.
\]

Its exact shell moments are
\[
I_f = \int_{-\infty}^{+\infty} d\xi\, f'(\xi)^2 = \frac13,
\qquad
I_g = \int_{-\infty}^{+\infty} d\xi\, f''(\xi)^2 = \frac{4}{15},
\]
so
\[
\frac{I_g}{I_f} = \frac45.
\]

Therefore the explicit branch coefficients become
\[
T_X = \frac{\pi a^2 \ell \hbar^2}{3m\rho_w},
\]
\[
K_X = \frac{4\pi a^2}{15\ell m\rho_w}\left(5m^2 c_{s,w}^2 \ell^2 + \hbar^2\right),
\]
\[
J_1 = \frac{1}{3H_w},
\]
and
\[
\kappa
=
4\left(\frac{m c_{s,w} L}{\hbar}\right)^2
+
\frac45\left(\frac{L}{\ell}\right)^2.
\]

Using the natural local mouth closure
\[
K_m = \frac{T_X}{\ell},
\]
the Robin ratio is
\[
\eta = \frac{K_m L}{T_X} = \frac{L}{\ell}.
\]

So the first canonical explicit parent branch is controlled by the three dimensionless variables
\[
\chi_s := \frac{m c_{s,w} L}{\hbar},
\qquad
\Lambda_\ell := \frac{L}{\ell},
\qquad
\Upsilon_w := \frac{4\rho_w^2 V_0^2}{\hbar^2 c_{s,w}^2},
\]
with exact branch map
\[
\kappa = 4\chi_s^2 + \frac45 \Lambda_\ell^2,
\qquad
\eta = \Lambda_\ell,
\qquad
W_{\rm wall} = \Upsilon_w \Lambda_\ell^2.
\]

So the first explicit parent branch is no longer described by a long symbolic list. It is a three-parameter branch-placement problem.

## Stage 246 — explicit branch threshold surfaces

Insert the canonical branch formulas into the universal matched-branch theorem:
\[
W_{\rm wall} \le \frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)}
\quad\Rightarrow\quad
\text{fail},
\qquad
W_{\rm wall} \ge \frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)}
\quad\Rightarrow\quad
\text{succeed}.
\]

With
\[
\kappa = 4\chi_s^2 + \frac45\Lambda_\ell^2,
\qquad
\eta = \Lambda_\ell,
\qquad
W_{\rm wall} = \Upsilon_w \Lambda_\ell^2,
\]
this becomes the exact explicit-branch surfaces
\[
\Upsilon_w \le \Upsilon_{\rm fail}(\chi_s,\Lambda_\ell)
\quad\Rightarrow\quad
\text{fail},
\]
\[
\Upsilon_w \ge \Upsilon_{\rm suff}(\chi_s,\Lambda_\ell)
\quad\Rightarrow\quad
\text{succeed},
\]
where
\[
\Upsilon_{\rm fail}
=
\frac{Pe_{\rm req}}{\Lambda_\ell^2\,
\Delta_\infty\!\left(4\chi_s^2+\frac45\Lambda_\ell^2,\Lambda_\ell\right)},
\]
\[
\Upsilon_{\rm suff}
=
\frac{Pe_{\rm req}}{\Lambda_\ell^2\,
\Delta_0\!\left(4\chi_s^2+\frac45\Lambda_\ell^2,\Lambda_\ell\right)}.
\]

Two asymptotic regimes are already exact:

### Shell-gradient dominated
If
\[
\frac45 \Lambda_\ell^2 \gg 4\chi_s^2,
\]
then
\[
\Upsilon_{\rm fail} \sim \frac{2 Pe_{\rm req}}{\sqrt5\,\Lambda_\ell},
\qquad
\Upsilon_{\rm suff} \to \frac45\left(1+\frac{2}{\sqrt5}\right) Pe_{\rm req}.
\]

### Compression dominated
If
\[
4\chi_s^2 \gg \frac45 \Lambda_\ell^2,
\]
then
\[
\Upsilon_{\rm fail} \sim \frac{2 Pe_{\rm req}\chi_s}{\Lambda_\ell^2},
\]
\[
\Upsilon_{\rm suff} \sim \frac{4 Pe_{\rm req}\chi_s^2(\Lambda_\ell+2\chi_s)}{\Lambda_\ell^3},
\]
and for `chi_s >> Lambda_ell`,
\[
\Upsilon_{\rm suff} \sim \frac{8 Pe_{\rm req}\chi_s^3}{\Lambda_\ell^3}.
\]

So compression-dominated branches are much harder to drive across the universal success threshold than shell-gradient dominated branches.

## Net result after Stage 246

The support/source theorem gap is now much smaller than it was after Stage 241.

1. The abstract loading amplitude `g_phi` has been replaced by the explicit thin-wall relation `g_phi = V0/ell`.
2. The abstract gain `G_eq` has collapsed to the wall figure of merit `W_wall`.
3. The matched thin-wall GNLS shell reduction fixes `T_X`, `K_X`, `kappa`, and shows `W_wall = Xi` exactly.
4. On the canonical `tanh` wall with the natural local mouth closure, the whole branch is controlled by the three parent dimensionless variables
   `(chi_s, Lambda_ell, Upsilon_w)`.
5. The remaining PDE-side theorem gate is now

> compute the actual moving-throat branch placement in `(chi_s, Lambda_ell, Upsilon_w)` and compare it directly to the explicit threshold surfaces `Upsilon_fail`, `Upsilon_suff`.

That is the cleanest continuation point after Stage 241.


## Stage 247 — Family-1 reference-branch geometry map

The next reduction is to stop treating the wall-geometry ratio `Lambda_ell = L/ell` as free and evaluate it on the first carried constructive throat branch.

Use the balanced thin-layer-consistent Family-1 wall width
\[
\epsilon_r = \frac{\ell}{a} = \frac1{20},
\]
together with the carried preferred throat aspect ratio
\[
\frac{L}{a} = \frac{37}{20}.
\]

Then the explicit branch fixes
\[
\Lambda_\ell = \frac{L}{\ell} = \frac{(L/a)}{(\ell/a)} = 37.
\]

On the natural local mouth closure from Stage 245,
\[
\eta = \frac{K_m L}{T_X} = \frac{L}{\ell} = \Lambda_\ell,
\]
so the same reference branch fixes
\[
\eta = 37.
\]

So the first explicit moving-throat support branch is no longer free in the Robin geometry variable. It already sits at one concrete large-`eta` point.

## Stage 248 — healing-length lock and exact support-scale reduction

Now identify the active shell width with the local GNLS healing width,
\[
\ell = \ell_h = \frac{\hbar}{2 m c_{s,w}}.
\]

Then the support-scale variable from Stage 245 becomes
\[
\chi_s = \frac{m c_{s,w} L}{\hbar} = \frac{L}{2\ell} = \frac{\Lambda_\ell}{2}.
\]

Since Stage 247 fixed `Lambda_ell = 37`, the same reference branch fixes
\[
\chi_s = \frac{37}{2} = 18.5.
\]

Therefore
\[
\kappa
=
4\chi_s^2 + \frac45 \Lambda_\ell^2
=
\frac95 \Lambda_\ell^2
=
\frac{12321}{5}
=
2464.2,
\]
and
\[
\alpha := \sqrt{\kappa}
=
\frac{111}{\sqrt5}
\approx 49.6407091.
\]

So after the healing lock the explicit branch has fixed
\[
\Lambda_\ell = 37,\qquad
\eta = 37,\qquad
\chi_s = \frac{37}{2},\qquad
\kappa = \frac{12321}{5}.
\]

The only remaining Stage-245 control is now the wall-loading amplitude `Upsilon_w`.

## Stage 249 — explicit Family-1 threshold window and reduction to one wall-depth datum

With
\[
\Lambda_\ell = 37,\qquad
\eta = 37,\qquad
\kappa = \frac{12321}{5},
\]
the exact operator threshold functions evaluate numerically to
\[
\Delta_0\!\left(\frac{12321}{5},37\right)
\approx
1.73302079021525\times 10^{-4},
\]
\[
\Delta_\infty\!\left(\frac{12321}{5},37\right)
\approx
2.01447565540522\times 10^{-2}.
\]

Therefore the explicit branch thresholds become
\[
\Upsilon_{\rm fail}
\approx
0.0362605617972939\,Pe_{\rm req},
\]
\[
\Upsilon_{\rm suff}
\approx
4.21495341569977\,Pe_{\rm req}.
\]

Equivalently, the fixed-point coupling window is
\[
\Xi_{\rm fail}
=
\frac{Pe_{\rm req}}{\Delta_\infty}
\approx
49.6407091004953\,Pe_{\rm req},
\]
\[
\Xi_{\rm suff}
=
\frac{Pe_{\rm req}}{\Delta_0}
\approx
5770.27122609299\,Pe_{\rm req}.
\]

Now write the Family-1 wall depth as
\[
V_0 = \alpha_r \mu_*,
\qquad
\alpha_r = 10,
\]
so
\[
\Upsilon_w = \alpha_r^2 \Theta_w = 100 \Theta_w,
\]
with the one remaining microscopic datum
\[
\Theta_w := \frac{4\rho_w^2 \mu_*^2}{\hbar^2 c_{s,w}^2}.
\]

The explicit threshold window is then
\[
\Theta_{\rm fail}
\approx
3.62605617972939\times 10^{-4}\,Pe_{\rm req},
\]
\[
\Theta_{\rm suff}
\approx
4.21495341569977\times 10^{-2}\,Pe_{\rm req}.
\]

So on the explicit Family-1 / healing-locked reference branch the reduced support/source placement problem is now almost finished:

- the geometry ratio is fixed,
- the support/healing ratio is fixed,
- the operator scales are fixed,
- and the only remaining explicit branch datum is the wall-depth amplitude `Theta_w`.

## Net result after Stage 249

The continuation past Stage 241 is now substantially sharper.

1. The abstract parent loading amplitude `g_phi` is replaced by the explicit thin-wall law `g_phi = V0/ell`.
2. The support theorem collapses to the dimensionless wall figure of merit `W_wall`.
3. The matched thin-wall GNLS shell reduction fixes `T_X`, `K_X`, `kappa`, and proves `W_wall = Xi`.
4. On the canonical `tanh` wall with the natural local mouth closure, the whole branch reduces to the three control variables
   \[
   (\chi_s,\Lambda_\ell,\Upsilon_w).
   \]
5. On the explicit Family-1 / healing-locked reference branch, even those collapse further:
   \[
   \Lambda_\ell = 37,\quad
   \eta = 37,\quad
   \chi_s = 37/2,\quad
   \kappa = 12321/5,
   \]
   so the only remaining branch datum is the wall-depth amplitude `Theta_w`.

That means the next honest theorem gate is no longer to derive more reduced support/source algebra.

It is:

> extract the actual Family-1 wall-depth datum `Theta_w` (or equivalently `Upsilon_w`) from the real moving-throat branch and compare it directly to the explicit threshold window above.
# 5PN continuation notes — Stages 250–253 exact `\Lambda_{\rm EM}` refresh of the explicit Family-1 branch

This session continued the explicit Family-1 support/source branch, but with one correction carried all the way through:

> do **not** use the shorthand reference freeze `L/a = 37/20`; use the exact derived EM-branch aspect ratio instead.

The carried exact 2PN geometry relation is
\[
\Lambda_{\rm EM}
=
\frac{\sqrt{2}\,\pi}{x_{01}},
\qquad
x_{01}=2.40482555769577\ldots
\]
so numerically
\[
\Lambda_{\rm EM}\approx 1.847486577120128.
\]

The older use of
\[
\frac{L}{a}=\frac{37}{20}=1.85
\]
was only a convenient **reference-branch numerical freeze**. This session replaced that freeze everywhere in the explicit Family-1 branch formulas.

---

## Stage 250 — exact `\Lambda_{\rm EM}` geometry refresh

With the explicit thin-wall Family-1 lock
\[
\frac{\ell}{a}=\frac1{20},
\]
the exact carried aspect ratio gives
\[
\Lambda_\ell:=\frac{L}{\ell}
=
20\,\Lambda_{\rm EM}
=
\frac{20\sqrt2\,\pi}{x_{01}}
\approx 36.94973154240256,
\]
so
\[
\eta=\Lambda_\ell\approx 36.94973154240256.
\]

With the healing-width lock
\[
\chi_s=\frac{L}{2\ell},
\]
we get
\[
\chi_s
=
10\,\Lambda_{\rm EM}
=
\frac{10\sqrt2\,\pi}{x_{01}}
\approx 18.47486577120128.
\]

Then the explicit support parameter becomes
\[
\kappa
=
4\chi_s^2+\frac45\Lambda_\ell^2
=
\frac95\Lambda_\ell^2
=
\frac{1440\pi^2}{x_{01}^2}
\approx 2457.508789900114.
\]

Relative to the shorthand reference freeze, this is only a small correction:
\[
\frac{\Lambda_{\rm EM}-(37/20)}{37/20}
\approx
-1.358606962\times 10^{-3}.
\]

So the branch geometry moves only slightly, but from this point onward the exact EM-branch formula is the one that should be propagated.

---

## Stage 251 — refreshed explicit threshold window

Using the exact refreshed values of `\eta` and `\kappa`, the explicit Family-1 support theorem gives
\[
\Delta_0 \approx 1.737739392346995\times 10^{-4},
\qquad
\Delta_\infty \approx 2.0172162594593645\times 10^{-2}.
\]

Therefore
\[
\Upsilon_{\rm fail}
\approx
3.630989267026856\times 10^{-2}\,Pe_{\rm req},
\]
\[
\Upsilon_{\rm suff}
\approx
4.214953415699773\,Pe_{\rm req},
\]
and hence
\[
\Theta_{\rm fail}
\approx
3.630989267026856\times 10^{-4}\,Pe_{\rm req},
\]
\[
\Theta_{\rm suff}
\approx
4.214953415699773\times 10^{-2}\,Pe_{\rm req}.
\]

The important observation is:

- the **fail** threshold shifts slightly upward relative to the `37/20` freeze,
- the **success** threshold is unchanged to the displayed precision.

So the explicit branch theorem is stable under the exact geometry refresh.

---

## Stage 252 — refreshed wall-depth verdict

The explicit wall-depth extraction from the Family-1 wall profile is unchanged:

\[
\Theta_w^{(\chi)} \approx 4.06863235008162\,\lambda_\mu^2,
\qquad
\Theta_w^{(J)} \approx 0.927552032539308\,\lambda_\mu^2.
\]

Comparing these to the refreshed threshold window gives

\[
Pe_{\rm suff}^{(\chi)}
=
\frac{\Theta_w^{(\chi)}}{\Theta_{\rm suff}}
\approx
96.5285247264385\,\lambda_\mu^2,
\]
\[
Pe_{\rm fail}^{(\chi)}
=
\frac{\Theta_w^{(\chi)}}{\Theta_{\rm fail}}
\approx
11205.2998532081\,\lambda_\mu^2.
\]

For the conservative lower envelope,
\[
Pe_{\rm suff}^{(J)}
\approx
22.0062226330754\,\lambda_\mu^2,
\qquad
Pe_{\rm fail}^{(J)}
\approx
2554.54358117343\,\lambda_\mu^2.
\]

So the qualitative branch-level verdict survives exactly as before:

> on the first explicit Family-1 branch, wall-depth is still **not** the natural bottleneck.

The wall-depth side still succeeds automatically for modest quadrupole demand and fails only for anomalously large demand.

---

## Stage 253 — refreshed quadrupole-demand window in `\zeta_{\rm req}` and `\Pi_{\rm tr}`

The exact Family-1 support-ratio map now uses the refreshed
\[
\eta=\Lambda_\ell\approx 36.94973154240256,
\qquad
\kappa\approx 2457.508789900114.
\]

Solving
\[
y\tan y=\eta
\]
on the principal branch gives
\[
y_{F1}\approx 1.5294278190457656,
\]
and therefore
\[
A_{F1}
=
\frac{\kappa+\pi^2/4}{\kappa+y_{F1}^2}
\approx
1.0000521380385143.
\]

The hard Family-1 support ceiling becomes
\[
\zeta_{\max}^{(F1)}
=
A_{F1}\frac{\pi^2}{4}
\approx
2.467529745725936.
\]

At `\lambda_\mu=1`, the refreshed explicit branch windows are

\[
\zeta_{\rm suff}^{(\chi)}(1)\approx 2.466223429475074,
\qquad
\zeta_{\rm fail}^{(\chi)}(1)\approx 2.467529648745268,
\]
\[
\zeta_{\rm suff}^{(J)}(1)\approx 2.442576225820804,
\qquad
\zeta_{\rm fail}^{(J)}(1)\approx 2.467527879753313.
\]

So the natural shell-weighted guaranteed-success threshold still sits only
\[
\zeta_{\max}^{(F1)}-\zeta_{\rm suff}^{(\chi)}(1)
\approx
1.306250899\times 10^{-3}
\]
below the hard Family-1 ceiling, while the guaranteed-failure threshold is essentially saturated at the ceiling.

In the unblocked product language `(\epsilon_{\rm blk}=0)`,
\[
\Pi_{\rm suff}^{(\chi)}/C_{\rm mix}\approx 3.466223429475074,
\]
\[
\Pi_{\rm fail}^{(\chi)}/C_{\rm mix}\approx 3.467529648745268,
\]
\[
\Pi_{\max}^{(F1)}/C_{\rm mix}\approx 3.467529745725936.
\]

The corresponding exact blocking ceiling is
\[
\epsilon_{\rm blk}^{\rm crit}
=
\frac{1}{\zeta_{\max}^{(F1)}}
\approx
0.405263604919909.
\]

So the explicit-branch conclusion is unchanged even in the exact `\Lambda_{\rm EM}` refresh:

> the Family-1 wall-depth side still pushes the guaranteed-success threshold extremely close to the hard support ceiling, so the unresolved issue remains the **selected quadrupole demand side**, not wall-depth supply.

---

## Net result after Stage 253

Replacing the shorthand
\[
L/a = 37/20
\]
by the exact carried relation
\[
L/a = \Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}
\]
does **not** change the qualitative support/source verdict.

It does three useful things:

1. it removes the last explicit use of a numerical freeze where an exact carried equation exists,
2. it refreshes the fail-side threshold numbers consistently,
3. and it confirms that the explicit Family-1 subprogram is still bottlenecked by the **quadrupole-demand branch**.

So the next honest theorem gate remains the same in substance, but now with the corrected geometry already folded in:

> compare the final selected quadrupole branch demand directly to the refreshed explicit Family-1 ceiling in `\zeta_{\rm req}` or `\Pi_{\rm tr}` language.
# 5PN continuation notes — Stages 254–257 exact Family-1 demand robustness after the `\Lambda_{\rm EM}` refresh

This session continued directly from the exact `\Lambda_{\rm EM}` refresh of the explicit Family-1 support/source branch.
The next clean move was **not** to re-assume the canonical `3/4 + 1/4` isotropic split immediately. Instead, the analysis kept a generic isotropic one-pole conservative quadrupole carrier all the way through the explicit Family-1 demand ceiling and only specialized to the actual isotropic branch at the end.

That avoids the main oversimplification risk that had been left open at the end of Stage 253.

The upshot is stronger than before:

> the refreshed Family-1 support/source side is not merely compatible with the canonical `c_{\rm pole}=1/4` branch. It remains safe on a **large exact interval** of isotropic one-pole conservative branches, and the actual isotropic branch stays automatic throughout the full admissible blocked regime.

---

## Stage 254 — exact Family-1 admissible one-pole window

Take a generic isotropic one-pole conservative quadrupole carrier
\[
\widehat Y_Q^{\rm cons}(\omega)
=
 c_{\rm geom}
+
\frac{c_{\rm pole}}{1-\omega^2/\Omega_Q^2},
\qquad
c_{\rm geom}+c_{\rm pole}=1.
\]
Then the exact reduced demand variables are
\[
\rho_\alpha = \frac{1}{1-c_{\rm pole}},
\]
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]
and the selected-branch product variable collapses exactly to
\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\rho_\alpha.
\]
So the exact Family-1 demand ceiling translates to one sharp conservative one-pole condition:
\[
 c_{\rm pole}
<
 c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})
=
\frac{\zeta_{\max}^{(F1)}(1-\epsilon_{\rm blk})}{1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}}.
\]
With the exact refreshed value
\[
\zeta_{\max}^{(F1)}\approx 2.4675297457259358,
\]
the unblocked (`\epsilon_{\rm blk}=0`) hard ceiling is
\[
 c_{\rm pole,max}^{(F1)}(0)
\approx 0.7116102605226109,
\qquad
c_{\rm geom,min}^{(F1)}(0)
\approx 0.2883897394773891.
\]
The shell-weighted guaranteed-success threshold is only slightly lower:
\[
 c_{\rm pole,suff}^{(\chi)}(0)
\approx 0.7115015750293280.
\]
So the explicit Family-1 branch leaves a very large exact admissible window in one-pole conservative language.

---

## Stage 255 — exact geometry-lane contamination ceiling

Now map the Family-1 one-pole ceiling back into the geometry-lane contamination variables using the Stage-75 exact one-pole relation
\[
 c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2}.
\]
Then the exact Family-1 safe region is
\[
1+\epsilon_4
<
4 c_\star (1+\epsilon_2)^2,
\]
with `c_\star` taken either as the hard ceiling `c_{\rm pole,max}` or the guaranteed-success ceiling `c_{\rm pole,suff}`.

Equivalently,
\[
\epsilon_{4,\max}(\epsilon_2;c_\star)
=
4 c_\star (1+\epsilon_2)^2 - 1.
\]
On the refreshed unblocked Family-1 branch this gives
\[
\epsilon_4 < 1.8464410420904435
\qquad (\epsilon_2=0,\; \text{hard ceiling}),
\]
while on the principal physical slice `1+\epsilon_2>0`
\[
\epsilon_2 > -0.4072809255389385
\qquad (\epsilon_4=0,\; \text{hard ceiling}).
\]
So the explicit support/source side remains tolerant to order-one geometry-lane contamination around the canonical point.
This is much stronger than merely checking the canonical `1/4` split by itself.

---

## Stage 256 — actual isotropic branch is automatic throughout the full blocked regime

Only after the generic window is established do we specialize to the actual isotropic one-pole branch
\[
 c_{\rm pole}=\frac14,
\qquad
c_{\rm geom}=\frac34.
\]
Then
\[
\rho_\alpha^{(\rm act)}=\frac43,
\qquad
\frac{\Pi_{\rm tr}^{(\rm act)}}{C_{\rm mix}}=\frac43,
\]
and the exact blocked-regime demand is
\[
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk})
=
\frac{1}{3-2\epsilon_{\rm blk}}.
\]
The admissible Family-1 blocked regime is
\[
0\le \epsilon_{\rm blk} < \epsilon_{\rm blk}^{\rm crit}
=
\frac{1}{\zeta_{\max}^{(F1)}}
\approx 0.40526360491990934.
\]
At the worst admissible blocking value,
\[
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk}^{\rm crit})
=
\frac{\zeta_{\max}^{(F1)}}{3\zeta_{\max}^{(F1)}-2}
\approx 0.45673095573242554,
\]
so the exact hard margin is still
\[
\zeta_{\max}^{(F1)}-
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk}^{\rm crit})
\approx 2.0107987899935102.
\]
Because `\zeta_{\max}^{(F1)} > 1`, the actual isotropic branch also stays in the symmetric-lowest-twin-safe regime throughout the full admissible blocked interval.
So the refreshed Family-1 support/source side is automatic not only at `\epsilon_{\rm blk}=0`, but on the whole blocked branch.

---

## Stage 257 — universal twin-safe strip and exact Family-1 extension

There is also a universal exact theorem independent of the explicit Family-1 details.
From
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]
solving `\zeta_{\rm req}=1` gives the exact, blocking-independent boundary
\[
 c_{\rm pole}=\frac12.
\]
So every isotropic one-pole conservative carrier with
\[
 c_{\rm pole}\le \frac12
\qquad\Longleftrightarrow\qquad
\rho_\alpha\le 2
\]
already lies in the universal symmetric-lowest-twin-safe strip.

The refreshed Family-1 branch extends this universal strip to the exact larger window
\[
 c_{\rm pole} < 0.7116102605226109,
\qquad
\rho_\alpha < 3.4675297457259358.
\]
So the extension beyond the universal twin-safe strip is
\[
\Delta c_{\rm pole}^{(F1)}
\approx 0.2116102605226109,
\qquad
\Delta \rho_\alpha^{(F1)}
\approx 1.4675297457259358.
\]
The actual isotropic branch sits at
\[
 c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\]
so it lies deeply inside both the universal strip and the larger exact Family-1 extension.

---

## Net result after Stage 257

The support/source side is now more robustly understood than before.

What is exact after this pass:

1. the exact `\Lambda_{\rm EM}`-refreshed Family-1 ceiling translated into generic isotropic one-pole variables,
2. the exact admissible window in `c_{\rm pole}` / `c_{\rm geom}` language,
3. the exact geometry-lane safe region in `(\epsilon_2,\epsilon_4)`,
4. the exact blocked-regime theorem for the actual isotropic `c_{\rm pole}=1/4` branch,
5. the exact universal twin-safe strip `c_{\rm pole}\le 1/2`,
6. and the exact extension supplied by the explicit Family-1 branch beyond that strip.

So the remaining theorem gap is **not** explicit support/source supply on the isotropic one-pole branch.
Even without re-assuming the canonical `1/4` split until the end, the refreshed Family-1 side stays safely non-bottlenecked over a large exact neighborhood of that branch.

The live work remains where the compact moving-throat ledger said it should be:

> branch realization and outgoing DtN selection/normalization on the actual moving-throat grouped-`P2` branch, not wall-depth or support-source supply on the explicit Family-1 side.
# 5PN continuation notes — Stages 258–261 exact one-pole regime bridge back to the selected moving-throat grouped-`P2` branch

This session continued directly from the exact `\Lambda_{\rm EM}` refresh and the caution not to re-freeze the canonical `c_{\rm pole}=1/4` split too early.

The next clean move was to propagate the exact support-regime theorems back into the **generic isotropic one-pole conservative grouped-`P2` carrier** and then into the geometry-lane contamination variables `\epsilon_2,\epsilon_4`, while keeping the blocking variable `\epsilon_{\rm blk}` explicit instead of silently setting it to zero again.

The payoff is sharper than before:

1. the entire mixed / twin / non-twin support classification is exact in one-pole language,
2. the same split becomes an exact phase portrait in `\bigl(\epsilon_2,\epsilon_4\bigr)`,
3. blocking changes the branch demand with a sign flip at the universal threshold `c_{\rm pole}=1/2`,
4. and on the exact `\Lambda_{\rm EM}`-refreshed Family-1 branch the admissible non-twin corridor actually **widens** with blocking even though the actual isotropic point stays twin-safe throughout.

So the continuation point is now much cleaner: if the true moving-throat grouped-`P2` branch ever leaves the universal twin-safe strip, we know the exact `c_{\rm pole}` and `\bigl(\epsilon_2,\epsilon_4\bigr)` corridor in which explicit Family-1 rescue is still possible.

---

## Stage 258 — exact one-pole selected-branch regime split

Take the generic isotropic one-pole conservative carrier
\[
\widehat Y_Q^{\rm cons}(\omega)
=
 c_{\rm geom}
+
\frac{c_{\rm pole}}{1-\omega^2/\Omega_Q^2},
\qquad
c_{\rm geom}+c_{\rm pole}=1.
\]
Then the exact selected-branch demand variables are
\[
\rho_\alpha = \frac{1}{1-c_{\rm pole}},
\qquad
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]
and, exactly as in Stage 254,
\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \rho_\alpha.
\]
So the support-regime split is independent of blocking in one-pole language:
\[
\Pi_{\rm tr}\le C_{\rm mix}
\iff
\rho_\alpha\le1
\iff
c_{\rm pole}\le0,
\]
\[
C_{\rm mix}<\Pi_{\rm tr}\le2C_{\rm mix}
\iff
1<\rho_\alpha\le2
\iff
0<c_{\rm pole}\le\frac12,
\]
\[
\Pi_{\rm tr}>2C_{\rm mix}
\iff
\rho_\alpha>2
\iff
c_{\rm pole}>\frac12.
\]
So mixed-only sufficiency is already impossible on every physical positive-pole branch.
The universal symmetric-lowest-twin strip is exactly
\[
0<c_{\rm pole}\le\frac12.
\]
The actual isotropic branch remains
\[
c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\]
so it lies strictly inside the twin-safe regime.

On the exact `\Lambda_{\rm EM}`-refreshed Family-1 branch at `\epsilon_{\rm blk}=0`, the hard ceiling is still
\[
c_{\rm pole}<c_{\rm pole,max}^{(F1)}(0)
\approx 0.7116102605226109,
\]
so there remains a finite non-twin corridor
\[
\frac12<c_{\rm pole}<c_{\rm pole,max}^{(F1)}(0).
\]

---

## Stage 259 — exact geometry-lane regime surface in `\bigl(\epsilon_2,\epsilon_4\bigr)`

Using the exact Stage-75 obstruction formula
\[
c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\]
the selected-branch demand becomes
\[
\rho_\alpha
=
\frac{4(1+\epsilon_2)^2}{4\epsilon_2^2+8\epsilon_2-\epsilon_4+3},
\qquad
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\rho_\alpha.
\]
So the exact regime surfaces are now explicit in geometry-contamination language:

### Mixed-only boundary
\[
c_{\rm pole}=0
\iff
1+\epsilon_4=0.
\]

### Universal lowest-twin boundary
\[
c_{\rm pole}=\frac12
\iff
1+\epsilon_4 = 2(1+\epsilon_2)^2.
\]
Equivalently,
\[
c_{\rm pole}-\frac12
=
\frac{(1+\epsilon_4)-2(1+\epsilon_2)^2}{4(1+\epsilon_2)^2},
\]
so the sign of
\[
(1+\epsilon_4)-2(1+\epsilon_2)^2
\]
exactly determines whether the branch is twin-safe or non-twin.

At `\epsilon_{\rm blk}=0`, the exact `\Lambda_{\rm EM}`-refreshed Family-1 upper ceiling becomes
\[
1+\epsilon_4
<
4c_{\rm pole,max}^{(F1)}(0)
(1+\epsilon_2)^2
\approx
2.8464410420904435
(1+\epsilon_2)^2.
\]
So the exact unblocked admissible non-twin strip is
\[
2(1+\epsilon_2)^2
<
1+\epsilon_4
<
4c_{\rm pole,max}^{(F1)}(0)(1+\epsilon_2)^2.
\]
The actual isotropic grouped-`P2` point
\[
\epsilon_2=\epsilon_4=0
\]
has exact twin margin
\[
2(1+0)^2-(1+0)=1,
\]
so it sits safely below the non-twin surface.

---

## Stage 260 — exact blocking monotonicity and asymmetry demand

Keeping `\epsilon_{\rm blk}` explicit, the exact support demand is
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}}.
\]
Its exact derivatives are
\[
\frac{\partial \zeta_{\rm req}}{\partial c_{\rm pole}}
=
\frac{1-\epsilon_{\rm blk}}{\bigl[1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}\bigr]^2}>0,
\]
\[
\frac{\partial \zeta_{\rm req}}{\partial \epsilon_{\rm blk}}
=
\frac{c_{\rm pole}(1-2c_{\rm pole})}{\bigl[1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}\bigr]^2}.
\]
So the branch ordering in `c_{\rm pole}` is exact, and blocking has a sign flip at the universal threshold `c_{\rm pole}=1/2`:
\[
c_{\rm pole}<\frac12
\Rightarrow
\partial_{\epsilon_{\rm blk}}\zeta_{\rm req}>0,
\]
\[
c_{\rm pole}=\frac12
\Rightarrow
\partial_{\epsilon_{\rm blk}}\zeta_{\rm req}=0,
\]
\[
c_{\rm pole}>\frac12
\Rightarrow
\partial_{\epsilon_{\rm blk}}\zeta_{\rm req}<0.
\]
So blocking hurts the twin-safe side but helps the non-twin side.

The exact excess beyond the symmetric-twin threshold is
\[
\Delta_\zeta:=\zeta_{\rm req}-1
=
\frac{(1-\epsilon_{\rm blk})(2c_{\rm pole}-1)}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}}.
\]
Therefore:
- `\Delta_\zeta<0` in the twin-safe strip,
- `\Delta_\zeta=0` on the universal boundary `c_{\rm pole}=1/2`,
- `\Delta_\zeta>0` on the non-twin side.

The actual isotropic branch stays
\[
\zeta_{\rm req}^{(act)}=\frac{1}{3-2\epsilon_{\rm blk}},
\]
so blocking increases its demand but never drives it out of the twin-safe regime on the admissible Family-1 branch.

---

## Stage 261 — exact blocked Family-1 non-twin corridor

Let `\zeta_{\max}^{(F1)}` be the exact refreshed Family-1 support ceiling in support-ratio language. Solving
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk}) = \zeta_{\max}^{(F1)}
\]
for `c_{\rm pole}` gives the exact blocked hard ceiling
\[
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})
=
\frac{\zeta_{\max}^{(F1)}(1-\epsilon_{\rm blk})}{1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}}.
\]
Because `\zeta_{\max}^{(F1)}>1`, its derivative is strictly positive:
\[
\frac{d c_{\rm pole,max}^{(F1)}}{d\epsilon_{\rm blk}}
=
\frac{\zeta_{\max}^{(F1)}(\zeta_{\max}^{(F1)}-1)}{\bigl[1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}\bigr]^2}>0.
\]
So the hard Family-1 pole ceiling grows with blocking.

Exact endpoint values:
\[
c_{\rm pole,max}^{(F1)}(0)
=
\frac{\zeta_{\max}^{(F1)}}{1+\zeta_{\max}^{(F1)}},
\qquad
\epsilon_{\rm blk}^{crit} = \frac{1}{\zeta_{\max}^{(F1)}},
\qquad
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk}^{crit})=1.
\]
So the admissible blocked non-twin corridor is exactly
\[
\frac12<c_{\rm pole}<c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk}),
\]
and its width
\[
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})-\frac12
\]
increases monotonically from
\[
0.2116102605226109\quad (\epsilon_{\rm blk}=0)
\]
to
\[
0.5\quad (\epsilon_{\rm blk}=\epsilon_{\rm blk}^{crit}).
\]
In geometry-contamination language the exact blocked non-twin strip becomes
\[
2(1+\epsilon_2)^2
<
1+\epsilon_4
<
4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})(1+\epsilon_2)^2.
\]
The upper coefficient therefore rises from
\[
4c_{\rm pole,max}^{(F1)}(0)
\approx 2.8464410420904435
\]
to
\[
4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk}^{crit})=4.
\]

A final exact invariance is that the maximal asymmetry demand on the blocked Family-1 ceiling stays fixed:
\[
\Delta_{\zeta,\max}
=
\zeta_{\rm req}(c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk});\epsilon_{\rm blk})-1
=
\zeta_{\max}^{(F1)}-1.
\]
So blocking widens the admissible `c_{\rm pole}` / `\bigl(\epsilon_2,\epsilon_4\bigr)` corridor, but it does **not** raise the maximal support asymmetry demanded at the Family-1 hard ceiling, because that ceiling is defined at fixed `\zeta_{\max}^{(F1)}`.

---

## Net result after Stage 261

The exact one-pole / geometry-lane branch picture is now much sharper.

What is exact after this pass:

1. the selected-branch support regime split is exact in `c_{\rm pole}` language,
2. the same split is an exact phase portrait in `\bigl(\epsilon_2,\epsilon_4\bigr)`,
3. blocking hurts the twin-safe side but helps the non-twin side, with the sign flip fixed at `c_{\rm pole}=1/2`,
4. the exact `\Lambda_{\rm EM}`-refreshed Family-1 non-twin corridor widens with blocking,
5. and the maximal asymmetry demand on the Family-1 ceiling stays fixed at `\zeta_{\max}^{(F1)}-1`.

So the live theorem gap is now even narrower than it was at Stage 257.
It is no longer about whether the explicit support/source side can tolerate one-pole deformation in the abstract.
It is:

> if the actual moving-throat grouped-`P2` branch ever leaves the universal twin-safe strip, does its exact `c_{\rm pole}` / `\bigl(\epsilon_2,\epsilon_4\bigr)` placement remain inside the widened blocked Family-1 corridor where explicit support rescue is still possible?
# 5PN continuation notes — Stages 262–264

This batch takes the exact blocked Family-1 corridor from Stages 258–261 and applies it directly to the **actual selected moving-throat grouped-`P2` branch** rather than to a generic one-pole carrier.

The key input from the earlier 5PN notes is the actual grouped-`P2` carrier formula

\[
c_{\rm pole}=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\]

with the actual isotropic branch at

\[
\epsilon_2=\epsilon_4=0
\quad\Longrightarrow\quad
c_{\rm pole}=\frac14.
\]

So the continuation point after Stage 261 is no longer “what happens for a generic one-pole branch?” It is:

> what does the **actual** moving-throat grouped-`P2` branch do inside the exact blocked Family-1 corridor, and how robust is that branch under the first weak-anisotropy contamination?

---

## Stage 262 — exact actual grouped-`P2` branch map

Using

\[
c_{\rm pole}=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\qquad
\zeta_{\rm req}
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]

we get the exact blocked demand map in the **actual** grouped-`P2` variables:

\[
\rho_\alpha
=
\frac{4(1+\epsilon_2)^2}{4\epsilon_2^2+8\epsilon_2-\epsilon_4+3},
\]

\[
\zeta_{\rm req}
=
\frac{1+\epsilon_4}
{4(1+\epsilon_2)^2(1-\epsilon_{\rm blk})-(1-2\epsilon_{\rm blk})(1+\epsilon_4)}.
\]

At the actual isotropic point,

\[
\epsilon_2=\epsilon_4=0,
\qquad
c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac{1}{3-2\epsilon_{\rm blk}}.
\]

So the exact actual isotropic grouped-`P2` branch carries precisely the same demand point already isolated earlier in generic one-pole language.

### Exact twin-safety numerator

The universal lowest-twin boundary becomes the exact branch numerator

\[
M_{\rm twin}
:=
2(1+\epsilon_2)^2-(1+\epsilon_4).
\]

It satisfies

\[
c_{\rm pole}-\frac12
=
-\frac{M_{\rm twin}}{4(1+\epsilon_2)^2},
\]

and

\[
\zeta_{\rm req}-1
=
-\frac{2(1-\epsilon_{\rm blk})M_{\rm twin}}
{4(1+\epsilon_2)^2(1-\epsilon_{\rm blk})-(1-2\epsilon_{\rm blk})(1+\epsilon_4)}.
\]

So the actual grouped-`P2` branch is

- twin-safe iff `M_twin >= 0`,
- exactly on the twin boundary iff `M_twin = 0`,
- non-twin iff `M_twin < 0`.

### Exact contamination monotonicity

The geometry-lane contamination acts asymmetrically:

\[
\frac{\partial c_{\rm pole}}{\partial \epsilon_2}
=
-\frac{1+\epsilon_4}{2(1+\epsilon_2)^3}<0,
\qquad
\frac{\partial c_{\rm pole}}{\partial \epsilon_4}
=
\frac{1}{4(1+\epsilon_2)^2}>0,
\]

and on every admissible blocked branch,

\[
\frac{\partial \zeta_{\rm req}}{\partial \epsilon_2}<0,
\qquad
\frac{\partial \zeta_{\rm req}}{\partial \epsilon_4}>0.
\]

So positive `epsilon_2` contamination **softens** the support demand, while positive `epsilon_4` contamination **hardens** it.

That is already more informative than the generic one-pole picture.

---

## Stage 263 — exact blocked Family-1 corridor on the actual branch

The exact blocked Family-1 ceiling from Stage 261 is

\[
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})
=
\frac{\zeta_{\max}^{(F1)}(1-\epsilon_{\rm blk})}
{1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}}.
\]

Rewriting directly in the actual grouped-`P2` variables gives the second exact numerator

\[
M_{F1}
:=
4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})(1+\epsilon_2)^2-(1+\epsilon_4).
\]

Then

- `M_F1 > 0` is exact Family-1 admissibility,
- `M_twin >= 0` is exact lowest-twin sufficiency.

So the whole blocked corridor is now expressed directly in the actual branch variables.

### Actual isotropic margins

At the physical isotropic point,

\[
M_{\rm twin}(0,0)=1,
\qquad
M_{F1}(0,0)=4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})-1.
\]

So the actual isotropic branch sits exactly **one full unit** below the universal twin/non-twin boundary.

And because `zeta_max^(F1) > 1`,

\[
\frac{d}{d\epsilon_{\rm blk}}M_{F1}(0,0)
=
\frac{4\zeta_{\max}^{(F1)}(\zeta_{\max}^{(F1)}-1)}{\bigl(1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}\bigr)^2}>0,
\]

the exact Family-1 admissibility margin at the isotropic point **grows** with blocking.

On the exact `Lambda_EM`-refreshed branch we get numerically:

- twin margin at the isotropic point: `1`,
- Family-1 admissibility margin at `epsilon_blk = 0`:
  `4 c_pole,max^(F1)(0) - 1 ≈ 1.84644104209044`,
- Family-1 admissibility margin at the critical blocking ceiling:
  `3`.

So the actual isotropic grouped-`P2` branch is not merely admissible. It sits comfortably inside the exact blocked corridor.

---

## Stage 264 — exact second-order weak-anisotropy tolerance theorem

Stage 78 in the notes says the first nonzero geometry contamination enters only at

\[
\epsilon_2,\epsilon_4 = O(\chi^2).
\]

Write this as

\[
\epsilon_2 = a_2 \chi^2,
\qquad
\epsilon_4 = a_4 \chi^2,
\qquad
y := \chi^2.
\]

Then the exact twin-safe and Family-1 conditions become explicit quadratics in `y`.

### Exact twin-safe quadratic

\[
M_{\rm twin}(y)
=
2a_2^2 y^2 + (4a_2-a_4)y + 1.
\]

So the actual branch remains twin-safe iff `M_twin(y) >= 0`.

### Exact Family-1 admissibility quadratic

\[
M_{F1}(y)
=
4c_{\rm pole,max}^{(F1)}a_2^2 y^2
+
(8c_{\rm pole,max}^{(F1)}a_2-a_4)y
+
(4c_{\rm pole,max}^{(F1)}-1).
\]

So the exact blocked Family-1 corridor persists iff `M_F1(y) > 0`.

### Three immediate theorem-level consequences

1. **Finite safe neighborhood around the isotropic point**

   Because
   \[
   M_{\rm twin}(0)=1,
   \qquad
   M_{F1}(0)=4c_{\rm pole,max}^{(F1)}-1>0,
   \]
   every weak-anisotropy branch has a finite safe neighborhood around `chi = 0`.

2. **Initial drift controlled only by `a4 - 2 a2`**

   The exact initial slopes are
   \[
   \left.\frac{dc_{\rm pole}}{d(\chi^2)}\right|_{\chi=0}
   =
   \frac{a_4-2a_2}{4},
   \]
   \[
   \left.\frac{d\zeta_{\rm req}}{d(\chi^2)}\right|_{\chi=0}
   =
   \frac{4(1-\epsilon_{\rm blk})(a_4-2a_2)}{(3-2\epsilon_{\rm blk})^2}.
   \]

   So:
   - `a4 - 2 a2 > 0` pushes the actual branch toward larger pole weight and larger support demand,
   - `a4 - 2 a2 < 0` pushes it toward smaller pole weight and smaller support demand.

3. **Exit can occur only through exact quadratic roots**

   Any eventual loss of twin-safety or Family-1 admissibility must occur through the corresponding positive roots of `M_twin(y)` or `M_F1(y)`. It cannot happen arbitrarily close to the actual isotropic point.

So the Stage-78 `O(chi^2)` statement is now turned into an explicit algebraic tolerance theorem.

---

## Net result after Stages 262–264

The “next move” after the generic one-pole corridor is now complete.

1. The exact blocked support corridor has been rewritten directly in the **actual** selected moving-throat grouped-`P2` branch variables `epsilon_2, epsilon_4`.
2. The actual isotropic branch sits at
   \[
   c_{\rm pole}=\frac14,
   \qquad
   M_{\rm twin}=1,
   \qquad
   M_{F1}=4c_{\rm pole,max}^{(F1)}-1>0,
   \]
   so it is safely inside both the universal twin-safe strip and the exact blocked Family-1 corridor.
3. The first weak-anisotropy contamination has been promoted from the qualitative statement `O(chi^2)` to exact quadratics in `chi^2`.

So the next honest theorem gate is now very sharp:

> extract the actual weak-anisotropy coefficients `(a_2,a_4)` — equivalently the first concrete `l=0 <-> l=2` induced geometry-mixing law — from the moving-throat branch, then test those coefficients against the exact quadratic twin-safe and Family-1 admissibility conditions above.
# 5PN continuation — Stages 265–267

## Goal

Take the Stage-78 weak-anisotropy statement

\[
\epsilon_2,\epsilon_4 = O(\chi^2)
\]

and turn it into the first concrete induced-geometry mixing law for the actual grouped-`P2` branch. The point is to stop talking about a generic `O(chi^2)` correction and instead extract the exact coefficients

\[
\epsilon_2 = a_2 \chi^2,
\qquad
\epsilon_4 = a_4 \chi^2,
\]

for the first explicit `l=0 \leftrightarrow l=2` mechanism, then test those coefficients against the exact twin-safe and Family-1 corridor conditions from the later grouped-`P2` support analysis.

---

## Stage 265 — exact coefficient extraction from the Stage-78 Schur complement

Start from the exact Stage-78 reduced mixed scalar-geometry / grouped-`P2` model,

\[
D_q(\omega)=K_{\rm stat}+\frac{K_{\rm pole}}{1-\omega^2/\Omega_Q^2},
\qquad
D_g(\omega)=G_0+G_2\omega^2+G_4\omega^4+O(\omega^6),
\]

with weak anisotropy generating a bilinear mixing term

\[
\chi M_0\,qg.
\]

Integrating out the scalar/geometry lane gives

\[
D_{\rm eff}(\omega)
=
D_q(\omega)-\frac{\chi^2 M_0^2}{D_g(\omega)}.
\]

Expanding through `O(omega^4)` yields the exact contamination moments

\[
K_{(g,2)}^{\rm eff} = \chi^2\frac{M_0^2 G_2}{G_0^2},
\qquad
K_{(g,4)}^{\rm eff} = \chi^2\frac{M_0^2(G_0G_4-G_2^2)}{G_0^3}.
\]

Therefore the grouped-`P2` obstruction coefficients are exactly

\[
a_2 = \frac{\Omega_Q^2 M_0^2 G_2}{K_{\rm pole} G_0^2},
\qquad
a_4 = \frac{\Omega_Q^4 M_0^2 (G_0G_4-G_2^2)}{K_{\rm pole} G_0^3}.
\]

So the Stage-78 `O(chi^2)` statement already has a fully explicit coefficient-level form.

### First concrete scalar-lane specialization

Take the scalar/geometry lane itself to be one effective contact plus one scalar pole,

\[
D_g(\omega)=G_c+\frac{G_p}{1-\omega^2/\Omega_g^2}.
\]

Define

\[
G_0 := G_c+G_p,
\qquad
c_g := \frac{G_p}{G_0},
\qquad
r := \frac{\Omega_Q^2}{\Omega_g^2},
\qquad
\mu_{\rm mix}:=\frac{M_0^2}{K_{\rm pole}G_0}.
\]

Then the exact induced-mixing coefficients collapse to

\[
a_2 = \mu_{\rm mix}\,r\,c_g,
\qquad
a_4 = \mu_{\rm mix}\,r^2 c_g(1-c_g).
\]

The initial grouped-`P2` support-demand drift is therefore controlled by the single combination

\[
a_4-2a_2
=
\mu_{\rm mix}\,r\,c_g\,[r(1-c_g)-2].
\]

So the first actual `l=0 \leftrightarrow l=2` mechanism is already much more rigid than a generic perturbative correction: it is governed by one amplitude `mu_mix`, one pole-ratio `r`, and one scalar-lane pole fraction `c_g`.

---

## Stage 266 — exact corridor theorem for the one-pole scalar lane

Introduce the compact control variable

\[
u := r(1-c_g).
\]

Then the exact weak-anisotropy corridor tests become

\[
M_{\rm twin}(y)=1+\alpha(4-u)y+2\alpha^2 y^2,
\qquad
\alpha:=\mu_{\rm mix}rc_g,
\qquad y:=\chi^2,
\]

and, for a general Family-1 ceiling `c_*`,

\[
M_{F1}(y)=(4c_*-1)+\alpha(8c_*-u)y+4c_*\alpha^2y^2.
\]

This gives a complete exact threshold ladder.

### Initial-drift thresholds

\[
u > 2
\quad\Longleftrightarrow\quad
\text{the actual grouped-`P2` support demand grows initially},
\]

\[
u > 4
\quad\Longleftrightarrow\quad
\text{the universal twin-safe margin shrinks initially},
\]

\[
u > 8c_*
\quad\Longleftrightarrow\quad
\text{the exact Family-1 margin shrinks initially}.
\]

### Actual failure thresholds

The exact twin discriminant is

\[
\Delta_{\rm twin}=\alpha^2[(u-4)^2-8],
\]

so a positive twin-boundary root exists only if

\[
u \ge 4+2\sqrt2.
\]

Likewise the exact Family-1 discriminant is

\[
\Delta_{F1}=\alpha^2[(u-8c_*)^2-16c_*(4c_*-1)],
\]

so an actual Family-1 failure requires

\[
u \ge 8c_* + 4\sqrt{c_*(4c_*-1)}.
\]

That is a strong theorem. The first induced geometry-mixing mechanism threatens the isotropic branch only when the scalar lane is both sufficiently pole-light and sufficiently faster than the grouped quadrupole pole.

### Pure-pole scalar lane

If `c_g=1`, then

\[
a_4=0,
\qquad
u=0,
\]

so the branch sits far below every drift/failure threshold. A pure scalar pole can soften the demand, but by itself it can never drive the actual grouped-`P2` branch out of either the universal twin-safe strip or the exact Family-1 corridor.

---

## Stage 267 — exact numerical thresholds on the Lambda_EM-refreshed Family-1 branch

Using the exact refreshed unblocked Family-1 ceiling,

\[
c_* = c_{\rm pole,max}^{(F1)}(0) \approx 0.7116102605,
\]

the hard threshold values are

\[
u_{\rm twin}^{\rm fail} = 4+2\sqrt2 \approx 6.8284271247,
\]

\[
u_{F1}^{\rm fail} = 8c_* + 4\sqrt{c_*(4c_*-1)} \approx 10.2779821110.
\]

Since `u = r(1-c_g)`, the required pole-ratio thresholds scale as

\[
r_{\rm crit} = \frac{u_{\rm crit}}{1-c_g}.
\]

Representative values on the exact Lambda_EM-refreshed hard ceiling are:

- `c_g = 0.25`:
  \[
  r_{\rm twin}^{\rm fail} \approx 9.10457,
  \qquad
  r_{F1}^{\rm fail} \approx 13.70398;
  \]
- `c_g = 0.50`:
  \[
  r_{\rm twin}^{\rm fail} \approx 13.65685,
  \qquad
  r_{F1}^{\rm fail} \approx 20.55596;
  \]
- `c_g = 0.75`:
  \[
  r_{\rm twin}^{\rm fail} \approx 27.31371,
  \qquad
  r_{F1}^{\rm fail} \approx 41.11193.
  \]

So unless the scalar geometry pole is *much* faster than the grouped quadrupole pole, the first actual `l=0 \leftrightarrow l=2` mixing mechanism does not by itself drive the real isotropic branch out of the universal twin-safe or exact Family-1 strips.

---

## Best current reading after Stages 265–267

The next theorem gate is now much narrower.

The first concrete induced geometry-mixing mechanism is no longer an unspecified `O(chi^2)` nuisance. It has one exact coefficient pair `(a_2,a_4)` and one exact danger variable

\[
u = r(1-c_g).
\]

The resulting message is surprisingly strong:

1. pure-pole scalar-lane mixing is automatically safe;
2. mixed contact+pole scalar lanes become dangerous only at rather large quadrupole/scalar pole ratios;
3. on the exact Lambda_EM-refreshed Family-1 branch, actual twin/F1 failure requires `u` above about `6.83` / `10.28`, respectively;
4. therefore this first induced `l=0 \leftrightarrow l=2` mechanism looks much more like a controlled correction than like the actual source of 5PN failure.

So the next honest continuation is no longer “some mixing exists.” It is:

> derive from the actual moving-throat reduced operator whether the first nontrivial scalar/geometry lane is closer to the safe pure-pole end, the safe contact-dominated end, or the genuinely dangerous fast mixed contact+pole regime.
# 5PN continuation — Stages 268–271

## Goal

Continue from the Stage-267 corridor theorem without oversimplifying the scalar/geometry lane.

The actual next question was:

> what does the first *known* scalar/geometry lane from the moving-throat / Family-1 program look like **at the level of the Stage-78 contamination coefficients** `(a_2,a_4)`, and does it really look dangerous enough to threaten the grouped-`P2` twin-safe / Family-1-safe corridor?

The answer from this batch is much sharper than a generic “maybe.” The explicit Family-1 scalar lane is a **coupled two-pole monopole breathing channel**. For the Stage-78 contamination coefficients, that two-pole channel has an **exact contamination-equivalent one-pole diagnostic** with

\[
\lambda_{\rm coeff}=\frac{G_2}{G_4}\approx 28.1365336267,
\qquad
c_{\rm eff}=\frac{G_2^2}{G_0G_4}\approx 0.6072383642.
\]

So the known scalar lane is neither pure pole nor pure contact at the contamination level. But on the actual isotropic moving-throat branch the geometry lane is dynamically inert through \(O(\omega^4)\), so this scalar lane only re-enters after an explicit \(l=0\leftrightarrow l=2\) mixing source is turned on. The actual isotropic branch itself remains exactly safe and the active reduced bottleneck remains the outgoing quadrupole normalization defect \(N_Q\).

---

## Stage 268 — exact effective one-pole diagnostic for arbitrary scalar lanes

Start from the exact Stage-78 contamination pair
\[
a_2=\frac{\Omega_Q^2M_0^2G_2}{K_{\rm pole}G_0^2},
\qquad
a_4=\frac{\Omega_Q^4M_0^2(G_0G_4-G_2^2)}{K_{\rm pole}G_0^3}.
\]

These depend only on the first three scalar-lane moments \((G_0,G_2,G_4)\). That means there is an exact “moment-equivalent” one-pole diagnostic for *any* scalar lane:

\[
\Omega_{\rm eff}^2=\frac{G_2}{G_4},
\qquad
c_{\rm eff}=\frac{G_2^2}{G_0G_4},
\qquad
G_{{\rm pole},\rm eff}=\frac{G_2^2}{G_4},
\qquad
G_{{\rm contact},\rm eff}=G_0-\frac{G_2^2}{G_4}.
\]

Then
\[
a_2
=
\left(\frac{M_0^2}{K_{\rm pole}G_0}\right)
\left(\frac{\Omega_Q^2}{\Omega_{\rm eff}^2}\right)c_{\rm eff},
\]
\[
a_4
=
\left(\frac{M_0^2}{K_{\rm pole}G_0}\right)
\left(\frac{\Omega_Q^2}{\Omega_{\rm eff}^2}\right)^2
c_{\rm eff}(1-c_{\rm eff}),
\]
**exactly**.

So the Stage-266 corridor theorem applies to any scalar/geometry lane through the exact effective one-pole diagnostic \((\Omega_{\rm eff}^2,c_{\rm eff})\). This is not a fit; it is an identity for the contamination pair \((a_2,a_4)\).

### Exact two-pole decomposition

For a genuine contact-plus-two-pole scalar lane
\[
D_g = G_c+\frac{R_1}{1-s/L_1}+\frac{R_2}{1-s/L_2},
\]
the obstruction term splits exactly as
\[
G_0G_4-G_2^2
=
G_c\!\left(\frac{R_1}{L_1^2}+\frac{R_2}{L_2^2}\right)
+
R_1R_2\!\left(\frac1{L_1}-\frac1{L_2}\right)^2.
\]

So a two-pole lane can pick up “effective contact” both from the literal contact slot and from pole separation.

---

## Stage 269 — the actual Family-1 monopole breathing channel

The 2PN constructive appendix already fixed the scalar monopole breathing channel as
\[
K_{00}(s)
=
-\frac{757}{2520}
+\frac{R_-}{1-s/\lambda_-}
+\frac{R_+}{1-s/\lambda_+},
\]
with
\[
\lambda_- \approx 6.405572392138922,
\qquad
\lambda_+ \approx 254.444968136936126,
\]
\[
R_- \approx 0.002552474771738,
\qquad
R_+ \approx 0.386733239513976,
\]
and exact static value
\[
K_{00}(0)=\frac{4}{45}.
\]

So the actual constructive scalar lane is explicitly **two-pole**.

From this we get the first three low-frequency moments
\[
G_0=\frac{4}{45},
\qquad
G_2=\frac{R_-}{\lambda_-}+\frac{R_+}{\lambda_+}
\approx 1.91838640121275\times10^{-3},
\]
\[
G_4=\frac{R_-}{\lambda_-^2}+\frac{R_+}{\lambda_+^2}
\approx 6.81813341566730\times10^{-5}.
\]

The exact contamination-equivalent one-pole diagnostic is therefore
\[
\lambda_{\rm coeff}=\frac{G_2}{G_4}\approx 28.1365336267,
\qquad
c_{\rm eff}=\frac{G_2^2}{G_0G_4}\approx 0.6072383642.
\]

So for the Stage-78 contamination coefficients, the actual Family-1 two-pole breathing lane behaves exactly like a **mixed contact-plus-pole** scalar lane with pole fraction about \(0.61\).

### Important non-equivalence

The 2PN note also quotes an interval-accuracy Padé reduction
\[
\lambda_{\rm pade}\approx 202.9235163675.
\]

This is **not** the same thing as \(\lambda_{\rm coeff}\). The Padé pole is useful for fitting \(K_{00}(s)\) over a finite \(s\)-interval. The coefficient-equivalent pole is the exact one-pole reduction that preserves \((a_2,a_4)\). Conflating them would be an unnecessary oversimplification.

---

## Stage 270 — exact corridor comparison for the actual Family-1 breathing lane

Using the exact coefficient-equivalent pole fraction \(c_{\rm eff}\), the scalar-lane danger variable becomes

\[
u_{\rm eff}=r_{\rm eff}(1-c_{\rm eff}),
\qquad
r_{\rm eff}:=\frac{\Omega_Q^2}{\lambda_{\rm coeff}},
\]
with
\[
1-c_{\rm eff}\approx 0.3927616358.
\]

So the Stage-266 thresholds become explicit pole-ratio thresholds:

### Initial-drift thresholds
\[
r_{\rm eff}>\frac{2}{1-c_{\rm eff}}
\approx 5.09214703737
\quad\Longrightarrow\quad
\text{support demand grows initially},
\]

\[
r_{\rm eff}>\frac{4}{1-c_{\rm eff}}
\approx 10.1842940747
\quad\Longrightarrow\quad
\text{the twin margin shrinks initially}.
\]

### Actual failure thresholds
\[
r_{\rm eff}\ge
\frac{4+2\sqrt2}{1-c_{\rm eff}}
\approx 17.3856774766
\quad\Longrightarrow\quad
\text{actual twin failure can occur},
\]

and, using the exact Lambda\(_{\rm EM}\)-refreshed hard Family-1 ceiling \(c_\*=0.7116102605\ldots\),
\[
r_{\rm eff}\ge
\frac{8c_\*+4\sqrt{c_\*(4c_\*-1)}}{1-c_{\rm eff}}
\approx 26.1684980784
\quad\Longrightarrow\quad
\text{actual Family-1 failure can occur}.
\]

So the first known constructive scalar lane is not dangerous by composition alone. It becomes dangerous only if the grouped quadrupole pole is **much faster** than the contamination-equivalent breathing pole in the same reduced spectral variable.

---

## Stage 271 — actual branch verdict

Now combine this with the later moving-throat theorem ledger.

On the actual isotropic branch the geometry lane is dynamically inert through \(O(\omega^4)\) with respect to the grouped real `P2` carrier:
\[
\epsilon_2=\epsilon_4=0.
\]

So the conservative grouped-`P2` branch is exactly
\[
c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13
\]
on the actual isotropic branch, and the explicit Family-1 support/source side is automatic there.

That means:

1. the actual isotropic branch is **exactly safe** against scalar/geometry contamination;
2. the first known scalar-lane correction enters only after an explicit \(l=0\leftrightarrow l=2\) mixing source is turned on;
3. that first correction is controlled by the **Family-1 two-pole monopole breathing lane**;
4. and even then it only threatens the corridor if the physical pole ratio \(r_{\rm eff}\) is large enough to exceed the thresholds above.

So the scalar/geometry lane is not the natural source of 5PN failure on the actual isotropic branch.

---

## Best current reading after Stages 268–271

The next theorem gate is now much smaller.

### What is now settled

- Any scalar/geometry lane has an exact contamination-equivalent one-pole diagnostic \((\Omega_{\rm eff}^2,c_{\rm eff})\) for the pair \((a_2,a_4)\).
- The actual constructive scalar lane from the 2PN Family-1 module is a **two-pole monopole breathing channel**, not a literal one-pole lane.
- For contamination purposes, that actual scalar lane is exactly equivalent to
  \[
  \lambda_{\rm coeff}\approx 28.1365,
  \qquad
  c_{\rm eff}\approx 0.60724.
  \]
- On the actual isotropic moving-throat branch,
  \[
  \epsilon_2=\epsilon_4=0
  \]
  exactly, so the scalar lane does not contaminate the grouped `P2` carrier there at all.
- The explicit Family-1 support/source side is already automatic on that same actual isotropic branch.

### What remains genuinely open

The scalar/geometry lane only re-enters after an explicit \(l=0\leftrightarrow l=2\) mixing source is turned on. At that point the only remaining scalar-lane uncertainty is the physical pole ratio
\[
r_{\rm eff}=\frac{\Omega_Q^2}{\lambda_{\rm coeff}}.
\]

So the next honest continuation is now:

> determine, on the actual weakly anisotropic moving-throat reduced operator, whether the physical grouped-quadrupole pole is anywhere near the explicit danger window
> \[
> r_{\rm eff}\gtrsim 17.39 \quad\text{(twin failure)},
> \qquad
> r_{\rm eff}\gtrsim 26.17 \quad\text{(Family-1 failure)}.
> \]

Unless it is, the first known scalar/geometry lane is a controlled correction, not the active 5PN bottleneck.
# 5PN / Moving-Throat Continuation — Stage 268: Bulk-PDE Firewall and Boundary-Response Reading

## Purpose

This note records a clarification that matters for the rest of the 5PN program:

> we should **not** try to repair the 5PN gap by changing the parent bulk GNLS equation.

The question is whether that statement is merely intuitive or whether it is actually aligned with the present derivation chain.

The verdict is:

- **yes**, the current program supports a strong **bulk-PDE firewall** reading;
- but the strongest honest theorem-level wording is **not** “5PN divergence is purely a boundary condition” in an absolute sense;
- it is the more careful statement that **within the present reduced hierarchy**, the remaining 5PN / 2.5PN / 4PN gap has been compressed to the **moving-throat branch data**, i.e. the conservative grouped-`P2` carrier, the outgoing quadrupole normalization, and the weak-axisymmetric orbit/selector defects, **not** a retuning of the parent GNLS medium.

---

## 1. What is already frozen and therefore should not be changed

The parent `4+1` bulk theory is already fixed as a gauged GNLS matter sector plus localized Maxwell sector, with the stiff-polytropic EOS
\[
P(\rho)=K\rho^5,
\]
so the exponent `n=5` is already part of the carried exact parent setup. The exact parent fields remain
\[
\psi(\mathbf X,t),\qquad A_M,\qquad \Sigma=r-R(\Omega,w,t),
\]
and the compact moving-throat program explicitly treats this as the correct parent-theory starting point rather than a tunable late-stage ansatz.

At the lower PN levels, the hierarchy has already frozen the `n=5` calibration and the conservative 1PN–4PN carry-forward ledgers. In particular:

- lower-order optical / 1PN matching fixed the EOS exponent to `n=5`;
- the 4PN summary says the full local 4PN sector is already closed within the declared hierarchy;
- and the only remaining full-4PN gap is **the same passive/outgoing quadrupole normalization** already isolated by the 2.5PN program.

So if we were to modify the parent GNLS medium now, we would no longer be continuing the already-calibrated hierarchy. We would be starting a different hierarchy.

---

## 2. What the actual 5PN reduced gap now depends on

The 5PN notes already say the explicit Family-1 support/source side is **not** the active bottleneck anymore. The live question is whether the actual grouped-`P2` / geometry branch realizes the minimal isotropic conservative quadrupole module, and then whether the passive/outgoing `l=2` branch carries the canonical outgoing normalization.

The same notes also show that on the natural isotropic branch:

- the grouped conservative carrier is the `3/4 + 1/4` module,
  \[
  \widehat Y_Q^{\rm cons}(\omega)
  =
  \frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
  \]
  so
  \[
  c_{\rm pole}=\frac14;
  \]
- the scalar/geometry `l=0` lane is orthogonal to the grouped real `l=2` bundle at linear isotropic order, so
  \[
  \epsilon_2=\epsilon_4=0
  \]
  on that branch;
- and the first nonzero geometry contamination requires explicit `l=0 \leftrightarrow l=2` mixing and enters only at
  \[
  O(\chi^2).
  \]

So the live 5PN ambiguity is already localized to the **moving-throat response branch**, not to the bulk EOS or bulk GNLS law.

---

## 3. Packet-level statement of the firewall

Stages 199–201 compressed the reduced closure problem to two finite data packets:

### Packet A — grouped branch packet
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad
A\in\{20,21,22\},
\]
plus the source-map factor `mhat_0`.

### Packet B — orbit/selector packet
any equivalent form of
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

Everything in the reduced 5PN / 2.5PN / 4PN endgame is an exact compiler output of those packets.

That means the remaining task is:

- compute the actual packet values from the moving-throat branch;
- test isotropy / `c_{\rm pole}=1/4` / outgoing normalization / weak-axisymmetric orbit lock;
- and **do not** reopen the parent GNLS action unless there is a direct contradiction with the already frozen lower-order hierarchy.

This is the precise sense in which the parent bulk PDE is “safe.”

---

## 4. The right way to phrase the user’s proposed update

The proposed reading is very close to the current theorem picture, but it needs one wording adjustment.

### Good and useful
It is correct and helpful to say:

1. the lower-order program already calibrated the background 4D medium;
2. the remaining 5PN problem is in the throat / response / boundary realization;
3. the active moving pieces are the restoring pole fraction and the first `l=0 \leftrightarrow l=2` mixing channel, not a new bulk EOS fit.

### Needs tightening
The strongest honest program-level wording is **not**
> “the 5PN divergence is purely a boundary condition”

because the compact moving-throat ledger is explicit that the remaining gap is a **realization** gap: whether the actual completed PDE realizes the reduced outgoing, mouth/core, coherent-support, and invariant structures strongly enough. So the right wording is:

> **Within the present reduced hierarchy, the unresolved 5PN data enter through the moving-throat branch / boundary-response sector, not through a retuning of the parent GNLS bulk medium.**

That is narrower, cleaner, and actually matches the current notes.

---

## 5. What to do with the black-hole / sink intuition

The “microscopic particle = subsonic / resonant interior survives” versus
“black hole = supersonic sink / resonance shredded” picture is a **reasonable physical intuition** for future branch classes, but it is **not yet** one of the program’s frozen theorem inputs.

So it should be used as:

- a motivation for why different compact objects might correspond to different throat boundary classes;

not yet as:

- a theorem we have already derived from the moving-throat PDE.

For the current 5PN derivation, the safe formal statement is simply:

- keep the exact parent bulk PDE fixed;
- classify the remaining freedom in the throat response/boundary branch.

---

## 6. Practical consequence for the next derivation move

The next derivation should therefore keep the bulk firewall explicit and continue only on the boundary/branch side:

1. keep the parent GNLS + localized Maxwell action frozen;
2. keep `n=5` and the 1PN–4PN carry-forward constants frozen;
3. treat the grouped `P2` carrier and outgoing normalization as branch data to be extracted, not bulk data to be refit;
4. continue deriving the actual grouped-`P2` / weak-axisymmetric branch formulas under that firewall.

This is the clean continuation point.

---

## 7. Recommended replacement text

A tighter version of the proposed context update is:

> **Bulk-PDE firewall for the 5PN moving-throat program.**  
> We keep the parent `4+1` GNLS + localized Maxwell bulk theory fixed, including the already-calibrated `P=K\rho^5` medium. The lower-order 1PN–4PN hierarchy has already frozen that background. Within the present reduced hierarchy, the remaining 5PN / 2.5PN / 4PN gap is localized to the moving-throat branch data: the isotropic grouped-`P2` conservative carrier, the passive/outgoing `l=2` normalization, and the weak-axisymmetric orbit/selector defects. On the natural isotropic branch this means `c_pole=1/4`, while the first geometry contamination enters only through explicit `l=0 \leftrightarrow l=2` mixing at `O(\chi^2)`. So the correct next move is to keep the bulk PDE fixed and continue extracting the throat boundary/response branch, rather than altering the parent GNLS equation.

That is the version I would carry forward.

# 5PN continuation notes — Stages 272–275

These stages pick up from the Stage-271 verdict that the actual isotropic grouped-`P2` / geometry branch is already conservatively clean and that the explicit Family-1 support/source side is automatic there. The point of the continuation is therefore **not** to reopen the bulk GNLS or to keep probing the scalar/geometry lane. It is to collapse the remaining reduced theorem gap all the way down to the outgoing quadrupole normalization factor.

The main new result is that the whole reduced 2.5PN / 4PN endgame can now be written as one outgoing scalar defect

a) `chi_Q - 1`,

with a clean canonical branch `chi_Q = 1`, and with an exact three-parameter isotropic DtN deformation algebra that shows which deformations can actually move it.

---

## Stage 272 — actual isotropic outgoing reduction

Start from the actual isotropic grouped-`P2` one-pole conservative carrier

`Y_Q^cons(omega) = N_Q [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2) ]`.

Expanding at low frequency gives

`K0 = N_Q`,

`K2 = N_Q / (4 Omega_Q^2)`,

`K4 = N_Q / (4 Omega_Q^4)`.

So once the actual isotropic one-pole branch is accepted, every conservative low-frequency coefficient is carried by one scalar defect `N_Q`.

Adding the retarded outgoing normalization factor gives

`Y_Q^ret(omega) = N_Q [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5) ] + O(omega^6)`.

The odd coefficient is therefore

`Gamma5 = N_Q chi_Q sigma_Q^can / 4`.

With the canonical target normalization

`K0_target = 54 G c_s^5 / (5 a^5 c^5)`,

`sigma_Q^can = 4 a^5 / (27 c_s^5)`,

one gets the exact bridge

`K0_target sigma_Q^can / 4 = 2 G / (5 c^5)`.

So the full odd normalization condition is

`mhat_0^2 chi_Q N_Q = 1`.

On the natural source-map branch `mhat_0 -> 1`, the remaining reduced obstruction is purely outgoing:

`N_Q = 1 / chi_Q`.

So after Stage 272 the live reduced theorem gap is no longer a mixed support/source/outgoing problem. It is one outgoing normalization number.

---

## Stage 273 — canonical outgoing `l=2` DtN match

Using the exact compact outgoing spherical `l=2` DtN fingerprint

`Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + O(z^8)`,

and the exact normalization relation

`Yhat_2^out(z) = -3 / Lambda_2^out(z)`,

one gets

`Yhat_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243 + O(z^8)`.

Now match this to the grouped-`P2` one-pole-plus-contact ansatz in the dimensionless variable `z = omega a / c_s`:

`Yhat_Q^ret(z) = 3/4 + (1/4)/(1 - alpha z^2 - i chi_Q B z^5) + O(z^6)`.

Exact coefficient matching gives

`alpha = 4/9`,

`B chi_Q = 4/27`.

Since

`alpha = c_s^2 / (a^2 Omega_Q^2)`,

the canonical pole scale is fixed to

`Omega_Q = 3 c_s / (2 a)`.

Then

`sigma_Q^can = 9 / (8 Omega_Q^5) = 4 a^5 / (27 c_s^5)`,

so the dimensionless odd coefficient is exactly `B = 4/27`, and therefore

`chi_Q = 1`.

This is the cleanest current proof that the canonical compact passive/outgoing `l=2` branch is the canonical reduced GR branch.

---

## Stage 274 — isotropic DtN deformation algebra

Now deform the isotropic outgoing DtN branch by

`Lambda_2^def(z)`
`= S Lambda_2^out(beta z) + Sigma0 + Sigma2 z^2 + Sigma4 z^4 + i Sigma5 z^5 + O(z^6)`.

Write again

`Yhat_2^def(z) = -3 / Lambda_2^def(z)`.

If the **normalized even moments** are required to stay canonical, i.e.

`y2 / y0 = 1/9`,

`y4 / y0 = 4/81`,

then the even deformation coefficients are fixed exactly:

`Sigma2 = - S beta^2 / 3 + S / 3 - Sigma0 / 9`,

`Sigma4 = - S beta^4 / 9 + S / 9 - Sigma0 / 27`.

After that elimination, the conservative and outgoing scalar defects are

`N_Q = y0 = 3 / (3 S - Sigma0)`,

`chi_Q = 27 y5 / y0 = 3 (S beta^5 + 9 Sigma5) / (3 S - Sigma0)`.

So once the even moments are pinned, only three isotropic DtN deformations can move the outgoing theorem:

- the scale deformation `beta`,
- the static DtN shift `Sigma0`,
- the odd `z^5` DtN shift `Sigma5`.

Everything else is already absorbed into maintaining the canonical even branch.

---

## Stage 275 — outgoing-defect linearization and exact no-shift family

The exact no-shift condition `chi_Q = 1` can be solved directly:

`Sigma0 = 3 S (1 - beta^5) - 27 Sigma5`.

So the canonical outgoing branch is **not unique**. Once the even moments are preserved, there is an exact two-parameter isotropic deformation family that still leaves `chi_Q = 1`.

Linearizing around the canonical compact branch,

`S = 1 + delta_S`,

`beta = 1 + delta_beta`,

`Sigma0 = delta_Sigma0`,

`Sigma5 = delta_Sigma5`,

one gets the first-order outgoing defect

`Delta_Q := chi_Q - 1`

`= 5 delta_beta + delta_Sigma0 / 3 + 9 delta_Sigma5 + O(2)`.

The overall amplitude deformation `delta_S` cancels out at first order. So at linear order the outgoing theorem is controlled only by:

1. the pole rescaling `beta`,
2. the static DtN shift `Sigma0`,
3. the odd DtN shift `Sigma5`.

A useful special slice is the already conservative-normalized branch `N_Q = 1`, which forces

`Sigma0 = 3 (S - 1)`.

On that slice the remaining outgoing defect reduces to

`chi_Q = S beta^5 + 9 Sigma5`.

So once the conservative isotropic branch is clean, the remaining reduced obstruction is entirely retarded.

---

## Best current reading after Stage 275

The reduced theorem gap is now smaller than it was at Stage 271.

1. The actual isotropic grouped-`P2` conservative branch carries one scalar defect `N_Q`.
2. On the natural source map, the remaining reduced odd obstruction is purely outgoing and is measured by `chi_Q`.
3. The canonical compact passive/outgoing `l=2` DtN branch gives `chi_Q = 1` exactly.
4. The most general isotropic DtN deformation that preserves the canonical even moments moves `chi_Q` only through `beta`, `Sigma0`, and `Sigma5`.
5. Near the canonical branch, the exact first-order outgoing defect is

   `Delta_Q = 5 delta_beta + delta_Sigma0 / 3 + 9 delta_Sigma5 + O(2)`.

So the next honest theorem gate is no longer on the scalar/geometry or explicit support/source side. It is to determine whether the actual passive/outgoing moving-throat DtN branch is the canonical one (`chi_Q = 1`) or, if not, which concrete isotropic DtN deformation data `beta`, `Sigma0`, `Sigma5` shift it away.
# 5PN continuation — Stages 276–279

## Goal

Continue directly from the Stage-272–275 outgoing-normalization collapse without reopening the bulk PDE.
The next clean move was to connect the **actual selected moving-throat outgoing branch** to the isotropic DtN deformation variables

\[
(\beta,\Sigma_0,\Sigma_5),
\]

so that `chi_Q` stops being an abstract leftover scalar and becomes a direct continuum-kernel observable.

The key point that emerged immediately is that this connection is **not unique**. The exact isotropic DtN deformation algebra only fixes a two-parameter **gauge family** of triples `(beta, Sigma_0, Sigma_5)` once the observable pair `(N_Q, chi_Q)` is specified. So the right result is not “the” unique triple, but an exact family together with physically useful gauge choices.

---

## Stage 276 — exact selected-branch DtN gauge family

Start from the exact isotropic deformation laws already fixed at Stage 274:

\[
N_Q = \frac{3}{3S-\Sigma_0},
\qquad
\chi_Q = \frac{3\bigl(S\beta^5+9\Sigma_5\bigr)}{3S-\Sigma_0}.
\]

Solving exactly for the deformation data gives

\[
\boxed{
\Sigma_0 = 3S - \frac{3}{N_Q},
\qquad
\Sigma_5 = \frac{\chi_Q}{9N_Q}-\frac{S\beta^5}{9}.
}
\]

So the observable pair `(N_Q, chi_Q)` fixes only those two combinations. The actual selected outgoing branch therefore determines a **two-parameter DtN gauge family**, labeled by `(S, beta)`.

Three exact gauges are especially useful.

### A. Core gauge
Set
\[
S=1,
\qquad
\beta=1.
\]
Then
\[
\boxed{
\Sigma_0 = 3\Bigl(1-\frac{1}{N_Q}\Bigr),
\qquad
\Sigma_5 = \frac{\chi_Q-N_Q}{9N_Q}.
}
\]
This packages the whole defect into a static core shift plus an odd core outlet.

### B. Scale gauge
Set
\[
\Sigma_0=0,
\qquad
\beta=1.
\]
Then
\[
\boxed{
S = \frac{1}{N_Q},
\qquad
\Sigma_5 = \frac{\chi_Q-1}{9N_Q}.
}
\]
This packages the same selected branch into pure mouth normalization plus an odd core outlet.

### C. Argument gauge
Set
\[
\Sigma_0=0,
\qquad
\Sigma_5=0,
\]
with the natural positive branch. Then
\[
\boxed{
S = \frac{1}{N_Q},
\qquad
\beta = \chi_Q^{1/5}.
}
\]
So the same selected branch can also be viewed as a pure effective-argument deformation.

---

## Stage 277 — exact compensated Robin–mixed outlet dictionary

Now take the explicit isotropic moving-throat outlet class

\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}.
\]

On the exact compensated canonical-even branch

\[
\boxed{
\rho_R = 4\sigma_W,
\qquad
\kappa_W = \frac13,
}
\]

one gets the exact outgoing normalization factor

\[
\boxed{
\chi_Q^{\rm hyb} = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
}
\]

This explicit outlet admits two exact DtN gauge embeddings.

### Core gauge embedding
Choose
\[
S=1,
\qquad
\beta=1,
\qquad
\Sigma_0 = 3\sigma_W,
\qquad
\Sigma_5 = -\sigma_W\gamma_W,
\]
with the required even slots
\[
\Sigma_2 = -\frac{\sigma_W}{3},
\qquad
\Sigma_4 = -\frac{\sigma_W}{9}.
\]
This reproduces the compensated hybrid outlet exactly.

### Scale gauge embedding
Choose
\[
S = 1-\sigma_W,
\qquad
\beta = 1,
\qquad
\Sigma_0 = 0,
\qquad
\Sigma_5 = \sigma_W\Bigl(\frac19-\gamma_W\Bigr),
\]
with
\[
\Sigma_2 = \Sigma_4 = 0.
\]
This also reproduces the same outlet exactly.

So the compensated Robin–mixed outlet is the first explicit moving-throat outgoing branch that can be embedded into the Stage-276 DtN gauge family in more than one natural way.

A very useful special case is the canonical compensated outgoing point
\[
\gamma_W = \frac19.
\]
Then in the **scale gauge**
\[
\Sigma_5 = 0,
\qquad
\Lambda_2^{\rm hyb}(z)=(1-\sigma_W)\Lambda_2^{\rm out}(z),
\]
so the whole compensated outlet reduces to the robust pure-scale class and therefore
\[
\boxed{\chi_Q=1.}
\]

---

## Stage 278 — `chi_Q` as a direct continuum-kernel observable

The exact compensated branch makes the outgoing factor explicit:

\[
\boxed{
\chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
\qquad
\Delta_Q:=\chi_Q-1 = \frac{\sigma_W(1-9\gamma_W)}{1-\sigma_W}.
}
\]

So the last reduced outgoing defect is already a direct continuum observable of the selected branch.
It depends only on

- the static mixed loading `sigma_W`, and
- the odd mixed outlet coefficient `gamma_W`.

The relation can be inverted exactly:

\[
\boxed{
\gamma_W = \frac{1-(1-\sigma_W)\chi_Q}{9\sigma_W}.
}
\]

If one additionally imposes the natural source-map odd normalization condition from Stages 83–84, then the **required** conservative normalization becomes

\[
\boxed{
N_Q^{\rm req} = \frac{1}{\chi_Q} = \frac{1-\sigma_W}{1-9\sigma_W\gamma_W}.
}
\]

So the actual selected moving-throat outgoing branch no longer hides behind an abstract `chi_Q`. Once `(sigma_W, gamma_W)` are computed on the physical branch, both `chi_Q` and the required `N_Q` are known immediately.

The exact DtN-gauge observables are then:

- **core gauge**:
  \[
  \Sigma_0 = 3\sigma_W,
  \qquad
  \Sigma_5 = -\sigma_W\gamma_W,
  \]
- **scale gauge**:
  \[
  S = 1-\sigma_W,
  \qquad
  \Sigma_5 = \sigma_W\Bigl(\frac19-\gamma_W\Bigr)=\frac{(\chi_Q-1)(1-\sigma_W)}{9}.
  \]

So in the scale gauge, the entire outgoing defect is seen as one odd core slot riding on top of a pure mouth renormalization.

---

## Stage 279 — linear Family-1 canonical-even projection into `(beta, Sigma_0, Sigma_5)`

The exact linear defect formulas on the compensated hybrid outlet are

\[
\delta E_2 = \frac{\delta\mathcal C - 9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\qquad
\delta E_4 = \frac{5\delta\mathcal C - 72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]

\[
\Delta_Q = \frac{\delta\mathcal C - 27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)}.
\]

Imposing exact first-order preservation of the canonical even `l=2` fingerprint,

\[
\delta E_2 = 0,
\qquad
\delta E_4 = 0,
\]

forces uniquely

\[
\boxed{
\delta\mathcal C = 0,
\qquad
\delta\kappa_W = 0.
}
\]

Using the explicit Family-1 transport law from Stage 141–142, that further implies

\[
\boxed{
\delta\mathfrak g = 0.
}
\]

So on the canonical-even compensated branch the only surviving linear mouth/core defect is the odd mixed-channel renormalization `delta gamma_W`, and the remaining outgoing defect is

\[
\boxed{
\Delta_Q = -\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

Now compare with the Stage-92 linear DtN law

\[
\Delta_Q = 5b + \frac{a_0}{3} + 9a_5.
\]

A natural compensated **core gauge** choice is therefore

\[
\boxed{
\delta\beta = 0,
\qquad
\delta\Sigma_0 = 0,
\qquad
\delta\Sigma_5 = -\frac{\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

Then indeed

\[
\boxed{
\Delta_Q = 9\,
\delta\Sigma_5.
}
\]

So after enforcing canonical-even preservation, the whole linearized Family-1 mouth/core defect projects into the isotropic DtN triple as a **pure odd core outlet**. That is the sharpest reduced bridge so far.

---

## Best current reading after Stage 279

The selected moving-throat outgoing branch is now connected to the isotropic DtN deformation variables in the strongest honest way available without solving the full PDE.

1. The exact isotropic deformation algebra determines only a **gauge family** of `(beta, Sigma_0, Sigma_5)` for a given `(N_Q, chi_Q)`.
2. The compensated Robin–mixed outlet gives the first explicit moving-throat branch inside that family.
3. On that explicit branch,
   
   \[
   \chi_Q = \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}
   \]
   
   is already a direct continuum-kernel observable.
4. On the canonical-even Family-1 branch, every first-order mouth/core defect except the odd mixed-channel renormalization is forced away, so the DtN projection becomes
   
   \[
   \delta\beta=0,
   \qquad
   \delta\Sigma_0=0,
   \qquad
   \delta\Sigma_5=-\frac{\sigma_*}{1-\sigma_*}\delta\gamma_W.
   \]

So the remaining reduced outgoing theorem gap is now as narrow as it can be in this language:

> compute the physical odd mixed-channel renormalization `gamma_W` of the selected passive/outgoing branch, because that one continuum observable already determines the surviving isotropic DtN defect.
# 5PN continuation — Stages 280–283: selected loading/product variables into compensated outlet observables

This session keeps the bulk-PDE firewall intact. Nothing here retunes the parent `4+1` GNLS/Maxwell medium. The work stays entirely on the selected branch / outlet side.

## Stage 280 — compensated hybrid outlet reproduces the minimal isotropic conservative module exactly

Start from the compensated Robin–mixed isotropic outlet branch

\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+4\sigma_W
-\frac{\sigma_W}{1-z^2/3-i\gamma_W z^5}
+O(z^6),
\]

with the compensated canonical-even conditions

\[
\rho_R = 4\sigma_W,
\qquad
\kappa_W = \frac13.
\]

Normalizing by the static outlet value gives

\[
\widehat Y_2^{\rm hyb}(z)
=
1+\frac{z^2}{9}+\frac{4z^4}{81}
+i\frac{1-9\sigma_W\gamma_W}{27(1-\sigma_W)}z^5+O(z^6).
\]

So the whole compensated hybrid outlet is exactly equivalent, through `O(z^5)`, to the minimal isotropic contact-plus-pole module

\[
\widehat Y_2^{\rm hyb}(z)
=
\frac34
+
\frac14\frac{1}{1-\frac49 z^2-i\frac{4}{27}\chi_Q z^5}
+O(z^6),
\qquad
\chi_Q=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]

Therefore the conservative selected-branch weights are frozen to

\[
c_0=\frac34,
\qquad
c_1=\frac14,
\qquad
\rho_\alpha=\frac{1}{c_0}=\frac43.
\]

Using the exact blocked-demand map

\[
\zeta_{\rm req}
=
\frac{\rho_\alpha-1}{1-\epsilon_{\rm blk}(2-\rho_\alpha)},
\]

this gives

\[
\zeta_{\rm req}(\epsilon_{\rm blk})=\frac{1}{3-2\epsilon_{\rm blk}}.
\]

But the selected loading/product ratio itself stays fixed:

\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
Q(\zeta_{\rm req},\epsilon_{\rm blk})
=
\frac43.
\]

So `epsilon_blk` changes the support-demand variable `zeta_req`, but on the minimal isotropic conservative branch it does **not** change the selected product ratio `Pi_tr/C_mix`.

## Stage 281 — exact outlet-observable transport on the compensated canonical-even branch

The compensated hybrid outgoing factor is

\[
\chi_Q(\sigma_W,\gamma_W)=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]

Its full first differential is

\[
d\chi_Q
=
\frac{1-9\gamma_W}{(1-\sigma_W)^2}\,d\sigma_W
-\frac{9\sigma_W}{1-\sigma_W}\,d\gamma_W.
\]

At the canonical outgoing point

\[
\gamma_W=\frac19,
\]

the `d sigma_W` term vanishes identically, so

\[
d\chi_Q\Big|_{\gamma_W=1/9}
=
-\frac{9\sigma_*}{1-\sigma_*}\,d\gamma_W.
\]

This matches the earlier outlet-projection formulas. Writing the first-order compensated-branch defects as

\[
\delta E_2 = \frac{\delta\mathcal C - 9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\]
\[
\delta E_4 = \frac{5\delta\mathcal C - 72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]
\[
\Delta_Q = \frac{\delta\mathcal C - 27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)},
\]

the canonical-even gate `delta E_2 = delta E_4 = 0` forces

\[
\delta\mathcal C = 0,
\qquad
\delta\kappa_W = 0,
\]

so the only surviving isotropic outlet defect is

\[
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}
\,\delta\gamma_W.
\]

Thus `sigma_W` survives only as a static amplification factor on the canonical-even tangent; the actual odd mismatch is carried purely by `delta gamma_W`.

## Stage 282 — exact loading/product to outlet-slippage bridge

The Family-1 loading side and the outlet side meet through the exact first-order identity

\[
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]

Comparing this with the canonical-even outlet formula

\[
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}
\,\delta\gamma_W,
\]

gives the exact tangent transport law

\[
\delta\gamma_W
=
\frac{1}{9}\,
\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]

On the weak-axisymmetric grouped branch, the same defect is also

\[
\Delta_Q = \frac{P_1}{P_0},
\]

because the compensated-branch slope is

\[
\gamma_1
=
-\frac{1-\sigma_*}{9\sigma_*}
\frac{P_1}{P_0}.
\]

So the remaining outgoing theorem gap is one scalar written in three exactly equivalent ways:

\[
\Delta_Q,
\qquad
\frac{P_1}{P_0},
\qquad
\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]

## Stage 283 — one-language end-to-end theorem ledger

After the bridge is made explicit, the selected-branch theorem data split cleanly into two packets.

### Conservative selected loading/product side

\[
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}(\epsilon_{\rm blk}) = \frac{1}{3-2\epsilon_{\rm blk}},
\qquad
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \frac43.
\]

### Compensated outlet side

\[
\chi_Q(\sigma_W,\gamma_W)=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
\]

and on the nontrivial compensated branch the canonical odd condition is simply

\[
\chi_Q=1
\iff
\gamma_W=\frac19.
\]

### Tangent bridge between them

\[
\delta\gamma_W
=
\frac{1}{9}\,
\Xi_{\rm slip}\,\delta\Pi_{\rm tan},
\qquad
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}
\Xi_{\rm slip}\,\delta\Pi_{\rm tan},
\qquad
\Delta_Q = \frac{P_1}{P_0}.
\]

So the final outgoing theorem gap is now written in one language end to end:

1. selected conservative product side: `Pi_tr/C_mix = 4/3`,
2. tangential branch slippage: `delta Pi_tan`,
3. outlet odd slippage: `delta gamma_W = Xi_slip delta Pi_tan / 9`,
4. final outgoing defect: `Delta_Q = chi_Q - 1 = P1/P0`.

That is the sharpest compact bridge we have so far between the selected branch, the explicit outlet observables, and the remaining 2.5PN / 4PN outgoing normalization defect.
# 5PN continuation — Stages 284–288: co-evolving Family-1 fixed-point bridge, on-family rigidity, and the off-family scalar

This session takes the end-to-end bridge from Stages 280–283

a) from the selected conservative product side,

b) through the compensated hybrid outlet,

c) down into the explicit co-evolving Family-1 fixed-point kernel.

The main result is that the final outgoing obstruction is now much sharper than
just “compute `deltaPi_tan`.”

On the explicit co-evolving Family-1 branch,

- `deltaPi_tan` is already explicit,
- the bare mixed-port odd slippage reduces to one D/N similarity-slippage scalar,
- that scalar vanishes identically on the exact lower parent compensation family,
- and any first-order failure is therefore carried by one off-family normal coordinate `delta_perp`.

So the next genuine PDE-facing question is no longer to guess an outlet defect. It is to determine whether the actual branch stays on the exact lower compensation family, and if not, what its single off-family scalar `delta_perp` is.

## Stage 284 — exact co-evolving Family-1 fixed-point transport

At the renormalized canonical point, the mouth/core load transport is

\[
\delta M_s = \delta\Sigma_0,
\]
\[
\delta M_q
=
-\frac14\,\delta\Sigma_0
+
\frac{\Sigma_0^{\rm can}}{\sqrt{1+\mathfrak r_{F1}^2}}\,\delta\mathfrak g,
\]
\[
\delta\Pi
=
\left(1-\frac14\mathcal S_{\rm can}\right)\delta\Sigma_0
-
\frac{\Sigma_0^{\rm can}}{4}\,\delta\mathcal S
+
\frac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_{F1}^2}}\,\delta\mathfrak g.
\]

On the canonical-even tangent `delta g = 0`, this collapses to

\[
\boxed{
\delta\Pi_{\rm tan}
=
\left(1-\frac14\mathcal S_{\rm can}\right)\delta\Sigma_0
-
\frac{\Sigma_0^{\rm can}}{4}\,\delta\mathcal S.
}
\]

Numerically,

\[
\boxed{
\delta\Pi_{\rm tan}
\approx
0.832409471081634\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S.
}
\]

Using

\[
\Sigma_0 = \frac{20}{9}\widehat T_m^2,
\qquad
\delta\Sigma_0 = \frac{40}{9}\widehat T_{m,\rm can}\,\delta\widehat T_m,
\]

this is equivalently

\[
\boxed{
\delta\Pi_{\rm tan}
\approx
5.35223887169622\,\delta\widehat T_m
-
1.16275838754222\,\delta\mathcal S.
}
\]

So the tangential defect transport is already explicit on the fixed-point branch.

## Stage 285 — bare mixed-port slippage and D/N similarity reduction

The compensated concrete core obeys

\[
\kappa_W = \frac{\kappa_0}{1+r_c},
\qquad
\gamma_W = \frac{\gamma_0}{1+r_c}.
\]

Linearizing on the compensated branch gives the exact identity

\[
\boxed{
\delta\gamma_W - \frac13\delta\kappa_W
=
\frac{1}{1+r_{c,*}}
\left(
\delta\gamma_0 - \frac13\delta\kappa_0
\right).
}
\]

Under the canonical-even gate `delta kappa_W = 0`, the odd outlet defect is therefore carried by the bare mixed-port slippage

\[
\delta\mathfrak B_W := \delta\gamma_0 - \frac13\delta\kappa_0,
\qquad
\boxed{
\delta\gamma_W = \frac{\delta\mathfrak B_W}{1+r_{c,*}}.
}
\]

Now use the D/N realization

\[
\kappa_0 = \frac{4L_W^2}{\pi^2 a^2},
\qquad
\gamma_0 = \frac{1+r_c}{9}.
\]

Then

\[
\delta\mathfrak B_W
=
\frac{1+r_{c,*}}{9}
\left[
\delta\ln\gamma_0 - 2\,\delta\ln\left(\frac{L_W}{a}\right)
\right].
\]

Define the D/N similarity-slippage scalar

\[
\Xi_{\rm slip} := \Xi_\gamma - 2\Xi_L.
\]

If `delta mathfrak B_W = Upsilon_Pi deltaPi_tan`, then

\[
\Upsilon_\Pi = \frac{1+r_{c,*}}{9}\,\Xi_{\rm slip},
\]

and the first-order outgoing defect collapses to

\[
\boxed{
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
}
\]

So the last effective odd defect is not a generic DtN susceptibility. It is one D/N similarity-strain mismatch.

## Stage 286 — parent compensation-surface rigidity and automatic similarity preservation

On the exact parent compensation family,

\[
1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2,
\qquad
\frac{L_W}{a} = \frac{\pi}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
\qquad
\gamma_0 = \frac{1+\mathfrak r^2}{9}.
\]

Therefore

\[
\delta\ln\gamma_0
-
2\,\delta\ln\left(\frac{L_W}{a}\right)=0
\qquad\Longrightarrow\qquad
\boxed{\Xi_{\rm slip}=0}
\]

identically along the family.

On the lower branch,

\[
\mathfrak g_-(\mathfrak r)=\mathfrak r-\frac12\sqrt{1+\mathfrak r^2},
\qquad
\mathfrak g_-'(\mathfrak r)
=
1-\frac{\mathfrak r}{2\sqrt{1+\mathfrak r^2}}>0.
\]

So the carried canonical-even condition `delta g = 0` forces

\[
\boxed{\delta\mathfrak r=0.}
\]

Hence all first-order D/N similarity defects freeze,

\[
\delta\ln\gamma_0=0,
\qquad
\delta\ln(L_W/a)=0,
\qquad
\delta\mathfrak B_W=0,
\qquad
\delta\gamma_W=0,
\]

and therefore

\[
\boxed{\Delta_Q=0,
\qquad
N_Q-1=0}
\]

at first order.

This is the strongest positive result of the batch:

> if the actual co-evolving moving-throat branch stays on the exact lower parent compensation family, the first-order reduced 2.5PN / 4PN outgoing defect disappears automatically.

## Stage 287 — off-family normal coordinate

To measure genuine departure from that family, define the exact parent compensation defect

\[
\mathcal F(\mathfrak g,\mathfrak r) := 1+\mathfrak r^2 - 4(\mathfrak g-\mathfrak r)^2,
\]

and the off-family normal coordinate

\[
\boxed{
\delta_\perp
:=
\delta\mathfrak g
-
\mathfrak g_-'(\mathfrak r_*)\,\delta\mathfrak r.
}
\]

Then on the lower branch,

\[
\boxed{
\delta\mathcal F = 4\sqrt{1+\mathfrak r_*^2}\,\delta_\perp,
}
\qquad
\boxed{
\delta R_q = -\frac{\delta_\perp}{\sqrt{1+\mathfrak r_*^2}}.
}
\]

The exact microscopic parent-variable formula is

\[
\boxed{
\delta_\perp
=
\mathfrak g_*\,\delta\ln\left(\frac{g_q K_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}
\,\delta\ln\left(\frac{K_s K_q}{\lambda^2}\right).
}
\]

So the whole first-order off-family defect is carried by one scalar only.

## Stage 288 — explicit microscopic log channels and lower-branch cancellation

Using the carried explicit throat-core formulas, the two imbalance channels become

\[
\delta\ln\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
\]

\[
\delta\ln\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
\]

So `delta_perp` is already an explicit linear combination of

\[
(\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W).
\]

Then inserting the exact lower compensated branch drift laws,

\[
\delta\ln L_W = \delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\left(
\delta\ln\mathcal Z_q - \delta\ln\rho_w + 3\delta\ln c_{s,w} + 2\delta\ln c_s - 5\delta\ln a
\right),
\]
\[
\delta\ln\mathcal T_m
=
\frac12\left(
\delta\ln\mathcal Z_q - \delta\ln\rho_w + 3\delta\ln c_{s,w} - 2\delta\ln c_s - 3\delta\ln a
\right),
\]

both imbalance channels vanish identically, hence

\[
\boxed{\delta_\perp=0}
\]

on the exact lower compensated branch.

So the theorem picture is now extremely tight:

- **on-family:** first-order outgoing defect vanishes automatically,
- **off-family:** the whole first-order defect is carried by the single scalar `delta_perp`.

## Net result after Stages 284–288

The remaining outgoing theorem gap is no longer the tangential load itself.
It is now the much smaller question:

> does the actual moving-throat branch stay on the exact lower parent compensation family to first order, or does it develop a nonzero off-family scalar `delta_perp`?

If it stays on-family, the first-order reduced 2.5PN / 4PN outgoing obstruction is gone. If not, the defect ledger is already explicit and one-dimensional.
# 5PN continuation — Stages 289–295: isotropic bundle tangency, off-bundle scalar slippage, and collapse to the physical weak-axisymmetric slopes

This batch continues directly from the Stage 284–288 fixed-point bridge.
At that point the first-order selected-branch outgoing defect had already been reduced to one off-family scalar
\[
\delta_\perp,
\]
and the live question was whether the actual moving-throat branch keeps that scalar at zero.

The new result is that the answer splits very sharply.

1. **Arbitrary first-order isotropic bundle drift cannot generate** `delta_perp`.
   It is tangent to the exact lower compensated parent family.
2. The first genuine first-order off-family correction is one weighted scalar slippage combination
   \[
   \varepsilon_\perp,
   \]
   built only from three exact transport-law failures.
3. A pure weak grouped real `P2` anisotropy cannot source that scalar at linear order.
4. Therefore the remaining linear grouped bottleneck is only in the direct outlet coefficients
   \[
   \delta\kappa_W,
   \qquad
   \delta\gamma_W,
   \]
   and the whole weak-axisymmetric outlet problem then collapses to the two physical slopes
   \[
   u_2^{(1)},
   \qquad
   P_1/P_0.
   \]

So the next theorem gate is no longer a broad “compute all first-order corrections” problem.
It is much smaller: derive the actual moving-throat values of the physical grouped slopes
`u_2^(1)` and `P_1/P_0`.

## Stage 289 — exact bundle transport compiler

The last four irreducible lower-branch microscopic drifts are exact algebraic images of
\[
(\Theta_w,
 K_s,
 K_q,
 P_0),
\qquad
P_0 = N_0/D_0.
\]

The exact inversion is
\[
\delta\ln\rho_w=
\frac12\,\delta\ln\Theta_w,
\]
\[
\delta\ln a=
\frac12\,\delta\ln K_s-
\frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln c_s=
\frac12\,\delta\ln K_s-
\frac14\,\delta\ln\Theta_w+
\frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q=
\delta\ln K_q-
\frac25\,\delta\ln P_0.
\]

All remaining mouth/background drifts are then co-transported:
\[
\delta\ln c_{s,w}=\delta\ln\Theta_w,
\qquad
\delta\ln\ell=-\delta\ln\Theta_w,
\qquad
\delta\ln L_W=
\frac12\,\delta\ln K_s-
\frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln v_{w0}=
-
\frac34\,\delta\ln K_s+
\frac12\,\delta\ln K_q+
\frac{13}{8}\,\delta\ln\Theta_w,
\]
\[
\delta\ln\mathcal T_m=
-
\frac54\,\delta\ln K_s+
\frac12\,\delta\ln K_q+
\frac{15}{8}\,\delta\ln\Theta_w-
\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_s=
-
\frac14\,\delta\ln K_s+
\frac12\,\delta\ln K_q+
\frac38\,\delta\ln\Theta_w-
\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_q=
-
\frac34\,\delta\ln K_s+
\delta\ln K_q+
\frac38\,\delta\ln\Theta_w-
\frac25\,\delta\ln P_0,
\]
\[
\delta\ln\lambda=
\frac12\bigl(\delta\ln K_s+
\delta\ln K_q\bigr).
\]

So the actual isotropic branch is no longer missing “many independent drifts.”
It is algebraically controlled by the four bundle observables above.

## Stage 290 — exact tangent-compensation theorem

Substituting the Stage-289 compiler into the exact parent invariants gives
\[
\delta\ln r_c = 2\,\delta\ln\lambda-
\delta\ln K_s-
\delta\ln K_q = 0,
\]
\[
\delta\ln\mathfrak r = \delta\ln\lambda-
\frac12(\delta\ln K_s+
\delta\ln K_q)=0,
\]
\[
\delta\ln\mathfrak g = \delta\ln g_q+
\frac12\,\delta\ln K_s-
\delta\ln g_s-
\frac12\,\delta\ln K_q=0.
\]

Hence both Stage-147 logarithmic imbalance channels vanish exactly:
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)=0,
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0,
\]
and therefore
\[
\boxed{\delta_\perp=0.}
\]

So arbitrary first-order isotropic bundle drift is tangent to the exact lower compensated parent family.
The mouth-bias law collapses to its tangential piece,
\[
\delta\Pi=\delta\Pi_{\rm tan},
\]
and the first genuine first-order danger must come from **off-bundle** structure.

## Stage 291 — exact off-bundle slippage bridge

The first off-bundle correction is carried by three exact slippages
\[
\varepsilon_L:=\delta\ln L_W-\delta\ln a,
\]
\[
\varepsilon_v:=
\delta\ln v_{w0}-
\left[
\frac12\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+
\frac32\delta\ln c_{s,w}+
\delta\ln c_s-
\frac52\delta\ln a
\right],
\]
\[
\varepsilon_T:=
\delta\ln\mathcal T_m-
\left[
\frac12\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+
\frac32\delta\ln c_{s,w}-
\delta\ln c_s-
\frac32\delta\ln a
\right].
\]

Substituting them into the Stage-147 normal-coordinate formula gives the exact scalar collapse
\[
\boxed{\delta_\perp=-\varepsilon_\perp,}
\]
with
\[
\boxed{
\varepsilon_\perp=
\mathfrak g_*\varepsilon_T+
\left(\mathfrak g_*+
\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_v+
\left(2\mathfrak g_*+
\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_L.
}
\]

Numerically on the renormalized Family-1 point,
\[
\delta_\perp
\approx
-0.758035078944663\,\varepsilon_T
-1.00314310113848\,\varepsilon_v
-1.88373219118005\,\varepsilon_L.
\]

The same scalar controls the mouth-bias and conservative-even outlet defects:
\[
\delta\Pi
=
\delta\Pi_{\rm tan}
-
\frac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_*^2}}\,\varepsilon_\perp,
\]
so numerically
\[
\delta\Pi
\approx
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
-
1.52843317823248\,\varepsilon_\perp,
\]
that is
\[
\delta\Pi
\approx
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
-
2.87915877990416\,\varepsilon_L
-
1.53323719829507\,\varepsilon_v
-
1.15860596492310\,\varepsilon_T.
\]

And the direct outlet ledger is
\[
\delta\mathcal C
=
-
\frac{16\sigma_*}{\sqrt{1+\mathfrak r_*^2}}\,\varepsilon_\perp,
\]
\[
\delta E_2
=
\frac{\sigma_*}{27(1-\sigma_*)}
\left[
-
\frac{16}{\sqrt{1+\mathfrak r_*^2}}\varepsilon_\perp
-
9\,\delta\kappa_W
\right],
\]
\[
\delta E_4
=
\frac{\sigma_*}{243(1-\sigma_*)}
\left[
-
\frac{80}{\sqrt{1+\mathfrak r_*^2}}\varepsilon_\perp
-
72\,\delta\kappa_W
\right],
\]
\[
\Delta_Q
=
\frac{\sigma_*}{3(1-\sigma_*)}
\left[
-
\frac{16}{\sqrt{1+\mathfrak r_*^2}}\varepsilon_\perp
-
27\,\delta\gamma_W
\right].
\]

So the first-order off-family defect is now only one weighted scalar, and the remaining odd normalization defect still needs `delta gamma_W` unless `eps_perp` and `delta kappa_W` are both killed.

## Stage 292 — no linear grouped-`P2` feed-down into the scalar channel

Every real `l=2` harmonic has zero sphere average,
\[
\int_{S^2}Y_{2m}^{\rm real}(\Omega)\,d\Omega=0,
\]
so every rotational scalar observable extracted from a pure grouped real `P2` perturbation has vanishing first variation.

On the weak-axisymmetric `Y_20` branch,
\[
 x_{20}=x^{(0)}+\epsilon x^{(1)},
\qquad
 x_{21}=x^{(0)}+\frac\epsilon2 x^{(1)},
\qquad
 x_{22}=x^{(0)}-\epsilon x^{(1)},
\]
so
\[
\boxed{b_x=3a_x,}
\qquad
\boxed{\mathcal I[X,Y]=\frac{7}{10}\epsilon^2 X^{(1)}Y^{(1)}.}
\]

Therefore the exact first-order scalar theorem is
\[
\boxed{
\varepsilon_L^{(1,P_2)}=
\varepsilon_v^{(1,P_2)}=
\varepsilon_T^{(1,P_2)}=0,
}
\]
and hence
\[
\boxed{
\varepsilon_\perp^{(1,P_2)}=
\delta_\perp^{(1,P_2)}=0.
}
\]

So a pure weak grouped-lane anisotropy cannot be the first **linear** source of the scalar off-family defect.
Its scalar feed-down begins only quadratically, through the grouped bilinears
\[
\mathcal I[X,Y]=4a_X a_Y+\frac45 b_X b_Y.
\]

That leaves the linear grouped bottleneck entirely in the **direct outlet coefficients**
\[
\delta\kappa_W,
\qquad
\delta\gamma_W,
\]
not in the scalar slippage channel.

## Stage 293 — exact linear grouped outlet map

For one grouped lane, linearizing the conservative response and outgoing prefactor gives
\[
\delta u_2 = -\frac{\delta D_2+\delta D_0/9}{D_0},
\qquad
\delta u_4 = -\frac{\delta D_4+\frac29\delta D_2+\frac{5}{81}\delta D_0}{D_0},
\]
\[
\delta P_0 = \frac{\delta N_0-P_0\delta D_0}{D_0}.
\]

Using the compensated-hybrid outlet identities, the direct outlet coefficients become
\[
\boxed{
\delta\kappa_W
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}
\left(\delta D_2+\frac19\delta D_0\right),
}
\]
\[
\boxed{
\delta\gamma_W
=
-
\frac{1-\sigma_*}{9\sigma_* N_0}
\left(\delta N_0-P_0\delta D_0\right).
}
\]

So the linear grouped-even problem has collapsed to the single combination
\[
\delta D_2+\frac19\delta D_0,
\]
while the linear grouped-odd problem has collapsed to
\[
\delta N_0-P_0\delta D_0.
\]

The exact hidden-even compatibility relation is
\[
\delta u_4=
\frac89\delta u_2
iff
\boxed{
\delta D_4=
\frac23\delta D_2+
\frac1{27}\delta D_0.
}
\]

## Stage 294 — exact microscopic grouped obstructions

Substituting the full bundle decomposition
\[
D_0=K-B_0-Z_0,
\qquad
D_2=-(M+B_2+Z_2),
\qquad
D_4=-(B_4+Z_4)
\]
into the two outlet obstructions gives
\[
\boxed{
\mathcal K_A=
\mathcal W_A-
\mathcal B_A-
\mathcal Z_A,
}
\]
with
\[
\mathcal W_A=
\frac19\delta K_A-
\delta M_A,
\qquad
\mathcal B_A=
\delta B_{A,2}+
\frac19\delta B_{A,0},
\qquad
\mathcal Z_A=
\delta Z_{A,2}+
\frac19\delta Z_{A,0},
\]
and
\[
\boxed{
\mathcal G_A=
-P_0\delta K_A+
P_0\delta B_{A,0}+
\mathcal N_A,
}
\]
with
\[
\mathcal N_A=
\delta N_{A,0}+
P_0\delta Z_{A,0}.
\]

For one BdG mode,
\[
\delta B_0=
\frac{2c}{\varpi^2}\delta c-
\frac{2c^2}{\varpi^3}\delta\varpi,
\qquad
\delta B_2=
\frac{2c}{\varpi^4}\delta c-
\frac{4c^2}{\varpi^5}\delta\varpi,
\]
so
\[
\mathcal B_A=
2c\left(\frac1{\varpi^4}+rac1{9\varpi^2}\right)\delta c_A
-
2c^2\left(\frac{2}{\varpi^5}+rac1{9\varpi^3}\right)\delta\varpi_A.
\]

For one Maxwell/mixed port,
\[
Z_0=Q/\Delta,
\qquad
Z_2=(QS-G\Delta)/\Delta^2,
\qquad
N_0=P^2/\Delta^2,
\]
so
\[
\delta Z_0=
\frac{\Delta\,\delta Q-Q\,\delta\Delta}{\Delta^2},
\]
\[
\delta Z_2=
\frac{S}{\Delta^2}\delta Q+
\frac{Q}{\Delta^2}\delta S-
\frac1\Delta\delta G+
\left(\frac{G}{\Delta^2}-\frac{2QS}{\Delta^3}\right)\delta\Delta,
\]
\[
\delta N_0=
\frac{2P}{\Delta^2}\delta P-
\frac{2P^2}{\Delta^3}\delta\Delta.
\]

So the conservative Maxwell/mixed outlet obstruction depends only on the port variations
\[
\delta Q,
\qquad
\delta S,
\qquad
\delta G,
\qquad
\delta\Delta,
\qquad
\delta P.
\]

## Stage 295 — collapse to the physical weak-axisymmetric slopes

On the weak-axisymmetric grouped branch,
\[
\delta X_A=
\epsilon\lambda_A X^{(1)},
\qquad
(\lambda_{20},\lambda_{21},\lambda_{22})=
\left(1,\frac12,-1\right).
\]

Then the whole linear outlet problem collapses to two scalar amplitudes only:
\[
\boxed{
\mathfrak K_1=-D_0 u_2^{(1)},
\qquad
\mathfrak G_1=D_0 P_1.
}
\]

The direct outlet deformations inherit the same grouped signature,
\[
\delta\kappa_W^{(20)}=\epsilon\kappa_1,
\qquad
\delta\kappa_W^{(21)}=\frac\epsilon2\kappa_1,
\qquad
\delta\kappa_W^{(22)}=-\epsilon\kappa_1,
\]
with
\[
\kappa_1=
-
\frac{3(1-\sigma_*)}{\sigma_*}
 u_2^{(1)},
\]
and
\[
\delta\gamma_W^{(20)}=\epsilon\gamma_1,
\qquad
\delta\gamma_W^{(21)}=\frac\epsilon2\gamma_1,
\qquad
\delta\gamma_W^{(22)}=-\epsilon\gamma_1,
\]
with
\[
\gamma_1=
-
\frac{1-\sigma_*}{9\sigma_*}
\frac{P_1}{P_0}.
\]

Their grouped trace/anomaly defects satisfy
\[
\boxed{b_\kappa=3a_\kappa,}
\qquad
\boxed{b_\gamma=3a_\gamma.}
\]

So after the full Stage 289–295 bridge, the remaining linear grouped-anisotropy problem is no longer the raw microscopic bundle and no longer the scalar off-bundle channel.
It is simply:

> compute the actual moving-throat physical slopes
> \[
> u_2^{(1)},
> \qquad
> P_1/P_0,
> \]
> and test whether they vanish on the physical branch.
# 5PN continuation notes — Stages 296–299

This continuation takes the Stage-295 collapse
\[
u_2^{(1)},\qquad \frac{P_1}{P_0},
\]
and pushes it all the way down to the actual coherent moving-throat branch variables and then to the exact microscopic monomial / similarity-orbit structure.

A deliberate firewall was kept throughout:

- **do not modify the parent GNLS bulk PDE,**
- continue only on the moving-throat / boundary-response side,
- and keep the live 5PN / 2.5PN / 4PN obstruction on the actual grouped branch rather than inventing a new bulk medium.

## Stage 296 — actual-branch static slope compiler

Files:
- `5pn_stage296_actual_branch_static_slope_compiler.py`
- `5pn_stage296_actual_branch_static_slope_compiler_output.txt`

Main result:

The physical weak-axisymmetric continuation point is now compiled directly in the sharpest static grouped form:
\[
D_{01}=K_1-B_0^{(1)}-Z_0^{(1)},
\qquad
D_{21}=-(M_1+B_2^{(1)}+Z_2^{(1)}),
\]
\[
u_2^{(1)}=-\frac{D_{21}+u_2D_{01}}{D_0},
\qquad
\frac{P_1}{P_0}=\frac{N_1}{N_0}-\frac{D_{01}}{D_0}.
\]

On the canonical even-preserving branch,
\[
D_{21}=-\frac{D_{01}}{9},
\qquad
u_2=\frac19,
\]
so the conservative slope vanishes exactly,
\[
u_2^{(1)}=0,
\]
and the whole remaining grouped defect becomes
\[
\Xi_1=\frac{P_1}{P_0}=\bar\nu_N-\kappa_1,
\]
with
\[
\bar\nu_N=\sum_r w_r \nu_r,
\qquad
\nu_r=\delta\ln N_{A,0}^{(r)}/(\epsilon\lambda_A).
\]

So the actual static continuation point is no longer “all grouped perturbations.”
It is

1. one conservative static slope `D01/D0`,
2. one outgoing-transfer static slope `N1/N0`,
3. and, on the even-preserving branch, only their difference `Xi_1`.

## Stage 297 — coherent tracking-branch transfer shape and support-blindness

Files:
- `5pn_stage297_coherent_tracking_branch_transfer_shape.py`
- `5pn_stage297_coherent_tracking_branch_transfer_shape_output.txt`

Main result:

On the actual coherent local D/N tracking branch, the exact transfer shape is
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
\epsilon
=
\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
with selected-branch identity
\[
R_{\mathrm{target}}\mathcal T^2=\Lambda_0(1-\epsilon_\eta).
\]

The exact weak-axisymmetric grouped defect is
\[
\Xi_1
=
\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}+\frac{2\epsilon_1}{1-\epsilon},
\]
where
\[
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\varepsilon_W
-
\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
\]

The strongest theorem from this stage is the exact support-blindness result:
\[
\partial_\zeta \ln \mathcal T^2
=
\partial_\zeta \ln R_{\mathrm{target}}
=
\partial_\zeta \Xi_1
=
0.
\]

So coherent support can move the steady normalization baseline, but it cannot repair or spoil the first weak-axisymmetric grouped defect at all. That defect is carried only by the mixed/outgoing placement variables.

## Stage 298 — microscopic coherent-kernel slippage normal form

Files:
- `5pn_stage298_microscopic_coherent_slippage_normal_form.py`
- `5pn_stage298_microscopic_coherent_slippage_normal_form_output.txt`

Main result:

The coherent weak-axisymmetric defect depends only on the microscopic slippages
\[
\Sigma_Z,\qquad
\Sigma_\chi,\qquad
\Sigma_\eta,\qquad
\Sigma_\epsilon,\qquad
\Sigma_\delta,
\]
with
\[
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W,
\]
\[
\Sigma_\chi=\gamma_1+c_1-\kappa_U,
\qquad
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta,
\]
\[
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
\qquad
\Sigma_\delta=\tau_1-\kappa_U.
\]

Then the exact branch-adapted coordinates are
\[
\Sigma_{\rm tr}=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
\]
\[
\Sigma_{\rm nt}
=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-
\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]

The observable drift ledger becomes the exact triangular normal form
\[
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
\]

So the full first weak-axisymmetric grouped normalization problem is no longer a five-slippage bookkeeping problem. It is exactly the vanishing of three branch-adapted scalars:
\[
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta.
\]

The inverse reconstruction is exact:
\[
\Sigma_{\rm tr}
=
-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1,
\]
\[
\Sigma_{\rm nt}
=
\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1,
\]
\[
\Sigma_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}(\mathcal R_1+\Xi_1).
\]

So the sharp zero-defect statement is now
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.
\]

## Stage 299 — direct microscopic monomials and similarity orbit

Files:
- `5pn_stage299_microscopic_monomial_orbit_bridge.py`
- `5pn_stage299_microscopic_monomial_orbit_bridge_output.txt`
- `fivepn_stage296_299_common.py`

Main result:

The three branch-adapted defect coordinates are exactly the first logarithmic drifts of three direct microscopic monomials:
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta=\Sigma_\eta.
\]

The exact direct monomials are
\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2 T_U}{L^2 K_U}\right)^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
\]

The exact monomial-drift map on the microscopic grouped drift vector
\[
\delta\mathbf x
=
(\delta\ln\lambda_W,\delta\ln c_{\eta U},\delta\ln\gamma,\delta\ln K_U,
\delta\ln K_\eta^{(\mathrm{eff})},\delta\ln K_W^{(\mathrm{eff})},
\delta\ln\mu_W,\delta\ln T_U)
\]
is the rank-3 matrix
\[
M_*
\delta\mathbf x
=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\
\delta\ln \mathfrak C_{{\rm nt},*}\\
\delta\ln \epsilon_\eta
\end{pmatrix},
\]
with exact minor
\[
\det M_*^{(\delta\ln T_U,\delta\ln K_\eta^{(\mathrm{eff})},\delta\ln\mu_W)}
=
1+\chi_{0,*}>0.
\]

So
\[
\mathrm{rank}(M_*)=3,
\qquad
\dim\ker M_* = 5.
\]

This is the exact tangent space of a five-parameter microscopic similarity orbit preserving
\[
\mathfrak C_{{\rm tr},*},\qquad
\mathfrak C_{{\rm nt},*},\qquad
\epsilon_\eta
\]
exactly.

So the final reduced weak-axisymmetric closure theorem is now:
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x \in \ker M_*.
\]

Equivalently, on the coherent local D/N tracking branch, the actual grouped weak-axisymmetric moving-throat tangent must lie inside the exact five-dimensional monomial-preserving similarity orbit.

## What changes after Stage 299

The continuation point is now much smaller than when we started this pass.

It is no longer:

> compute all grouped-lane perturbations of
> \(
> K_A,\ M_A,\ c_{A\alpha},\ \varpi_{A\alpha},\ \Omega_{U/W,A,r},\ R_{A,r},\ g_{U/W,A,r}
> \)
> and then somehow infer \(u_2^{(1)}\) and \(P_1/P_0\).

It is now:

1. compile the actual branch static slopes into
   \[
   \kappa_1=\frac{D_{01}}{D_0},
   \qquad
   \bar\nu_N=\frac{N_1}{N_0},
   \qquad
   \Xi_1=\bar\nu_N-\kappa_1,
   \]
   on the even-preserving branch;
2. compute the three direct weak-axisymmetric monomial drifts
   \[
   \delta\ln \mathfrak C_{{\rm tr},*},\qquad
   \delta\ln \mathfrak C_{{\rm nt},*},\qquad
   \delta\ln \epsilon_\eta,
   \]
   on the actual coherent moving-throat branch;
3. test whether the physical branch tangent lies in the exact five-dimensional similarity kernel.

So the next real theorem gate is not more algebra. It is an **actual branch-tangent extraction problem** on the moving-throat PDE side.
# 5PN continuation notes — Stages 300–302

These stages continue directly from Stages 296–299.

The earlier reduction already did two key things:

1. Stage 296–298 showed that the physical first-order grouped weak-axisymmetric defect is carried by the three-scalar normal form
   
   - `Theta_1`,
   - `Xi_1 = P1/P0`,
   - `R_1`,
   
   with exact microscopic coordinates `(Sigma_tr, Sigma_nt, Sigma_eta)`.

2. Stage 299 showed that these are exactly the first logarithmic drifts of the three direct microscopic monomials
   
   - `C_tr,*`,
   - `C_nt,*`,
   - `epsilon_eta`,
   
   and that the zero-defect branch is the tangent space of the exact five-parameter monomial-preserving similarity orbit.

What was still missing was the **actual branch tester**:

> given an actual moving-throat microscopic branch state or branch tangent, how do we decide whether it is on-orbit or off-orbit without re-solving the whole drift system every time?

Stages 300–302 supply that missing tester.

---

## Stage 300 — exact finite quotient coordinates for the actual branch

Instead of working only infinitesimally, introduce the exact finite quotient coordinates between an actual microscopic branch state and a reference branch state:

\[
Q_{\rm tr}=\ln\!\frac{\mathfrak C_{{\rm tr},*}}{\mathfrak C_{{\rm tr},*,\rm ref}},
\qquad
Q_{\rm nt}=\ln\!\frac{\mathfrak C_{{\rm nt},*}}{\mathfrak C_{{\rm nt},*,\rm ref}},
\qquad
Q_{\eta}=\ln\!\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

With the direct microscopic monomials from Stage 299,

\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2 T_U}{L^2 K_U}\right)^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{\rm (eff)}(K_W^{\rm (eff)})^2}
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{\rm (eff)}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
\]

\[
\epsilon_\eta=
\frac{c_{\eta U}^2}{K_U K_\eta^{\rm (eff)}},
\]

these quotient coordinates are exact on the positive microscopic state space.

The key Stage-300 theorem is that their first linearization reproduces the Stage-299 monomial-drift vector exactly:

\[
\delta Q = M_*\,\delta x,
\]

where `M_*` is the exact rank-3 monomial-drift matrix from Stage 299.

So the finite quotient language and the infinitesimal defect language are now stitched together exactly.

---

## Stage 301 — exact affine decomposition of an actual branch tangent

Let the microscopic drift vector be ordered as

\[
\delta x =
(d\ln\lambda_W,
 d\ln c_{\eta U},
 d\ln\gamma,
 d\ln K_U,
 d\ln K_\eta^{\rm (eff)},
 d\ln K_W^{\rm (eff)},
 d\ln\mu_W,
 d\ln T_U).
\]

The Stage-299 compatibility equations can be solved exactly for the three dependent drifts

- `d ln K_eta^(eff)`,
- `d ln T_U`,
- `d ln mu_W`,

in terms of the five free drifts

- `d ln lambda_W`,
- `d ln c_etaU`,
- `d ln gamma`,
- `d ln K_U`,
- `d ln K_W^(eff)`,

plus the three quotient residuals `(q_tr, q_nt, q_eta)`.

The exact affine solve is

\[
d\ln K_\eta^{\rm (eff)} = 2\,d\ln c_{\eta U} - d\ln K_U - q_\eta,
\]

\[
d\ln T_U = d\ln K_U
+\frac{q_{\rm tr}-(1+\delta_{U,*})(d\ln\gamma+d\ln c_{\eta U}-d\ln K_U)}{1+\chi_{0,*}},
\]

\[
d\ln\mu_W
=
q_{\rm nt}-q_\eta - d\ln K_U + 2d\ln K_W^{\rm (eff)} + 2d\ln c_{\eta U} - 2d\ln\lambda_W
\]
\[
\qquad
+E_*\bigl(d\ln K_U + d\ln K_W^{\rm (eff)} - 2d\ln\gamma - 2d\ln\lambda_W\bigr)
\]
\[
\qquad
+\frac{F_*\bigl((1+\delta_{U,*})(d\ln K_U-d\ln c_{\eta U}-d\ln\gamma)+q_{\rm tr}\bigr)}{1+\chi_{0,*}}.
\]

This produces an exact decomposition

\[
\delta x_{\rm actual} = \delta x_{\rm orbit} + N\,q,
\qquad
q=(q_{\rm tr},q_{\rm nt},q_\eta)^T,
\]

where

- `delta x_orbit` lies in the five-dimensional similarity-orbit tangent space,
- `N` is a convenient exact `8 x 3` normal-basis matrix satisfying
  \[
  M_* N = I_3.
  \]

A particularly simple exact normal basis is

\[
n_{\rm tr} = \left(0,0,0,0,0,0,\frac{F_*}{1+\chi_{0,*}},\frac{1}{1+\chi_{0,*}}\right)^T,
\]

\[
n_{\rm nt} = (0,0,0,0,0,0,1,0)^T,
\]

\[
n_{\eta} = (0,0,0,0,-1,0,-1,0)^T.
\]

So any actual branch tangent now has a unique exact split:

- five tangent directions that preserve the three monomials,
- plus three normal coordinates `(q_tr, q_nt, q_eta)` measuring the off-orbit failure.

---

## Stage 302 — exact first-order defect compiler from quotient residuals

The three quotient residuals are not just convenient coordinates. They map exactly into the physical first-order defect triplet.

Using the reference-branch coefficients

\[
C_{\rm tr}
=
\frac{\chi_{0,*}\delta_{U,*}}{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
\qquad
A_{\rm tr}
=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})},
\]

and `epsilon_eta,*`, the exact first-order map is

\[
\Theta_1 = - C_{\rm tr}\, q_{\rm tr},
\]

\[
\Xi_1 = A_{\rm tr}\, q_{\rm tr} + q_{\rm nt},
\qquad
\Xi_1 = \frac{P_1}{P_0},
\]

\[
\mathcal R_1 = -\Xi_1 - \frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\, q_\eta.
\]

So in matrix form,

\[
\begin{pmatrix}
\Theta_1\\[3pt]
\Xi_1\\[3pt]
\mathcal R_1
\end{pmatrix}
=
\begin{pmatrix}
-C_{\rm tr} & 0 & 0\\
A_{\rm tr} & 1 & 0\\
-A_{\rm tr} & -1 & -\dfrac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
\end{pmatrix}
\begin{pmatrix}
q_{\rm tr}\\[3pt]
q_{\rm nt}\\[3pt]
q_\eta
\end{pmatrix}.
\]

The determinant of this defect matrix is nonzero on the physical reference branch,
so the map is exactly invertible. The inverse is

\[
q_{\rm tr} = -\frac{\Theta_1}{C_{\rm tr}},
\qquad
q_{\rm nt} = \Xi_1 + \frac{A_{\rm tr}}{C_{\rm tr}}\Theta_1,
\qquad
q_\eta = -\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}(\mathcal R_1+\Xi_1).
\]

So the actual moving-throat branch can now be tested in either of two exactly equivalent first-order languages:

1. finite quotient residuals `(q_tr, q_nt, q_eta)`,
2. physical defect triplet `(Theta_1, Xi_1, R_1)`.

And the zero-defect theorem is now the cleanest it has been:

\[
\Theta_1 = \Xi_1 = \mathcal R_1 = 0
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

---

## Net result after Stage 302

The continuation point is now smaller again.

We no longer need to test an actual moving-throat weak-axisymmetric branch against the raw eight-dimensional drift space.

The branch only needs to provide either:

- the three finite quotient residuals
  \[
  (Q_{\rm tr},Q_{\rm nt},Q_\eta),
  \]
  whose first linearization is `(q_tr,q_nt,q_eta)`,

or equivalently

- the three physical first-order defects
  \[
  (\Theta_1,\Xi_1,\mathcal R_1).
  \]

So the next honest theorem gate is no longer algebraic compression. It is finally a real branch test:

> compute the actual moving-throat branch state (or branch tangent), form its three exact quotient residuals against the reference coherent branch, and see whether those residuals vanish.

That is the cleanest actual-branch projector we have so far.
# 5PN continuation notes — Stages 303–305

These stages take the Stage 300–302 continuum-ratio orbit and push it one step closer to the actual moving-throat branch observables.

The conceptual simplification is stronger than I expected.
After using the exact coherent-branch identities already established in the moving-throat program,
the first-order weak-axisymmetric orbit-lock test no longer needs the full five-variable reduced ratio state

a) `(chi_0, delta_U, epsilon_eta, Z_W, Lambda)`,

and it does not even need the finite quotient packet `(q_tr, q_nt, q_eta)` as a primary object.

It can be written directly in the three exact branch observables

a) `R_tr`,

b) `R_target`,

c) `epsilon_eta`.

That is now the cleanest actual-branch chart for the coherent weak-axisymmetric endgame.

## Stage 303 — exact finite quotient coordinates from direct branch observables

Using the exact branch-invariant definitions

a) `T_* = R_tr^(-C_*)`,

b) `N_* = Lambda_0 (1 - epsilon_eta) R_tr^(B_*) / R_target`,

with

a) `C_* = ((1+chi_{0,*})(1+delta_{U,*})(1+chi_{0,*}+delta_{U,*}))/(chi_{0,*} delta_{U,*})`,

b) `B_* = 2(1+chi_{0,*}+delta_{U,*})/delta_{U,*}`,

the finite quotient coordinates relative to a reference branch become exactly

`q_tr = - C_* ln(R_tr / R_tr,ref)`,

`q_nt = B_* ln(R_tr / R_tr,ref)
        + ln[(1-epsilon_eta)/(1-epsilon_eta,ref)]
        - ln(R_target / R_target,ref)`,

`q_eta = ln(epsilon_eta / epsilon_eta,ref)`.

So the Stage 300 finite quotient packet has an exact direct branch-observable form.
The inverse map is also exact:

`R_tr = R_tr,ref exp(-q_tr / C_*)`,

`epsilon_eta = epsilon_eta,ref exp(q_eta)`,

`R_target = R_target,ref
            exp[-q_nt - (B_*/C_*) q_tr]
            * (1-epsilon_eta)/(1-epsilon_eta,ref)`.

This proves that the coherent branch can be charted equivalently by either

- `(R_tr, R_target, epsilon_eta)`, or
- `(q_tr, q_nt, q_eta)`.

## Stage 304 — exact first-order defect compiler in direct branch language

Linearizing the Stage 303 finite quotient map gives

`q_tr = - C_* delta ln R_tr`,

`q_nt = B_* delta ln R_tr
        - [epsilon_eta,* / (1-epsilon_eta,*)] delta ln epsilon_eta
        - delta ln R_target`,

`q_eta = delta ln epsilon_eta`.

Composing this with the exact Stage 302 quotient-to-defect compiler yields a completely triangular first-order map:

`Theta_1 = delta ln R_tr`,

`Xi_1 = - delta ln R_target
        - [epsilon_eta,* / (1-epsilon_eta,*)] delta ln epsilon_eta`,

`R_1 = delta ln R_target`.

The inverse map is equally simple:

`delta ln R_tr = Theta_1`,

`delta ln R_target = R_1`,

`delta ln epsilon_eta = -[(1-epsilon_eta,*)/epsilon_eta,*] (R_1 + Xi_1)`.

So the physical first-order defect packet on the actual coherent branch is now exactly equivalent to the three direct observable drifts

`(delta ln R_tr, delta ln R_target, delta ln epsilon_eta)`.

## Stage 305 — exact three-observable and support-blind theorem

The coherent support-compensation side from Stages 30–35 only rescales the total baseline through

`M_tr = M_mix S(zeta; epsilon)`,

`S(zeta; epsilon) = 1 + zeta(1-epsilon)/(1-zeta epsilon)`.

But Stage 303 shows that the finite quotient packet depends only on

- `R_tr`,
- `R_target`,
- `epsilon_eta`,

and Stage 304 shows the same is true for the first-order defect packet.
Therefore both packets are exactly blind to the total baseline and to the coherent support-enhancement ratio:

`partial_(M_tr) (q_tr, q_nt, q_eta) = 0`,

`partial_(M_tr) (Theta_1, Xi_1, R_1) = 0`,

`partial_(zeta) (q_tr, q_nt, q_eta) = 0`,

`partial_(zeta) (Theta_1, Xi_1, R_1) = 0`.

So the support-compensation branch and the weak-axisymmetric orbit-lock branch are exact but distinct pieces of the 5PN endgame.

### Final structural consequence of Stages 303–305

The actual coherent weak-axisymmetric zero-defect theorem now has its sharpest direct form so far:

`Theta_1 = Xi_1 = R_1 = 0`

if and only if

`delta ln R_tr = 0`,

`delta ln R_target = 0`,

`delta ln epsilon_eta = 0`.

So the actual moving-throat branch no longer needs to be tested against the older eight-variable microscopic ledger, or even explicitly against the five continuum ratios, to answer the first-order orbit-lock question.
It only has to answer three exact branch-observable questions:

1. is `R_tr` invariant?
2. is `R_target` invariant?
3. is `epsilon_eta` invariant?

And importantly, coherent support enhancement does not enter that test.
It matters for whether the selected quadrupole branch can hit the normalization target, but not for whether the branch stays on the exact first-order similarity orbit.
# 5PN continuation notes — Stages 306–308

These stages continue the Stage 300–305 projector program by eliminating the last abstract step between the quotient/defect packet and the actual coherent-kernel branch variables.

The main new result is that the first-order orbit-lock test is now written directly in the physical coherent variables

- `chi_0`
- `delta_U`
- `epsilon_eta`
- `Z_W`
- `Omega_W^2`
- `epsilon`

and is exactly factorized away from the support-compensation baseline sector.

## Stage 306 — exact coherent-branch observable map

The direct coherent-branch observables are now explicit:


a) tracking observable

`R_tr = (1 + chi_0/(1+delta_U)) / (1 + chi_0)`

which simplifies to

`R_tr = (1 + chi_0 + delta_U) / [ (1 + chi_0)(1 + delta_U) ]`.

b) transfer-shape / target observable

`T^2 = Z_W (1 + chi_0)^2 / [ Omega_W^2 (1 - epsilon)^2 ]`

`R_target = Lambda_0 Omega_W^2 (1 - epsilon_eta) (1 - epsilon)^2 / [ Z_W (1 + chi_0)^2 ]`

so the exact selected-branch identity is

`R_target * T^2 = Lambda_0 (1 - epsilon_eta)`.

c) dressing observable

`epsilon_eta = epsilon_eta`.

So the Stage 303 direct quotient packet is no longer abstract. The actual coherent branch itself supplies the finite observable packet

`(R_tr, R_target, epsilon_eta)`

from which `(q_tr, q_nt, q_eta)` follow exactly by the Stage-303 formulas.

## Stage 307 — exact physical branch drift compiler

Linearizing the Stage-306 observables gives the exact direct drift laws

`d ln R_tr = - C_tr [ (1 + delta_U) dlnchi_0 + (1 + chi_0) dlndelta_U ]`

with

`C_tr = chi_0 delta_U / [ (1 + chi_0)(1 + delta_U)(1 + chi_0 + delta_U) ]`,

`d ln R_target = dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1 + chi_0)`
`               - 2 epsilon dlnepsilon/(1 - epsilon)`
`               - epsilon_eta dlnepsilon_eta/(1 - epsilon_eta)`,

`d ln epsilon_eta = dlnepsilon_eta`.

Composing these with the Stage-304 defect map reproduces the earlier physical-branch formulas exactly:

`Theta_1 = d ln R_tr`,

`Xi_1 = - d ln R_target - [epsilon_eta/(1 - epsilon_eta)] d ln epsilon_eta`
`     = dlnZ_W - dlnOmega_W^2 + 2 chi_0 dlnchi_0/(1 + chi_0) + 2 epsilon dlnepsilon/(1 - epsilon)`,

`R_1 = d ln R_target`.

In the older slope notation this is the exact same triplet as

`Theta_1 = -( chi_0 (1 + chi_0) delta_U1 + delta_U (1 + delta_U) chi_1 )`
`          / [ (1 + chi_0)(1 + delta_U)(1 + chi_0 + delta_U) ]`,

`Xi_1 = zeta_Z - omega_W + 2 chi_1/(1 + chi_0) + 2 epsilon_1/(1 - epsilon)`,

`R_1 = omega_W - eta_1/(1 - epsilon_eta) - zeta_Z - 2 chi_1/(1 + chi_0) - 2 epsilon_1/(1 - epsilon)`.

So the actual first-order defect packet is now carried directly by the physical coherent variables rather than only by the abstract quotient packet.

## Stage 308 — exact physical-variable factorization theorem

The coherent physical branch factors exactly into two sectors:

1. orbit-lock observables

`(R_tr, R_target, epsilon_eta)`

2. support-compensation baseline

`M_tr = M_mix * S(zeta;epsilon)`

with

`S(zeta;epsilon) = 1 + zeta (1 - epsilon)/(1 - zeta epsilon)`.

The exact support-blindness identities are

`partial_zeta ln R_tr = 0`,
`partial_zeta ln R_target = 0`,
`partial_zeta ln epsilon_eta = 0`,

`partial_Mmix ln R_tr = 0`,
`partial_Mmix ln R_target = 0`,
`partial_Mmix ln epsilon_eta = 0`.

So the support lane enters only through `M_tr`; it is completely invisible to the first-order orbit-lock packet.

That turns the zero-defect theorem into one explicit system on the physical coherent variables.
The exact first-order orbit-lock conditions are

1. tracking:

`(1 + delta_U) dlnchi_0 + (1 + chi_0) dlndelta_U = 0`

2. target-ratio invariance:

`dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1 + chi_0)`
` - 2 epsilon dlnepsilon/(1 - epsilon)`
` - epsilon_eta dlnepsilon_eta/(1 - epsilon_eta) = 0`

3. dressing invariance:

`dlnepsilon_eta = 0`.

With condition 3 imposed, condition 2 simplifies to

`dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1 + chi_0)`
` - 2 epsilon dlnepsilon/(1 - epsilon) = 0`.

So the weak-axisymmetric first-order orbit-lock theorem is now a direct differential statement on the physical coherent branch variables, completely separated from the support-compensation baseline problem.

## What this changes

After Stage 308, the next theorem gate is no longer “how do we interpret the quotient packet?”
It is much sharper:

1. extract the actual branch drifts of
   `chi_0, delta_U, epsilon_eta, Z_W, Omega_W^2, epsilon`
   from the moving-throat PDE,
2. test the three direct conditions above,
3. separately test whether the support baseline sector can hit the Stage 31 / 33 / 34 support threshold.

So the 5PN endgame has now split cleanly into

- an orbit-lock / branch-invariance test, and
- a support-compensation / normalization test,

both written directly in physical branch variables.
# 5PN continuation notes — Stages 309–311

These stages continue directly from the Stage-306–308 physical branch observable compiler.
The key question there was whether the actual physical coherent branch can satisfy the nontrivial orbit-lock condition
\[
 d\ln \Omega_W^2-d\ln Z_W-
 \frac{2\chi_0}{1+\chi_0}d\ln\chi_0-
 \frac{2\epsilon}{1-\epsilon}d\epsilon = 0
\]
on the same branch where the support baseline is strong enough.

The answer is sharper than a generic “maybe.”

1. The first-order grouped weak-axisymmetric defect is carried only by the **mixed/outgoing placement variables**, not by the coherent support-enhancement baseline.
2. Exact tracking rigidity is **necessary but not sufficient**.
3. The full coherent first-order problem collapses to an exact **triangular three-scalar normal form**.
4. Equivalently, first-order zero defect is just invariance of three exact branch composites.

This is exactly the direction the later moving-throat notes had already pointed toward: the remaining weak-axisymmetric theorem gap is not support amplitude, but the drift of the physical mixed/outgoing placement variables and the corresponding exact branch composites. fileciteturn80file6turn80file13

## Stage 309 — physical tracking-branch transfer-shape defect law

On the coherent local D/N tracking branch,
\[
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\qquad
\epsilon
=\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
\[
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]
\[
R_{\rm target}
=\Lambda_0\,\Omega_W^2\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\Lambda_0=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
\]
The exact transfer shape is therefore
\[
\mathcal T^2
=\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}
=\Lambda_0\,\frac{1-\epsilon_\eta}{R_{\rm target}}.
\]

The first exact theorem is support-blindness:
\[
\partial_\zeta \ln \mathcal T^2=0,
\qquad
\partial_\zeta \ln R_{\rm target}=0,
\qquad
\partial_\zeta \Xi_1=0.
\]
So coherent support can raise the baseline through \(M_{\rm tr}\), but it cannot repair or spoil the first-order grouped defect.

Writing the weak-axisymmetric drifts as
\[
\delta\ln Z_W=\epsilon\lambda_A\zeta_Z,
\qquad
\delta\ln \Omega_W^2=\epsilon\lambda_A\omega_W,
\]
\[
\delta\chi_0=\epsilon\lambda_A\chi_1,
\qquad
\delta\epsilon=\epsilon\lambda_A\epsilon_1,
\qquad
\delta\epsilon_\eta=\epsilon\lambda_A\eta_1,
\]
one gets the exact physical defect laws
\[
\Xi_1
=\zeta_Z-\omega_W+\frac{2\chi_1}{1+\chi_0}+\frac{2\epsilon_1}{1-\epsilon},
\]
\[
\mathcal R_1
=\omega_W-\zeta_Z-\frac{2\chi_1}{1+\chi_0}-\frac{2\epsilon_1}{1-\epsilon}-\frac{\eta_1}{1-\epsilon_\eta},
\]
so that
\[
\mathcal R_1+\Xi_1=-\frac{\eta_1}{1-\epsilon_\eta}.
\]

The exact split-blocking drift is
\[
\epsilon_1
=
\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)\epsilon_{W,1}
-\frac{2\epsilon_W}{11(1+\delta_U)^2}\,\delta_{U,1}.
\]

## Stage 310 — microscopic coherent-kernel slippages

The Stage-309 physical drift variables are exact functions of the microscopic coherent-kernel drifts
\[
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1).
\]
The exact port-side translations are
\[
\zeta_Z=2\lambda_1-\kappa_\eta-\kappa_W,
\qquad
\omega_W=\kappa_W-\mu_1,
\]
\[
\chi_1=\chi_0(\gamma_1+c_1-\kappa_U),
\qquad
\eta_1=\epsilon_\eta(2c_1-\kappa_U-\kappa_\eta),
\]
\[
\epsilon_{W,1}=\epsilon_W(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W),
\qquad
\delta_{U,1}=\delta_U(\tau_1-\kappa_U).
\]

This motivates the exact microscopic slippage variables
\[
\Sigma_Z:=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W=\zeta_Z-\omega_W,
\]
\[
\Sigma_\chi:=\gamma_1+c_1-\kappa_U=\frac{\chi_1}{\chi_0},
\qquad
\Sigma_\eta:=2c_1-\kappa_U-\kappa_\eta=\frac{\eta_1}{\epsilon_\eta},
\]
\[
\Sigma_\epsilon:=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W=\frac{\epsilon_{W,1}}{\epsilon_W},
\qquad
\Sigma_\delta:=\tau_1-\kappa_U=\frac{\delta_{U,1}}{\delta_U}.
\]
Then the grouped weak-axisymmetric defect becomes
\[
\Xi_1
=\Sigma_Z
+\frac{2\chi_0}{1+\chi_0}\Sigma_\chi
+\frac{2\epsilon_W}{1-\epsilon}
\left[
\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-\frac{2\delta_U}{11(1+\delta_U)^2}\Sigma_\delta
\right].
\]

The exact tracking combination is already isolated:
\[
\Sigma_{\rm tr}:=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
\]
with
\[
\Theta_1
=-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr}.
\]
So exact tracking rigidity kills \(\Theta_1\), but it does **not** kill \(\Xi_1\) unless the remaining nontracking slippages also cooperate. That is the explicit microscopic version of the note-side warning that exact tracking-factor rigidity is not sufficient by itself. fileciteturn80file6

## Stage 311 — exact triangular normal form and direct branch composites

The exact branch-adapted nontracking slippage is
\[
\Sigma_{\rm nt}:=
\Sigma_Z
+\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\Sigma_\epsilon
-
\left[
\frac{2\chi_0}{1+\delta_U}
+\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]
Then the three observables collapse to the exact triangular system
\[
\Theta_1=-C_{\rm tr}\,\Sigma_{\rm tr},
\qquad
C_{\rm tr}:=\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)},
\]
\[
\Xi_1=A_{\rm tr}\,\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
A_{\rm tr}:=\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)},
\]
\[
\mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
\]
This inverts exactly to
\[
\Sigma_{\rm tr}
=-\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U}\,\Theta_1,
\]
\[
\Sigma_{\rm nt}
=\Xi_1+\frac{2(1+\chi_0+\delta_U)}{\delta_U}\,\Theta_1,
\]
\[
\Sigma_\eta
=-\frac{1-\epsilon_\eta}{\epsilon_\eta}(\mathcal R_1+\Xi_1).
\]

So the coherent first-order problem is no longer a five-slippage bookkeeping problem. It is an exact three-scalar normal form.

There is also an exact direct branch-composite theorem. Define
\[
C_*:=\frac{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}{\chi_0\delta_U},
\qquad
B_*:=\frac{2(1+\chi_0+\delta_U)}{\delta_U},
\]
\[
\mathfrak T_*:=R_{\rm tr}^{-C_*},
\qquad
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*}.
\]
Then their logarithmic drifts are exactly
\[
\delta\ln \mathfrak T_* = \Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak N_* = \Sigma_{\rm nt},
\qquad
\delta\ln \epsilon_\eta = \Sigma_\eta.
\]
Therefore
\[
\Theta_1=\Xi_1=\mathcal R_1+\Xi_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0
\iff
\delta\ln R_{\rm tr}=0,
\ \delta\ln(\mathcal T^2 R_{\rm tr}^{B_*})=0,
\ \delta\ln\epsilon_\eta=0.
\]
So the actual first-order coherent theorem gate is now exactly the invariance of three branch composites, not the adjustment of the support baseline. That matches the compact continuation point already recorded in the moving-throat ledger. fileciteturn80file13turn80file15turn80file11

## Bottom line after Stages 309–311

The program did **not** hit a dead end here.
It hit a sharper separation of roles.

- The coherent support lane can still rescue the **steady normalization baseline**.
- But it is exactly absent from the **first-order grouped weak-axisymmetric defect**.
- The first-order defect is carried only by the mixed/outgoing placement variables and their exact branch composites.

So the next honest theorem gate is now very small:

> compute the actual first grouped weak-axisymmetric drifts of
> \[
> R_{\rm tr},
> \qquad
> \mathcal T^2 R_{\rm tr}^{B_*},
> \qquad
> \epsilon_\eta
> \]
> on the physical moving-throat branch.

If all three are invariant, the coherent first-order grouped defect vanishes automatically. If not, the failure is localized immediately into tracking, nontracking transfer shape, or dressing. fileciteturn80file13turn80file15
# 5PN continuation notes — Stages 312–314

These stages continue directly from the Stage-309–311 coherent branch defect normal form.

At Stage 311 the first coherent weak-axisymmetric theorem had already collapsed to
three exact branch-adapted scalars,
\[
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta),
\]
and equivalently to invariance of three direct branch composites.

The next honest step was to remove the remaining dependence on those intermediate
composites and rewrite the whole zero-defect test directly in the microscopic
coherent-kernel variables. That is what Stages 312–314 do.

## Stage 312 — direct microscopic monomial coordinates

Using the coherent-kernel definitions
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2 T_U}{L^2K_U},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\]
\[
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\qquad
\frac{Z_W}{\Omega_W^2}=\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2},
\]
with frozen reference-branch exponents
\[
E_*=
\frac{2\epsilon_{W,*}(11+9\delta_{U,*})}{11(1-\epsilon_*)(1+\delta_{U,*})},
\qquad
F_*=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}{11(1-\epsilon_*)(1+\delta_{U,*})^2},
\]
the three direct microscopic monomials are
\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}}
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*}
}
\]
\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
}
\]

With microscopic grouped weak-axisymmetric log-drifts
\[
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1),
\]
the exact logarithmic derivatives are
\[
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
\]
So the coherent first-order zero-defect theorem becomes
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln \mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln \mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]

That is the smallest exact microscopic formulation reached so far.

## Stage 313 — exact microscopic compatibility ledger

Once the three direct monomials are used, the zero-defect theorem is equivalent
to the explicit microscopic compatibility system
\[
\boxed{
(1+\chi_{0,*})(\tau_1-\kappa_U)
+
(1+\delta_{U,*})(\gamma_1+c_1-\kappa_U)=0,
}
\]
\[
\boxed{
2c_1-\kappa_U-\kappa_\eta=0,
}
\]
\[
\boxed{
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
+
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*(\tau_1-\kappa_U)=0.
}
\]

These solve exactly for the three dependent drifts:
\[
\boxed{
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
(\gamma_1+c_1-\kappa_U),
}
\]
\[
\boxed{
\kappa_\eta = 2c_1-\kappa_U,
}
\]
\[
\boxed{
\mu_1
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\gamma_1+c_1-\kappa_U).
}
\]

So the zero-defect branch is no longer an abstract slippage statement. It is a
codimension-3 microscopic rigidity ledger inside the eight-drift kernel space.

## Stage 314 — exact microscopic similarity orbit

The compatibility ledger is not an isolated fine-tuning condition. It is the
logarithmic tangent form of a finite five-parameter multiplicative similarity orbit.

Take free logarithmic orbit parameters
\[
(u_\lambda,u_c,u_\gamma,u_{K_U},u_{K_W}).
\]
The dependent exponents are fixed by exact monomial preservation:
\[
\boxed{
u_{K_\eta}=2u_c-u_{K_U},}
\]
\[
\boxed{
u_{T_U}=u_{K_U}-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(u_\gamma+u_c-u_{K_U}),}
\]
\[
\boxed{
u_{\mu_W}=2u_c-u_{K_U}+2u_{K_W}-2u_\lambda
-
E_*(2u_\gamma+2u_\lambda-u_{K_U}-u_{K_W})
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(u_\gamma+u_c-u_{K_U}).}
\]

The finite orbit acts multiplicatively by
\[
\lambda_W\mapsto e^{u_\lambda}\lambda_W,
\qquad
c_{\eta U}\mapsto e^{u_c}c_{\eta U},
\qquad
\gamma\mapsto e^{u_\gamma}\gamma,
\]
\[
K_U\mapsto e^{u_{K_U}}K_U,
\qquad
K_\eta^{(\mathrm{eff})}\mapsto e^{u_{K_\eta}}K_\eta^{(\mathrm{eff})},
\qquad
K_W^{(\mathrm{eff})}\mapsto e^{u_{K_W}}K_W^{(\mathrm{eff})},
\]
\[
\mu_W\mapsto e^{u_{\mu_W}}\mu_W,
\qquad
T_U\mapsto e^{u_{T_U}}T_U.
\]

On this finite family,
\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta
\]
are preserved **exactly**.

The tangent matrix of this orbit has rank `5`, while the monomial-drift matrix has
rank `3`, so
\[
\dim\ker M_*=5,
\]
and the orbit tangent space is precisely the zero-defect tangent space.

Therefore the coherent weak-axisymmetric theorem can now be read geometrically:
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x\in T_{\mathrm{id}}\mathcal G_*,
}
\]
where `\mathcal G_*` is the exact five-parameter monomial-preserving similarity orbit.

## Net result after Stage 314

The algebraic side of the coherent weak-axisymmetric branch is now effectively finished.

1. The remaining defect is measured by the drift of three exact microscopic monomials.
2. The zero-defect conditions are an explicit three-equation compatibility ledger.
3. That ledger is exactly the tangent-space condition for a finite five-parameter similarity orbit.

So the next honest theorem gate is no longer algebraic compression.
It is the actual branch-realization question:

> does the true moving-throat grouped weak-axisymmetric branch move tangent to the exact monomial-preserving similarity orbit `\mathcal G_*`?

That is the cleanest continuation point we have right now.
# 5PN continuation notes — Stages 315–318

These stages finally turn the Stage 312–314 microscopic similarity-orbit theorem into a **usable actual-branch tester** in the reduced continuum variables that the moving-throat branch is expected to hand us.

The key reduction is that the actual coherent branch does **not** need the full eight-variable microscopic kernel state to test weak-axisymmetric orbit lock. It only needs the five-state packet

\[
(\chi_0,\ \delta_U,\ \widehat Z_W,\ \epsilon_W,\ \epsilon_\eta),
\]

where

\[
\widehat Z_W := \frac{\lambda_W^2 \mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
= \frac{Z_W}{\Omega_W^2}.
\]

So the actual-branch weak-axisymmetric theorem gate is now formulated directly in the language of the coherent continuum branch rather than in the older full microscopic kernel ledger.

## Stage 315 — exact continuum monomial drift map

The direct Stage-312 monomials become, in the reduced continuum variables,

\[
\mathfrak C_{{\rm tr},*}
=
\chi_0^{1+\delta_{U,*}}\,\delta_U^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\widehat Z_W\,\epsilon_W^{E_*}\,\delta_U^{-F_*},
\]

\[
\epsilon_\eta = \epsilon_\eta.
\]

Therefore the finite quotient packet is exactly

\[
q_{\rm tr}
=
(1+\delta_{U,*})\ln\!\frac{\chi_0}{\chi_{0,\rm ref}}
+
(1+\chi_{0,*})\ln\!\frac{\delta_U}{\delta_{U,\rm ref}},
\]

\[
q_{\rm nt}
=
\ln\!\frac{\widehat Z_W}{\widehat Z_{W,\rm ref}}
+
E_*\ln\!\frac{\epsilon_W}{\epsilon_{W,\rm ref}}
-
F_*\ln\!\frac{\delta_U}{\delta_{U,\rm ref}},
\]

\[
q_\eta
=
\ln\!\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

At infinitesimal level, the branch-adapted monomial drifts are

\[
\Sigma_{\rm tr}
=
(1+\delta_{U,*})\,d\ln\chi_0
+
(1+\chi_{0,*})\,d\ln\delta_U,
\]

\[
\Sigma_{\rm nt}
=
 d\ln\widehat Z_W + E_*\,d\ln\epsilon_W - F_*\,d\ln\delta_U,
\]

\[
\Sigma_\eta = d\ln\epsilon_\eta.
\]

So the actual orbit-lock tester now acts on the reduced five-drift packet

\[
(d\ln\chi_0,\ d\ln\delta_U,\ d\ln\widehat Z_W,\ d\ln\epsilon_W,\ d\ln\epsilon_\eta).
\]

## Stage 316 — exact support-blindness of the continuum orbit tester

The coherent support-enhancement branch still has

\[
M_{\rm tr}=M_{\rm mix}\,S(\zeta;\epsilon),
\qquad
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]

and the explicit twin-family support tower still determines the physical support ratio

\[
\zeta_n^{(\rm twin)}
=
\frac{1}{(2n+1)^2\bigl(1+x\,n(n+1)\bigr)}.
\]

But the actual continuum quotient packet and the actual defect packet are exactly blind to

- the support-enhancement ratio \(\zeta\),
- the total baseline \(M_{\rm tr}\),
- and the explicit twin-family stiffness \(x\).

Symbolically,

\[
\partial_\zeta(q_{\rm tr},q_{\rm nt},q_\eta)=0,
\qquad
\partial_{M_{\rm tr}}(q_{\rm tr},q_{\rm nt},q_\eta)=0,
\qquad
\partial_x(q_{\rm tr},q_{\rm nt},q_\eta)=0,
\]

and likewise for the physical defect triplet

\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]

So support compensation and support-harmonic selection belong entirely to the **isotropic normalization** side of the 5PN endgame. They do **not** move the coherent weak-axisymmetric branch on or off the exact similarity orbit.

## Stage 317 — exact reduced actual-branch orbit tester

The five-drift packet maps linearly to the monomial-drift triple through the exact reduced projector

\[
M_{\rm red}:
(d\ln\chi_0,\ d\ln\delta_U,\ d\ln\widehat Z_W,\ d\ln\epsilon_W,\ d\ln\epsilon_\eta)
\mapsto
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta).
\]

Then the physical defect triplet is

\[
\Theta_1 = -C_{\rm tr}\,\Sigma_{\rm tr},
\]

\[
\Xi_1 = A_{\rm tr}\,\Sigma_{\rm tr}+\Sigma_{\rm nt},
\]

\[
\mathcal R_1 = -\Xi_1 - \frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\Sigma_\eta,
\]

with

\[
C_{\rm tr}=
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
\qquad
A_{\rm tr}=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}.
\]

The reduced tester has exact rank `3`, so its zero-defect kernel has dimension `2` inside the five-drift continuum space.

A convenient exact basis is

\[
v_\chi
=
\left(
1,
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}},
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}},
0,
0
\right),
\]

\[
v_{\epsilon_W}
=
(0,0,-E_*,1,0).
\]

Equivalently, exact orbit lock on the actual coherent branch is the codrift system

\[
d\ln\delta_U
=
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}\,d\ln\chi_0,
\]

\[
d\ln\widehat Z_W
=
- E_*\,d\ln\epsilon_W + F_*\,d\ln\delta_U,
\]

\[
d\ln\epsilon_\eta = 0.
\]

So the actual coherent branch lies on the exact similarity orbit iff its reduced drift vector lies in the 2-plane spanned by \(v_\chi\) and \(v_{\epsilon_W}\).

## Stage 318 — packaged actual-branch tester API

The new tester is now ready in two equivalent forms.

### Finite tester
Input:

\[
(\chi_0,\delta_U,\widehat Z_W,\epsilon_W,\epsilon_\eta)
\]

plus a reference branch state and the frozen reference exponents.

Output:

\[
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

### Infinitesimal tester
Input:

\[
(d\ln\chi_0,\ d\ln\delta_U,\ d\ln\widehat Z_W,\ d\ln\epsilon_W,\ d\ln\epsilon_\eta).
\]

Output:

either

\[
(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)
\]

or

\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]

The exact zero-defect packet is the 2-plane above, so the physical branch no longer needs to be checked against a larger microscopic state space.

## Net result after Stages 315–318

These stages close the exact step we had been missing after Stage 312–314:

1. the microscopic similarity-orbit theorem is now expressed directly in the reduced continuum variables of the actual coherent branch;
2. support enhancement and twin/non-twin support selection are proven to be orthogonal to the weak-axisymmetric orbit-lock problem;
3. the actual branch tester is now a concrete five-variable input / three-variable output object;
4. and the exact zero-defect branch inside that five-variable continuum state space is a 2-plane.

So the next theorem gate is no longer “how do we test the real branch?”
It is narrower:

> extract the actual moving-throat branch drifts of
> \(\chi_0,\delta_U,\widehat Z_W,\epsilon_W,\epsilon_\eta\),
> feed them into the new tester,
> and check whether the branch lies in the exact reduced orbit plane while simultaneously satisfying the isotropic normalization branch conditions.
# 5PN continuation notes — Stages 319–322

These stages take the Stage-315–318 reduced continuum orbit tester and complete the missing bridge back to the actual microscopic coherent-kernel variables.

The net effect is that the weak-axisymmetric orbit-lock problem is no longer formulated in an abstract five-drift packet. It is now directly readable from the microscopic grouped drifts of

- `lambda_W`
- `c_etaU`
- `gamma`
- `K_U`
- `K_eta^(eff)`
- `K_W^(eff)`
- `mu_W`
- `T_U`

while the explicit support variables `lambda_phi` and `K_phi^(eff)` are shown to live only on the separate isotropic support-compensation branch.

## Stage 319 — exact microscopic -> reduced drift extractor

Files:
- `fivepn_stage319_322_common.py`
- `5pn_stage319_microscopic_reduced_drift_extractor.py`
- `5pn_stage319_microscopic_reduced_drift_extractor_output.txt`

The exact reduced Stage-318 tester inputs are

- `dln chi0`
- `dln deltaU`
- `dln ZhatW`
- `dln epsilonW`
- `dln epsilon_eta`

with

`chi0 = gamma c_etaU / K_U`,

`deltaU = pi^2 T_U / (L^2 K_U)`,

`ZhatW = lambda_W^2 mu_W / (K_eta^(eff) (K_W^(eff))^2)`,

`epsilonW = gamma^2 lambda_W^2 sigma / (K_U K_W^(eff))`,

`epsilon_eta = c_etaU^2 / (K_U K_eta^(eff))`.

The exact microscopic grouped-drift extractor is therefore

`dln chi0     = gamma1 + c1 - kappaU`,

`dln deltaU   = tau1 - kappaU`,

`dln ZhatW    = 2 lambda1 + mu1 - kappaEta - 2 kappaW`,

`dln epsilonW = 2 gamma1 + 2 lambda1 - kappaU - kappaW`,

`dln epsilon_eta = 2 c1 - kappaU - kappaEta`.

So the reduced five-drift packet is an exact linear image of the eight microscopic grouped drifts. No extra reduced closure assumptions are needed to pass from the microscopic branch to the Stage-318 tester.

## Stage 320 — exact microscopic actual-branch tester

Files:
- `5pn_stage320_microscopic_actual_branch_tester.py`
- `5pn_stage320_microscopic_actual_branch_tester_output.txt`

Composing the Stage-319 extractor with the Stage-317 reduced tester gives the direct microscopic monomial-drift packet

`Sigma_tr`
`= (1+deltaU_*) (gamma1 + c1 - kappaU) + (1+chi0_*) (tau1 - kappaU)`,

`Sigma_nt`
`= (2 lambda1 + mu1 - kappaEta - 2 kappaW)`
`  + E_* (2 gamma1 + 2 lambda1 - kappaU - kappaW)`
`  - F_* (tau1 - kappaU)`,

`Sigma_eta = 2 c1 - kappaU - kappaEta`.

So the physical defect packet is still the same triangular normal form,

`Theta_1 = - C_tr,* Sigma_tr`,

`Xi_1    =   A_tr,* Sigma_tr + Sigma_nt`,

`R_1     = - Xi_1 - epsilon_eta,* Sigma_eta/(1-epsilon_eta,*)`.

The direct microscopic tester has exact rank `3` inside the eight-drift kernel space, so its zero-defect kernel is `5`-dimensional.

A convenient exact solve is

`tau1 = kappaU - (1+deltaU_*) (gamma1 + c1 - kappaU)/(1+chi0_*)`,

`kappaEta = 2 c1 - kappaU`,

`mu1 = 2 c1 - kappaU + 2 kappaW - 2 lambda1`
`      - E_* (2 gamma1 + 2 lambda1 - kappaU - kappaW)`
`      - F_* (1+deltaU_*) (gamma1 + c1 - kappaU)/(1+chi0_*)`.

So the codimension-3 similarity-orbit closure from Stages 312–314 now appears directly inside the Stage-318 actual-branch tester language.

## Stage 321 — finite microscopic branch packet

Files:
- `5pn_stage321_finite_microscopic_branch_packet.py`
- `5pn_stage321_finite_microscopic_branch_packet_output.txt`

The finite Stage-318 quotient packet is now evaluated directly on the microscopic coherent-kernel state

`(lambda_W, c_etaU, gamma, K_U, K_eta^(eff), K_W^(eff), mu_W, T_U)`.

The three direct microscopic monomials are

`C_tr,*`
`= (gamma c_etaU / K_U)^(1+deltaU_*)`
`  (pi^2 T_U / (L^2 K_U))^(1+chi0_*)`,

`C_nt,*`
`= [lambda_W^2 mu_W / (K_eta^(eff) (K_W^(eff))^2)]`
`  [gamma^2 lambda_W^2 sigma / (K_U K_W^(eff))]^(E_*)`
`  [pi^2 T_U / (L^2 K_U)]^(-F_*)`,

`epsilon_eta = c_etaU^2 / (K_U K_eta^(eff))`.

And the finite quotient packet is exactly

`q_tr  = ln(C_tr,* / C_tr,*,ref)`,

`q_nt  = ln(C_nt,* / C_nt,*,ref)`,

`q_eta = ln(epsilon_eta / epsilon_eta,ref)`.

The direct physical branch observables are also explicit in microscopic variables:

`R_tr = (1 + chi0/(1+deltaU)) / (1 + chi0)`,

`R_target = Lambda_0 (1-epsilon_eta) (1-epsilon)^2 / [ ZhatW (1+chi0)^2 ]`,

with

`epsilon = epsilonW [ 1 - (2/11) deltaU/(1+deltaU) ]`.

So the actual moving-throat branch no longer needs a separate reduced-state postulate before the finite quotient packet can be compiled: it is already a direct log-ratio packet of three microscopic monomials and three direct microscopic branch observables.

## Stage 322 — exact support-blind orbit lock vs support-sensitive isotropic split

Files:
- `5pn_stage322_support_blind_orbit_split.py`
- `5pn_stage322_support_blind_orbit_split_output.txt`

The final bridge of this batch is the exact separation theorem.

Introduce the explicit support variables

- `lambda_phi`
- `K_phi^(eff)`

through the coherent support ratio

`zeta = lambda_phi^2 K_W^(eff) / (lambda_W^2 K_phi^(eff))`.

Then:

1. the orbit-lock observables

   `R_tr`, `R_target`, `epsilon_eta`

   are exactly independent of `lambda_phi` and `K_phi^(eff)`;

2. the finite quotient packet `(q_tr, q_nt, q_eta)` is also exactly independent of them;

3. but the isotropic support lane depends on them through the single enhancement ratio

   `zeta = lambda_phi^2 K_W^(eff) / (lambda_W^2 K_phi^(eff))`,

   with exact support-enhancement factor

   `S(zeta;epsilon) = 1 + zeta (1-epsilon)/(1-zeta epsilon)`.

The exact logarithmic support sensitivities are

`d ln zeta / d ln lambda_phi = 2`,

`d ln zeta / d ln K_phi^(eff) = -1`,

and therefore

`d M_tr / d ln lambda_phi > 0`,

`d M_tr / d ln K_phi^(eff) < 0`

on the physical branch.

So the weak-axisymmetric orbit-lock test and the isotropic support-compensation / normalization test are rigorously separated at the microscopic level.

## Bottom line after Stage 322

The continuation point has tightened again.

1. The weak-axisymmetric orbit-lock problem is now a direct tester on the microscopic drift vector of
   `(lambda_W, c_etaU, gamma, K_U, K_eta^(eff), K_W^(eff), mu_W, T_U)`.
2. Its exact zero-defect kernel is still the codimension-3 similarity-orbit closure already isolated in the notes.
3. The explicit support variables `(lambda_phi, K_phi^(eff))` do **not** enter that orbit-lock packet at all.
4. They enter only through the separate isotropic support-enhancement lane `zeta`.

So the next honest theorem gate is now even smaller:

> extract the actual microscopic branch drifts and finite kernel states from the completed moving-throat operator, feed them into the Stage-320/321 orbit packet, and test the support-enhancement lane separately through `zeta`.
# 5PN continuation notes — Stages 323–325

These stages take the Stage-319–322 microscopic orbit tester and splice it directly into the coherent local-kernel placement map that the moving-throat compact program already uses for the physical branch.

The net effect is that the weak-axisymmetric orbit-lock packet no longer needs the full microscopic kernel tuple

- `lambda_W`
- `c_etaU`
- `gamma`
- `K_U`
- `K_eta^(eff)`
- `K_W^(eff)`
- `mu_W`
- `T_U`

and it no longer needs the intermediate reduced state variable `ZhatW` as an independent input either.

On the coherent branch, everything now factors through the six physical placement variables

- `chi0`
- `deltaU`
- `Z_W`
- `epsilon_W`
- `epsilon_eta`
- `Lambda`

plus the separate isotropic support ratio

- `zeta`

which lives only on the support-compensation lane.

## Stage 323 — coherent placement map already determines the reduced Stage-318 state

The exact coherent placement variables are

\[
\chi_0,
\qquad
\delta_U,
\qquad
Z_W,
\qquad
\epsilon_W,
\qquad
\epsilon_\eta,
\qquad
\Lambda,
\qquad
\zeta.
\]

The key bridge is

\[
\Lambda_0 := \frac{27\pi^2 G c_s^5}{20 a^5 c^5},
\qquad
\widehat Z_W = Z_W\frac{\Lambda_0}{\Lambda}.
\]

So the reduced Stage-318 state is already reconstructed from the coherent placement map by

\[
(\chi_0,\delta_U,\widehat Z_W,\epsilon_W,\epsilon_\eta)
=
\left(
\chi_0,
\delta_U,
Z_W\frac{\Lambda_0}{\Lambda},
\epsilon_W,
\epsilon_\eta
\right).
\]

The direct coherent-branch observables are

\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]

\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]

\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]

And because

\[
R_{\rm target}
=
\Lambda_0\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{\widehat Z_W(1+\chi_0)^2},
\]

the placement-map and reduced-state versions of the physical branch are exactly equivalent.

Most importantly, `zeta` does **not** enter

- `ZhatW`,
- `epsilon`,
- `R_tr`,
- or `R_target`.

So the coherent support ratio is absent from the extracted orbit state.

## Stage 324 — exact finite and infinitesimal orbit packet in coherent placement variables

The finite quotient packet is now

\[
q_{\rm tr}
=
(1+\delta_{U,*})\ln\frac{\chi_0}{\chi_{0,\rm ref}}
+
(1+\chi_{0,*})\ln\frac{\delta_U}{\delta_{U,\rm ref}},
\]

\[
q_{\rm nt}
=
\ln\frac{Z_W}{Z_{W,\rm ref}}
-
\ln\frac{\Lambda}{\Lambda_{\rm ref}}
+
E_*\ln\frac{\epsilon_W}{\epsilon_{W,\rm ref}}
-
F_*\ln\frac{\delta_U}{\delta_{U,\rm ref}},
\]

\[
q_\eta = \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

So the only change relative to the older reduced-state formula is the exact replacement

\[
\ln \widehat Z_W = \ln Z_W - \ln \Lambda + \ln \Lambda_0,
\]

with the constant `Lambda_0` dropping out of the quotient packet.

At infinitesimal level, the exact coherent monomial-drift packet becomes

\[
\Sigma_{\rm tr}
=
(1+\delta_{U,*})\,d\ln\chi_0
+
(1+\chi_{0,*})\,d\ln\delta_U,
\]

\[
\Sigma_{\rm nt}
=
(d\ln Z_W-d\ln\Lambda)
+
E_*\,d\ln\epsilon_W
-
F_*\,d\ln\delta_U,
\]

\[
\Sigma_\eta = d\ln\epsilon_\eta.
\]

Then the physical defect packet is still the same triangular normal form,

\[
\Theta_1 = -C_{\rm tr,*}\Sigma_{\rm tr},
\qquad
\Xi_1 = A_{\rm tr,*}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1 = -\Xi_1 - \frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\Sigma_\eta.
\]

The strongest exact bridge in this stage is the direct-observable drift form:

\[
\Theta_1 = d\ln R_{\rm tr},
\]

\[
\Xi_1 = -d\ln R_{\rm target} - \frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}d\ln\epsilon_\eta,
\]

\[
\mathcal R_1 = d\ln R_{\rm target}.
\]

So the coherent orbit-lock theorem on the physical branch is now simply

\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
 d\ln R_{\rm tr}=0,
\quad
 d\ln R_{\rm target}=0,
\quad
 d\ln\epsilon_\eta=0.
\]

## Stage 325 — exact two-packet split on the coherent branch

The coherent branch now separates exactly into:

### Orbit-lock packet
Finite:

\[
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

Infinitesimal:

\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]

This packet depends only on

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda)
\]

and is exactly blind to `zeta`.

### Support / normalization packet
The split-support variables are

\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]

\[
M_{\rm mix}
=
\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
\]

\[
S(\zeta;\epsilon)
=
1+
\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\qquad
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\]

\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]

The exact mixed-only product law is

\[
R_{\rm target}M_{\rm mix} = \frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

The crucial separation statement is therefore:

- `zeta` changes only the available baseline through `S(zeta;epsilon)`,
- `zeta` leaves the finite orbit packet unchanged,
- and `zeta` leaves the infinitesimal orbit defect packet unchanged.

So support compensation cannot rescue orbit-lock failure, and orbit lock does not determine support enhancement. The completed moving-throat PDE has to satisfy the two packets separately on the same physical branch.

## Bottom line after Stage 325

The next honest theorem gate is now smaller again.

The completed moving-throat operator no longer needs to be interpreted through the full microscopic 8-tuple before the orbit-lock problem can even be asked. It only needs to provide the six coherent placement variables

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda),
\]

and the separate support ratio

\[
\zeta.
\]

Then:

1. the weak-axisymmetric orbit-lock packet is compiled exactly from the first six,
2. the support / normalization packet is compiled from all seven,
3. and the two tests are rigorously separated.

So the remaining gap is no longer algebraic compression. It is actual branch realization of these placement variables on the completed moving-throat operator.
# 5PN continuation notes — Stages 326–328

These stages perform the extraction that Stage 323–325 was still waiting for.

Before this batch, the coherent branch was already split into

- an **orbit-lock packet** built from
  \[(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda),\]
- and a separate **support/normalization packet** that adds only
  \[\zeta.\]

What was still missing was the exact bridge from the **microscopic kernel state and its grouped weak-axisymmetric drifts** into those actual branch variables.

These stages close that gap.

---

## Stage 326 — exact microscopic coherent placement-state extractor

Files:
- `fivepn_stage326_328_common.py`
- `5pn_stage326_microscopic_coherent_placement_state.py`
- `5pn_stage326_microscopic_coherent_placement_state_output.txt`

The actual coherent placement state is now extracted directly from the microscopic kernel state
\[(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U,\lambda_\phi,K_\phi^{(\mathrm{eff})}).\]

The exact formulas are

\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\qquad
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\]
\[
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\Lambda=\Lambda_0\frac{K_W^{(\mathrm{eff})}}{\mu_W},
\]
\[
\zeta=\frac{\lambda_\phi^2K_W^{(\mathrm{eff})}}{\lambda_W^2K_\phi^{(\mathrm{eff})}},
\qquad
\Lambda_0=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
\]

The older reduced Stage-318 variable is recovered exactly as
\[
\widehat Z_W = Z_W\frac{\Lambda_0}{\Lambda}
=\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}.
\]
So the earlier reduced state was not wrong; it was just hiding the quotient structure of two actual branch observables.

The exact coherent observables then become
\[
R_{\rm tr}=\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]
\[
\epsilon=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
\[
R_{\rm target}=\Lambda\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]
The support packet also compiles directly in microscopic variables through
\[
M_{\rm mix}=\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
\qquad
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\]
with
\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\qquad
R_{\rm target}M_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

So after Stage 326, the actual coherent branch is no longer expressed by an abstract reduced packet. It is a direct quotient state of the microscopic kernel.

---

## Stage 327 — exact microscopic coherent placement-drift extractor

Files:
- `5pn_stage327_microscopic_coherent_placement_drifts.py`
- `5pn_stage327_microscopic_coherent_placement_drifts_output.txt`

The grouped weak-axisymmetric drifts of the actual placement variables are now explicit.
Writing the microscopic drift variables as
\[(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1,\phi_1,\kappa_\phi),\]
we get

\[
\delta\ln\chi_0=\gamma_1+c_1-\kappa_U,
\qquad
\delta\ln\delta_U=\tau_1-\kappa_U,
\]
\[
\delta\ln Z_W=2\lambda_1-\kappa_\eta-\kappa_W,
\qquad
\delta\ln\epsilon_W=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
\]
\[
\delta\ln\epsilon_\eta=2c_1-\kappa_U-\kappa_\eta,
\qquad
\delta\ln\Lambda=\kappa_W-\mu_1,
\]
\[
\delta\ln\zeta=2\phi_1-2\lambda_1+\kappa_W-\kappa_\phi.
\]

The older reduced drift is recovered exactly:
\[
\delta\ln\widehat Z_W
=
\delta\ln Z_W-
\delta\ln\Lambda
=
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W.
\]

This stage also makes the packet split completely sharp at drift level:

- the **orbit side** uses only
  \[(\delta\ln\chi_0,\delta\ln\delta_U,\delta\ln Z_W,\delta\ln\epsilon_W,\delta\ln\epsilon_\eta,\delta\ln\Lambda),\]
- the **support side** adds only
  \[\delta\ln\zeta.\]

So support transport cannot hide inside the orbit packet even infinitesimally.

---

## Stage 328 — exact microscopic actual-branch packet compiler

Files:
- `5pn_stage328_microscopic_actual_branch_packet_compiler.py`
- `5pn_stage328_microscopic_actual_branch_packet_compiler_output.txt`

This stage closes the loop between the older microscopic monomial theorem and the newer actual coherent placement-map packet.

### 1. Finite orbit packet

The finite quotient packet
\[(q_{\rm tr},q_{\rm nt},q_\eta)\]
compiled from the actual placement state is exactly the finite log-ratio packet of the three direct microscopic monomials. There is no mismatch between the two languages.

### 2. Infinitesimal orbit defect packet

Feeding the actual placement drifts into the Stage-324 packet compiler gives
\[
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta,
\]
and from them
\[
\Theta_1,
\qquad
\Xi_1,
\qquad
\mathcal R_1.
\]
The exact theorem proved in the script is that these are **identical** to the older Stage-313/314 microscopic compatibility ledger.

So the actual coherent placement map introduces no hidden extra weak-axisymmetric obstruction. The packet is the same packet in two equivalent coordinate systems.

### 3. Direct observable zero-defect form

On the actual microscopic branch the first-order orbit-lock theorem is exactly
\[
\delta\ln R_{\rm tr}=0,
\qquad
\delta\ln R_{\rm target}=0,
\qquad
\delta\ln\epsilon_\eta=0.
\]

### 4. Final two-packet split

The actual coherent moving-throat branch now ends at the exact microscopic two-packet decomposition

- **orbit packet**
  \[(q_{\rm tr},q_{\rm nt},q_\eta)\]
  or equivalently
  \[(\Theta_1,\Xi_1,\mathcal R_1),\]
- **support packet**
  \[(\zeta;M_{\rm mix},S,M_{\rm tr}).\]

The orbit packet is blind to `lambda_phi`, `K_phi^(eff)`, `phi_1`, and `kappa_phi`; the support packet carries them all.

---

## Where this leaves the 5PN / moving-throat program

The extraction step is now done in the sharpest exact form we have had so far.

1. The actual coherent branch state is a **7-scalar microscopic quotient state**:
   \[(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta).\]
2. The weak-axisymmetric orbit-lock theorem depends only on the first six.
3. The support/normalization theorem adds only `\zeta`.
4. The actual orbit packet compiled from those six scalars is **exactly the same** packet as the older microscopic monomial-rigidity ledger.

So there is no further hidden algebraic bottleneck between the microscopic kernel variables and the actual coherent branch observables.

The remaining honest theorem gate is now purely PDE-side:

> compute the actual finite state and grouped weak-axisymmetric drifts of
> \[(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta)\]
> on the completed moving-throat operator and test them against the exact orbit packet and support packet above.
# 5PN continuation notes — Stages 329–332

These stages take the Stage 323–328 orbit/support packet compiler and splice it directly into the first concrete coherent local D/N kernel family from the moving-throat support-compensation program.

The main point is that the abstract branch packet is no longer floating over an unspecified support closure. On the concrete coherent local D/N family, the physical support ratio is explicit, the orbit packet stays completely support-blind, and the support theorem collapses to one exact regime test on the same branch.

## Stage 329 — concrete coherent local D/N branch state

Start from the microscopic coherent kernel family in which the mixed lane `W` and the support lane `phi_n` are sourced by the same local throat density. Then

\[
\lambda_W = \lambda_* I_0,
\qquad
\lambda_{\phi,n} = \lambda_* I_n,
\]
with the exact coherent placement-state map
\[
\chi_0 = \frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U = \frac{\pi^2 T_U}{L^2 K_U},
\qquad
Z_W = \frac{\lambda_W^2}{K_{\eta}^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\]
\[
\epsilon_W = \frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\qquad
\epsilon_\eta = \frac{c_{\eta U}^2}{K_U K_{\eta}^{(\mathrm{eff})}},
\qquad
\Lambda = \frac{27\pi^2 G c_s^5}{20 a^5 c^5}\frac{K_W^{(\mathrm{eff})}}{\mu_W}.
\]

On the same branch the physical coherent support ratio is no longer free. It is
\[
\zeta_n^{(\mathrm{phys})}
=
\frac{K_W^{(\mathrm{eff})}}{K_{\phi,n}^{(\mathrm{eff})}}
\left(\frac{I_n}{I_0}\right)^2.
\]
For the first uniform local source density,
\[
\frac{I_n}{I_0} = \frac{1}{2n+1},
\]
so
\[
\zeta_n^{(\mathrm{uniform})}
=
\frac{K_W^{(\mathrm{eff})}}{K_{\phi,n}^{(\mathrm{eff})}}\frac{1}{(2n+1)^2}.
\]
On the same-operator twin family,
\[
\zeta_n^{(\mathrm{twin})}=
\frac{1}{(2n+1)^2\bigl(1+x\,n(n+1)\bigr)}.
\]
The lowest symmetric twin lane is exact:
\[
\zeta_0^{(\mathrm{twin})}=1,
\qquad
S_0=2.
\]
So the lowest twin branch is the universal doubling branch of the concrete coherent support sector.

## Stage 330 — actual coherent local D/N orbit packet

On the actual coherent local D/N branch, the orbit packet is carried entirely by the three exact branch observables
\[
(R_{\rm tr},\,R_{\rm target},\,\epsilon_\eta),
\]
with
\[
R_{\rm tr}=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\qquad
R_{\rm target}=
\Lambda\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\epsilon=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right).
\]

The exact finite quotient packet is therefore
\[
q_{\rm tr} = -C_*\ln\frac{R_{\rm tr}}{R_{{\rm tr},{\rm ref}}},
\]
\[
q_{\rm nt} = B_*\ln\frac{R_{\rm tr}}{R_{{\rm tr},{\rm ref}}}
+ \ln\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}
- \ln\frac{R_{\rm target}}{R_{{\rm target},{\rm ref}}},
\]
\[
q_\eta = \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]
The infinitesimal defect packet is still
\[
\Theta_1=d\ln R_{\rm tr},
\qquad
\Xi_1=-d\ln R_{\rm target}-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}d\ln\epsilon_\eta,
\qquad
\mathcal R_1=d\ln R_{\rm target}.
\]

Most importantly, every one of these objects is exactly support-blind:
\[
\partial_\zeta q_{\rm tr}=
\partial_\zeta q_{\rm nt}=
\partial_\zeta q_\eta=
\partial_\zeta\Theta_1=
\partial_\zeta\Xi_1=
\partial_\zeta\mathcal R_1=0.
\]
So choosing the physical D/N support harmonic does not move the branch on or off the weak-axisymmetric similarity orbit.

## Stage 331 — actual coherent local D/N support threshold

On the same physical branch, the support theorem is controlled by the tracking-branch product
\[
\Pi_{\rm tr}(\xi,\delta;R)=F_{\rm tr}(\xi,\delta;R)\,G_{\rm tr}(\xi,\delta;R),
\]
and the mixed-only product scale
\[
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]
Then
\[
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}},
\qquad
\zeta_{\rm req}=
\frac{\Pi_{\rm tr}-C_{\rm mix}}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]

The exact regime split is:

1. mixed-only enough iff
   \[
   \Pi_{\rm tr} \le C_{\rm mix};
   \]
2. the symmetric lowest twin lane is enough iff
   \[
   C_{\rm mix} < \Pi_{\rm tr} \le 2C_{\rm mix};
   \]
3. non-twin asymmetry is required iff
   \[
   \Pi_{\rm tr} > 2C_{\rm mix}.
   \]

So on the concrete coherent local D/N family the support problem is no longer vague. The first physical support lane is the universal doubling branch, and it succeeds or fails according to one exact branch-product inequality.

Higher D/N harmonics are strongly suppressed. For the same-operator twin family,
\[
\zeta_n^{(\mathrm{twin})}
=
\frac{1}{(2n+1)^2(1+x n(n+1))}
<
\frac{1}{(2n+1)^2}
\qquad (n\ge 1),
\]
so an exact necessary condition for the `n`th twin harmonic to work is
\[
\zeta_{\rm req} \le \frac{1}{(2n+1)^2}.
\]
If that is satisfied, the exact softness threshold is
\[
x \le x_{\max}(n;\zeta_{\rm req})
=
\frac{1/((2n+1)^2\zeta_{\rm req})-1}{n(n+1)}.
\]

## Stage 332 — actual coherent local D/N branch theorem gate

Putting the two packets together, the actual coherent local D/N branch now ends at one exact theorem gate:

### Orbit packet
\[
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
\]
or equivalently
\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln\epsilon_\eta=0.
\]

### Support packet
With
\[
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\]
we have
- mixed-only enough iff `Π_tr <= C_mix`,
- lowest symmetric twin enough iff `C_mix < Π_tr <= 2 C_mix`,
- non-twin asymmetry required iff `Π_tr > 2 C_mix`.

On the lowest twin branch `zeta=1`,
\[
M_{\rm tr}=2M_{\rm mix},
\qquad
R_{\rm target}M_{\rm tr}=\frac{16\Lambda(1-\epsilon)}{\pi^2}.
\]

## Where we end up

The abstract packet compiler is now fully tied to the first concrete coherent local D/N operator family.

1. The actual weak-axisymmetric orbit-lock test is carried exactly by
   \[
   (R_{\rm tr},R_{\rm target},\epsilon_\eta),
   \]
   and it is rigorously blind to the support lane.
2. The actual coherent support theorem is carried exactly by
   \[
   \Pi_{\rm tr},\quad C_{\rm mix},\quad \zeta_n^{(\mathrm{phys})},
   \]
   with the lowest symmetric twin branch giving the universal doubling law.
3. So the remaining PDE-side question is not another algebraic compression. It is simply:

> does the completed moving-throat operator place the physical branch in the exact orbit-lock locus, and if so, is its support sector already in the mixed-only regime, the lowest-twin regime, or the genuinely non-twin regime?
# 5PN / Moving-Throat continuation — Stages 333–335

This session pushed the coherent D/N support theorem one step further down to an explicit regime classifier and then tied it to the natural minimal isotropic conservative quadrupole precursor.

The immediate continuation point before this session was:

> determine whether the actual passive/outgoing branch lands in the mixed-only, lowest-twin, or genuinely non-twin support regime.

The result is that this classification is now exact in one scalar ratio, and the canonical minimal isotropic branch lands comfortably in the symmetric-lowest-twin window.

---

## Stage 333 — exact branch-product support regimes

Let
\[
C_{\rm mix}:=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\qquad
S_{\rm req}:=\frac{\Pi_{\rm tr}}{C_{\rm mix}},
\]
and write the exact coherent support demand as
\[
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}
=
\frac{\Pi_{\rm tr}-C_{\rm mix}}
{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]

The exact thresholds are:

- at \(\Pi_{\rm tr}=C_{\rm mix}\),
  \[
  \zeta_{\rm req}=0;
  \]
- at \(\Pi_{\rm tr}=2C_{\rm mix}\),
  \[
  \zeta_{\rm req}=1.
  \]

The exact derivative is
\[
\frac{d\zeta_{\rm req}}{d\Pi_{\rm tr}}>0
\]
throughout the blocked branch \(0<\epsilon<1\), so the support demand increases strictly with the tracking-branch product.

Therefore the support problem splits into three exact regimes:

\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\Longrightarrow\quad
\text{mixed-only already enough},
\]

\[
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\Longrightarrow\quad
\text{symmetric lowest twin enough},
\]

\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\Longrightarrow\quad
\text{non-twin asymmetry required}.
\]

So the next branch decision no longer depends on many reduced parameters separately. It depends only on where the actual branch lands relative to the two product thresholds \(C_{\rm mix}\) and \(2C_{\rm mix}\).

---

## Stage 334 — exact loading-ratio compiler and canonical isotropic conservative precursor

On the explicit support/source branch, the natural contact-plus-pole conservative precursor may be written as
\[
Y_Q^{\rm cons}(\omega)
=
c_0+\frac{c_1}{1-\omega^2/\Omega_Q^2},
\qquad
c_0+c_1=1.
\]

Its exact support/source loading-ratio dictionary is
\[
c_0=\frac{1}{\rho_\alpha},
\qquad
c_1=\frac{\rho_\alpha-1}{\rho_\alpha},
\qquad
\rho_\alpha=\frac{1}{c_0}=\frac{1}{1-c_1},
\]
with
\[
\zeta_{\rm req}=\frac{c_1}{c_0}.
\]

The minimal isotropic conservative quadrupole module already fixed earlier is
\[
c_0=\frac34,
\qquad
c_1=\frac14.
\]
So it implies exactly
\[
\rho_\alpha^{(\min)}=\frac43,
\qquad
\zeta_{\rm req}^{(\min)}=\frac13.
\]

Since
\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\rho_\alpha,
\]
the same branch gives
\[
\Pi_{\rm tr}^{(\min)}=\frac43\,C_{\rm mix}.
\]

With blocking retained,
\[
\zeta_{\rm req}^{(\min)}(\epsilon)
=
\frac{\rho_\alpha^{(\min)}-1}{1-\epsilon(2-\rho_\alpha^{(\min)})}
=
\frac{1}{3-2\epsilon}.
\]

This is still strictly below \(1\) for every blocked branch with \(0\le \epsilon<1\). So the canonical isotropic passive/outgoing branch is always twin-safe before any non-twin asymmetry is needed.

---

## Stage 335 — explicit Family-1 verdict for the canonical isotropic branch

The explicit Family-1 reduced theorem window already frozen earlier is

\[
\rho_\alpha \le 3.46622291347846
\quad\Longrightarrow\quad
\text{guaranteed success},
\]
\[
\rho_\alpha \ge 3.46752913273870
\quad\Longrightarrow\quad
\text{guaranteed failure},
\]
with hard constructive ceiling
\[
\rho_\alpha < 3.46752922945601.
\]

Comparing the canonical isotropic value
\[
\rho_\alpha^{(\min)}=\frac43
\]
to that window gives very large exact margins:
\[
\Delta_{\rm suff}=3.46622291347846-\frac43\approx 2.13288958014513,
\]
\[
\Delta_{\rm fail}=3.46752913273870-\frac43\approx 2.13419579940537,
\]
\[
\Delta_{\rm max}=3.46752922945601-\frac43\approx 2.13419589612268.
\]

Likewise, on the support-ratio side
\[
\zeta_{\rm req}^{(\min)}=\frac13,
\qquad
\zeta_{\max}^{(F1)}\approx 2.46752922945601,
\]
so
\[
\zeta_{\max}^{(F1)}-\zeta_{\rm req}^{(\min)}
\approx 2.13419589612268.
\]

This branch also satisfies the explicit zero-bias Family-1 criterion, because the exact trigger is
\[
\zeta_{\rm req}<A_{F1},
\qquad
A_{F1}\approx 1.00005192880220,
\]
and
\[
\frac13 < 1.00005192880220.
\]
Therefore
\[
Pe_{\rm req}=0
\]
on the explicit Family-1 branch for the canonical isotropic passive/outgoing quadrupole module.

So the reduced support/source side is no longer the live bottleneck on that branch.

---

## Net result after Stages 333–335

The moving-throat support/source theorem has now been reduced to one exact scalar classifier:

\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\rho_\alpha.
\]

That one number tells us everything:

- \(\rho_\alpha\le 1\): mixed-only enough,
- \(1<\rho_\alpha\le 2\): symmetric lowest twin enough,
- \(\rho_\alpha>2\): non-twin asymmetry required.

The canonical minimal isotropic conservative quadrupole precursor gives
\[
\rho_\alpha=\frac43,
\]
so it lies strictly in the symmetric-lowest-twin regime and, on the explicit Family-1 branch, already succeeds at zero transport bias.

So the reduced theorem gap is now even sharper than before:

> the explicit support/source side is no longer the active uncertainty on the canonical isotropic branch.  
> The remaining reduced question is exactly what loading ratio \(\rho_\alpha\) the actual passive/outgoing moving-throat quadrupole branch selects.
# 5PN / Moving-Throat continuation — Stages 336–339

This pass takes the exact support-regime classifier from Stages 333–335 and ties it back into the physical selected-branch softening variable `xi` and the coherent-kernel placement map. That was the cleanest open step left on the support/source side.

The continuation point before this session was:

> express the actual support regime directly on the physical selected branch, instead of only in the abstract ratio `rho_alpha = Pi_tr / C_mix`.

The result is that the support/source side is now controlled by one exact scalar branch ratio
\[
\rho_\alpha^{(\mathrm{phys})}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{G_{\rm tr}}{M_{\rm mix}},
\]
and every regime boundary is an exact algebraic surface in the same selected-branch variable `xi`.

---

## Stage 336 — exact selected-branch loading ratio

On the tracking branch,
\[
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]
\[
F_{\rm tr}(\xi,\delta;R)
=
\frac{[9\delta+(9+2R^2)\xi]^2[9\delta+(9+2R)\xi]^2}
{81(1-\xi)\,[9\delta^2+18\delta\xi+(9+2R^2)\xi^2]^2},
\]
so the exact branch product is
\[
\Pi_{\rm tr}=F_{\rm tr}G_{\rm tr}.
\]

On the physical target branch, where
\[
R_{\rm target}=F_{\rm tr},
\qquad
C_{\rm mix}=R_{\rm target}M_{\rm mix},
\]
the loading ratio becomes exactly
\[
\boxed{
\rho_\alpha^{(\mathrm{phys})}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{G_{\rm tr}}{M_{\rm mix}}.
}
\]

So the support regime no longer needs the full product language once the branch is on-shell. It is just the ratio of the required tracking load to the mixed-only baseline.

The regime split is therefore

\[
G_{\rm tr}\le M_{\rm mix}
\quad\Longrightarrow\quad
\text{mixed-only enough},
\]

\[
M_{\rm mix}<G_{\rm tr}\le 2M_{\rm mix}
\quad\Longrightarrow\quad
\text{lowest symmetric twin enough},
\]

\[
G_{\rm tr}>2M_{\rm mix}
\quad\Longrightarrow\quad
\text{non-twin asymmetry required}.
\]

Because
\[
\frac{dG_{\rm tr}}{d\xi}>0,
\]
one also gets
\[
\frac{d\rho_\alpha^{(\mathrm{phys})}}{d\xi}>0.
\]
So the required support regime becomes strictly harder as the selected branch softens more deeply.

---

## Stage 337 — exact saturation depths `xi_(1x)` and `xi_(2x)`

The mixed-only threshold is the exact positive root of
\[
G_{\rm tr}=M_{\rm mix},
\]
namely
\[
\boxed{
\xi_{(1x)}
=
\frac{
M_{\rm mix}(9+2R^2)-9\delta
+
\sqrt{[M_{\rm mix}(9+2R^2)-9\delta]^2+324M_{\rm mix}\delta}
}{18}.
}
\]

The lowest-twin threshold is the exact positive root of
\[
G_{\rm tr}=2M_{\rm mix},
\]
namely
\[
\boxed{
\xi_{(2x)}
=
\frac{
2M_{\rm mix}(9+2R^2)-9\delta
+
\sqrt{[2M_{\rm mix}(9+2R^2)-9\delta]^2+648M_{\rm mix}\delta}
}{18}.
}
\]

Because `G_tr` is strictly increasing,
\[
\xi_{(1x)}<\xi_{(2x)}.
\]

So the support/source regime can now be read directly from the selected-branch depth:

- `xi_phys <= xi_(1x)` means mixed-only is enough,
- `xi_(1x) < xi_phys <= xi_(2x)` means the lowest symmetric twin is enough,
- `xi_phys > xi_(2x)` means non-twin asymmetry is required.

That is the cleanest direct classifier we have had so far on the support side.

---

## Stage 338 — coherent-kernel regime map

On the coherent local D/N branch,
\[
M_{\rm mix}
=
\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)}.
\]

Substituting this into the selected-branch loading ratio gives
\[
\boxed{
\rho_\alpha^{(\mathrm{coh})}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)}{8Z_W(1+\chi_0)^2}\,
G_{\rm tr}(\xi,\delta;R).
}
\]

So at fixed selected point, increasing `Z_W` lowers the required regime ratio and makes support success easier.

This gives two exact overlap thresholds:

### Mixed-only threshold
\[
\boxed{
Z_W^{(\mathrm{mix})}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)}{8(1+\chi_0)^2}
\,G_{\rm tr}(\xi,\delta;R).
}
\]

### Lowest-twin threshold
\[
\boxed{
Z_W^{(\mathrm{twin})}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)}{16(1+\chi_0)^2}
\,G_{\rm tr}(\xi,\delta;R).
}
\]

Equivalently, using
\[
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\]
the exact radiative-demand thresholds are

\[
\boxed{
\Lambda_{\rm mix}
=
\frac{\pi^2}{8(1-\epsilon)}\,
\Pi_{\rm tr},
}
\qquad
\boxed{
\Lambda_{\rm twin}
=
\frac{\pi^2}{16(1-\epsilon)}\,
\Pi_{\rm tr}.
}
\]

So the support theorem is now directly tied to the same coherent-kernel variables that place the actual moving-throat branch.

---

## Stage 339 — exact non-twin asymmetry requirement

The support-ratio demand is still
\[
\zeta_{\rm req}
=
\frac{\rho_\alpha-1}{1-\epsilon(2-\rho_\alpha)}.
\]

The exact excess over the symmetric-twin value is
\[
\boxed{
\zeta_{\rm req}-1
=
\frac{(1-\epsilon)(\rho_\alpha-2)}{1+\epsilon(\rho_\alpha-2)}.
}
\]

So the three support regimes become completely explicit:

- `rho_alpha <= 1`: mixed-only enough,
- `1 < rho_alpha <= 2`: lowest symmetric twin enough, with `0 < zeta_req <= 1`,
- `rho_alpha > 2`: non-twin asymmetry required, with `zeta_req > 1`.

On the physical selected branch,
\[
\rho_\alpha^{(\mathrm{phys})}=\frac{G_{\rm tr}}{M_{\rm mix}},
\]
so the non-twin condition is exactly
\[
G_{\rm tr}(\xi,\delta;R) > 2 M_{\rm mix},
\]
equivalently
\[
\xi_{\rm phys} > \xi_{(2x)}.
\]

That means the first true asymmetry requirement is no longer vague. It is the exact amount by which the selected branch has crossed above the twin threshold.

---

## Net result after Stages 336–339

The support/source side is now effectively solved at the reduced-theorem level.

1. The exact regime classifier has been pulled back to the physical selected-branch variable `xi`.
2. The mixed-only and lowest-twin thresholds are exact algebraic depths `xi_(1x)` and `xi_(2x)`.
3. The coherent-kernel placement map rewrites the classifier directly in terms of the actual branch variables `(chi_0, eps_eta, eps, Z_W)`.
4. The non-twin asymmetry requirement is now an exact excess formula rather than a qualitative statement.

So the remaining reduced question is no longer “what support regime might the branch be in?”
It is narrower:

> once the actual PDE gives the physical branch point `(xi_phys, delta, R, chi_0, eps_eta, eps, Z_W)`, the support regime and any required asymmetry are read off immediately from these exact formulas.  
> The only thing still open is where the completed moving-throat branch actually lands.

# 5PN continuation notes — Stages 340–342: reduced finish-line compression

This session does not introduce a new microscopic branch family. It compresses the
current theorem ledger to the narrowest reduced finish line reached so far.

The guiding idea is simple:

- the Stage-200 endgame compiler still carries an 8-component reduced branch packet
  together with the 3-component orbit packet;
- the later moving-throat branch work shows that, on the **actual isotropic**
  static-geometry one-pole contact-plus-pole branch, almost all of that structure
  has already collapsed.

So the right final move is to prove that collapse cleanly and to state the exact
remaining theorem gate in the smallest possible language.

---

## Stage 340 — the Stage-200 branch packet collapses to one scalar on the actual isotropic branch

The exact Stage-200 reduced branch packet is

\[
\Delta_{\rm branch}
=
(a_2,\ b_2,\ a_4,\ b_4,\ a(P_0),\ b(P_0),\ \Delta_{\rm pole},\ \Delta_{\rm norm}).
\]

On the later actual isotropic one-pole branch, the following hold simultaneously:

1. grouped-lane isotropy kills
   \[
   a_2=b_2=a_4=b_4=a(P_0)=b(P_0)=0,
   \]
2. the static-geometry one-pole conservative branch kills
   \[
   \Delta_{\rm pole}=0,
   \]
3. the only surviving reduced branch mismatch is the outgoing-normalization scalar
   \[
   N_Q-1.
   \]

So the full 8-component reduced branch packet collapses exactly to

\[
\Delta_{\rm branch}
\longrightarrow
(0,0,0,0,0,0,0,N_Q-1).
\]

Equivalently,
\[
\Delta_{\rm branch}=0
\iff
N_Q=1
\]
once the actual isotropic one-pole branch is accepted.

That is the cleanest bridge from the Stage-200 general endgame to the later
moving-throat finish-line statement.

---

## Stage 341 — the orbit packet is exactly a three-observable invariance test

The weak-axisymmetric orbit-lock side can now be written directly in the three
physical branch observables

\[
R_{\rm tr},\qquad R_{\rm target},\qquad \epsilon_\eta.
\]

At first order the defect triplet is

\[
\Theta_1 = \delta\ln R_{\rm tr},
\]
\[
\Xi_1 = -\,\delta\ln R_{\rm target}
        -\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\delta\ln \epsilon_\eta,
\]
\[
\mathcal R_1 = \delta\ln R_{\rm target}.
\]

And the inverse map is exact:

\[
\delta\ln R_{\rm tr} = \Theta_1,
\]
\[
\delta\ln R_{\rm target} = \mathcal R_1,
\]
\[
\delta\ln \epsilon_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}\,(\mathcal R_1+\Xi_1).
\]

So the full weak-axisymmetric orbit-lock problem is now exactly equivalent to the
three-observable invariance test

\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln R_{\rm tr}=0,\quad
\delta\ln R_{\rm target}=0,\quad
\delta\ln\epsilon_\eta=0.
\]

This is the direct-observable form of the older similarity-orbit / monomial-rigidity
theorem.

---

## Stage 342 — reduced completion theorem

There are now three exact pieces to combine.

### 1. Support/source theorem on the actual isotropic branch

On the actual isotropic branch,
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk})
=
\frac{1}{3-2\epsilon_{\rm blk}}.
\]

If a constructive explicit support/source family has ceiling `zeta_max > 1`, then
through the full admissible blocked window
\[
0\le \epsilon_{\rm blk}<\frac1{\zeta_{\max}},
\]
one has
\[
\zeta_{\rm req}^{(\rm act)} < \zeta_{\max}.
\]

Indeed, at the worst admissible edge
\[
\zeta_{\rm req,worst}
=
\frac{\zeta_{\max}}{3\zeta_{\max}-2},
\]
and
\[
\zeta_{\max}-\zeta_{\rm req,worst}
=
\frac{3\zeta_{\max}(\zeta_{\max}-1)}{3\zeta_{\max}-2}>0
\quad\text{for}\quad
\zeta_{\max}>1.
\]

So the explicit support/source side is automatic on the actual isotropic branch.
It is no longer an active reduced bottleneck.

### 2. Outgoing normalization on the natural source-map branch

On the natural source-map branch,
\[
\hat m_0\to 1,
\]
so the exact outgoing factorization reduces to
\[
N_Q=\frac{1}{\chi_Q}.
\]

Therefore
\[
N_Q=1
\iff
\chi_Q=1.
\]

So on the natural source-map branch the last outgoing defect is purely the scalar
outgoing renormalization `chi_Q - 1`.

### 3. Final reduced finish line

Inside the present reduced hierarchy, once the actual isotropic conservative precursor
is accepted, the full reduced 5PN / 2.5PN / 4PN closure problem has collapsed to

\[
\delta\ln R_{\rm tr}=0,
\qquad
\delta\ln R_{\rm target}=0,
\qquad
\delta\ln\epsilon_\eta=0,
\qquad
N_Q=1,
\]
or, on the natural source-map branch,
\[
\delta\ln R_{\rm tr}=0,
\qquad
\delta\ln R_{\rm target}=0,
\qquad
\delta\ln\epsilon_\eta=0,
\qquad
\chi_Q=1.
\]

So the reduced theorem gate is no longer:

- a large grouped branch packet,
- a separate support/source branch phase,
- and a diffuse outgoing normalization problem.

It is now exactly:

1. a **three-observable orbit-lock test**, and
2. a **one-scalar outgoing-normalization test**.

That is the narrowest honest reduced finish line reached so far.

---

## Best current verdict

Within the present reduced hierarchy:

- the explicit Family-1 support/source side is no longer an independent obstacle,
- the Stage-200 branch packet has collapsed to `N_Q - 1` on the actual isotropic branch,
- the weak-axisymmetric side is exactly the invariance of
  \[
  (R_{\rm tr},R_{\rm target},\epsilon_\eta),
  \]
- and the only remaining reduced theorem gap is whether the completed moving-throat PDE
  actually realizes

\[
\delta\ln R_{\rm tr}
=
\delta\ln R_{\rm target}
=
\delta\ln\epsilon_\eta
=
0,
\qquad
N_Q=1
\]
(or equivalently `chi_Q = 1` on the natural source-map branch).

That is as close to “finished” as the reduced program can honestly get without the
final PDE branch realization.
# 5PN continuation — Stages 343–345: landing on the four reduced finish-line conditions

This session takes the reduced finish line from Stages 340–342 and asks the next obvious question:

> can the current reduced moving-throat hierarchy actually land on the four surviving conditions simultaneously, or is there still a hidden algebraic incompatibility?

The answer is better than expected.

Within the present reduced hierarchy there is **no algebraic contradiction** among the four surviving conditions. The whole problem splits cleanly into:

1. a three-equation **orbit-lock landing surface** for
   - `delta ln R_tr = 0`,
   - `delta ln R_target = 0`,
   - `delta ln epsilon_eta = 0`,
2. and a one-equation **outgoing landing surface** for
   - `N_Q = 1`, equivalently `chi_Q = 1` on the natural source-map branch.

So the remaining obstruction is now purely PDE-side branch realization, not another reduced-sector inconsistency.

---

## Stage 343 — exact orbit-lock landing surface

On the coherent branch,

a) the direct tracking observable is

a) `R_tr = (1 + chi_0/(1+delta_U)) / (1 + chi_0)`,

b) the target observable is

b) `R_target = Lambda (1-epsilon_eta)(1-epsilon)^2 / [ Z_W (1+chi_0)^2 ]`,

c) and the third direct observable is simply

c) `epsilon_eta`.

The exact first-order landing conditions are therefore

1. `delta ln R_tr = 0`,
2. `delta ln R_target = 0`,
3. `delta ln epsilon_eta = 0`.

The script shows these are exactly equivalent to the solved co-drift system

`dln delta_U = - (1 + delta_U)/(1 + chi_0) * dln chi_0`,

`dln Lambda = dln Z_W + 2 chi_0 dln chi_0/(1 + chi_0) + 2 epsilon dln epsilon/(1 - epsilon)`,

`dln epsilon_eta = 0`.

So the first three finish-line conditions live on one explicit orbit-lock surface in the coherent branch variables.

A useful structural fact is that this whole orbit packet is support-blind at the reduced level:

`d_zeta ln R_tr = d_zeta ln R_target = d_zeta ln epsilon_eta = 0`.

So support enhancement still matters for the isotropic normalization baseline, but it does **not** move the branch on or off the weak-axisymmetric orbit-lock surface.

---

## Stage 344 — exact outgoing landing surface

On the natural source-map branch,

`N_Q = 1 / chi_Q`.

So the fourth finish-line condition is exactly

`N_Q = 1  <=>  chi_Q = 1`.

That already shows one exact outgoing landing surface: the canonical outgoing branch with `chi_Q = 1`.

The session also rewrote the first-order Family-1 outgoing defect in the form

`Delta_Q = - [sigma_* / (1 - sigma_*)] Xi_slip deltaPi_tan`,

so that

`N_Q - 1 = 1/(1 + Delta_Q) - 1`.

Therefore, on the exact lower parent compensation family, the first-order outgoing defect vanishes automatically whenever

`Xi_slip = 0`.

Then

`Delta_Q = 0`,

`N_Q - 1 = 0`.

So there are two exact reduced ways to land the fourth condition:

1. the canonical outgoing branch `chi_Q = 1`,
2. or, at first order, the lower parent compensation family with `Xi_slip = 0`.

---

## Stage 345 — the combined four-condition landing theorem

Putting the two pieces together gives the current cleanest reduced theorem:

### Orbit-lock surface

`dln delta_U = - (1 + delta_U)/(1 + chi_0) * dln chi_0`,

`dln Lambda = dln Z_W + 2 chi_0 dln chi_0/(1 + chi_0) + 2 epsilon dln epsilon/(1 - epsilon)`,

`dln epsilon_eta = 0`.

### Outgoing surface

Either

`chi_Q = 1`,

or, on the first-order lower compensation family,

`Xi_slip = 0`.

### Combined reduced finish line

If the branch lies on the orbit-lock surface and also lies on either outgoing surface above, then all four reduced finish-line conditions vanish simultaneously:

`delta ln R_tr = 0`,

`delta ln R_target = 0`,

`delta ln epsilon_eta = 0`,

`N_Q - 1 = 0`.

So the four-condition finish line is now shown to be **algebraically reachable** inside the current reduced hierarchy.

---

## Best current verdict

This session removes the last serious worry that there might still be some hidden reduced-sector contradiction among the four finish-line conditions.

There is not.

What remains open is exactly what the compact moving-throat master has been saying:

- whether the completed PDE actually realizes the orbit-lock surface,
- and whether it actually realizes the canonical / lower-compensation outgoing surface.

So the open problem is now purely one of **branch realization**, not another reduced-algebra obstruction.
# 5PN continuation notes — Stages 346–349

These stages stop trying to compress the reduced theorem any further and instead extract the **actual coherent-branch finish packet** as far as the present file stack allows.

The goal was the one left at Stage 345:

\[
(R_{\rm tr},\ R_{\rm target},\ \epsilon_\eta,\ N_Q)
\]

for the actual coherent local D/N branch, and then the exact surfaces on which the four reduced finish-line conditionals vanish.

The main outcome is sharp:

1. the current notes/scripts do fix the **actual coherent-branch formulas** for
   \[
   R_{\rm tr},\quad R_{\rm target},\quad \epsilon_\eta,
   \]
   and the natural-source-map relation
   \[
   N_Q=\chi_Q^{-1},
   \]
2. the first three finish-line conditions reduce to exact codrift surfaces in the physical branch variables,
3. the outgoing condition is the exact canonical branch condition
   \[
   \chi_Q=1
   \quad\Longleftrightarrow\quad
   N_Q=1,
   \]
   with the lower-parent compensation family giving the linear sufficient condition
   \[
   \Xi_{\rm slip}\,\delta\Pi_{\rm tan}=0,
   \]
4. but the file stack still does **not** contain a completed numerical PDE-selected branch point.  
   It contains the exact symbolic branch packet and the exact landing surfaces.

So the reduced-theorem side is finished, but the PDE-realization side is still open.

---

## Stage 346 — actual coherent local D/N branch values

On the actual coherent tracking branch,

\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{\chi_0+\delta_U+1}{(1+\chi_0)(1+\delta_U)},
\]

\[
\epsilon
=
\epsilon_W^{(\mathrm{split})}
=
\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)
=
\epsilon_W\frac{11+9\delta_U}{11(1+\delta_U)},
\]

\[
R_{\rm target}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]

\[
M_{\rm mix}
=
\frac{8 Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
\qquad
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\]

\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]

and on the natural source-map branch,

\[
N_Q=\frac{1}{\chi_Q}.
\]

The exact support-blindness check is explicit:
\[
\partial_\zeta R_{\rm tr}=0,
\qquad
\partial_\zeta \epsilon_\eta=0,
\qquad
\partial_\zeta R_{\rm target}=0.
\]

So the coherent support lane changes only the baseline through \(M_{\rm tr}=M_{\rm mix}S\); it does not move the orbit packet.

There is also an exact product law:
\[
R_{\rm target}M_{\rm mix}
=
\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

---

## Stage 347 — the first three finish-line conditions in physical branch variables

Using logarithmic branch drifts
\[
d\ln\chi_0,\quad d\ln\delta_U,\quad d\ln\epsilon_W,\quad d\ln\epsilon_\eta,\quad d\ln Z_W,\quad d\ln\Lambda,
\]
the exact split-blocking drift is
\[
d\ln\epsilon
=
d\ln\epsilon_W
-
\frac{2\delta_U}{(1+\delta_U)(11+9\delta_U)}\,d\ln\delta_U.
\]

The three reduced finish-line conditionals are

\[
d\ln R_{\rm tr}
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}
\Big[(1+\delta_U)d\ln\chi_0+(1+\chi_0)d\ln\delta_U\Big],
\]

\[
d\ln R_{\rm target}
=
d\ln\Lambda-d\ln Z_W
-\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
-\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
-\frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
\]

\[
d\ln\epsilon_\eta=d\ln\epsilon_\eta.
\]

So the exact orbit-lock / target-ratio / dressing landing surfaces are

\[
(1+\delta_U)d\ln\chi_0+(1+\chi_0)d\ln\delta_U=0,
\]

\[
d\ln\Lambda-d\ln Z_W
-\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
-\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
-\frac{2\epsilon}{1-\epsilon}d\ln\epsilon=0,
\]

\[
d\ln\epsilon_\eta=0.
\]

After imposing the first and third, the second becomes the actual coherent-branch codrift law
\[
d\ln\Lambda
=
d\ln Z_W
+\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
+\frac{2\epsilon}{1-\epsilon}d\ln\epsilon.
\]

---

## Stage 348 — the outgoing normalization finish surface

The exact outgoing condition on the natural source-map branch is
\[
N_Q=\chi_Q^{-1},
\]
so the fourth finish-line condition is simply
\[
N_Q=1
\quad\Longleftrightarrow\quad
\chi_Q=1.
\]

On the canonical compact passive/outgoing branch this is exactly the Stage-87/89 condition.

On the lower-parent compensation family at first order,
\[
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan},
\qquad
N_Q=\frac{1}{1+\Delta_Q},
\]
so either
\[
\Xi_{\rm slip}=0
\qquad\text{or}\qquad
\delta\Pi_{\rm tan}=0
\]
is sufficient to land on
\[
N_Q=1.
\]

So the outgoing side is not another orbit condition. It is a separate exact surface.

---

## Stage 349 — actual four-condition extractor and honest verdict

The actual coherent branch packet is therefore

\[
R_{\rm tr}
=
\frac{\chi_0+\delta_U+1}{(1+\chi_0)(1+\delta_U)},
\qquad
R_{\rm target}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]

\[
\epsilon_\eta=\epsilon_\eta,
\qquad
N_Q=\chi_Q^{-1}.
\]

The exact four finish-line conditionals are

1. \(d\ln R_{\rm tr}=0\),
2. \(d\ln R_{\rm target}=0\),
3. \(d\ln\epsilon_\eta=0\),
4. \(N_Q-1=0\).

Their exact combined landing compiler is:

- orbit lock
  \[
  d\ln\delta_U
  =
  -\frac{1+\delta_U}{1+\chi_0}\,d\ln\chi_0,
  \]
- dressing lock
  \[
  d\ln\epsilon_\eta=0,
  \]
- target-ratio lock
  \[
  d\ln\Lambda
  =
  d\ln Z_W
  +\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
  +\frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
  \]
- outgoing lock
  \[
  \chi_Q=1
  \quad\text{(equivalently \(N_Q=1\)).}
  \]

This is the sharpest actual-branch extraction the present file stack supports.

### Honest end-state

What the current files **do** contain:

- the actual coherent-branch symbolic packet,
- the exact support-blind orbit packet,
- the exact outgoing normalization surface,
- and the exact combined four-condition landing surfaces.

What they **do not** yet contain:

- a completed numerical or closed-form PDE-selected point
  \[
  (R_{\rm tr}^{\rm phys},R_{\rm target}^{\rm phys},\epsilon_\eta^{\rm phys},N_Q^{\rm phys}).
  \]

So the reduced program has reached the point where **all remaining uncertainty is branch realization**.

That means the next real move is no longer another reduced compiler.
It is to solve or numerically locate the completed moving-throat branch strongly enough to read off those four actual values.
# 5PN Stage 350–351 — Exact Family-1 branch location on the refreshed `\Lambda_{EM}` geometry

This pass stops compressing the support/source side and actually solves the exact operator-selected support/source fixed-point equation on the first explicit Family-1 branch, using the **refreshed**
\[
\Lambda_{EM}=\frac{\sqrt2\,\pi}{x_{01}}
\]
geometry rather than the old `37/20` shorthand.

## What was numerically located

Using
\[
\Lambda_\ell = \frac{L}{\ell}=20\Lambda_{EM}\approx 36.94973154240256,
\qquad
\eta=\Lambda_\ell,
\qquad
\kappa=\frac{1440\pi^2}{x_{01}^2}\approx 2457.5087899001137,
\]
the exact Robin/support constants are

\[
y\tan y=\eta
\quad\Longrightarrow\quad
y\approx 1.5294278190457656,
\]
\[
A_K=\frac{\kappa+\pi^2/4}{\kappa+y^2}\approx 1.0000521380385143,
\]
\[
\Delta_0\approx 1.7377393923469950\times 10^{-4},
\qquad
\Delta_\infty\approx 2.0172162594593645\times 10^{-2},
\]
\[
\zeta_{\max}=A_K\frac{\pi^2}{4}\approx 2.4675297457259358.
\]

Then, using the two explicit wall-depth extractions already carried in the notes,
\[
\Theta_w^{(\chi)}\approx 4.06863235008162\,\lambda_\mu^2,
\qquad
\Theta_w^{(J)}\approx 0.927552032539308\,\lambda_\mu^2,
\]
and setting the benchmark \(\lambda_\mu=1\), the exact wall/source figure of merit is
\[
\Xi = W_{\rm wall}=100\,\Theta_w\,\Lambda_\ell^2.
\]

This gives:

### χ-weighted extraction
\[
\Xi_\chi \approx 5.5548332017764099\times 10^5,
\]
and the exact fixed-point equation
\[
Pe=\Xi_\chi\,\Delta(Pe;\kappa,\eta)
\]
has the numerically located branch point
\[
Pe_*^{(\chi)} \approx 11155.7265863205869.
\]

At that root,
\[
\zeta_{\rm phys}^{(\chi)} \approx 2.4675296478814376,
\qquad
\rho_{\alpha,\max}^{(\chi)} = 1+\zeta_{\rm phys}^{(\chi)} \approx 3.4675296478814376.
\]

### J-weighted extraction
\[
\Xi_J \approx 1.2663707072528143\times 10^5,
\]
and the exact fixed-point branch point is
\[
Pe_*^{(J)} \approx 2504.9703142859238.
\]

At that root,
\[
\zeta_{\rm phys}^{(J)} \approx 2.4675278051675084,
\qquad
\rho_{\alpha,\max}^{(J)} = 1+\zeta_{\rm phys}^{(J)} \approx 3.4675278051675084.
\]

So both explicit operator-selected Family-1 branches sit extremely close to the exact Family-1 ceiling
\[
\zeta_{\max}\approx 2.4675297457259358.
\]

## What this means for the actual isotropic branch

The natural isotropic passive/outgoing grouped-`P2` branch still requires only
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13
\]
in the unblocked case, and it stays lowest-twin-safe throughout the admissible blocked interval.

Comparing against the numerically located Family-1 roots gives very large margins:

### χ branch
\[
\zeta_{\rm phys}^{(\chi)}-\zeta_{\rm req} \approx 2.1341963145481043,
\]
\[
\rho_{\alpha,\max}^{(\chi)}-\frac43 \approx 2.1341963145481043.
\]

### J branch
\[
\zeta_{\rm phys}^{(J)}-\zeta_{\rm req} \approx 2.1341944718341751,
\]
\[
\rho_{\alpha,\max}^{(J)}-\frac43 \approx 2.1341944718341751.
\]

So the support/source side is no longer just “probably okay.” On the refreshed exact `\Lambda_{EM}` Family-1 branch it is numerically located and **deep inside** the lowest-symmetric-twin-safe region.

## Finish-line status after the numerical location

This pass closes the strongest remaining **numerically accessible** part of the problem from the current stack:

- the explicit support/source branch has now been located numerically;
- it is strongly non-bottlenecked;
- on the canonical compact passive/outgoing branch one still has
  \[
  \chi_Q=1,\qquad N_Q=1
  \]
  exactly.

What is **still not numerically present in the files** is a completed PDE-selected point for the orbit-lock packet:
\[
d\ln R_{\rm tr}=0,\qquad d\ln R_{\rm target}=0,\qquad d\ln\epsilon_\eta=0.
\]

So the honest finish-line verdict is now:

1. **support/source** — numerically located and safe;
2. **outgoing normalization** — exactly fixed on the canonical passive/outgoing branch;
3. **orbit lock** — still only given as exact drift surfaces, not as a numerically selected PDE point.

That is the sharpest landing available from the current notes and scripts without a stronger direct solve of the completed moving-throat operator.
# Stage 352–353 — Final Numerical Object Assessment

## Goal

Decide whether the current file stack already contains the **last missing numerical object**: the PDE-selected orbit-lock point, rather than only its exact theorem surface.

## What is already numerically located

The explicit Family-1 support/source branch **is** numerically located on the refreshed exact `\Lambda_{EM}` geometry branch. The current reduced operator family gives

- `\Lambda_\ell = L/\ell \approx 36.9497315424`
- `\kappa \approx 2457.5087899001`
- `\zeta_{\max} \approx 2.4675297457`
- on the `\chi`-weighted wall-depth extraction,
  - `Pe_* \approx 11155.7265863`
  - `\zeta_{\rm phys}(Pe_*) \approx 2.4675296479`
  - `\rho_{\alpha,\max} \approx 3.4675296479`

So the explicit support/source side is not only algebraically under control; it is numerically located and safely above the canonical isotropic demand.

## What the current files still do **not** numerically locate

The coherent-branch finish packet is still only present **symbolically** in the current scripts/notes. The exact actual-branch quantities are already compiled as

- `R_tr = (\chi_0 + \delta_U + 1) / ((\chi_0+1)(\delta_U+1))`
- `\epsilon = \epsilon_W (9\delta_U+11)/(11(\delta_U+1))`
- `R_target = \Lambda (1-\epsilon_\eta)(1-\epsilon)^2 / (Z_W(1+\chi_0)^2)`
- `N_Q = 1/\chi_Q` on the natural source-map branch

but these remain formulas in the branch variables, not a numerically selected branch point.

The decisive script verdict from the current stack is still:

> *“No numerical PDE-selected point is present yet in the notes/scripts.”*

So the final unlocated object is specifically the **orbit-lock / coherent placement point** — equivalently the realized values of `(R_tr, R_target, \epsilon_\eta)` or the realized drift verdict `(\delta\ln R_tr, \delta\ln R_target, \delta\ln \epsilon_\eta)` on the completed moving-throat branch.

## Bottom line

We have already numerically located the explicit support/source branch, but **we have not numerically located the final PDE-selected orbit-lock point** from the present file stack alone.

So the remaining gap is not more reduced algebra. It is the actual numerical or analytic solve of the completed moving-throat branch strongly enough to read off the coherent placement variables.

## Deliverables in this turn

- `5pn_stage352_final_object_status.py`
- `5pn_stage352_final_object_status_output.txt`
- `5pn_stage352_353_final_object_assessment.md`
# 5PN computational handoff — what still has to be done to finish the program

## 0. Executive status

The reduced theorem side is now essentially closed.

The current chain has already reduced the problem to four finish-line conditions on the actual PDE-selected branch:

1. `d ln R_tr = 0`
2. `d ln R_target = 0`
3. `d ln epsilon_eta = 0`
4. `N_Q = 1`  (equivalently `chi_Q = 1` on the natural source-map branch)

The explicit Family-1 support/source side has already been located on the refreshed exact `Lambda_EM` geometry and sits safely inside the lowest-symmetric-twin regime. So the support/source side is no longer the active bottleneck in the reduced hierarchy.

What is still missing is the **actual PDE-selected orbit-lock point** on the coherent local D/N branch, plus confirmation that the realized outgoing branch is the canonical passive/outgoing one rather than a nearby deformation.

So the remaining computational job is not “more reduced algebra.” It is:

- solve or continue the completed moving-throat branch strongly enough to extract the actual coherent placement state,
- compute the weak-axisymmetric tangent of that branch,
- evaluate the orbit packet and outgoing normalization on the same branch,
- and return a final four-condition verdict.

---

## 1. Non-negotiable computational firewall

These rules are part of the handoff. They should not be changed casually.

### 1.1 Keep the parent bulk PDE fixed

Do **not** retune the parent GNLS / bulk medium to force 5PN.

Keep:

- the `4+1` parent field theory,
- the `P = K rho^5` medium,
- the corrected charge / ontology firewall,
- and the moving-throat geometry lift.

The remaining gap is branch realization, not a bulk-theory refit.

### 1.2 Keep the exact geometry relation

Do **not** revert to the old `L/a = 37/20` shorthand.

Use the carried exact relation

`Lambda_EM = sqrt(2) pi / x_01`

and any downstream geometry derived from it.

### 1.3 Keep support/source and orbit-lock separate

The support ratio `zeta` and the support enhancement factor `S(zeta; epsilon)` affect the isotropic normalization / baseline lane.

They do **not** enter the coherent weak-axisymmetric orbit packet.

So the code path must keep two packets separate:

- **support / normalization packet**
- **orbit / similarity packet**

---

## 2. Exact remaining data packet to extract from the PDE branch

There are really two extraction layers.

### 2.1 Stationary isotropic coherent placement state

The first required packet is the actual stationary coherent-branch state

`(chi0, delta_U, Z_W, epsilon_W, epsilon_eta, Lambda, zeta)`

with derived quantities

`epsilon = epsilon_W * (1 - (2/11) * delta_U/(1 + delta_U))`

`R_tr = (1 + chi0/(1 + delta_U)) / (1 + chi0)`

`R_target = Lambda * (1 - epsilon_eta) * (1 - epsilon)^2 / ( Z_W * (1 + chi0)^2 )`

`M_mix = 8 * Z_W * (1 + chi0)^2 / ( pi^2 * (1 - epsilon_eta) * (1 - epsilon) )`

`M_tr = M_mix * S(zeta; epsilon)`

and, on the outgoing side,

- `chi_Q` and/or
- `N_Q`

with the natural-source-map relation `N_Q = 1/chi_Q` if that branch is used.

### 2.2 Weak-axisymmetric branch tangent / drift packet

The second required packet is the actual first grouped weak-axisymmetric drift of the same branch.

The clean reduced form is

`(d ln chi0, d ln delta_U, d ln Z_W, d ln epsilon_W, d ln epsilon_eta, d ln Lambda)`

Optionally, the microscopic eight-drift packet is also acceptable:

`(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1)`

because that can be pushed through the Stage-320 / Stage-321 / Stage-323 machinery.

### 2.3 Optional but valuable regime-classification packet

For the support regime classifier, it is also useful to extract one of the following:

- `Pi_tr`
- `rho_alpha = Pi_tr / C_mix`
- `xi_phys` together with the selected-branch parameters used in the coherent support notes

This is lower priority than the orbit packet, because the explicit Family-1 support/source side is already known to be non-bottlenecked on the canonical isotropic branch.

---

## 3. Work package sequence

## WP0 — Freeze conventions and branch choice

Before any solve:

1. fix the exact `Lambda_EM` geometry,
2. fix the sign conventions for the outgoing DtN branch,
3. fix the natural source-map convention,
4. fix the coherent local D/N branch definitions,
5. fix whether `chi_Q` or `N_Q` is the primary outgoing observable in the solver output.

Deliverable:

- one short conventions file or JSON manifest used by every downstream script.

## WP1 — Solve / continue the stationary isotropic moving-throat branch

Goal:

Compute the actual stationary isotropic branch state from the completed operator, not from the reduced prototype.

Minimum acceptable output:

`chi0, delta_U, Z_W, epsilon_W, epsilon_eta, Lambda, zeta, chi_Q or N_Q`

Recommended method:

- finite throat in axial coordinate,
- wall profile + support + localized Maxwell/mixed sectors retained,
- continuation from the explicit Family-1 / D/N reduced branch,
- Newton or Newton-Krylov solve on the stationary nonlinear system,
- pseudo-arclength continuation if the branch folds.

Numerical tasks:

1. choose the stationary throat profile family / discretization,
2. solve the isotropic stationary branch,
3. postprocess the solution into the coherent placement variables.

Stop condition:

- stationary residuals small enough that the placement variables are stable under refinement.

## WP2 — Extract isotropic support / normalization packet on the actual branch

Goal:

Evaluate the isotropic branch packet on the realized branch.

Compute:

- `epsilon`
- `R_tr`
- `R_target`
- `M_mix`
- `M_tr`
- `C_mix = 8 Lambda (1 - epsilon) / pi^2`
- `rho_alpha = Pi_tr / C_mix` if `Pi_tr` is available
- `chi_Q` and `N_Q`

Check:

1. support regime:
   - mixed-only,
   - lowest symmetric twin,
   - non-twin.
2. outgoing normalization:
   - `N_Q - 1`, or equivalently `chi_Q - 1`.

This is where the explicit Family-1 branch-location results become a benchmark, not the final answer.

## WP3 — Solve the first weak-axisymmetric tangent problem on that branch

Goal:

Compute the grouped weak-axisymmetric tangent of the **actual** realized stationary branch.

Minimum acceptable output:

`d ln chi0, d ln delta_U, d ln Z_W, d ln epsilon_W, d ln epsilon_eta, d ln Lambda`

Recommended method:

- linearize the completed stationary operator in the grouped real `P2` weak-axisymmetric sector,
- solve the tangent / response problem at zero or low frequency,
- extract the coherent log drifts by projection.

Optional richer output:

- microscopic drift vector `(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1)`.

## WP4 — Evaluate the orbit-lock packet on the realized branch

Goal:

Test whether the actual branch is tangent to the exact similarity orbit.

Compute either the finite quotient packet or its linearization:

### Finite packet

`q_tr, q_nt, q_eta`

or

`Q_tr = ln( C_tr / C_tr,ref )`

`Q_nt = ln( C_nt / C_nt,ref )`

`Q_eta = ln( epsilon_eta / epsilon_eta,ref )`

### Linear packet

`Theta_1, Xi_1, R_1`

with the direct observable formulas

`Theta_1 = d ln R_tr`

`Xi_1 = - d ln R_target - (epsilon_eta / (1 - epsilon_eta)) d ln epsilon_eta`

`R_1 = d ln R_target`

and inverse relations if needed.

Success condition:

`d ln R_tr = 0`

`d ln R_target = 0`

`d ln epsilon_eta = 0`

or equivalently

`q_tr = q_nt = q_eta = 0`.

## WP5 — Final four-condition verdict

Goal:

Issue the final reduced completion verdict on the actual realized branch.

The final conditions are

1. `d ln R_tr = 0`
2. `d ln R_target = 0`
3. `d ln epsilon_eta = 0`
4. `N_Q = 1`  (equivalently `chi_Q = 1` on the natural source-map branch)

Decision tree:

- if 1–3 fail and 4 passes: outgoing normalization is fine, orbit lock fails;
- if 1–3 pass and 4 fails: coherent branch is on-orbit, outgoing normalization fails;
- if support/source fails its regime test but 1–4 pass: check convention mismatch, because that would be inconsistent with the reduced theorem chain;
- only if all four pass do we have the reduced 5PN / 2.5PN / 4PN closure on the realized branch.

---

## 4. Recommended solver / post-processing architecture

### 4.1 Solver side

Use one code path for the stationary isotropic solve and one code path for the weak-axisymmetric tangent solve.

Suggested split:

- `solver_stationary.py` or equivalent:
  solves the isotropic stationary branch.

- `solver_tangent_p2.py` or equivalent:
  linearizes around the stored stationary branch and solves the grouped weak-axisymmetric tangent.

### 4.2 Post-processing side

Re-use the reduced compilers already built in the session chain instead of rewriting formulas by hand.

The most useful existing script anchors are:

- `5pn_stage201_home_stretch_theorem.py`
- `5pn_stage318_actual_branch_orbit_tester_api.py`
- `5pn_stage325_coherent_branch_two_packet_compiler.py`
- `5pn_stage332_actual_branch_theorem_gate.py`
- `5pn_stage339_nontwin_asymmetry_requirement.py`
- `5pn_stage349_actual_four_condition_extractor.py`
- `5pn_stage350_family1_exact_branch_locator.py`

These are best treated as validators / packet compilers once real branch data exist.

### 4.3 Suggested data format

Save one isotropic branch JSON and one tangent JSON.

Recommended isotropic JSON fields:

```json
{
  "chi0": 0.0,
  "delta_U": 0.0,
  "Z_W": 0.0,
  "epsilon_W": 0.0,
  "epsilon_eta": 0.0,
  "Lambda": 0.0,
  "zeta": 0.0,
  "chi_Q": 0.0,
  "N_Q": 0.0,
  "Pi_tr": null,
  "mhat0": 1.0
}
```

Recommended tangent JSON fields:

```json
{
  "dln_chi0": 0.0,
  "dln_delta_U": 0.0,
  "dln_Z_W": 0.0,
  "dln_epsilon_W": 0.0,
  "dln_epsilon_eta": 0.0,
  "dln_Lambda": 0.0
}
```

If the microscopic drift packet is easier to extract, also save it.

---

## 5. Practical numerical warnings

### 5.1 Do not force canonical outgoing normalization

If the realized PDE branch gives `chi_Q != 1`, record it.

Do **not** manually project it back onto the canonical outgoing branch just to make the theorem pass.

### 5.2 Do not mix support and orbit packets

Because `zeta` drops out of the orbit packet exactly, a change in support enhancement should never be used to “explain” an orbit-lock miss.

If the orbit packet misses, it is a real branch-realization miss.

### 5.3 Track refinement sensitivity

For every extracted quantity, compare at least two discretizations / resolutions.

The handoff is not complete unless the final packet is numerically stable under refinement.

### 5.4 Save the full branch state before reduction

Do not only save reduced numbers. Save the actual solved stationary profiles / operator data too, so later work can revisit the reduction if the verdict is surprising.

---

## 6. Definition of done

The computational side is complete when all of the following exist:

1. one converged stationary isotropic branch solution,
2. one converged grouped weak-axisymmetric tangent about that same branch,
3. extracted coherent placement state,
4. extracted orbit packet,
5. extracted outgoing normalization scalar,
6. one final four-condition verdict.

At that point the reduced program is either:

- **closed on the realized PDE branch**, or
- **falsified at one explicit condition**.

Either outcome is useful.
