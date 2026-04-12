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
