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
