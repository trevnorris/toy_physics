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
