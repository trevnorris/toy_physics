
# Moving-Throat PDE — Stage 116: Finite D/N Mixed-Tube Realization

## Goal

Give the concrete core-balance theorem a geometric throat-core realization rather than leaving the bare mixed-channel data `(\kappa_0,\gamma_0)` abstract.

## Bare mixed D/N tube

Take the mixed side-channel to live on a finite auxiliary tube of length `L_W` with the first D/N half-wave:
\[
q''+k^2 q=0,
\qquad
q(0)=0,
\qquad
q'(L_W)=0.
\]
Then
\[
\boxed{
k_W=\frac{\pi}{2L_W},
\qquad
\Omega_W=\frac{\pi c_s}{2L_W}.
}
\]
Writing the side-channel pole denominator in the outlet variable
\[
z=\frac{a\omega}{c_s},
\]
the bare even coefficient is
\[
\boxed{
\kappa_0=\frac{\omega^2}{\Omega_W^2}\Big/\! z^2
=\frac{4L_W^2}{\pi^2 a^2}.
}
\]

## Compensation-selected tube length

Stage 115 requires
\[
\kappa_0=\frac{1+r_c}{3},
\qquad
r_c=\frac{\lambda^2}{K_sK_q}.
\]
So the auxiliary mixed-tube length is fixed to
\[
\boxed{
L_W
=
\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
}
\]
Thus the static/mixed hybridization `r_c` has a direct geometric meaning: it sets the ratio between the auxiliary mixed-tube length and the exterior mouth radius.

## Bare outgoing normalization

The same compensation theorem requires
\[
\gamma_0=\frac{1+r_c}{9}.
\]
A simple concrete realization is to take the bare mixed outlet to be a pure-scale deformation of the canonical compact outgoing `l=2` branch:
\[
\boxed{
D_W^{\rm bare}(z)
=
(1+r_c)\left(1-\frac{z^2}{3}-i\frac{z^5}{9}\right)+O(z^6).
}
\]
Then the hybridization factor `(1+r_c)` is removed exactly by the Stage-114 denominator renormalization, leaving the canonical final coefficients
\[
\kappa_c=\frac13,
\qquad
\gamma_c=\frac19.
\]

## Consequence

The outlet program is now geometrically concrete.

A fully compensated canonical branch exists whenever:

1. the core couplings lie on the exact balance surface
   \[
   g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
   \]
2. the auxiliary mixed side-channel is the first D/N half-wave of length
   \[
   L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}},
   \]
3. and its bare outgoing port is a pure-scale deformation of the canonical compact outgoing `l=2` branch.

So the remaining PDE-side question is no longer “some unknown outlet deformation.”
It is whether the actual moving-throat core realizes this specific D/N-tube + coupling-balance structure.
