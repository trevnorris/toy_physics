# Moving-Throat PDE — Stage 93: Explicit Isotropic Robin Outlet Model

## Goal

Replace the abstract static deformation slot `a_0` by the first explicit isotropic geometric outlet model and compute its exact effect on the outgoing quadrupole normalization.

## Raw Robin DtN deformation

Take the canonical outgoing `l=2` branch
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6),
\qquad z:=\frac{a\omega}{c_s},
\]
and add a dimensionless isotropic Robin throat-core shift
\[
\boxed{\Lambda_2^{\rm R}(z)=\Lambda_2^{\rm out}(z)+\rho_R,}
\qquad \rho_R:=a h_R.
\]
This is the natural reduced form of an isotropic Robin-type mouth law.

Normalizing by the actual static slot,
\[
\widehat Y_2^{\rm R}(z)
:=
\frac{-3+\rho_R}{\Lambda_2^{\rm R}(z)},
\]
one finds the exact low-frequency expansion
\[
\boxed{
\widehat Y_2^{\rm R}(z)
=
1+
\frac{z^2}{9-3\rho_R}
+
\frac{4-\rho_R}{9(3-\rho_R)^2}z^4
+
 i\frac{z^5}{27-9\rho_R}
+O(z^6).
}
\]

So the raw isotropic Robin outlet changes both the even branch and the odd outgoing normalization.

## Effective quadrupole normalization factor

Relative to the canonical outgoing `l=2` fingerprint,
\[
\widehat Y_2^{\rm out}(z)=1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6),
\]
the raw Robin outlet carries the exact normalization factor
\[
\boxed{
\chi_Q^{\rm R}=rac{3}{3-\rho_R}.
}
\]
Expanding around the canonical branch,
\[
\chi_Q^{\rm R}=1+\frac{\rho_R}{3}+\frac{\rho_R^2}{9}+O(\rho_R^3).
\]
So a pure isotropic Robin core generically pushes the branch away from `chi_Q=1`.

## Branch-selection triple

In the Stage-92 linearized notation, the raw Robin outlet is the pure static deformation
\[
\boxed{(b,a_0,a_5)=(0,\rho_R,0).}
\]
Therefore the linearized outgoing-normalization shift is
\[
\delta\chi_Q^{\rm R}=\frac{\rho_R}{3}+O(\rho_R^2).
\]

## Consequence

A pure isotropic geometric Robin outlet is **not** automatically harmless. By itself it deforms both the canonical even fingerprint and the odd normalization. So if the already-fixed conservative grouped-`P_2` branch is to survive, a Robin core must either be negligible or be compensated by additional outlet structure.
