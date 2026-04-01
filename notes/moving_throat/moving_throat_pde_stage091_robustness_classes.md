# Moving-Throat PDE — Stage 91: Robustness Classes for `chi_Q`

## Goal

Classify which explicit isotropic DtN perturbations leave the canonical outgoing normalization untouched and which ones genuinely shift it.

## Class A — pure scale deformation

Take
\[
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(z).
\]
Then
\[
\widehat Y_2^{\rm def}(z)=\widehat Y_2^{\rm out}(z),
\]
exactly, so
\[
\boxed{\chi_Q=1.}
\]
Thus pure mouth normalization is invisible to the normalized outgoing quadrupole fingerprint.

## Class B — pure scale+argument deformation

Take
\[
\Lambda_2^{\rm def}(z)=S\Lambda_2^{\rm out}(\beta z).
\]
Then
\[
\widehat Y_2^{\rm def}(z)
=
1+\frac{\beta^2 z^2}{9}+\frac{4\beta^4 z^4}{81}+i\frac{\beta^5 z^5}{27}+O(z^6).
\]
So if the conservative even fingerprint is kept canonical,
\[
\beta^2=1,
\qquad
\beta^4=1,
\]
which on the natural positive branch forces
\[
\boxed{\beta=1.}
\]
Hence
\[
\boxed{\chi_Q=1.}
\]
So a pure effective radius/sound-speed rescaling cannot move the canonical outgoing normalization without simultaneously spoiling the already fixed even moments.

## Class C — additive isotropic throat-core channel

Set `beta = 1` but allow a genuine throat-core self-energy,
\[
\Lambda_2^{\rm def}(z)
=
S\Lambda_2^{\rm out}(z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5.
\]
If the even moments are held canonical, then
\[
\Sigma_2=-\frac{\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{\Sigma_0}{27},
\]
and the odd normalization is
\[
\boxed{
\chi_Q=\frac{3(S+9\Sigma_5)}{3S-\Sigma_0}.
}
\]
So an additive throat-core channel can move `chi_Q` even while leaving the lower even moments unchanged.

A particularly clean special case is a purely even additive core with `Sigma_5 = 0`:
\[
\boxed{
\chi_Q=\frac{3S}{3S-\Sigma_0}.
}
\]
So a static additive core shift is a genuine candidate branch-selection effect.

## Exact preservation submanifold

The deformation preserves the canonical outgoing normalization iff
\[
\chi_Q=1.
\]
The exact condition is
\[
\boxed{
\Sigma_5=\frac{S(1-\beta^5)}{9}-\frac{\Sigma_0}{27}.
}
\]
This includes as special cases:

- pure scale deformation,
- pure scale+argument deformation with `beta = 1`,
- additive core deformations whose odd slot is locked to the static shift.

## Consequence

The canonical value `chi_Q = 1` is robust against:

1. overall mouth normalization,
2. pure effective radius/sound-speed rescaling once the conservative even fingerprint is fixed.

It is shifted only by a genuine isotropic throat-core self-energy that is not on the exact preservation submanifold above.
