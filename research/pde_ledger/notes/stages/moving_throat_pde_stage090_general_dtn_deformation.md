# Moving-Throat PDE — Stage 90: General Isotropic `l=2` DtN Deformation Algebra

## Goal

Replace the symbolic branch parameter `xi_Q` by an explicit low-frequency moving-throat DtN deformation model and derive the exact map from deformation data to the retarded quadrupole normalization factor `chi_Q`.

## Canonical outgoing branch

From the exact outgoing spherical `l=2` DtN model,
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6),
\qquad z:=\frac{a\omega}{c_s}.
\]

The corresponding normalized branch is
\[
\widehat Y_2^{\rm out}(z)
=1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6),
\]
so the canonical outgoing normalization is `chi_Q = 1`.

## General isotropic deformation model

Take the first explicit isotropic moving-throat deformation family
\[
\boxed{
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6).
}
\]
Here:

- `S` is an overall mouth normalization,
- `beta` rescales the effective outgoing argument,
- `Sigma_0, Sigma_2, Sigma_4` are isotropic throat-core even self-energy data,
- `Sigma_5` is the first extra isotropic odd `l=2` core outlet.

Expanding,
\[
\Lambda_2^{\rm def}(z)=L_0+L_2 z^2+L_4 z^4+iL_5 z^5+O(z^6),
\]
with
\[
L_0=-3S+\Sigma_0,
\qquad
L_2=\frac{S\beta^2}{3}+\Sigma_2,
\qquad
L_4=\frac{S\beta^4}{9}+\Sigma_4,
\qquad
L_5=\frac{S\beta^5}{9}+\Sigma_5.
\]

## Normalized deformation law

Normalize by the actual static slot,
\[
\widehat Y_2^{\rm def}(z):=\frac{L_0}{L_0+L_2 z^2+L_4 z^4+iL_5 z^5}+O(z^6).
\]
Then the exact low-frequency coefficients are
\[
\boxed{
\widehat Y_2^{\rm def}(z)
=
1-\frac{L_2}{L_0}z^2
+\left(\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}\right)z^4
-i\frac{L_5}{L_0}z^5
+O(z^6).
}
\]

So the deformation-normalized quadrupole factor is
\[
\chi_Q=\frac{-L_5/L_0}{1/27}.
\]

## Exact canonical-even matching conditions

Demand that the deformed branch preserve the canonical conservative even fingerprint,
\[
\frac{z^2}{9},\qquad \frac{4z^4}{81}.
\]
Then
\[
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81}.
\]
Solving for the even throat-core coefficients gives
\[
\boxed{
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27}.
}
\]

With those imposed, the exact odd normalization becomes
\[
\boxed{
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
}
\]
Equivalently,
\[
\boxed{
\chi_Q-1=
\frac{3S(\beta^5-1)+\Sigma_0+27\Sigma_5}{3S-\Sigma_0}.
}
\]

## Consequence

This is the first explicit moving-throat DtN deformation model for the last reduced 2.5PN scalar. It shows exactly which isotropic branch data can move the canonical value `chi_Q = 1`:

- argument deformation `beta`,
- static additive throat-core shift `Sigma_0`,
- odd `l=2` throat-core outlet `Sigma_5`.

Overall scale `S` is not itself an independent obstruction; it only enters through the ratios above.
