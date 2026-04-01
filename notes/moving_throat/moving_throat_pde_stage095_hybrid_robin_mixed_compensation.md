# Moving-Throat PDE — Stage 95: Exact Robin–Mixed Compensation Law

## Goal

Combine the geometric Robin core and the hidden mixed side-channel and determine whether an explicit compensated moving-throat outlet can preserve the canonical outgoing quadrupole branch.

## Hybrid outlet model

Take
\[
\boxed{
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+ho_R
-
\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}
+O(z^6).
}
\]
Expanding,
\[
L_0=-3+\rho_R-\sigma_W,
\qquad
L_2=\frac13-\sigma_W\kappa_W,
\qquad
L_4=\frac19-\sigma_W\kappa_W^2,
\qquad
L_5=\frac19-\sigma_W\gamma_W.
\]

## Exact canonical-even solutions

Imposing the canonical conservative even fingerprint,
\[
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81},
\]
yields exactly two branches:
\[
\boxed{\rho_R=\sigma_W,\qquad \kappa_W=0,}
\]
or
\[
\boxed{\rho_R=4\sigma_W,\qquad \kappa_W=\frac13.}
\]

The first is a trivial cancellation branch: the static Robin shift simply cancels the static mixed loading.
The second is the nontrivial compensated branch.

## Nontrivial compensated branch

On
\[
\rho_R=4\sigma_W,
\qquad
\kappa_W=\frac13,
\]
one finds
\[
\Lambda_2^{\rm hyb}(z)
=
-3(1-\sigma_W)
+
\frac{1-\sigma_W}{3}z^2
+
\frac{1-\sigma_W}{9}z^4
+
 i\left(\frac19-\sigma_W\gamma_W\right)z^5
+O(z^6).
\]
So the exact outgoing-normalization factor is
\[
\boxed{
\chi_Q^{\rm hyb}
=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
}
\]
Canonical outgoing normalization is preserved iff
\[
\boxed{\gamma_W=\frac19.}
\]
With that value,
\[
\boxed{
\Lambda_2^{\rm hyb}(z)=(1-\sigma_W)\Lambda_2^{\rm out}(z)+O(z^6),
}
\]
so the whole hybrid outlet collapses to the harmless pure-scale deformation class.

## Branch-selection data

On the nontrivial compensated branch, the net Stage-92 deformation data are
\[
\boxed{(b,a_0,a_5)=\left(0,\,3\sigma_W,\,-\sigma_W\gamma_W\right).}
\]
Hence the linearized preservation condition becomes
\[
\frac{a_0}{3}+9a_5=\sigma_W(1-9\gamma_W)=0,
\]
again giving
\[
\boxed{\gamma_W=\frac19.}
\]

## Consequence

This is the first explicit compensated moving-throat outlet model that preserves the canonical outgoing quadrupole branch. It shows that neither a pure geometric Robin core nor a naive hidden mixed pole is sufficient by itself, but a specific Robin–mixed balance law can reduce the whole deformation to a pure mouth renormalization, which is exactly the robust class already identified in Stages 90–92.
