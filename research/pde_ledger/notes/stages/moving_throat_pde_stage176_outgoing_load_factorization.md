# Moving-Throat PDE — Stage 176: Outgoing Load-Factor Factorization and the Square-Root Mixed-Leg Law

## Purpose

Stage 243 reduced the remaining linear grouped weak-axisymmetric defect to the outgoing-load theorem
\[
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r=\delta_K,
\qquad
\Lambda_r=\frac{P_r}{\Delta_r},
\]
on conservative-shape-preserving branches. The next honest step is therefore not to reopen the whole static self-similarity story, but to **factor the outgoing load itself** into the smallest set of microscopic slippage variables.

That is what this stage does.

The main result is that the full outgoing-load defect splits exactly into three pieces:

1. a wall-normalized **raw mixed-leg load**,
2. a dimensionless **interference ratio**,
3. a dimensionless **hybridization ratio**.

On branches where the last two are rigid, the entire remaining linear grouped `2.5`PN defect collapses to a single square-root wall-loading law for the raw mixed leg.

---

## 1. Starting point

From the earlier grouped Maxwell/mixed reduction, the outgoing load factor in port `r` is
\[
\Lambda_r=\frac{P_r}{\Delta_r},
\]
with
\[
P_r=\Omega_{U,r}^2 G_{W,r}+R_r G_{U,r},
\qquad
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2.
\]
For consistency with Stage 243, the outgoing-leg couplings are written here as
\(G_{W,r}\) and \(G_{U,r}\); these are the same quantities previously denoted
\(g_{W,r}\) and \(g_{U,r}\).
The Stage 243 defect field was
\[
\Sigma_r^{(N)}:=\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right),
\]
so that on conservative-shape-preserving branches
\[
\Xi_{\rm load}=\sum_r \rho_r^{(N)}\,\Sigma_r^{(N)}.
\]
The remaining theorem gate was exactly to compute the grouped weak-axisymmetric drift of \(\Lambda_r\) on the actual moving-throat branch.

---

## 2. Exact factorization of the outgoing load

Define the three microscopic dimensionless objects
\[
\mathcal M_r:=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
\qquad
\mathcal I_r:=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
\qquad
\mathcal H_r:=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}.
\]
Then the outgoing load factor obeys the exact identity
\[
\boxed{
\frac{\Lambda_r^2}{K}
=
\mathcal M_r^2\,
\frac{(1+\mathcal I_r)^2}{(1-\mathcal H_r)^2}.
}
\]
This is just algebra, but it is the right algebra.

It shows that the outgoing load is not one opaque object. It is the product of:

- the wall-normalized raw mixed leg \(\mathcal M_r\),
- the interference correction \((1+\mathcal I_r)\),
- the hybridization denominator \((1-\mathcal H_r)^{-1}\).

---

## 3. Exact decomposition of the outgoing defect field

Taking the logarithmic variation gives
\[
\boxed{
\Sigma_r^{(N)}
=
2\,\delta\ln\mathcal M_r
+2\,\delta\ln(1+\mathcal I_r)
-2\,\delta\ln(1-\mathcal H_r).
}
\]
So the Stage 243 grouped defect becomes
\[
\boxed{
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}
\Big[
2\,\delta\ln\mathcal M_r
+2\,\delta\ln(1+\mathcal I_r)
-2\,\delta\ln(1-\mathcal H_r)
\Big]
}
\]
on conservative-shape-preserving branches.

This is the cleanest exact form reached so far. The last linear grouped `2.5`PN bottleneck is now the weighted failure of three microscopic loading factors to track one another.

---

## 4. First-order transport of the three microscopic slippages

The three logarithmic variations are
\[
\delta\ln\mathcal M_r
=
\delta\ln G_{W,r}-\delta\ln\Omega_{W,r}^2-\frac12\,\delta_K,
\]
\[
\delta\ln\mathcal I_r
=
\delta\ln R_r+\delta\ln G_{U,r}-\delta\ln\Omega_{U,r}^2-\delta\ln G_{W,r},
\]
\[
\delta\ln\mathcal H_r
=
2\,\delta\ln R_r-\delta\ln\Omega_{U,r}^2-\delta\ln\Omega_{W,r}^2.
\]
Therefore the first-order defect field can be written exactly as
\[
\boxed{
\Sigma_r^{(N)}
=
2\,\delta\ln\mathcal M_r
+
\frac{2\mathcal I_r}{1+\mathcal I_r}\,\delta\ln\mathcal I_r
+
\frac{2\mathcal H_r}{1-\mathcal H_r}\,\delta\ln\mathcal H_r.
}
\]
This is the sharpest microscopic decomposition of the outgoing-load theorem.

The three channels have a clear interpretation:

- \(\delta\ln\mathcal M_r\): raw mixed-leg wall-loading mismatch,
- \(\delta\ln\mathcal I_r\): interference-ratio slippage,
- \(\delta\ln\mathcal H_r\): hybridization-ratio slippage.

---

## 5. Fully expanded first-order transport in the primitive port variables

Expanding the previous formula gives
\[
\boxed{
\Sigma_r^{(N)}
=
-\delta_K
+\frac{2}{1+\mathcal I_r}\,\delta\ln G_{W,r}
+\frac{2\mathcal I_r}{1+\mathcal I_r}\,\delta\ln G_{U,r}
}
\]
\[
\boxed{
\qquad
+2\left(\frac{\mathcal I_r}{1+\mathcal I_r}+\frac{2\mathcal H_r}{1-\mathcal H_r}\right)\delta\ln R_r
-2\left(\frac{\mathcal I_r}{1+\mathcal I_r}+\frac{\mathcal H_r}{1-\mathcal H_r}\right)\delta\ln\Omega_{U,r}^2
-\frac{2}{1-\mathcal H_r}\,\delta\ln\Omega_{W,r}^2.
}
\]
So the remaining linear grouped defect has now been reduced to a concrete weighted drift law in the primitive outgoing-port variables.

---

## 6. Conservative-shape theorem in the new variables

Stage 243 already showed that if the conservative shapes are preserved, then
\[
\Xi_{\rm load}=\sum_r \rho_r^{(N)}\Sigma_r^{(N)}.
\]
The new factorization says more.

If, in addition, the interference and hybridization ratios are rigid,
\[
\delta\ln\mathcal I_r=0,
\qquad
\delta\ln\mathcal H_r=0,
\]
then
\[
\boxed{
\Sigma_r^{(N)}=2\,\delta\ln\mathcal M_r
=2\,\delta\ln\!\left(\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K}\right).
}
\]
So on conservative-shape-preserving, interference-rigid, hybridization-rigid branches,
\[
\boxed{
\Xi_{\rm load}
=2\sum_r \rho_r^{(N)}\,\delta\ln\!\left(\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K}\right).
}
\]
This is a much sharper theorem than Stage 243.

---

## 7. The square-root mixed-leg law

The previous formula immediately gives a strong sufficient condition for the full linear grouped defect to vanish.

If every outgoing port obeys the lane-by-lane wall-loading law
\[
\boxed{
\delta\ln\mathcal I_r=0,
\qquad
\delta\ln\mathcal H_r=0,
\qquad
\delta\ln\!\left(\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K}\right)=0,
}
\]
then
\[
\Sigma_r^{(N)}=0
\qquad \text{for every port }r,
\]
and therefore
\[
\boxed{\Xi_{\rm load}=0.}
\]
Equivalently, the raw mixed leg must scale as
\[
\boxed{
\frac{G_{W,r}}{\Omega_{W,r}^2}\propto \sqrt K.
}
\]
This is the stage’s main positive theorem.

It is the exact replacement for the old vague statement that the outgoing bundle must “track the wall.”
The sharp statement is:

> once the interference and hybridization geometry are rigid, the remaining linear grouped `2.5`PN defect vanishes exactly when the raw mixed leg obeys a square-root wall-loading law.

---

## 8. Dominant-port corollary

If a single outgoing port `r_*` dominates the weight distribution,
\[
\rho_{r_*}^{(N)}\approx1,
\qquad
\rho_{r\neq r_*}^{(N)}\approx0,
\]
then
\[
\boxed{
\Xi_{\rm load}\approx \Sigma_{r_*}^{(N)}.
}
\]
Under interference/hybridization rigidity this collapses further to
\[
\boxed{
\Xi_{\rm load}
\approx
2\,\delta\ln\!\left(\frac{G_{W,*}}{\Omega_{W,*}^2\sqrt K}\right).
}
\]
So in the dominant-port limit the whole remaining linear grouped defect is just a measurement of how the leading mixed leg fails to satisfy the square-root wall-loading law.

---

## 9. What Stage 244 changes

Stage 243 reduced the remaining linear grouped defect to the outgoing-load theorem
\[
2\sum_r \rho_r^{(N)}\delta\ln\Lambda_r=\delta_K.
\]
Stage 244 sharpens that result substantially.

It shows that the remaining outgoing-load theorem is not one indivisible condition. It splits into:

1. raw mixed-leg wall loading,
2. interference-ratio rigidity,
3. hybridization-ratio rigidity.

So the next honest theorem gate is no longer simply “compute \(\delta\ln\Lambda_r\).”
It is:

> compute the grouped weak-axisymmetric drifts of
> \[
> \mathcal M_r=\frac{G_{W,r}}{\Omega_{W,r}^2\sqrt K},
> \qquad
> \mathcal I_r=\frac{R_r G_{U,r}}{\Omega_{U,r}^2 G_{W,r}},
> \qquad
> \mathcal H_r=\frac{R_r^2}{\Omega_{U,r}^2\Omega_{W,r}^2}
> \]
> on the actual moving-throat branch.

Once those three microscopic slippages are known, the remaining linear grouped `2.5`PN defect is already explicit.
