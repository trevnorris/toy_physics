# Moving-Throat PDE — Stage 140 / v59 Status

## New closure result

The explicit Family-1 GNLS + localized-Maxwell program now has a **fully co-evolving
core–mouth closure**.

For a positive normalized source \(\Sigma\), the source moment
\[
\mathfrak g[\Sigma]
=
\int_0^1 \Sigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx
\]
feeds the Family-1 core loading ratio
\[
\mathcal R[\Sigma]
=
\frac{(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\qquad
\mathfrak r_{F1}\approx 1.77799353547498,
\]
which then re-enters the full mouth potential
\[
\Phi_{\Sigma_0}[\Sigma](x)
=
\Sigma_0\Bigl[
\mathcal T_s[\Sigma](x)-\mathcal R[\Sigma]\mathcal T_q[\Sigma](x)
\Bigr].
\]

So the explicit Family-1 problem is now an honest nonlinear fixed point, not a
frozen-core correction.

## What the full solve says

The numerics below are carried directly from Stage 138 for the frozen canonical
traction and Stage 139 for the renormalized canonical branch.

At the old canonical traction from Stage 138
\[
\Sigma_0^*\approx 1.80594111095636
\quad
\left(
\widehat T_{m,*}\approx 0.901484054174204
\right),
\]
the co-evolving fixed point lands at
\[
\mathfrak g_{\rm fp}\approx 0.693352419668063,
\qquad
\mathcal R_{\rm fp}\approx 0.2827139049082381,
\qquad
\Pi_{\rm fp}\approx 1.4885734438300713.
\]
So the branch survives and stays close in bias, but it is no longer exactly
compensated.

Restoring the exact lower compensated branch requires the unique renormalized
traction from Stage 139
\[
\Sigma_0^{\rm can}\approx 4.651033550168876,
\qquad
\widehat T_{m,\rm can}\approx 1.446708366456762,
\qquad
\Pi_{\rm can}\approx 3.871564377479009.
\]

## Interpretation

The mouth/core ambiguity is now essentially gone inside the explicit Family-1
closure:

- the upper compensated branch is impossible,
- the equal-normalized branch is singular,
- the lower compensated branch remains the only regular canonical branch,
- and full self-consistency promotes it to a renormalized finite-traction,
  finite-bias fixed point.

## Next serious step

The next PDE-facing task is no longer branch selection. It is to derive the actual
deviation of the moving-throat mouth/core system from this explicit Family-1
co-evolving fixed point, and then translate that deviation into the remaining
outgoing quadrupole-normalization defect.
