# Moving-Throat PDE — Stage 157 / v59 Status

## New closure result

The explicit Family-1 GNLS + localized-Maxwell program now has a **fully co-evolving
core–mouth closure** on the analyzed positive branch window.

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

So the explicit Family-1 problem is now an honest nonlinear fixed-point problem
inside this reduced closure, not merely a frozen-core correction.

## What the full solve says

The numerics below are carried directly from Stage 155 for the frozen canonical
traction on the analyzed fixed-point branch and from Stage 156 for the
renormalized canonical branch.

At the old canonical traction from Stage 155, the solved positive fixed point on the
analyzed branch lands at
\[
\Sigma_0^*\approx 1.80594111095636
\quad
\left(
\widehat T_{m,*}\approx 0.901484054174204
\right),
\qquad
\mathfrak g_{\rm fp}\approx 0.693352419668063,
\qquad
\mathcal R_{\rm fp}\approx 0.2827139049082381,
\qquad
\Pi_{\rm fp}\approx 1.4885734438300713.
\]
So the branch survives and stays close in bias, but it is no longer exactly
compensated.

Restoring the exact lower compensated branch requires the numerically located
renormalized traction from Stage 156 on the analyzed monotone bracket
\[
\Sigma_0^{\rm can}\approx 4.651033550168867,
\qquad
\widehat T_{m,\rm can}\approx 1.4467083664567613,
\qquad
\Pi_{\rm can}\approx 3.871564377479002.
\]

## Interpretation

The mouth/core ambiguity is now essentially gone inside the explicit Family-1
closure on the analyzed positive branch window:

- the upper compensated branch is impossible,
- the equal-normalized branch is singular,
- the lower compensated branch remains the only regular canonical branch inside the explicit Family-1 closure,
- and full self-consistency promotes it to a renormalized finite-traction,
  finite-bias fixed point.

This stage identifies the reduced co-evolving target. It does **not** yet prove that
the full moving-throat PDE branch dynamically realizes that target without additional
deviation analysis.

## Next serious step

The next PDE-facing task is no longer branch selection. It is to derive the actual
deviation of the moving-throat mouth/core system from this explicit Family-1
co-evolving fixed point, and then translate that deviation into the remaining
outgoing quadrupole-normalization defect.
