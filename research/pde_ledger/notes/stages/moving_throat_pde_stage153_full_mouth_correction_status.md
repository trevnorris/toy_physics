# Moving-Throat PDE — Stage 153 / v58 Status

## New closure result

The full explicit GNLS + localized-Maxwell Family-1 mouth boundary layer now selects
a **definite first-order non-exponential correction** around the unique regular lower
compensated branch in the mouth-only closure.

The exact compensated mouth potential
\[
\Phi_*(x)=4\Sigma_m^*T_s(x;\Pi_*)-\Sigma_m^*T_q(x;\Pi_*)
\]
defines a residual
\[
R_*(x)=\Phi_*(x)-\Pi_*x
\]
with
\[
R_*(0)=R_*'(0)=0,
\qquad
R_*''(0)=-3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0.
\]
So the full mouth profile is tangent-matched but sublinear at the mouth, and therefore
broadens the actual source relative to the tangent exponential branch.

## Actual first-order Family-1 correction

Projecting the exact residual against the Stage 249 rigidity kernel gives
\[
\delta\mathfrak g_{\rm act}\approx -0.0648069688,
\qquad
\delta\mathcal S_{\rm act}\approx -0.0388718369,
\]
\[
\delta\Pi_{\rm act}\approx 0.9070844148,
\qquad
\delta\widehat T_{m,\rm act}\approx 0.2716539795.
\]

So the corrected canonical point is
\[
\Pi_{\rm corr}\approx 2.4159139283,
\qquad
\widehat T_{m,\rm corr}\approx 1.1731380336.
\]

## One-step nonlinear check

Using the one-step fully exponentiated profile
\[
\Sigma_1(x)\propto e^{-\Pi_*x-R_*(x)}
\]
gives
\[
\mathfrak g_1\approx 0.6844235741,
\qquad
\mathcal S_1\approx 0.6163331306,
\]
and a retuned point
\[
\Pi_1\approx 2.5391484761,
\qquad
\widehat T_{m,1}\approx 1.2103694208.
\]

So the nonlinear correction is slightly stronger than the linear estimate but follows
the same direction and scale.
It is a one-step consistency check, not a proof of full nonlinear convergence for the
mouth-only Picard iteration.

## Interpretation

Inside the explicit Family-1 closure, and before full core-mouth co-evolution is
turned back on:

- branch selection is finished;
- the lower compensated branch remains the unique regular branch in this mouth-only correction analysis;
- the full coupled mouth profile does **not** destroy it;
- but it shifts the preferred point upward in both bias and normalized traction.

Numerical stress tests confirm that the linear correction is quantitatively good for
small residual loads and remains directionally correct at the full one-step nonlinear
load, but it should not be read as an exact finite-amplitude law.

The mouth-side problem is therefore no longer one of branch ambiguity. It is now a
finite correction problem around a unique regular branch.

## Next serious step

The next PDE-facing step is to let the **core outlet coefficients and the mouth profile
co-evolve** self-consistently, instead of holding the compensated core branch fixed while
correcting only the mouth source law.
