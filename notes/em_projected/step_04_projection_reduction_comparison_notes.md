# Projection-First Versus Reduction-First Coupling Notes

## Purpose

This note records the zero-mode coupling comparison audited by
`step_04_projection_reduction_comparison_sympy.py`.

The goal is to show exactly when projection-first Maxwell reproduces the
reduction-first effective coupling, and when it remains observer dependent.

## Zero-Mode Projection

Assume

\[
F^{\mu\nu}(x,w)=f^{\mu\nu}(x),
\qquad
F^{w\nu}=0,
\qquad
J^\nu(x,w)=j^\nu(x)S(w).
\]

Then the projected inhomogeneous equation reduces to

\[
I_{WZ}\,\partial_\mu f^{\mu\nu}
=
\mu_0 I_{WS}j^\nu,
\]

with

\[
I_{WZ}=\int WZ\,dw,
\qquad
I_{WS}=\int WS\,dw.
\]

After dividing by `I_WZ`, the projection-first effective coupling is

\[
\mu_{0,\rm eff}^{\rm proj}
=
\mu_0\frac{I_{WS}}{I_{WZ}}.
\]

Reduction-first gives

\[
\mu_{0,\rm eff}^{\rm red}
=
\frac{\mu_0}{Z_{\rm int}},
\qquad
Z_{\rm int}=\int Z\,dw.
\]

So the two agree only if

\[
\frac{I_{WS}}{I_{WZ}}=\frac{1}{Z_{\rm int}}.
\]

The audit script no longer treats this as a bare factorization claim. It also
checks a concrete smooth projected residual with

\[
f(x)=\sin(kx)+x^2,\qquad j(x)=\cos(kx)+x,
\]

Gaussian observer/source profiles, and the same Gaussian localizer \(Z(w)\).
It verifies directly that

\[
\int W\left(Z\,\partial_x f-\mu_0 S\,j\right)dw
=
I_{WZ}\partial_x f-\mu_0 I_{WS}j.
\]

Then it mutates the assumptions by adding \(w^2\)-dependent field and source
pieces. Those mutations produce nonzero extra projected moments, so the
zero-mode/factorized-source reduction fails as it should.

## Gaussian Matched-Kernel Example

Take

\[
Z(w)=e^{-w^2/\lambda^2},
\qquad
W(w)=\frac{Z(w)}{Z_{\rm int}}.
\]

Then

\[
Z_{\rm int}=\sqrt{\pi}\lambda,
\qquad
\int Z^2dw=\frac{\sqrt{2\pi}\lambda}{2},
\qquad
I_{WZ}=\frac{1}{\sqrt2}.
\]

For a delta-localized source, `S(w)=delta(w)`,

\[
I_{WS}=W(0)=\frac{1}{\sqrt{\pi}\lambda},
\]

so

\[
\mu_{0,\rm eff}^{\rm proj}
=
\frac{\sqrt2\,\mu_0}{\sqrt{\pi}\lambda}
=
\sqrt2\,\mu_{0,\rm eff}^{\rm red}.
\]

Thus projection-first does not automatically reproduce the reduced Maxwell
coupling.

## Regularized Sharp-Layer Check

The audit script also replaces the distribution-undefined sharp-slice
`delta(w) * delta(w)` case with a normalized Gaussian regulator:

\[
W_\epsilon(w)=S_\epsilon(w)
=
\frac{e^{-w^2/\epsilon^2}}{\sqrt{\pi}\epsilon}.
\]

For the same Gaussian localizer,

\[
I_{WZ}(\epsilon)
=
\frac{\lambda}{\sqrt{\epsilon^2+\lambda^2}}
\longrightarrow 1
\qquad (\epsilon\to0),
\]

while

\[
I_{WS}(\epsilon)
=
\frac{\sqrt2}{2\sqrt{\pi}\epsilon}
\]

diverges as the layer is squeezed. So an exact delta observer measuring an exact
delta source is not a finite coupling prescription; it must be regulated or
replaced by a matched distributed source profile.

## Script

- `step_04_projection_reduction_comparison_sympy.py`
