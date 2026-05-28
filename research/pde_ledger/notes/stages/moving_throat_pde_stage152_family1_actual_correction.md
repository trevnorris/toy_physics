# Moving-Throat PDE — Stage 152: Actual Family-1 Mouth Correction and One-Step Nonlinear Check

## Goal

Evaluate the exact full-profile residual \(R_*(x)\) on the explicit Family-1
canonical branch and project it against the rigidity formulas.

---

## 1. Actual first-order Family-1 correction

On the compensated Family-1 branch
\[
\Pi_* \approx 1.50882951349316,
\qquad
\Sigma_m^*\approx 0.451485277739090,
\qquad
\mathfrak g_* \approx 0.758035078944663,
\qquad
\mathcal S_* \approx 0.658075937605429,
\]
the exact weighted residual covariances are
\[
\boxed{
\operatorname{Cov}_*(c,R_*)\approx 0.0648069687666328,
}
\]
\[
\boxed{
\operatorname{Cov}_*(K_q,R_*)\approx 0.0388718368650403.
}
\]

Therefore the actual first-order moment shifts are
\[
\boxed{
\delta\mathfrak g_{\rm act}\approx -0.0648069687666328,
\qquad
\delta\mathcal S_{\rm act}\approx -0.0388718368650403.
}
\]

So the full mouth profile broadens the source relative to the tangent exponential:
the overlap factor and the mixed D/N response factor both move downward.

---

## 2. Retuned canonical point

Using the previously frozen canonical derivatives and rigidity coefficients,
\[
\mathfrak g_*' \approx 0.0714453558083195,
\qquad
A_T\approx -4.27263956256927,
\qquad
B_T\approx 0.134875005736706,
\]
the actual retuning is
\[
\boxed{
\delta\Pi_{\rm act}
\approx 0.907084414842908,
}
\]
\[
\boxed{
\delta\widehat T_{m,{\rm act}}
\approx 0.271653979462338.
}
\]

So the corrected canonical Family-1 point is
\[
\boxed{
\Pi_{\rm corr}
=
\Pi_*+\delta\Pi_{\rm act}
\approx 2.41591392833607,
}
\]
\[
\boxed{
\widehat T_{m,\rm corr}
=
\widehat T_{m,*}+\delta\widehat T_{m,\rm act}
\approx 1.17313803363654.
}
\]

This is the first actual non-exponential correction selected by the full explicit
mouth-layer model.

---

## 3. One-step nonlinear Picard check

A useful nonlinear check is to replace the exponential source by the one-step
full-profile iterate
\[
\boxed{
\Sigma_1(x)=
\frac{e^{-\Pi_*x-R_*(x)}}{\int_0^1 e^{-\Pi_*y-R_*(y)}\,dy}.
}
\]

This is not yet the full nonlinear fixed point, but it keeps the entire finite
residual \(R_*(x)\) rather than only its linearized projection.

Its actual moments are
\[
\boxed{
\mathfrak g_1\approx 0.684423574065325,
\qquad
\mathcal S_1\approx 0.616333130570251.
}
\]

Retuning the canonical branch with those exact one-step moments gives
\[
\boxed{
\Pi_1\approx 2.53914847609768,
\qquad
\widehat T_{m,1}\approx 1.21036942084359.
}
\]

So the one-step nonlinear correction shifts the canonical point slightly more than
the linearized estimate, but in the same direction and by the same overall scale.

---

## 4. Effective positive-family interpretation

Comparing the actual first-order correction to the explicit positive-family
interpolation from Stage 148 gives
\[
\lambda_{\rm eff}^{(\Pi)}\approx 0.380487632771110,
\qquad
\lambda_{\rm eff}^{(T)}\approx 0.378939241176339.
\]

So the full coupled mouth-layer correction is well approximated by a point about
\[
\boxed{
\lambda_{\rm eff}\approx 0.38
}
\]
on the positive interpolation line between the uniform family and the self-matched
derivative family.

Equivalently, the actual selected correction corresponds to a **broadening fraction**
\[
1-\lambda_{\rm eff}\approx 0.62.
\]

So the full mouth profile behaves much more like a moderate broadening toward
uniformity than like a sharpening toward the self-matched derivative branch.