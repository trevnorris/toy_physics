
# Moving-Throat PDE — Stage 131: Representative Non-Exponential Positive Mouth Families

## Goal

Evaluate the first-order Family-1 correction formulas on two explicit positive
non-exponential mouth families and on their convex interpolation.

This converts the abstract rigidity kernel of Stage 130 into concrete scale estimates.

---

## 1. Uniform broadening family

Take the broadest positive normalized profile on the mouth interval,
\[
\varsigma_u(x)=1.
\]

Then
\[
\bar g_u=\int_0^1 \cos\!\left(\frac{\pi x}{2}\right)\,dx=\frac{2}{\pi},
\]
and
\[
\bar S_u
=
\int_0^1 K_q(x)\,dx
=
\frac{2\tanh(\pi/2)}{\pi}.
\]

Numerically,
\[
\bar g_u\approx 0.636619772367581,
\qquad
\bar S_u\approx 0.583877311158896.
\]

The canonical retuning shifts are therefore
\[
\boxed{
\frac{\delta\Pi_u}{\epsilon}\approx +1.699414961314297,
\qquad
\frac{\delta\widehat T_{m,u}}{\epsilon}\approx +0.508756302215084.
}
\]

So broadening the source toward uniformity forces the canonical branch to move to
**larger** bias and **larger** traction.

---

## 2. Self-matched derivative family

Take the positive self-matched derivative profile,
\[
\varsigma_d(x)=\frac{\pi}{2}\cos\!\left(\frac{\pi x}{2}\right).
\]

Then
\[
\bar g_d=\frac{\pi}{4},
\]
and
\[
\bar S_d
=
\frac{1+\sinh(\pi/2)}{2\cosh(\pi/2)}
\approx 0.657844575502831.
\]

The canonical retuning shifts are
\[
\boxed{
\frac{\delta\Pi_d}{\epsilon}\approx -0.382993186095928,
\qquad
\frac{\delta\widehat T_{m,d}}{\epsilon}\approx -0.116943802151811.
}
\]

So sharpening the source toward the self-matched derivative profile moves the
canonical branch to **smaller** bias and **smaller** traction.

---

## 3. Convex interpolation between the two positive non-exponential families

Now interpolate the two explicit positive families:
\[
\varsigma_\lambda(x)=(1-\lambda)\varsigma_u(x)+\lambda \varsigma_d(x),
\qquad 0\le \lambda\le 1.
\]

Because the first-order formulas are affine in \((\bar g_\varsigma,\bar S_\varsigma)\),
the canonical first-order shifts are exactly
\[
\boxed{
\frac{\delta\Pi_\lambda}{\epsilon}
=
1.699414961314297
-
2.082408147410224\,\lambda,
}
\]
\[
\boxed{
\frac{\delta\widehat T_{m,\lambda}}{\epsilon}
=
0.508756302215084
-
0.625700104366895\,\lambda.
}
\]

The **bias-neutral** interpolation point is
\[
\boxed{
\lambda_{\Pi,0}\approx 0.816081594488460.
}
\]
Equivalently,
\[
1-\lambda_{\Pi,0}\approx 0.183918405511540,
\]
which agrees to numerical precision with the earlier Stage-109 broadening fraction
\[
\xi_*
=
\frac{\frac{\pi}{4}-\mathfrak g_-^{F1}}{\frac{\pi}{4}-\frac{2}{\pi}}
=
\frac{-37\sqrt3-5\pi^2+2\sqrt{4107-100\pi^2}}{5(8-\pi^2)}
\approx 0.183918405511538
\]
for the positive family that hit the Family-1 lower compensated branch.

The **traction-neutral** interpolation point is
\[
\boxed{
\lambda_{T,0}\approx 0.813099276577333.
}
\]

So the first-order finite-correction theory is perfectly consistent with the
earlier exact positive-family compensation result.

---

## 4. Meaning

For the two most natural explicit positive non-exponential mouth families:

- broadening toward uniformity changes the canonical point upward,
- sharpening toward the self-matched derivative profile changes it downward,
- and the zero-shift mixture lies very close to the earlier exact compensation
  broadening fraction.

This is the best first-order evidence so far that the canonical Family-1 mouth branch is
**rigid but not brittle**: moderate positive non-exponential corrections move it in a
controlled, almost one-parameter way.
