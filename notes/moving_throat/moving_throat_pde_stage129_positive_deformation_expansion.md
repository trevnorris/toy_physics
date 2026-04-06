
# Moving-Throat PDE — Stage 129: First-Order Expansion for Positive Mouth-Layer Deformations

## Goal

Test how rigid the unique regular Family-1 canonical mouth branch is under **non-exponential but still positive localized** mouth-source deformations.

The canonical exponential branch already fixed
\[
\Pi_* \approx 1.50882951349316,
\qquad
\widehat T_{m,*}\approx 0.901484054174205,
\qquad
\mathfrak g_*=\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathcal S_*=\mathcal S_q(\Pi_*)\approx 0.658075937605429.
\]

The question now is: if the actual mouth layer is **not exactly exponential**, how much must the canonical point move?

---

## 1. Positive normalized deformation family

Work on the dimensionless mouth interval
\[
x=\frac{z}{L}\in[0,1].
\]
Let the canonical exponential source be
\[
\Sigma_*(x)=\frac{\Pi_* e^{-\Pi_* x}}{1-e^{-\Pi_*}}.
\]

Take any other positive normalized mouth profile
\[
\varsigma(x)\ge 0,
\qquad
\int_0^1 \varsigma(x)\,dx=1,
\]
and define the exact convex positive deformation family
\[
\boxed{
\Sigma_\epsilon(x)=(1-\epsilon)\Sigma_*(x)+\epsilon\,\varsigma(x),
\qquad
0\le \epsilon\ll 1.
}
\]

This preserves positivity and normalization exactly.

---

## 2. The only two source moments that matter at first order

Define the Family-1 overlap and mixed-channel kernels
\[
c(x)=\cos\!\left(\frac{\pi x}{2}\right),
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

Then the source shape enters the mouth-core system only through the two averages
\[
\bar g[\sigma]=\int_0^1 \sigma(x)\,c(x)\,dx,
\qquad
\bar S[\sigma]=\int_0^1 \sigma(x)\,K_q(x)\,dx.
\]

For the positive convex family this is exact:
\[
\boxed{
\bar g_\epsilon
=
\mathfrak g_*+\epsilon(\bar g_\varsigma-\mathfrak g_*),
\qquad
\bar S_\epsilon
=
\mathcal S_*+\epsilon(\bar S_\varsigma-\mathcal S_*),
}
\]
where
\[
\bar g_\varsigma=\int_0^1\varsigma(x)c(x)\,dx,
\qquad
\bar S_\varsigma=\int_0^1\varsigma(x)K_q(x)\,dx.
\]

So the first-order mouth-shape problem is **two-dimensional**:
only \((\bar g_\varsigma,\bar S_\varsigma)\) matter.

---

## 3. Retuning the electrochemical bias to stay on the canonical lower branch

To remain on the canonical outgoing-preserving lower branch, the overlap must stay fixed at
\[
\bar g = \mathfrak g_*.
\]
So if the source is deformed, the exponential bias must shift
\[
\Pi=\Pi_*+\delta\Pi.
\]

Let
\[
\mathfrak g'_*=
\left.\frac{d\mathfrak g_\Pi}{d\Pi}\right|_{\Pi_*},
\qquad
\mathcal S'_*
=
\left.\frac{d\mathcal S_q}{d\Pi}\right|_{\Pi_*}.
\]
Then the first-order compensation law is
\[
\boxed{
\delta\Pi
=
-\epsilon\,\frac{\bar g_\varsigma-\mathfrak g_*}{\mathfrak g'_*}.
}
\]

Numerically,
\[
\boxed{
\mathfrak g'_* \approx 0.0714453558083195,
\qquad
\mathcal S'_* \approx 0.0483709542125041.
}
\]

So the compensated mixed-channel response changes by
\[
\boxed{
\delta \mathcal S_q
=
\epsilon\left[
(\bar S_\varsigma-\mathcal S_*)
-
\frac{\mathcal S'_*}{\mathfrak g'_*}\,(\bar g_\varsigma-\mathfrak g_*)
\right].
}
\]

This is the first non-exponential correction formula at `O(epsilon)`.
