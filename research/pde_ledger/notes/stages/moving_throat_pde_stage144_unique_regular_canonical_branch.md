# Moving-Throat PDE — Stage 144: Unique Regular Canonical Mouth Branch

## Goal

Combine the positive-source theorem, the explicit mouth-layer bias family, and the
self-consistent core-to-mouth gain law to decide which compensated branch survives as a
**regular finite-bias / finite-traction** Family-1 solution.

---

## 1. Upper branch remains impossible

The explicit Family-1 compensated branches are
\[
\mathfrak g_\pm^{F1}
=
\mathfrak r_{F1}\pm\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
with
\[
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
\]

But every positive normalized source profile on the first D/N interval satisfies
\[
0\le \mathfrak g_c\le 1.
\]
So
\[
\boxed{
\mathfrak g_+^{F1}>1
\quad\Rightarrow\quad
\text{the upper compensated branch is impossible.}
}
\]

---

## 2. Lower branch is uniquely reachable at finite bias

From Stage 232, \(\mathfrak g_\Pi\) is strictly monotone increasing in \(\Pi\), with range
\[
\frac{2}{\pi}<\mathfrak g_\Pi<1
\qquad (\Pi>0).
\]
Since
\[
\frac{2}{\pi}\approx 0.636619772367581
<
\mathfrak g_-^{F1}\approx 0.758035078944663
<
1,
\]
there exists a unique positive finite \(\Pi_*\) such that
\[
\mathfrak g_{\Pi_*}=\mathfrak g_-^{F1}.
\]
Numerically,
\[
\boxed{
\Pi_*\approx 1.50882951349316.
}
\]

So the lower compensated branch is not merely allowed. Inside the explicit exponential
positive mouth family it is the **unique finite-bias** compensated branch.

---

## 3. Finite regular traction at the canonical point

At \(\Pi=\Pi_*\), Stage 244 gives
\[
\Sigma_0(\Pi_*)\approx 1.80594111095636,
\qquad
\widehat T_m(\Pi_*)\approx 0.901484054174205.
\]
These are finite and moderate.

For comparison, the self-matched derivative-profile point \(\mathfrak g=\pi/4\) occurs at
\[
\Pi_{\rm match}\approx 1.90848600654854,
\qquad
\widehat T_m(\Pi_{\rm match})\approx 1.01132972803599.
\]
So the canonical compensated branch is reached **before** the self-matched derivative point
and far before the singular point-source limit.

---

## 4. Final branch-selection theorem

Within the explicit Family-1 core + positive exponential mouth-layer closure:

- the upper compensated branch is impossible,
- the equal-normalized branch is a singular limit,
- and the lower compensated branch is the unique regular finite-bias / finite-traction
  branch that preserves the canonical outgoing quadrupole fingerprint.

Equivalently,
\[
\boxed{
\text{inside the explicit Family-1 positive exponential mouth-layer closure, the canonical mouth branch is uniquely selected as the regular finite positive-source branch.}
}
\]

---

## Meaning

This is a real narrowing.

The remaining mouth ambiguity is no longer a branch-choice ambiguity at all.
On the explicit Family-1 positive mouth-layer closure, the canonical outgoing-preserving
branch is the only regular finite one.
