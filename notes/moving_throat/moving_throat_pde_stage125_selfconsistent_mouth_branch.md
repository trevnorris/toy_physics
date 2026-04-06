# Moving-Throat PDE — Stage 125: Self-Consistent Mouth-Branch Law

## Goal

Close the loop between:

1. the explicit positive mouth-layer source family
   \[
   \mathfrak g_\Pi
   =
   \frac{2\Pi(2\Pi e^\Pi+\pi)}
   {(4\Pi^2+\pi^2)(e^\Pi-1)},
   \]
2. the explicit core-to-mouth gain map
   \[
   M_q=-M_s\,\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2},
   \]
3. and the coupled Family-1 mouth fixed-point law
   \[
   \Pi=M_s+M_q\,\mathcal S_q(\Pi).
   \]

This removes the last free branch label \(\mathfrak g_c\) from the mouth problem.

---

## 1. Self-consistent Family-1 gain law

On the explicit positive mouth-layer branch we must identify
\[
\mathfrak g_c=\mathfrak g_\Pi.
\]
So the Family-1 mixed-to-shell gain ratio becomes
\[
\boxed{
R_q(\Pi):=-\frac{M_q}{M_s}
=
\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
}
\]

The coupled Family-1 mouth equation then closes to
\[
\boxed{
\Pi
=
\Sigma_0\Big[1-R_q(\Pi)\,\mathcal S_q(\Pi)\Big],
}
\]
where
\[
\Sigma_0:=M_s,
\qquad
\mathcal S_q(\Pi)=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
\]

Equivalently, the shell gain required to realize a given bias \(\Pi\) is
\[
\boxed{
\Sigma_0(\Pi)
=
\frac{\Pi}{1-R_q(\Pi)\mathcal S_q(\Pi)}.
}
\]

So the explicit mouth branch is now one scalar function \(\Sigma_0(\Pi)\), not a free gain pair.

---

## 2. Self-matched mouth susceptibility closure

Under the self-matched closure from Stage 123,
\[
\Sigma_0=\frac{20}{9}\,\widehat T_m^2.
\]
Therefore the normalized mouth-traction branch is
\[
\boxed{
\widehat T_m(\Pi)
=
\sqrt{\frac{9\Pi}{20\left[1-R_q(\Pi)\mathcal S_q(\Pi)\right]}}.
}
\]

This formula is the first explicit self-consistent map from the parent mouth bias \(\Pi\)
to the physical normalized traction \(\widehat T_m\).

---

## 3. Canonical Family-1 point

The already-frozen lower compensated Family-1 mouth point is
\[
\mathfrak g_-^{F1}
=
\mathfrak r_{F1}-\frac12\sqrt{1+\mathfrak r_{F1}^2}
\approx 0.758035078944663,
\]
and the unique mouth-layer bias solving
\[
\mathfrak g_\Pi=\mathfrak g_-^{F1}
\]
is
\[
\boxed{
\Pi_*\approx 1.50882951349316.
}
\]

At that point
\[
R_q(\Pi_*)=\frac14,
\qquad
\mathcal S_q(\Pi_*)\approx 0.658075937605429,
\]
so
\[
\boxed{
\Sigma_0(\Pi_*)\approx 1.80594111095636,
\qquad
\widehat T_m(\Pi_*)\approx 0.901484054174205.
}
\]

This reproduces the earlier outlet-consistent canonical branch exactly, but now from the
fully self-consistent positive mouth-layer + core map inside this explicit closure.

---

## 4. Meaning

The mouth problem is now much sharper:

- the source shape is no longer free,
- the gain ratio is no longer free,
- and the shell traction is no longer free.

Once this explicit positive mouth layer is assumed, the branch is governed by one scalar variable
\(\Pi\), with a unique canonical compensation point at moderate finite bias and moderate
finite traction inside that closure.
