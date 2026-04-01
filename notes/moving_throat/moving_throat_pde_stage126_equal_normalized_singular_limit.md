# Moving-Throat PDE — Stage 126: Equal-Normalized Branch Is a Singular Limit

## Goal

Determine whether the simple equal-normalized mouth-source branch
\[
\mathfrak g_c=1
\]
can occur at any **finite** positive mouth-layer bias \(\Pi\), and whether it corresponds
to a regular finite-traction mouth state.

---

## 1. Exact strict inequality for the exponential mouth layer

For the explicit positive mouth-layer family,
\[
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}
{(4\Pi^2+\pi^2)(e^\Pi-1)}.
\]
A direct rearrangement gives
\[
1-\mathfrak g_\Pi
=
\frac{\pi^2(e^\Pi-1)-2\pi\Pi-4\Pi^2}
{(4\Pi^2+\pi^2)(e^\Pi-1)}.
\]

The numerator splits exactly as
\[
\pi^2\!\left(e^\Pi-1-\Pi-\frac{\Pi^2}{2}\right)
+\Pi(\pi^2-2\pi)
+\Pi^2\!\left(\frac{\pi^2}{2}-4\right).
\]
Every term is positive for \(\Pi>0\), since
\[
e^\Pi-1-\Pi-\frac{\Pi^2}{2}>0,
\qquad
\pi^2-2\pi>0,
\qquad
\frac{\pi^2}{2}-4>0.
\]

Therefore
\[
\boxed{
0<\mathfrak g_\Pi<1
\qquad \text{for every finite } \Pi>0.
}
\]

So the equal-normalized branch \(\mathfrak g_c=1\) does **not** occur at finite positive bias.

---

## 2. Singular point-source limit

The same family obeys
\[
\lim_{\Pi\to\infty}\mathfrak g_\Pi=1.
\]
So the equal-normalized branch is recovered only in the singular point-source limit
\[
\boxed{
\Pi\to+\infty.
}
\]

In source-profile language, this is the delta-function mouth limit.

---

## 3. Traction divergence on the equal-normalized branch

From Stage 125,
\[
\widehat T_m(\Pi)
=
\sqrt{\frac{9\Pi}{20\left[1-R_q(\Pi)\mathcal S_q(\Pi)\right]}},
\qquad
R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

As \(\Pi\to\infty\),
\[
\mathfrak g_\Pi\to 1,
\qquad
\mathcal S_q(\Pi)\to 1,
\qquad
R_q(\Pi)\to R_\infty
=
\frac{(1-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}
\approx 0.145454452260420.
\]
Hence
\[
\Sigma_0(\Pi)
=
\frac{\Pi}{1-R_q(\Pi)\mathcal S_q(\Pi)}
\sim
\frac{\Pi}{1-R_\infty},
\]
and therefore
\[
\boxed{
\widehat T_m(\Pi)\sim
\sqrt{\frac{9}{20(1-R_\infty)}}\,\Pi^{1/2}
\approx 0.725669130700713\,\Pi^{1/2}.
}
\]

So the equal-normalized branch requires
\[
\boxed{
\widehat T_m\to\infty.
}
\]

---

## Result

The explicit positive exponential mouth-layer family proves that:

1. \(\mathfrak g_c=1\) is not a finite-bias branch,
2. it is reached only as \(\Pi\to\infty\),
3. and in that same limit the normalized mouth traction diverges.

So the naive equal-normalized mouth-source branch is **not** a regular finite branch of the
explicit mouth-layer dynamics. It is a singular point-source limit.
