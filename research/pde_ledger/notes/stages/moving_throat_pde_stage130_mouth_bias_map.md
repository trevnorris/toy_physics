
# Moving-Throat PDE — Stage 232: Exact Mouth-Bias Map and Family-1 Compensation Point

## Goal

Compute the exact Family-1 mouth-bias factor for the explicit boundary-layer
profile
\[
\sigma_{\Pi}(z)=\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})},
\]
and solve the canonical lower-branch compensation condition on the explicit
Family-1 geometry.

---

## 1. Exact mouth-bias factor

The positive mouth-source branch is measured against the first D/N derivative
shape
\[
\cos\!\left(\frac{\pi z}{2L}\right).
\]
So the explicit mouth-bias factor is
\[
\mathfrak g_{\Pi}
=
\int_0^L dz\,
\sigma_{\Pi}(z)\,
\cos\!\left(\frac{\pi z}{2L}\right).
\]

This evaluates exactly to
\[
\boxed{
\mathfrak g_{\Pi}
=
\frac{2\Pi\bigl(2\Pi e^{\Pi}+\pi\bigr)}
{\left(4\Pi^2+\pi^2\right)\left(e^{\Pi}-1\right)}.
}
\]

Equivalent \(x=1/\Pi\) form:
\[
\mathfrak g_{\Pi}
=
\frac{2\left(2+\pi x\,e^{-1/x}\right)}
{\left(4+\pi^2x^2\right)\left(1-e^{-1/x}\right)},
\qquad x=\frac1\Pi,
\]
which matches the earlier truncated-exponential penetration family exactly.

---

## 2. Exact monotonicity law

Because \(\sigma_\Pi(z)\) is an exponential family,
\[
\partial_\Pi \sigma_\Pi
=
-\frac{1}{L}\,\sigma_\Pi
\Bigl(z-\langle z\rangle_\Pi\Bigr),
\]
so
\[
\boxed{
\frac{d\mathfrak g_\Pi}{d\Pi}
=
-\frac{1}{L}\,
\mathrm{Cov}_\Pi\!\left(
\cos\!\frac{\pi z}{2L},
\,z
\right).
}
\]

Since \(\cos(\pi z/2L)\) is strictly decreasing on \([0,L]\), the covariance is
negative, hence
\[
\boxed{
\frac{d\mathfrak g_\Pi}{d\Pi}>0.
}
\]

So the explicit GNLS + localized-Maxwell boundary-layer family is strictly
monotone:
\[
\mathfrak g_{\Pi}: \frac{2}{\pi}\longrightarrow 1
\qquad (\Pi:0^+\to+\infty).
\]

---

## 3. Family-1 compensation point

The explicit Family-1 lower compensated branch is already fixed from the earlier
core balance analysis:
\[
\mathfrak g_-^{F1}\approx 0.758035078944663.
\]

Solving
\[
\mathfrak g_{\Pi}=\mathfrak g_-^{F1}
\]
gives the unique explicit compensation point
\[
\boxed{
\Pi_* \approx 1.50882951349316.
}
\]

Equivalently,
\[
\boxed{
x_*=\frac{1}{\Pi_*}\approx 0.662765402623160,
}
\]
which is exactly the penetration-depth value already found earlier for the
truncated exponential family.

So the compensated Family-1 mouth source is selected by a **moderate**
electrochemical bias:
\[
\Pi_*\sim O(1).
\]

---

## 4. Interpretation

This is a sharper statement than the earlier positive-family result.

It now says:

- the lower compensated branch is not only positive and reachable;
- it is the unique target of a concrete monotone GNLS + localized-Maxwell
  boundary-layer family;
- and the required bias is moderate rather than extreme.

So the remaining branch-selection ambiguity has collapsed to one explicit number:

\[
\boxed{
\Pi_m \stackrel{?}{\approx} 1.50882951349.
}
\]
