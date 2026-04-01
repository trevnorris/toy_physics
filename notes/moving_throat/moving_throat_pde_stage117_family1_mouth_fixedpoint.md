
# Moving-Throat PDE — Stage 117: Family-1 Shell + First Mixed D/N Tube Reduction

## Goal

Specialize the general Stage-116 fixed-point law to the first explicit
moving-throat mouth-layer branch compatible with the carried Family-1 geometry:

- one **static shell-compliance** lane,
- one **first mixed D/N half-wave** lane.

This is the cleanest reduction of the coupled mouth-layer solve that still knows
about the localized-Maxwell mixed tube.

---

## 1. Channel identifications

Take the two diagonal D/N channels to be

- shell/compliance lane:
  \[
  \boxed{\kappa_s=0,}
  \]
  so it contributes the exact static limit \(\mathcal S(\Pi,0)=1\);

- mixed localized-Maxwell lane:
  \[
  \boxed{\kappa_q=\frac{\pi}{2},}
  \]
  the first D/N half-wave on the actual throat span.

Define the corresponding dimensionless gains
\[
M_s:=\frac{L\,H_sG_s}{\Theta_\sigma},
\qquad
M_q:=\frac{L\,H_qG_q}{\Theta_\sigma}.
\]

Then the coupled mouth-layer law collapses to
\[
\boxed{
\Pi = M_s + M_q\,\mathcal S_q(\Pi),
\qquad
\mathcal S_q(\Pi):=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
}
\]

Using the exact Stage-116 kernel,
\[
\boxed{
\mathcal S_q(\Pi)
=
\frac{
\Pi\left[\frac{\pi}{2}\tanh\!\frac{\pi}{2}
+\Pi\left(e^{-\Pi}\operatorname{sech}\!\frac{\pi}{2}-1\right)\right]
}{
(1-e^{-\Pi})\left(\frac{\pi^2}{4}-\Pi^2\right)
}.
}
\]

So the Family-1 mouth-layer branch is now an explicit one-dimensional fixed-point
problem in the two gains \(M_s\) and \(M_q\).

---

## 2. Canonical compensation line

Stage 114 already fixed the exact source-shape compensation point
\[
\Pi_* \approx 1.50882951349316.
\]

Inserting this into the Family-1 fixed-point law gives the exact gain-line
condition
\[
\boxed{
M_s = \Pi_* - M_q\,\mathcal S_q(\Pi_*).
}
\]

Numerically,
\[
\boxed{
\mathcal S_q(\Pi_*) \approx 0.658075937605428.
}
\]
So the canonical Family-1 compensation line is
\[
\boxed{
M_s \approx 1.50882951349316 - 0.658075937605428\,M_q.
}
\]

Interpretation:

- if the mixed localized-Maxwell lane is mouth-localizing (\(M_q>0\)),
  the shell traction demand is reduced;
- if the mixed lane is de-localizing (\(M_q<0\)),
  the shell traction must be larger than the Stage-114 pure-mechanical threshold.

---

## 3. What has improved

Compared to Stage 114, the threshold is no longer written in terms of the effective
combination
\[
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_zA_0\right|_{\rm m}
\]
alone.

It is now resolved into an explicit two-channel boundary-layer law with:

- a static shell lane,
- a first mixed D/N lane,
- and one exact nontrivial response factor
  \(\mathcal S_q(\Pi)\).

So the “actual \(\Pi_m\)” problem has narrowed from an effective slope datum to a
simple two-gain fixed-point problem.

---

## Result

On the first explicit Family-1 coupled mouth-layer branch, the actual source-shape
bias is determined by
\[
\boxed{
\Pi = M_s + M_q\,\mathcal S_q(\Pi),
}
\]
and the canonical compensated branch sits on the straight gain line
\[
\boxed{
M_s \approx 1.50882951349316 - 0.658075937605429\,M_q.
}
\]

The remaining ambiguity is therefore no longer profile selection and no longer a
free mouth slope: it is just the signed gain pair \((M_s,M_q)\).
