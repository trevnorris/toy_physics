# Parent Throat Action — Moderate Branch-B Sector

## Purpose

`step_28` established that the exact Stage-95 branch-B amplitude coordinate
helps primarily by reducing the required normalization, not by improving the
best sampled isotropic defect.

That still leaves an important quantitative question:

> is there a **moderate** branch-B sector where the normalization burden drops
> sharply while the isotropic packet barely moves?

This step answers that by extracting the sampled low-defect Pareto sector under
moderate \(|\sigma|\) budgets.

---

## 1. Low-defect sampled frontiers

The script keeps the same frozen sampled grid as `step_24` and `step_28`, but
restricts attention to the low-defect sector

\[
Q_{\rm iso}\le 1.5.
\]

It then computes the sampled Pareto frontier in the two objectives

\[
(\widehat m_0^{\rm req},Q_{\rm iso})
\]

for the moderate branch-B budgets

\[
|\sigma|\le 19,
\qquad
|\sigma|\le 49,
\]

with

\[
\lambda_{\rm out}=1-\sigma.
\]

The frontiers confirm that the most useful moderate points sit on the
high-scale edge `scale = 0.092`.

---

## 2. The \(\sigma=-19\) point

The undeformed low-defect baseline is

\[
(\text{scale},\lambda_{\rm out},\sigma)
=(0.092,1,0),
\]
\[
\widehat m_0^{\rm req}\approx 214.2709506975902,
\qquad
Q_{\rm iso}\approx 0.26705543084121786.
\]

The same-scale moderate branch-B point is

\[
(\text{scale},\lambda_{\rm out},\sigma)
=(0.092,20,-19),
\]
\[
\widehat m_0^{\rm req}\approx 47.912441136331765,
\qquad
Q_{\rm iso}\approx 0.2670781822626949.
\]

So:

\[
\frac{\widehat m_{0,\rm base}^{\rm req}}{\widehat m_{0,\sigma=-19}^{\rm req}}
=\sqrt{20}\approx 4.47213595499958,
\]

while

\[
\Delta Q_{\rm iso}
\approx 2.2751421477039345\times 10^{-5},
\]
\[
\frac{\Delta Q_{\rm iso}}{Q_{\rm iso}^{\rm base}}
\approx 8.519643312775214\times 10^{-5}.
\]

So in this sampled reduced family, \(\sigma=-19\) cuts the normalization burden
by a factor of `sqrt(20)` while changing the isotropic defect only in the fifth
decimal place.

The corresponding exact transfer-shape amplitude shift is

\[
\delta\ln T_{\rm eff}^2=\ln\lambda_{\rm out}=\ln 20\approx 2.995732273553991.
\]

---

## 3. The \(\sigma=-49\) point

The next moderate point on the same high-scale edge is

\[
(\text{scale},\lambda_{\rm out},\sigma)
=(0.092,50,-49),
\]
\[
\widehat m_0^{\rm req}\approx 30.302488449910886,
\qquad
Q_{\rm iso}\approx 0.26719789465663646.
\]

Now

\[
\frac{\widehat m_{0,\rm base}^{\rm req}}{\widehat m_{0,\sigma=-49}^{\rm req}}
=\sqrt{50}\approx 7.0710678118654755,
\]

while

\[
\Delta Q_{\rm iso}
\approx 1.4246381541860577\times 10^{-4},
\]
\[
\frac{\Delta Q_{\rm iso}}{Q_{\rm iso}^{\rm base}}
\approx 5.334557393551715\times 10^{-4}.
\]

So even `sigma = -49` still leaves the sampled isotropic defect nearly
unchanged while producing a much larger normalization reduction.

The corresponding transfer-shape amplitude shift is

\[
\delta\ln T_{\rm eff}^2=\ln 50\approx 3.912023005428146.
\]

---

## 4. Interpretation

This is the sharpest reduced-branch result so far.

1. The strong `step_24` `lambda_out = 2000` corner is not the only useful
   region.
2. The exact Stage-95 branch-B amplitude channel already has a very favorable
   **moderate** sector.
3. On the sampled frontier, `sigma = -19` to `sigma = -49` gives order-`1`
   logarithmic transfer-shape shifts while leaving the isotropic defect almost
   unchanged and sharply reducing the required normalization.

So the next realization target should no longer be phrased as “can we realize
huge `lambda_out`?”. The sharper question is:

> can a realized branch generate a moderate exact branch-B amplitude deformation
> in the corridor \(\sigma\in[-49,-19]\)?

That is now the most plausible route to a physically acceptable reduced branch.
