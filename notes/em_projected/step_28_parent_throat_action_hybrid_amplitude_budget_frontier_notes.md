# Parent Throat Action — Hybrid-Amplitude Budget Frontier

## Purpose

`step_27` showed that the exact PDE algebra leaves room for a genuine
transfer/outlet amplitude coordinate: the Stage-95 hybrid branch-B family, with

\[
\lambda_{\rm out}=1-\sigma,
\qquad
\chi_Q=1.
\]

That removes the **algebraic** objection to `step_24`.

But it leaves the quantitative one:

> do we need an extremely large branch-B deformation to see any benefit, or do
> modest \(|\sigma|\) budgets already capture most of the reduced frontier gain?

This step answers that by replaying the `step_24` sampled grid under explicit
\(|\sigma|\) budgets.

---

## 1. Frozen branch-B budgets

The script keeps exactly the `step_24` sampled grid

\[
\text{scale}\in\{0,0.03,0.05,0.08,0.088,0.09,0.092\},
\]
\[
\lambda_{\rm out}\in\{1,5,20,50,100,200,500,1000,2000\},
\]

but reinterprets

\[
\lambda_{\rm out}=1-\sigma
\]

through the exact Stage-95 branch-B amplitude coordinate.

It then scans the budget set

\[
|\sigma|\le 4,\ 19,\ 49,\ 199,\ 1999.
\]

So these budgets correspond exactly to allowing

\[
\lambda_{\rm out}\le 5,\ 20,\ 50,\ 200,\ 2000.
\]

---

## 2. Best sampled isotropic defect by \(|\sigma|\) budget

This is the most important surprise in the scan.

For every tested amplitude budget, the best sampled isotropic defect is the
same undeformed point:

\[
(\text{scale},\lambda_{\rm out},\sigma)
=(0.092,1,0),
\]
\[
\widehat m_0^{\rm req}\approx 214.2709506975902,
\qquad
Q_{\rm iso}\approx 0.26705543084121786.
\]

So on this sampled grid, the exact Stage-95 branch-B amplitude channel does
**not** improve the isotropic defect itself.

The best-\(Q_{\rm iso}\) point stays on the undeformed branch for every tested
budget.

---

## 3. Best sampled normalization at \(Q_{\rm iso}\le 1\)

What the amplitude channel *does* improve is the required normalization.

The best sampled normalization among points satisfying \(Q_{\rm iso}\le 1\) is:

\[
|\sigma|\le 4
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out},\sigma)=(0.09,5,-4),
\]
\[
\widehat m_0^{\rm req}\approx 87.03141955826572,
\qquad
Q_{\rm iso}\approx 0.4513189656743675.
\]

\[
|\sigma|\le 19
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out},\sigma)=(0.09,20,-19),
\]
\[
\widehat m_0^{\rm req}\approx 43.51570977913286,
\qquad
Q_{\rm iso}\approx 0.45133756597801233.
\]

\[
|\sigma|\le 49
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out},\sigma)=(0.09,50,-49),
\]
\[
\widehat m_0^{\rm req}\approx 27.521751380184508,
\qquad
Q_{\rm iso}\approx 0.4514417135159099.
\]

\[
|\sigma|\le 199
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out},\sigma)=(0.09,200,-199),
\]
\[
\widehat m_0^{\rm req}\approx 13.760875690092254,
\qquad
Q_{\rm iso}\approx 0.45329746222514433.
\]

\[
|\sigma|\le 1999
\quad\Rightarrow\quad
(\text{scale},\lambda_{\rm out},\sigma)=(0.09,2000,-1999),
\]
\[
\widehat m_0^{\rm req}\approx 4.351570977913287,
\qquad
Q_{\rm iso}\approx 0.618690285150578.
\]

So the exact branch-B amplitude channel buys normalization relief steadily, and
already cuts the required normalization in half by the moderate budget
\(|\sigma|\le 19\), while barely moving \(Q_{\rm iso}\).

---

## 4. Interpretation

This budgeted scan resolves the next ambiguity.

1. The exact Stage-95 branch-B amplitude coordinate does **not** help by
   reducing the sampled isotropic defect.
2. It helps by reducing the required normalization at essentially fixed
   \(Q_{\rm iso}\), at least through moderate budgets.
3. Large budgets continue to reduce \(\widehat m_0^{\rm req}\), but eventually
   start to degrade \(Q_{\rm iso}\) more visibly.

So the branch-realization question is now sharper again:

> can a realized branch produce a moderate branch-B amplitude deformation,
> roughly on the order of \(\sigma\sim -10\) to \(-20\), which already gives a
> substantial normalization gain without materially worsening the isotropic
> packet?

That is the next physically meaningful regime to target.
