# Parent Throat Action — Natural-Source Outgoing Burden

## Purpose

`step_25` showed that if the `step_24` outgoing-amplitude coordinate is simply
the PDE scalar \(N_Q\), then the Packet-A finish line is not rescued.

This step sharpens that statement by carrying the exact PDE natural-source-map
law explicitly:

\[
m_{\hat 0}\to 1,
\qquad
N_Q=\chi_Q^{-1}.
\]

The question here is not whether the reduced frontier looks better after
introducing \(\lambda_{\rm out}\). It does.

The question is:

> does \(\lambda_{\rm out}\) reduce the **exact natural-source outgoing
> burden**, or does it merely shift that burden between
> \(\lambda_{\rm out}\) and \(\widehat m_0^{\rm req}\)?

---

## 1. Exact bridge

From the reduced static-normalized packet one has

\[
\widehat m_0^{\,2}\lambda_{\rm out}P_0^{\rm base}=P_0^{\rm target}.
\]

On the exact PDE natural-source branch,

\[
N_Q=\chi_Q^{-1},
\qquad
m_{\hat 0}=1.
\]

So the same base branch would require

\[
\chi_Q=\frac{P_0^{\rm base}}{P_0^{\rm target}}
=\frac{1}{\lambda_{\rm out}\widehat m_0^{\,2}},
\]
\[
N_Q^{\rm nat}=\chi_Q^{-1}
=\lambda_{\rm out}\widehat m_0^{\,2}.
\]

This is the key identity.

At fixed frozen scale, the exact natural-source outgoing burden depends only on
the underlying base branch through \(P_0^{\rm base}\), not on how that burden is
split between \(\lambda_{\rm out}\) and \(\widehat m_0\).

---

## 2. Same-scale invariance

The script compares two already-audited points at the same frozen scale
\(0.09\):

- the `step_23` point
  \[
  (\lambda_{\rm out},\widehat m_0^{\rm req})=(1,194.6081703105869),
  \]
- the `step_24` best point with \(Q_{\rm iso}\le 1\)
  \[
  (\lambda_{\rm out},\widehat m_0^{\rm req})=(2000,4.351570977913287).
  \]

Despite the dramatic change in the reduced bookkeeping variables, the exact
natural-source burden is identical:

\[
N_Q^{\rm nat}
=\lambda_{\rm out}\widehat m_0^{\,2}
\approx 37872.3399516344.
\]

Equivalently,

\[
\chi_Q^{\rm nat}\approx 2.6404494712422554\times 10^{-5}.
\]

So the `step_24` outgoing-amplitude win does **not** reduce the exact
natural-source outgoing burden at that scale. It only repartitions it.

At the `step_24` best sampled point with \(Q_{\rm iso}\le 0.5\),

\[
N_Q^{\rm nat}\approx 45912.04031284913,
\qquad
\chi_Q^{\rm nat}\approx 2.1780778923914126\times 10^{-5}.
\]

That is worse, not better.

---

## 3. Minimal Robin outlet interpretation

The PDE ledger’s minimal Robin outlet model gives

\[
\chi_Q^{R}=\frac{3}{3-\rho},
\qquad
\rho=3-3N_Q.
\]

Feeding in the natural-source burdens above yields

\[
\rho\approx -113614.01985490319
\]

for the scale-\(0.09\) points, and

\[
\rho\approx -137733.1209385474
\]

for the best sampled \(Q_{\rm iso}\le 0.5\) point.

So in the simplest exact outlet-deformation family, realizing these tiny
\(\chi_Q\) values would demand an enormous negative Robin shift, far outside any
plausible perturbative or near-canonical interpretation.

---

## 4. Interpretation

This is the stronger admissibility verdict.

1. The `step_24` outgoing-amplitude frontier improves the reduced
   normalization tradeoff.
2. But under the exact PDE natural-source map, that improvement does not reduce
   the actual outgoing burden.
3. At fixed frozen scale it only trades a large \(\widehat m_0\) burden for a
   large \(\lambda_{\rm out}\) burden.
4. In the minimal Robin outlet family, the required \(\chi_Q\) values correspond
   to enormous negative \(\rho\).

So the current `step_24` coordinate is still not a realized outgoing-branch
success. To make it physical, we would need a genuinely independent
branch-derived transfer/outlet coordinate, not just a rescaling that the exact
PDE source map reabsorbs.
