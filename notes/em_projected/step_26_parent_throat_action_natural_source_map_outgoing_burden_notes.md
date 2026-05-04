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
  (\lambda_{\rm out},\widehat m_0^{\rm req})=(1,194.60817031038846),
  \]
- the `step_24` best point with \(Q_{\rm iso}\le 1\)
  \[
  (\lambda_{\rm out},\widehat m_0^{\rm req})=(2000,4.3515709779088505).
  \]

These records are imported from `step_24`, so this burden check does not keep
its own independent copy of the frontier points.

Despite the dramatic change in the reduced bookkeeping variables, the exact
natural-source burden is identical:

\[
N_Q^{\rm nat}
=\lambda_{\rm out}\widehat m_0^{\,2}
\approx 37872.33995155716.
\]

Equivalently,

\[
\chi_Q^{\rm nat}\approx 2.6404494712476408\times 10^{-5}.
\]

So the `step_24` outgoing-amplitude win does **not** reduce the exact
natural-source outgoing burden at that scale. It only repartitions it.

As a mutation guard, lowering the \(Q_{\rm iso}\le 1\) point from
\(\lambda_{\rm out}=2000\) to \(1999\) changes the natural-source burden by
about \(-18.936169975779194\), so the outgoing load is active in the check.
The symbolic guard is upstream of that concrete point check: the script
perturbs the static-normalized target equation
\[
\widehat m_0^{\,2}\lambda_{\rm out}P_0^{\rm base}=P_0^{\rm target}
\]
to
\[
\widehat m_0^{\,2}\lambda_{\rm out}P_0^{\rm base}
=P_0^{\rm target}+\epsilon,
\]
re-solves for \(P_0^{\rm base}\), rebuilds \(\chi_Q\), and then rebuilds
\(N_Q^{\rm nat}\). The resulting residual is
\[
\Delta N_Q
=-\frac{\epsilon\,\lambda_{\rm out}\widehat m_0^{\,2}}
{P_0^{\rm target}+\epsilon},
\]
so the check probes the upstream solve rather than a cosmetic rearrangement of
\(N_Q=\lambda_{\rm out}\widehat m_0^2\).

At the `step_24` best sampled point with \(Q_{\rm iso}\le 0.5\),

\[
N_Q^{\rm nat}\approx 45912.04031275389,
\qquad
\chi_Q^{\rm nat}\approx 2.1780778923959307\times 10^{-5}.
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
\rho\approx -113614.01985467147
\]

for the scale-\(0.09\) points, and

\[
\rho\approx -137733.12093826168
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
