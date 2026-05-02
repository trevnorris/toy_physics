# Parent Throat Action — Outgoing-Amplitude Admissibility

## Purpose

`step_24` showed that a target-blind outgoing amplitude coordinate
\(\lambda_{\rm out}\) can move the reduced sampled frontier to much smaller
required static normalization \(\widehat m_0^{\rm req}\).

That is only useful if \(\lambda_{\rm out}\) is a physically admissible
branch-derived coordinate.

This step tests the sharpest conservative interpretation:

> suppose the reduced outgoing amplitude is really the PDE outgoing
> normalization scalar \(N_Q\).

---

## 1. Exact elimination

The script combines the two exact laws

\[
\widehat m_0^{\,2} N_Q P_0^{\rm base}=P_0^{\rm target},
\qquad
\widehat m_0^{\,2}\chi_QN_Q=1.
\]

Writing \(X=\widehat m_0^{\,2}N_Q\), eliminating \(X\) gives

\[
X=\frac{P_0^{\rm target}}{P_0^{\rm base}},
\qquad
\chi_Q=\frac{P_0^{\rm base}}{P_0^{\rm target}},
\qquad
N_Q=\frac{P_0^{\rm target}}{\widehat m_0^{\,2}P_0^{\rm base}}.
\]

The script now checks the elimination by direct residual substitution and pins
the elimination Jacobian determinant to \(P_0^{\rm target}\). It also includes
mutation guards: changing either \(X\) or \(\chi_Q\) breaks the corresponding
law instead of silently passing.

So once the static target is imposed, \(\chi_Q\) is fixed by the unscaled base
branch and is **independent** of the outgoing amplitude bookkeeping parameter.

---

## 2. Concrete frontier points

The script reuses three already-audited frontier points:

- the `step_23` target-blind ray point with
  \[
  (\text{scale},\lambda_{\rm out},\widehat m_0^{\rm req})
  =(0.09,1,194.6081703105869),
  \]
- the `step_24` best sampled point with \(Q_{\rm iso}\le 1\),
  \[
  (0.09,2000,4.351570977913287),
  \]
- the `step_24` best sampled point with \(Q_{\rm iso}\le 0.5\),
  \[
  (0.092,2000,4.7912441136331765).
  \]

From these, the inferred outgoing scalar is

\[
\chi_Q=\frac{1}{\widehat m_0^{\,2}\lambda_{\rm out}}.
\]

The first two points lie at the same frozen scale, and the script verifies

\[
\chi_Q(0.09,\lambda_{\rm out}=1)
=
\chi_Q(0.09,\lambda_{\rm out}=2000)
\approx 2.6404494712422554\times 10^{-5}.
\]

At the smaller sampled defect point one gets

\[
\chi_Q(0.092,\lambda_{\rm out}=2000)
\approx 2.1780778923914126\times 10^{-5}.
\]

If \(\lambda_{\rm out}\) is identified with \(N_Q\), then the outgoing
finish-line defect is simply

\[
N_Q-1=1999
\]

at both large-amplitude `step_24` points.

The same-scale static-load ratio is also checked:

\[
\frac{\widehat m_0(0.09,\lambda=2000)^2\,2000}
{\widehat m_0(0.09,\lambda=1)^2}
=1.
\]

As a mutation guard, replacing \(\lambda_{\rm out}=2000\) by \(1999\) at the
same point leaves an odd-closure residual of about
\(-5.0\times10^{-4}\), so the large outgoing normalization is load-bearing in
the test.

---

## 3. Interpretation

This closes the easy loophole.

If the reduced outgoing-amplitude coordinate is really the PDE
outgoing-normalization scalar \(N_Q\), then:

1. the `step_24` frontier improvement does **not** move \(\chi_Q\) toward the
   Packet-A finish line;
2. it leaves \(\chi_Q\) at order \(10^{-5}\), far from the canonical
   outgoing value \(\chi_Q=1\);
3. the large-amplitude points carry the explicit outgoing defect
   \[
   N_Q-1=1999.
   \]

So `step_24` should be read as an upper-bound reduced-family diagnostic unless
the outgoing amplitude is shown to be a genuinely independent branch-derived
coordinate rather than the PDE scalar \(N_Q\).
