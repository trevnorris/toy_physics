# Moving-Throat PDE — Stage 121: Geometric Selection of the Core Hybridization Ratio

## Goal

Turn the normalized core-hybridization parameter
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}}
\]
from a free reduced number into an explicit branch value.

## Parent motivation

The parent 4D stack already keeps the throat geometry explicit through the collective throat variables `(a,L)` in the confinement / geometry sector, while the localized-Maxwell/plasma papers keep the mixed channels
\[
A_w,\qquad F_{\mu w},\qquad J^w
\]
in the microscopic ontology outside the strict far-field brane limit. fileciteturn33file4 fileciteturn33file1turn33file3
The carried constructive hierarchy also keeps the preferred aspect ratio
\[
L/a \approx 1.85.
\]
fileciteturn35file0

So the first concrete branch test is to identify the auxiliary mixed D/N tube of Stage 99 with the **actual throat axial span**:
\[
\boxed{L_W=L.}
\]

This is the simplest honest moving-throat core identification: the mixed side-channel is not a detached auxiliary cavity but the same axial throat corridor already selected by the carried geometry branch.

## 1. Exact geometric branch law

Stage 99 already fixed the compensation-selected D/N tube length to
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
\]
Imposing \(L_W=L\) gives the exact geometric hybridization branch
\[
\boxed{
\mathfrak r_{\rm geom}\!\left(\frac{L}{a}\right)
=
\sqrt{\frac{12}{\pi^2}\left(\frac{L}{a}\right)^2-1}.
}
\]

So the normalized shell/mixed hybridization is no longer free once the mixed tube is identified with the actual throat length.

The existence condition is
\[
\boxed{
\frac{L}{a}\ge \frac{\pi}{2\sqrt3}\approx 0.9069,
}
\]
which is comfortably satisfied on the carried throat branch.

## 2. Family-1 / preferred-aspect-ratio value

Using the carried preferred aspect ratio
\[
\frac{L}{a}=\frac{37}{20},
\]
the geometric branch value is
\[
\boxed{
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
=
\frac{\sqrt{4107-168\pi^2}}{10\pi}
\approx 1.77799353547498.
}
\]

Equivalently,
\[
\boxed{
r_c^{F1}=\mathfrak r_{F1}^2\approx 3.16126101219081.
}
\]

So the first concrete moving-throat core solve already fixes the static/mixed hybridization to an `O(1)` value rather than a tunable small parameter.

## 3. Immediate mixed-tube consequence

With \(L_W=L\), the auxiliary mixed-tube pole is simply
\[
\boxed{
\Omega_W=\frac{\pi c_s}{2L},
}
\]
the first D/N half-wave of the actual throat corridor. So the mixed side-channel is now tied directly to the real throat geometry, not to a detached bookkeeping cavity.

## Result

The surviving outlet-core problem has already shrunk:

- before Stage 223, the normalized hybridization \(\mathfrak r\) was free;
- after identifying the mixed D/N tube with the actual throat length, it is fixed to
  \[
  \mathfrak r_{F1}\approx 1.778
  \]
  on the preferred aspect-ratio branch.

So the next question is no longer “what is \(\mathfrak r\)?” It is “what normalized mouth-coupling ratio \(\mathfrak g\) does the actual core pick at this geometrically fixed \(\mathfrak r\)?”
