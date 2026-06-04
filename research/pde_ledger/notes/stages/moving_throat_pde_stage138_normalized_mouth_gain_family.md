# Moving-Throat PDE — Stage 138: Normalized Mouth-Gain Family and Compensation Ratio

## Goal

Rewrite the explicit gain map of Stage 137 in the normalized parent variables already
used for the throat-core compensation family.

This turns the actual mouth gains into one overall amplitude times one exact ratio.

---

## 1. Normalized parent variables

Use the same core-normalized quantities as in Stages 119–123:
\[
\mathfrak r:=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g_c:=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\qquad
\Sigma_0:=\frac{L g_s^2}{K_s\Theta_\sigma}.
\]

Then Stage 137 gives immediately
\[
\boxed{M_s=\Sigma_0.}
\]

For the mixed gain,
\[
K_sg_q-\lambda g_s = g_s\sqrt{K_sK_q}(\mathfrak g_c-\mathfrak r),
\]
so
\[
\boxed{
M_q
=
-\Sigma_0\,
\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

Define the exact mixed-to-shell gain ratio
\[
\boxed{
R_q:=-\frac{M_q}{M_s}=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

So the full coupled Family-1 mouth law is
\[
\boxed{
\Pi = \Sigma_0\Big[1-R_q\,\mathcal S_q(\Pi)\Big].
}
\]

---

## 2. Exact compensation family

Stage 115 already fixed the core-balance family
\[
\mathfrak g_c
=
\mathfrak r\pm \frac12\sqrt{1+\mathfrak r^2}.
\]

Inserting this into the exact ratio gives
\[
\boxed{R_q=\frac14.}
\]

Therefore on the exact compensated branch,
\[
\boxed{
M_s=\Sigma_0,
\qquad
M_q=-\frac{\Sigma_0}{4}.
}
\]
If one defines
\[
\Sigma_m:=\frac{\Sigma_0}{4},
\]
then this is exactly the Stage 135 one-parameter closure
\[
\boxed{M_s=4\Sigma_m,\qquad M_q=-\Sigma_m.}
\]

So the Stage 135 closure is not independent bookkeeping. It is the normalized image of
Stage 137 plus the exact core-balance compensation surface.

---

## Result

The actual mouth gains are completely controlled by one overall amplitude `\Sigma_0`
and one exact ratio
\[
\boxed{R_q=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.}
\]

On the exact compensated branch this collapses to
\[
\boxed{R_q=\frac14,}
\]
so the outlet-consistent mouth closure is derived rather than assumed.
