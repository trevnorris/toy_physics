# Moving-Throat PDE — Stage 124: Mouth-Gain Status Update

The coupled mouth-layer problem is now much tighter than it was at Stage 116.

## What is now explicit

1. The actual mouth gains are derived from the explicit throat-core ansatz:
   \[
   M_s=\frac{L g_s^2}{K_s\Theta_\sigma},
   \qquad
   M_q=-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
   \]
   This statement is within the explicit throat-core plus mouth-layer closure of
   Stages 120–123; it is not yet a theorem that the full PDE admits no other
   microscopic realization of the mouth source law.

2. In normalized core variables,
   \[
   M_q=-M_s\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
   \]

3. On the exact core-balance compensation family,
   \[
   M_q=-\frac{M_s}{4},
   \qquad
   M_s=4\Sigma_m,
   \qquad
   M_q=-\Sigma_m.
   \]

4. On the explicit Family-1 branch,
   the natural equal-normalized source law gives
   \[
   M_s\approx 1.66854,
   \qquad
   M_q\approx -0.24270,
   \]
   while the exact outlet-consistent canonical branch gives
   \[
   M_s\approx 1.80594,
   \qquad
   M_q\approx -0.45149.
   \]

5. Under the self-matched mouth susceptibility closure,
   \[
   M_s=\frac{20}{9}\widehat T_m^2,
   \]
   so the natural and canonical branches differ by only about `4%` in the normalized
   traction amplitude.

## What remains open

Within the explicit mouth-layer closure, the ambiguity is no longer a profile ambiguity
and no longer a free gain pair.
It has shrunk to two very concrete microscopic questions:

1. does the real moving-throat mouth source choose the natural equal-normalized branch
   \(\mathfrak g_c\approx1\) or the nearby lower compensated branch,
2. and does the actual mouth traction land at the corresponding `\widehat T_m` value?

That is a much smaller target than the original abstract `\Pi_m` problem.
