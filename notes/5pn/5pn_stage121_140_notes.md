# 5PN continuation — Stages 121–140

This batch covers the next explicit mouth/core chain after Stage 120.

## What is now fixed

Stages 121–123 rewrite the mouth gains in normalized parent variables and collapse the explicit core-to-mouth map to

\[
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2},
\qquad
R_q:=-\frac{M_q}{M_s}=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
\]

On the exact compensation family
\[
\mathfrak g_c=\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2},
\]
that ratio collapses to
\[
R_q=\frac14.
\]
So the Stage-118 outlet-consistent closure is derived rather than assumed. On the self-matched mouth-susceptibility closure the overall shell gain becomes
\[
\Sigma_0=M_s=\frac{20}{9}\,\widehat T_m^2.
\]

Stage 122 evaluates this on the explicit Family-1 branch:

- natural equal-normalized branch:
  \[
  M_s\approx 1.66854,
  \qquad
  M_q\approx -0.24270;
  \]
- exact compensated branch:
  \[
  M_s\approx 1.80594,
  \qquad
  M_q\approx -0.45149.
  \]

The corresponding normalized traction amplitudes differ by only about 4.04%.

## Branch selection

Stages 125–128 close the self-consistent Family-1 mouth branch law,
\[
\Pi=\Sigma_0\bigl[1-R_q(\Pi)\,\mathcal S_q(\Pi)\bigr],
\qquad
R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

The positive exponential mouth family obeys
\[
0<\mathfrak g_\Pi<1
\quad\text{for every finite }\Pi>0,
\]
so the equal-normalized branch \(\mathfrak g_c=1\) is only a singular point-source limit. The upper compensated branch is impossible because \(\mathfrak g_+^{F1}>1\). The lower compensated branch is the unique regular finite branch, reached at
\[
\Pi_*\approx 1.50882951349316,
\qquad
\widehat T_{m,*}\approx 0.901484054174204.
\]

## Mouth rigidity under positive non-exponential deformations

Stages 129–132 derive the first-order rigidity kernel around the canonical point. For any positive normalized deformation, only two source moments matter,
\[
\bar g_\varsigma=\int_0^1\varsigma(x)\cos\!\left(\frac{\pi x}{2}\right)dx,
\qquad
\bar S_\varsigma=\int_0^1\varsigma(x)K_q(x)dx,
\]
and the traction shift is
\[
\delta\widehat T_m
=
\epsilon\Bigl[
A_T(\bar g_\varsigma-\mathfrak g_*)
+
B_T(\bar S_\varsigma-\mathcal S_*)
\Bigr].
\]
The explicit coefficients are
\[
A_T\approx -4.27263956256927,
\qquad
B_T\approx 0.134875005736706,
\]
so overlap changes dominate mixed-kernel changes by a factor
\[
|A_T|/B_T\approx 31.68.
\]

Representative families:

- uniform broadening raises the canonical point,
  \[
  \frac{\delta\Pi_u}{\epsilon}\approx 1.69941,
  \qquad
  \frac{\delta\widehat T_{m,u}}{\epsilon}\approx 0.508756;
  \]
- self-matched derivative sharpening lowers it,
  \[
  \frac{\delta\Pi_d}{\epsilon}\approx -0.382993,
  \qquad
  \frac{\delta\widehat T_{m,d}}{\epsilon}\approx -0.116944.
  \]

So the mouth-side ambiguity is no longer branch choice; it is one explicit rigidity-kernel problem.

## Full-profile mouth correction

Stages 133–136 replace the tangent exponential potential by the full coupled shell + mixed D/N mouth profile. The exact residual
\[
R_*(x)=\Phi_*(x)-\Pi_*x
\]
is tangent matched,
\[
R_*(0)=0,
\qquad
R_*'(0)=0,
\]
and has negative curvature at the mouth,
\[
R_*''(0)=-3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0,
\]
so the actual full profile broadens the source relative to the tangent exponential branch.

Projecting that residual onto the Stage-130 rigidity kernel gives the actual first-order correction:
\[
\delta\mathfrak g_{\rm act}\approx -0.06480697,
\qquad
\delta\mathcal S_{\rm act}\approx -0.03887184,
\]
\[
\delta\Pi_{\rm act}\approx 0.907084,
\qquad
\delta\widehat T_{m,\rm act}\approx 0.271654.
\]
So the corrected canonical point is approximately
\[
\Pi_{\rm corr}\approx 2.41591,
\qquad
\widehat T_{m,\rm corr}\approx 1.17314.
\]
A one-step nonlinear iterate shifts slightly further but in the same direction.

## Full core–mouth co-evolution

Stages 137–140 promote the corrected mouth layer to a fully co-evolving fixed point,
\[
\Sigma(x)=\frac{e^{-\Phi_{\Sigma_0}[\Sigma](x)}}{\int_0^1 e^{-\Phi_{\Sigma_0}[\Sigma](y)}dy},
\qquad
\Phi_{\Sigma_0}[\Sigma](x)=\Sigma_0\bigl[\mathcal T_s[\Sigma](x)-\mathcal R[\Sigma]\mathcal T_q[\Sigma](x)\bigr],
\]
with
\[
\mathcal R[\Sigma]=\frac{(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
\]

At the old canonical traction \(\Sigma_0^*\), the fixed point stays close in bias but drifts off exact compensation,
\[
\mathfrak g_{\rm fp}\approx 0.69336,
\qquad
\mathcal R_{\rm fp}\approx 0.28271,
\qquad
\Pi_{\rm fp}\approx 1.48857.
\]

Restoring exact lower-branch compensation requires a unique renormalized traction. With the reduced numerical fixed-point solver used in this batch, that root is
\[
\Sigma_0^{\rm can}\approx 4.65077,
\qquad
\widehat T_{m,\rm can}\approx 1.44667,
\qquad
\Pi_{\rm can}\approx 3.87134.
\]
These values are very close to the handoff’s quoted ones; the small differences come from the reduced collocation/iteration resolution used in the executable scripts.

## Bottom line after Stage 140

Inside the explicit Family-1 closure:

1. the upper compensated branch is impossible;
2. the equal-normalized branch is singular;
3. the lower compensated branch remains the unique regular canonical branch;
4. full self-consistency preserves that branch only after a finite traction/bias renormalization.

So the mouth/core side is no longer blocked by branch ambiguity. The remaining 2.5PN/5PN theorem gap is the projection of the residual deviation from this renormalized canonical Family-1 point into the outgoing quadrupole-normalization defect.
