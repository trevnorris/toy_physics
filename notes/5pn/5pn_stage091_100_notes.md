# 5PN / Moving-Throat Continuation — Stages 91–100

This batch continues directly from Stage 90, where the isotropic outgoing normalization factor was reduced to
\[
\chi_Q=\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
\]
The goal of Stages 91–100 is to stop treating \((\beta,\Sigma_0,\Sigma_5)\) as abstract branch labels and push them through explicit outlet/core models.

## Stage 91
Classifies the robustness classes of \(\chi_Q\).

- Pure overall scaling is exactly invisible:
  \[
  \Lambda_2^{\rm def}=S\Lambda_2^{\rm out}
  \quad\Longrightarrow\quad
  \chi_Q=1.
  \]
- Pure scale+argument deformation preserves the already-fixed even moments only on the natural positive branch \(\beta=1\), so it is also harmless.
- A genuine additive isotropic throat-core channel can move \(\chi_Q\) while leaving the lower even moments canonical:
  \[
  \chi_Q=\frac{3(S+9\Sigma_5)}{3S-\Sigma_0}
  \qquad(\beta=1).
  \]
- The exact preservation submanifold is
  \[
  \Sigma_5=\frac{S(1-\beta^5)}{9}-\frac{\Sigma_0}{27}.
  \]

## Stage 92
Linearizes around the canonical outgoing branch:
\[
S=1+\varepsilon s,\qquad
\beta=1+\varepsilon b,\qquad
\Sigma_0=\varepsilon a_0,\qquad
\Sigma_5=\varepsilon a_5.
\]
Then
\[
\chi_Q
=
1+\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)+O(\varepsilon^2).
\]
So the minimal isotropic branch-selection data are the triple
\[
(b,a_0,a_5),
\]
and first-order preservation requires
\[
5b+\frac{a_0}{3}+9a_5=0.
\]

## Stage 93
Introduces the first explicit isotropic geometric outlet:
\[
\Lambda_2^{\rm R}(z)=\Lambda_2^{\rm out}(z)+\rho_R.
\]
The normalized response is
\[
\widehat Y_2^{\rm R}(z)
=
1+\frac{z^2}{9-3\rho_R}
+\frac{4-\rho_R}{9(3-\rho_R)^2}z^4
+i\frac{z^5}{27-9\rho_R}+O(z^6),
\]
so
\[
\chi_Q^{\rm R}=\frac{3}{3-\rho_R}.
\]
A pure Robin core therefore generically shifts both the even branch and the odd normalization.

## Stage 94
Adds the first explicit isotropic hidden mixed side-channel:
\[
\Lambda_2^{\rm mix}(z)
=
\Lambda_2^{\rm out}(z)
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}+O(z^6).
\]
The even-preserving conditions force
\[
\kappa_W=-\frac19,
\qquad
\sigma_W=0.
\]
So a standalone isotropic hidden pole cannot sit on the canonical even branch unless it is absent. Its normalization factor is
\[
\chi_Q^{\rm mix}
=
\frac{3(1-9\sigma_W\gamma_W)}{3+\sigma_W}.
\]

## Stage 95
Combines the Robin core and the mixed side-channel:
\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5}+O(z^6).
\]
Imposing the canonical even branch yields exactly two solutions:
\[
(\rho_R,\kappa_W)=(\sigma_W,0)
\quad\text{or}\quad
(\rho_R,\kappa_W)=\left(4\sigma_W,\frac13\right).
\]
The second is the nontrivial compensated branch. On it,
\[
\chi_Q^{\rm hyb}=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
\]
and canonical odd normalization is preserved iff
\[
\gamma_W=\frac19.
\]
With that value the whole outlet collapses to a pure harmless scale deformation.

## Stage 96
Freezes the outlet-model classification:

1. pure Robin loading generically shifts \(\chi_Q\);
2. a standalone isotropic mixed pole is too rigid and cannot preserve the canonical even branch unless it vanishes;
3. the hybrid Robin–mixed outlet admits one nontrivial compensated canonical-even branch, with
   \[
   \rho_R=4\sigma_W,\qquad \kappa_W=\frac13,\qquad \gamma_W=\frac19
   \]
   on the fully preserved canonical branch.

## Stage 97
Replaces the reduced outlet coefficients by a concrete two-channel core model with internal variables

- \(s(\omega)\): static shell/compliance mode,
- \(q(\omega)\): mixed \(A_w/F_{\mu w}\)-type side-channel mode.

The linear core system
\[
\begin{pmatrix}
K_s & \lambda\\
\lambda & -K_q D_W^{\rm bare}(z)
\end{pmatrix}
\binom{s}{q}
=
u\binom{g_s}{g_q}
\]
gives the exact Schur-complement outlet
\[
\delta\Lambda_{\rm core}(z)
=
\frac{g_s^2}{K_s}
-
\frac{(K_s g_q-\lambda g_s)^2}
{K_s\big(K_sK_q D_W^{\rm bare}(z)+\lambda^2\big)}.
\]
Defining
\[
r_c=\frac{\lambda^2}{K_sK_q},
\]
the reduced Robin–mixed parameters are
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=\frac{(K_s g_q-\lambda g_s)^2}{K_s^2K_q(1+r_c)},
\qquad
\kappa_c=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_c=\frac{\gamma_0}{1+r_c}.
\]

## Stage 98
Determines exactly when the concrete core lands on the compensated canonical branch. The nontrivial compensation condition is
\[
\rho_c=4\sigma_c,\qquad \kappa_c=\frac13,\qquad \gamma_c=\frac19.
\]
This is equivalent to the exact coupling-balance law
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2,
\]
so
\[
g_q=
\frac{g_s}{2K_s}\left(2\lambda\pm\sqrt{K_sK_q+\lambda^2}\right).
\]
The bare mixed channel must itself be a scale-deformed copy of the canonical compact outgoing branch:
\[
\kappa_0=\frac{1+r_c}{3},\qquad
\gamma_0=\frac{1+r_c}{9}.
\]
On that surface,
\[
\delta\Lambda_{\rm core}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9},
\qquad
\sigma_*=\frac{g_s^2}{4K_s},
\]
and the normalized outgoing fingerprint stays exactly canonical.

## Stage 99
Realizes the bare mixed side-channel geometrically as the first D/N half-wave on an auxiliary tube of length \(L_W\):
\[
k_W=\frac{\pi}{2L_W},
\qquad
\Omega_W=\frac{\pi c_s}{2L_W}.
\]
In outlet variables \(z=a\omega/c_s\),
\[
\kappa_0=\frac{4L_W^2}{\pi^2 a^2}.
\]
The compensation condition \(\kappa_0=(1+r_c)/3\) fixes
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+r_c}{3}}.
\]
If the bare mixed outlet is a pure-scale deformation of the canonical compact outgoing \(l=2\) branch,
\[
D_W^{\rm bare}(z)=(1+r_c)\left(1-\frac{z^2}{3}-i\frac{z^5}{9}\right)+O(z^6),
\]
then the hybridization factor is removed exactly and
\[
\kappa_c=\frac13,\qquad \gamma_c=\frac19.
\]

## Stage 100
Freezes the concrete outlet-core status. The surviving PDE-facing question is no longer “is there some deformed outlet?” It is much sharper:

> Does the actual moving-throat core realize the concrete coupling-balance surface together with the auxiliary D/N-tube normalization?

On that surface the effective outlet preserves the canonical normalized outgoing fingerprint exactly.

## Net result of Stages 91–100

The isotropic outgoing-branch ambiguity is no longer open-ended. The first explicit outlet audit shows:

- pure scale deformations are harmless,
- pure Robin loading and a standalone isotropic mixed pole do not by themselves preserve the canonical branch,
- a specific compensated Robin–mixed outlet does preserve it,
- that compensated branch is realized by a concrete two-channel throat-core model,
- and the mixed side-channel can be given a direct finite D/N tube realization.

So the next honest step is Stage 101: extract the concrete core parameters \((K_s,K_q,\lambda,g_s,g_q)\) from a parent-action throat-core ansatz rather than leaving them as reduced variables.
