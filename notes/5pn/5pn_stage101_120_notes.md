
# 5PN stages 101–120 notes

This bundle continues the chain immediately after the compensated outlet/core result at Stage 100.
The focus of this block is the next honest microscopic gate: replace the reduced outlet/core
coefficients by explicit parent-action overlaps, then carry that data all the way into the
mouth-layer fixed-point problem.

## Stage 101 — parent-action extraction of the core parameters

The reduced two-channel core variables from Stages 97–100 are replaced by explicit overlap
formulas from one concrete GNLS + localized-Maxwell throat-core ansatz.

The shell/compliance mode gives
\[
K_s
=
4\pi a^2\left(
\frac{H_w\ell}{3}
+
\frac{\hbar^2}{15m_\psi\rho_w\,\ell}
\right),
\]
and on the healing-locked shell branch
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\qquad
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
\]

The mixed D/N half-wave gives
\[
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi}{\sqrt2\,L_W^{3/2}}.
\]

The shell–mixed hybridization and shell mouth coupling become
\[
\lambda=-q_* v_{w0}\mathcal I_{sq},
\qquad
\mathcal I_{sq}=\frac{8\sqrt2}{3}a^2\ell\sqrt{L_W},
\qquad
g_s=\mathcal T_m \frac{4\pi a^2\ell}{3}.
\]

So the Stage-97 core matrix is no longer abstract; every entry now has an explicit parent
overlap meaning. This is exactly the next gate identified in the moving-throat notes. fileciteturn30file16

## Stages 102–103 — collapse to normalized parent ratios

The exact core-balance condition
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2
\]
collapses to the two dimensionless ratios
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\]
through
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.
\]

The D/N mixed-tube length becomes
\[
L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
\]

So the surviving outlet/core theorem gate is no longer “what are the reduced coefficients?”
It is only: which branch point \((\mathfrak r,\mathfrak g)\) the actual GNLS + localized-Maxwell
core selects. fileciteturn30file16turn30file5

## Stages 104–107 — explicit Family-1 bridge

To keep the executable chain sequential, I added the missing Family-1 bridge audits that the
later notes assume implicitly.

Using the carried Family-1 geometry \(L/a=37/20\) together with the D/N length law gives
\[
\mathfrak r_{F1}
=
\sqrt{\frac{12}{\pi^2}\left(\frac{37}{20}\right)^2-1}
\approx 1.77799353547498.
\]

The two compensated coupling branches are
\[
\mathfrak g_\pm^{F1}
=
\mathfrak r_{F1}\pm\frac12\sqrt{1+\mathfrak r_{F1}^2},
\]
numerically
\[
\mathfrak g_-^{F1}\approx 0.758035078944663,
\qquad
\mathfrak g_+^{F1}\approx 2.79795199200529.
\]

These bridge scripts also verify the useful ordering
\[
\frac{2}{\pi}<\mathfrak g_-^{F1}<\frac{\pi}{4}<1<\mathfrak g_+^{F1},
\]
which is exactly the window the later positive-source theorem needs.

## Stages 108–111 — positive mouth-source selection

For any positive normalized axial mouth source profile \(\sigma(z)\) on the first D/N interval,
the mouth-bias factor is
\[
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz,
\]
so positivity immediately forces
\[
0\le \mathfrak g[\sigma]\le 1.
\]

Since \(\mathfrak g_+^{F1}>1\), the upper compensated branch is impossible under any positive
localized mouth source, while \(\mathfrak g_-^{F1}\in(0,1)\), so the lower branch is the unique
physically admissible canonical candidate.

The first explicit positive families then show that the lower branch is easy to reach:

- self-matched derivative source:
  \[
  \mathfrak g_{\rm match}=\frac{\pi}{4},
  \]
  only \(3.61\%\) away in traction from exact lower-branch compensation;

- convex derivative/uniform family:
  \[
  \sigma_\xi=(1-\xi)k\cos(kz)+\xi/L,
  \]
  reaches the exact lower branch at
  \[
  \xi_*\approx 0.183918405511540;
  \]

- slab and truncated-exponential penetration laws reach the same branch at
  \[
  x_*^{\rm slab}\approx 0.797839360904564,
  \qquad
  x_*^{\exp}\approx 0.662765402623161.
  \]

So by Stage 111 the branch-selection ambiguity has collapsed: the lower compensated branch is
the unique positive-source branch and is reached with moderate penetration, not by exotic
sign-changing mouth forcing.

## Stages 112–115 — explicit GNLS + localized-Maxwell mouth boundary layer

The abstract positive-source family is replaced by the first explicit boundary-layer law.
With the mouth free energy
\[
F_{\rm mouth}[\sigma]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\!\left(\ln\frac{\sigma}{\sigma_*}-1\right)
+
V_m(z)\sigma
\Big],
\qquad
V_m(z)\approx V_1 z,
\]
and zero-flux Onsager current, the exact normalized source law is
\[
\sigma_\Pi(z)=\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})},
\qquad
\Pi=\frac{V_1L}{\Theta_\sigma}.
\]

Its exact mouth-bias factor is
\[
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}
{(4\Pi^2+\pi^2)(e^\Pi-1)},
\]
with range \(2/\pi\to1\) as \(\Pi:0^+\to\infty\).

Solving \(\mathfrak g_\Pi=\mathfrak g_-^{F1}\) gives the unique canonical Family-1 point
\[
\Pi_* \approx 1.50882951349316.
\]

So the parent threshold is now concrete:
\[
\Pi_m=\frac{L}{\Theta_\sigma}(T_m-q_*A_0')=\Pi_*,
\]
with local sensitivity
\[
\mathfrak g_*' \approx 0.0714453558083196.
\]

At this point the mouth problem is no longer branch sign or family choice. It is the
single microscopic bias question: does the real mouth layer select \(\Pi_m\approx1.51\)?

## Stages 116–120 — coupled mouth fixed point and explicit core-to-mouth gain map

The next honest step is to couple the shell/compliance and mixed Maxwell channels directly in
the mouth layer. The exact scalar D/N response kernel to the exponential source is
\[
\mathcal S(\Pi,\kappa)
=
\frac{
\Pi\left[\kappa\tanh\kappa+\Pi\left(e^{-\Pi}\operatorname{sech}\kappa-1\right)\right]
}{
(1-e^{-\Pi})(\kappa^2-\Pi^2)
},
\]
with \(\mathcal S(\Pi,0)=1\).

So the fully coupled mouth bias obeys
\[
\Pi = \sum_\alpha M_\alpha \mathcal S(\Pi,\kappa_\alpha).
\]

On the first explicit Family-1 branch with one static shell lane and one mixed D/N half-wave,
\[
\kappa_s=0,
\qquad
\kappa_q=\frac{\pi}{2},
\]
this reduces to
\[
\Pi = M_s + M_q \mathcal S_q(\Pi),
\qquad
\mathcal S_q(\Pi)=\mathcal S\!\left(\Pi,\frac{\pi}{2}\right).
\]

At the canonical point,
\[
\mathcal S_q(\Pi_*)\approx 0.658075937605429,
\]
so the exact gain line is
\[
M_s=\Pi_* - M_q \mathcal S_q(\Pi_*).
\]

Imposing the outlet-consistent shell:mixed ratio \(4:-1\) collapses the mouth problem to
\[
\Pi = \Sigma_m[4-\mathcal S_q(\Pi)].
\]
At \(\Pi_*\),
\[
\Sigma_m^*\approx 0.451485277739088,
\qquad
M_s^*\approx 1.80594111095635,
\qquad
M_q^*\approx -0.451485277739088.
\]

Finally, Stage 120 removes the last abstract gain pair entirely. Using the exact Stage-97
Schur complement,
\[
\rho_c=\frac{g_s^2}{K_s},
\qquad
\sigma_c=\frac{(K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)},
\]
the actual mouth-layer gains are
\[
M_s=\frac{L g_s^2}{K_s\Theta_\sigma},
\qquad
M_q=
-\frac{L (K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
\]

So by the end of Stage 120 the mouth fixed-point ambiguity has collapsed from a free source
profile and a free gain pair to one explicit set of parent core quantities. The next clean step
is to normalize these gains and carry them into the self-consistent branch law beyond 120. fileciteturn30file11turn30file16
