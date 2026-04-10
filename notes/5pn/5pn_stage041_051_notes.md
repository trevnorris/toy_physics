# 5PN stages 41–51 notes

This bundle continues the post-Stage-40 support/normalization program by turning the
physical placement map into an operator-selected branch law and then projecting the
remaining theorem gap back into explicit parent-action overlap data.

## Stage 41 — coupled support/source operator and exact `Pe` branch equation

The first coupled axial operator is

a) source transport
\[
\partial_t\sigma + \partial_s J = 0,
\qquad
J = -D_\sigma \partial_s\sigma + v_\sigma \sigma,
\]

b) support field
\[
-T_X \partial_s^2 \phi + K_X \phi = \Lambda_\phi \sigma,
\]
with support boundary conditions
\[
T_X \phi_s(0)=K_m\phi(0),
\qquad
\phi_s(L)=0.
\]

After nondimensionalization,
\[
\kappa = K_X L^2/T_X,
\qquad
\eta = K_m L/T_X,
\qquad
Pe = v_\sigma L/D_\sigma,
\qquad
\Xi = \mu_\sigma \Lambda_\phi^2 L^2/(D_\sigma T_X).
\]

On the stationary zero-flux branch,
\[
\Sigma_{Pe}(x)=\frac{Pe\,e^{Pe x}}{e^{Pe}-1},
\qquad x=s/L\in[0,1].
\]
The exact support-drop kernel is
\[
K_{\kappa,\eta}(x)=\frac{\cosh(\alpha x)+(\eta/\alpha)\sinh(\alpha x)-\cosh(\alpha(1-x))}{\alpha\sinh\alpha+\eta\cosh\alpha},
\qquad \alpha=\sqrt\kappa,
\]
with derivative
\[
\frac{dK_{\kappa,\eta}}{dx}=
\frac{\alpha\sinh(\alpha x)+\eta\cosh(\alpha x)+\alpha\sinh(\alpha(1-x))}{\alpha\sinh\alpha+\eta\cosh\alpha}>0.
\]

So the exact dimensionless support drop is
\[
\Delta(Pe;\kappa,\eta)=\int_0^1 dx\;K_{\kappa,\eta}(x)\Sigma_{Pe}(x),
\]
with endpoint values
\[
\Delta_0(\kappa,\eta)=
\frac{\eta(\cosh\alpha-1)}{\alpha^2(\alpha\sinh\alpha+\eta\cosh\alpha)},
\]
\[
\Delta_\infty(\kappa,\eta)=
\frac{\cosh\alpha+(\eta/\alpha)\sinh\alpha-1}{\alpha\sinh\alpha+\eta\cosh\alpha}.
\]

The branch point is therefore selected by the exact fixed-point equation
\[
Pe = \Xi\,\Delta(Pe;\kappa,\eta),
\]
and every constructive branch root obeys
\[
\Xi\Delta_0(\kappa,\eta)
\le Pe_* \le
\Xi\Delta_\infty(\kappa,\eta).
\]
At weak coupling,
\[
Pe_* = \Xi \Delta_0(\kappa,\eta)+O(\Xi^2).
\]

## Stage 42 — exact residual bounds on the operator-selected branch

Evaluating the Stage-40 physical support ratio on the branch root gives
\[
\zeta_{\rm phys}(\Xi,\eta;\kappa)=
\Omega_{Pe_*}^2\,\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2},
\qquad y\tan y=\eta.
\]
Since \(\Omega_{Pe}\) is strictly increasing, the Stage-41 branch interval gives the exact support bracket
\[
\zeta_-(\Xi,\eta;\kappa)
\le \zeta_{\rm phys}(\Xi,\eta;\kappa)
\le \zeta_+(\Xi,\eta;\kappa),
\]
where
\[
\zeta_-=
\Omega_{\Xi\Delta_0}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2},
\qquad
\zeta_+=
\Omega_{\Xi\Delta_\infty}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2}.
\]

Inside the Stage-40 reachability window, define the unique constructive point \(Pe_{\rm req}\) by
\[
\Omega_{Pe_{\rm req}}^2=
\zeta_{\rm req}\,\frac{\kappa+y^2}{\kappa+\pi^2/4}.
\]
Then the exact coupling thresholds are
\[
\Xi_{\rm fail}=\frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)},
\qquad
\Xi_{\rm suff}=\frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\qquad
\Xi_{\rm fail}\le \Xi_{\rm suff}.
\]
So the branch has a sharp three-zone structure:

- \(\Xi\le\Xi_{\rm fail}\): impossible,
- \(\Xi\ge\Xi_{\rm suff}\): guaranteed,
- \(\Xi_{\rm fail}<\Xi<\Xi_{\rm suff}\): only here is the full root solve needed.

The exact residual envelope is
\[
R_- \le R_{\rm phys} \le R_+,
\qquad
R_{\rm phys}=\zeta_{\rm req}-\zeta_{\rm phys},
\]
with
\[
R_- = \zeta_{\rm req}-\zeta_+,
\qquad
R_+ = \zeta_{\rm req}-\zeta_-.
\]

Using
\[
\Omega_{Pe}^2 = 1 + \frac{4-\pi}{\pi}Pe + O(Pe^2),
\]
one gets the weak-coupling branch law
\[
\zeta_{\rm phys}=
A_K(\eta;\kappa)
\Bigl[1+\frac{4-\pi}{\pi}\,\Xi\Delta_0(\kappa,\eta)+O(\Xi^2)\Bigr].
\]

## Stage 43 — entropic source microclosure and microscopic gain

The first explicit source/support free energy is
\[
F[\sigma,\phi]=
\int_0^L ds\Bigl[
\Theta_\sigma\sigma(\log(\sigma/\sigma_*)-1)-\Lambda_\phi\sigma\phi
+\frac{T_X}{2}\phi_s^2+\frac{K_X}{2}\phi^2
\Bigr]+
\frac{K_m}{2}\phi(0)^2.
\]
Its exact variations are
\[
\mu_\sigma^{\rm chem}=
\frac{\delta F}{\delta\sigma}=
\Theta_\sigma\log(\sigma/\sigma_*)-\Lambda_\phi\phi,
\]
\[
-T_X\phi_{ss}+K_X\phi=\Lambda_\phi\sigma,
\qquad
T_X\phi_s(0)=K_m\phi(0),
\qquad
\phi_s(L)=0.
\]

The minimal positive-density Onsager current is
\[
J=-M_\sigma\sigma\partial_s\mu_\sigma^{\rm chem}
  =-D_\sigma\partial_s\sigma + M_\sigma\Lambda_\phi\sigma\partial_s\phi,
\]
with exact Einstein relation
\[
D_\sigma=M_\sigma\Theta_\sigma.
\]

Under the affine-drop reduction
\[
\phi(s)\approx \phi(0)+[\Delta\phi]s/L,
\]
the stationary zero-flux branch gives the exact exponential family
\[
\sigma(s)=C\exp\!igl[(\Lambda_\phi\Delta\phi)/(\Theta_\sigma L)\,s\bigr],
\]
so
\[
Pe=(\Lambda_\phi/\Theta_\sigma)\Delta\phi.
\]
Using
\[
\Delta\phi=(\Lambda_\phi L^2/T_X)\Delta(Pe;\kappa,\eta),
\]
one gets the exact microscopic coupling
\[
\Xi_{\rm micro}=
\frac{\Lambda_\phi^2 L^2}{\Theta_\sigma T_X}
=
\chi_\sigma\frac{\Lambda_\phi^2 L^2}{T_X},
\qquad
\chi_\sigma:=1/\Theta_\sigma.
\]

The closure is automatically passive:
\[
\frac{dF}{dt}=-\int_0^L ds\;\frac{J^2}{M_\sigma\sigma}\le 0
\]
under no-flux boundaries.

## Stage 44 — microscopic gain thresholds and exact phase diagram

Using \(\kappa=K_XL^2/T_X\), the Stage-43 coupling becomes
\[
\Xi_{\rm micro}=\kappa G_{\rm micro},
\qquad
G_{\rm micro}:=\chi_\sigma\Lambda_\phi^2/K_X.
\]
So the exact microscopic thresholds are
\[
G_{\rm fail}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_\infty(\kappa,\eta)},
\qquad
G_{\rm suff}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_0(\kappa,\eta)}.
\]
The exact phase diagram is:

- \(G_{\rm micro}\le G_{\rm fail}\): impossible,
- \(G_{\rm micro}\ge G_{\rm suff}\): guaranteed,
- only the bounded interval in between needs the full root solve.

Equivalent threshold surfaces are
\[
\chi_\sigma \le \frac{T_X Pe_{\rm req}}{\Lambda_\phi^2 L^2 \Delta_\infty}\quad\Rightarrow\quad\text{fail},
\]
\[
\chi_\sigma \ge \frac{T_X Pe_{\rm req}}{\Lambda_\phi^2 L^2 \Delta_0}\quad\Rightarrow\quad\text{succeed},
\]
and similarly for \(\Lambda_\phi^2\).

Soft-support limit:
\[
\Delta_0\to \frac12,
\qquad
\Delta_\infty\to 1,
\qquad
G_{\rm fail}\sim \frac{Pe_{\rm req}}{\kappa},
\qquad
G_{\rm suff}\sim \frac{2Pe_{\rm req}}{\kappa}.
\]
So very soft support is strongly disfavored.

Highly compliant-mouth limit:
\[
\Delta_0^{(\infty)}=
\frac{1-\operatorname{sech}(\sqrt\kappa)}{\kappa},
\qquad
\Delta_\infty^{(\infty)}=
\frac{\tanh(\sqrt\kappa)}{\sqrt\kappa},
\]
so
\[
G_{\rm fail}^{(\infty)}=
\frac{Pe_{\rm req}}{\sqrt\kappa\tanh(\sqrt\kappa)},
\qquad
G_{\rm suff}^{(\infty)}=
\frac{Pe_{\rm req}}{1-\operatorname{sech}(\sqrt\kappa)}.
\]
For \(\kappa\gg1\), these reduce to
\[
G_{\rm fail}^{(\infty)}\sim \frac{Pe_{\rm req}}{\sqrt\kappa},
\qquad
G_{\rm suff}^{(\infty)}\sim Pe_{\rm req}.
\]

## Stage 45 — parent-action projection of the microscopic gain

Starting from the parent matter energy
\[
H_\psi=
\int d^4X\left[\frac{\hbar^2}{2m}|D_i\psi|^2+V_{\rm conf}\rho+U(\rho)\right],
\]
with frozen EOS
\[
U(\rho)=K\rho^5/4,
\qquad
h(\rho)=\frac{5K}{4}\rho^4,
\qquad
h'(\rho)=5K\rho^3=\frac{m c_s^2(\rho)}{\rho},
\]
the local compressional quadratic energy is
\[
\delta H_{\rm comp}=
\frac12\int d^4X\;h'(\rho_*)(\delta\rho)^2.
\]

Project one source channel
\[
\delta\rho(s,y)=\sigma(s)\chi_\sigma(y)
\]
and one support channel entering the confinement as
\[
\delta V_{\rm conf}(s,y)=-g_\phi\chi_\phi(y)\phi(s).
\]
Then the exact reduced coefficients are
\[
\Theta_\sigma=h'(\rho_*)N_{\sigma\sigma},
\qquad
\Lambda_\phi=g_\phi O_{\sigma\phi},
\]
with parent overlap invariants
\[
N_{\sigma\sigma}=\int d^3y\,\chi_\sigma^2,
\qquad
N_{\phi\phi}=\int d^3y\,\chi_\phi^2,
\qquad
O_{\sigma\phi}=\int d^3y\,\chi_\sigma\chi_\phi.
\]

So the effective source susceptibility is
\[
\chi_\sigma^{\rm eff}=
\frac{1}{\Theta_\sigma}=
\frac{\rho_*}{m c_{s,*}^2 N_{\sigma\sigma}},
\]
and the microscopic gain becomes the explicit parent quantity
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 O_{\sigma\phi}^2}{m c_{s,*}^2 K_X N_{\sigma\sigma}}.
\]
Introducing the coherence factor
\[
C_{\sigma\phi}^2=
\frac{O_{\sigma\phi}^2}{N_{\sigma\sigma}N_{\phi\phi}},

aud via Cauchy–Schwarz, one gets the exact factorization
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
C_{\sigma\phi}^2,
\qquad 0\le C_{\sigma\phi}^2\le 1.
\]

## Stage 46 — parent-overlap threshold theorem

Combining the parent gain with the Stage-44 phase diagram gives exact parent thresholds
\[
g_{\phi,\rm fail}^2=
\frac{m c_{s,*}^2 K_X N_{\sigma\sigma} G_{\rm fail}}{\rho_* O_{\sigma\phi}^2},
\qquad
g_{\phi,\rm suff}^2=
\frac{m c_{s,*}^2 K_X N_{\sigma\sigma} G_{\rm suff}}{\rho_* O_{\sigma\phi}^2}.
\]
Equivalently, in coherence form,
\[
C_{\rm fail}^2=
\frac{m c_{s,*}^2 K_X G_{\rm fail}}{\rho_* g_\phi^2 N_{\phi\phi}},
\qquad
C_{\rm suff}^2=
\frac{m c_{s,*}^2 K_X G_{\rm suff}}{\rho_* g_\phi^2 N_{\phi\phi}}.
\]

There is an exact Cauchy no-go theorem:
if
\[
G_{\max}(g_\phi):=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
< G_{\rm fail}(\kappa,\eta),
\]
then no profile engineering of \(\chi_\sigma\) can rescue the branch.

Inserting
\[
G_{\rm fail}=
\frac{Pe_{\rm req}}{\kappa\Delta_\infty},
\qquad
G_{\rm suff}=
\frac{Pe_{\rm req}}{\kappa\Delta_0},
\qquad
\kappa=K_XL^2/T_X,
\]
one finds exact amplitude thresholds
\[
g_{\phi,\rm fail}^2=
\frac{m c_{s,*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2 \Delta_\infty},
\qquad
g_{\phi,\rm suff}^2=
\frac{m c_{s,*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2 \Delta_0}.
\]
So \(K_X\) cancels from the explicit prefactor and survives only through the geometry-shape functions \(\Delta_0,\Delta_\infty\).

## Stage 47 — parent equilibrium source/support alignment

The parent equilibrium law
\[
H(y)\,\delta\rho(s,y)+\delta V_{\rm conf}(s,y)=0,
\qquad H(y):=h'(\rho_*(y)),
\]
forces the aligned source profile
\[
\chi_\sigma(y)=g_\phi\chi_\phi(y)/H(y).
\]
So the overlap invariants become
\[
O_{\sigma\phi}=g_\phi I_1,
\qquad
N_{\sigma\sigma}=g_\phi^2 I_2,
\]
with
\[
I_1=\int d^3y\;\frac{\chi_\phi(y)^2}{H(y)},
\qquad
I_2=\int d^3y\;\frac{\chi_\phi(y)^2}{H(y)^2}.
\]
Therefore
\[
C_{\sigma\phi}^2=
\frac{I_1^2}{N_{\phi\phi} I_2}\le 1.
\]
In the thin active layer where \(H(y)\approx H_w\) is nearly constant,
\[
I_1=N_{\phi\phi}/H_w,
\qquad
I_2=N_{\phi\phi}/H_w^2,
\qquad
C_{\sigma\phi}^2=1.
\]
So the matched-layer branch is not arbitrary; it is the natural thin-layer limit of the parent equilibrium branch.

The exact eliminated-source support softening is
\[
\Delta K_X^{\rm (eq)}=g_\phi^2 I_1,
\qquad
G_{\rm eq}=\Delta K_X^{\rm(eq)}/K_X = g_\phi^2 I_1/K_X.
\]

## Stage 48 — explicit thin-wall confinement branch

For the explicit wall family
\[
V_{\rm conf}(r;a)=V_0 f\bigl((r-a)/\ell\bigr),
\]
with wall coordinate \(\xi=(r-a)/\ell\), a support displacement \(a\to a+\phi(s)\) gives
\[
\delta V_{\rm conf} = +\frac{V_0}{\ell} f'(\xi)\phi(s),
\]
so the parent loading amplitude is exactly
\[
g_\phi=V_0/\ell.
\]

The shell integral entering the equilibrium gain is
\[
I_1 = 4\pi\ell\bigl[a^2J_1+2a\ell J_2+\ell^2J_3\bigr],
\]
where
\[
J_n := \int d\xi\;\frac{\xi^n f'(\xi)^2}{H(\xi)}.
\]
For a centered symmetric wall layer, \(J_2=0\), so
\[
I_1 = 4\pi\ell\bigl[a^2J_1+\ell^2J_3\bigr].
\]
The exact equilibrium gain is
\[
G_{\rm eq}=4\pi V_0^2\Bigl[\frac{a^2J_1}{\ell}+2aJ_2+\ell J_3\Bigr]/K_X.
\]
In the thin-wall limit \(\ell\ll a\), the leading gain is
\[
G_{\rm eq}^{\rm(tw)}=\frac{4\pi a^2 V_0^2 J_1}{K_X\ell}.
\]

Comparing with the Stage-44 thresholds gives wall-amplitude surfaces
\[
V_{0,\rm fail}^2=\frac{K_X\ell G_{\rm fail}}{4\pi a^2J_1},
\qquad
V_{0,\rm suff}^2=\frac{K_X\ell G_{\rm suff}}{4\pi a^2J_1}.
\]
After inserting
\[
G_{\rm fail}=
\frac{Pe_{\rm req}}{\kappa\Delta_\infty},
\qquad
G_{\rm suff}=
\frac{Pe_{\rm req}}{\kappa\Delta_0},
\qquad
\kappa=K_XL^2/T_X,
\]
the explicit \(K_X\) factor cancels:
\[
V_{0,\rm fail}^2=
\frac{T_X\ell Pe_{\rm req}}{4\pi a^2L^2J_1\Delta_\infty},
\qquad
V_{0,\rm suff}^2=
\frac{T_X\ell Pe_{\rm req}}{4\pi a^2L^2J_1\Delta_0}.
\]

For an almost constant-compressibility active wall layer, \(H(\xi)\approx H_w\),
\[
J_1 = I_f/H_w,
qquad I_f:=\int d\xi\;f'(\xi)^2,
\]
so the thresholds reduce to
\[
V_{0,\rm fail}^2=
\frac{H_w T_X\ell Pe_{\rm req}}{4\pi a^2L^2I_f\Delta_\infty},
\qquad
V_{0,\rm suff}^2=
\frac{H_w T_X\ell Pe_{\rm req}}{4\pi a^2L^2I_f\Delta_0}.
\]

## Stage 49 — dimensionless wall figure of merit

For the same thin-wall matched branch, define
\[
W_{\rm wall}:=
\frac{4\pi a^2L^2J_1V_0^2}{T_X\ell}.
\]
Since \(\kappa=K_XL^2/T_X\), this is exactly
\[
W_{\rm wall}=\kappa G_{\rm eq}^{\rm(tw)}.
\]
The Stage-44 operator theorem therefore becomes
\[
W_{\rm fail}=\frac{Pe_{\rm req}}{\Delta_\infty(\kappa,\eta)},
\qquad
W_{\rm suff}=\frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\]
with exact wall-form theorem:

- \(W_{\rm wall}\le W_{\rm fail}\): fail,
- \(W_{\rm wall}\ge W_{\rm suff}\): succeed,
- only the narrow intermediate band still needs the full root solve.

If \(H(\xi)\approx H_w\), the wall control variable becomes
\[
W_H=
\frac{4\pi a^2L^2I_fV_0^2}{H_w T_X\ell}.
\]
So the explicit parent branch is now controlled by one dimensionless figure of merit rather than a diffuse set of amplitudes.

## Stage 50 — sech–Gaussian coherence resonance benchmark

For the explicit independent-profile benchmark
\[
\chi_\sigma(y)=\operatorname{sech}(y/w_f),
\qquad
\chi_\phi(y)=e^{-y^2/w_g^2},
\qquad
r:=w_g/w_f,
\]
the exact norms are
\[
N_{\sigma\sigma}=2w_f,
\qquad
N_{\phi\phi}=w_g\sqrt{\pi/2}.
\]
The overlap is
\[
O_{\sigma\phi}=w_f I(r),
\qquad
I(r):=\int_{-\infty}^{\infty} dx\;\operatorname{sech}(x)e^{-x^2/r^2}.
\]
So the coherence is
\[
C^2(r)=\frac{I(r)^2}{r\sqrt{2\pi}}.
\]

Parseval/Fourier-transform arguments give the exact duality
\[
I(r)=\frac{r}{\sqrt\pi}I(\pi/r),
\qquad
C^2(r)=C^2(\pi/r).
\]
Hence the self-dual stationary point is
\[
r_* = \sqrt\pi.
\]
Numerically,
\[
C_{\rm res}^2:=C^2(\sqrt\pi)=0.9944188364515293487\ldots,
\]
so the resonance penalty factor is
\[
P_{\rm res}:=1/C_{\rm res}^2=1.0056124877605762169\ldots.
\]
Thus the best independent sech–Gaussian mismatch branch misses the ideal matched-layer coherence only by about \(0.56\%\).

## Stage 51 — resonance-corrected thresholds

The general parent gain is
\[
G_{\rm micro}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}
C_{\sigma\phi}^2.
\]
On the Stage-47 matched equilibrium branch, \(C_{\sigma\phi}^2=1\), so the matched gain is
\[
G_{\rm match}=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s,*}^2 K_X}.
\]
Stage 49 repackaged this as the wall figure of merit
\[
W_{\rm wall}=
\frac{4\pi a^2L^2J_1V_0^2}{T_X\ell}.
\]
For the independent sech–Gaussian family,
\[
G_{\rm res}(r)=C^2(r)G_{\rm match},
\qquad
W_{\rm res}(r)=C^2(r)W_{\rm wall}.
\]
Therefore the exact profile-family thresholds are
\[
W_{\rm wall}\le \frac{Pe_{\rm req}}{C^2(r)\Delta_\infty}
\quad\Rightarrow\quad \text{fail},
\]
\[
W_{\rm wall}\ge \frac{Pe_{\rm req}}{C^2(r)\Delta_0}
\quad\Rightarrow\quad \text{succeed}.
\]
At resonance \(r=\sqrt\pi\), this becomes
\[
W_{\rm wall}\le \frac{P_{\rm res} Pe_{\rm req}}{\Delta_\infty}
\quad\Rightarrow\quad \text{fail on the resonance family},
\]
\[
W_{\rm wall}\ge \frac{P_{\rm res} Pe_{\rm req}}{\Delta_0}
\quad\Rightarrow\quad \text{succeed on the resonance family}.
\]
So the explicit independent-profile benchmark modifies the Stage-49 wall thresholds by only the tiny factor \(P_{\rm res}\approx1.00561249\).
