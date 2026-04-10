# 5PN stages 34–40 notes

This bundle continues the post-Stage-33 support/normalization program by turning the
`zeta_req` threshold into exact operator criteria on the moving-throat branch.

## Stage 34 — exact lowest-twin sufficiency criterion

Using the tracking-branch functions
\[
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]
\[
F_{\rm tr}(\xi,\delta;R)
=
\frac{\bigl[9\delta+(9+2R^2)\xi\bigr]^2\bigl[9\delta+(9+2R)\xi\bigr]^2}
{81(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2},
\]
the exact product is
\[
\Pi_{\rm tr}=F_{\rm tr}G_{\rm tr}
=
\frac{\xi(\xi+\delta)\bigl[9\delta+(9+2R)\xi\bigr]^2\bigl[9\delta+(9+2R^2)\xi\bigr]}
{9(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2}.
\]

With
\[
C_{\rm mix}:=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\qquad
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}},
\]
the symmetric lowest twin lane is sufficient iff
\[
\Pi_{\rm tr}(\xi_{\rm req},\delta;R_{\rm tr})
\le
\frac{16\Lambda(1-\epsilon)}{\pi^2}.
\]

Equivalent threshold scales:
\[
\Lambda_{\rm twin,req}=\frac{\pi^2}{16(1-\epsilon)}\Pi_{\rm tr},
\qquad
M_{\rm mix}^{(\rm twin,req)}=\frac{G_{\rm tr}}{2},
\]
\[
Z_W^{(\rm twin,req)}
=
\frac{\pi^2(1-\epsilon_\eta)(1-\epsilon)\,G_{\rm tr}}
{16(1+\chi_0)^2}.
\]

The exact twin-saturation depth at fixed mixed baseline is the unique root of
\[
G_{\rm tr}(\xi_{2\times},\delta;R)=2M_{\rm mix},
\]
namely
\[
\xi_{2\times}
=
\frac{2M_{\rm mix}(9+2R^2)-9\delta+\sqrt{(2M_{\rm mix}(9+2R^2)-9\delta)^2+648M_{\rm mix}\delta}}{18}.
\]

## Stage 35 — exact non-twin asymmetry requirement

Define
\[
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}
=
\frac{\Pi_{\rm tr}-C_{\rm mix}}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]
Then
\[
\frac{d\zeta_{\rm req}}{d\Pi_{\rm tr}}
=
\frac{C_{\rm mix}(1-\epsilon)}{\bigl[C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})\bigr]^2}>0,
\]
so the required coherent support ratio grows monotonically with the branch product.

Exact regime split:
\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\Rightarrow\quad
\text{mixed-only enough},
\]
\[
C_{\rm mix}<\Pi_{\rm tr}\le2C_{\rm mix}
\quad\Rightarrow\quad
\text{symmetric lowest twin enough},
\]
\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\Rightarrow\quad
\text{non-twin asymmetry required}.
\]

The exact excess beyond the symmetric twin branch is
\[
\Delta_\zeta:=\zeta_{\rm req}-1
=
\frac{(1-\epsilon)(\Pi_{\rm tr}-2C_{\rm mix})}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]

For a general lowest support lane,
\[
\zeta_0^{(\rm phys)}
=
\frac{K_W^{(\rm eff)}}{K_{\phi,0}^{(\rm eff)}}\Omega_0^2.
\]
So the two equivalent exact rescue thresholds are
\[
\Omega_0^2 \ge \zeta_{\rm req}\frac{K_{\phi,0}^{(\rm eff)}}{K_W^{(\rm eff)}},
\qquad
K_{\phi,0}^{(\rm eff)} \le K_W^{(\rm eff)}\frac{\Omega_0^2}{\zeta_{\rm req}}.
\]

## Stage 36 — exact overlap-boost window

For the D/N lowest support mode
\[
\chi_0(s)=\sqrt{\frac{2}{L}}\sin\!\frac{\pi s}{2L},
\qquad
I_W=\int_0^L\chi_0(s)\,ds=\frac{2\sqrt{2L}}{\pi},
\]
and the normalized exponential source family
\[
\sigma_\alpha(s)=\frac{\alpha e^{\alpha s/L}}{e^\alpha-1},
\qquad
\int_0^L \sigma_\alpha(s)\,ds=L,
\]
the overlap boost is
\[
\Omega_{\exp}(\alpha)
=
\frac{\int_0^L \sigma_\alpha(s)\chi_0(s)\,ds}{I_W}
=
\frac{\pi\alpha\bigl(2\alpha e^\alpha+\pi\bigr)}
{(4\alpha^2+\pi^2)(e^\alpha-1)}.
\]

Exact endpoint values:
\[
\Omega_{\exp}(0)=1,
\qquad
\lim_{\alpha\to+\infty}\Omega_{\exp}(\alpha)=\frac{\pi}{2}.
\]
Therefore
\[
0\le \Omega_0\le\frac{\pi}{2},
\qquad
A_I:=\Omega_0^2\le\frac{\pi^2}{4}.
\]

So pure overlap rescue alone is possible only if
\[
\zeta_{\rm req}\le\frac{\pi^2}{4}.
\]

## Stage 37 — Robin-compliance softening

Replacing the Dirichlet mouth by a Robin compliance gives the lowest-lane eigenvalue condition
\[
y\tan y=\eta,
\qquad
0<y<\frac{\pi}{2},
\]
with
\[
\eta:=hL.
\]

If
\[
x:=\frac{\pi^2 T_X}{L^2K_W^{(\rm eff)}},
\qquad
0<x<4,
\]
then the exact support-softening factor is
\[
A_K(\eta)
=
\frac{K_W^{(\rm eff)}}{K_{\phi,0}^{(\rm eff)}}
=
\frac{1}{1-x/4+xy^2/\pi^2}.
\]

Endpoint window:
\[
A_K=\;1\;\text{at}\;y=\frac{\pi}{2},
\qquad
A_{K,\max}=\frac{4}{4-x}\;\text{at}\;y\to0.
\]

So pure support softening alone can rescue the Stage-35 threshold only if
\[
\zeta_{\rm req}\le\frac{4}{4-x}.
\]

At fixed \(\zeta_{\rm req}\), the exact eigenvalue and Robin thresholds are
\[
y_{\rm req}^2=\frac{\pi^2}{x}\left(\frac{1}{\zeta_{\rm req}}-1+\frac{x}{4}\right),
\qquad
\eta_{\rm req}=y_{\rm req}\tan y_{\rm req}.
\]

## Stage 38 — explicit non-twin lowest-lane reachability

Combining the Stage-36 overlap family with Stage-37 Robin softening gives
\[
\zeta_0^{(\exp+R)}(\alpha,\eta)
=
\frac{\Omega_{\exp}(\alpha)^2}{1-x/4+xy(\eta)^2/\pi^2}.
\]

Its exact closure range is
\[
1\le \zeta_0^{(\exp+R)} \le \frac{\pi^2}{4-x}.
\]

So the explicit family reaches the Stage-35 threshold iff
\[
\zeta_{\rm req}\le\frac{\pi^2}{4-x}.
\]

This gives the exact three-regime split:
\[
\zeta_{\rm req}\le\frac{\pi^2}{4}
\quad\Rightarrow\quad
\text{overlap alone enough},
\]
\[
\frac{\pi^2}{4}<\zeta_{\rm req}\le\frac{\pi^2}{4-x}
\quad\Rightarrow\quad
\text{overlap + softening enough},
\]
\[
\zeta_{\rm req}>\frac{\pi^2}{4-x}
\quad\Rightarrow\quad
\text{even this explicit family fails}.
\]

## Stage 39 — transport origin of source asymmetry

The Stage-36 exponential family is exactly the stationary zero-flux branch of
\[
\partial_t \sigma + \partial_s J = 0,
\qquad
J=-D_\sigma \partial_s \sigma + v_\sigma \sigma.
\]

On the stationary recirculating branch \(J=0\),
\[
\sigma_{Pe}(s)
=
\frac{Pe\,e^{Pe s/L}}{e^{Pe}-1},
\qquad
Pe:=\frac{v_\sigma L}{D_\sigma}.
\]

The corresponding overlap boost is
\[
\Omega_{Pe}
=
\frac{\pi Pe\bigl(2Pe\,e^{Pe}+\pi\bigr)}
{(4Pe^2+\pi^2)(e^{Pe}-1)}.
\]

Exact identities:
\[
\Omega_{Pe}(0)=1,
\qquad
\lim_{Pe\to+\infty}\Omega_{Pe}=\frac{\pi}{2},
\]
and the score-function identity
\[
\partial_{Pe}p_{Pe}(x)=(x-\mathbb E_{Pe}[x])p_{Pe}(x)
\]
implies
\[
\frac{d}{dPe}\mathbb E_{Pe}[\chi_0]=\operatorname{Cov}_{Pe}(\chi_0,x),
\]
so the constructive branch is monotone increasing because \(\chi_0\) is increasing on \([0,1]\).

## Stage 40 — physical \((Pe,\kappa,\eta)\) placement map

Define the physical support ratios
\[
\kappa:=\frac{K_XL^2}{T_X},
\qquad
\eta:=hL=\frac{K_mL}{T_X}.
\]
Then
\[
x=\frac{\pi^2}{\kappa+\pi^2/4},
\]
and the Robin softening factor becomes
\[
A_K(\eta;\kappa)=\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2}.
\]

The explicit physical lowest-lane family is therefore
\[
\zeta_0^{(Pe+R)}(Pe,\eta;\kappa)
=
\Omega_{Pe}^2\frac{\kappa+\pi^2/4}{\kappa+y(\eta)^2}.
\]

This map is monotone:
- increasing in \(Pe\),
- decreasing in \(\eta\),
- decreasing in \(\kappa\).

Its exact constructive-branch ceiling is
\[
\zeta_{\max}(\kappa)=\frac{\pi^2}{4}\frac{\kappa+\pi^2/4}{\kappa}.
\]

So the Stage-35 demand is reachable on this first physical family iff
\[
\zeta_{\rm req}\le \zeta_{\max}(\kappa),
\]
equivalently, whenever \(\zeta_{\rm req}>\pi^2/4\),
\[
\kappa \le \kappa_{\max}(\zeta_{\rm req})
:= \frac{\pi^4}{4(4\zeta_{\rm req}-\pi^2)}.
\]

The exact physical threshold surfaces are
\[
\Omega_{\rm req}^2
=
\zeta_{\rm req}\frac{\kappa+y(\eta)^2}{\kappa+\pi^2/4},
\]
\[
y_{\rm req}^2
=
\frac{\Omega_{Pe}^2}{\zeta_{\rm req}}(\kappa+\pi^2/4)-\kappa,
\]
\[
\kappa_{\rm req}
=
\frac{\Omega_{Pe}^2\pi^2/4-\zeta_{\rm req}y(\eta)^2}{\zeta_{\rm req}-\Omega_{Pe}^2}.
\]

## Where the program stands after Stage 40

The support/normalization problem is no longer phrased in abstract deformation variables.
It has collapsed to three physical moving-throat operator ratios:

- axial source Peclet number `Pe`,
- mouth compliance `eta`,
- baseline support stiffness ratio `kappa`.

That makes the next clean move very sharp: derive the coupled support/source branch equation that selects `Pe`
from the same operator that already fixes `eta` and `kappa`.
