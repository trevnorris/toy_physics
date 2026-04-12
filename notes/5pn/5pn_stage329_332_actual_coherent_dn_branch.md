# 5PN continuation notes — Stages 329–332

These stages take the Stage 323–328 orbit/support packet compiler and splice it directly into the first concrete coherent local D/N kernel family from the moving-throat support-compensation program.

The main point is that the abstract branch packet is no longer floating over an unspecified support closure. On the concrete coherent local D/N family, the physical support ratio is explicit, the orbit packet stays completely support-blind, and the support theorem collapses to one exact regime test on the same branch.

## Stage 329 — concrete coherent local D/N branch state

Start from the microscopic coherent kernel family in which the mixed lane `W` and the support lane `phi_n` are sourced by the same local throat density. Then

\[
\lambda_W = \lambda_* I_0,
\qquad
\lambda_{\phi,n} = \lambda_* I_n,
\]
with the exact coherent placement-state map
\[
\chi_0 = \frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U = \frac{\pi^2 T_U}{L^2 K_U},
\qquad
Z_W = \frac{\lambda_W^2}{K_{\eta}^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\]
\[
\epsilon_W = \frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}},
\qquad
\epsilon_\eta = \frac{c_{\eta U}^2}{K_U K_{\eta}^{(\mathrm{eff})}},
\qquad
\Lambda = \frac{27\pi^2 G c_s^5}{20 a^5 c^5}\frac{K_W^{(\mathrm{eff})}}{\mu_W}.
\]

On the same branch the physical coherent support ratio is no longer free. It is
\[
\zeta_n^{(\mathrm{phys})}
=
\frac{K_W^{(\mathrm{eff})}}{K_{\phi,n}^{(\mathrm{eff})}}
\left(\frac{I_n}{I_0}\right)^2.
\]
For the first uniform local source density,
\[
\frac{I_n}{I_0} = \frac{1}{2n+1},
\]
so
\[
\zeta_n^{(\mathrm{uniform})}
=
\frac{K_W^{(\mathrm{eff})}}{K_{\phi,n}^{(\mathrm{eff})}}\frac{1}{(2n+1)^2}.
\]
On the same-operator twin family,
\[
\zeta_n^{(\mathrm{twin})}=
\frac{1}{(2n+1)^2\bigl(1+x\,n(n+1)\bigr)}.
\]
The lowest symmetric twin lane is exact:
\[
\zeta_0^{(\mathrm{twin})}=1,
\qquad
S_0=2.
\]
So the lowest twin branch is the universal doubling branch of the concrete coherent support sector.

## Stage 330 — actual coherent local D/N orbit packet

On the actual coherent local D/N branch, the orbit packet is carried entirely by the three exact branch observables
\[
(R_{\rm tr},\,R_{\rm target},\,\epsilon_\eta),
\]
with
\[
R_{\rm tr}=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\qquad
R_{\rm target}=
\Lambda\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\epsilon=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right).
\]

The exact finite quotient packet is therefore
\[
q_{\rm tr} = -C_*\ln\frac{R_{\rm tr}}{R_{{\rm tr},{\rm ref}}},
\]
\[
q_{\rm nt} = B_*\ln\frac{R_{\rm tr}}{R_{{\rm tr},{\rm ref}}}
+ \ln\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}
- \ln\frac{R_{\rm target}}{R_{{\rm target},{\rm ref}}},
\]
\[
q_\eta = \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]
The infinitesimal defect packet is still
\[
\Theta_1=d\ln R_{\rm tr},
\qquad
\Xi_1=-d\ln R_{\rm target}-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}d\ln\epsilon_\eta,
\qquad
\mathcal R_1=d\ln R_{\rm target}.
\]

Most importantly, every one of these objects is exactly support-blind:
\[
\partial_\zeta q_{\rm tr}=
\partial_\zeta q_{\rm nt}=
\partial_\zeta q_\eta=
\partial_\zeta\Theta_1=
\partial_\zeta\Xi_1=
\partial_\zeta\mathcal R_1=0.
\]
So choosing the physical D/N support harmonic does not move the branch on or off the weak-axisymmetric similarity orbit.

## Stage 331 — actual coherent local D/N support threshold

On the same physical branch, the support theorem is controlled by the tracking-branch product
\[
\Pi_{\rm tr}(\xi,\delta;R)=F_{\rm tr}(\xi,\delta;R)\,G_{\rm tr}(\xi,\delta;R),
\]
and the mixed-only product scale
\[
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]
Then
\[
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}},
\qquad
\zeta_{\rm req}=
\frac{\Pi_{\rm tr}-C_{\rm mix}}{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]

The exact regime split is:

1. mixed-only enough iff
   \[
   \Pi_{\rm tr} \le C_{\rm mix};
   \]
2. the symmetric lowest twin lane is enough iff
   \[
   C_{\rm mix} < \Pi_{\rm tr} \le 2C_{\rm mix};
   \]
3. non-twin asymmetry is required iff
   \[
   \Pi_{\rm tr} > 2C_{\rm mix}.
   \]

So on the concrete coherent local D/N family the support problem is no longer vague. The first physical support lane is the universal doubling branch, and it succeeds or fails according to one exact branch-product inequality.

Higher D/N harmonics are strongly suppressed. For the same-operator twin family,
\[
\zeta_n^{(\mathrm{twin})}
=
\frac{1}{(2n+1)^2(1+x n(n+1))}
<
\frac{1}{(2n+1)^2}
\qquad (n\ge 1),
\]
so an exact necessary condition for the `n`th twin harmonic to work is
\[
\zeta_{\rm req} \le \frac{1}{(2n+1)^2}.
\]
If that is satisfied, the exact softness threshold is
\[
x \le x_{\max}(n;\zeta_{\rm req})
=
\frac{1/((2n+1)^2\zeta_{\rm req})-1}{n(n+1)}.
\]

## Stage 332 — actual coherent local D/N branch theorem gate

Putting the two packets together, the actual coherent local D/N branch now ends at one exact theorem gate:

### Orbit packet
\[
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
\]
or equivalently
\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln\epsilon_\eta=0.
\]

### Support packet
With
\[
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\]
we have
- mixed-only enough iff `Π_tr <= C_mix`,
- lowest symmetric twin enough iff `C_mix < Π_tr <= 2 C_mix`,
- non-twin asymmetry required iff `Π_tr > 2 C_mix`.

On the lowest twin branch `zeta=1`,
\[
M_{\rm tr}=2M_{\rm mix},
\qquad
R_{\rm target}M_{\rm tr}=\frac{16\Lambda(1-\epsilon)}{\pi^2}.
\]

## Where we end up

The abstract packet compiler is now fully tied to the first concrete coherent local D/N operator family.

1. The actual weak-axisymmetric orbit-lock test is carried exactly by
   \[
   (R_{\rm tr},R_{\rm target},\epsilon_\eta),
   \]
   and it is rigorously blind to the support lane.
2. The actual coherent support theorem is carried exactly by
   \[
   \Pi_{\rm tr},\quad C_{\rm mix},\quad \zeta_n^{(\mathrm{phys})},
   \]
   with the lowest symmetric twin branch giving the universal doubling law.
3. So the remaining PDE-side question is not another algebraic compression. It is simply:

> does the completed moving-throat operator place the physical branch in the exact orbit-lock locus, and if so, is its support sector already in the mixed-only regime, the lowest-twin regime, or the genuinely non-twin regime?
