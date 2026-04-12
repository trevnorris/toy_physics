# 5PN continuation notes — stages 237 through 241

These stages fold the Stage-236 twin-selection theorem into the explicit moving-throat support/source operator from the compact PDE program.

The main structural change is that the coherent support problem is no longer phrased in abstract overlap variables.
It is now organized in three nested layers:

1. the **twin theorem** from Stage 236,
2. the first explicit **physical non-twin lowest-lane family** in terms of transport/compliance/stiffness variables,
3. the **operator-selected** and then **parent-projected** microscopic gain thresholds.

So after Stage 241 the next theorem gate is not “derive more support formulas somehow.”
It is:

> compute the actual parent confinement-loading amplitude and the equilibrium source/support alignment on the real moving-throat branch, then compare the resulting microscopic gain directly against the exact fail/succeed thresholds.

## Stage 237 — physical lowest-lane placement map

Starting from the first explicit non-twin lowest support family from the moving-throat notes, the physical support ratio is now

\[
\zeta_0^{(Pe+R)}(Pe,y;\kappa)
=
\Omega_{Pe}^2\,\frac{\kappa+\pi^2/4}{\kappa+y^2},
\]

with

\[
\Omega_{Pe}
=
\frac{\pi Pe\bigl(2Pe\,e^{Pe}+\pi\bigr)}{(4Pe^2+\pi^2)(e^{Pe}-1)}.
\]

Here:

- `Pe` is the physical axial source-drift Peclet number,
- `y` is the Robin eigenvalue parameter for the lowest support lane,
- `kappa = K_X L^2/T_X` is the baseline support stiffness ratio.

Exact identities proved in the script:

\[
\Omega_{Pe}(0)=1,
\qquad
\lim_{Pe\to +\infty}\Omega_{Pe}=\frac{\pi}{2}.
\]

So the Stage-236 symmetric twin baseline is recovered exactly at

\[
Pe=0,
\qquad
y=\frac{\pi}{2},
\qquad
\zeta_0^{(Pe+R)}=1.
\]

The exact closure ceiling of the first explicit non-twin family is

\[
\zeta_{\max}(\kappa)
=
\frac{\pi^2}{4}\,\frac{\kappa+\pi^2/4}{\kappa}.
\]

Therefore the Stage-236 twin theorem expands to the first physical phase split:

\[
\zeta_{\rm req}\le 1
\quad\Rightarrow\quad
\text{symmetric lowest twin enough},
\]

\[
1<\zeta_{\rm req}\le \zeta_{\max}(\kappa)
\quad\Rightarrow\quad
\text{first explicit non-twin family can in principle rescue},
\]

\[
\zeta_{\rm req}>\zeta_{\max}(\kappa)
\quad\Rightarrow\quad
\text{first explicit non-twin family fails identically}.
\]

Equivalently, when `\zeta_req > \pi^2/4`, the exact stiffness ceiling is

\[
\kappa \le \kappa_{\max}(\zeta_{\rm req})
:=
\frac{\pi^4}{4(4\zeta_{\rm req}-\pi^2)}.
\]

## Stage 238 — operator-selected support branch and exact `Xi` thresholds

The Stage-237 placement map is then evaluated on the first explicit coupled support/source branch.
The transport bias is no longer free; it is the root of the exact fixed-point equation

\[
Pe_* = \Xi\,\Delta(Pe_*;\kappa,\eta),
\]

with support-drop kernel endpoints

\[
\Delta_0(\kappa,\eta)
=
\frac{\eta(\cosh\sqrt{\kappa}-1)}{\kappa\bigl(\eta\cosh\sqrt{\kappa}+\sqrt{\kappa}\sinh\sqrt{\kappa}\bigr)},
\]

\[
\Delta_{\infty}(\kappa,\eta)
=
\frac{\eta\sinh\sqrt{\kappa}+\sqrt{\kappa}(\cosh\sqrt{\kappa}-1)}{\sqrt{\kappa}\bigl(\eta\cosh\sqrt{\kappa}+\sqrt{\kappa}\sinh\sqrt{\kappa}\bigr)}.
\]

Every constructive branch point obeys the exact bracket

\[
\Xi\Delta_0\le Pe_*\le \Xi\Delta_{\infty}.
\]

This immediately gives the exact support thresholds

\[
\Xi_{\rm fail} = \frac{Pe_{\rm req}}{\Delta_{\infty}(\kappa,\eta)},
\qquad
\Xi_{\rm suff} = \frac{Pe_{\rm req}}{\Delta_0(\kappa,\eta)},
\]

where `Pe_req` is the unique constructive transport bias that hits the demanded support ratio.
So the operator-selected branch has the exact three-zone structure:

- `Xi <= Xi_fail` : impossible in this operator family,
- `Xi >= Xi_suff` : guaranteed success,
- `Xi_fail < Xi < Xi_suff` : only then is the full root solve still needed.

Useful exact endpoint data verified in the script:

\[
\lim_{\kappa\to 0^+}\Delta_0 = \frac12,
\qquad
\lim_{\kappa\to 0^+}\Delta_{\infty}=1,
\]

and in the highly compliant mouth limit,

\[
\lim_{\eta\to +\infty}\Delta_0 = \frac{1-\operatorname{sech}(\sqrt{\kappa})}{\kappa},
\qquad
\lim_{\eta\to +\infty}\Delta_{\infty} = \frac{\tanh(\sqrt{\kappa})}{\sqrt{\kappa}}.
\]

## Stage 239 — parent-action microscopic gain

The phenomenological operator strength is then pushed back to the parent 4D action.
Projecting the `n=5` GNLS/compressional sector onto one source channel and one support channel gives the exact microscopic gain

\[
G_{\rm micro}
=
\frac{\rho_* g_\phi^2 O_{\sigma\phi}^2}{m c_{s*}^2 K_X N_{\sigma\sigma}}.
\]

Equivalently, introducing the source/support coherence factor

\[
C_{\sigma\phi}^2 := \frac{O_{\sigma\phi}^2}{N_{\sigma\sigma}N_{\phi\phi}},
\]

one gets the exact factorization

\[
G_{\rm micro}
=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s*}^2 K_X}
\,C_{\sigma\phi}^2,
\qquad
0\le C_{\sigma\phi}^2\le 1.
\]

So the best-case gain at fixed confinement loading is

\[
G_{\max}(g_\phi)
=
\frac{\rho_* g_\phi^2 N_{\phi\phi}}{m c_{s*}^2 K_X},
\]

reached only at perfect source/support alignment.

Using `\kappa = K_X L^2/T_X`, the exact fixed-point strength becomes

\[
\Xi_{\rm micro}=\kappa G_{\rm micro}
=
\frac{\rho_* g_\phi^2 O_{\sigma\phi}^2 L^2}{m c_{s*}^2 T_X N_{\sigma\sigma}}.
\]

So the support/source theorem gap is no longer an abstract gain problem.
It is a parent-overlap problem.

## Stage 240 — exact parent threshold decision theorem

Combining Stages 237–239 gives the exact microscopic phase diagram

\[
G_{\rm fail}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_{\infty}(\kappa,\eta)},
\qquad
G_{\rm suff}(\kappa,\eta)=\frac{Pe_{\rm req}}{\kappa\Delta_0(\kappa,\eta)}.
\]

So inside the non-twin reachability window:

- `G_micro <= G_fail` : fail,
- `G_micro >= G_suff` : success,
- only the band `G_fail < G_micro < G_suff` requires the full fixed-point solve.

In parent variables, this becomes exact threshold surfaces on the confinement-loading amplitude:

\[
g_{\phi,\rm fail}^2
=
\frac{m c_{s*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2\Delta_{\infty}(\kappa,\eta)},
\]

\[
g_{\phi,\rm suff}^2
=
\frac{m c_{s*}^2 T_X N_{\sigma\sigma} Pe_{\rm req}}{\rho_* L^2 O_{\sigma\phi}^2\Delta_0(\kappa,\eta)},
\]

and on the coherence factor:

\[
C_{\rm fail}^2
=
\frac{m c_{s*}^2 K_X G_{\rm fail}}{\rho_* g_\phi^2 N_{\phi\phi}},
\qquad
C_{\rm suff}^2
=
\frac{m c_{s*}^2 K_X G_{\rm suff}}{\rho_* g_\phi^2 N_{\phi\phi}}.
\]

This also yields the first exact parent-overlap no-go theorem:
if

\[
G_{\max}(g_\phi) < G_{\rm fail}(\kappa,\eta),
\]

then the branch fails for **every possible** source profile in the chosen support channel.

So the full support decision now has a strict order:

1. check the Stage-236 twin theorem,
2. if needed, check the Stage-237 non-twin reachability ceiling,
3. only then compare the actual parent branch against the exact microscopic gain thresholds.

## Stage 241 — parent equilibrium alignment removes the free coherence datum

The parent equilibrium equations sharpen things one more step.
On the local static branch,

\[
H(y)\,\delta\rho(s,y)+\delta V_{\rm conf}(s,y)=0,
\qquad
H(y):=h'(\rho_*(y)),
\]

so the source profile induced by a support displacement is not free:

\[
\chi_\sigma(y)=\frac{g_\phi\chi_\phi(y)}{H(y)}.
\]

This forces the overlap invariants to become

\[
O_{\sigma\phi}=g_\phi I_1,
\qquad
N_{\sigma\sigma}=g_\phi^2 I_2,
\]

with

\[
I_1=\int d^3y\,\frac{\chi_\phi(y)^2}{H(y)},
\qquad
I_2=\int d^3y\,\frac{\chi_\phi(y)^2}{H(y)^2}.
\]

Therefore the coherence factor is now derived, not free:

\[
C_{\sigma\phi}^2=
\frac{I_1^2}{N_{\phi\phi} I_2}\le 1.
\]

The exact equilibrium gain becomes

\[
G_{\rm eq} = \frac{g_\phi^2 I_1}{K_X}.
\]

In the thin matched layer where `H(y)` is nearly constant,

\[
C_{\sigma\phi}^2=1,
\qquad
G_{\rm eq} = \frac{g_\phi^2 N_{\phi\phi}}{K_X H_w}
= \frac{\rho_w g_\phi^2 N_{\phi\phi}}{m c_{s,w}^2 K_X}.
\]

So the best-alignment formulas used above are not ad hoc. They are the natural thin-layer limit of the parent equilibrium branch.

## Net result after Stage 241

The support theorem is now much narrower than it was at Stage 236.

1. The symmetric lowest twin lane still succeeds exactly iff `zeta_req <= 1`.
2. If `zeta_req > 1`, the first explicit non-twin lowest-lane family is controlled by the physical placement map `zeta_0^(Pe+R)(Pe,y;kappa)`.
3. The transport bias `Pe` is not free; it is selected by the exact support/source fixed-point equation `Pe = Xi Delta(Pe;kappa,eta)`.
4. The microscopic gain is no longer phenomenological; it is a parent-action overlap quantity.
5. On the equilibrium branch, the source/support coherence is no longer independent either.

So the remaining theorem gap is now almost completely localized:

> compute the actual parent confinement-loading amplitude `g_phi`, the active-layer compressional stiffness profile `H(y)`, and the support profile `chi_phi(y)` on the real moving-throat branch, form the resulting equilibrium gain `G_eq`, and compare it directly to the exact fail/succeed thresholds inside the first explicit non-twin support family.
