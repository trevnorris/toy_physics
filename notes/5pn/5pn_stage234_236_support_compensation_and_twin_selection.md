# 5PN continuation notes — stages 234 through 236

These stages push the coherent tracking branch into the first explicit support-compensation and harmonic-selection theorem.

## Stage 234 — exact support-compensation theorem on the coherent tracking branch

Using the tracking-branch functions

\[
G_{\rm tr}(\xi,\delta;R)=\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]

\[
F_{\rm tr}(\xi,\delta;R)=
\frac{\bigl[9\delta+(9+2R^2)\xi\bigr]^2\bigl[9\delta+(9+2R)\xi\bigr]^2}
{81(1-\xi)\bigl(9\delta^2+18\delta\xi+(9+2R^2)\xi^2\bigr)^2},
\]

and the coherent support-enhancement factor

\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]

the exact support threshold is

\[
\zeta_{\rm req}=\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)},
\qquad
S_{\rm req}=\frac{M_{\rm req}}{M_{\rm mix}},
\qquad
M_{\rm req}=G_{\rm tr}(\xi_{\rm req},\delta;R).
\]

The derivative

\[
\frac{dS}{d\zeta}=\frac{1-\epsilon}{(1-\zeta\epsilon)^2}>0
\]

shows strict monotonicity, so the support threshold is unique once the stable-side branch point \(\xi_{\rm req}\) is fixed.

Because

\[
\frac1\epsilon-\zeta_{\rm req}=\frac{1-\epsilon}{\epsilon\,[1+\epsilon(S_{\rm req}-2)]}>0,
\]

there is no reduced-level support no-go for any finite \(S_{\rm req}>1\) on \(0<\epsilon<1\).

## Stage 235 — explicit D/N overlap extraction of the coherent support ratio

For the exact finite-throat D/N family

\[
\chi_n(s)=\sqrt{\frac{2}{L}}\sin\!\frac{(n+1/2)\pi s}{L},
\]

with the first uniform local source density \(\sigma(s)=1\), the exact overlap law is

\[
I_n=\int_0^L \chi_n(s)\,ds
=\frac{2\sqrt{2L}}{\pi(2n+1)},
\qquad
\frac{I_n}{I_0}=\frac{1}{2n+1}.
\]

So the physical coherent support ratio becomes

\[
\zeta_n^{(\rm phys)}=
\frac{K_W^{(\rm eff)}}{K_{\phi,n}^{(\rm eff)}}\,\frac{1}{(2n+1)^2}.
\]

On the same-operator twin family,

\[
K_{\phi,n}^{(\rm eff)}=K_W^{(\rm eff)}\bigl(1+x n(n+1)\bigr),
\]

this collapses to

\[
\zeta_n^{(\rm twin)}=
\frac{1}{(2n+1)^2\bigl(1+x n(n+1)\bigr)}.
\]

The exact lowest-twin value is therefore

\[
\zeta_0^{(\rm twin)}=1.
\]

## Stage 236 — exact lowest-twin sufficiency and higher-harmonic selection rules

Because \(\zeta_0^{(\rm twin)}=1\), the lowest symmetric twin lane has exact enhancement

\[
S_0=S(1;\epsilon)=2,
\]

independent of \(\epsilon\).

Therefore the lowest symmetric twin lane succeeds exactly when

\[
\zeta_{\rm req}\le 1,
\qquad\text{equivalently}\qquad
S_{\rm req}\le 2.
\]

For higher D/N harmonics,

\[
\zeta_n^{(\rm twin)}\le \frac{1}{(2n+1)^2},
\]

so harmonic \(n\) is ruled out immediately if

\[
\zeta_{\rm req}>\frac{1}{(2n+1)^2}.
\]

When that immediate impossibility bound is not violated, the exact twin softness threshold is

\[
x\le x_{\max}(n;\zeta_{\rm req})
:=
\frac{1/\bigl((2n+1)^2\zeta_{\rm req}\bigr)-1}{n(n+1)}.
\]

The exact enhancement ceiling at fixed \(\epsilon\) is

\[
S_n^{(\max)}
=
1+\frac{1-\epsilon}{(2n+1)^2-\epsilon}.
\]

So the explicit D/N support tower is strongly biased toward the lowest symmetric twin lane.

## Net result after stage 236

The coherent support question is no longer vague.

1. The tracking-branch support threshold is one exact scalar \(\zeta_{\rm req}\).
2. The explicit D/N overlap tower makes that threshold microscopic.
3. The lowest symmetric twin lane is an exact universal doubling branch with \(S_0=2\).
4. Every higher support harmonic is rapidly ruled out by exact overlap and stiffness suppression.

So the next clean theorem gate is now extremely narrow:

> Does the actual moving-throat PDE place the physical coherent support sector on the lowest twin D/N branch with \(\zeta_{\rm req}\le 1\), or is stronger-than-twin asymmetry required?
