# 5PN / Moving-Throat continuation — Stages 333–335

This session pushed the coherent D/N support theorem one step further down to an explicit regime classifier and then tied it to the natural minimal isotropic conservative quadrupole precursor.

The immediate continuation point before this session was:

> determine whether the actual passive/outgoing branch lands in the mixed-only, lowest-twin, or genuinely non-twin support regime.

The result is that this classification is now exact in one scalar ratio, and the canonical minimal isotropic branch lands comfortably in the symmetric-lowest-twin window.

---

## Stage 333 — exact branch-product support regimes

Let
\[
C_{\rm mix}:=\frac{8\Lambda(1-\epsilon)}{\pi^2},
\qquad
S_{\rm req}:=\frac{\Pi_{\rm tr}}{C_{\rm mix}},
\]
and write the exact coherent support demand as
\[
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)}
=
\frac{\Pi_{\rm tr}-C_{\rm mix}}
{C_{\rm mix}-\epsilon(2C_{\rm mix}-\Pi_{\rm tr})}.
\]

The exact thresholds are:

- at \(\Pi_{\rm tr}=C_{\rm mix}\),
  \[
  \zeta_{\rm req}=0;
  \]
- at \(\Pi_{\rm tr}=2C_{\rm mix}\),
  \[
  \zeta_{\rm req}=1.
  \]

The exact derivative is
\[
\frac{d\zeta_{\rm req}}{d\Pi_{\rm tr}}>0
\]
throughout the blocked branch \(0<\epsilon<1\), so the support demand increases strictly with the tracking-branch product.

Therefore the support problem splits into three exact regimes:

\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\Longrightarrow\quad
\text{mixed-only already enough},
\]

\[
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\Longrightarrow\quad
\text{symmetric lowest twin enough},
\]

\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\Longrightarrow\quad
\text{non-twin asymmetry required}.
\]

So the next branch decision no longer depends on many reduced parameters separately. It depends only on where the actual branch lands relative to the two product thresholds \(C_{\rm mix}\) and \(2C_{\rm mix}\).

---

## Stage 334 — exact loading-ratio compiler and canonical isotropic conservative precursor

On the explicit support/source branch, the natural contact-plus-pole conservative precursor may be written as
\[
Y_Q^{\rm cons}(\omega)
=
c_0+\frac{c_1}{1-\omega^2/\Omega_Q^2},
\qquad
c_0+c_1=1.
\]

Its exact support/source loading-ratio dictionary is
\[
c_0=\frac{1}{\rho_\alpha},
\qquad
c_1=\frac{\rho_\alpha-1}{\rho_\alpha},
\qquad
\rho_\alpha=\frac{1}{c_0}=\frac{1}{1-c_1},
\]
with
\[
\zeta_{\rm req}=\frac{c_1}{c_0}.
\]

The minimal isotropic conservative quadrupole module already fixed earlier is
\[
c_0=\frac34,
\qquad
c_1=\frac14.
\]
So it implies exactly
\[
\rho_\alpha^{(\min)}=\frac43,
\qquad
\zeta_{\rm req}^{(\min)}=\frac13.
\]

Since
\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\rho_\alpha,
\]
the same branch gives
\[
\Pi_{\rm tr}^{(\min)}=\frac43\,C_{\rm mix}.
\]

With blocking retained,
\[
\zeta_{\rm req}^{(\min)}(\epsilon)
=
\frac{\rho_\alpha^{(\min)}-1}{1-\epsilon(2-\rho_\alpha^{(\min)})}
=
\frac{1}{3-2\epsilon}.
\]

This is still strictly below \(1\) for every blocked branch with \(0\le \epsilon<1\). So the canonical isotropic passive/outgoing branch is always twin-safe before any non-twin asymmetry is needed.

---

## Stage 335 — explicit Family-1 verdict for the canonical isotropic branch

The explicit Family-1 reduced theorem window already frozen earlier is

\[
\rho_\alpha \le 3.46622291347846
\quad\Longrightarrow\quad
\text{guaranteed success},
\]
\[
\rho_\alpha \ge 3.46752913273870
\quad\Longrightarrow\quad
\text{guaranteed failure},
\]
with hard constructive ceiling
\[
\rho_\alpha < 3.46752922945601.
\]

Comparing the canonical isotropic value
\[
\rho_\alpha^{(\min)}=\frac43
\]
to that window gives very large exact margins:
\[
\Delta_{\rm suff}=3.46622291347846-\frac43\approx 2.13288958014513,
\]
\[
\Delta_{\rm fail}=3.46752913273870-\frac43\approx 2.13419579940537,
\]
\[
\Delta_{\rm max}=3.46752922945601-\frac43\approx 2.13419589612268.
\]

Likewise, on the support-ratio side
\[
\zeta_{\rm req}^{(\min)}=\frac13,
\qquad
\zeta_{\max}^{(F1)}\approx 2.46752922945601,
\]
so
\[
\zeta_{\max}^{(F1)}-\zeta_{\rm req}^{(\min)}
\approx 2.13419589612268.
\]

This branch also satisfies the explicit zero-bias Family-1 criterion, because the exact trigger is
\[
\zeta_{\rm req}<A_{F1},
\qquad
A_{F1}\approx 1.00005192880220,
\]
and
\[
\frac13 < 1.00005192880220.
\]
Therefore
\[
Pe_{\rm req}=0
\]
on the explicit Family-1 branch for the canonical isotropic passive/outgoing quadrupole module.

So the reduced support/source side is no longer the live bottleneck on that branch.

---

## Net result after Stages 333–335

The moving-throat support/source theorem has now been reduced to one exact scalar classifier:

\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\rho_\alpha.
\]

That one number tells us everything:

- \(\rho_\alpha\le 1\): mixed-only enough,
- \(1<\rho_\alpha\le 2\): symmetric lowest twin enough,
- \(\rho_\alpha>2\): non-twin asymmetry required.

The canonical minimal isotropic conservative quadrupole precursor gives
\[
\rho_\alpha=\frac43,
\]
so it lies strictly in the symmetric-lowest-twin regime and, on the explicit Family-1 branch, already succeeds at zero transport bias.

So the reduced theorem gap is now even sharper than before:

> the explicit support/source side is no longer the active uncertainty on the canonical isotropic branch.  
> The remaining reduced question is exactly what loading ratio \(\rho_\alpha\) the actual passive/outgoing moving-throat quadrupole branch selects.
