# 5PN continuation notes — Stages 250–253 exact `\Lambda_{\rm EM}` refresh of the explicit Family-1 branch

This session continued the explicit Family-1 support/source branch, but with one correction carried all the way through:

> do **not** use the shorthand reference freeze `L/a = 37/20`; use the exact derived EM-branch aspect ratio instead.

The carried exact 2PN geometry relation is
\[
\Lambda_{\rm EM}
=
\frac{\sqrt{2}\,\pi}{x_{01}},
\qquad
x_{01}=2.40482555769577\ldots
\]
so numerically
\[
\Lambda_{\rm EM}\approx 1.847486577120128.
\]

The older use of
\[
\frac{L}{a}=\frac{37}{20}=1.85
\]
was only a convenient **reference-branch numerical freeze**. This session replaced that freeze everywhere in the explicit Family-1 branch formulas.

---

## Stage 250 — exact `\Lambda_{\rm EM}` geometry refresh

With the explicit thin-wall Family-1 lock
\[
\frac{\ell}{a}=\frac1{20},
\]
the exact carried aspect ratio gives
\[
\Lambda_\ell:=\frac{L}{\ell}
=
20\,\Lambda_{\rm EM}
=
\frac{20\sqrt2\,\pi}{x_{01}}
\approx 36.94973154240256,
\]
so
\[
\eta=\Lambda_\ell\approx 36.94973154240256.
\]

With the healing-width lock
\[
\chi_s=\frac{L}{2\ell},
\]
we get
\[
\chi_s
=
10\,\Lambda_{\rm EM}
=
\frac{10\sqrt2\,\pi}{x_{01}}
\approx 18.47486577120128.
\]

Then the explicit support parameter becomes
\[
\kappa
=
4\chi_s^2+\frac45\Lambda_\ell^2
=
\frac95\Lambda_\ell^2
=
\frac{1440\pi^2}{x_{01}^2}
\approx 2457.508789900114.
\]

Relative to the shorthand reference freeze, this is only a small correction:
\[
\frac{\Lambda_{\rm EM}-(37/20)}{37/20}
\approx
-1.358606962\times 10^{-3}.
\]

So the branch geometry moves only slightly, but from this point onward the exact EM-branch formula is the one that should be propagated.

---

## Stage 251 — refreshed explicit threshold window

Using the exact refreshed values of `\eta` and `\kappa`, the explicit Family-1 support theorem gives
\[
\Delta_0 \approx 1.737739392346995\times 10^{-4},
\qquad
\Delta_\infty \approx 2.0172162594593645\times 10^{-2}.
\]

Therefore
\[
\Upsilon_{\rm fail}
\approx
3.630989267026856\times 10^{-2}\,Pe_{\rm req},
\]
\[
\Upsilon_{\rm suff}
\approx
4.214953415699773\,Pe_{\rm req},
\]
and hence
\[
\Theta_{\rm fail}
\approx
3.630989267026856\times 10^{-4}\,Pe_{\rm req},
\]
\[
\Theta_{\rm suff}
\approx
4.214953415699773\times 10^{-2}\,Pe_{\rm req}.
\]

The important observation is:

- the **fail** threshold shifts slightly upward relative to the `37/20` freeze,
- the **success** threshold is unchanged to the displayed precision.

So the explicit branch theorem is stable under the exact geometry refresh.

---

## Stage 252 — refreshed wall-depth verdict

The explicit wall-depth extraction from the Family-1 wall profile is unchanged:

\[
\Theta_w^{(\chi)} \approx 4.06863235008162\,\lambda_\mu^2,
\qquad
\Theta_w^{(J)} \approx 0.927552032539308\,\lambda_\mu^2.
\]

Comparing these to the refreshed threshold window gives

\[
Pe_{\rm suff}^{(\chi)}
=
\frac{\Theta_w^{(\chi)}}{\Theta_{\rm suff}}
\approx
96.5285247264385\,\lambda_\mu^2,
\]
\[
Pe_{\rm fail}^{(\chi)}
=
\frac{\Theta_w^{(\chi)}}{\Theta_{\rm fail}}
\approx
11205.2998532081\,\lambda_\mu^2.
\]

For the conservative lower envelope,
\[
Pe_{\rm suff}^{(J)}
\approx
22.0062226330754\,\lambda_\mu^2,
\qquad
Pe_{\rm fail}^{(J)}
\approx
2554.54358117343\,\lambda_\mu^2.
\]

So the qualitative branch-level verdict survives exactly as before:

> on the first explicit Family-1 branch, wall-depth is still **not** the natural bottleneck.

The wall-depth side still succeeds automatically for modest quadrupole demand and fails only for anomalously large demand.

---

## Stage 253 — refreshed quadrupole-demand window in `\zeta_{\rm req}` and `\Pi_{\rm tr}`

The exact Family-1 support-ratio map now uses the refreshed
\[
\eta=\Lambda_\ell\approx 36.94973154240256,
\qquad
\kappa\approx 2457.508789900114.
\]

Solving
\[
y\tan y=\eta
\]
on the principal branch gives
\[
y_{F1}\approx 1.5294278190457656,
\]
and therefore
\[
A_{F1}
=
\frac{\kappa+\pi^2/4}{\kappa+y_{F1}^2}
\approx
1.0000521380385143.
\]

The hard Family-1 support ceiling becomes
\[
\zeta_{\max}^{(F1)}
=
A_{F1}\frac{\pi^2}{4}
\approx
2.467529745725936.
\]

At `\lambda_\mu=1`, the refreshed explicit branch windows are

\[
\zeta_{\rm suff}^{(\chi)}(1)\approx 2.466223429475074,
\qquad
\zeta_{\rm fail}^{(\chi)}(1)\approx 2.467529648745268,
\]
\[
\zeta_{\rm suff}^{(J)}(1)\approx 2.442576225820804,
\qquad
\zeta_{\rm fail}^{(J)}(1)\approx 2.467527879753313.
\]

So the natural shell-weighted guaranteed-success threshold still sits only
\[
\zeta_{\max}^{(F1)}-\zeta_{\rm suff}^{(\chi)}(1)
\approx
1.306250899\times 10^{-3}
\]
below the hard Family-1 ceiling, while the guaranteed-failure threshold is essentially saturated at the ceiling.

In the unblocked product language `(\epsilon_{\rm blk}=0)`,
\[
\Pi_{\rm suff}^{(\chi)}/C_{\rm mix}\approx 3.466223429475074,
\]
\[
\Pi_{\rm fail}^{(\chi)}/C_{\rm mix}\approx 3.467529648745268,
\]
\[
\Pi_{\max}^{(F1)}/C_{\rm mix}\approx 3.467529745725936.
\]

The corresponding exact blocking ceiling is
\[
\epsilon_{\rm blk}^{\rm crit}
=
\frac{1}{\zeta_{\max}^{(F1)}}
\approx
0.405263604919909.
\]

So the explicit-branch conclusion is unchanged even in the exact `\Lambda_{\rm EM}` refresh:

> the Family-1 wall-depth side still pushes the guaranteed-success threshold extremely close to the hard support ceiling, so the unresolved issue remains the **selected quadrupole demand side**, not wall-depth supply.

---

## Net result after Stage 253

Replacing the shorthand
\[
L/a = 37/20
\]
by the exact carried relation
\[
L/a = \Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}
\]
does **not** change the qualitative support/source verdict.

It does three useful things:

1. it removes the last explicit use of a numerical freeze where an exact carried equation exists,
2. it refreshes the fail-side threshold numbers consistently,
3. and it confirms that the explicit Family-1 subprogram is still bottlenecked by the **quadrupole-demand branch**.

So the next honest theorem gate remains the same in substance, but now with the corrected geometry already folded in:

> compare the final selected quadrupole branch demand directly to the refreshed explicit Family-1 ceiling in `\zeta_{\rm req}` or `\Pi_{\rm tr}` language.
