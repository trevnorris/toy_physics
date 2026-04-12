# 5PN continuation notes — Stages 254–257 exact Family-1 demand robustness after the `\Lambda_{\rm EM}` refresh

This session continued directly from the exact `\Lambda_{\rm EM}` refresh of the explicit Family-1 support/source branch.
The next clean move was **not** to re-assume the canonical `3/4 + 1/4` isotropic split immediately. Instead, the analysis kept a generic isotropic one-pole conservative quadrupole carrier all the way through the explicit Family-1 demand ceiling and only specialized to the actual isotropic branch at the end.

That avoids the main oversimplification risk that had been left open at the end of Stage 253.

The upshot is stronger than before:

> the refreshed Family-1 support/source side is not merely compatible with the canonical `c_{\rm pole}=1/4` branch. It remains safe on a **large exact interval** of isotropic one-pole conservative branches, and the actual isotropic branch stays automatic throughout the full admissible blocked regime.

---

## Stage 254 — exact Family-1 admissible one-pole window

Take a generic isotropic one-pole conservative quadrupole carrier
\[
\widehat Y_Q^{\rm cons}(\omega)
=
 c_{\rm geom}
+
\frac{c_{\rm pole}}{1-\omega^2/\Omega_Q^2},
\qquad
c_{\rm geom}+c_{\rm pole}=1.
\]
Then the exact reduced demand variables are
\[
\rho_\alpha = \frac{1}{1-c_{\rm pole}},
\]
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]
and the selected-branch product variable collapses exactly to
\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\rho_\alpha.
\]
So the exact Family-1 demand ceiling translates to one sharp conservative one-pole condition:
\[
 c_{\rm pole}
<
 c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})
=
\frac{\zeta_{\max}^{(F1)}(1-\epsilon_{\rm blk})}{1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}}.
\]
With the exact refreshed value
\[
\zeta_{\max}^{(F1)}\approx 2.4675297457259358,
\]
the unblocked (`\epsilon_{\rm blk}=0`) hard ceiling is
\[
 c_{\rm pole,max}^{(F1)}(0)
\approx 0.7116102605226109,
\qquad
c_{\rm geom,min}^{(F1)}(0)
\approx 0.2883897394773891.
\]
The shell-weighted guaranteed-success threshold is only slightly lower:
\[
 c_{\rm pole,suff}^{(\chi)}(0)
\approx 0.7115015750293280.
\]
So the explicit Family-1 branch leaves a very large exact admissible window in one-pole conservative language.

---

## Stage 255 — exact geometry-lane contamination ceiling

Now map the Family-1 one-pole ceiling back into the geometry-lane contamination variables using the Stage-75 exact one-pole relation
\[
 c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2}.
\]
Then the exact Family-1 safe region is
\[
1+\epsilon_4
<
4 c_\star (1+\epsilon_2)^2,
\]
with `c_\star` taken either as the hard ceiling `c_{\rm pole,max}` or the guaranteed-success ceiling `c_{\rm pole,suff}`.

Equivalently,
\[
\epsilon_{4,\max}(\epsilon_2;c_\star)
=
4 c_\star (1+\epsilon_2)^2 - 1.
\]
On the refreshed unblocked Family-1 branch this gives
\[
\epsilon_4 < 1.8464410420904435
\qquad (\epsilon_2=0,\; \text{hard ceiling}),
\]
while on the principal physical slice `1+\epsilon_2>0`
\[
\epsilon_2 > -0.4072809255389385
\qquad (\epsilon_4=0,\; \text{hard ceiling}).
\]
So the explicit support/source side remains tolerant to order-one geometry-lane contamination around the canonical point.
This is much stronger than merely checking the canonical `1/4` split by itself.

---

## Stage 256 — actual isotropic branch is automatic throughout the full blocked regime

Only after the generic window is established do we specialize to the actual isotropic one-pole branch
\[
 c_{\rm pole}=\frac14,
\qquad
c_{\rm geom}=\frac34.
\]
Then
\[
\rho_\alpha^{(\rm act)}=\frac43,
\qquad
\frac{\Pi_{\rm tr}^{(\rm act)}}{C_{\rm mix}}=\frac43,
\]
and the exact blocked-regime demand is
\[
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk})
=
\frac{1}{3-2\epsilon_{\rm blk}}.
\]
The admissible Family-1 blocked regime is
\[
0\le \epsilon_{\rm blk} < \epsilon_{\rm blk}^{\rm crit}
=
\frac{1}{\zeta_{\max}^{(F1)}}
\approx 0.40526360491990934.
\]
At the worst admissible blocking value,
\[
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk}^{\rm crit})
=
\frac{\zeta_{\max}^{(F1)}}{3\zeta_{\max}^{(F1)}-2}
\approx 0.45673095573242554,
\]
so the exact hard margin is still
\[
\zeta_{\max}^{(F1)}-
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk}^{\rm crit})
\approx 2.0107987899935102.
\]
Because `\zeta_{\max}^{(F1)} > 1`, the actual isotropic branch also stays in the symmetric-lowest-twin-safe regime throughout the full admissible blocked interval.
So the refreshed Family-1 support/source side is automatic not only at `\epsilon_{\rm blk}=0`, but on the whole blocked branch.

---

## Stage 257 — universal twin-safe strip and exact Family-1 extension

There is also a universal exact theorem independent of the explicit Family-1 details.
From
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]
solving `\zeta_{\rm req}=1` gives the exact, blocking-independent boundary
\[
 c_{\rm pole}=\frac12.
\]
So every isotropic one-pole conservative carrier with
\[
 c_{\rm pole}\le \frac12
\qquad\Longleftrightarrow\qquad
\rho_\alpha\le 2
\]
already lies in the universal symmetric-lowest-twin-safe strip.

The refreshed Family-1 branch extends this universal strip to the exact larger window
\[
 c_{\rm pole} < 0.7116102605226109,
\qquad
\rho_\alpha < 3.4675297457259358.
\]
So the extension beyond the universal twin-safe strip is
\[
\Delta c_{\rm pole}^{(F1)}
\approx 0.2116102605226109,
\qquad
\Delta \rho_\alpha^{(F1)}
\approx 1.4675297457259358.
\]
The actual isotropic branch sits at
\[
 c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\]
so it lies deeply inside both the universal strip and the larger exact Family-1 extension.

---

## Net result after Stage 257

The support/source side is now more robustly understood than before.

What is exact after this pass:

1. the exact `\Lambda_{\rm EM}`-refreshed Family-1 ceiling translated into generic isotropic one-pole variables,
2. the exact admissible window in `c_{\rm pole}` / `c_{\rm geom}` language,
3. the exact geometry-lane safe region in `(\epsilon_2,\epsilon_4)`,
4. the exact blocked-regime theorem for the actual isotropic `c_{\rm pole}=1/4` branch,
5. the exact universal twin-safe strip `c_{\rm pole}\le 1/2`,
6. and the exact extension supplied by the explicit Family-1 branch beyond that strip.

So the remaining theorem gap is **not** explicit support/source supply on the isotropic one-pole branch.
Even without re-assuming the canonical `1/4` split until the end, the refreshed Family-1 side stays safely non-bottlenecked over a large exact neighborhood of that branch.

The live work remains where the compact moving-throat ledger said it should be:

> branch realization and outgoing DtN selection/normalization on the actual moving-throat grouped-`P2` branch, not wall-depth or support-source supply on the explicit Family-1 side.
