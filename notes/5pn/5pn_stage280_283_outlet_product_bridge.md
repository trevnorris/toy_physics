# 5PN continuation — Stages 280–283: selected loading/product variables into compensated outlet observables

This session keeps the bulk-PDE firewall intact. Nothing here retunes the parent `4+1` GNLS/Maxwell medium. The work stays entirely on the selected branch / outlet side.

## Stage 280 — compensated hybrid outlet reproduces the minimal isotropic conservative module exactly

Start from the compensated Robin–mixed isotropic outlet branch

\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+4\sigma_W
-\frac{\sigma_W}{1-z^2/3-i\gamma_W z^5}
+O(z^6),
\]

with the compensated canonical-even conditions

\[
\rho_R = 4\sigma_W,
\qquad
\kappa_W = \frac13.
\]

Normalizing by the static outlet value gives

\[
\widehat Y_2^{\rm hyb}(z)
=
1+\frac{z^2}{9}+\frac{4z^4}{81}
+i\frac{1-9\sigma_W\gamma_W}{27(1-\sigma_W)}z^5+O(z^6).
\]

So the whole compensated hybrid outlet is exactly equivalent, through `O(z^5)`, to the minimal isotropic contact-plus-pole module

\[
\widehat Y_2^{\rm hyb}(z)
=
\frac34
+
\frac14\frac{1}{1-\frac49 z^2-i\frac{4}{27}\chi_Q z^5}
+O(z^6),
\qquad
\chi_Q=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]

Therefore the conservative selected-branch weights are frozen to

\[
c_0=\frac34,
\qquad
c_1=\frac14,
\qquad
\rho_\alpha=\frac{1}{c_0}=\frac43.
\]

Using the exact blocked-demand map

\[
\zeta_{\rm req}
=
\frac{\rho_\alpha-1}{1-\epsilon_{\rm blk}(2-\rho_\alpha)},
\]

this gives

\[
\zeta_{\rm req}(\epsilon_{\rm blk})=\frac{1}{3-2\epsilon_{\rm blk}}.
\]

But the selected loading/product ratio itself stays fixed:

\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
Q(\zeta_{\rm req},\epsilon_{\rm blk})
=
\frac43.
\]

So `epsilon_blk` changes the support-demand variable `zeta_req`, but on the minimal isotropic conservative branch it does **not** change the selected product ratio `Pi_tr/C_mix`.

## Stage 281 — exact outlet-observable transport on the compensated canonical-even branch

The compensated hybrid outgoing factor is

\[
\chi_Q(\sigma_W,\gamma_W)=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]

Its full first differential is

\[
d\chi_Q
=
\frac{1-9\gamma_W}{(1-\sigma_W)^2}\,d\sigma_W
-\frac{9\sigma_W}{1-\sigma_W}\,d\gamma_W.
\]

At the canonical outgoing point

\[
\gamma_W=\frac19,
\]

the `d sigma_W` term vanishes identically, so

\[
d\chi_Q\Big|_{\gamma_W=1/9}
=
-\frac{9\sigma_*}{1-\sigma_*}\,d\gamma_W.
\]

This matches the earlier outlet-projection formulas. Writing the first-order compensated-branch defects as

\[
\delta E_2 = \frac{\delta\mathcal C - 9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\]
\[
\delta E_4 = \frac{5\delta\mathcal C - 72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]
\[
\Delta_Q = \frac{\delta\mathcal C - 27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)},
\]

the canonical-even gate `delta E_2 = delta E_4 = 0` forces

\[
\delta\mathcal C = 0,
\qquad
\delta\kappa_W = 0,
\]

so the only surviving isotropic outlet defect is

\[
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}
\,\delta\gamma_W.
\]

Thus `sigma_W` survives only as a static amplification factor on the canonical-even tangent; the actual odd mismatch is carried purely by `delta gamma_W`.

## Stage 282 — exact loading/product to outlet-slippage bridge

The Family-1 loading side and the outlet side meet through the exact first-order identity

\[
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]

Comparing this with the canonical-even outlet formula

\[
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}
\,\delta\gamma_W,
\]

gives the exact tangent transport law

\[
\delta\gamma_W
=
\frac{1}{9}\,
\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]

On the weak-axisymmetric grouped branch, the same defect is also

\[
\Delta_Q = \frac{P_1}{P_0},
\]

because the compensated-branch slope is

\[
\gamma_1
=
-\frac{1-\sigma_*}{9\sigma_*}
\frac{P_1}{P_0}.
\]

So the remaining outgoing theorem gap is one scalar written in three exactly equivalent ways:

\[
\Delta_Q,
\qquad
\frac{P_1}{P_0},
\qquad
\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]

## Stage 283 — one-language end-to-end theorem ledger

After the bridge is made explicit, the selected-branch theorem data split cleanly into two packets.

### Conservative selected loading/product side

\[
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}(\epsilon_{\rm blk}) = \frac{1}{3-2\epsilon_{\rm blk}},
\qquad
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \frac43.
\]

### Compensated outlet side

\[
\chi_Q(\sigma_W,\gamma_W)=\frac{1-9\sigma_W\gamma_W}{1-\sigma_W},
\]

and on the nontrivial compensated branch the canonical odd condition is simply

\[
\chi_Q=1
\iff
\gamma_W=\frac19.
\]

### Tangent bridge between them

\[
\delta\gamma_W
=
\frac{1}{9}\,
\Xi_{\rm slip}\,\delta\Pi_{\rm tan},
\qquad
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}
\Xi_{\rm slip}\,\delta\Pi_{\rm tan},
\qquad
\Delta_Q = \frac{P_1}{P_0}.
\]

So the final outgoing theorem gap is now written in one language end to end:

1. selected conservative product side: `Pi_tr/C_mix = 4/3`,
2. tangential branch slippage: `delta Pi_tan`,
3. outlet odd slippage: `delta gamma_W = Xi_slip delta Pi_tan / 9`,
4. final outgoing defect: `Delta_Q = chi_Q - 1 = P1/P0`.

That is the sharpest compact bridge we have so far between the selected branch, the explicit outlet observables, and the remaining 2.5PN / 4PN outgoing normalization defect.
