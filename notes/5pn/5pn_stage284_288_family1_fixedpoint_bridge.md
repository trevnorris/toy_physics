# 5PN continuation — Stages 284–288: co-evolving Family-1 fixed-point bridge, on-family rigidity, and the off-family scalar

This session takes the end-to-end bridge from Stages 280–283

a) from the selected conservative product side,

b) through the compensated hybrid outlet,

c) down into the explicit co-evolving Family-1 fixed-point kernel.

The main result is that the final outgoing obstruction is now much sharper than
just “compute `deltaPi_tan`.”

On the explicit co-evolving Family-1 branch,

- `deltaPi_tan` is already explicit,
- the bare mixed-port odd slippage reduces to one D/N similarity-slippage scalar,
- that scalar vanishes identically on the exact lower parent compensation family,
- and any first-order failure is therefore carried by one off-family normal coordinate `delta_perp`.

So the next genuine PDE-facing question is no longer to guess an outlet defect. It is to determine whether the actual branch stays on the exact lower compensation family, and if not, what its single off-family scalar `delta_perp` is.

## Stage 284 — exact co-evolving Family-1 fixed-point transport

At the renormalized canonical point, the mouth/core load transport is

\[
\delta M_s = \delta\Sigma_0,
\]
\[
\delta M_q
=
-\frac14\,\delta\Sigma_0
+
\frac{\Sigma_0^{\rm can}}{\sqrt{1+\mathfrak r_{F1}^2}}\,\delta\mathfrak g,
\]
\[
\delta\Pi
=
\left(1-\frac14\mathcal S_{\rm can}\right)\delta\Sigma_0
-
\frac{\Sigma_0^{\rm can}}{4}\,\delta\mathcal S
+
\frac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_{F1}^2}}\,\delta\mathfrak g.
\]

On the canonical-even tangent `delta g = 0`, this collapses to

\[
\boxed{
\delta\Pi_{\rm tan}
=
\left(1-\frac14\mathcal S_{\rm can}\right)\delta\Sigma_0
-
\frac{\Sigma_0^{\rm can}}{4}\,\delta\mathcal S.
}
\]

Numerically,

\[
\boxed{
\delta\Pi_{\rm tan}
\approx
0.832409471081634\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S.
}
\]

Using

\[
\Sigma_0 = \frac{20}{9}\widehat T_m^2,
\qquad
\delta\Sigma_0 = \frac{40}{9}\widehat T_{m,\rm can}\,\delta\widehat T_m,
\]

this is equivalently

\[
\boxed{
\delta\Pi_{\rm tan}
\approx
5.35223887169622\,\delta\widehat T_m
-
1.16275838754222\,\delta\mathcal S.
}
\]

So the tangential defect transport is already explicit on the fixed-point branch.

## Stage 285 — bare mixed-port slippage and D/N similarity reduction

The compensated concrete core obeys

\[
\kappa_W = \frac{\kappa_0}{1+r_c},
\qquad
\gamma_W = \frac{\gamma_0}{1+r_c}.
\]

Linearizing on the compensated branch gives the exact identity

\[
\boxed{
\delta\gamma_W - \frac13\delta\kappa_W
=
\frac{1}{1+r_{c,*}}
\left(
\delta\gamma_0 - \frac13\delta\kappa_0
\right).
}
\]

Under the canonical-even gate `delta kappa_W = 0`, the odd outlet defect is therefore carried by the bare mixed-port slippage

\[
\delta\mathfrak B_W := \delta\gamma_0 - \frac13\delta\kappa_0,
\qquad
\boxed{
\delta\gamma_W = \frac{\delta\mathfrak B_W}{1+r_{c,*}}.
}
\]

Now use the D/N realization

\[
\kappa_0 = \frac{4L_W^2}{\pi^2 a^2},
\qquad
\gamma_0 = \frac{1+r_c}{9}.
\]

Then

\[
\delta\mathfrak B_W
=
\frac{1+r_{c,*}}{9}
\left[
\delta\ln\gamma_0 - 2\,\delta\ln\left(\frac{L_W}{a}\right)
\right].
\]

Define the D/N similarity-slippage scalar

\[
\Xi_{\rm slip} := \Xi_\gamma - 2\Xi_L.
\]

If `delta mathfrak B_W = Upsilon_Pi deltaPi_tan`, then

\[
\Upsilon_\Pi = \frac{1+r_{c,*}}{9}\,\Xi_{\rm slip},
\]

and the first-order outgoing defect collapses to

\[
\boxed{
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
}
\]

So the last effective odd defect is not a generic DtN susceptibility. It is one D/N similarity-strain mismatch.

## Stage 286 — parent compensation-surface rigidity and automatic similarity preservation

On the exact parent compensation family,

\[
1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2,
\qquad
\frac{L_W}{a} = \frac{\pi}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
\qquad
\gamma_0 = \frac{1+\mathfrak r^2}{9}.
\]

Therefore

\[
\delta\ln\gamma_0
-
2\,\delta\ln\left(\frac{L_W}{a}\right)=0
\qquad\Longrightarrow\qquad
\boxed{\Xi_{\rm slip}=0}
\]

identically along the family.

On the lower branch,

\[
\mathfrak g_-(\mathfrak r)=\mathfrak r-\frac12\sqrt{1+\mathfrak r^2},
\qquad
\mathfrak g_-'(\mathfrak r)
=
1-\frac{\mathfrak r}{2\sqrt{1+\mathfrak r^2}}>0.
\]

So the carried canonical-even condition `delta g = 0` forces

\[
\boxed{\delta\mathfrak r=0.}
\]

Hence all first-order D/N similarity defects freeze,

\[
\delta\ln\gamma_0=0,
\qquad
\delta\ln(L_W/a)=0,
\qquad
\delta\mathfrak B_W=0,
\qquad
\delta\gamma_W=0,
\]

and therefore

\[
\boxed{\Delta_Q=0,
\qquad
N_Q-1=0}
\]

at first order.

This is the strongest positive result of the batch:

> if the actual co-evolving moving-throat branch stays on the exact lower parent compensation family, the first-order reduced 2.5PN / 4PN outgoing defect disappears automatically.

## Stage 287 — off-family normal coordinate

To measure genuine departure from that family, define the exact parent compensation defect

\[
\mathcal F(\mathfrak g,\mathfrak r) := 1+\mathfrak r^2 - 4(\mathfrak g-\mathfrak r)^2,
\]

and the off-family normal coordinate

\[
\boxed{
\delta_\perp
:=
\delta\mathfrak g
-
\mathfrak g_-'(\mathfrak r_*)\,\delta\mathfrak r.
}
\]

Then on the lower branch,

\[
\boxed{
\delta\mathcal F = 4\sqrt{1+\mathfrak r_*^2}\,\delta_\perp,
}
\qquad
\boxed{
\delta R_q = -\frac{\delta_\perp}{\sqrt{1+\mathfrak r_*^2}}.
}
\]

The exact microscopic parent-variable formula is

\[
\boxed{
\delta_\perp
=
\mathfrak g_*\,\delta\ln\left(\frac{g_q K_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}
\,\delta\ln\left(\frac{K_s K_q}{\lambda^2}\right).
}
\]

So the whole first-order off-family defect is carried by one scalar only.

## Stage 288 — explicit microscopic log channels and lower-branch cancellation

Using the carried explicit throat-core formulas, the two imbalance channels become

\[
\delta\ln\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
\]

\[
\delta\ln\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
\]

So `delta_perp` is already an explicit linear combination of

\[
(\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W).
\]

Then inserting the exact lower compensated branch drift laws,

\[
\delta\ln L_W = \delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\left(
\delta\ln\mathcal Z_q - \delta\ln\rho_w + 3\delta\ln c_{s,w} + 2\delta\ln c_s - 5\delta\ln a
\right),
\]
\[
\delta\ln\mathcal T_m
=
\frac12\left(
\delta\ln\mathcal Z_q - \delta\ln\rho_w + 3\delta\ln c_{s,w} - 2\delta\ln c_s - 3\delta\ln a
\right),
\]

both imbalance channels vanish identically, hence

\[
\boxed{\delta_\perp=0}
\]

on the exact lower compensated branch.

So the theorem picture is now extremely tight:

- **on-family:** first-order outgoing defect vanishes automatically,
- **off-family:** the whole first-order defect is carried by the single scalar `delta_perp`.

## Net result after Stages 284–288

The remaining outgoing theorem gap is no longer the tangential load itself.
It is now the much smaller question:

> does the actual moving-throat branch stay on the exact lower parent compensation family to first order, or does it develop a nonzero off-family scalar `delta_perp`?

If it stays on-family, the first-order reduced 2.5PN / 4PN outgoing obstruction is gone. If not, the defect ledger is already explicit and one-dimensional.
