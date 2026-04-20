# Moving-Throat PDE — Stage 150: Bundle Transport of the Remaining Mouth Variables and the Tangent-Compensation Theorem

## Goal

Stage 149 reduced the unresolved first-order branch transport to the four bundle observables
\[
(\Theta_w,\;K_s,\;K_q,\;P_0),
\qquad
P_0=\frac{N_0}{D_0}.
\]
That was already a major simplification, but it still left one natural question open:

> once these four bundle drifts are given, what happens to the remaining microscopic mouth/background variables
> \[
> (L_W,\;c_{s,w},\;\ell,\;\mathcal T_m,\;v_{w0},\;g_s,\;g_q,\;\lambda)
> \]
> and does the actual branch stay on the exact compensated parent family?

This note answers that question exactly.

The result is sharper than I expected:

1. all of the remaining first-order mouth/background drifts are explicit algebraic images of
   \((\delta\ln\Theta_w,\delta\ln K_s,\delta\ln K_q,\delta\ln P_0)\);
2. the outgoing-normalization drift \(\delta\ln P_0\) feeds only the traction side, not the background mixed velocity;
3. and, most importantly, the exact parent ratios
   \[
   r_c=\frac{\lambda^2}{K_sK_q},
   \qquad
   \mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
   \qquad
   \mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}
   \]
   are all **first-order invariants** under arbitrary bundle drifts.

So the Stage-149 bundle motion is not merely “close” to the compensated Family-1 branch.
It is tangent to it exactly at first order.

---

## 1. Carry-forward identities

From Stage 149,
\[
\delta\ln \rho_w
=
\frac12\,\delta\ln \Theta_w,
\]
\[
\delta\ln a
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln \Theta_w,
\]
\[
\delta\ln c_s
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln \Theta_w
+
\frac15\,\delta\ln P_0,
\]
\[
\delta\ln \mathcal Z_q
=
\delta\ln K_q
-
\frac25\,\delta\ln P_0.
\]

On the frozen `n=5` wall-EOS branch,
\[
c_{s,w}^2 \propto \rho_w^4,
\]
so
\[
\boxed{
\delta\ln c_{s,w}=2\,\delta\ln\rho_w=\delta\ln\Theta_w.
}
\]

On the healing-locked shell branch,
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\]
hence
\[
\boxed{
\delta\ln \ell=-\delta\ln c_{s,w}=-\delta\ln\Theta_w.
}
\]

And from the exact lower compensated branch,
\[
\delta\ln L_W=\delta\ln a.
\]
So
\[
\boxed{
\delta\ln L_W
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln \Theta_w.
}
\]

---

## 2. Exact bundle transport of \(\mathcal T_m\) and \(v_{w0}\)

Stage 148 already fixed the lower-branch transport laws
\[
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
+\delta\ln c_s
-\frac52\,\delta\ln a,
\]
\[
\delta\ln \mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
-\delta\ln c_s
-\frac32\,\delta\ln a.
\]

Substituting the Stage-149 inversion formulas gives
\[
\boxed{
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{13}{8}\,\delta\ln\Theta_w.
}
\]
Remarkably, the outgoing-normalization drift cancels:
\[
\boxed{
\partial_{\ln P_0}\,\delta\ln v_{w0}=0.
}
\]

For the mouth traction amplitude,
\[
\boxed{
\delta\ln \mathcal T_m
=
-\frac54\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{15}{8}\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0.
}
\]
So the first-order outgoing-normalization drift is carried entirely by the traction side, not by the background mixed velocity.

It is also useful to keep the sum/difference forms:
\[
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln\Theta_w
+
\frac25\,\delta\ln P_0,
}
\]
\[
\boxed{
\delta\ln(v_{w0}\mathcal T_m)
=
-2\,\delta\ln K_s
+\delta\ln K_q
+\frac72\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0.
}
\]

---

## 3. Exact bundle transport of the parent couplings \(g_s,g_q,\lambda\)

The explicit parent-action formulas are
\[
g_s=\mathcal T_m J_s,
\qquad
J_s\propto a^2\ell,
\]
\[
g_q\propto \mathcal Z_q\,L_W^{-3/2},
\]
\[
\lambda\propto v_{w0}\,a^2\ell\,L_W^{1/2}.
\]

Therefore
\[
\delta\ln g_s
=
\delta\ln \mathcal T_m + 2\,\delta\ln a + \delta\ln\ell,
\]
\[
\delta\ln g_q
=
\delta\ln\mathcal Z_q-\frac32\,\delta\ln L_W,
\]
\[
\delta\ln \lambda
=
\delta\ln v_{w0}+2\,\delta\ln a+\delta\ln\ell+\frac12\,\delta\ln L_W.
\]

Substituting the bundle transport laws yields
\[
\boxed{
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
}
\]
\[
\boxed{
\delta\ln g_q
=
-\frac34\,\delta\ln K_s
+\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
}
\]
and the especially clean hybridization law
\[
\boxed{
\delta\ln \lambda
=
\frac12\,\delta\ln K_s
+
\frac12\,\delta\ln K_q.
}
\]

So the bilinear shell/mixed coupling is blind, at first order, to both the wall-depth drift and the outgoing-normalization drift. It tracks only the geometric mean of the two stiffness drifts.

---

## 4. Tangent-compensation theorem

Now compute the three parent invariants.

### 4.1 Hybridization ratio \(r_c\)

Because
\[
r_c=\frac{\lambda^2}{K_sK_q},
\]
we get
\[
\delta\ln r_c
=
2\,\delta\ln\lambda
-\delta\ln K_s
-\delta\ln K_q.
\]
Using the result above,
\[
\boxed{
\delta\ln r_c=0.
}
\]

So the static/mixed hybridization ratio is fixed at first order under arbitrary bundle drifts.

### 4.2 Normalized mixed background flow \(\mathfrak r\)

Since
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\]
we obtain
\[
\delta\ln\mathfrak r
=
\delta\ln\lambda
-\frac12(\delta\ln K_s+\delta\ln K_q)
=
0.
\]
Thus
\[
\boxed{
\delta\ln\mathfrak r=0.
}
\]

### 4.3 Normalized mouth-coupling ratio \(\mathfrak g\)

Since
\[
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\]
we have
\[
\delta\ln\mathfrak g
=
\delta\ln g_q
+\frac12\,\delta\ln K_s
-\delta\ln g_s
-\frac12\,\delta\ln K_q.
\]
Substituting the explicit bundle formulas gives
\[
\boxed{
\delta\ln\mathfrak g=0.
}
\]

So both parent compensation coordinates \((\mathfrak r,\mathfrak g)\) are exact first-order invariants.

---

## 5. Exact vanishing of the off-family channels

The two Stage-146 logarithmic channels are
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathfrak g-\delta\ln\mathfrak r,
\]
\[
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
-\delta\ln r_c.
\]
Therefore
\[
\boxed{
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)=0,
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0.
}
\]

And since Stage 146 defined the off-family normal coordinate by
\[
\delta_\perp
=
\mathfrak g_*\,
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}\,
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right),
\]
we get the exact first-order result
\[
\boxed{
\delta_\perp=0.
}
\]

This is the strongest geometric interpretation so far:

> arbitrary first-order drift of the bundle observables \((\Theta_w,K_s,K_q,P_0)\) moves the actual lower compensated branch only **tangentially** inside the exact parent compensation family.

No first-order normal/off-family departure is generated inside this reduced bundle.

---

## 6. Interpretation

This stage changes the status of the remaining first-order moving-throat problem.

Before Stage 149, the unresolved statement looked like:

> compute the last four microscopic drifts.

After Stage 149, it became:

> compute the four bundle observables \((\Theta_w,K_s,K_q,P_0)\).

After the present stage, the sharper statement is:

> first-order drift of those four bundle observables is **not enough** to leave the exact compensated Family-1 branch.

That means the remaining first-order danger cannot come from generic isotropic bundle transport alone.
It must come from genuinely **off-bundle** structure, for example:

1. grouped-lane anisotropy \((20/21/22)\),
2. breakdown of the healing-lock / exact lower-branch transport identities,
3. deformation of the wall-depth extraction away from the explicit Family-1 shell-weighted law,
4. or higher-order/nonlinear corrections beyond the present first-order algebraic bundle.

So the live PDE-facing bottleneck is now narrower again:

- not first-order isotropic bundle drift,
- but the first true off-bundle correction that spoils the tangent-compensation theorem.

---

## 7. What is now fully explicit

Collecting Stages 149–150, the actual first-order branch transport of the moving-throat variables is now completely explicit in terms of the four bundle observables:
\[
(\Theta_w,\;K_s,\;K_q,\;P_0)
\longrightarrow
(\rho_w,\;c_{s,w},\;\ell,\;a,\;L_W,\;c_s,\;\mathcal Z_q,\;\mathcal T_m,\;v_{w0},\;g_s,\;g_q,\;\lambda).
\]

The exact bundle-side transport laws are:
\[
\delta\ln\rho_w=\frac12\,\delta\ln\Theta_w,
\qquad
\delta\ln c_{s,w}=\delta\ln\Theta_w,
\qquad
\delta\ln\ell=-\delta\ln\Theta_w,
\]
\[
\delta\ln a=\delta\ln L_W=\frac12\,\delta\ln K_s-\frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln c_s=\frac12\,\delta\ln K_s-\frac14\,\delta\ln\Theta_w+\frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q=\delta\ln K_q-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln v_{w0}
=
-\frac34\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{13}{8}\,\delta\ln\Theta_w,
\]
\[
\delta\ln \mathcal T_m
=
-\frac54\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac{15}{8}\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_s
=
-\frac14\,\delta\ln K_s
+\frac12\,\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln g_q
=
-\frac34\,\delta\ln K_s
+\delta\ln K_q
+\frac38\,\delta\ln\Theta_w
-\frac25\,\delta\ln P_0,
\]
\[
\delta\ln\lambda
=
\frac12(\delta\ln K_s+\delta\ln K_q).
\]

So the first-order isotropic branch transport is no longer missing at all.
What is missing is the first correction that escapes this closed bundle algebra.
