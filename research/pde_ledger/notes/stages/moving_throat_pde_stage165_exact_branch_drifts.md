# Moving-Throat PDE — Stage 165: Exact Lower-Branch Drift Laws for `L_W`, `v_{w0}`, and `\mathcal T_m`

## Goal

Stage 164 reduced the first off-family defect to the explicit microscopic logarithmic channels
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right),
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right),
\]
but it still left the actual branch drifts of
\[
(\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W)
\]
looking broader than they really are.

This stage computes the part that the explicit coupled throat-core branch already fixes.

The result is that on the exact lower compensated Family-1 branch:

1. the D/N tube length co-transports with the mouth radius,
   \[
   \delta\ln L_W=\delta\ln a,
   \]
2. the background mixed flow and the mouth traction are not independent,
3. their product and ratio split into two exact transport channels,
4. and after the frozen `n=5` wall-EOS reduction the whole actual branch-drift problem collapses to only four independent microscopic drifts:
   \[
   (\delta\ln\mathcal Z_q,\ \delta\ln\rho_w,\ \delta\ln c_s,\ \delta\ln a).
   \]

So this really does compute the “actual drifts” that Stage 164 asked for, up to the irreducible microscopic data that only the full PDE can still choose.

---

## 1. Lower-branch conditions already fixed

From the earlier lower compensated Family-1 analysis, the exact parent ratios
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}
\]
are pinned at first order on the exact lower branch:
\[
\delta\mathfrak g=0,
\qquad
\delta\mathfrak r=0.
\]

Equivalently, the two Stage 164 logarithmic imbalance channels vanish on the true lower branch:
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)=0,
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0.
\]

So the remaining task is not to hunt for new branch constraints.
It is to solve these exact constraints in terms of the explicit microscopic variables.

---

## 2. Exact D/N geometric co-transport

The exact D/N mixed-tube selection law is
\[
L_W=\frac{\pi a}2\sqrt{\frac{1+\mathfrak r^2}3}.
\]
Differentiate at fixed \(\mathfrak r\):
\[
\boxed{
\delta\ln L_W=\delta\ln a.
}
\]

So the first actual branch drift is immediate:
the auxiliary D/N tube length is not an independent linearized variable on the lower branch.

---

## 3. Exact background-flow transport from fixed \(\mathfrak r\)

On the explicit throat-core branch,
\[
\lambda=-q_*v_{w0}\mathcal I_{sq},
\qquad
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\]
and under the carried uniform-overlap + healing-lock closure,
\[
\mathcal I_{sq}\propto a^2\ell L_W^{1/2},
\qquad
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\qquad
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\ell}.
\]

Using \(\delta\mathfrak r=0\) gives
\[
0
=
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
\]

Solving,
\[
\boxed{
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
+\delta\ln c_s
-\delta\ln a
-\frac32\,\delta\ln L_W.
}
\]

Using the D/N co-transport law \(\delta\ln L_W=\delta\ln a\),
\[
\boxed{
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_s^2 c_{s,w}^3}{\rho_w a^5}\right).
}
\]

So the actual background mixed-flow drift is already determined by the remaining shell/gauge/geometry drifts.

---

## 4. Exact traction transport from fixed \(\mathfrak g\)

The exact mouth traction law on the lower compensated branch is
\[
\mathcal T_m
=
\frac{\sqrt{2\mathcal Z_q K_s}}{J_s c_s \sqrt{\mu_0 L_W}}
\frac{1}{\mathfrak g_*},
\qquad
J_s=\frac{4\pi a^2\ell}3.
\]

Using the same healing-lock shell reduction and \(\delta\mathfrak g=0\),
\[
0
=
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W.
\]

Solving directly for the mouth traction,
\[
\boxed{
\delta\ln\mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
-\delta\ln c_s
-\delta\ln a
-\frac12\,\delta\ln L_W.
}
\]

Then the D/N co-transport law gives
\[
\boxed{
\delta\ln\mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_{s,w}^3}{\rho_w c_s^2 a^3}\right).
}
\]

So the actual traction drift is also already fixed by the same reduced microscopic data.

---

## 5. Exact branch formulas before differentiation

The previous two sections integrate to exact lower-branch amplitude laws.

Using the fixed lower compensated Family-1 values
\[
\mathfrak r_* \approx 1.77799353547498,
\qquad
\mathfrak g_* \approx 0.758035078944663,
\]
the D/N law gives
\[
L_W=\frac{\pi a}2\sqrt{\frac{1+\mathfrak r_*^2}3}.
\]

Substituting this into the explicit throat-core formulas yields
\[
\boxed{
\mathcal T_m
=
\frac{3\sqrt{10}\,3^{3/4}}{5\pi\mathfrak g_*(1+\mathfrak r_*^2)^{1/4}}
\frac{m_\psi}{\sqrt{\hbar\mu_0}}
\frac{\sqrt{\mathcal Z_q}\,c_{s,w}^{3/2}}{a^{3/2}c_s\sqrt{\rho_w}},
}
\]
\[
\boxed{
v_{w0}
=
-\frac{9\sqrt{10}\,3^{1/4}\mathfrak r_*}{20(1+\mathfrak r_*^2)^{3/4}}
\frac{m_\psi}{q_*\sqrt{\hbar\mu_0}}
\frac{\sqrt{\mathcal Z_q}\,c_s c_{s,w}^{3/2}}{a^{5/2}\sqrt{\rho_w}}.
}
\]

Numerically, the pure branch constants are
\[
\boxed{
\mathcal T_m
\approx
1.2715890393387603\,
\frac{m_\psi}{\sqrt{\hbar\mu_0}}
\frac{\sqrt{\mathcal Z_q}\,c_{s,w}^{3/2}}{a^{3/2}c_s\sqrt{\rho_w}},
}
\]
\[
\boxed{
v_{w0}
\approx
-1.1428896163056477\,
\frac{m_\psi}{q_*\sqrt{\hbar\mu_0}}
\frac{\sqrt{\mathcal Z_q}\,c_s c_{s,w}^{3/2}}{a^{5/2}\sqrt{\rho_w}}.
}
\]

So the actual lower-branch amplitudes are now explicit rather than merely logarithmic.

---

## 6. Exact product/ratio factorization

A very useful simplification appears when the two branch formulas are combined.

Their ratio is exactly
\[
\boxed{
\frac{v_{w0}}{\mathcal T_m}
=
-\frac{\sqrt3\,\pi\,\mathfrak g_* \mathfrak r_*}{4q_*\sqrt{1+\mathfrak r_*^2}}
\frac{c_s^2}{a}.
}
\]

So the **ratio** of background mixed flow to mouth traction carries only the mixed-tube wave-speed/geometry information.

Its logarithmic drift is therefore
\[
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
2\,\delta\ln c_s
-
\delta\ln a.
}
\]

Likewise the product is exactly
\[
\boxed{
v_{w0}\mathcal T_m
=
-\frac{81\mathfrak r_*}{10\pi\mathfrak g_*(1+\mathfrak r_*^2)}
\frac{m_\psi^2}{\hbar\mu_0 q_*}
\frac{\mathcal Z_q c_{s,w}^3}{\rho_w a^4}.
}
\]

So the **product** carries only the shell/localization sector.
Its logarithmic drift is
\[
\boxed{
\delta\ln(v_{w0}\mathcal T_m)
=
\delta\ln\mathcal Z_q
+
3\,\delta\ln c_{s,w}
-
\delta\ln\rho_w
-
4\,\delta\ln a.
}
\]

Numerically, the pure lower-branch constants are
\[
\boxed{
\frac{v_{w0}}{\mathcal T_m}
\approx
-\frac{0.8987885086678338}{q_*}
\frac{c_s^2}{a},
}
\]
\[
\boxed{
v_{w0}\mathcal T_m
\approx
-\frac{1.4532859092683434}{q_*}
\frac{m_\psi^2}{\hbar\mu_0}
\frac{\mathcal Z_q c_{s,w}^3}{\rho_w a^4}.
}
\]

This split is one of the cleanest new results of the stage.

---

## 7. Frozen `n=5` wall-EOS reduction

On the frozen GNLS branch,
\[
c_{s,w}^2 \propto \rho_w^4,
\]
so
\[
\boxed{
\delta\ln c_{s,w} = 2\,\delta\ln\rho_w.
}
\]

Substituting into the exact lower-branch drift laws gives
\[
\boxed{
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_s^2 \rho_w^5}{a^5}\right),
}
\]
\[
\boxed{
\delta\ln\mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q \rho_w^5}{c_s^2 a^3}\right).
}
\]

And the product/ratio channels reduce to
\[
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
2\,\delta\ln c_s-\delta\ln a,
}
\]
\[
\boxed{
\delta\ln(v_{w0}\mathcal T_m)
=
\delta\ln\mathcal Z_q
+
5\,\delta\ln\rho_w
-
4\,\delta\ln a.
}
\]

So after the exact lower-branch and `n=5` reductions, the apparently large drift space has collapsed to four irreducible microscopic coordinates:
\[
(\delta\ln\mathcal Z_q,\ \delta\ln\rho_w,\ \delta\ln c_s,\ \delta\ln a).
\]

Everything else is co-transported.

---

## 8. Collapse of the Stage 164 off-family defect

Stage 164 expressed the microscopic off-family normal coordinate in terms of the two imbalance channels
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right),
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right).
\]

But these are exactly the lower-branch fixed-\(\mathfrak g\) and fixed-\(\mathfrak r\) conditions used above.
So after substituting the branch drift laws derived here,
\[
\boxed{
\delta_\perp=0
}
\]
identically.

That is the honest reason the Stage 164 normal coordinate disappears on the exact lower compensated branch:
the branch drift constraints themselves already compute the cancellation.

---

## 9. Best current theorem statement after Stage 165

The explicit lower compensated Family-1 branch now determines the linearized drifts of

- `L_W`,
- `v_{w0}`,
- and `\mathcal T_m`

exactly in terms of the remaining microscopic shell/gauge/geometry variables.

After also using the frozen `n=5` wall-EOS identity, the actual branch-drift problem is no longer eight-dimensional. It is only four-dimensional:

\[
(\mathcal Z_q,\rho_w,c_s,a).
\]

So the next PDE-facing question is now as small as it has ever been:

> what does the true moving-throat solution do to the four irreducible lower-branch variables
> \[
> (\mathcal Z_q,\rho_w,c_s,a),
> \]
> given that the exact compensated branch already co-transports
> `L_W`, `v_{w0}`, `\mathcal T_m`, and `c_{s,w}` for us?

That is a much sharper continuation target than Stage 164 left us with.
