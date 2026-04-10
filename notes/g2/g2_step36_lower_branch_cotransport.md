# Step 36 — Exact lower-branch co-transport and the collapse to four irreducible drifts

## Goal

Step 35 showed that the compensated Family-1 outlet is not a delicate fit: the
lower branch is the unique regular positive-source branch, and it sits close to
several very natural mouth-source laws.

The next natural question is then sharper:

> once we sit on that exact lower compensated branch, how many microscopic drifts
> are still genuinely free?

This step answers that by solving the exact lower-branch transport constraints.
The main result is that the branch already co-transports most of the apparent
microscopic freedom, and the old first-order off-family defect vanishes
identically on the exact lower branch.

---

## Step 36A — The two exact lower-branch channel constraints

On the explicit lower compensated branch, the parent ratios

```math
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\qquad
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}}
```

are fixed to first order. So the two logarithmic imbalance channels must vanish:

```math
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)=0,
\qquad
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)=0.
```

Using the explicit throat-core branch formulas from the moving-throat notes, the
corresponding drift equations are

```math
0=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
```

```math
0=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
```

So the lower branch is already imposing two exact co-transport conditions before
we solve any further PDE detail.

---

## Step 36B — The D/N geometry removes `L_W` as an independent drift

The exact D/N side-tube law is

```math
L_W=
\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r_*^2}{3}}.
```

Since `\mathfrak r_*` is fixed on the exact lower branch, differentiation gives

```math
\boxed{
\delta\ln L_W = \delta\ln a.
}
```

So the side-tube length is not an independent linearized variable on the lower
branch. It co-transports with the mouth radius.

---

## Step 36C — Exact drift laws for `v_{w0}` and `\mathcal T_m`

Substituting `\delta\ln L_W=\delta\ln a` into the two branch constraints and
solving gives

```math
\boxed{
\delta\ln v_{w0}
=
\frac12
\Big(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
+2\,\delta\ln c_s
-5\,\delta\ln a
\Big),
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
\frac12
\Big(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
-2\,\delta\ln c_s
-3\,\delta\ln a
\Big).
}
```

So the background mixed flow and the mouth traction are not independent branch
variables at first order. Once the shell/gauge/background drifts are fixed, the
branch already fixes both of them.

---

## Step 36D — Exact product/ratio factorization

The two drift laws split into a particularly clean product/ratio form.
Subtracting and adding them gives

```math
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
2\,\delta\ln c_s - \delta\ln a,
}
```

```math
\boxed{
\delta\ln\!(v_{w0}\mathcal T_m)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-4\,\delta\ln a.
}
```

This is very informative.

- The **ratio** `v_{w0}/\mathcal T_m` only remembers the wave-speed/geometry side.
- The **product** `v_{w0}\mathcal T_m` only remembers the shell/localization side.

So the branch is already separating the two physical channels for us.

---

## Step 36E — Frozen `n=5` wall-EOS reduction

On the frozen wall GNLS branch,

```math
c_{s,w}^2\propto \rho_w^4,
```

so

```math
\delta\ln c_{s,w}=2\,\delta\ln\rho_w.
```

Then the exact lower-branch drift laws reduce to

```math
\boxed{
\delta\ln v_{w0}
=
\frac12\Big(
\delta\ln\mathcal Z_q
+5\,\delta\ln\rho_w
+2\,\delta\ln c_s
-5\,\delta\ln a
\Big),
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
\frac12\Big(
\delta\ln\mathcal Z_q
+5\,\delta\ln\rho_w
-2\,\delta\ln c_s
-3\,\delta\ln a
\Big).
}
```

And the product/ratio channels become

```math
\boxed{
\delta\ln\!\left(\frac{v_{w0}}{\mathcal T_m}\right)
=
2\,\delta\ln c_s-\delta\ln a,
}
```

```math
\boxed{
\delta\ln\!(v_{w0}\mathcal T_m)
=
\delta\ln\mathcal Z_q+5\,\delta\ln\rho_w-4\,\delta\ln a.
}
```

So after the exact branch constraints and the frozen `n=5` wall-EOS law are both
used, the apparently large drift space collapses to only four irreducible
microscopic directions:

```math
\boxed{
(\delta\ln\mathcal Z_q,\;\delta\ln\rho_w,\;\delta\ln c_s,\;\delta\ln a).
}
```

Everything else is already co-transported by the exact lower branch.

---

## Step 36F — The old off-family scalar vanishes identically on the exact lower branch

The earlier moving-throat notes reduced the first off-family defect to the scalar

```math
\delta_\perp
=
\mathfrak g_*\,\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right).
```

But on the exact lower compensated branch, those two channels are precisely the
constraints we solved above, so

```math
\boxed{
\delta_\perp=0
}
```

identically.

That is the honest reason the old first-order defect disappears on the exact
lower branch: the lower-branch transport laws already enforce its cancellation.

---

## Main result of the step

The compensated Family-1 outlet is now much sharper than it looked even one step
ago.

The exact lower branch already forces

```math
\boxed{
\delta\ln L_W=\delta\ln a,
}
```

and fixes the mixed-flow and traction drifts by

```math
\boxed{
\delta\ln v_{w0}
=
\frac12(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
+2\,\delta\ln c_s
-5\,\delta\ln a),
}
```

```math
\boxed{
\delta\ln\mathcal T_m
=
\frac12(
\delta\ln\mathcal Z_q
-\delta\ln\rho_w
+3\,\delta\ln c_{s,w}
-2\,\delta\ln c_s
-3\,\delta\ln a).
}
```

After the frozen `n=5` wall-EOS reduction, the entire actual branch-drift problem
collapses to only

```math
\boxed{
(\delta\ln\mathcal Z_q,\;\delta\ln\rho_w,\;\delta\ln c_s,\;\delta\ln a),
}
```

and the old first-order off-family scalar vanishes identically.

So the next question is no longer “is there still a first-order branch defect?”
There is not. The sharper PDE-facing question is:

> what does the true moving-throat solution do to those **four** residual lower-branch variables?
