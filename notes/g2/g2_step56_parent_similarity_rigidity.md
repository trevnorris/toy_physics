# Step 56 — The exact parent compensation family preserves D/N similarity automatically

## Goal

Step 55 ended with the strongest reduced statement so far:

- the natural compensated branch does **not** generate the electron sliver at first order,
- but that result was still phrased in terms of the reduced similarity-slippage scalar
  `Xi_slip`.

The later moving-throat notes go one level deeper.
They show that on the exact parent compensation family,

```math
\gamma_0 = \frac{1+r^2}{9},
\qquad
\frac{L_W}{a} = \frac{\pi}{2}\sqrt{\frac{1+r^2}{3}},
\qquad
g_-(r)=r-\frac12\sqrt{1+r^2},
```

so the D/N similarity law is not an extra assumption at all.
It is an **identity** along the full parent family.

This step folds that exact rigidity back into the g-2 chain.

---

## Step 56A — Exact parent-family similarity identity

Differentiate the exact parent-family formulas:

```math
\ln \gamma_0 = \ln(1+r^2)-\ln 9,
```

```math
\ln\left(\frac{L_W}{a}\right)
=
\ln\frac{\pi}{2}-\frac12\ln 3 + \frac12\ln(1+r^2).
```

Then

```math
\delta\ln\gamma_0
=
\frac{2r}{1+r^2}\,\delta r,
\qquad
\delta\ln\left(\frac{L_W}{a}\right)
=
\frac{r}{1+r^2}\,\delta r.
```

So exactly

```math
\boxed{
\delta\ln\gamma_0 - 2\,\delta\ln\left(\frac{L_W}{a}\right)=0.
}
```

Equivalently,

```math
\boxed{\Xi_{\rm slip}=\Xi_\gamma-2\Xi_L=0}
```

on the entire exact parent compensation family.

So the Stage-55 “similarity preservation” is not extra structure. It is already built into the exact parent family.

---

## Step 56B — Lower-branch rigidity under the canonical-even gate

On the lower compensated branch,

```math
g_-(r)=r-\frac12\sqrt{1+r^2},
```

and

```math
\frac{dg_-}{dr}
=
1-\frac{r}{2\sqrt{1+r^2}}
=
\frac{4+3r^2}{2\sqrt{1+r^2}(2\sqrt{1+r^2}+r)}
>0.
```

So if the carried canonical-even gate still gives

```math
\delta g = 0,
```

then the lower branch is first-order rigid:

```math
\boxed{\delta r = 0.}
```

At the Family-1 point `r ≈ 1.77799353547498`, the branch slope is

```math
\frac{dg_-}{dr} \approx 0.564199521046343,
```

so the rigidity is numerically strong rather than marginal.

---

## Step 56C — First-order outgoing defect collapses to zero

Once `Xi_slip = 0`, the reduced first-order defect law becomes

```math
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}
=0.
```

And since the bare mixed-port slippage

```math
\delta\mathfrak B_W
:=
\delta\gamma_0-\frac13\delta\kappa_0
```

also vanishes identically on the parent family, the renormalized odd outlet shift

```math
\delta\gamma_W = \frac{\delta\mathfrak B_W}{1+r^2}
```

vanishes too.

So at first order on the natural compensated lower branch,

```math
\boxed{
\delta\mathfrak B_W = 0,
\qquad
\delta\gamma_W = 0,
\qquad
\Delta_Q = 0,
\qquad
N_Q-1 = 0.
}
```

---

## Main result of the step

The natural compensated lower branch does not merely *allow* zero outgoing slippage.
It **forces** it at first order because the exact parent compensation family already preserves the D/N similarity law automatically.

So the g-2 picture is now sharper:

- the canonical outgoing normalization is a natural exact first-order consequence of the parent family,
- and any nonzero electron sliver must come from a genuine departure from that family, not from the natural branch itself.
