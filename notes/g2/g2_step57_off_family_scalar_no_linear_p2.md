# Step 57 — The first off-family defect is one scalar, and linear grouped-`P2` anisotropy cannot generate it

## Goal

Step 56 showed that the exact parent compensation family forces

```math
\Xi_{\rm slip}=0,
\qquad
\delta\mathfrak B_W=0,
\qquad
\delta\gamma_W=0,
\qquad
\Delta_Q=0
```

at first order on the natural compensated lower branch.

The next honest question is:

> if a first-order defect does survive somewhere, what is the **smallest** object that can carry it?

The later moving-throat notes answer that cleanly.
They show that all first-order off-family structure collapses to one weighted scalar slippage combination, and then show that pure grouped real `P2` anisotropy cannot feed that scalar at linear order.

So the first-order theorem gate narrows again.

---

## Step 57A — One exact scalar carries the off-family drift

The three first-order transport slippages are

- `\varepsilon_L` for D/N geometric similarity,
- `\varepsilon_v` for the mixed background-flow transport law,
- `\varepsilon_T` for the mouth-traction transport law.

But they matter only through the one weighted combination

```math
\boxed{
\varepsilon_\perp
=
\mathfrak g_*\,\varepsilon_T
+
\left(\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_v
+
\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_L.
}
```

and the exact off-family normal coordinate is simply

```math
\boxed{\delta_\perp = -\varepsilon_\perp.}
```

So the first off-family defect is not a large vector of drifts anymore.
It is one scalar.

At the Family-1 point,

```math
\mathfrak r_* \approx 1.77799353547498,
\qquad
\mathfrak g_* \approx 0.758035078944663,
```

so numerically

```math
\boxed{
\delta_\perp
\approx
-0.758035078944663\,\varepsilon_T
-1.00314310113848\,\varepsilon_v
-1.88373219118005\,\varepsilon_L.
}
```

The strongest weight falls on the geometric slippage `\varepsilon_L`, then the mixed-speed slippage `\varepsilon_v`, then the traction slippage `\varepsilon_T`.

---

## Step 57B — The outgoing defect depends only on `\varepsilon_\perp` and `\delta\gamma_W`

The exact first-order outgoing defect ledger becomes

```math
\boxed{
\Delta_Q
=
\frac{\sigma_*}{3(1-\sigma_*)}
\left[
-\frac{16}{\sqrt{1+\mathfrak r_*^2}}\,\varepsilon_\perp
-27\,\delta\gamma_W
\right].
}
```

So there are only two linear routes left:

1. a true off-family scalar slippage `\varepsilon_\perp`,
2. a direct odd mixed-port renormalization `\delta\gamma_W`.

If the canonical even branch is preserved, then the same moving-throat theorem again gives

```math
\varepsilon_\perp=0,
```

leaving only the odd mixed-port channel.

---

## Step 57C — Weak axisymmetric grouped-`P2` anisotropy is quadratic as a scalar invariant

For a weak axisymmetric grouped signature

```math
x_{20}=x^{(0)}+\epsilon x^{(1)},
\qquad
x_{21}=x^{(0)}+\frac\epsilon2 x^{(1)},
\qquad
x_{22}=x^{(0)}-\epsilon x^{(1)},
```

one gets the exact grouped variables

```math
\bar x = x^{(0)},
\qquad
a_x = \frac\epsilon4 x^{(1)},
\qquad
b_x = \frac{3\epsilon}{4} x^{(1)}.
```

So the weak axisymmetric signature always satisfies

```math
\boxed{b_x = 3 a_x.}
```

The first scalar grouped invariant is

```math
\mathcal A_x^2 = 4 a_x^2 + \frac45 b_x^2,
```

which reduces exactly to

```math
\boxed{
\mathcal A_x^2 = \frac{7}{10}\,\epsilon^2 (x^{(1)})^2.
}
```

So the first scalar invariant is **quadratic**, not linear.

That means pure grouped real `P2` anisotropy cannot generate any scalar feed-down at linear order.
In particular,

```math
\boxed{
\varepsilon_L^{(1,P_2)}=0,
\qquad
\varepsilon_v^{(1,P_2)}=0,
\qquad
\varepsilon_T^{(1,P_2)}=0,
\qquad
\varepsilon_\perp^{(1,P_2)}=0.
}
```

---

## Main result of the step

The first-order g-2 bottleneck is now smaller than before.

A nonzero first-order defect cannot be blamed on a generic grouped-`P2` anisotropy, because pure grouped real `P2` anisotropy has no linear scalar feed-down.

So the only remaining linear possibilities are:

- a genuine off-family scalar slippage `\varepsilon_\perp`,
- or a direct odd mixed-port/output renormalization `\delta\gamma_W`.

That is a much sharper continuation target than “look for any anisotropy somewhere.”
