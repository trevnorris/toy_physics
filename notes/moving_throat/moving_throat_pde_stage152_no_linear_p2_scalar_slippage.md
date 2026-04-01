# Moving-Throat PDE — Stage 152: No Linear Grouped-`P2` Feed-Down into the Scalar Off-Bundle Slippages

## Purpose

Stage 151 reduced the first-order off-bundle defect to three scalar slippages
\[
\varepsilon_L,
\qquad
\varepsilon_v,
\qquad
\varepsilon_T,
\]
and showed that they matter only through the weighted combination
\[
\varepsilon_\perp
=
\mathfrak g_*\,\varepsilon_T
+\left(\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_v
+\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_L.
\]

But Stage 151 deliberately left open the microscopic source of these slippages.
One of the first candidate sources was the weak grouped-lane anisotropy of the real
`P2` bundle already isolated much earlier.

This stage answers that candidate cleanly:

> a **pure grouped real `P2` anisotropy cannot generate the scalar off-bundle slippages at linear order**.

So grouped-lane anisotropy is **not** the first linear source of
\(\varepsilon_L,\varepsilon_v,\varepsilon_T\).
Its scalar feed-down starts only at quadratic order, through the grouped invariant bilinears.

That narrows the first-order theorem gate again.
The linear grouped-anisotropy problem moves entirely into the direct outlet coefficients
\(\delta\kappa_W\) and \(\delta\gamma_W\), not into the scalar slippage channel.

---

## 1. Carry-forward grouped-lane anisotropy data

From the grouped real `P2` bundle, any grouped triple
\[
x=(x_{20},x_{21},x_{22})^T
\]
is naturally measured with the grouped metric
\[
G_{\rm grp}=\operatorname{diag}(1,2,2).
\]
The isotropic trace/anomaly variables are
\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\qquad
 a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
 b_x=\frac{x_{21}-x_{22}}{2}.
\]

The exact grouped-deviation vector is
\[
\delta x:=x-\bar x\,(1,1,1)^T,
\]
and its natural quadratic invariant is
\[
\boxed{
\mathcal I[x,y]
:=
\frac{1}{5}\,\delta x^T G_{\rm grp}\,\delta y
=
4a_x a_y+\frac{4}{5}b_x b_y.
}
\]
For one bundle this reduces to the anisotropy norm
\[
\boxed{
\mathcal A_x^2:=\mathcal I[x,x]=4a_x^2+\frac45 b_x^2.
}
\]

On the weak axisymmetric `Y_20` branch, the grouped signature is
\[
x_{20}=x^{(0)}+\epsilon x^{(1)},
\qquad
x_{21}=x^{(0)}+\frac\epsilon2 x^{(1)},
\qquad
x_{22}=x^{(0)}-\epsilon x^{(1)},
\]
so that
\[
\boxed{b_x=3a_x,}
\qquad
\boxed{
\mathcal A_x^2=\frac{7}{10}\,\epsilon^2\,(x^{(1)})^2.
}
\]
Thus the first weak axisymmetric grouped anisotropy is already quadratic as a scalar invariant.

---

## 2. Monopole-selection theorem

Let \(\mathcal S[\cdot]\) be any smooth **rotational scalar observable** extracted from the grouped-
`P2` moving-throat data around an isotropic reference branch. Examples include the mouth-average
radius, the mean throat depth, the scalar background flow speed, the scalar mouth traction, and the
scalar bundle observables that enter the lower-branch transport map.

Then the first variation of \(\mathcal S\) at the isotropic branch is an
\(O(3)\)-invariant linear functional on the `l=2` perturbation space.
But there is no nonzero invariant linear map
\[
V_{l=2}\to V_{l=0}\cong\mathbb R.
\]
So the linear functional must vanish.

Equivalently, in the explicit harmonic language, the sphere average of every real `P2` harmonic is zero:
\[
\int_{S^2}Y_{2m}^{\rm real}(\Omega)\,d\Omega=0.
\]
Therefore
\[
\boxed{
\delta^{(1)}_{P_2}\mathcal S=0
}
\]
for every scalar observable \(\mathcal S\) extracted from the isotropic branch.

This is the exact reduced statement that pure grouped real `P2` anisotropy has **no linear scalar feed-down**.

---

## 3. Immediate consequences for the moving-throat scalar observables

Apply the theorem to the scalar observables already used in the lower-branch transport map:
\[
a,
\quad
L_W,
\quad
v_{w0},
\quad
\mathcal T_m,
\quad
\rho_w,
\quad
c_{s,w},
\quad
c_s,
\quad
\mathcal Z_q.
\]
Then on a pure grouped real `P2` anisotropy,
\[
\boxed{
\delta^{(1)}_{P_2}\ln a
=
\delta^{(1)}_{P_2}\ln L_W
=
\delta^{(1)}_{P_2}\ln v_{w0}
=
\delta^{(1)}_{P_2}\ln \mathcal T_m
=
0,
}
\]
\[
\boxed{
\delta^{(1)}_{P_2}\ln\rho_w
=
\delta^{(1)}_{P_2}\ln c_{s,w}
=
\delta^{(1)}_{P_2}\ln c_s
=
\delta^{(1)}_{P_2}\ln \mathcal Z_q
=
0.
}
\]

Therefore the three Stage-151 slippages obey
\[
\boxed{
\varepsilon_L^{(1,P_2)}=0,
\qquad
\varepsilon_v^{(1,P_2)}=0,
\qquad
\varepsilon_T^{(1,P_2)}=0.
}
\]
Hence
\[
\boxed{
\varepsilon_\perp^{(1,P_2)}=0,
\qquad
\delta_\perp^{(1,P_2)}=0.
}
\]

So a pure weak grouped-lane anisotropy cannot be the first linear source of the scalar normal/off-family defect.

---

## 4. Quadratic invariant feed-down law

Once the linear term is forbidden, the first nonzero scalar feed-down from grouped-lane anisotropy is quadratic.
The natural quadratic scalar basis is built from the grouped bilinears
\[
\mathcal I[x,y]
=
\frac15\,\delta x^T G_{\rm grp}\,\delta y
=
4a_x a_y+\frac45 b_x b_y.
\]

Let
\[
\mathcal B:=\{D_0,D_2,D_4,N_0,N_2,N_4\}
\]
be the reduced grouped-lane families already isolated by the coupled wall/BdG/Maxwell/mixed bundle.
Then the most general first nonzero grouped-anisotropy contribution to the three scalar slippages has the form
\[
\boxed{
\varepsilon_L
=
\sum_{X,Y\in\mathcal B}\Xi_L^{(XY)}\,\mathcal I[X,Y]
+O(3),
}
\]
\[
\boxed{
\varepsilon_v
=
\sum_{X,Y\in\mathcal B}\Xi_v^{(XY)}\,\mathcal I[X,Y]
+O(3),
}
\]
\[
\boxed{
\varepsilon_T
=
\sum_{X,Y\in\mathcal B}\Xi_T^{(XY)}\,\mathcal I[X,Y]
+O(3),
}
\]
where the omitted terms are cubic and higher in the grouped-anisotropy amplitudes.

If one wants the more microscopic decomposition, each `D_n` can then be expanded further into the wall, BdG, and conservative Maxwell/mixed families from the full grouped bundle.

On a one-parameter weak axisymmetric branch,
\[
\mathcal I[X,Y]
=
\frac{7}{10}\,\epsilon^2 X^{(1)}Y^{(1)}.
\]
So the scalar feed-down is automatically quadratic in the weak axisymmetric anisotropy amplitude.

---

## 5. Exact consequence for the Stage-151 normal defect ledger

Using the Stage-151 weighting,
\[
\varepsilon_\perp
=
\mathfrak g_*\,\varepsilon_T
+\left(\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_v
+\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_L,
\]
we get the quadratic grouped-anisotropy transport law
\[
\boxed{
\varepsilon_\perp
=
\sum_{X,Y\in\mathcal B}\Xi_\perp^{(XY)}\,\mathcal I[X,Y]
+O(3),
}
\]
with
\[
\boxed{
\Xi_\perp^{(XY)}
=
\mathfrak g_*\,\Xi_T^{(XY)}
+\left(\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\Xi_v^{(XY)}
+\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\Xi_L^{(XY)}.
}
\]

Numerically, on the renormalized Family-1 point,
\[
\boxed{
\Xi_\perp^{(XY)}
\approx
0.758035078944663\,\Xi_T^{(XY)}
+1.00314310113848\,\Xi_v^{(XY)}
+1.88373219118005\,\Xi_L^{(XY)}.
}
\]

Therefore
\[
\boxed{
\delta_\perp
=
-\sum_{X,Y\in\mathcal B}\Xi_\perp^{(XY)}\,\mathcal I[X,Y]
+O(3).
}
\]

So the grouped-anisotropy contribution to the scalar off-family defect is necessarily **quadratic**.

---

## 6. Mouth-bias and conservative-even consequences

The Stage-151 mouth-bias law becomes
\[
\delta\Pi
=
\delta\Pi_{\rm tan}
-
\frac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_*^2}}
\sum_{X,Y\in\mathcal B}\Xi_\perp^{(XY)}\,\mathcal I[X,Y]
+O(3).
\]

Likewise the conservative-even outlet defect becomes
\[
\delta\mathcal C
=
-
\frac{16\sigma_*}{\sqrt{1+\mathfrak r_*^2}}
\sum_{X,Y\in\mathcal B}\Xi_\perp^{(XY)}\,\mathcal I[X,Y]
+O(3).
\]

So grouped-lane anisotropy does not enter the scalar slippage channel until quadratic order.
At linear order, the only grouped-anisotropy dangers that remain are the **direct outlet coefficients**
\[
\delta\kappa_W,
\qquad
\delta\gamma_W,
\]
not the scalar slippages.

This is the key narrowing achieved by the stage.

---

## 7. What Stage 152 changes

Before this step, Stage 151 still allowed the possibility that grouped real `P2` anisotropy might directly generate the scalar slippages
\((\varepsilon_L,\varepsilon_v,\varepsilon_T)\) at first order.

After this step, that possibility is gone.

The theorem status is now:

1. **pure grouped-lane anisotropy cannot linearly source the scalar off-bundle slippages;**
2. its scalar feed-down begins only through the quadratic grouped invariants \(\mathcal I[X,Y]\);
3. therefore the remaining **linear** grouped-anisotropy bottleneck is entirely in the direct outlet coefficients
   \(\delta\kappa_W\) and \(\delta\gamma_W\),
   not in \(\varepsilon_\perp\).

So the next smallest theorem gate is no longer

> derive \((\varepsilon_L,\varepsilon_v,\varepsilon_T)\) from weak grouped-lane anisotropy.

That part is finished.

The next honest gate is

> derive the direct weak grouped-anisotropy contribution to the outlet coefficients
> \(\delta\kappa_W\) and \(\delta\gamma_W\),
> and only after that return to the quadratic scalar feed-down coefficients
> \(\Xi_I^{(XY)}\).

That is a genuine narrowing of the moving-throat PDE bottleneck.
