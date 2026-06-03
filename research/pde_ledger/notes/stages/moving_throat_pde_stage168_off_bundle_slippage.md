# Moving-Throat PDE — Stage 168: Exact Off-Bundle Slippage Decomposition

## Purpose

Stage 167 proved that arbitrary **first-order isotropic bundle drift**
does not move the co-evolving Family-1 system off the exact parent compensation family.
So the next live question is no longer

> what do the bundle drifts do?

That part is finished.

The real next question is sharper:

> what is the first correction that escapes the closed bundle algebra and actually generates the normal/off-family defect?

This stage answers that at first order.

The main result is that the entire first-order off-bundle defect collapses to **three scalar slippage variables** measuring failure of the exact lower-branch transport laws:

1. axial-length slippage,
2. mixed-background-speed slippage,
3. mouth-traction slippage.

Those three slippages enter the normal coordinate only through one weighted scalar combination.
So the first-order off-bundle problem is smaller than “all microscopic drifts.”
It is a three-parameter defect ledger.

---

## 1. Carry-forward formulas

From Stage 164, the off-family normal coordinate is
\[
\delta_\perp
=
A_*\,
\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+3A_*\,\delta\ln c_{s,w}
+B_*\,\delta\ln c_s
-\mathfrak g_*\,\delta\ln\mathcal T_m
-(\mathfrak g_*+B_*)\,\delta\ln v_{w0}
-2A_*\,\delta\ln a
-C_*\,\delta\ln L_W,
\]
with
\[
A_*=\mathfrak g_*+\frac{1}{4\sqrt{1+\mathfrak r_*^2}},
\qquad
B_*=\frac{1}{2\sqrt{1+\mathfrak r_*^2}},
\qquad
C_*=2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}.
\]

From Stage 165, the exact lower-branch transport laws are
\[
\delta\ln L_W=\delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
+\delta\ln c_s
-\frac52\,\delta\ln a,
\]
\[
\delta\ln\mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
-\delta\ln c_s
-\frac32\,\delta\ln a.
\]

Stage 167 then showed that substituting those laws into the Stage 164 normal coordinate gives
\[
\delta_\perp=0.
\]

So any nonzero first-order \(\delta_\perp\) must come from **violations** of the lower-branch transport laws themselves.

---

## 2. Define the three off-bundle slippage variables

Introduce the exact first-order slippages
\[
\boxed{
\varepsilon_L
:=
\delta\ln L_W-\delta\ln a,
}
\]
\[
\boxed{
\varepsilon_v
:=
\delta\ln v_{w0}
-
\left[
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
+\delta\ln c_s
-\frac52\,\delta\ln a
\right],
}
\]
\[
\boxed{
\varepsilon_T
:=
\delta\ln\mathcal T_m
-
\left[
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+\frac32\,\delta\ln c_{s,w}
-\delta\ln c_s
-\frac32\,\delta\ln a
\right].
}
\]

Interpretation:

- \(\varepsilon_L\) measures departure from the exact lower-branch D/N geometric similarity law \(L_W\propto a\);
- \(\varepsilon_v\) measures slippage of the mixed background flow away from its exact lower-branch transport law;
- \(\varepsilon_T\) measures slippage of the mouth traction away from its exact lower-branch transport law.

These are the first genuinely off-bundle first-order variables.

---

## 3. Exact collapse of the normal coordinate

Substituting
\[
\delta\ln L_W=\delta\ln a+\varepsilon_L,
\]
\[
\delta\ln v_{w0}
=
\left(\delta\ln v_{w0}\right)_{\rm br}
+\varepsilon_v,
\qquad
\delta\ln\mathcal T_m
=
\left(\delta\ln\mathcal T_m\right)_{\rm br}
+\varepsilon_T,
\]
into the Stage 164 formula gives the exact identity
\[
\boxed{
\delta_\perp
=
-\mathfrak g_*\,\varepsilon_T
-\left(\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_v
-\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_L.
}
\]

So all first-order off-bundle normal motion is carried by one weighted scalar combination of the three slippages.

Define
\[
\boxed{
\varepsilon_\perp
:=
\mathfrak g_*\,\varepsilon_T
+\left(\mathfrak g_*+\frac{1}{2\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_v
+\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\varepsilon_L.
}
\]
Then the exact first-order normal coordinate is simply
\[
\boxed{
\delta_\perp=-\varepsilon_\perp.
}
\]

This is the cleanest reduced result of the stage:

> the first-order off-bundle defect is not a large vector of microscopic drifts;
> it is one scalar normal-slippage combination \(\varepsilon_\perp\).

---

## 4. Numerical Family-1 coefficients

On the renormalized co-evolving canonical point,
\[
\mathfrak r_*=\mathfrak r_{F1}\approx 1.77799353547498,
\qquad
\mathfrak g_*\approx 0.758035078944663.
\]

Therefore
\[
\boxed{
\delta_\perp
\approx
-0.758035078944663\,\varepsilon_T
-1.00314310113848\,\varepsilon_v
-1.88373219118005\,\varepsilon_L.
}
\]

So the first-order off-family defect is weighted most strongly by the axial-length slippage,
next by mixed-speed slippage,
and least by mouth-traction slippage.

---

## 5. Mouth-bias transport with off-bundle slippage

Stage 163 already split the mouth-bias variation as
\[
\delta\Pi
=
\delta\Pi_{\rm tan}
+
\frac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_*^2}}\,\delta_\perp,
\]
with
\[
\delta\Pi_{\rm tan}
=
\left(1-\frac14\mathcal S_{\rm can}\right)\delta\Sigma_0
-\frac{\Sigma_0^{\rm can}}{4}\delta\mathcal S.
\]

So using \(\delta_\perp=-\varepsilon_\perp\),
\[
\boxed{
\delta\Pi
=
\delta\Pi_{\rm tan}
-
\frac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_*^2}}\,
\varepsilon_\perp.
}
\]

Numerically,
\[
\delta\Pi_{\rm tan}
=
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S,
\]
and
\[
\boxed{
\delta\Pi
\approx
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
-
1.52843317823248\,\varepsilon_\perp.
}
\]

Equivalently,
\[
\delta\Pi
\approx
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
-
1.15860596492310\,\varepsilon_T
-
1.53323719829507\,\varepsilon_v
-
2.87915877990416\,\varepsilon_L.
\]

So the same scalar off-bundle combination \(\varepsilon_\perp\) that drives the normal family defect also drives the non-tangential mouth-bias defect.

---

## 6. Outlet defects in terms of the slippage scalar

Stage 163 already gave
\[
\delta\mathcal C
=
\frac{16\sigma_*}{\sqrt{1+\mathfrak r_*^2}}\,\delta_\perp,
\]
\[
\delta E_2
=
\frac{\sigma_*}{27(1-\sigma_*)}
\left[
\frac{16}{\sqrt{1+\mathfrak r_*^2}}\delta_\perp
-
9\,\delta\kappa_W
\right],
\]
\[
\delta E_4
=
\frac{\sigma_*}{243(1-\sigma_*)}
\left[
\frac{80}{\sqrt{1+\mathfrak r_*^2}}\delta_\perp
-
72\,\delta\kappa_W
\right],
\]
\[
\Delta_Q
=
\frac{\sigma_*}{3(1-\sigma_*)}
\left[
\frac{16}{\sqrt{1+\mathfrak r_*^2}}\delta_\perp
-
27\,\delta\gamma_W
\right].
\]

Using \(\delta_\perp=-\varepsilon_\perp\), these become
\[
\boxed{
\delta\mathcal C
=
-\frac{16\sigma_*}{\sqrt{1+\mathfrak r_*^2}}\,\varepsilon_\perp,
}
\]
\[
\boxed{
\delta E_2
=
\frac{\sigma_*}{27(1-\sigma_*)}
\left[
-\frac{16}{\sqrt{1+\mathfrak r_*^2}}\varepsilon_\perp
-
9\,\delta\kappa_W
\right],
}
\]
\[
\boxed{
\delta E_4
=
\frac{\sigma_*}{243(1-\sigma_*)}
\left[
-\frac{80}{\sqrt{1+\mathfrak r_*^2}}\varepsilon_\perp
-
72\,\delta\kappa_W
\right],
}
\]
\[
\boxed{
\Delta_Q
=
\frac{\sigma_*}{3(1-\sigma_*)}
\left[
-\frac{16}{\sqrt{1+\mathfrak r_*^2}}\varepsilon_\perp
-
27\,\delta\gamma_W
\right].
}
\]

So the first-order off-bundle conservative-even defect is carried by exactly one slippage scalar \(\varepsilon_\perp\), while the odd quadrupole defect still needs the mixed-port/output coefficient \(\delta\gamma_W\).

---

## 7. First-order preservation theorem in the new variables

Requiring preservation of the canonical conservative even \(l=2\) fingerprint means
\[
\delta E_2=0,
\qquad
\delta E_4=0.
\]

Using the formulas above, the same determinant argument as in Stage 159 gives
\[
\boxed{
\varepsilon_\perp=0,
\qquad
\delta\kappa_W=0.
}
\]

So the even-preservation theorem now has a microscopic off-bundle reading:

> exact first-order preservation of the canonical conservative even branch kills the weighted scalar slippage combination \(\varepsilon_\perp\) and the mixed-lane scale defect \(\delta\kappa_W\).

When that happens, the remaining quadrupole-normalization defect again collapses to
\[
\boxed{
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

So the Stage 159 result survives intact, but we now know exactly what off-bundle microscopic scalar defect must vanish for it to hold.

---

## 8. Interpretation

This stage sharpens the Stage 167 conclusion in the cleanest possible way.

Stage 167 said:

- first-order isotropic bundle drift is tangent to the exact lower compensated branch,
- so the first genuine defect must come from off-bundle structure.

Stage 168 shows that, at first order, this off-bundle structure is not arbitrary.
It collapses to three slippages
\[
(\varepsilon_L,\varepsilon_v,\varepsilon_T),
\]
and those three slippages only matter through the one scalar combination
\[
\varepsilon_\perp.
\]

So the live first-order bottleneck is now:

1. derive the grouped-lane anisotropy defects separately,
2. derive the three scalar transport-law slippages \((\varepsilon_L,\varepsilon_v,\varepsilon_T)\),
3. compute the one weighted combination \(\varepsilon_\perp\).

That is much smaller than “solve all first-order corrections.”

---

## 9. Immediate next theorem gate

The next clean pass is now obvious.

Do **not** reopen generic isotropic bundle transport; that theorem is finished.

Instead derive the microscopic source of the three slippages:
\[
\varepsilon_L,
\qquad
\varepsilon_v,
\qquad
\varepsilon_T,
\]
for example from:

- grouped-lane anisotropy and mouth/worldtube splitting,
- breakdown of the exact D/N geometric similarity law,
- departure from the lower-branch traction/background transport identities,
- or the first nonlinear correction to the explicit Family-1 co-evolving branch.

That is the smallest next derivation step that can actually move the theorem forward.
