
# Moving-Throat PDE — Stage 142: Linear Projection of the Co-Evolving Family-1 Defect onto the Compensated Robin–Mixed Outlet

## Goal

Stage 141 reduced the explicit co-evolving Family-1 mouth/core system to the linear
defect ledger
\[
(\delta\Sigma_0,\delta\mathfrak g,\delta\mathcal S)
\longrightarrow
(\delta M_s,\delta M_q,\delta\Pi).
\]

The clean next pass was to test the only explicit isotropic outlet class already known
to preserve the canonical even `l=2` fingerprint while still allowing a nontrivial odd
renormalization: the **compensated Robin–mixed outlet** of Stages 95–98.

This note carries out that projection.

The main result is that, on the compensated hybrid outlet:

1. the direct mouth/core departure from the compensated outlet surface is measured by
   one scalar only,
   \[
   \delta\mathcal C := \delta\rho_R - 4\,\delta\sigma_W,
   \]
2. in the explicit Family-1 closure that scalar depends **only** on the source-moment
   defect `\delta\mathfrak g`, not on `\delta\Sigma_0` or `\delta\mathcal S`,
3. exact preservation of the canonical even `l=2` fingerprint at first order forces
   \[
   \delta\mathfrak g = 0,
   \qquad
   \delta\kappa_W = 0,
   \]
4. and then the last reduced 2.5PN defect collapses to a single odd quantity:
   \[
   \Delta_Q = \chi_Q-1
   =
   -\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W
   +O(2).
   \]

So the next serious PDE-facing task is no longer a broad
`(\delta\Sigma_0,\delta\mathfrak g,\delta\mathcal S)\to(b,a_0,a_5)` projection.
Inside the explicit compensated hybrid class, the even branch already tells us that
the only surviving linear obstruction is the odd mixed-channel renormalization
`\delta\gamma_W`.

---

## 1. Core-outlet loadings from the actual mouth gains

Stages 120–121 already fixed the exact relation between the isotropic hybrid outlet
loads and the mouth gains:
\[
M_s=\frac{L}{\Theta_\sigma}\,\rho_R,
\qquad
M_q=-\frac{L}{\Theta_\sigma}\,\sigma_W.
\]

So with
\[
\Xi:=\frac{\Theta_\sigma}{L},
\]
one has
\[
\boxed{
\rho_R=\Xi M_s,
\qquad
\sigma_W=-\Xi M_q.
}
\]

At the renormalized canonical Family-1 point of Stage 139,
\[
M_{s,*}=\Sigma_0^{\rm can},
\qquad
M_{q,*}=-\frac{\Sigma_0^{\rm can}}{4},
\]
so the compensated hybrid branch parameter is
\[
\boxed{
\sigma_*=\frac{\Xi\Sigma_0^{\rm can}}{4},
\qquad
\rho_{R,*}=4\sigma_*.
}
\]

---

## 2. Direct mouth/core transport into the hybrid loading defect

Stage 141 gave
\[
\delta M_s=\delta\Sigma_0,
\qquad
\delta M_q
=
-\frac14\,\delta\Sigma_0
-\Sigma_0^{\rm can}\,\delta\mathcal R,
\]
with
\[
\delta\mathcal R
=
-\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+O(\delta\mathfrak g^2),
\qquad
\mathfrak r_{F1}\approx 1.77799353547498.
\]

Therefore
\[
\delta\rho_R=\Xi\,\delta\Sigma_0,
\]
\[
\delta\sigma_W
=
-\Xi\,\delta M_q
=
\Xi\left(
\frac14\,\delta\Sigma_0
+
\Sigma_0^{\rm can}\,\delta\mathcal R
\right).
\]

So the direct departure from the compensated outlet surface is
\[
\boxed{
\delta\mathcal C
:=
\delta\rho_R-4\,\delta\sigma_W
=
-4\Xi\,\Sigma_0^{\rm can}\,\delta\mathcal R.
}
\]

Using
\[
\sigma_*=\frac{\Xi\Sigma_0^{\rm can}}{4},
\]
this becomes
\[
\boxed{
\delta\mathcal C
=
-16\sigma_*\,\delta\mathcal R.
}
\]

And substituting the Family-1 transport law,
\[
\boxed{
\delta\mathcal C
=
\frac{16\sigma_*}{\sqrt{1+\mathfrak r_{F1}^2}}\,
\delta\mathfrak g
+O(\delta\mathfrak g^2).
}
\]

Numerically,
\[
\sqrt{1+\mathfrak r_{F1}^2}\approx 2.039916913060632,
\]
so
\[
\boxed{
\delta\mathcal C
\approx
7.84345671020202\,\sigma_*\,\delta\mathfrak g.
}
\]

The important structural point is that the `\delta\Sigma_0` pieces cancel exactly.
So traction changes are tangent to the compensated hybrid family at first order,
whereas `\delta\mathfrak g` is the genuine normal/off-compensation defect.

---

## 3. Exact linear outlet algebra around the compensated hybrid branch

Take the full hybrid outlet
\[
\Lambda_2^{\rm hyb}(z)
=
\Lambda_2^{\rm out}(z)
+\rho_R
-
\frac{\sigma_W}{1-\kappa_W z^2-i\gamma_W z^5},
\]
and linearize around the compensated canonical branch
\[
\rho_{R,*}=4\sigma_*,
\qquad
\kappa_{W,*}=\frac13,
\qquad
\gamma_{W,*}=\frac19.
\]

Define the even-fingerprint defects
\[
\delta E_2
:=
-\frac{L_2}{L_0}-\frac19,
\qquad
\delta E_4
:=
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}-\frac{4}{81},
\]
and the outgoing-normalization defect
\[
\Delta_Q:=\chi_Q-1.
\]

Then the exact first-order formulas are
\[
\boxed{
\delta E_2
=
\frac{\delta\mathcal C-9\sigma_*\,\delta\kappa_W}
{27(1-\sigma_*)},
}
\]
\[
\boxed{
\delta E_4
=
\frac{5\,\delta\mathcal C-72\sigma_*\,\delta\kappa_W}
{243(1-\sigma_*)},
}
\]
\[
\boxed{
\Delta_Q
=
\frac{\delta\mathcal C-27\sigma_*\,\delta\gamma_W}
{3(1-\sigma_*)}.
}
\]

So the compensated hybrid outlet separates the remaining linear data into:

- one **even/off-compensation** scalar `\delta\mathcal C`,
- one **even dispersion** scalar `\delta\kappa_W`,
- one **odd** scalar `\delta\gamma_W`.

---

## 4. Substituting the explicit Family-1 mouth/core transport

Using the Stage-141/Stage-142 transport
\[
\delta\mathcal C=-16\sigma_*\,\delta\mathcal R
=
\frac{16\sigma_*}{\sqrt{1+\mathfrak r_{F1}^2}}\,
\delta\mathfrak g,
\]
the outlet defects become
\[
\boxed{
\delta E_2
=
-\frac{16\sigma_*}{27(1-\sigma_*)}\,\delta\mathcal R
-
\frac{\sigma_*}{3(1-\sigma_*)}\,\delta\kappa_W,
}
\]
\[
\boxed{
\delta E_4
=
-\frac{80\sigma_*}{243(1-\sigma_*)}\,\delta\mathcal R
-
\frac{8\sigma_*}{27(1-\sigma_*)}\,\delta\kappa_W,
}
\]
\[
\boxed{
\Delta_Q
=
-\frac{16\sigma_*}{3(1-\sigma_*)}\,\delta\mathcal R
-
\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

In terms of `\delta\mathfrak g`,
\[
\boxed{
\Delta_Q
=
\frac{16\sigma_*}{3\sqrt{1+\mathfrak r_{F1}^2}\,(1-\sigma_*)}\,
\delta\mathfrak g
-
\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

Numerically,
\[
\boxed{
\Delta_Q
\approx
\frac{2.61448557006734\,\sigma_*}{1-\sigma_*}\,\delta\mathfrak g
-
\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
}
\]

So if the outlet is not protected to stay on the canonical-even branch, the
direct mouth/core broadening `\delta\mathfrak g` feeds linearly into the final
outgoing-normalization defect.

---

## 5. Exact linear canonical-even gate

Now impose the actual condition we want:

> preserve the canonical conservative even `l=2` fingerprint at first order.

That means
\[
\delta E_2=0,
\qquad
\delta E_4=0.
\]

Using the formulas above, this is the linear system
\[
\delta\mathcal C-9\sigma_*\,\delta\kappa_W=0,
\]
\[
5\delta\mathcal C-72\sigma_*\,\delta\kappa_W=0.
\]

Its determinant is
\[
-72\sigma_*+45\sigma_*=-27\sigma_*,
\]
so on a nontrivial loaded branch `\sigma_*\neq 0` the unique solution is
\[
\boxed{
\delta\mathcal C=0,
\qquad
\delta\kappa_W=0.
}
\]

Because
\[
\delta\mathcal C
=
-16\sigma_*\,\delta\mathcal R
=
\frac{16\sigma_*}{\sqrt{1+\mathfrak r_{F1}^2}}\,
\delta\mathfrak g,
\]
this immediately implies
\[
\boxed{
\delta\mathcal R=0,
\qquad
\delta\mathfrak g=0.
}
\]

So the compensated hybrid outlet does something very strong:

> exact first-order preservation of the canonical even `l=2` branch forces the
> explicit Family-1 mouth/core system to stay on the compensated source moment
> `\mathfrak g=\mathfrak g_*`, and it simultaneously forbids a first-order shift
> in the hidden even pole parameter `\kappa_W`.

This is the sharpened bridge that Stage 141 was pointing toward.

---

## 6. Collapse of the last linear 2.5PN defect

Once the canonical-even gate is imposed,
\[
\delta\mathfrak g=0,
\qquad
\delta\kappa_W=0,
\]
so the Stage-142 defect law collapses to
\[
\boxed{
\Delta_Q
=
-\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W
+O(2).
}
\]

On the natural source-map branch,
\[
N_Q-1=-\Delta_Q+O(2),
\]
so
\[
\boxed{
N_Q-1
=
\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W
+O(2).
}
\]

This is the cleanest reduced statement reached so far.

Inside the explicit compensated Robin–mixed outlet class:

- the linear mouth/core off-compensation defect is exactly `\delta\mathfrak g`,
- canonical-even preservation kills it,
- traction shifts `\delta\Sigma_0` are tangent and drop out,
- the direct mouth-layer transport coefficient `\delta\mathcal S` also drops out of the
  outlet loading defect,
- and the only surviving first-order 2.5PN obstruction is the odd mixed-channel
  renormalization `\delta\gamma_W`.

---

## 7. Immediate PDE-facing consequence

Stage 141 ended with the unresolved bridge
\[
(\delta M_s,\delta M_q,\delta\Pi)\to(b,a_0,a_5).
\]

Inside the compensated hybrid outlet class, Stage 142 sharpens that to:

1. the even fingerprint already fixes the direct loading projection:
   \[
   \delta\mathfrak g=0,
   \qquad
   \delta\kappa_W=0,
   \]
2. so the only remaining live projection is the **odd tangential one**
   \[
   (\delta\Sigma_0,\delta\mathcal S)\longrightarrow \delta\Pi_{\rm tan}
   \longrightarrow \delta\gamma_W,
   \]
   with
   \[
   \delta\Pi_{\rm tan}
   =
   0.832409471081635\,\delta\Sigma_0
   -
   1.16275838754222\,\delta\mathcal S
   \qquad
   (\delta\mathfrak g=0).
   \]
3. The last reduced 2.5PN scalar is then
   \[
   \boxed{
   \Delta_Q
   =
   -\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
   }
   \]

So the next serious calculation is no longer a broad three-variable mouth/core
projection. It is the much narrower DtN problem:

> compute the odd mixed-channel renormalization `\delta\gamma_W` induced by the
> tangential co-evolving mouth deformation on the compensated hybrid outlet.

That is now the exact next theorem gate.
