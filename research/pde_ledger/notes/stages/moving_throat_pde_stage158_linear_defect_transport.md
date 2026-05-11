# Moving-Throat PDE — Stage 158: Linear Defect Transport from the Renormalized Family-1 Canonical Point

## Goal

Stage 242 left the program at a very narrow finish line:

> derive the actual deviation of the moving-throat mouth/core system from the
> explicit co-evolving Family-1 canonical fixed point, and translate that
> deviation into the remaining outgoing quadrupole-normalization defect.

This note does the first half of that job exactly.

It does **not** yet solve the full DtN projection into the outgoing branch data.
What it does is reduce the mouth/core side to a small linear defect ledger around
`(Sigma0_can, g_*, S_can)`.

The remaining unsolved step is then only the projection

`(delta Sigma0, delta g, delta S) -> (b, a0, a5)`

into the isotropic outgoing-branch deformation triple of Stages 90–92.

---

## 1. Renormalized canonical base point

From Stage 241, the exact co-evolving compensated Family-1 point is

\[
\Sigma_0^{\rm can}
\approx 4.651033550168876,
\qquad
\widehat T_{m,\rm can}
\approx 1.4467083664567624,
\]
\[
\mathfrak g_{\rm can}=\mathfrak g_*
\approx 0.758035078944663,
\qquad
\mathcal R_{\rm can}=\frac14,
\qquad
\mathcal S_{\rm can}
\approx 0.6703621156734617,
\]
\[
\Pi_{\rm can}
\approx 3.8715643774790087.
\]

The carried Family-1 core parameter is

\[
\mathfrak r_{F1}\approx 1.77799353547498,
\qquad
\sqrt{1+\mathfrak r_{F1}^2}
\approx 2.039916913060632.
\]

---

## 2. Exact first transport step: `delta g -> delta R`

Stage 239 already fixed the exact Family-1 core law

\[
\mathcal R(\mathfrak g)
=
\frac{(\mathfrak g-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\]

and on the lower compensated branch

\[
\mathfrak g=\mathfrak g_*+\delta\mathfrak g
\qquad\Longrightarrow\qquad
\mathcal R
=
\frac14
-
\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+
\frac{(\delta\mathfrak g)^2}{1+\mathfrak r_{F1}^2}.
\]

So the exact linear defect transport is

\[
\boxed{
\delta\mathcal R
=
-
\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+O(\delta\mathfrak g^2)
}
\]

and numerically

\[
\boxed{
\delta\mathcal R
\approx
-0.490216044387626\,\delta\mathfrak g.
}
\]

So broadening of the mouth source (`delta g < 0`) still pushes the shell/mixed
ratio upward (`delta R > 0`).

---

## 3. Mouth-gain transport at the co-evolving canonical point

From Stages 188–189, the actual coupled mouth gains are

\[
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\,\mathcal R.
\]

Therefore around the co-evolving canonical point,

\[
\delta M_s=\delta\Sigma_0,
\]
\[
\delta M_q
=
-\mathcal R_*\,\delta\Sigma_0
-\Sigma_{0,*}\,\delta\mathcal R
+O(2),
\qquad
\mathcal R_* = \frac14.
\]

Using the exact `delta R` transport above,

\[
\boxed{
\delta M_q
=
-
\frac14\,\delta\Sigma_0
+
\frac{\Sigma_0^{\rm can}}{\sqrt{1+\mathfrak r_{F1}^2}}\,\delta\mathfrak g
+O(2)
}
\]

with numerical coefficient

\[
\boxed{
\delta M_q
\approx
-0.25\,\delta\Sigma_0
+
2.28001126927792\,\delta\mathfrak g.
}
\]

Because

\[
\Sigma_0 = \frac{20}{9}\widehat T_m^2,
\]

one also has the traction-level transport

\[
\boxed{
\delta\Sigma_0
=
\frac{40}{9}\widehat T_{m,\rm can}\,\delta\widehat T_m
+O(\delta\widehat T_m^2)
\approx
6.42981496203006\,\delta\widehat T_m.
}
\]

So in the traction variable,

\[
\delta M_q
\approx
-1.60745374050751\,\delta\widehat T_m
+
2.28001126927792\,\delta\mathfrak g.
\]

---

## 4. Slope / mouth-bias transport

The exact slope identity on the coupled branch is

\[
\Pi = \Sigma_0\bigl(1-\mathcal R\mathcal S\bigr)=M_s+M_q\mathcal S.
\]

Linearizing about the renormalized canonical point gives

\[
\delta\Pi
=
(1-\mathcal R_*\mathcal S_*)\,\delta\Sigma_0
-
\Sigma_{0,*}\bigl(\mathcal R_*\,\delta\mathcal S + \mathcal S_*\,\delta\mathcal R\bigr)
+O(2).
\]

Substituting `R_* = 1/4`, `S_* = S_can`, and the exact `delta R` transport,

\[
\boxed{
\delta\Pi
=
\left(1-\frac14\mathcal S_{\rm can}\right)\delta\Sigma_0
-
\frac{\Sigma_0^{\rm can}}{4}\,\delta\mathcal S
+
\frac{\Sigma_0^{\rm can}\mathcal S_{\rm can}}{\sqrt{1+\mathfrak r_{F1}^2}}\,\delta\mathfrak g
+O(2)
}
\]

or numerically

\[
\boxed{
\delta\Pi
\approx
0.832409471081635\,\delta\Sigma_0
-
1.16275838754222\,\delta\mathcal S
+
1.52843317823248\,\delta\mathfrak g.
}
\]

In terms of the traction deviation,

\[
\boxed{
\delta\Pi
\approx
5.35223887169622\,\delta\widehat T_m
-
1.16275838754222\,\delta\mathcal S
+
1.52843317823248\,\delta\mathfrak g.
}
\]

So the co-evolving canonical point is now surrounded by a fully explicit linear
mouth/core defect map.

---

## 5. Reduced 2.5PN bridge on the natural source-map branch

Stages 83–84 already reduced the remaining point-particle 2.5PN obstruction to the
outgoing-normalization factor `chi_Q`.

On the natural source-map branch,

\[
N_Q = \frac{1}{\chi_Q},
\qquad
\Delta_Q := \chi_Q-1,
\qquad
N_Q-1 = -\Delta_Q + O(\Delta_Q^2).
\]

Stages 90–92 then showed that the most general **small isotropic outgoing-branch**
deformation is controlled by the triple

\[
(b,a_0,a_5),
\]

with linearized normalization shift

\[
\boxed{
\Delta_Q
=
5b + \frac{a_0}{3} + 9a_5 + O(2).
}
\]

So the whole reduced theorem gap after Stage 242 is now exactly the Jacobian

\[
\boxed{
(\delta\Sigma_0,\delta\mathfrak g,\delta\mathcal S)
\longrightarrow
(b,a_0,a_5).
}
\]

Everything on the mouth/core side upstream of that projection is already explicit.

---

## 6. Three immediate projection models

Even before the full DtN projection is solved, the earlier outlet audit already
provides three useful reduced subcases.

### A. Pure argument deformation

If the mouth/core defect appears only as an effective outgoing-argument shift,

\[
a_0=a_5=0,
\]

then

\[
\boxed{\Delta_Q \approx 5b.}
\]

### B. Pure Robin/static core deformation

If the defect is dominantly a static isotropic Robin outlet,

\[
b=a_5=0,
\qquad
a_0=\rho_R,
\]

then

\[
\boxed{\Delta_Q \approx \frac{\rho_R}{3}.}
\]

### C. Even-preserving compensated Robin–mixed branch

On the exact compensated hybrid branch of Stage 95,

\[
\rho_R=4\sigma_W,
\qquad
\kappa_W=\frac13,
\qquad
\chi_Q^{\rm hyb}=
\frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]

So if the actual branch stays close to the canonical outgoing value

\[
\gamma_W = \frac19 + \delta\gamma_W,
\]

then

\[
\boxed{
\Delta_Q
=
-\frac{9\sigma_W}{1-\sigma_W}\,\delta\gamma_W
+O(\delta\gamma_W^2).
}
\]

In this subcase the even branch stays canonical and the whole residual defect is
purely the odd mixed-channel renormalization.

---

## 7. Best current reduced statement after Stage 242

The program no longer needs to discover a broad new mouth/core structure.
That structure is already explicit.

What remains is much narrower:

1. the exact co-evolving Family-1 canonical point is known,
2. the linear defect transport
   
   `delta Sigma0, delta g, delta S -> delta M_s, delta M_q, delta Pi`
   
   is known,
3. the remaining reduced 2.5PN scalar is
   
   `Delta_Q = 5 b + a0/3 + 9 a5 + O(2)`,
4. so the only unsolved bridge is the DtN projection
   
   `(delta M_s, delta M_q, delta Pi) -> (b,a0,a5)`.

That is the exact next theorem gate.

## Immediate next step

The clean next derivation is now:

1. choose the first DtN projection class to test,
2. express `(b,a0,a5)` in terms of the actual mouth/core defect variables,
3. substitute the transport formulas above,
4. and evaluate the resulting `Delta_Q`.

The most economical first pass is to test the **compensated Robin–mixed branch**,
because it is the only explicit isotropic outlet class already known to preserve the
canonical even `l=2` fingerprint while allowing a nontrivial odd renormalization.
