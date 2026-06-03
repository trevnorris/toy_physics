# Moving-Throat PDE — Stage 164: Explicit Microscopic Log-Imbalance Channels on the Linearized Wall/BdG/Maxwell/Mixed Branch

## Goal

Stage 163 reduced the first off-family defect of the co-evolving Family-1 branch to the single scalar
\[
\delta_\perp
=
\mathfrak g_*\,
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}\,
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right).
\]

So the next real question is no longer broad:

> what do the two logarithmic imbalance channels
> \[
> \delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right),
> \qquad
> \delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
> \]
> actually become on the explicit linearized moving-throat wall/BdG/Maxwell/mixed branch?

This note computes them.

The main result is that, inside the explicit throat-core closure already fixed at Stages 118–119,
the two channels collapse to concrete logarithmic drifts of the microscopic variables
\[
(a,L_W,\rho_w,c_{s,w},c_s,\mathcal Z_q,\mathcal T_m,v_{w0}),
\]
and the full off-family scalar becomes an explicit linear co-transport law on that branch.

So the remaining PDE-facing task is now much narrower than Stage 163 made it sound:
it is no longer “compute the channels somehow,” but only to determine the actual linearized
values of these specific branch-variable drifts on the true moving-throat solution.

---

## 1. Two exact identities before any further closure

The Stage 119 parent ratios are
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\qquad
r_c=\frac{\lambda^2}{K_sK_q}=\mathfrak r^2.
\]

So the two Stage 163 logarithmic channels already have an exact microscopic meaning:
\[
\boxed{
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\!\left(\frac{\mathfrak g}{\mathfrak r}\right),
}
\]
\[
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
-\delta\ln r_c.
}
\]

This is a very strong simplification.

- The first channel is the logarithmic drift of the **relative mouth-coupling ratio**
  \(\mathfrak g/\mathfrak r\).
- The second channel is exactly the negative logarithmic drift of the
  **static/mixed hybridization ratio** \(r_c\).

So Stage 163’s microscopic normal coordinate can already be rewritten as
\[
\boxed{
\delta_\perp
=
\mathfrak g_*\,
\delta\ln\!\left(\frac{\mathfrak g}{\mathfrak r}\right)
-
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}\,
\delta\ln r_c.
}
\]

---

## 2. Exact parent-action formulas on the explicit throat-core branch

Stage 118 fixed
\[
g_s=\mathcal T_m J_s,
\qquad
J_s=\frac{4\pi a^2\ell}{3},
\]
\[
g_q=\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi}{\sqrt2\,L_W^{3/2}},
\qquad
K_q=\frac{\mathcal Z_q}{\mu_0}\,\frac{\pi^2 c_s^2}{4L_W^2},
\]
\[
\lambda=-q_*v_{w0}\,\mathcal I_{sq}.
\]

So before using any overlap simplification, the two logarithmic channels are exactly
\[
\boxed{
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln K_s
+\delta\ln\mathcal Z_q
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-\delta\ln J_s
-\delta\ln\mathcal I_{sq}
-\frac32\,\delta\ln L_W,
}
\]
\[
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln K_s
+\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
-2\,\delta\ln v_{w0}
-2\,\delta\ln\mathcal I_{sq}
-2\,\delta\ln L_W.
}
\]

This is already an honest moving-throat branch statement.

At this level the two channels are controlled by four microscopic sectors:

1. the scalar/wall shell stiffness \(K_s\),
2. the mixed-channel localization norm \(\mathcal Z_q\),
3. the mouth-source amplitude \(\mathcal T_m\),
4. and the background mixed overlap \(v_{w0}\mathcal I_{sq}\).

---

## 3. Uniform-overlap + D/N simplification

Now impose the same explicit overlap closure already used at Stage 118:
\[
\mathcal I_{sq}=J_s I_q,
\qquad
J_s=\frac{4\pi a^2\ell}{3},
\qquad
I_q=\frac{2\sqrt{2L_W}}{\pi}.
\]

Then
\[
\mathcal I_{sq}\propto a^2\ell\,L_W^{1/2},
\]
so
\[
\delta\ln J_s = 2\,\delta\ln a + \delta\ln\ell,
\qquad
\delta\ln\mathcal I_{sq}
=
2\,\delta\ln a+\delta\ln\ell+\frac12\,\delta\ln L_W.
\]

Substituting into the exact general formulas gives
\[
\boxed{
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln K_s
+\delta\ln\mathcal Z_q
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-4\,\delta\ln a
-2\,\delta\ln\ell
-2\,\delta\ln L_W,
}
\]
\[
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln K_s
+\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
-2\,\delta\ln v_{w0}
-4\,\delta\ln a
-2\,\delta\ln\ell
-3\,\delta\ln L_W.
}
\]

So under the first explicit throat-core closure, the Stage 163 channels are already
concrete logarithmic drifts of the actual branch variables
\[
(K_s,a,\ell,L_W,\mathcal Z_q,\mathcal T_m,v_{w0},c_s).
\]

---

## 4. Healing-locked shell simplification

Stage 118 also recorded the carried healing-lock shell branch
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\qquad
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
\]

So on that branch
\[
K_s
=
\frac{6\pi\hbar}{5}\,
\frac{a^2 c_{s,w}}{\rho_w},
\]
and the two products themselves become
\[
\boxed{
\frac{g_qK_s}{g_s\lambda}
=
-\frac{27\pi m_\psi^2}{40\hbar\mu_0 q_*}\,
\frac{\mathcal Z_q\,c_{s,w}^3}
{\rho_w\,\mathcal T_m\,v_{w0}\,a^2L_W^2},
}
\]
\[
\boxed{
\frac{K_sK_q}{\lambda^2}
=
\frac{27\pi^3 m_\psi^2}{320\hbar\mu_0 q_*^2}\,
\frac{\mathcal Z_q\,c_s^2 c_{s,w}^3}
{\rho_w\,v_{w0}^2\,a^2L_W^3}.
}
\]

Therefore the Stage 163 logarithmic channels become
\[
\boxed{
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
}
\]
\[
\boxed{
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
=
\delta\ln\mathcal Z_q
+2\,\delta\ln c_s
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-2\,\delta\ln v_{w0}
-2\,\delta\ln a
-3\,\delta\ln L_W.
}
\]

This is the first point where the two Stage 163 channels have become a literal
linearized moving-throat wall/BdG/Maxwell/mixed-branch formula.

Interpretation:

- \((a,L_W)\) are the wall-geometry coordinates,
- \((\rho_w,c_{s,w})\) are the scalar/BdG shell-background variables,
- \(\mathcal Z_q\) is the localized-Maxwell mixed-channel norm,
- \(v_{w0}\) is the static mixed background flow,
- \(\mathcal T_m\) is the mouth traction amplitude,
- \(c_s\) is the mixed-tube propagation speed.

---

## 5. Explicit off-family scalar on the actual branch variables

Now insert these Stage 164 formulas into the exact Stage 163 identity
\[
\delta_\perp
=
\mathfrak g_*\,
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}\,
\delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right).
\]

Define
\[
A_*:=\mathfrak g_*+\frac{1}{4\sqrt{1+\mathfrak r_*^2}},
\qquad
B_*:=\frac{1}{2\sqrt{1+\mathfrak r_*^2}}.
\]
Then the off-family scalar becomes
\[
\boxed{
\delta_\perp
=
A_*\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+3A_*\,\delta\ln c_{s,w}
+B_*\,\delta\ln c_s
-\mathfrak g_*\,\delta\ln\mathcal T_m
-\left(\mathfrak g_*+B_*\right)\delta\ln v_{w0}
-2A_*\,\delta\ln a
-\left(2\mathfrak g_*+\frac{3}{4\sqrt{1+\mathfrak r_*^2}}\right)\delta\ln L_W.
}
\]

On the renormalized Family-1 canonical point,
\[
\mathfrak g_* \approx 0.758035078944663,
\qquad
\mathfrak r_* \approx 1.77799353547498,
\qquad
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}
\approx 0.122554011096907,
\]
so
\[
A_*\approx 0.880589090041570,
\qquad
B_*\approx 0.245108022193814.
\]

Therefore
\[
\boxed{
\delta_\perp
\approx
0.880589090041570\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+
2.64176727012471\,\delta\ln c_{s,w}
+
0.245108022193814\,\delta\ln c_s
-
0.758035078944663\,\delta\ln\mathcal T_m
}
\]
\[
\boxed{
\qquad
-
1.00314310113848\,\delta\ln v_{w0}
-
1.76117818008314\,\delta\ln a
-
1.88373219118005\,\delta\ln L_W.
}
\]

This is the explicit microscopic normal-coordinate formula that Stage 163 was asking for.

So the first true off-family defect is now a weighted sum of exactly seven logarithmic branch drifts.

---

## 6. Exact tangency law on the explicit branch

Setting \(\delta_\perp=0\) gives the exact first-order condition for the true moving-throat
branch to stay tangent to the lower compensated parent family.

Solving for the mouth-traction drift gives
\[
\boxed{
\delta\ln\mathcal T_m
=
\left(1+\frac{1}{4\mathfrak g_*\sqrt{1+\mathfrak r_*^2}}\right)
\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+
3\left(1+\frac{1}{4\mathfrak g_*\sqrt{1+\mathfrak r_*^2}}\right)\delta\ln c_{s,w}
}
\]
\[
\boxed{
\qquad
+
\frac{1}{2\mathfrak g_*\sqrt{1+\mathfrak r_*^2}}\,\delta\ln c_s
-
\left(1+\frac{1}{2\mathfrak g_*\sqrt{1+\mathfrak r_*^2}}\right)\delta\ln v_{w0}
-
2\left(1+\frac{1}{4\mathfrak g_*\sqrt{1+\mathfrak r_*^2}}\right)\delta\ln a
}
\]
\[
\boxed{
\qquad
-
\left(
2+\frac{3}{4\mathfrak g_*\sqrt{1+\mathfrak r_*^2}}
\right)\delta\ln L_W.
}
\]

Numerically,
\[
\boxed{
\delta\ln\mathcal T_m
\approx
1.16167271174229\,\delta\ln\!\left(\frac{\mathcal Z_q}{\rho_w}\right)
+
3.48501813522686\,\delta\ln c_{s,w}
+
0.323345423484574\,\delta\ln c_s
}
\]
\[
\boxed{
\qquad
-
1.32334542348457\,\delta\ln v_{w0}
-
2.32334542348457\,\delta\ln a
-
2.48499999999999\,\delta\ln L_W.
}
\]

So the “stay on the parent family” condition is now an explicit co-transport law for the actual microscopic branch variables.

---

## 7. Best current theorem statement after Stage 164

The two Stage 163 microscopic logarithmic imbalance channels are no longer abstract.

Inside the explicit first throat-core closure:

1. one channel is exactly the logarithmic drift of the relative mouth-coupling ratio,
   \[
   \delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
   =
   \delta\ln\!\left(\frac{\mathfrak g}{\mathfrak r}\right),
   \]
2. the other is exactly the negative logarithmic drift of the hybridization ratio,
   \[
   \delta\ln\!\left(\frac{K_sK_q}{\lambda^2}\right)
   =-\delta\ln r_c,
   \]
3. and on the carried healing-locked wall/BdG/Maxwell/mixed branch they reduce to the explicit microscopic drifts of
   \[
   (\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W).
   \]

So the remaining PDE-facing task is no longer to discover what the channels are.
It is only to determine the actual linearized values of these branch-variable drifts on the true moving-throat solution.

If they satisfy the tangency law above, then \(\delta_\perp=0\) and the first-order reduced
\(2.5\)PN obstruction disappears automatically.
If they do not, the whole first-order off-family defect is now already explicit.
