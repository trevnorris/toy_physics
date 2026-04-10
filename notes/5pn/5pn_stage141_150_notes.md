
# 5PN / Moving-Throat continuation — Stages 141–150

This batch codifies the next handwritten moving-throat block after Stage 140: the linear defect-transport, hybrid-outlet projection, bare mixed-port slippage, D/N similarity decomposition, parent-family rigidity, microscopic off-family normal coordinate, explicit log-channel reduction, exact lower-branch drift laws, four-observable bundle inversion, and the tangent-compensation theorem.

## What is newly frozen in executable form

### Stage 141
The co-evolving Family-1 point is reduced to the linear mouth/core defect ledger
\[
\delta\mathcal R = -\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_*^2}},
\]
\[
\delta M_s = \delta\Sigma_0,
\qquad
\delta M_q = -\frac14\,\delta\Sigma_0 + \frac{\Sigma_{0,*}}{\sqrt{1+\mathfrak r_*^2}}\,\delta\mathfrak g,
\]
\[
\delta\Pi
=
\left(1-\frac14\mathcal S_*\right)\delta\Sigma_0
-\frac{\Sigma_{0,*}}{4}\,\delta\mathcal S
+\frac{\Sigma_{0,*}\mathcal S_*}{\sqrt{1+\mathfrak r_*^2}}\delta\mathfrak g.
\]

### Stage 142
Projecting the defect to the compensated Robin–mixed outlet gives the exact first-order outlet defects
\[
\delta E_2 = \frac{\delta\mathcal C - 9\sigma_*\,\delta\kappa_W}{27(1-\sigma_*)},
\]
\[
\delta E_4 = \frac{5\delta\mathcal C - 72\sigma_*\,\delta\kappa_W}{243(1-\sigma_*)},
\]
\[
\Delta_Q = \frac{\delta\mathcal C - 27\sigma_*\,\delta\gamma_W}{3(1-\sigma_*)}.
\]
Imposing the canonical-even gate \(\delta E_2=\delta E_4=0\) yields
\[
\delta\mathcal C = 0,
\qquad
\delta\kappa_W = 0,
\qquad
\delta\mathfrak g = 0,
\]
and therefore
\[
\Delta_Q = -\frac{9\sigma_*}{1-\sigma_*}\,\delta\gamma_W.
\]

### Stage 143
The concrete core algebra collapses the remaining odd defect to the bare mixed-port slippage scalar
\[
\delta\mathfrak B_W := \delta\gamma_0 - \frac13\,\delta\kappa_0,
\qquad
\delta\gamma_W = \frac{\delta\mathfrak B_W}{1+r_{c,*}}.
\]
With the tangential susceptibility \(\delta\mathfrak B_W = \Upsilon_\Pi\,\delta\Pi_{\rm tan}\),
\[
\Delta_Q
=
-\frac{9\sigma_*\,\Upsilon_\Pi}{(1-\sigma_*)(1+r_{c,*})}\,\delta\Pi_{\rm tan}.
\]

### Stage 144
The black-box susceptibility is decomposed into the D/N similarity-slippage scalar
\[
\Xi_{\rm slip}:=\Xi_\gamma - 2\Xi_L,
\qquad
\Upsilon_\Pi = \frac{1+r_{c,*}}{9}\,\Xi_{\rm slip},
\]
so the reduced defect becomes
\[
\Delta_Q = -\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]
If the D/N similarity law
\[
\gamma_0 = \frac{4L_W^2}{3\pi^2 a^2}
\]
is preserved to first order, then \(\Xi_{\rm slip}=0\).

### Stage 145
On the exact parent compensation family
\[
1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2,
\qquad
\frac{L_W}{a} = \frac{\pi}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
\qquad
\gamma_0 = \frac{1+\mathfrak r^2}{9},
\]
automatic similarity preservation is exact:
\[
\delta\ln\gamma_0 - 2\,\delta\ln(L_W/a) = 0.
\]
On the lower branch, \(\delta\mathfrak g=0\) implies \(\delta\mathfrak r=0\), so every first-order similarity defect vanishes and
\[
\Delta_Q = 0.
\]

### Stage 146
The exact off-family normal coordinate is
\[
\delta_\perp := \delta\mathfrak g - \mathfrak g'_-(\mathfrak r_*)\,\delta\mathfrak r,
\]
with
\[
\delta\mathcal F = 4\sqrt{1+\mathfrak r_*^2}\,\delta_\perp,
\qquad
\delta R_q = -\frac{\delta_\perp}{\sqrt{1+\mathfrak r_*^2}}.
\]
Its explicit microscopic form is
\[
\delta_\perp
=
\mathfrak g_*\,
\delta\ln\!\left(\frac{g_q K_s}{g_s\lambda}\right)
+
\frac{1}{4\sqrt{1+\mathfrak r_*^2}}\,
\delta\ln\!\left(\frac{K_s K_q}{\lambda^2}\right).
\]

### Stage 147
Those two log channels are reduced to explicit wall/BdG/Maxwell/mixed drifts. Under the carried overlap and healing-lock closures:
\[
\delta\ln\!\left(\frac{g_qK_s}{g_s\lambda}\right)
=
\delta\ln\mathcal Z_q
+3\,\delta\ln c_{s,w}
-\delta\ln\rho_w
-\delta\ln\mathcal T_m
-\delta\ln v_{w0}
-2\,\delta\ln a
-2\,\delta\ln L_W,
\]
\[
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
So \(\delta_\perp\) becomes an explicit linear combination of
\[
(\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W).
\]

### Stage 148
The exact lower compensated branch fixes
\[
\delta\ln L_W = \delta\ln a,
\]
\[
\delta\ln v_{w0}
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_s^2 c_{s,w}^3}{\rho_w a^5}\right),
\]
\[
\delta\ln\mathcal T_m
=
\frac12\,\delta\ln\!\left(\frac{\mathcal Z_q c_{s,w}^3}{\rho_w c_s^2 a^3}\right).
\]
With the frozen \(n=5\) wall-EOS branch, the irreducible microscopic drift space collapses to
\[
(\delta\ln\mathcal Z_q,\ \delta\ln\rho_w,\ \delta\ln c_s,\ \delta\ln a).
\]

### Stage 149
Those four drifts are exactly inverted by the bundle observables
\[
(\Theta_w,\ K_s,\ K_q,\ P_0),
\qquad
P_0=\frac{N_0}{D_0}.
\]
The exact inversion is
\[
\delta\ln\rho_w = \frac12\,\delta\ln\Theta_w,
\]
\[
\delta\ln a = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w,
\]
\[
\delta\ln c_s = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w + \frac15\,\delta\ln P_0,
\]
\[
\delta\ln\mathcal Z_q = \delta\ln K_q - \frac25\,\delta\ln P_0.
\]

### Stage 150
Every remaining first-order mouth/background drift is an explicit algebraic image of \((\Theta_w,K_s,K_q,P_0)\):
\[
\delta\ln c_{s,w} = \delta\ln\Theta_w,
\qquad
\delta\ln\ell = -\delta\ln\Theta_w,
\qquad
\delta\ln L_W = \frac12\,\delta\ln K_s - \frac14\,\delta\ln\Theta_w,
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
\delta\ln\lambda = \frac12(\delta\ln K_s+\delta\ln K_q).
\]
The tangent-compensation theorem then holds exactly:
\[
\delta\ln r_c = 0,
\qquad
\delta\ln\mathfrak r = 0,
\qquad
\delta\ln\mathfrak g = 0,
\qquad
\delta_\perp = 0.
\]

## What this means

The remaining first-order isotropic problem is now no longer “general branch drift.” It has collapsed to a bundle-observable calculation:

\[
(\Theta_w,\ K_s,\ K_q,\ P_0)
\longrightarrow
\text{all first-order isotropic mouth/core/background drifts}.
\]

And the executable result is stronger than expected: **arbitrary first-order isotropic bundle drift is tangent to the exact compensated Family-1 parent family.**

So the next live theorem gate after Stage 150 is not first-order isotropic transport anymore. It is the first correction that escapes this closed algebra, i.e. the first **off-bundle** slippage.
