
# Moving-Throat PDE — Stage 119: One-Parameter Parent Compensation Family

## Goal

Use the explicit parent-action formulas from Stage 118 to rewrite the compensated canonical outlet branch as a small set of dimensionless microscopic balance conditions.

## 1. Dimensionless parent ratios

Define the two dimensionless parent ratios
\[
\boxed{
\mathfrak r := \frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g := \frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
\]

Then the exact Stage-115 core-balance theorem
\[
g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2
\]
collapses to
\[
\boxed{
1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2.
}
\]

So the compensated canonical branch is not a five-parameter tuning problem. It is a relation between just two normalized parent ratios.

## 2. Exact mouth-coupling balance law

Solving for \(\mathfrak g\) gives the exact two-branch law
\[
\boxed{
\mathfrak g
=
\mathfrak r
\pm
\frac12\sqrt{1+\mathfrak r^2}.
}
\]

This is the parent-action version of the Stage-115 coupling-balance surface.

Interpretation:

- \(\mathfrak r\) is the normalized static/mixed hybridization,
- \(\mathfrak g\) is the normalized relative mouth coupling of the mixed tube to the shell mode.

Once \(\mathfrak r\) is fixed by the core background, the required mouth-coupling ratio is fixed exactly.

## 3. D/N tube selection in parent variables

The Stage-115 D/N-tube condition
\[
\kappa_c=\frac13,
\qquad
r_c=\frac{\lambda^2}{K_sK_q}
\]
is simply
\[
r_c=\mathfrak r^2.
\]
So the auxiliary mixed-tube length becomes
\[
\boxed{
L_W
=
\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
}
\]

This is a strong reduction of the geometric freedom: the D/N mixed-tube length is controlled by the **same** normalized hybridization that enters the coupling-balance law.

## 4. Parent formulas for \(\mathfrak r\) and \(\mathfrak g\)

From Stage 118,
\[
\mathfrak r
=
\frac{-q_* v_{w0}\mathcal I_{sq}}{\sqrt{K_sK_q}}.
\]

Under the simplest uniform-core overlap closure,
\[
\mathcal I_{sq}
=
\frac{8\sqrt2}{3}a^2\ell\sqrt{L_W},
\]
so
\[
\boxed{
\mathfrak r
=
-\frac{8\sqrt2}{3}\,
\frac{q_* v_{w0}\,a^2\ell\sqrt{L_W}}{\sqrt{K_sK_q}}.
}
\]

Likewise,
\[
g_s=\mathcal T_m J_s,
\qquad
J_s=\frac{4\pi a^2\ell}{3},
\qquad
g_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi}{\sqrt2\,L_W^{3/2}}.
\]
Using
\[
K_q=\frac{\mathcal Z_q}{\mu_0}\frac{\pi^2 c_s^2}{4L_W^2},
\]
the normalized mouth ratio becomes
\[
\boxed{
\mathfrak g
=
\frac{\sqrt{2\mathcal Z_q K_s}}
{\mathcal T_m J_s\,c_s\sqrt{\mu_0 L_W}}.
}
\]

So the compensated canonical outlet is selected by just two physical background controls:

- the normalized mixed background flow \(\mathfrak r\),
- the normalized mouth-traction strength \(\mathfrak g\).

## 5. Exact traction law

Solving the balance relation for the required traction amplitude gives
\[
\boxed{
\mathcal T_m
=
\frac{\sqrt{2\mathcal Z_q K_s}}
{J_s\,c_s\sqrt{\mu_0 L_W}}
\,
\frac{1}{\mathfrak r \pm \frac12\sqrt{1+\mathfrak r^2}}.
}
\]

So once the background mixed flow fixes \(\mathfrak r\), the traction amplitude required to preserve the canonical outgoing quadrupole fingerprint is fixed exactly.

## 6. Family-1 / healing-lock shell simplification

If the shell is on the carried healing-locked GNLS branch,
\[
\ell=\frac{\hbar}{2m_\psi c_{s,w}},
\]
then Stage 118 gave
\[
K_s=\frac{3\pi a^2\hbar^2}{5m_\psi\rho_w\,\ell}.
\]
So the only shell-side microscopic inputs left are \((a,\ell,\rho_w)\). All other shell data collapse into fixed prefactors.

## Result

The compensated canonical outlet has now been reduced to a **one-parameter parent family**:

1. choose the normalized mixed background flow \(\mathfrak r\),
2. then the D/N mixed-tube length is fixed by
   \[
   L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}},
   \]
3. and the required normalized mouth coupling is fixed by
   \[
   \mathfrak g=\mathfrak r\pm \frac12\sqrt{1+\mathfrak r^2}.
   \]

So the remaining PDE-facing question is no longer a large outlet-coefficient search. It is whether the actual GNLS + localized-Maxwell throat core picks one of these parent-balance branches.
