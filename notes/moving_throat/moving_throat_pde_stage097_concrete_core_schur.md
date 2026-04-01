
# Moving-Throat PDE — Stage 97: Concrete Two-Channel Core Outlet Model

## Goal

Replace the reduced outlet coefficients `(\rho_R,\sigma_W,\kappa_W,\gamma_W)` by a concrete throat-core response model with explicit internal variables.

## Core variables

Let `u(\omega)` be the isotropic `l=2` mouth amplitude seen by the exterior branch. Introduce two internal core variables:

- `s(\omega)`: a static geometric compliance coordinate,
- `q(\omega)`: a mixed `A_w/F_{\mu w}`-type side-channel coordinate.

Take the linear core system
\[
\begin{pmatrix}
K_s & \lambda\\[2pt]
\lambda & -K_q D_W^{\rm bare}(z)
\end{pmatrix}
\binom{s}{q}
=
u\binom{g_s}{g_q},
\qquad
z:=\frac{a\omega}{c_s},
\]
with bare mixed denominator
\[
D_W^{\rm bare}(z)=1-\kappa_0 z^2-i\gamma_0 z^5+O(z^6).
\]
The mouth feedback is defined by
\[
\delta\Lambda_{\rm core}(z)\,u = g_s s + g_q q.
\]

This is the simplest concrete isotropic core model that contains:
- a static geometric compliance,
- a dynamic mixed side-channel,
- and a nontrivial static/mixed hybridization `\lambda`.

## Exact Schur-complement outlet

Eliminating `(s,q)` gives the exact core correction
\[
\boxed{
\delta\Lambda_{\rm core}(z)
=
\frac{g_s^2}{K_s}
-
\frac{(K_s g_q-\lambda g_s)^2}
{K_s\big(K_sK_q D_W^{\rm bare}(z)+\lambda^2\big)}.
}
\]

Define the dimensionless hybridization ratio
\[
r_c:=\frac{\lambda^2}{K_sK_q}.
\]
Then the outlet takes the reduced Stage-95 form
\[
\boxed{
\delta\Lambda_{\rm core}(z)
=
\rho_c
-
\frac{\sigma_c}{1-\kappa_c z^2-i\gamma_c z^5}
+O(z^6),
}
\]
with exact core-level identifications
\[
\boxed{
\rho_c=\frac{g_s^2}{K_s},
}
\]
\[
\boxed{
\sigma_c=
\frac{(K_s g_q-\lambda g_s)^2}
{K_s^2K_q(1+r_c)},
}
\]
\[
\boxed{
\kappa_c=\frac{\kappa_0}{1+r_c},
\qquad
\gamma_c=\frac{\gamma_0}{1+r_c}.
}
\]

So the concrete core model reproduces the full reduced Robin–mixed hybrid outlet with no extra assumptions.

## Interpretation

This is already a significant narrowing. The outlet is no longer described by four unrelated reduced coefficients. It is controlled by:

- one static shell stiffness `K_s`,
- one mixed-channel stiffness `K_q`,
- one static/mixed hybridization `\lambda`,
- and two mouth couplings `(g_s,g_q)`,
- plus the bare mixed low-frequency pair `(\kappa_0,\gamma_0)`.

The next question is whether this concrete core model can *naturally* land on the compensated canonical branch found algebraically in Stage 95.
