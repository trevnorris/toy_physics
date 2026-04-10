
# Stages 199–201 — Bringing the 5PN endgame home

This block stops the long exploratory chain and compresses the remaining theorem gap
into the smallest exact packets the completed moving-throat PDE must still supply.

## Stage 199 — exact final branch residual packet

Take the actual grouped-lane low-frequency bundle data
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad A\in\{20,21,22\},
\]
together with the source-map factor \(m_{\hat 0}\).

Compile the normalized grouped response moments
\[
u_2^{(A)}=-\frac{D_{A2}}{D_{A0}},
\qquad
u_4^{(A)}=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2},
\]
and the outgoing prefactor moments
\[
P_0^{(A)}=\frac{N_{A0}}{D_{A0}},
\]
\[
P_2^{(A)}=\frac{D_{A0}N_{A2}-2D_{A2}N_{A0}}{D_{A0}^2},
\]
\[
P_4^{(A)}=
\frac{D_{A0}^2N_{A4}-2D_{A0}(D_{A2}N_{A2}+D_{A4}N_{A0})+3D_{A2}^2N_{A0}}
{D_{A0}^3}.
\]

Then extract the grouped trace/anomaly data
\[
(\bar u_2,a_2,b_2),\qquad
(\bar u_4,a_4,b_4),\qquad
(\bar P_0,a_{P_0},b_{P_0}).
\]

The exact final branch residual packet is
\[
\Delta_{\rm branch}
=
\bigl(
a_2,\ b_2,\ a_4,\ b_4,\ a_{P_0},\ b_{P_0},\ \Delta_{\rm pole},\ \Delta_{\rm norm}
\bigr),
\]
with
\[
\Delta_{\rm pole}=\bar u_4-4\bar u_2^{\,2},
\]
\[
\Delta_{\rm norm}=m_{\hat 0}^{\,2}\bar P_0-\frac{54Gc_s^5}{5a^5c^5}.
\]

So the completed PDE no longer has to “show 5PN somehow.” It has to drive one exact
finite-dimensional residual packet to zero.

## Stage 200 — exact endgame compiler

The orbit side was already reduced in Stages 181–198 to any one of the equivalent packets
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

This stage combines that orbit packet with \(\Delta_{\rm branch}\) and shows that the
whole reduced 5PN / 2.5PN / 4PN closure problem depends only on

1. the grouped branch packet \(\Delta_{\rm branch}\),
2. the orbit-lock packet \(\Delta_{\rm orbit}\).

It also records one useful practical simplification: on the minimal isotropic conservative
module, the explicit Family‑1 support/source side is already above the required threshold,
because
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
A_{F1}\approx 1.00005192880220 > \frac13.
\]
So the support/source side is no longer the active bottleneck inside the current hierarchy.

The reduced closure problem is therefore exactly
\[
\Delta_{\rm branch}=0,
\qquad
\Delta_{\rm orbit}=0.
\]

## Stage 201 — home-stretch theorem and minimal PDE data packet

The final theorem gate is now stated as a minimal-data problem.

### Packet A — grouped bundle data
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad
A\in\{20,21,22\},
\]
plus \(m_{\hat 0}\).

### Packet B — orbit/invariant data
Any one of
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

Everything else in the reduced closure test is an exact compiler output of these packets.

So the remaining theorem gap is no longer diffuse. The completed moving-throat PDE only
has to supply the data needed to evaluate
\[
\Delta_{\rm branch}
\quad\text{and}\quad
\Delta_{\rm orbit}.
\]

If \(\Delta_{\rm branch}\neq0\), the branch fails the reduced GR test.
If \(\Delta_{\rm branch}=0\) but \(\Delta_{\rm orbit}\neq0\), the branch is isotropic
but off the exact similarity orbit.
Only if both vanish is the reduced closure complete inside the current hierarchy.
