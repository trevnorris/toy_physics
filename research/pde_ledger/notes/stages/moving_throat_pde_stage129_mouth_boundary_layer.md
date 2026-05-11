
# Moving-Throat PDE — Stage 231: Explicit GNLS + Localized-Maxwell Mouth Boundary Layer

## Goal

Replace the remaining abstract positive mouth-source family by the first explicit
boundary-layer law derived from a GNLS + localized-Maxwell electrochemical balance
at the mouth.

The point is not to solve the full moving-throat core. The point is to derive an
honest **local mouth source law** from the same parent structures already frozen
in the 4D stack:
- entropic/compressible source thermodynamics from the GNLS matter sector,
- a localized Maxwell scalar potential \(A_0\),
- and the mechanical mouth-traction/confinement channel already carried by
  \(\delta V_{\rm conf}\).

---

## 1. Local mouth free energy and electrochemical potential

Let \(z\in[0,L]\) denote the axial coordinate measured inward from the mouth.

Use the minimal positive source-density free energy
\[
\boxed{
F_{\rm mouth}[\sigma]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\!\left(\ln\frac{\sigma}{\sigma_*}-1\right)
+
V_m(z)\,\sigma
\Big],
}
\]
with
\[
\Theta_\sigma>0
\]
the effective source compressibility/entropic scale.
The normalization imposed below factors out the overall source amplitude, so this
stage is about the **shape law** of the stationary mouth profile rather than the
total amount of source injected into the reduced throat model.

Near the mouth, linearize the parent potential load as
\[
\boxed{
V_m(z)\approx V_1 z,
\qquad
V_1:=\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_z A_0\right|_{\rm m}.
}
\]
So \(V_1\) is the **net mouth-localizing slope** coming from the mechanical mouth
traction and the localized-Maxwell scalar-potential drop.
Using the linear form \(V_m(z)\approx V_1 z\) across the active interval
\([0,L]\) is the reduced mouth-layer closure of this stage. It should be read as
the lowest-order effective potential on the interval where the source profile is
being modeled, not as a claim that the full throat potential is globally linear.

The corresponding electrochemical potential is
\[
\mu_\sigma^{\rm chem}(z)
=
\Theta_\sigma\ln\!\frac{\sigma}{\sigma_*}
+
V_1 z.
\]

---

## 2. Onsager current and stationary zero-flux branch

Take the standard positive Onsager current
\[
J_\sigma
=
-M_\sigma\,\sigma\,\partial_z\mu_\sigma^{\rm chem},
\qquad
M_\sigma>0.
\]
Then
\[
J_\sigma
=
-M_\sigma\left(\Theta_\sigma\,\sigma'(z)+V_1\sigma(z)\right).
\]

On the stationary recirculating mouth branch,
\[
\partial_t \sigma + \partial_z J_\sigma = 0,
\qquad
J_\sigma=0,
\]
so the source law is exactly
\[
\Theta_\sigma\,\sigma'(z)+V_1\sigma(z)=0.
\]

Therefore
\[
\sigma(z)=C\,e^{-V_1 z/\Theta_\sigma}.
\]

Normalizing to one total mouth source on \([0,L]\),
\[
\int_0^L dz\,\sigma(z)=1,
\]
gives the exact positive family
\[
\boxed{
\sigma_{\Pi}(z)
=
\frac{\Pi\,e^{-\Pi z/L}}{L\left(1-e^{-\Pi}\right)},
\qquad
\Pi:=\frac{V_1L}{\Theta_\sigma}>0.
}
\]

So the previously ad hoc truncated exponential family is now the exact zero-flux
equilibrium branch of a GNLS + localized-Maxwell mouth boundary layer.

---

## 3. Physical meaning of \(\Pi\)

The single remaining source-shape parameter is now

\[
\boxed{
\Pi_m=\frac{V_1L}{\Theta_\sigma}
=
\frac{L}{\Theta_\sigma}
\left(
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_z A_0\right|_{\rm m}
\right).
}
\]

Interpretation:

- \(\Pi_m\to 0\): almost uniform source along the throat;
- \(\Pi_m\to\infty\): point-like source concentrated at the mouth;
- finite \(\Pi_m\): positive localized penetration into the throat interior.

So the mouth-source ambiguity has collapsed from a free profile \(\sigma(z)\) to a
single **dimensionless electrochemical bias** \(\Pi_m\).

---

## 4. What this changes

The branch-selection problem is now sharper:

1. positivity already killed the upper compensated branch;
2. the lower compensated branch was already known to be reachable by truncated
   exponentials;
3. now that truncated exponential is no longer just a convenient family — it is
   the exact equilibrium law of a concrete mouth boundary-layer model.

So the next question is no longer “what source family should we try?” It is simply:

\[
\boxed{
\text{what value of }\Pi_m\text{ does the actual moving-throat mouth layer select?}
}
