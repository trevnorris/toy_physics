# Moving-Throat PDE — Stage 109: Linearized Branch-Selection Law Near the Canonical Outgoing Branch

## Goal

Extract the first-order deformation law around the canonical outgoing branch and isolate the minimal PDE-facing branch-selection data.

## Linearized deformation ansatz

Expand around the canonical branch by writing
\[
S=1+\varepsilon s,
\qquad
\beta=1+\varepsilon b,
\qquad
\Sigma_0=\varepsilon a_0,
\qquad
\Sigma_5=\varepsilon a_5,
\]
with the even slots adjusted to preserve the canonical conservative fingerprint.

Then the exact Stage 107 deformation formula
\[
\chi_Q=\frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0}
\]
expands to
\[
\boxed{
\chi_Q
=
1+
\varepsilon\left(5b+\frac{a_0}{3}+9a_5\right)
+O(\varepsilon^2).
}
\]

## Immediate implications

### 1. Pure overall scaling cancels
The coefficient `s` drops out to first order. This matches the exact Class-A invariance: overall mouth renormalization does not move the normalized outgoing quadrupole fingerprint.

### 2. Effective argument deformation matters only through `b`
A small outgoing-argument deformation contributes as
\[
\delta\chi_Q^{(\beta)}=5b.
\]
So a genuine effective radius/sound-speed shift that is not removed by the even-moment matching will move the outgoing normalization linearly.

### 3. Static additive slot and odd core slot are the remaining direct branch data
The first-order additive contributions are
\[
\delta\chi_Q^{(0)}=\frac{a_0}{3},
\qquad
\delta\chi_Q^{(5)}=9a_5.
\]
So the minimal isotropic branch-selection triple is
\[
\boxed{(b,\,a_0,\,a_5).}
\]

## Linearized preservation condition

To keep the canonical outgoing normalization at first order, the deformation must satisfy
\[
\boxed{
5b+\frac{a_0}{3}+9a_5=0.
}
\]
Equivalently,
\[
\boxed{
a_5=-\frac{5b}{9}-\frac{a_0}{27}.
}
\]

## Consequence for the reduced theorem program

After the grouped-`P2` conservative split, geometry cleaning, and Family-1 support sufficiency results, the final reduced 2.5PN branch-selection problem can now be stated very sharply:

> compute the isotropic moving-throat DtN branch data `(b, a_0, a_5)` and test whether they satisfy the exact nonlinear condition of Stages 107--108, or at least the linearized condition above.

So the remaining PDE-facing ambiguity is no longer an open-ended “deformed branch somehow.” It is a small explicit set of outgoing-branch deformation scalars.
