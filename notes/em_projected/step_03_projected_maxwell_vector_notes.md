# Projected Maxwell Vector Notes

## Purpose

This note records the vector-form projected Maxwell derivation audited by
`step_03_projected_maxwell_vector_sympy.py`.

The main point is that projection-first electrodynamics naturally separates
measured fields from source-coupled flux fields.

## Field Layers

Measured brane fields are built from projected field strength components:

\[
E_{\rm meas}=\operatorname{Proj}_W[F_{a0}],
\qquad
B_{\rm meas}=\operatorname{Proj}_W[F_{ab}].
\]

Source-coupled flux fields are built from projected localized flux components:

\[
D_{\rm flux}=\operatorname{Proj}_W[ZF^{a0}],
\qquad
H_{\rm flux}=\operatorname{Proj}_W[ZF^{ab}].
\]

In general these are not the same objects.

## Homogeneous Sector

The projected Bianchi identities give the usual measured-field homogeneous laws:

\[
\nabla\cdot B_{\rm meas}=0,
\]

\[
\nabla\times E_{\rm meas}+\partial_tB_{\rm meas}=0.
\]

No leakage term appears here because the homogeneous identities project cleanly.

## Inhomogeneous Sector

The projected inhomogeneous equations instead use the flux fields:

\[
\nabla\cdot D_{\rm flux}
+
\operatorname{Leak}_0
+
\operatorname{Gauge}_0
=
\mu_0\rho_{\rm proj},
\]

\[
\nabla\times H_{\rm flux}
-
\partial_tD_{\rm flux}
+
\operatorname{Leak}_{\rm vec}
+
\operatorname{Gauge}_{\rm vec}
=
\mu_0J_{\rm proj}.
\]

The leakage terms come from the projected transverse components `ZF^{w\nu}`.
The gauge terms come from the projected gauge-driver sector.

The updated concrete audit also carries the transverse integration-by-parts
boundary term explicitly and only drops it after verifying that the localized
Gaussian profiles make that boundary flux vanish.

The updated concrete audit now instantiates this with
\[
W(w)=\frac{e^{-w^2}}{\sqrt{\pi}},
\qquad
Z(w)=e^{-w^2},
\qquad
F^{w1}(w)=w,
\]
and checks the genuinely localized leakage
\[
-\operatorname{Proj}_{W'}[ZF^{w1}]
=
\operatorname{Proj}_{W}[\partial_w(ZF^{w1})]
=
\frac{\sqrt2}{4}\neq 0.
\]

For this Gaussian example the accompanying boundary term
\(\operatorname{Boundary}[WZF^{w1}]\) is checked separately and found to vanish,
so the displayed leakage value is not coming from a hidden decay assumption.

It also now checks the homogeneous and inhomogeneous projected vector laws on
an explicit bulk potential
\[
A_0=xz(1+w^2),\qquad
A_1=ty(1+w^2),\qquad
A_2=tz(1+w^2),\qquad
A_3=xy(1+w^2),
\]
with \(A_w=0\). From this concrete field, the script projects the measured
\((E,B)\), the flux fields \((D,H)\), and the projected source terms
\(\rho_{\rm proj},J_{\rm proj}\), and verifies
\[
\nabla\!\cdot B_{\rm meas}=0,\qquad
\nabla\times E_{\rm meas}+\partial_t B_{\rm meas}=0,
\]
\[
\nabla\!\cdot D_{\rm flux}+\operatorname{Leak}_0=\mu_0\rho_{\rm proj},
\qquad
\nabla\times H_{\rm flux}-\partial_t D_{\rm flux}+\operatorname{Leak}_{\rm vec}
=\mu_0J_{\rm proj}
\]
by direct symbolic projection, with the corresponding transverse boundary fluxes
checked to vanish before the compact vector laws are asserted.

## Takeaway

The projected theory looks like an open-system effective medium:

- measured fields satisfy the homogeneous equations,
- source-coupled flux fields enter the inhomogeneous equations,
- hidden-direction leakage controls the failure of the brane theory to close by
  itself.

## Script

- `step_03_projected_maxwell_vector_sympy.py`
