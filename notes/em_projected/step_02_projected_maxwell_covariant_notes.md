# Projected Maxwell Covariant Projection Notes

## Purpose

This note records the covariant projection-first derivation audited by
`step_02_projected_maxwell_covariant_sympy.py`.

The starting point is the localized `4+1` inhomogeneous Maxwell equation in a
generic brane component:

\[
\partial_M\bigl(ZF^{M\nu}\bigr)+\Gamma^\nu=\mu_0J^\nu,
\]

where the script keeps the gauge-driver term generic as `Gamma^nu`.

## Projection Identity

For a normalized observation kernel `W(w)`,

\[
\operatorname{Proj}_W[Q]=\int W(w)Q(w)\,dw.
\]

Projection commutes with brane derivatives because `W` depends only on `w`.
The transverse derivative instead gives

\[
\operatorname{Proj}_W[\partial_w Q]
=
[WQ]_{\partial w}
-
\operatorname{Proj}_{W'}[Q].
\]

Under boundary decay or compact support this becomes

\[
\operatorname{Proj}_W[\partial_w Q]
=
-
\operatorname{Proj}_{W'}[Q].
\]

The updated script now keeps the boundary contribution explicit in the symbolic
identity and only drops it in the concrete Gaussian audit after evaluating the
actual boundary value.

## Projected Inhomogeneous Law

Writing `G^{M\nu}=ZF^{M\nu}`, projection gives the exact brane-component law

\[
\partial_\mu\operatorname{Proj}_W[G^{\mu\nu}]
+
\operatorname{Boundary}[WG^{w\nu}]
-
\operatorname{Proj}_{W'}[G^{w\nu}]
+
\operatorname{Proj}_W[\Gamma^\nu]
=
\mu_0\operatorname{Proj}_W[J^\nu].
\]

With boundary decay:

\[
\partial_\mu\operatorname{Proj}_W[ZF^{\mu\nu}]
-
\operatorname{Proj}_{W'}[ZF^{w\nu}]
+
\operatorname{Proj}_W[\Gamma^\nu]
=
\mu_0\operatorname{Proj}_W[J^\nu].
\]

This is the first clean audit point showing that projection-first Maxwell is an
open-system brane theory unless the transverse leakage term is suppressed.

## Projected Charge Continuity

Bulk current conservation projects to

\[
\partial_t\operatorname{Proj}_W[J^0]
+
\partial_a\operatorname{Proj}_W[J^a]
=
\operatorname{Proj}_{W'}[J^w]
\]

after boundary decay.

So apparent brane charge nonconservation is precisely transverse current
exchange with the hidden `w` direction.

## Script

- `step_02_projected_maxwell_covariant_sympy.py`
