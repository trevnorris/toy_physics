# Parent Throat Action — Full-Bundle Bridge Master Note

## What this continuation adds

This continuation picks up after the parent-throat-action promotion and pushes it into the same bundle language already used by the moving-throat program:

- isotropic full-bundle moments `D0,D2,D4` and `N0,N2,N4`,
- normalized response moments `u2,u4`,
- outgoing prefactor moments `P0,P2,P4`,
- the weak-axisymmetric gate packet `Xi1, K1, H_even`.

The point is to stop treating the wall entries `K` and `M` as abstract closure knobs.  After parent promotion they become explicit branch integrals of the throat action.

## Files in this bundle

- `step_17_parent_throat_action_isotropic_bundle_notes.md`
- `step_18_parent_throat_action_weak_axisym_packet_notes.md`
- `step_19_parent_throat_action_actual_branch_export_notes.md`
- `step_16_parent_throat_action_bundle_master_sympy.py`
- `step_17_parent_throat_action_isotropic_bundle_sympy.py`
- `step_18_parent_throat_action_weak_axisym_packet_sympy.py`
- `step_19_parent_throat_action_actual_branch_export_sympy.py`

## Executive result

Let the promoted throat action be the nonlinear parent wall block

\[
S_\Sigma[R]
=
\int dt\,dw\,d\Omega\;\mathcal L_\Sigma
\]

with quadratic limit determined by

\[
\mu_\eta(w),\qquad T_w(w),\qquad T_\Omega(w),\qquad K_\eta(w).
\]

Then on a chosen grouped-`P2` wall profile `beta_2(w)` the wall enters the isotropic full bundle only through the two parent branch integrals

\[
M_\Sigma = \int dw\,\mu_\eta\,\beta_2^2,
\]

\[
K_\Sigma = \int dw\,\Big[T_w(\beta_2')^2 + (K_\eta+6T_\Omega)\beta_2^2\Big].
\]

With those substitutions the full isotropic bundle is

\[
D_0 = K_\Sigma - B_0 - Z_0,
\qquad
D_2 = -(M_\Sigma + B_2 + Z_2),
\qquad
D_4 = -(B_4 + Z_4),
\]

and the existing one-pole / normalization targets become exact equations for the parent wall block rather than for abstract wall knobs.

The most useful new algebraic consequences are:

### 1. Parent-complete isotropic compatibility surface

The one-pole condition and the outgoing normalization target fix `K_\Sigma` in two exact ways,

\[
K_\Sigma
=
B_0+Z_0
+
\frac{3\,(M_\Sigma+B_2+Z_2)^2}{B_4+Z_4},
\]

\[
K_\Sigma
=
B_0+Z_0
+
\frac{N_0}{P_{0,\rm target}},
\qquad
P_{0,\rm target}
=
\frac{54Gc_s^5}{5a^5c^5\,\widehat m_0^2}.
\]

So the parent-complete branch must satisfy the single compatibility equation

\[
\boxed{
\frac{N_0}{P_{0,\rm target}}
=
\frac{3\,(M_\Sigma+B_2+Z_2)^2}{B_4+Z_4}.
}
\]

### 2. Exact wall-slope solve for the weak-axisymmetric even gates

Write the first-order parent wall slopes as

\[
D_{01}=\delta K_\Sigma-B_{01}-Z_{01},
\qquad
D_{21}=-(\delta M_\Sigma+B_{21}+Z_{21}),
\qquad
D_{41}=-(B_{41}+Z_{41}).
\]

Imposing the live even gates

\[
K_1 = D_{21}+\frac{D_{01}}{9}=0,
\qquad
H_{\rm even}=D_{41}-\frac{2}{3}D_{21}-\frac{D_{01}}{27}=0
\]

solves the wall slopes **uniquely** as

\[
\boxed{
\delta K_\Sigma = B_{01}+Z_{01}+27(B_{41}+Z_{41}),
}
\]

\[
\boxed{
\delta M_\Sigma = -\,(B_{21}+Z_{21}) + 3(B_{41}+Z_{41}).
}
\]

So once the support and conservative Maxwell/mixed anisotropy slopes are known, the parent wall block has no remaining first-order freedom if one insists on the canonical even-compensated branch.

### 3. Residual weak-axisymmetric normalization scalar

After that even-gate compensation, the transported prefactor slope becomes

\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
-
\frac{27(B_{41}+Z_{41})}{K_\Sigma-B_0-Z_0}.
}
\]

So the parent-complete weak-axisymmetric packet again collapses to one scalar normalization defect, but now with the wall contribution written entirely in parent-action branch data.

## Why this matters

This continuation shows exactly what the promoted throat action buys us and what it does not.

What it buys us:

- `K` and `M` are no longer abstract wall coefficients; they are explicit branch integrals of `S_\Sigma`.
- The isotropic one-pole and outgoing normalization targets can be written directly as equations for those branch integrals.
- The weak-axisymmetric wall slopes are no longer free tuning knobs once the even gates are imposed.

What it does **not** buy us:

- It does not determine the support moments `B_n`.
- It does not determine the conservative mixed moments `Z_n`.
- It does not determine the outgoing slopes `N_{0},N_{01}`.

So the remaining gap is exactly what the later moving-throat notes already suggest: actual branch realization, not more algebraic repackaging.
