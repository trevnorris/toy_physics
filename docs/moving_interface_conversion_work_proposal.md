# Moving-interface conversion and work

Status: new paper derivation proposed for Claude-only assessment, October 4, 2026. No scientific worker, restored symbolic objects or numerical calculation. This is a small bridge between existing material balances and face conventions; it supplies no conversion law or throat solution.

The physical question is: **when disordered material becomes brane again, what crosses the moving boundary, and what work could reach a trapped shear wave?** Return already means re-ordering. We are not choosing a new return loop, pump or reservoir. Particle orientation and the desired like-repulsion/opposite-attraction remain requirements of the eventual shared core model, not conclusions of this balance.

## Existing inputs

The [V3 plan, S12](../research/pde_ledger_v3/V3_STEP_PLAN.md) commits dynamical order conversion and rejects a total-density sink with a frozen wall. The amended [stage-006 material record](../research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md), P1, P6, P8 and A1, supplies

\[
\partial_t n+\nabla_4\!\cdot(nv)=0,\qquad
\partial_t(n\chi)+\nabla_4\!\cdot(n\chi v+J_\chi)=n\Gamma.
\]

Here \(v\) denotes the original material velocity \(u\), \(\chi=\chi_B\), and \(\Gamma=\Gamma_{\rm return}-\Gamma_{\rm drain}\) is the **total rate in this balance**. The separately labelled relaxation adjunct is not added a second time. No constitutive choice for \(\Gamma\) or \(J_\chi\) is made.

For constant constituent mass \(m=m_{\rm GNLS}\), put \(\rho=mn>0\) and \(K=mJ_\chi\). This converts number to mass units without changing the original balance. \(K\) is an order-transport current relative to material advection, not an additional total-mass source.

## 1. Keep boundary motion before taking a front limit

For a smooth moving four-dimensional control volume \(\Omega(t)\), outward normal \(N\), and boundary velocity \(b\), the proposed transport identities are

\[
\frac{d}{dt}\int_{\Omega(t)}\rho\,d^4X
=-\int_{\partial\Omega(t)}\rho(v-b)\!\cdot N\,dA,
\]
\[
\frac{d}{dt}\int_{\Omega(t)}\rho\chi\,d^4X
=-\int_{\partial\Omega(t)}[\rho\chi(v-b)+K]\!\cdot N\,dA
+\int_{\Omega(t)}\rho\Gamma\,d^4X.
\]

These retain accumulation and transport through every boundary. The source integral alone is not generally the outward face flux. In particular, re-ordering can advance an interface through material; it need not mean that the material follows a closed circulating trajectory.

Only as an explicit **planar travelling-front limit**, take a profile steady in a frame moving at normal speed \(V\), no tangential variation, ordered interior \(\chi_{\rm in}=1\), and disordered exterior \(\chi_{\rm out}=0\). The normal coordinate \(z\) points from ordered to disordered. Total continuity gives the constant relative mass flux \(j=\rho(v_n-V)\). Integrating the order balance across the front gives

\[
[\chi j+K_n]_{\rm in}^{\rm out}
=\int_{\rm in}^{\rm out}\rho\Gamma\,dz=:G,
\qquad j=K_{n,\rm out}-K_{n,\rm in}-G.
\]

With zero endpoint order-current, \(j=-G\): positive net re-ordering means negative outward relative mass flux. This is a proposed sign consequence under the stated front hypotheses, **not a stationary-front existence result** or an imposed S11c background. A time-dependent or curved front returns to the full balance above. Density may vary through this front; equal endpoint density is not needed for this identity.

The native [S11b convention](../research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md) is \(J_f=\rho_m(v_{{\rm bulk},n_f}-V_f)\), with each face's own outward normal. Identifying a limit of \(j\) with \(J_f\) requires the same material velocity, \(\rho_{\rm out}=\rho_m\), boundary motion and true face measure. The original [S11c-a geometry source](../research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py) exposes these relative-flux and area operands. This proposal does not perform that native curved-geometry join. Both faces use the same outward sign rule; the face label \(f=\pm\) is not the particle's electric orientation label.

There is also an order-of-perturbation distinction: \(j\) above is the total travelling-front flux, whereas the carried S11b face closure describes wave perturbations. A future bridge must first join a background and then linearize consistently; it cannot identify the DC integral \(G\) with the harmonic \(J_f\) by notation.

## 2. Isolate conversion work without calling it mode power

The corrected source convention, including the [handoff erratum](../notes/brane_bulk_handoff.md), is \(\mu_\chi=\delta F/\delta\chi\) and \(P_\chi=\int\mu_\chi D_t\chi\,d^4X\). \(\mu_\chi\) is energy per four-volume, not specific chemical energy. There is no extra factor of number density in this power.

The two balances imply \(D_t\chi=\Gamma-\rho^{-1}\nabla_4\cdot K\). For a smooth finite region at an instant, integration by parts proposes

\[
P_\chi=\int_\Omega
\left[\mu_\chi\Gamma+K\!\cdot\nabla_4(\mu_\chi/\rho)\right]d^4X
-\int_{\partial\Omega}(\mu_\chi/\rho)K\!\cdot N\,dA.
\]

This rewrites the stated **order-work contribution**, not the derivative of all energy inside a moving volume. Material energy transport, kinetic/compressional work, interface stress and any gradient-energy boundary work still need their own terms. The expression does not close a thermodynamic or mode-energy balance. Although \(\mu_\chi/\rho\) has specific-energy units, dimensions do not identify it with native \(\mu_s\) or the affinity \(A_f=\mu_s-p_f/\rho_m\).

There is already a useful face-power identity in [the S11b source, task_b2d](../research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py): with \(v_{\rm bulk}=V+J/\rho_m\),

\[
(p+\Lambda_X A)\overline V+\mu_s\overline J
=p\overline{v_{\rm bulk}}+A\overline J+\Lambda_X A\overline V.
\]

Its real part divided by two is the recorded harmonic mean-power convention. It is an accounting identity with native normalization, not an assumed equality to \(P_\chi\). Frequency-dependent memory/storage and the actual nonlinear background must be joined before interpreting a sustained energy transfer. The earlier strict-rest-bulk wave calculation did not carry a finite background conversion flow.

## 3. What the existing shear gate can and cannot supply

The amended material energy already contains

\[
\chi f_{\rm sh},\qquad
f_{\rm sh}=\tfrac12\mu_R^{(4)}|\operatorname{curl}_4 u_d|^2.
\]

The displacement \(u_d\) is distinct from material velocity \(v\). With the shear configuration and coefficient held fixed, and no extra hidden \(\chi\) dependence assigned to them, this term contributes \(f_{\rm sh}\) to \(\mu_\chi\). Positive shear stiffness makes this particular stored-energy contribution increase when order increases at fixed displacement. That is not a conclusion that re-ordering supplies free energy to a wave, or that it cannot support one once the full state adjusts.

For example, on a material volume moving with \(v\), the product rule gives the proposed bookkeeping identity

\[
\frac{d}{dt}\int_{\Omega_m(t)}\chi f_{\rm sh}\,d^4X
=\int_{\Omega_m(t)}
\left[f_{\rm sh}D_t\chi+\chi D_tf_{\rm sh}
+\chi f_{\rm sh}\nabla_4\!\cdot v\right]d^4X.
\]

Keeping all three terms prevents the first from being renamed net wave power. A mode reduction also needs its kinetic normalization and the actual relation of \(u_d\) to the material/core dynamics. The original bare order well has equal pure-phase minima at fixed density; a net sustaining energy difference is not supplied by that well alone. \(f_{\rm throat}\) and \(f_{\rm mix}\) remain explicit placeholders in the inspected material specification. The historical Maxwell content of \(f_{\rm mix}\) is not a derived native electric coupling.

## Bounded outcome and next physics step

This paper should settle the conditional conversion-to-relative-flux sign, preserve interface motion, and identify the work terms that must be joined. It should not trigger a general verification framework or a numerical front solve. The next physical deliverable, if these balances withstand assessment, is a minimal core/mode work specification from the existing order/shear energy: which displacement stores the support energy, which stress it exerts on the mouth, and how conversion and material transport exchange that energy. If the existing sources cannot supply one of those functions, name that particular missing input and explain concrete physical options to the user.

No fixed-amplitude reservoir, phase-energy bias, mode inertia gate or conversion kinetics is adopted here. Nothing selects the electric force sign, a trapped equilibrium, an expansion rate or a leakage factor. The September [native interpretation](native_light_em_and_vortex_throat_interpretation.md) remains the interpretation baseline. Completed S11c calculations remain preserved; leakage workers, centre-drive implementation and mixed numerical recovery stay parked. A conceptual verdict is not authorization for an unspecified scientific worker.
