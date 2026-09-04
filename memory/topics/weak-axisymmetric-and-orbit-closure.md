---
title: Weak-axisymmetric and orbit closure
type: topic
status: current
sources:
  - research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex
  - research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex
  - research/pde_audit/scripts/stage_v2_17_weak_axisymmetric_splitting_sympy_audit.py
  - research/pde_audit/notes/stage_v2_18_monomial_quotient_similarity_orbit_derivation.md
  - research/pde_audit/notes/stage_v2_26_program_status_after_audit.md
  - software/stage1_solver/README.md
  - software/stage1_solver/derivations/pathA_01_return_source_and_balance.md
  - software/stage1_solver/STAGE1_VERDICT.md
last_updated: 2026-09-03
---

# Weak-axisymmetric and orbit closure

## Scope

This page connects three distinct uses of “orbit”: the older 1PN test-body
orbit, translation/spin response in the N-body extension, and the later
similarity orbit of microscopic response parameters. They are related pieces
of the program but not one derivation.

There is no completed 5PN or orbit-lock result in the permitted sources. The
weak-axisymmetric and similarity formulas are criteria that a future actual
branch must satisfy.

## Older 1PN orbital result

For the static late-time scalar lag, the orbital paper recovers a Poisson
sector under its boundary and regularity assumptions, but finds zero 1PN
precession from that static scalar alone. A velocity-dependent kinetic factor
gives

$$
\Delta\phi=\frac{2\pi\beta GM}{c_s^2a(1-e^2)}.
$$

The choice $\beta=3$ is fixed by Schwarzschild matching. Its decomposition,
including the required $\kappa_{\rm PV}=3/2$ remainder, was not derived from
throat microphysics in this source.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — labels `sec:static_limit_theorem`, `sec:retarded_scalar_1pn`, `eq:Delta_phi_scalar_result`, `eq:Delta_phi_total_beta`, `sec:beta`, and `sec:hydro_status`

## Spin and N-body lane

The spin/N-body paper treats rotation-driven dyon flow and
translation-driven EIH wake as different structures. The spin flow is
phenomenologically calibrated to weak Kerr. Given a long-range translational
wake with the declared tensor form, the N-body construction reproduces the
EIH cross tensor; the existence and microscopic source of that wake remain
PDE-facing open questions.

Its full isotropic-projector treatment also corrects an earlier truncated
wake-basis interpretation: the viable mixing is real with
$\alpha^2=3/4$, not an imaginary/negative-weight mode.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — labels `subsec:spin-dyon`, `eq:D-calibration`, `subsec:spin-mismatch`, `eq:wake-ansatz`, `subsec:nbody-derivation`, and `subsec:nbody-signature`

## Weak-axisymmetric response diagnostic

Inside the grouped reduced response, a first-order pure axisymmetric
quadrupole has lane signature

$$
(\lambda_{20},\lambda_{21},\lambda_{22})=(1,\tfrac12,-1),
\qquad b=3a,
$$

and prefactor slope

$$
\Xi_1=\frac{P_1}{P_0}
=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
$$

The audit verifies the corresponding first-order splitting and conditional
even-sector relations. This is exact algebra within the declared reduced
carrier; it does not show that a nonlinear throat realizes the axisymmetric
tangent or its required scalar loading.

Sources:

- `research/pde_audit/scripts/stage_v2_17_weak_axisymmetric_splitting_sympy_audit.py` — module description and function `main`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — summary row `V2-17`

## Similarity-orbit criterion

For the declared monomial response map, the eight-dimensional microscopic
tangent space splits into a five-dimensional similarity orbit plus three
quotient directions
$q=(q_{\rm tr},q_{\rm nt},q_\eta)$. Under the stated nondegeneracy conditions,
the compiler from these quotient coordinates to
$(\Theta_1,\Xi_1,\mathcal R_1)$ is invertible. Consequently, zero first-order
defects are equivalent to the tangent lying in $\ker M_*$.

This is a classification theorem for a fixed monomial map. It does not prove
that an actual moving-throat branch stays on the zero-defect similarity orbit;
that tangent must be extracted from the branch.

Source:

- `research/pde_audit/notes/stage_v2_18_monomial_quotient_similarity_orbit_derivation.md` — headings `Exact normal basis and similarity-orbit split`, `Exact physical defect compiler`, `Bridge back to the prefactor slope`, and `Stage verdict`

## Branch-realization status

The original Stage-1 validation harness established numerical machinery, not
a physical result. A later frozen effective-closure branch was tested and
robustly missed its target. The return-source derivation closes part of the
reciprocal structure but does not solve the branch. The actual exporter and
weak-axisymmetric tangent therefore remain open.

Sources:

- `software/stage1_solver/README.md` — headings `Build status` and `What the numbers mean — and what they do not`
- `software/stage1_solver/STAGE1_VERDICT.md` — headings `The result — a robust MISS` and `Interpretation — what it means (and does NOT mean)`
- `software/stage1_solver/derivations/pathA_01_return_source_and_balance.md` — heading `Posited or Open`
- `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` — heading `Current status statement`

Related pages:

- `memory/topics/moving-throat-dynamics.md`
- `memory/topics/post-newtonian-ladder.md`
