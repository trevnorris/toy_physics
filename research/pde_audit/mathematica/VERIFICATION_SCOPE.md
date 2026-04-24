# PDE V2 Mathematica Mirror Scope

Snapshot date: 2026-04-24

The Mathematica files in this directory are load-bearing execution mirrors for
the V2 audit stack. They are intentionally treated as secondary CAS coverage,
not as independent derivations.

Mirrored coverage:

- V2-04 open finite-radius exit / impedance patch.
- V2-13 grouped normalization and constant-prefactor identities.
- V2-16 branch-freeze / no-refit target-slot Jacobian.
- V2-17 weak-axisymmetric grouped-P2 splitting.
- V2-19 isotropic full-bundle target surface.
- V2-21 through V2-22C fixture-backed extraction and handoff checks.
- V2-23 formula audit plus solver residual-packet sanity checks.

The Python scripts remain the primary implementation for the heavier profile
integration, schema conversion, and reduced FEM branch solve. The Mathematica
mirrors verify the algebraic gates and fixture contracts that are most likely to
hide a false positive if they drift.
