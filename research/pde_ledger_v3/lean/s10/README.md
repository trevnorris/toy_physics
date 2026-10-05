# S10 formalization

Follow [the Lean scope policy](../FORMALIZATION_POLICY.md) before resuming work.
The compact S10 Lean contract is complete: both independent fidelity reviews
returned CLEAR, and the proof/build and mutation evidence pass. See
[FIDELITY_REVIEW.md](FIDELITY_REVIEW.md) for the reviewed revision, dispositions
and stopping point. The committed CAS bridge is retained without further
systematic expansion.
See [COVERAGE.md](COVERAGE.md) for the current Lean and production obligations.

S10 proofs and reports live here. The Lean toolchain, dependency pins, and build
cache are shared at [../](../README.md).

- `S10Pilot/` and `S10Pilot.lean`: the arbitrary-dimensional baseline, integrated
  variational principle, phase average, mode census, and exact S9 specialization.
  See [RESULT.md](RESULT.md) and [VERIFICATION.txt](VERIFICATION.txt).
- `S10Controls/`: full-gradient and divergence-only stiffness actions and their
  variational and spectral proofs. See [CONTROLS_RESULT.md](CONTROLS_RESULT.md)
  and [CONTROLS_VERIFICATION.txt](CONTROLS_VERIFICATION.txt).
  Coefficient rescaling and sign reversal are covered in [SCALAR_RESULT.md](SCALAR_RESULT.md).
- `S10Anisotropic/`: one-axis inertia action, integrated variation, split spectrum,
  and oblique, perpendicular, and parallel direction counts. See
  [ANISOTROPIC_RESULT.md](ANISOTROPIC_RESULT.md) and
  [ANISOTROPIC_VERIFICATION.txt](ANISOTROPIC_VERIFICATION.txt).
- `S10Audit/`: PhysLean expression-tree dimensions, action-derived coefficient
  units, root dimensions, matrices, mixed-unit minors, complete basis families,
  residual dimensions and Levi-Civita Q7 comparisons for all six packages.
  See [Q6_Q7_RESULT.md](Q6_Q7_RESULT.md), [MATRIX_RESULT.md](MATRIX_RESULT.md),
  and [AUDIT_VERIFICATION.txt](AUDIT_VERIFICATION.txt).
- `S10Audit/CAS/`: only `PY`, `WL`, `Support`, `BasisCompletion` and `Bindings` remain,
  preserving the limited D3 imported-matrix identity used by the compact C2 fidelity contract.
  The transcript-wide audit, exceptional reruns, count/root/locus records and their bindings
  are at `archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/`.
  Their historical reports are not claims that those modules remain in the current build.

Run from `research/pde_ledger_v3/lean/`:

```sh
lake build S10Pilot
lake build S10Controls
lake build S10Anisotropic
LAKE_CACHE_DIR=.lake/cache lake build S10Audit
```

The analytic setting is flat R^(D+1), smooth backgrounds, and smooth compact test
variations. The supplied action, field content, constant coefficients, and
physical dimension remain premises. The reports distinguish the proved
statements from the remaining S9 and S10 exclusions.
See [COVERAGE.md](COVERAGE.md) for the remaining S10 obligations.
