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
- `S10Audit/CAS/`: generated trees from both engines' generic and exceptional anisotropic D3
  transcripts, checked values and units, route normalization, explicit
  denominator domains and complete displayed bases. See
  [CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md) and
  [CAS_BRIDGE_VERIFICATION.txt](CAS_BRIDGE_VERIFICATION.txt) for the current
  increment; earlier verification files retain the committed checkpoint.
  The [minor/locus extension](MINOR_LOCUS_RESULT.md) additionally certifies all
  printed D3 rank-drop minors, complete selection maps, guarded exceptional
  locus predicates and the four targeted points. The
  [exceptional rerun extension](EXCEPTIONAL_RERUN_RESULT.md) connects their
  matrices, roots and complete bases, including both parallel basis vectors
  and the restoration of physical units after numerical specialization. The
  [count extension](COUNT_RESULT.md) connects all 112 generic and exceptional N2/N3/N4/N7
  records to matrix ranks, kernel dimensions, basis cardinalities and signed
  residuals, with explicit generic chart assumptions and preservation of the
  exceptional counts under nonzero physical scaling.
  The [root-list extension](ROOT_RESULT.md) certifies complete solution lists,
  distinct-root counts, algebraic multiplicities and syntactic filter records.
  The [coincidence extension](COINCIDENCE_RESULT.md) connects primary emitted
  root differences, guarded loci, allowed regions, decisions and witnesses,
  including the full exceptional parallel axis. The
  [metadata extension](METADATA_RESULT.md) checks all coincidence aggregate/Q8
  fields, root signs, spectrum solve operands/statuses, empty root-condition
  lists and retained/skipped stratum dispositions.

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
