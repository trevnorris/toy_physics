# S9/S10 Lean checkpoint — 2026-09-11

This records the original checkpoint. For current work, follow
[FORMALIZATION_POLICY.md](FORMALIZATION_POLICY.md) and
[s10/COVERAGE.md](s10/COVERAGE.md); systematic CAS bridge expansion is no longer
part of the Lean completion plan. Historical evidence below is preserved.

This checkpoint includes the pinned Lean/PhysLean environment, S9 and S10 proof
sources, verification records, and the S10 changes outside `lean/`: the shared
rank-stratum requirement, both CAS engines, focused transcripts/comparator,
validation instruments and ledger updates.

## Completion status

**S9's original formalization pilot is complete.** It proves the supplied
three-dimensional action's integrated first variation, equivalence of compact-test
stationarity to the local PDE, plane-wave reduction, and complete mode census
within that ansatz. The later S10 library proves exact D=3 agreement with these
definitions. The [S9 coverage map](s9/COVERAGE.md) distinguishes that completed
pilot from the full ledger step: the scalar-GNLS no-transverse-mode argument and
certification of the original S9 CAS/export chain remain open.

**S10 is not complete end to end.** Its Lean mathematical core now covers all
six action families, arbitrary-dimensional mode classifications, anisotropic
exceptional directions, coefficient/sign controls, dimensions, Levi-Civita Q7,
modal matrices, minors, complete basis constructions and residual dimensions.
The remaining work is to connect actual CAS expressions/emissions to the proved
objects, align production Q7 and route normalizations, rerun and integrate the
full comparator/export chain, and reconcile ledger/paper claims. See the
[S10 coverage map](s10/COVERAGE.md).

The given action, constant coefficients, physical field interpretation and
choice of physical dimension remain premises. Deriving those premises from a
microscopic medium is outside the completed formalization. Neither pilot proves
general Fourier completeness, finite-boundary/interface physics or weak-solution
extensions.

## Evidence included

- Five Lean libraries build with warnings treated as errors: 250 audited
  declarations across 60 canonical `.lean` files. No proof admissions occur;
  audited dependencies use only `propext`, `Classical.choice` and `Quot.sound`.
  The [audit record](s10/AUDIT_VERIFICATION.txt) captures the successful 3779-job
  build and source/dependency hashes.
- The focused SymPy/Wolfram anisotropic repair discovers both the parallel and
  perpendicular exceptional directions. Its [comparator](../scripts/out/S10_anisotropic_strata_comparator.json)
  passes ten root cases; eight [validation checks](../_measurements/S10_anisotropic_strata_checks.json)
  include deliberate corruptions, stratum renumbering and the export guard.
- The [Q6/Q7](../_measurements/S10_lean_q6_q7_checks.json) and
  [matrix/minor/basis](../_measurements/S10_lean_matrix_checks.json) instruments
  record five rejected source mutations, using isolated copies. Earlier S9
  and anisotropic mutation checks are documented in their proof reports.

Before checkpointing, the verification hashes were checked against the current
sources and dependency pins, and all eight focused CAS validation checks were
rerun successfully. The focused transcripts are sampled implementation evidence;
the general anisotropic classification comes from Lean.

The original full-sweep CAS outputs and `S10_exports.py` are retained.
`S10_exports.py` has SHA-256
`bc8de16bae05dcf6caa71d82184f5aa95e2a9d6fd157fdbb674b88f185ed34c9`.
The checkpoint contains no S11 changes. Build caches, temporary mutation copies
and transient logs are ignored; durable evidence is in the committed reports
and manifests.
