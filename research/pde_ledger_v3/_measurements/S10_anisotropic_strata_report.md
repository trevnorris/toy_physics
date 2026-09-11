# S10 anisotropic exceptional-stratum repair

2026-09-11. Focused SymPy/Wolfram reruns and a Lean-derived coverage contract.
Included in the [S9/S10 checkpoint](../lean/CHECKPOINT.md). The original full-sweep transcripts, broad
comparator output, and `S10_exports.py` have not been regenerated.

## Finding and repair

The [Lean anisotropic classification](../lean/s10/ANISOTROPIC_RESULT.md)
identified a missing exceptional direction in the original CAS discovery.
With positive one-axis inertia scale `sRho != 1` and nonzero real wavevector,
the extra propagating mode becomes exactly transverse when `k1=0`. Its
unstacked matrix nullity does not change there, and its frequency stays
distinct from the ordinary branch. Looking only for unstacked rank drops and
root coincidences therefore misses this change in N3.

Q8a in the shared specification now requires the same minor-based rank-locus
calculation for the transverse stack `[M_r; k^T]`, using that stack's own
computed generic rank. Both engines implement this second discovery family.
The resulting allowed regions feed the existing point substitution and full
matrix/nullspace reruns; the perpendicular root formula is not inserted as an
answer.

The Wolfram implementation also partitions the discovered regions by their
Boolean membership combinations before sampling. A witness for a union can
otherwise land on the parallel branch and miss the perpendicular branch.
The new partition separates the currently discovered memberships and rejects
an undecided partition operationally. It is not a general irreducible
decomposition or a proof that every remaining region has constant rank.
One witness per component would not by itself establish that either.

## Focused runs and comparison

Commands run from the repository root, both with exit status 0:

```sh
python3 -u research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py --package XFORM_ANISO --no-export > research/pde_ledger_v3/scripts/out/S10_anisotropic_strata_sympy_audit.out
S10_PACKAGES=XFORM_ANISO WolframKernel -script research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl > research/pde_ledger_v3/mathematica/out/S10_anisotropic_strata_mathematica_audit.out
python3 research/pde_ledger_v3/scripts/S10_anisotropic_strata_comparator.py research/pde_ledger_v3/scripts/out/S10_anisotropic_strata_sympy_audit.out research/pde_ledger_v3/mathematica/out/S10_anisotropic_strata_mathematica_audit.out > research/pde_ledger_v3/scripts/out/S10_anisotropic_strata_comparator.json
```

The comparator exits 0 and records PASS for ten root cases. Each engine
discovers and reruns both directions at D=3 and D=4. The engines choose
different nonzero witness magnitudes; comparison uses `r/|k|²` and matches
direction geometry, never the incidental stratum numbering.

| D | Direction | `r/|k|²` | N2 / N3 |
|---:|---|---|---|
| 3 | parallel | 0 | 1 / 0 |
| 3 | parallel | `muR/rhoBr` | 2 / 2 |
| 3 | perpendicular | 0 | 1 / 0 |
| 3 | perpendicular | `muR/rhoBr` | 1 / 1 |
| 3 | perpendicular | `muR/(rhoBr*sRho)` | 1 / 1 |
| 4 | parallel | 0 | 1 / 0 |
| 4 | parallel | `muR/rhoBr` | 3 / 3 |
| 4 | perpendicular | 0 | 1 / 0 |
| 4 | perpendicular | `muR/rhoBr` | 2 / 2 |
| 4 | perpendicular | `muR/(rhoBr*sRho)` | 1 / 1 |

The [focused comparator result](../scripts/out/S10_anisotropic_strata_comparator.json)
also retains each actual wavevector and unnormalized root. It independently
checks the emitted determinant's complete root list at each witness, matrix
and stack dimensions, ranks/nullities, N4 differences, and N6 basis rank and
annihilation; the cross-engine join compares normalized roots and N2/N3/N4/N7.
It requires the four direction/dimension cases supplied by the Lean
classification. It does not infer a global completeness claim from those
finite samples or prove equivalence of the engines' general Boolean loci.

Wolfram prints 20 `Solve::svars` diagnostics for underdetermined solutions.
The comparator accounts for that exact diagnostic and rejects other untagged
output; local solver diagnostics remain in the transcript.

The [validation instrument](S10_anisotropic_strata_check.py) has eight
[recorded checks](S10_anisotropic_strata_checks.json): the canonical pair
passes; swapping Wolfram stratum numbers still passes; corrupting a transverse
count, omitting the perpendicular rerun, and corrupting a root list are each
rejected. All 48 selected generic rows in each engine match the original
full-sweep transcript byte-for-byte. A subset invocation without `--no-export`
is rejected before touching the export. This avoids replacing a full export
with the focused subset.

## Proof and integration boundaries

The general direction classification is proved in Lean; the CAS result here
is a sampled implementation check against that classification. The scalar
coefficient/sign controls and anisotropic scaling extension are described in
[SCALAR_RESULT.md](../lean/s10/SCALAR_RESULT.md). The full Lean build succeeds
with 250 audited declarations after the later Q6/Q7 and matrix/basis extensions,
and only the standard logical axioms; pinned
source hashes are recorded in its verification files.

`S10_exports.py` retains SHA-256
`bc8de16bae05dcf6caa71d82184f5aa95e2a9d6fd157fdbb674b88f185ed34c9`.
The [run manifest](S10_anisotropic_strata_runs.json) records the focused
sources, outputs and pre-change digests. This increment has not modified the
frozen S11c-b/c1/c2 inputs or the other session's S11c-d files.

The same sampling risk can apply to S11 spectral/sheet constructions even
without a shared-file change. This S10 result does not certify those
constructions. Thresholds, branch intersections, defective roots and bound
poles require explicit domain accounting in their own stage.

The Lean Q6 dimensional-analysis and Q7 Levi-Civita extensions are now complete
for the objects described in the [formal reports](../lean/s10/Q6_Q7_RESULT.md).
The actual CAS expression bridge, production Q7 alignment, full CAS rerun and
broad-comparator/export integration remain in the
[S10 coverage map](../lean/s10/COVERAGE.md).
