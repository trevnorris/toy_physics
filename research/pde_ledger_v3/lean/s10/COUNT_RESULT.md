# S10 generic and exceptional rank, nullity and count bridge

This extension certifies all **112 printed N2/N3/N4/N7 count records** in the
focused anisotropic D3 transcripts: 42 generic records and 70 parallel and
perpendicular rerun records, across eight root cases per engine. The
subsequent [root-list](ROOT_RESULT.md), [coincidence](COINCIDENCE_RESULT.md), and
[metadata](METADATA_RESULT.md) extensions bring the combined arithmetic bridge
to **289 records and 916 scalar expressions**, plus 12 exceptional-locus
predicates, four targeted points, three empty discarded-root lists, 56 primary
coincidence logical/container records and 48 additional metadata payloads.

## Meaning of each count

| Emission | Records | Kernel-checked meaning |
|---|---:|---|
| N2 rank | 16 | `Matrix.rank` of the actual imported N1 matrix |
| N2 nullity | 16 | Dimension of that matrix map's kernel |
| N3 stacked rank | 16 | `Matrix.rank` of the actual imported N3 stack |
| N3 transverse nullity | 16 | Dimension of the stacked matrix's kernel |
| N4 nullity difference | 16 | Signed N2 nullity minus N3 transverse nullity |
| N7 basis count | 16 | Cardinality of the actual printed basis family |
| N7 count residual | 16 | Signed basis count minus N2 nullity |

The expression parser checks exact integral constants and dimensionless units.
Each record then has a separate semantic theorem identifying its value with
the matrix or basis quantity above. The manifest records both stages and the
coefficient domain. Agreement with an expected integer alone is not used as
proof of the mathematical count.

CountSupport.lean (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/CountSupport.lean`) proves that appending a
wavevector row gives exactly the intersection of the modal kernel with the
Euclidean transverse space. The generated CountReference.lean (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/CountReference.lean`)
connects both imported matrix maps to the previously classified spaces, then
uses their proved dimensions. Rank follows from rank-nullity with domain
dimension three. This applies to the 4-by-3 N3 matrix as well: four rows do
not change the dimension of its input space.

## Results on the generic chart

| Case | N2 rank | Nullity | N3 rank | Transverse nullity | N4 difference | Basis count | N7 residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| Static root | 2 | 1 | 3 | 0 | 1 | 1 | 0 |
| Ordinary root | 2 | 1 | 2 | 1 | 0 | 1 | 0 |
| Extra root | 2 | 1 | 3 | 0 | 1 | 1 | 0 |

Each generic semantic binding assumes `rho != 0`, `mu != 0`, and
`GenericChart sigma k`: `0 < sigma`, `sigma != 1`, and all three components
of `k` nonzero. These are explicit hypotheses, also recorded in the manifest.
The chart supports the displayed basis normalizations. It is not claimed to
be the largest domain on which an individual count holds.

GenericCountReference.lean (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/GenericCountReference.lean`) assembles
each engine's complete N1 and N3 matrices directly from the already imported
expression trees and proves their whole-matrix reference equalities. It
identifies their kernels with the classified modal and transverse subspaces,
then derives all ranks, nullities and signed counts. Each imported basis is
also proved independent and complete for its corresponding N1 kernel.
GenericCountSupport.lean (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/GenericCountSupport.lean`) derives the
generic dimensions from the existing static and oblique census theorems.
These formulas use the physical wavevector and squared-frequency convention
of the generic expressions; no numerical point substitution is involved.

The extra-root transverse nullity is zero on this chart and one on the
perpendicular stratum below. The boundary regression proves that the same
generic N3 matrix formulas, evaluated at the perpendicular target, give
nullity one; an attempted extension of the generic zero count is rejected.
This tests the domain boundary independently of changing an emitted integer.

## Results on the exceptional strata

The counts agree across both engines, although the targeted points and some
basis normalizations differ.

| Case | N2 rank | Nullity | N3 rank | Transverse nullity | N4 difference | Basis count | N7 residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| Parallel static root | 2 | 1 | 3 | 0 | 1 | 1 | 0 |
| Parallel nonzero root | 1 | 2 | 1 | 2 | 0 | 2 | 0 |
| Perpendicular static root | 2 | 1 | 3 | 0 | 1 | 1 | 0 |
| Perpendicular ordinary root | 2 | 1 | 2 | 1 | 0 | 1 | 0 |
| Perpendicular extra root | 2 | 1 | 2 | 1 | 0 | 1 | 0 |

The semantic theorems assume `rho != 0`, `mu != 0`, `0 < sigma`, and
`sigma != 1`. They quantify over that coefficient domain at each exact printed
point. The earlier [complete-basis bridge](EXCEPTIONAL_RERUN_RESULT.md) proves
independence and completeness of the displayed families, including both
vectors at each parallel nonzero root.

SymPy obtains its N7 count from a separate `DomainMatrix` nullspace computation;
Wolfram uses the length of its displayed basis. The new certificate proves that
each printed number equals both the displayed family's cardinality and the
mathematical nullity on the stated domain. It does not verify either CAS
nullspace implementation.

## Signed residuals and physical scaling

N4 and N7 subtraction takes place in the integers. In particular, retaining
one vector from a two-dimensional parallel kernel produces N7 residual `-1`.
Natural-number subtraction would truncate this to zero and conceal the error.
The regression suite checks this case explicitly.

The numerical specialization uses the existing convention `k = κ p` and
`omegaSquared = κ² z`. The modal block scales by `κ²`; the appended constraint
row scales by `κ`. For `κ != 0`, Lean proves equality of both kernels before
and after scaling, hence equality of their nullities and ranks. Every exceptional
case has a `*_physical_counts` theorem carrying all four N2/N3 counts to the
physical variables. Basis cardinality remains the same, and the signed
residuals consequently have the same values.

## Verification and remaining work

S10_lean_cas_counts.py (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/S10_lean_cas_counts.py`) and
S10_lean_cas_generic_counts.py (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/S10_lean_cas_generic_counts.py`)
generate the exceptional and generic arithmetic and semantic bindings through
the main bridge generator. The
manifest (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_manifest.json`) preserves tags,
line numbers, exact payloads, raw trees, units, source hashes, mathematical
meanings and proof names for every record.

The regression instrument compares all 916 arithmetic expressions using the
independently implemented existing parsers. New mutations change a printed
modal nullity, confuse static modal and transverse nullities, use four as the
N3 input dimension, and truncate a negative N7 residual. A positive control
proves that the incomplete basis produces `-1`. The restricted count grammar
also rejects fractional, symbolic, out-of-bound and container payloads.
The generic extension rejects altered SymPy and Wolfram count records and
an invalid extension across the perpendicular boundary; its positive control
proves the correct count at that boundary.
The subsequent root-list tests also cover incomplete multisets, incorrect
multiplicities, distinct-root counts and raw-tree filtering. The complete
regression suite passes 916 arithmetic comparisons, 324 locus sign-pattern
checks, four point comparisons, three filter-list comparisons and 56 complete
coincidence payload comparisons, plus 48 complete metadata payload comparisons.
It rejects all 80 malformed
inputs and all 36 mathematical mutations; all 14 positive controls pass.
Canonical source and input hashes remain unchanged by the tests.

The full build passes **1962 selected axiom audits**, including **1712 CAS
audits**, across **115 canonical Lean files** and 3834 build jobs. All 112
semantic count bindings are included in that audit. Only `propext`,
`Classical.choice` and `Quot.sound` occur; warnings are errors and there are
no proof admissions. Full-build and test results are recorded in
CAS_BRIDGE_VERIFICATION.txt (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/CAS_BRIDGE_VERIFICATION.txt`) and
the check record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_checks.json`).
Reproduction commands are in [CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md).

The [root-list extension](ROOT_RESULT.md) now certifies the raw solution lists,
distinct-root counts, algebraic multiplicities and Wolfram's syntactic candidate
filter records, including three candidates and two distinct roots at the
parallel point. The [coincidence extension](COINCIDENCE_RESULT.md) now connects
primary root differences, guarded loci, allowed regions, decisions and witnesses.
The [metadata extension](METADATA_RESULT.md) connects Q8/aggregate fields, root
signs, root-condition lists, spectrum solve operands/statuses and stratum
dispositions. Reality-filter traces and remaining Q5/Q6/period-average metadata
still need integration. Certification of the original solver implementations remains separate.
D4, other action packages, Q7 production alignment, broader
comparator/export integration and ledger/paper reconciliation remain open in
[COVERAGE.md](COVERAGE.md).

These proofs certify the selected emissions and their declared mathematical
meaning. The transcript parser remains tested software. The CAS engines,
frozen exports and S11 operands are unchanged by this extension.
