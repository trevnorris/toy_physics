# S10 emitted minors and exceptional loci

This extension connects the actual D3 rank-drop minor lists and exceptional
locus descriptions from both focused CAS transcripts to Lean. The checked
conditions range over every real wavevector, including the parallel and
perpendicular directions excluded by the earlier generic basis chart.

## Complete minor families

The new arithmetic bridge covers **83 scalar minor expressions in 12 records**:
30 entries from SymPy and 53 from Wolfram. Each retains its raw arithmetic,
has a checked unit, and is proved equal to its stated determinant.

| Root | Modal minor order | Modal selections | Stacked minor order | Stacked selections |
|---|---:|---:|---:|---:|
| Static | 2 | 9 | 3 | 4 |
| Ordinary | 2 | 9 | 2 | 18 |
| Extra | 2 | 9 | 3 | 4 |

SymPy removes duplicate expressions; Wolfram retains every row/column
selection. The certificate records the index maps explicitly, proves every
mapped determinant identity, and proves that the maps cover their complete
canonical families. Lean also checks that the finite selection tables enumerate
every strictly increasing row and column selection. Thus vanishing of each
printed list is equivalent to vanishing of **all ordered minors of that size**
of the root-substituted matrix or N3 stack. This equivalence holds for both
engines despite their different list lengths.

All root matrices are half the previously verified action matrix, evaluated
at the static, ordinary or extra squared-frequency root. The existing generic
N1/N3 expression bridge supplies the connection to the emitted matrices.
The proof does not infer completeness from the number of printed entries.

Units follow the selected rows: a 2x2 modal minor has dimensions
`L^-6 T^-4 M^2`; a 2x2 minor containing the wavevector row instead has
`L^-4 T^-2 M`. The analogous 3x3 stack dimensions are `L^-7 T^-4 M^2`, while
the all-modal 3x3 selection has `L^-9 T^-6 M^3`. Literal zeros retain explicit
slot annotations, including when one deduplicated zero represents selections
with different units. The earlier general minor-tree theorem covers the unit
of every ordered selection.

## Exceptional loci, including both branches

The bridge additionally parses **12 emitted locus records**. A locus is a union
of branches, and each branch is a conjunction of coordinate equations and any
printed guards. Wolfram's `ConditionalExpression` guards are retained.
Lean proves equivalence between each printed predicate, the complete minor-zero
condition and the following geometry. Coordinates below use the transcript's
one-based notation; the distinguished axis is `k1`.

| Root | Modal rank-drop minor locus | Transverse rank-drop minor locus |
|---|---|---|
| Static | `k = 0` | `k = 0` |
| Ordinary | `k2 = k3 = 0` | `k2 = k3 = 0` |
| Extra | `k2 = k3 = 0` | `k1 = 0` **or** `k2 = k3 = 0` |

The shared matrix/locus certificates assume `rho != 0`, `mu != 0`,
`sigma != 0` and `sigma != 1`; the canonical geometric lemmas state their
weaker individual assumptions. These are exact implications, not a finite
sample of wavevectors.

On `k != 0`, the extra transverse locus splits into two disjoint, nonempty
parts: the perpendicular plane with its origin removed, and the parallel axis
with its origin removed. Further Lean theorems show:

- The perpendicular part makes the transverse minors vanish while the modal
  rank-drop minors do not all vanish.
- The parallel part makes both families vanish.
- Neither static minor family has a rank-drop point with nonzero wavevector.

These facts connect the actual CAS rank tests to the earlier general spectral
classification. They do not rely on a basis vector chosen on a generic chart.

## Targeted records and verification boundary

All four emitted D3 stratum points are parsed and checked. SymPy uses
`(1,0,0)` and `(0,1,0)`; Wolfram uses `(27,0,0)` and `(0,27,-1/2)`.
Lean proves that each point is nonzero, lies on the extra transverse locus,
and belongs to the specified parallel or perpendicular part. The two engines
need not choose the same point.

The proof of locus coverage is universal. The point checks separately verify
the placement of the actual targeted records. They do **not** certify every
matrix, root list, rank integer or multidimensional basis printed during those
exceptional reruns. The subsequent [rerun extension](EXCEPTIONAL_RERUN_RESULT.md)
connects their determinants, roots, matrices, stacks, full bases and residuals,
retaining both parallel vectors and restoring physical units through an
explicit coordinate scale. The subsequent [count bridge](COUNT_RESULT.md)
connects all 70 exceptional rank/count integers. Generic count records and
root metadata remain open.

The coordinate parser is restricted to the printed coordinate/rational grammar
and rejects unsupported forms. Like the arithmetic parser, it remains tested
software; Lean checks the generated predicates and mathematical equivalences.
The bridge does not prove a general CAS stratum-discovery algorithm or certify
the original solver status labels.

The manifest (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_manifest.json`) retains
payloads, hashes, minor selections, complete index maps, predicate trees and
point coordinates. The regression instrument (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_check.py`)
compares the new arithmetic to the existing parsers and checks all 27 coordinate
sign patterns for each locus. This exhausts the distinct truth patterns of the
restricted coordinate-to-zero comparison grammar. Lean independently proves
the geometric equivalences without sampling.

Targeted Lean mutations swap two minor selections, omit either exceptional
branch, and corrupt a conditional guard. The earlier coefficient, unit,
normalization and denominator mutations remain in the same instrument. See
CAS_BRIDGE_VERIFICATION.txt (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/CAS_BRIDGE_VERIFICATION.txt`) and the
check record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_checks.json`) for current
build evidence and outcomes.

At this minor/locus checkpoint, the full build passed 536 selected axiom audits
across 75 canonical Lean files,
using only the standard logical axioms and no proof admissions. That regression
suite rejected all nine mathematical mutations and passed both
positive controls. It included 333 independent arithmetic comparisons,
324 locus sign-pattern comparisons, four point comparisons and 34 rejected
malformed inputs.

Reproduce with the commands in [CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md).
This increment leaves the CAS engines, frozen exports and S11 operands unchanged.
