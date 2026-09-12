# S10 exceptional rerun matrices and complete bases

This extension connects the actual parallel and perpendicular D3 rerun
arithmetic in both focused CAS transcripts to the Lean action and spectral
classification. It adds **68 tagged records and 338 scalar expressions**.
This matrix/basis increment brought the bridge to **128 arithmetic records
and 671 expressions**, plus 12 locus predicates and four targeted points. The
subsequent [count extension](COUNT_RESULT.md) adds 70 exceptional and 42 generic
count records. The [root-list](ROOT_RESULT.md) and
[coincidence](COINCIDENCE_RESULT.md) extensions, followed by the
[metadata extension](METADATA_RESULT.md), bring the current total to
**289 arithmetic records and 916 expressions**, plus three empty discarded-root
lists, 56 primary coincidence logical/container records and 48 metadata payloads.

## Imported expressions

| Emission | Records | Scalar expressions |
|---|---:|---:|
| Q3 determinant at each point | 4 | 4 |
| Ordered root formulas | 4 | 10 |
| N1 matrices at each root | 10 | 90 |
| N3 matrices with the coordinate row appended | 10 | 120 |
| N5 matrix times coordinate vector | 10 | 30 |
| N6 complete printed basis families | 10 | 36 |
| N6 basis dot products | 10 | 12 |
| N6 longitudinality residuals | 10 | 36 |

Every imported scalar has a proof of its value and units. The strict parser
retains the emitted arithmetic and denominator conditions. Whole-matrix,
whole-stack and whole-basis binding theorems connect those cells to their
independently defined mathematical references.

SymPy's points are `(1,0,0)` and `(0,1,0)`. Wolfram's are `(27,0,0)` and
`(0,27,-1/2)`. Their nonvanishing and stratum placement were proved in the
[minor/locus extension](MINOR_LOCUS_RESULT.md).

## Both vectors at the repeated root

For each of the ten engine/root cases, Lean proves that every reference basis
vector belongs to the kernel, that the full family is linearly independent,
and that its cardinality equals the kernel dimension from the earlier
spectral classification. These three facts prove equality between its span
and the entire kernel. The emitted family is then connected to that proof.

At the parallel nonzero root, SymPy prints `(e2,e3)` and Wolfram prints
`(e3,e2)`, using one-based coordinate labels. Both orders are retained. Each
family has two independent vectors and spans the two-dimensional transverse
kernel. Checking one vector would not establish this result.

At either perpendicular point, there are three distinct reference root
formulas: the static root, the ordinary root and the extra root. Each kernel
has dimension one. The static kernel has zero transverse dimension; the two
nonzero kernels each have transverse dimension one. Wolfram's printed static
vector `(0,-54,1)` and ordinary vector `(0,1/54,1)` retain their exact
normalizations.

The complete-kernel theorems assume `rho != 0`, `mu != 0`, `0 < sigma`, and
`sigma != 1`. Individual arithmetic theorems record their weaker denominator
conditions. No generic chart denominator is imposed on these exceptional
basis families.

## Units after numerical specialization

A fixed numerical coordinate does not retain the units of a symbolic
wavevector. The bridge therefore declares an explicit convention:

```
k = κ p,       omegaSquared = κ² z,
[κ] = L^-1,    [p] = 1,      [z] = L² T^-2.
```

Here `p` is the printed rational point. The transcript's frequency symbol is
interpreted as the reduced variable `z` in this part of the bridge. Density,
stiffness and anisotropy retain their original units. Lean checks all raw
arithmetic under this coordinate convention and proves how to restore the
physical quantities:

| Printed coordinate quantity | Physical quantity |
|---|---|
| Root `z` | `κ² z` |
| Modal matrix `B` | `κ² B` |
| Determinant | `κ⁶ det B` |
| Appended coordinate row `p` | `κ p` |
| N5 product `B p` | `κ³ B p` |
| Basis vector | Same dimensionless direction |
| Dot product `p · a` | `κ (p · a)` |
| Residual `|p|² a - (p · a)p` | `κ²` times that residual |

These are polynomial identities with the previously proved physical reference
matrix, not an assignment of physical units to nonzero literals. Dedicated
unit theorems check the scale factors against the original dimensions.
Literal zeros retain the existing explicit typed-zero adapter policy.
For `κ != 0`, the physical matrix has precisely the same complete kernel as
the coordinate matrix. The generated `*_physical_complete` theorems connect
each emitted basis to the physical matrix at `κ p` and its scaled root.

## Verification and remaining boundary

The generator is [S10_lean_cas_reruns.py](../../scripts/S10_lean_cas_reruns.py),
invoked through the existing main bridge generator. The
[manifest](../../_measurements/S10_lean_cas_bridge_manifest.json) records raw
payloads, parsed trees, source hashes, units, domains and reference equalities,
including the coordinate convention. The parser-to-transcript correspondence
remains tested software; Lean checks the generated mathematical claims.

The regression suite compares all 916 arithmetic cells with the independent
existing parsers. New Lean mutations duplicate a parallel basis vector,
replace its kernel dimension by one, use an incorrect frequency scale, and
use an incorrect scale-unit exponent. Positive controls check the original
basis and scaling proofs. The parser also rejects removal of a parallel basis
vector from the actual transcript. The [count extension](COUNT_RESULT.md)
also tests altered counts, signed residuals and a generic count extended
beyond its chart. The root-list extension tests algebraic multiplicity,
complete root multisets, distinct counts and the raw-tree filter. The coincidence
extension checks all 56 selected logical payloads and tests guards, allowed
regions, decisions and witnesses. The metadata extension compares 48 complete
payloads and tests repeated fields, aggregate equations, root signs and skipped
strata. All 36 Lean mutations are
rejected and all 14 positive controls pass, alongside 80 malformed-input
rejections.
The full build passes **1962 selected axiom audits** across **115 canonical Lean
files** and 3834 build jobs; only the standard logical axioms occur, and there
are no proof admissions. See [CAS_BRIDGE_VERIFICATION.txt](CAS_BRIDGE_VERIFICATION.txt)
for the current full build and [the check record](../../_measurements/S10_lean_cas_bridge_checks.json)
for exact regression outcomes. Reproduction commands are in
[CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md).

The [count extension](COUNT_RESULT.md) connects all 70 printed N2/N3/N4/N7
records across these ten engine/root cases to the actual matrices and bases.
It also certifies the 42 generic count records on their explicit chart.
The [root-list extension](ROOT_RESULT.md) certifies solution lists, distinct-root
counts, algebraic multiplicities and the syntactic candidate filter. The
[coincidence extension](COINCIDENCE_RESULT.md) connects primary equations,
guarded loci, decisions and witnesses. The [metadata extension](METADATA_RESULT.md)
connects Q8/aggregate fields, root signs, root-condition lists, spectrum solve
operands/statuses and stratum dispositions. Reality-filter traces and remaining
Q5/Q6/period-average metadata remain open. D4, other action packages,
Q7 production alignment and broader comparator/export integration also remain
open in [COVERAGE.md](COVERAGE.md).

This certifies the printed points and their nonzero scale families. The earlier
spectral and locus proofs provide universal stratum classification; these
point bindings do not prove a general CAS discovery or sampling algorithm.
The CAS engines, frozen exports and S11 files are unchanged by this extension.
