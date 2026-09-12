# S10 CAS expression bridge: anisotropic D=3 pilot

This increment connects actual arithmetic printed by both CAS engines to the
existing Lean action, spectral and dimensional proofs. It covers the generic
`XFORM_ANISO_D3` records in the two checked-in focused transcripts. The subsequent
[minor and locus extension](MINOR_LOCUS_RESULT.md) adds complete emitted minor
families, exceptional predicates and targeted point records. The
[rerun extension](EXCEPTIONAL_RERUN_RESULT.md) adds exceptional matrices, roots,
complete bases and physical scale restoration. The
[count extension](COUNT_RESULT.md) connects the generic and exceptional rank, nullity and
count records to those objects. The
[root-list extension](ROOT_RESULT.md) certifies the raw solution lists,
distinct-root counts, algebraic multiplicities and syntactic candidate filters.
The [coincidence extension](COINCIDENCE_RESULT.md) connects primary emitted
root differences, guarded loci, allowed regions, decisions and witnesses.
The [metadata extension](METADATA_RESULT.md) checks aggregate/Q8 copies, reported
root signs, spectrum operands and statuses, root solver-condition lists, and
retained/skipped-stratum dispositions.
The [manifest](../../_measurements/S10_lean_cas_bridge_manifest.json) records each
input file hash, tag, line number, exact payload, parsed tree, denominator
condition, unit and generated theorem.

## Certified expressions

Each engine contributes 24 records and 125 scalar expressions, for **48 records
and 250 scalar expressions** in the initial arithmetic bridge:

| Emission | Expressions per engine | Reference checked in Lean |
|---|---:|---|
| Q2 route matrices, difference and ratio | 28 | Actual entry values and route normalization |
| Q3 determinant | 1 | Determinant of the emitted route-B matrix |
| Ordered squared-frequency roots | 3 | Static, ordinary and extra root formulas |
| N1 root-substituted matrices | 27 | Route B evaluated at each of the three roots |
| N3 stacked matrices | 36 | N1 matrix with the wavevector row appended |
| N5 matrix times wavevector | 9 | Full matrix-vector product |
| N6 displayed bases | 9 | Normalized static, ordinary and extra basis vectors |
| N6 dot products and longitudinality residuals | 12 | Every component of the displayed expressions |

The minor extension adds 12 records and 83 scalar expressions, bringing the
minor/locus checkpoint to **60 records and 333 expressions**. A separate
predicate bridge certifies 12 locus records and four targeted points. See
[MINOR_LOCUS_RESULT.md](MINOR_LOCUS_RESULT.md) for the complete selection and
stratum-coverage proofs. The exceptional reruns add 68 records and 338
expressions; their 70 N2/N3/N4/N7 count records and the 42 generic count records
bring the count checkpoint to **240 records and 783 expressions**. The generic
count bindings retain the explicit basis chart and connect complete imported
N1/N3 matrices to the classified modal and transverse subspaces.
The root-list extension adds 27 records and 51 scalar expressions, bringing
the root-list checkpoint to **267 records and 834 expressions**, plus three empty
discarded-root list records. It distinguishes the parallel double root from
the two-element distinct-root list and ties every distinct list to the actual
emitted determinant's complete zero set. The coincidence extension adds 13
arithmetic records and 28 expressions, bringing the current arithmetic total
to **280 records and 862 expressions** at that checkpoint. It also certifies 56 logical/container
records; seven Wolfram operand records receive both arithmetic and logical
checks. Its locus proofs cover arbitrary wavevectors on the declared positive
coefficient domain, including the exceptional parallel axis. The metadata
extension adds nine arithmetic records and 54 expressions, bringing the current
arithmetic total to **289 records and 916 expressions**. It adds 48 complete
metadata payload checks; six Wolfram aggregates receive both arithmetic and
logical checks. All 16 root sign observations are preserved, with independent
proofs of six zero and ten positive mathematical signs.

Every scalar expression has a kernel-checked equality theorem and a
`Expr.HasDim` proof. Arithmetic trees retain the emitted additions,
subtractions, products, quotients and natural powers. Sharing identical
subtrees does not simplify their arithmetic.

The reference matrix is polynomial in an independent real squared-frequency
variable `z`. Both engines satisfy, for every real assignment,

```
M_A = -referenceMatrix
M_B = (1/2) referenceMatrix
M_A - M_B = (-3/2) referenceMatrix
M_A = -2 M_B
```

At `z = omega^2`, the reference matrix equals the previously verified action
matrix. The generated matrix objects from both engines agree entry by entry,
and route B has exactly the action matrix's kernel. No equality is inferred
merely from proportionality or a zero/nonzero residual flag.

## Denominators and complete bases

The generator records nonzero factors from every raw denominator and adds any
conditions needed by the independent reference expression. These appear as
hypotheses in the generated equalities. Root substitutions retain `rho != 0`
and, for the extra root, `sigma != 0`; basis normalizations retain their chart
conditions. Numeric zero denominators are rejected before code generation.

Using the transcript's one-based wavevector coordinates, the displayed bases
are:

| Root | Displayed vector | Chart condition |
|---|---|---|
| Static | `(k1/k3, k2/k3, 1)` | `k3 != 0` |
| Ordinary | `(0, -k3/k2, 1)` | `k2 != 0` |
| Extra | `(-(k2^2+k3^2)/(sigma*k1*k3), k2/k3, 1)` | `sigma*k1*k3 != 0` |

[BasisCompletion.lean](S10Audit/CAS/BasisCompletion.lean) proves nonvanishing,
linear independence and preservation of the relevant spans. The generated
`PY.basis_complete` and `WL.basis_complete` theorems identify the **full kernel**
of each emitted route-B matrix with the span of its displayed vector on the
common generic chart:

```
rho != 0, mu != 0,
0 < sigma, sigma != 1,
k1 != 0, k2 != 0, k3 != 0.
```

These conditions describe the pilot's domain. A zero chart denominator does
not imply absence of a physical mode. Parallel, perpendicular, coalescent and
other exceptional cases have earlier Lean classifications; their actual CAS
rerun matrices and bases are connected separately in the
[exceptional extension](EXCEPTIONAL_RERUN_RESULT.md). It uses the targeted
coordinate points and complete basis families with explicit scale restoration.
The generic chart above does not cover an exceptional stratum.

## Dimensional and trust boundaries

The symbol units specialize the already inferred coefficient dimensions to
`D=3` and a displacement with length units. Lean proves the coefficient-unit
identifications. Basis coordinates are dimensionless; the separately chosen
mode amplitude supplies the field's units. Fixed-point reruns use dimensionless
coordinates and a reduced frequency, with `k = κ p` and `omegaSquared = κ² z`.
Their raw units and proved physical scale restoration are detailed in
[EXCEPTIONAL_RERUN_RESULT.md](EXCEPTIONAL_RERUN_RESULT.md).

A printed literal zero carries no recoverable dimensional annotation. Only
such a scalar slot may receive its declared unit through an explicit typed-zero
adapter. The manifest marks each adapter. Nonzero expressions and internal
additive branches must type-check as emitted; algebraic cancellation cannot
hide a dimensionally invalid branch. Adapter values remain zero without new
denominator assumptions.

The Python parser and transcript-to-tree correspondence are tested software,
not a formally proved lexer/parser. The parser uses a restricted arithmetic
grammar and never evaluates input code. It rejects unknown symbols/functions,
unsupported syntax, floats, unconsumed tokens, malformed containers, duplicate
tags and missing selected records. The generated tree semantics, reference
equalities, units and kernel/basis statements are checked by Lean's kernel.
The bridge certifies the selected printed expressions; it does not prove that
the CAS engines derived or exhaustively emitted every required object.

## Reproduction and verification

From `research/pde_ledger_v3/lean`:

```sh
python3 ../scripts/S10_lean_cas_bridge.py --check
LAKE_CACHE_DIR=.lake/cache lake build
python3 ../_measurements/S10_lean_cas_bridge_check.py
python3 ../_measurements/S10_lean_cas_verify.py
```

Regeneration uses the same script without `--check`. The check mode compares
the generated sources and manifest byte for byte without writing them.
The full build now passes **1962 selected axiom audits**, including 1712 in the
CAS bridge, across 115 canonical Lean source files and 3834 build jobs.
The exceptional count extension adds 316 audits to the earlier 788-audit
rerun build; the generic count extension adds another 205 and the root-list
extension adds 137. The coincidence extension adds 182, and the aggregate/metadata
extension adds 334. All 112 matrix/basis count records and all 27 root-list/count
records have their semantic bindings included in the axiom audit, along with
the three empty discarded-root list records and every selected coincidence
record's semantic binding. All 48 metadata records and nine associated arithmetic
records also have their semantic bindings audited. Thirty declarations require
no axioms; the others use only `propext`, `Classical.choice` and `Quot.sound`.
Warnings are errors, and there are no proof admissions.
[CAS_BRIDGE_VERIFICATION.txt](CAS_BRIDGE_VERIFICATION.txt) records the current
build, axiom audits and source hashes. The earlier verification files describe
the [committed 250-audit checkpoint](../CHECKPOINT.md).

The [regression instrument](../../_measurements/S10_lean_cas_bridge_check.py)
compares all 916 cells against the independently implemented existing parsers,
tests malformed inputs, and compiles isolated Lean mutations. It changes a
dimensionless coefficient in each engine's matrix, removes a chart assumption,
changes a basis normalization and assigns a wrong slot unit. Positive controls
distinguish intended proof failures from a broken test environment. The
[check record](../../_measurements/S10_lean_cas_bridge_checks.json) retains
outcomes and verifies canonical hashes before and after the tests.
The extended checks also cover locus guards, targeted points, row/column
selection errors, omission of either exceptional branch, incomplete/dependent
parallel bases, incorrect coordinate scaling, altered rank/nullity records,
the N3 input dimension, signed residuals and extension of a generic count
beyond its chart. Root-list checks cover omitted roots, incorrect algebraic
multiplicity, candidate/distinct-count confusion and the raw-tree frequency
filter. Coincidence tests cover guarded branches, both allowed half-axes,
Boolean decisions, actual witnesses and the orientation of root differences.
Metadata tests cover repeated fields, aggregate pair indices and equations,
root signs, the nonzero wavevector premise and skipped-origin decisions.
All 916 arithmetic comparisons, 56 complete coincidence payload comparisons,
48 complete metadata payload comparisons,
324 locus sign-pattern comparisons,
four point comparisons, three empty filter-list comparisons, 80 malformed-input
rejections, 36 Lean mutation rejections and 14 positive Lean controls pass.
Exact outcomes are retained in the check record.

## Remaining integration

Extend the bridge to reality-filter traces and remaining Q5/Q6/period-average
metadata, then D=4 and
other action packages. Q7 production construction, broad comparator
and export integration, and ledger/paper reconciliation remain open in
[COVERAGE.md](COVERAGE.md). This increment changes no CAS engine, frozen S10
export or S11 operand.
