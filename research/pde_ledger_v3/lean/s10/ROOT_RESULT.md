# S10 D3 root lists, multiplicities and filter counts

This extension connects the solution lists and distinct-root lists in both
focused anisotropic D3 transcripts to the actual emitted determinants. It
adds **27 arithmetic records containing 51 scalar expressions**, plus three
empty discarded-root list records. At this root-list checkpoint, the combined bridge covered **267 arithmetic
records and 834 scalar expressions**, 12 exceptional-locus predicates, four
targeted points and these three filter lists.

| Added emission | Records | Scalar expressions |
|---|---:|---:|
| Raw solution rules | 6 | 17 |
| Distinct-root lists | 6 | 16 |
| Distinct-root counts | 6 | 6 |
| Wolfram candidate counts before and after filtering | 6 | 6 |
| Wolfram candidate/distinct count associations | 3 | 6 |

## Complete spectrum and algebraic multiplicity

[RootSupport.lean](S10Audit/CAS/RootSupport.lean) defines the cubic polynomial

```
(rho^3 sigma / 8) z (z - ordinaryRoot) (z - extraRoot).
```

Lean proves that its evaluation equals the determinant of the reference
route-B matrix for every real squared-frequency parameter `z`. The coefficient
is nonzero when `rho != 0` and `sigma != 0`. Its root multiset is therefore
exactly `[0, ordinaryRoot, extraRoot]`, retaining repeated roots. Its total
algebraic multiplicity is three. These are polynomial-root multiplicities in
squared frequency; they are separate from the modal kernel dimensions proved
in the [count bridge](COUNT_RESULT.md).

[RootBindings.lean](S10Audit/CAS/RootBindings.lean) connects the imported
solution and distinct-root lists to this multiset. For each of the six
engine/case combinations, a squared frequency belongs to the imported
distinct-root list if and only if the **actual emitted determinant** vanishes.
The distinct lists have no duplicates, and each printed distinct-root count
equals the cardinality of the corresponding finite set.

| Case | SymPy solution entries | Wolfram solution entries | Distinct roots | Algebraic multiplicities |
|---|---:|---:|---:|---|
| Generic chart | 3 | 3 | 3 | Static 1, ordinary 1, extra 1 |
| Parallel target | 2 | 3 | 2 | Static 1, nonzero 2 |
| Perpendicular target | 3 | 3 | 3 | Static 1, ordinary 1, extra 1 |

The parallel case is deliberately asymmetric. SymPy prints the nonzero root
once; Wolfram prints it twice. Both lists contain exactly the determinant's
distinct zeros. Wolfram's candidate multiset also equals the polynomial's
root multiset, including the repeated entry. SymPy's two-entry parallel list
does not encode algebraic multiplicity; the separate polynomial theorem
proves the double root.

The semantic bindings assume `rho != 0`, `mu != 0`, `0 < sigma`, and
`sigma != 1`. Generic distinctness and counts additionally retain
`GenericChart sigma k`, including all three nonzero coordinate assumptions.
The exceptional cases use each engine's exact previously certified point,
without substituting coefficient values. Generic roots use physical squared
frequency; fixed-point roots retain the reduced convention
`omegaSquared = kappa^2 z` from the [rerun bridge](EXCEPTIONAL_RERUN_RESULT.md).

## The actual candidate filter

The Wolfram audit selects candidates with `FreeQ[expression, omegaSquared]`.
This is a syntactic test for an unresolved squared-frequency symbol, not a
positivity or reality test. The Lean `frequencyFree` function performs that
test on the **raw imported expression trees**. The generated proofs establish
that all candidates pass, that the rejected list is empty, and that the
printed before/after counts equal the lengths of the actual lists. Both
fields of each printed count association receive semantic bindings.

Filtering must precede the literal-zero unit adapter. A printed zero is
frequency-free, while its dimension-carrying representation may introduce a
frequency symbol. The manifest records the raw tree separately from the tree
used for dimensional checking. The regression suite tests this distinction,
and also checks that `z - z` remains syntactically dependent on `z` even though
its numerical value cancels.

## Parsing, provenance and verification

[S10_lean_cas_roots.py](../../scripts/S10_lean_cas_roots.py) accepts only
single-variable solution rules with the exact `omegaSquared` key and scalar
arithmetic values. It preserves every candidate and its order, including
duplicates. The association parser accepts exactly the two named count fields.
Unsupported conditions, extra assignments, wrong variables, duplicate fields,
fractional counts, malformed containers and trailing text are rejected.
Original payloads, line numbers, source hashes, arithmetic trees, units,
domains and semantic proof names remain in the
[manifest](../../_measurements/S10_lean_cas_bridge_manifest.json).

The independent comparison uses the existing SymPy and Wolfram parsers to
extract each solution value and association field separately. It also checks
all three empty filter lists. New Lean regressions omit a root from the cubic
multiset, collapse the parallel multiplicity to one, confuse candidate count
with distinct-root count, and filter a unit-adapted zero. Positive controls
prove the double root and the intended raw-tree filter behavior.

The root-list checkpoint build passed **1446 selected axiom audits**, including **1196 CAS
audits**, across **98 canonical Lean files** and 3817 build jobs. Only
`propext`, `Classical.choice` and `Quot.sound` occur; warnings are errors and
there are no proof admissions. All 27 new arithmetic records and all three
filter-list records have their semantic bindings included in the audit.

That checkpoint regression suite passed 834 independent arithmetic comparisons,
324 locus sign-pattern checks, four point checks and three filter-list checks.
It rejected all 53 malformed inputs and all 24 mathematical mutations; all eight
positive controls passed. Canonical source, input and export hashes remain
unchanged by the tests. The latest full-build and regression results are recorded in
[CAS_BRIDGE_VERIFICATION.txt](CAS_BRIDGE_VERIFICATION.txt) and
[the check record](../../_measurements/S10_lean_cas_bridge_checks.json).
Reproduction commands are in [CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md).

## Remaining work

The subsequent [coincidence extension](COINCIDENCE_RESULT.md) now certifies
the primary emitted coincidence equations, loci, allowed regions, decisions
and witnesses against this classification. The [metadata extension](METADATA_RESULT.md)
also connects Q8/aggregate fields, reported root signs, root-condition lists,
spectrum solve operands/statuses and retained/skipped-stratum records.
Reality-filter traces, further Q5/Q6/period-average metadata, D4 and other
packages remain outside these extensions. The original solver algorithms,
full comparator/export chain, Q7 production alignment and paper reconciliation
remain separate obligations in [COVERAGE.md](COVERAGE.md).

The parser and transcript translation remain tested software. The generated
mathematical claims are checked by Lean's kernel. This extension changes no
CAS engine, frozen export, input transcript or S11 operand.
