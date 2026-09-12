# S10 D3 emitted coincidence loci, decisions and witnesses

This extension connects the generic and targeted anisotropic D3 coincidence
records from both focused CAS transcripts to the proved root classification.
It adds 13 arithmetic records containing 28 scalar expressions and 56
logical/container records containing 84 translated predicates. Seven Wolfram
operand records occur in both categories: the selection contains 62 distinct
transcript tags.

The arithmetic selection contains all 14 emitted root differences across the
six engine/case combinations, plus all 14 SymPy pair-index entries. Every
scalar has checked units and an equality to its reference operand. The logical
selection covers the primary loci, allowed operands/regions, Boolean decisions,
Wolfram decided-empty/nonempty outcomes, and SymPy witness and branch records.
Exact source payloads, line numbers, hashes, trees and semantic theorem names
are retained in the [manifest](../../_measurements/S10_lean_cas_bridge_manifest.json).

## Full loci on the stated domain

[CoincidenceSupport.lean](S10Audit/CAS/CoincidenceSupport.lean) proves the
classification at arbitrary real wavevectors, under

```
rho > 0, mu > 0, sigma > 0, sigma != 1, lambdaScale > 0.
```

The dimension is specialized to D=3. Coordinate and amplitude realness are
represented by their Lean types. No generic basis chart is assumed: the
parallel axis is included in the proof's domain. Using the transcript's
one-based coordinate names:

| Root pair | Complete coincidence locus | Allowed part, with k != 0 |
|---|---|---|
| Static / ordinary | k1 = k2 = k3 = 0 | Empty |
| Static / extra | k1 = k2 = k3 = 0 | Empty |
| Ordinary / extra | k2 = k3 = 0 | k1 != 0, k2 = k3 = 0 |

[CoincidenceBindings.lean](S10Audit/CAS/CoincidenceBindings.lean) proves each
translated locus equivalent to this geometry. It also proves that each actual
emitted root-difference expression vanishes if and only if its printed locus
holds. Each allowed-region predicate has the same complete geometry, and the
printed Boolean and decided-empty/nonempty records agree with existence in
that printed allowed region. Cross-engine theorems identify the generic loci,
allowed regions and Boolean decisions pair by pair.

The targeted reruns use each engine's own previously certified parallel and
perpendicular point. The parallel rerun has two distinct roots; its sole pair
cannot coincide on the coefficient domain. All three distinct-root pairs at
each perpendicular target are likewise separated. Wolfram's parameter loci,
including `mu = 0` and `sigma = 1` branches, are retained and proved inadmissible
under the declared coefficients. These are statements about the distinct-root
lists; the parallel nonzero root still has algebraic multiplicity two, as
proved in [ROOT_RESULT.md](ROOT_RESULT.md).

## Conditional branches and witnesses

Wolfram's static/extra locus contains four square-root branches guarded by
`sigma < 0`, followed by the positive-sigma origin branch. The translation
preserves every branch, square-root expression and conditional guard. Lean
eliminates the negative-sigma branches using the stated domain. Square roots
are represented by `Real.sqrt` only inside those excluded branches; the proof
uses their false guards and makes no claim about Wolfram's complex-square-root
semantics or those loci at negative sigma. The allowed ordinary/extra region retains
both `k1 < 0` and `k1 > 0` banks and both admissible sigma intervals.

The SymPy witness `(1,0,0)` has a proof that it satisfies its translated witness
predicate and lies in the allowed coincidence locus. All emitted empty witness
lists and branch groups retain their empty interpretation. The witness is an
existence certificate; the universal locus equivalence proves coverage of the
whole allowed axis. The result is not inferred from sampling that point.

## Translation and verification

[S10_lean_cas_coincidence.py](../../scripts/S10_lean_cas_coincidence.py) uses a
restricted, non-evaluating arithmetic/logical grammar. It preserves unions,
conjunctions, inequalities, equality rules and conditional expressions.
Dimension membership is specialized only for D=3; unsupported integer
membership is rejected. Denominators are restricted to nonzero constants and
the declared nonzero coefficient factors. Unsupported functions, wrong
arities/operators, oversized powers, zero numeric denominators, coordinate
denominators, trailing syntax and unknown symbols are rejected.

The regression instrument independently parses all 56 complete logical
payloads. It uses the existing SymPy parser and SymPy's Mathematica tokenizer
and full-form parser. Wolfram membership expressions remain opaque Boolean
atoms in this comparison because the older comparator's final Boolean
conversion does not accept them inside conjunctions. No production comparator
code is changed. The independent arithmetic comparison additionally covers
all 28 new scalar expressions.

New Lean mutation checks reverse a root difference, admit the negative-sigma
branches, remove the positive half-axis, flip the nonempty decision, replace
the witness with zero wavevector, and omit an axis constraint. Three positive
controls exercise the guarded locus, both allowed banks and the actual witness.
The full build passes **1962 selected axiom audits**, including **1712 CAS
audits**, across **115 canonical Lean files** and 3834 build jobs. Only
`propext`, `Classical.choice` and `Quot.sound` occur. Warnings are errors and
there are no proof admissions. The complete bridge now contains 289 arithmetic
records and 916 scalar expressions. All 916 arithmetic comparisons, 56 complete
coincidence payload comparisons, 48 complete metadata payload comparisons,
324 coordinate-locus sign patterns, four point
checks and three filter-list checks pass. All 80 malformed inputs and 36
mathematical mutations are rejected; all 14 positive controls pass. Canonical
source, input and frozen export hashes remain unchanged by the tests.

Full-build and regression evidence is in
[CAS_BRIDGE_VERIFICATION.txt](CAS_BRIDGE_VERIFICATION.txt) and the
[check record](../../_measurements/S10_lean_cas_bridge_checks.json).

## Remaining S10 integration

This certifies the selected primary records on the explicit coefficient
domain. It does not prove that the CAS solvers discover all strata for arbitrary
inputs. The [metadata extension](METADATA_RESULT.md) now connects duplicated
Q8 fields, Wolfram aggregate associations, root sign labels, root solver-condition
lists, spectrum solve operands/statuses and retained/skipped-stratum records.
SymPy reality-filter traces and remaining Q5/Q6/period-average metadata still
need translation or reconciliation.
D4 and the other emitted action packages, Q7 production construction, the broad
comparator/export chain and paper alignment remain listed in
[COVERAGE.md](COVERAGE.md).

The parsers and transcript translation are tested software; the generated
mathematical claims are checked by Lean's kernel. No CAS engine, input
transcript, frozen export or S11 operand is changed by this extension.
