# S10 cross-engine name-join comparator directive

## Authority and boundary

Build one comparator that reads exactly two `.out` paths supplied on its command line. One transcript must
contain `PY_` tags and the other `WL_` tags. The comparator has no config, YAML, manifest, per-step pair
file, runner, or authored name-to-name correspondence table. It does not import either engine and does not
run either CAS.

The committed inputs for this build are:

```text
research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out
research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out
```

Do not change those transcripts or either engine. The build adds the comparator, its controlled ablation,
their literal stdout, and this directive only.

## Transcript reader and join

Every nonempty transcript line must have one engine prefix, one emitted name, `: `, and one payload. A
malformed line and a duplicate emitted name are uncomparable findings. Print the offending line or all
duplicate operands, the reason no residual was computed, and a failing guard. Never select one duplicate
or discard an unparsed line.

Strip only `PY_` or `WL_` from each tag and join by exact equality of the remaining string. Names present
on only one side are inventory: print every name and the computed count for each side, but do not fail the
comparison for absence alone. Do not use similarity, payload equality, question numbering, or a special
case to turn two unequal emitted names into a joined row.

## Payload correspondences

Parse both payload syntaxes into a common exact symbolic representation. Apply snake_case to lowerCamel
mechanically to SymPy payload identifiers. The conversion has no exceptions. Normalize the engines'
power notation and ordered sequence delimiters mechanically; preserve sequence order and shape. Treat a
SymPy one-column matrix as the corresponding ordered sequence used by Wolfram. Parse the two derivative
surface syntaxes into the derivative's named variable/order pairs.

Do not rename a nonmechanical symbol to make a row pass. If a failed row becomes identical under a
bijective symbol-only substitution discovered from that row's structure, print the SymPy spelling, its
mechanical spelling, the Wolfram spelling, and every shared name that supplies evidence for the pair.
This is a D12 work-list entry, not a comparator correspondence: compute and guard the original residual.

For scalar expressions compute the exact simplified difference. For ordered structures compute a
componentwise residual. Preserve relational operators and operands; a type, shape, length, operator, or
non-subtractable structural difference is an explicit nonzero residual, not a parse success and not an
agreement. Print the raw SymPy operand, raw Wolfram operand, computed residual, and then the row guard for
every parsed shared name.

If either payload cannot be parsed or normalized, print both raw operands and the exception reason in a
separate uncomparable category, print `NOT_COMPUTED` as the residual, and fail the row. An unparsed shared
row never contributes to agreement.

## Noncanonical representations of canonical objects

A direct nullspace basis is representation, but its span is a canonical object. For every shared name
ending in `N6_NULLSPACE_BASIS`, interpret the emitted vectors as rows and reduce them to a canonical
row-space representative. Compare those representatives componentwise. Also obtain each engine's
same-root emitted `N1_MATRIX` and compute its action on the transpose of that engine's canonical span.
Print the raw operands, canonical spans, span residual, both matrix-action residuals, the complete row
residual, and the guard.

This residual is invariant under permutations, nonzero rescalings, and general invertible basis changes,
but remains sensitive to a changed subspace, rank loss, and vectors outside the emitted matrix's
nullspace. A missing, duplicate, malformed, or shape-incompatible basis or paired matrix is uncomparable
and failing. No non-canonical object is excluded merely because direct representation equality is not
meaningful.

## Final guard and emitted record

After all row records, inventories, format findings, and D12 symbol work-list entries, emit the measured
tallies. Decompose disagreements into naming-only, representational, and content counts without changing
the original residual or guard. A naming-only row is one made identical by a row-discovered bijective
symbol substitution; route-tag rows ending in `_ROUTE` or `_ROUTES` are content and never symbol-spelling
evidence. A representational row is one made identical by mechanically converting equality-to-zero leaves
to left-minus-right expressions, either directly or after the same row-discovered symbol substitution;
classify this before naming-only. All other disagreements are content divergences. Also decompose agreements
into bare integers, empty containers, and symbolic or structured objects.

The process exits nonzero if any shared row disagrees, any shared payload is uncomparable, any emitted name
is duplicated, or any transcript line is malformed. One-sided inventory does not make the process fail.

Do not encode an expected tally, agreement total, or acceptance count in the script or this directive.
The committed stdout is the measurement.

## Controlled re-point ablation

The saved ablation operates on the parsed transcript records in memory and leaves both committed inputs
untouched. Re-point the SymPy standard name `S10_MAIN_D3_Q2_MATRIX_A` to the neighbouring live payload
emitted as `S10_MAIN_D3_Q2_MATRIX_B`. Compare the target row before and after re-pointing against the same
Wolfram target. Print both complete comparator row records, whether the residual moved, whether the row
guard moved, and an ablation guard. The ablation succeeds only when the original row passes and the
re-pointed row fails with a different residual.

Generate the committed records without starting either engine:

```bash
python research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py \
  research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out \
  research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out \
  > research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out

python research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py \
  research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out \
  research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out \
  > research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
```
