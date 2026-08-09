# S10 cross-engine comparator fix round 4 directive

## Authority and boundary

Repair only the existing cross-engine comparator and regenerate only its comparator-side stdout records
from the committed PY and WL transcripts. Treat both transcripts as read-only inputs. Do not modify or run
either engine, do not start `wolframscript`, and do not modify any Wolfram file, step record, or TeX file.
Add no configuration, runner, authored name correspondence, expected tally, or acceptance population.

The live repository changes for this round are limited to:

```text
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_comparator_fix_round4_directive.md
```

Every tally remains a computed emission. Neither this directive nor the comparator may state its measured
value as an expected result.

## Classify only after transcript refutation

Retain row-local symbol substitution as diagnostic evidence. It does not change the operands, residual, or
row guard. Build the candidate evidence and the transcript-wide distinction evidence before assigning the
authoritative disagreement kind. A symbol pair that the transcripts refute supplies no naming evidence to
that assignment: a row whose substitution depends on such a pair is content, not naming-only. Apply the
same exclusion when a projected equation-system comparison considers symbol substitution while deciding a
representational kind.

Continue to emit the rejected pair, its candidate rows, and its per-engine, per-spelling distinction
evidence. Continue to omit it from the rename worklist. A candidate pair for which the transcripts provide
no contradiction remains explicitly `UNREFUTED` and may still support a naming-only classification. That
status is not acceptance of a rename and creates no comparator correspondence.

## Delete non-measuring residuals

Delete `CONTENT_DECOMPOSITION_RESIDUAL` and `SHARED_ACCOUNTING_RESIDUAL`: their operands are constructed as
complete disjoint case splits, so neither residual can move. Remove both terms from `FINAL_GUARD` and do
not replace them. Keep emitting the content parent and sub-populations and the shared-name accounting
pieces so their arithmetic remains visible in the stdout record.

## Recorded limits; do not repair here

The PY payload path is protected against bare names resolving to SymPy callables, classes, registries, or
namespace objects: `_sympy_payload_locals` turns such a would-be shadow into a free payload symbol. The WL
payload path through `parse_mathematica` has no corresponding shadow protection. In particular, the bare
spellings `beta` and `gamma` resolve to library callables on that path. If a future shared row uses either
as a symbol, the comparison is loud but incomplete: the row emits `RESIDUAL: NOT_COMPUTED` and fails rather
than producing a false agreement. Do not change either parser in this round.

The content route-token and numeric-or-algebraic selectors inspect the emitted name suffix rather than the
residual. The committed input routes the selected rows consistently, but a future row can therefore carry
real algebra into the route population or route data into the algebra population. Record this limit and do
not change either selector in this round.

## Preserve earlier repairs

Keep the object-based constant/shadow distinction, mechanical symbol normalization, derivative dependence
sets, duplicate accounting, mapping-order handling, and computed disagreement kinds. Preserve both parts
of the nullspace comparison: canonical-span disagreement and each engine's emitted-matrix membership
residual. A symbolic invertible basis change must preserve the complete nullspace residual, while rank loss
and matrix perturbation remain distinguishable by the span and membership halves.

## Acceptance and saved evidence

Use constructed in-memory records and disposable transcript copies outside the repository. Save every
ablation script and its literal stdout at named absolute paths outside the repository.

- On the committed transcripts, alter one emitted matrix only in a disposable copy by re-pointing its shear
  stiffness symbol to the brane-density symbol. Show the pair is transcript-refuted, the affected row is
  content, the naming-only population loses the row, and the content population gains it. Show separately
  that a genuine unrefuted spelling difference remains naming-only.
- Exercise each deleted residual's former construction with duplicate names, a duplicate on only one side,
  and a nullspace-basis row whose supporting matrix is absent. From the emitted parent and component
  populations, compute the former residuals externally and show they cannot move. Confirm neither deleted
  label is emitted and neither term participates in the final guard.
- Re-run the nullspace residual on committed payloads with a wrong subspace and with a symbolic invertible
  change of basis. Demonstrate independently that rank loss moves the span half and an unchanged-span matrix
  perturbation moves the membership half.
- Re-run the constant/shadow split, derivative-dependence, mapping-order, duplicate-accounting, and computed
  disagreement-kind controls from the earlier rounds.

Regenerate the comparator and existing re-point ablation stdout without running either engine. Save the
complete before/after diff for every regenerated `.out` outside the repository and account for every
changed line. Literal stdout is the measurement; do not transcribe its measured tallies into this directive.
