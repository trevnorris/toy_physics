# S10 cross-engine comparator fix round 1 directive

## Authority, inputs, and boundary

Repair the existing name-join comparator using the two committed S10 transcripts as read-only inputs. Do
not modify or re-run either engine, do not modify either transcript, and do not start `wolframscript`.
The live build is limited to:

```text
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_cross_engine_comparator_directive.md
research/pde_ledger_v3/directives/S10_comparator_fix_round1_directive.md
```

Keep the exact-name join, parsing rules, raw operands, original per-row residuals, one-sided inventories,
format and duplicate failures, and final guard. Add no config, name-pair table, per-step file, runner, or
expected count.

## Compare a nullspace basis as a subspace

A basis is non-canonical, but its span is an object. A shared name ending in `N6_NULLSPACE_BASIS` must be
compared; it must not be excluded. Interpret each emitted ordered collection of vectors as the rows of a
basis matrix `B` and compute its reduced row-echelon form `C = rref(B)`. This is the canonical row-space
representative and is unchanged by a permutation, a nonzero rescale, or a general invertible change of
basis.

For the same root and engine, obtain the already-emitted `N1_MATRIX` named by replacing the basis suffix
with `N1_MATRIX`. Print both raw basis operands, both canonical spans, the componentwise residual between
the spans, and each engine's matrix action `N1_MATRIX C^T`. The row residual contains the span residual and
both matrix-action residuals. Guard the row on that complete residual. A missing, duplicate, malformed, or
dimensionally incompatible paired matrix makes the basis row uncomparable and failing; never select an
ambiguous matrix.

This comparison must be insensitive to every invertible basis change while remaining sensitive to a
changed span, lost basis rank, or a vector that is not in the emitted matrix's nullspace. If this
construction cannot be computed for an emitted basis shape, stop with the operands and explicit reason;
do not restore an exclusion.

## Decompose disagreements without absorbing them

Every parsed shared row still prints and guards its original residual. A diagnostic classification does
not turn a failed row into an agreement. Print one classification on every disagreement and emit separate
summary counts:

- `REPRESENTATIONAL`: mechanically converting equality-to-zero leaves to their left-minus-right
  expressions makes the row residual zero, either directly or after the same row-discovered symbol-only
  substitution. Apply this classification first: a row with both surface-form and spelling differences is
  not naming-only.
- `NAMING_ONLY`: absent that equality surface-form difference, a bijective substitution of unmatched
  payload symbols makes the row residual zero. Discover this from the row; do not author a spelling pair.
  Permit an algebraic single-symbol check so expanded and factored forms supply the same evidence.
- `CONTENT`: neither diagnostic projection makes the residual zero.

An emitted route tag is content, not an algebraic symbol spelling. Rows whose standard name ends in
`_ROUTE` or `_ROUTES` do not supply symbol-rename candidates. In particular, never bind a route label to a
description merely because each parsed as one symbol. Keep every naming-only row in the D12 worklist with
its evidence names, and let the algebraic inference see symbol pairs that occur inside differently
factored expressions.

Also emit agreement-form counts for bare integers, empty containers, and symbolic or structured objects.
These are a decomposition of computed agreements, not expected targets.

## Recorded baseline facts, not regression targets

Before this repair, the comparator's `FINAL_GUARD` was already `FAIL` and its process exit code was already
`1`. Neither can show whether a value mutation was detected; the per-row residual is the ablation signal.
Do not use the final exit status as the mutation regression bar.

Of the 362 agreements reported by that baseline, 215 were bare integers, 14 were empty containers, and
133 were symbolic or structured. This records the support of the baseline agreement claim. It is not an
acceptance count for the repaired build.

## Acceptance and saved evidence

Regenerate only the comparator stdout and its existing repoint-ablation stdout from the committed
transcripts. Save each complete before/after `.out` diff outside the repository.

Use disposable, in-memory operand mutations and save every ablation script plus its literal stdout at
named absolute paths outside the repository. Demonstrate separately that the basis residual moves when:

1. a root-one basis is replaced by the neighbouring root-two basis;
2. one vector component is scaled without scaling the rest of its vector;
3. a root-one basis is replaced by a coordinate vector.

Demonstrate that the residual does not move under both a legitimate nonzero whole-vector rescale and a
legitimate general invertible basis change. Print the baseline and mutated residuals, whether each moved,
and an ablation guard. The ablation passes only if all three invalid mutations move the residual and both
legitimate changes leave it unchanged.
