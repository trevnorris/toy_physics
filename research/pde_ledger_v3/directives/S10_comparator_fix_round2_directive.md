# S10 cross-engine comparator fix round 2 directive

## Authority and boundary

Repair only the existing cross-engine comparator and regenerate only its comparator-side stdout records
from the committed PY and WL transcripts. The transcripts are read-only inputs. Do not modify or run either
engine, do not start `wolframscript`, and do not modify any Wolfram file, step record, or TeX file. Add no
configuration, runner, authored name correspondence, expected tally, or acceptance population.

The live repository changes for this round are limited to:

```text
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_comparator_fix_round2_directive.md
```

Every tally remains a computed emission. Neither this directive nor the comparator may state its measured
value as an expected result.

## Preserve the differentiated object

The common derivative representation must contain the field name, the field's dependence set, and the
differentiated variable/order pairs. Canonicalize the dependence set without preserving argument slots,
because the engines may order the same field arguments differently. Keep differentiated variables paired
with their own orders, so changing the differentiation variable remains observable.

A dependency added, removed, or replaced must move the derivative residual. A permutation of the same
dependency set must not move it. The field name alone is not an adequate representation of the
differentiated object.

## Parse payload symbols as symbols

When a bare identifier in a SymPy transcript collides with a SymPy global object, supply an explicit local
symbol for that identifier. Determine whether an identifier is bare from its token context: a callable
head and an attribute namespace remain syntax, not payload symbols. Thus predicate forms and constructors
continue to parse through their callable objects while a standalone colliding spelling remains a free
symbol. The transcript spelling `pi` remains the genuine constant π.

This is parser behavior, not a cross-engine rename. Do not add a name-pair table or normalize one payload
symbol into a differently named symbol.

## Preserve mapping shape and order

Delete the dictionary branch from the normalizer. Do not sort mappings, turn them into rule sequences, or
replace that branch with another default unordered comparison. If a future object requires mapping or set
semantics, that correspondence requires its own stated decision.

## Make shared-name accounting total

Partition shared names into ordinary comparison names and shared duplicate names. A duplicate standard
name receives one duplicate verdict even when either or both transcripts repeat it; print every involved
operand under that single verdict. Continue to fail on every duplicate row, including a duplicate name
that is present on only one side.

Emit the number of shared-name verdicts and the residual between that number and the shared-name
population. The shared verdict population consists of parsed comparisons, unparsed comparisons, and
shared duplicate-name verdicts. Include the accounting residual in the final guard. Compute and emit all
of these quantities; encode no target population.

## Compute every disagreement kind

A failing nullspace-basis comparison obtains its disagreement kind through the same classifier used by
other parsed comparisons. Do not type a content label into the nullspace comparison path. Classification
remains diagnostic: it never changes the complete residual or its row guard.

The nullspace span residual and the matrix-membership residuals are distinct and load-bearing. Moving an
emitted span out of its paired matrix kernel moves membership, while loss of basis rank can move only the
span comparison. Retain both. Neither may be simplified away merely because the other detects some of the
same mutations.

## Acceptance and saved evidence

Use constructed in-memory records and disposable transcripts outside the repository. Save every ablation
script and its literal stdout at named absolute paths outside the repository. Each check must reject the
corresponding weaker repair:

- For derivative dependence, show that the former name-and-order representation accepts a changed
  dependence set while the repaired comparator rejects it. Also show that dependency permutation remains
  accepted and changing the differentiated variable is rejected.
- For symbol capture, show an input where the former no-locals parse collapses a free payload symbol onto a
  CAS object and therefore agrees, while the repaired comparator produces a nonzero residual. Also show
  that π and callable namespace syntax retain their intended parses.
- For mapping normalization, show that the former branch sorts and changes the object shape while the
  repaired normalizer returns the mapping with its insertion order intact.
- For accounting, construct a shared name duplicated by both transcripts. Show that row-based duplicate
  arithmetic does not equal the shared-name population, while the repaired run prints one verdict for the
  duplicated name and a zero accounting residual.
- For nullspace classification, replace the common classifier with a probe result during a failing
  nullspace comparison. The comparison must carry the probe result; a typed label fails this check.
- As a regression control, show that the complete nullspace residual fires on a wrong subspace and remains
  silent under a general invertible change of basis.

Regenerate the comparator and existing repoint-ablation stdout without running either engine. Save the
complete before/after diff for every regenerated `.out` outside the repository, and account for every
changed line. Literal stdout is the measurement; do not transcribe its measured tallies into this
directive.
