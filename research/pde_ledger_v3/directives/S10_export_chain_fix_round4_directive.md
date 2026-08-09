# S9/S10 export chain — round-four repair directive

## Authority and scope

Apply the round-four record-integrity repair to the S9 and S10 SymPy export chain. Regenerate S9 before
S10. Change only the two SymPy engines, `scripts/extract_knob_inventory.py`, both generated export
modules, both literal stdout captures, and this directive. Do not change or run a Wolfram engine, and do
not change a step record, TeX artifact, action, ansatz, derivation, or computed physics value.

This repair closes ways that an otherwise successful ledger can misstate its construction or the kind of
object it carries. It adds measurements at the generated-module and record boundaries; it does not add a
new physics conclusion.

## Failed-attempt invalidation and build inputs

Invalidate each generated-module path before any fallible project or upstream import. In particular, an
S10 attempt that cannot import S9 must already have removed the previous S10 product, and an S9 attempt
that cannot import its shared extractor must already have removed the previous S9 product. Keep the
existing validate-before-write order for failures reached after import.

Each generated module publishes an immutable mapping from every source input it consumed to the digest of
the exact bytes consumed. The S9 inputs include its engine and the shared extractor. The S10 inputs also
include the generated S9 module it imported. Reconstruct and compare this mapping before publication.
Thus invalidation identifies the most recent failed attempt, while the input mapping identifies the run
that produced a surviving module.

## Whole-value kind

Give every ledger record a structural `value_kind`. Determine it from the normalized whole value: a
top-level SymPy `Str` is an `AUTHORED_WORD`; every other stored SymPy object is a `COMPUTED_OBJECT`. Apply
the same rule to Python text after export-boundary normalization has converted it to `Str`. Do not use an
allowed-word list or infer the kind from record names, class tags, or spellings.

A `Str` nested inside a row does not make the whole row an authored word. Branch outcomes whose operands
are emitted alongside them remain honest authored words and must not be deleted. Emit the complete
population of whole-value authored-word records and its measured size. Include `value_kind` in exact
round-trip reconstruction, so changing only the kind of a reconstructed record is rejected.

## Assumption-channel comparison

For the MAIN records exported by each engine, collect the unary `Q.*` predicates on symbols from the
assumption objects under which the engine refined its results. Resolve each predicate's symbol through the
reconstructed ledger binding and ask the same predicate using only the assumptions intrinsic to that
exported symbol. Emit the ledger binding reconstruction, the reasoning predicate, the bound predicate,
the inherited result, and a per-row residual, followed by the total residual. Guard the total.

The comparison is source-driven over predicate objects. Do not maintain a symbol/predicate pair table.
Predicates on compound expressions are outside this particular two-channel comparison because a consumer
does not inherit them from one exported `Symbol`; keep emitting those premises in their existing records.
Control packages that do not enter the MAIN export population likewise remain outside this consumer
boundary.

## Dimension-link coverage

Keep `dimension_key` as the single discoverability mechanism for exported units. Emit every owner-to-
dimension-record link and the measured size of that population. This replaces the deleted counter over a
field that no longer exists: the readout now measures the live mechanism and changes when a link changes.
Do not restore the former computed/absent-dimension counters or write an expected coverage size into
source or prose.

## Immutable mapping boundary

Publish `LEDGER` as a read-only mapping whose record mappings are also read-only. Continue normalizing all
stored values to immutable SymPy objects. Emit the outer and per-record mapping operands and a residual,
and guard the residual. A consumer must be unable to replace a ledger entry or edit a record's metadata
through the public mapping.

## Name-split limit

The present S9/S10 ledger is clean of quantity-to-name splits by inspection, not by
`symbol_binding_residual`. That residual counts object variants independently under each printed name and
cannot detect one quantity published under two different names. S11 adding names is therefore unpoliced
by this instrument. Do not add a hand-maintained name-to-quantity map or claim that the existing residual
closes this gap.

## Acceptance evidence

Keep each ablation program, its literal stdout, and any child-engine stdout at a named absolute path
outside the repository. Exercise the weakest implementation that each property must reject:

- move invalidation back below the upstream import and show that an import failure preserves a stale
  marker only in the weakened copy; remove build-input metadata and show that a consumer cannot validate
  the weakened module;
- mark one top-level `Str` as a computed object and show the structural record audit rejects it without a
  spelling list, while nested `Str` values remain outside the whole-value population;
- widen an exported symbol's intrinsic assumption while leaving the engine's `Q.*` predicate unchanged,
  verify the mutation remains in the copied source, and show the assumption-channel residual rejects the
  run;
- remove one live `dimension_key` and show the emitted live-link population changes and a consumer audit
  rejects the weakened record;
- replace the public outer and inner read-only mappings with ordinary dictionaries and show mutation
  succeeds only in the weakened copy.

Apply a separate FORM ablation that changes the structure of the load-bearing main action rather than
rescaling a coefficient. Verify the source mutation remains present, run the ablated engine, and compare
its generated values with the control module to demonstrate that the export chain carries the changed
objects.

Run S9 and capture literal stdout. Then run S10 under `timeout 1800` and capture literal stdout. Preserve
complete diffs of both engines, the extractor, generated modules, and transcripts against the pre-round
tree. Verify every S9 record not overwritten by S10 is field-for-field identical across the generated
modules. Report from the diffs whether any transcript change extends beyond record/provenance bookkeeping,
without placing an evaluated count, tally, partition size, digest, or expected physics value in this
directive.
