# S9/S10 export chain — round-two repair directive

## Authority and scope

Apply the round-two export repair to the S9 and S10 SymPy engines. Regenerate S9 before running S10 so
the downstream engine imports the repaired upstream module. Change only the two SymPy engines,
`scripts/extract_knob_inventory.py`, the two generated export modules, the two literal stdout captures,
and this directive. Do not change a Wolfram engine, step record, or TeX artifact, and do not change either
engine's action, ansatz, derivation, or computed physics.

The generated modules are a public object boundary. A consumer must be able to bind every algebraic
symbol found in a stored value from a ledger record instead of recreating it by printed name.

## Structural symbol boundary

Use the existing declaration annotations as the source of class metadata. Extend the declaration-class
extractor so annotated assignments inside construction scopes are available to the export path without
making unannotated local working variables into inventory declarations. Collect symbols directly held by
an annotated declaration, including nested homogeneous symbol collections, without walking through a
constructed expression and reclassifying its dependencies.

Authored field names in structured export rows are string atoms, not algebraic symbols. Separate these
populations by SymPy type. Do not enumerate field-name spellings, symbolic variable spellings, or a list
of names to exclude. Build symbol records from the symbol atoms actually present in the records that will
be written. Prefer a symbol's declaration class; for a symbol that is itself a computed symbolic status,
derive its class from the record population that carries it. Preserve an upstream binding when the same
symbol reaches S10 from S9.

Every stored symbol name receives a ledger record whose value is that symbol. The key, assumptions, and
class therefore travel with the object. Symbols absent from stored values do not become records merely
because an engine constructed them on a control or eliminated path.

## Overwrite provenance

An S10 construction explicitly marks whether it is an upstream re-derivation. For each key collision,
emit the upstream and downstream values, their exact residual, both class tags and their residual, both
step tags, the overwrite authorization, and its residual. Guard after emitting the complete collection.
An unmarked collision fails even when its value is equal. A marked overwrite must preserve class and must
describe the S9-to-S10 provenance transition; only its last-setting step changes.

Imported records that do not pass through a marked overwrite remain exact copies of their S9 records.
Do not add a separate in-run preservation check; the generated-module diff establishes that property.

## One representation for dimensions

Keep named standalone dimension records and remove `dim` metadata from records. Do this at the record
writer, not with a list of dimension keys. A consumer retains directly addressable dimension objects but
loses the inline association between a value record and its dimension. The material-symbol records also
lose their metadata-only dimensions because those objects have no separate named record at this boundary.
This is an intentional information loss chosen to prevent downstream comparisons between two exports of
one constructed object.

## Traversable stored values

Every stored `value` is a SymPy object supporting `free_symbols` and `atoms()`. Convert Python container
trees generically to SymPy tuples, dictionaries, sets, and string atoms at the export boundary; do not
special-case a ledger key. In addition, make each dimension solver construct its returned solution
collection from SymPy containers so the live solve result remains traversable when nested in another
SymPy object.

Exact reconstruction continues to compare Python type and `srepr` after this normalization. Display text
is derived from the normalized stored object.

## Written-module symbol check

After reconstructing the complete generated module, walk every stored value. For each symbol name, emit
the ledger binding and the distinct symbol objects found under that name, followed by a residual. The
single guard rejects a missing or wrong ledger binding and rejects more than one assumption-bearing object
under one printed name. This check reads the written module; it is the only new writer check authorized by
this repair.

Do not add any other writer check. Retain the existing round-trip, population, class-tally, production
partition, and overwrite checks, adjusted only for the new record representation.

## Acceptance evidence

Keep every ablation program and its literal stdout at a named absolute path outside the repository. For
each item, exercise the weakest implementation that should be rejected:

- Replace the source-driven symbol population with a partial population and show that a consumer audit
  reports the missing binding, even though the weakened writer can validate what it chose to write.
- Change one indexed field-name atom back to an algebraic symbol and show that it enters the symbol
  population; demonstrate that no field-name spelling table participates in the separation.
- Remove overwrite authorization while colliding with an equal-valued upstream record, and separately
  alter only its class and only its source step. Retain the emitted operands, residuals, and failing exits.
- Restore one duplicate `dim` field and show an external object-duplication audit rejects it. Separately
  remove a standalone dimension record and show the same audit identifies the lost representation.
- Inject a same-named symbol with different assumptions into a stored value and show the written-module
  symbol guard rejects the run.
- Revert only the solver result to a Python list and show the traversal audit fails. Also place a Python
  container in a different export value and show generic normalization, rather than a key-specific rule,
  makes the stored value traversable.
- Apply an action-form ablation that changes the stiffness structure. Preserve the literal exported-value
  diff proving that the export path still carries moved physics objects.

Run the unmodified S9 engine and capture its literal stdout, then run S10 under `timeout 1800` and capture
its literal stdout. Compare every S9 record that S10 does not explicitly overwrite field for field across
the generated modules. Preserve complete diffs of both generated modules and both stdout captures against
the pre-repair working tree, and report the added and removed line totals for this round.

Do not write an evaluated count, tally, or partition size into this directive or a source comment. Emit
computed readouts and leave their interpretation to the step record.
