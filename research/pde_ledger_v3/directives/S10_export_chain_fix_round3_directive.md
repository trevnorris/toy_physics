# S9/S10 export chain — round-three repair directive

## Authority and scope

Apply the round-three export-chain repair to the S9 and S10 SymPy engines. Regenerate S9 before S10 so
the downstream engine imports the repaired upstream module. Change only the two SymPy engines,
`scripts/extract_knob_inventory.py`, both generated export modules, both literal stdout captures, and this
directive. Do not change a Wolfram engine, step record, TeX artifact, action, ansatz, derivation, or
computed physics value.

The ledger is an object boundary. One physical or structural quantity has one printed name and one SymPy
object throughout the S9-to-S10 chain; authored words remain typed data rather than algebraic variables.

## Canonical shared objects

Use `omegaSquared` and `lambdaScale` as the shared names. This choice does not defer to the newer SymPy
engine: the independently written Wolfram engine already uses those names. Bind both names in S10 by
importing their S9 ledger records instead of constructing replacements.

The scale's assumptions are not an open choice here. S9 declares it positive, while both S10 engines put
the same object under an explicit positivity assumption. Preserve that positive object when consolidating
the binding; do not weaken it to the newer SymPy declaration's merely-real assumption.

Use the upstream `dim_<coefficient>_<axis>` alphabet for coefficient-dimension solver unknowns in both
steps. The Wolfram engine independently uses the same dimension-prefix and axis-suffix structure. Do not
retain S10's newer `<coefficient>_dim_<axis-word>` aliases, and do not publish aliases for any of these
shared quantities. A consumer substitution through the one ledger binding must reach every stored record
that contains the quantity.

## Corroborated overwrite provenance

Keep `step` as the record's last-setting step. When S10 independently re-derives an S9 record and the
overwrite operands, classes, and provenance guards agree, also write a `corroborated_steps` tuple naming
both computing steps into that ledger record. Do not put this field on a merely inherited record.

The overwrite authorization remains authored metadata on `ExportRecord`. The emitted overwrite rows can
check that authorization, equality, class preservation, and the S9-to-S10 transition; no programmatic
guard can establish that the downstream operand was derived rather than copied from S9. Record this limit
in the implementation evidence rather than presenting the flag as proof of independent derivation.

## Whole-tree traversability and immutable values

Normalize Python mappings, sequences, sets, and strings at the export boundary. Convert every matrix
handed to the ledger to an immutable matrix. Then inspect the entire stored tree, including every argument
of a SymPy object and every matrix entry. Reject a raw Python container or a mutable matrix found inside a
SymPy object; accepting the top-level type is insufficient.

The generated ledger must expose no value that a consumer can edit in place. Preserve exact reconstruction
comparison after normalization, so container mutability changes are explicit export-boundary changes.

## Authored words

Represent status markers, branch dispositions, route names, solver classifications, field names, class
labels, production-partition member names, run-pair names, emission-name inventories, and other authored
words with `Str`. This applies to unconditional readouts as well as branch-selected statuses. Do not use a
spelling list to keep these words out of the algebraic-symbol population; their SymPy type is the
separation.

## Discoverable dimensions

Do not restore an inline copy of a dimension object. For every export record that has dimension metadata,
write a `dimension_key` field whose value is the ledger key of the record carrying that dimension object.
A consumer obtains the units as `LEDGER[record["dimension_key"]]["value"]`; it neither recreates the
object nor guesses a sibling spelling.

Resolve the link from exact reconstructed object identity among the records being exported, preferring
dimension-named records structurally and never maintaining a per-object spelling table. Preserve inherited
links unchanged through the S10 merge. Remove the former computed/absent-dimension counters entirely; the
links themselves are the inventory.

## Validate before publishing

Construct and reconstruct each generated module in memory, emit the existing operands and residuals, and
run all authorized guards before writing the module path. Remove the prior generated module at engine
startup so any failure, including a failure before the writer is reached, leaves no rejected or stale
module behind. Do not add another in-run writer check beyond repairing the whole-tree symbol/traversability
check already in scope.

## Acceptance evidence

Keep every ablation program and its literal stdout at a named absolute path outside the repository. For
each repair, make the weakest mutation that should be rejected and record the result:

- reintroduce one alternate shared name and show that binding the canonical object no longer reaches the
  mutated record;
- remove corroboration metadata from an authorized equal overwrite and show the ledger audit rejects it,
  while separately recording that copied derivations cannot be detected;
- hide a Python list inside a SymPy object's argument tree and show the whole-tree traversal rejects it;
- turn one authored `Str` back into a `Symbol` and show it enters the bindable-symbol population;
- remove one dimension link and show a consumer-side discoverability audit rejects it;
- move S9 publication before its writer guards, force a guard failure, and show the rejected module is
  left behind only by the weakened copy;
- turn one immutable matrix back into a mutable matrix and show an attempted edit succeeds only on the
  weakened copy.

Apply a separate FORM ablation that changes the structure of a load-bearing action. Compare exported
values from the control and ablated copies to demonstrate that the export path carries the moved objects.
Do not use a coefficient rescale as this control.

Run S9 and capture literal stdout, then run S10 under `timeout 1800` and capture literal stdout. Preserve
complete source, generated-module, and transcript diffs against the pre-round working tree. Verify every
S9 record not explicitly overwritten is field-for-field identical in S10. Report whether any transcript
change extends beyond names, containers, and export bookkeeping. Emit counts and partitions for readers;
do not put their evaluated values into this directive or an acceptance target.
