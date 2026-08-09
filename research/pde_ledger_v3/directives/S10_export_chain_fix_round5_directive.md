# S9/S10 export chain — round five. A removal round.

## Authority and scope

Round five deletes two residuals from the S9 and S10 SymPy export chain, repairs one regression that
round four introduced, and records six limits that were measured and are being left in place. It adds
nothing: no new check, residual, guard, count or emitted tag exists after this round that did not exist
before it.

Changed: both SymPy engines, both generated export modules, both literal stdout captures, and this
directive. `scripts/extract_knob_inventory.py` was in scope and needed no change. No Wolfram engine, step
record, TeX artifact, action, ansatz, derivation or computed physics value was read, run or changed.

Every physics value in both transcripts is untouched. Both generated modules differ from their round-four
predecessors in exactly one line, and that line is the build-input digest map, which is a function of the
engine's own bytes and therefore moves whenever the engine is edited at all. Every ledger record in both
modules is field-identical to round four, and every S9 record that S10 does not restate is field-identical
between the two modules.

## R1 · Two residuals are deleted because neither could ever have been nonzero

**`value_kind` round trip.** Each engine wrote the `value_kind` field with `export_value_kind` and then
compared the written field against a fresh call to `export_value_kind` on the same object. Writer and
checker were one function, so the comparison was an identity on that function's output rather than a
measurement of anything about the record.

Deleted: the `value_kind` term in the export round-trip residual, in both engines. Kept: the `value_kind`
field itself, its classifier, and the emitted authored-word population. The field is a carrier a consumer
reads; the comparison was decoration on it.

**`BUILD_INPUT_DIGEST_RESIDUAL`.** Each engine built a `str → str` digest map, rendered it into the
generated module with `repr`, executed that module, and compared the executed map against the map it had
just rendered. The two operands are related by a `repr` → `eval` round trip on a mapping of plain strings,
so they are the same value by construction; the residual, its emission and its guard are deleted in both
engines. Kept: the digest map itself, in the generated module and in the emitted operand tag. The map is
load-bearing — it is how a stale module was detected — and it is only the self-comparison that was inert.

**How this was established: by running it, not by argument.** Four form ablations were run against the
round-four engines, each in a private tree built from byte-identical copies, and each preceded by a
control run confirming the tree reproduces the repository transcript byte for byte.

- Pin `export_value_kind` to each of its two constants in turn, in each engine. Every run exits zero, the
  export round-trip residual reads zero, and the module publishes. Where the pinned constant is not the
  one a record deserves, the ablation stdout enumerates the published records whose `value_kind` then
  contradicts the object it labels.
- Pin `file_sha256` to a constant, in each engine. The engine exits zero, the digest residual reads zero,
  the module publishes, and every published digest entry disagrees with the file on disk.

A residual that reads zero while the property it names is false in every row is not measuring that
property. Both are gone rather than repaired, because every repair available to them is another
comparison of the writer against itself.

**One thing the round measured without looking for it.** In S9 the `value_kind` field takes one value
across the entire ledger, so pinning the classifier to that same constant changes nothing at all; the
mislabelling only appears under the other constant. The field distinguishes records in S10 only. Recorded,
not acted on.

## R2 · Importing an engine now has no effect on any file

Round four moved each engine's generated-module `unlink` to module scope so that a failed run could not
leave a stale module behind. That made merely importing an engine delete the corresponding generated
module — for S10 silently, on an import that then succeeds.

Two properties have to hold at once:

- **P1** importing an engine, by any route, changes no file;
- **P2** the chain's own failure — an upstream input missing — leaves no stale generated module.

Both engines now guard the module-scope `unlink` with `__name__ == "__main__"`. Because the guarded
statement still stands above the upstream imports, a script run invalidates the previous product before
anything fallible is attempted, and an import performs no file operation at all.

**S9 needed a second change that R2 did not name.** S9's `unlink` is not its only module-scope file
effect: S9 calls `write_exports()` at module scope, so importing S9 rewrote the generated module even with
the `unlink` guarded. That call is now guarded the same way. S10 needed no second guard; its entry point
was already behind `__main__`.

**Measured, over four placements of the same statement**, each built from pristine sources, each measured
by snapshotting the whole tree on existence, content digest and modification time, so that a rewrite with
identical content still registers as an effect:

| placement | P1 import inert | P2 no stale module |
| --- | --- | --- |
| A · as built (module scope, unconditional) | no — deletes or rewrites | yes |
| B · weakest change (moved into the entry function) | no — S9 still rewrites on import | no |
| C− · guard the `unlink` only (S9) | no — S9 still rewrites on import | yes |
| C · guard the `unlink`, and S9's `write_exports()` call | yes | yes |

Two import routes were exercised for every placement: loading the engine file by path without its
directory on `sys.path`, which is how a scanner or doc generator loads a module and which makes the
engine's own upstream import fail; and importing by name with the directory on `sys.path`, which is how a
comparator sitting beside the engines would load it and which lets the import succeed. Under placement A
the second route is the dangerous one: the S10 module vanishes and the import reports success.

The shipped sources were then measured directly, both in a private tree and against the repository
directory itself: four imports, no file effect; and, in the tree, an upstream-missing run of each engine
leaves no stale module.

Placement B and placement C− are each a smaller change than the one shipped, and each fails a property the
shipped one holds. The content change B and C− produce on import is the generated module being rewritten;
its bytes differ from the baseline module only because a variant engine has a different digest of its own
source. The finding is the write, not the difference.

## R3 · Six limits, measured, recorded, and deliberately not fixed

Each of these is a real limit on an instrument. Every repair proposed for any of them is either a denylist
or another comparison of a writer against itself, which is what R1 just deleted twice. They are recorded
here so that a consumer knows what the chain does not establish.

1. **`value_kind` marks one carrier only.** A sentence can reach the ledger as a `Symbol` name, inside a
   `Function` head, inside a `Tuple` of `Str`, on either side of an `Eq` of `Str`, or in the `display`
   field — which is a raw `str` on every record and appears in no residual. A `Symbol` whose name is a
   sentence would acquire an auto-created record keyed by that sentence; that particular route is measured
   empty in both published modules today, but nothing in the chain closes it. The field marks a top-level
   `Str` and nothing else.
2. **The authored-word count does not bound the population.** The emitted authored-word record count, the
   number of records carrying authored text anywhere in the value, and the number of distinct authored
   strings are three different quantities in the S10 module, in that increasing order, and only the first
   is emitted. In the S9 module all three are empty, which is also why S9's `value_kind` field takes one
   value everywhere.
3. **The assumption-channel residual is coverage-free.** Emptying the entire `Q` channel scores zero and
   publishes, in both engines. Strengthening a `Symbol` that no `Q` predicate mentions scores zero and is
   inherited downstream. The channel reaches part of the symbols each ledger binds, from the `MAIN`
   packages only, and the residual says nothing about the rest.
4. **A parse failure leaves a stale module.** Python compiles a module before executing it, so a syntax
   error anywhere in an engine stops the guarded `unlink` from ever running and the previous generated
   module survives, in both engines. No module-scope statement inside the engine can prevent this. It is
   detectable from outside, via `BUILD_INPUT_DIGESTS`.
5. **Every export guard is an `assert`.** Under `PYTHONOPTIMIZE=1` a violated guard publishes: with the
   count residual forced nonzero, the engine exits zero and writes the module. The transcript still
   carries the nonzero residual, because the operands and the residual are emitted before the guard runs —
   which is the only reason this failure is visible at all.
6. **`BUILD_INPUT_DIGESTS` omits the CAS.** It names files only — no SymPy or Python version — so two
   modules with identical digests can hold different values if the CAS changed underneath them.

Limits 1, 2 and 6 were read off the published modules. Limits 3, 4 and 5 were established by running the
shipped engines under an ablation: emptying the `Q` channel entirely, appending a syntax error to an
engine, and forcing a guarded residual nonzero under each setting of `PYTHONOPTIMIZE`. None of them was
taken on report.

**What a consumer must therefore check for itself**, because the chain cannot check it from inside:

- that each digest in `BUILD_INPUT_DIGESTS` matches the file of that name on disk — this is the check
  deleted from the engines in R1, and it is only informative when performed by something that is not the
  writer (limits 4 and 6);
- that the CAS the consumer is running is the one the module was built under, since no field records it
  (limit 6);
- that the engines were run without `PYTHONOPTIMIZE`, since a guard is an `assert` (limit 5);
- that any authored text a consumer cares about is sought in `display`, in `Symbol` names, and inside
  composite values, and not only in the `value_kind` field or the authored-word population (limits 1
  and 2);
- that a symbol's assumptions are those the consumer needs, read off the exported `Symbol` itself, rather
  than inferred from the assumption-channel residual (limit 3).

## What this round did not build

Nothing here was blocked for want of a new check, so nothing was named and stopped on. Every property
above that remains unverified is unverified because verifying it from inside the writer is exactly the
defect this round removed; the list immediately above is where that work belongs, and it belongs to a
consumer, not to an engine.

## What moved in the transcripts

Both `.out` files differ from round four by one deleted line — the digest residual — and by the
build-input digest line, whose hex values change because the engine sources changed. Nothing else moved in
either transcript. The digest line will move on every future round that touches an engine at all; a
comparator diffing transcripts across rounds should expect it.

## Where the measurements are

Every ablation script and its literal stdout is outside the repository, under
`/tmp/claude-1000/builder_round5/`: the scripts in `ablations/`, their stdout and every engine transcript
they produced in `stdout/`, the round-four sources and products they were built from in `baseline/`, and
the private trees each ran in under `ablation_tree/`. The control run that shows a private tree reproduces
the repository transcript byte for byte is `stdout/control_S9.out` and `stdout/control_S10.out`.
