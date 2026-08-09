# S9 export — structural inputs, without in-run preservation claims

**User decisions, 2026-08-08.** Read `CLAUDE.md` and
`research/pde_ledger_v3/S9_REWRITE_PLAN.md` first. The rewrite plan remains closed. This directive
authorizes only the structural-input export repair and check removal described here.

Repository root: `/var/projects/toy_physics`. Ledger root: `research/pde_ledger_v3`.

## Export boundary and `identity3`

The declared-input export boundary admits both `KNOB` and `STRUCTURAL` declarations. The writer reads
declaration names, live values, and class tags from the existing annotation extractor and handles every
admitted declaration through the same generated-record and exact-`srepr` reconstruction path. Do not name
an individual declaration in the export path, inspect a payload for a symbol, or add an exclusion list.

`identity3` is a constructed object, not a declared parameter. Construct it by reading the component count
from `spatial_coordinates`, rather than typing that count again, and classify it as `DERIVED`. The class
boundary then keeps this fixed-count matrix out of a ledger whose structural inputs can be consumed at a
different component count. This is a classification correction, not a name-based export exception.

## Engine machinery is `CONTROL`

Classify the emitted-tag uniqueness state, standard-name substitution table, standard-name collection,
posited output-class map, production-dimension metadata key, action package table, and the emitted-name
string built while traversing packages as `CONTROL`. They are Python engine machinery, not structural
physics inputs or computed physics objects. Do not change their values or behavior.

`CONTROL` now has two meanings in this engine: an ablation control coefficient and engine machinery. A
reader of the knob inventory must not interpret every `CONTROL` declaration as a tunable physics
coefficient.

## Remove the in-run check battery

Delete the previous round's existing-record preservation emissions, classified-versus-selected input
emissions, consumer-symbol binding emissions, all residuals supporting those emissions, and their asserts.
Delete supporting bookkeeping and refactoring that no longer serves the pre-existing checks. Do not repair
or replace any part of that battery, and do not add a check, guard, helper, registry, runner, or name-based
rule.

The reason is structural: the export writer's operands are its own selected inputs. An in-run comparison of
those operands cannot establish that selection or preservation is right. Selection and preservation are
claims made by an external diff against the committed generated export.

Keep the pre-existing full-ledger round-trip, count, computed/absent-dimension count, class-tally, and
production-dimension partition emissions. Keep the pre-existing guards on the round-trip and count and
class-tally residuals. The exact-`srepr` round trip continues to traverse every generated record.

## Preservation and verification

Run S9 and capture literal stdout. Diff the regenerated `S9_exports.py` against the committed file. Every
committed record must remain byte-identical; the only export additions are the declarations admitted by the
structural class boundary, and `identity3` must not be among them. Also compare the new export with the
previous-round generated file to show that the `identity3` record is the only previous-round record removed
or changed.

Report the complete `.out` diff against the committed baseline. Report added and removed line totals for
this round.

Outside the repository, weaken the structural export gate back to `KNOB` only and run that copy. Report its
exit status and diff its generated export against the correct generated export. This ablation demonstrates
the boundary explicitly: the run checks exact reconstruction, population consistency, and the class tally
for the records it actually writes; only the external generated-file diff catches omitted structural
records and establishes preservation against the committed baseline.

## Files and invariants

Change only:

- `scripts/S9_light_requires_shear_sympy_audit.py`;
- `scripts/S9_exports.py`, only by running the engine;
- `scripts/out/S9_light_requires_shear_sympy_audit.out`, only as literal stdout from that run;
- this directive.

Do not change the derivation, action, ansatz, assumptions, premises, computed values, production-dimension
metadata, key formatter, or Wolfram engine. Do not state in this directive or a source comment what an
emitted count or tally evaluates to.
