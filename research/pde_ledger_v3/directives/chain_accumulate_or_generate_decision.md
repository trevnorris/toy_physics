# Decision — generate the whole-model ledger from per-step exports

**Status: AUTHOR DECISION, 2026-08-11.** This decides one question only: whether the whole-model list is
stored repeatedly or generated from step-owned exports. It does not decide names, row identity, or a build.

## Decision

**Each step exports only the rows it defines or re-derives; the one authoritative whole-model list and the
true knob frontier are GENERATED on demand from those per-step exports, not accumulated into every write.**

The per-step files are the authorities for their own rows. At one repository revision, their generated fold
is the sole authoritative whole-model view: every distinct object defined anywhere is present, and a later
re-derivation supplies its current classification and comparison evidence without republishing untouched
upstream rows. The frontier is exactly D7's set difference: the union of step-local `KNOB` rows minus every
object a later step marks `DERIVED`. Thus this choice delivers both required outputs from the same sources.

## Measurement

The literal measurement is `/tmp/measure_export_copying.stdout`; its script is
`/tmp/measure_export_copying.py`.

| committed export | total rows | exact upstream copies | changed shared rows | new rows |
|---|---:|---:|---:|---:|
| `S9_exports.py` | 165 | 0 | — | 165 |
| `S10_exports.py` | 617 | 162 (26.26%) | 3 | 452 |

S10 contains every S9 key: no upstream row is missing, all 162 untouched rows are exact copies, and the
three non-copies are deliberate S10 re-derivations. The accumulated artifact is correct today. A local-only
S10 export would contain its 455 S10-set rows; generation over S9 and S10 would still yield all 617 rows.

The handoff and step plan leave 18 named builds: S11, S11b-A/B/C, S12, S13, S14a, S14–S20, S20a, and
S21–S23. Exact future row counts do not exist. A mechanical projection that repeats only the observed S10
transition (452 new, 3 re-derived per step) puts 614 copied rows in the next 1,069-row artifact and 8,298
copied rows in the final 8,753-row artifact: 80,208 copied-row instances across the 18 writes. Even the
no-growth/no-overwrite scenario repeats today's 617 rows 18 times, or 11,106 copies. These are storage
projections, not expected physics results.

## The value/name distinction survives measurement

S10 directly binds 11 S9 rows as computation inputs and independently re-derives 3 more. Of its 162 exact
copies, 151 are not directly bound by S10 at all. A name or reference edit confined to any of those rows
cannot move S10's derivation, yet the accumulated S10 file changes because it republishes the row; its
whole-file S9 digest changes too. A changed value that S10 actually consumes is different: S10 computed
with that object, so its physics must be recomputed.

Generation removes the recursive copy coupling: a producer's spelling/reference edit changes that producer
and the generated view, not every later step-owned row set. It does not erase real direct dependencies.
S10's 11 live lookups, source spellings, and whole-file input digest remain coupled today; changing that
provenance contract is outside this decision. Therefore the distinction is real and measured, but this
storage choice does not pretend to solve every source-level invalidation.

## Reliability and losses

Accumulation gives four real benefits: one already-materialized import gives a consumer the complete state;
each committed file is a self-contained snapshot of what that step saw; a diff between adjacent exports
shows an overwrite; and `KNOB → DERIVED` plus the moving `step` field marks retirement at the write.
Those are useful audit affordances, and generation gives them up: a local export is not the whole model,
adjacent local-file diffs are not state diffs, and retirement becomes visible only in the generated fold.

The cost is an extra correctness premise at every later step. The final accumulated list is correct only if
every intermediate writer preserved every untouched row and reference. The current S10 transition satisfies
that premise, but the obligation grows once per row per downstream artifact. The generated view has no
intermediate copies to become stale or be re-pointed: given the same per-step rows and fold semantics, it is
fresh whenever generated. Both choices still fail on a wrong or omitted producer row, a wrong class, or a
wrong decision that two rows denote the same object. Generation is better for the requirement because it
removes a failure mode while retaining the producer evidence from which both the full list and frontier are
derived.

## Consumers that exist today

The literal inventory is `/tmp/export_consumer_inventory.stdout`.

| consumer | measured state | consequence of this decision |
|---|---|---|
| `scripts/S10_brane_mode_spectrum_sympy_audit.py` | the only executable importer; reads S9, binds 11 rows, compares 3, then copies the entire ledger | physics imports/comparisons remain; only its accumulated writer/output shape changes |
| `scripts/S11_stray_longitudinal_sympy_audit.py` and `_asbuilt/…` | zero S10-ledger references | nothing breaks today |
| `directives/S11_SHARED_PHYSICS.md#Q6r` | specifies a future import of the complete `S10_exports.LEDGER`; no such rebuilt engine exists | that future consumer would read the generated whole-model view, not S10's local slice |
| `scripts/extract_knob_inventory.py` | scans engine annotations; it does not read either ledger or compute the global frontier | its present per-file report is unchanged and is not the required frontier |
| S10 comparator, Wolfram engine, record, and card | the executable paths read tag streams, not a ledger | no computation changes; the S10 record's export line citations would move when S10 shrinks |
| `steps/S9_light_requires_shear.md` | cites only S9's own first-step export | no architecture change |

Keeping accumulation would break none of these today, but it would preserve the measured false coupling and
make every future consumer inherit the growing copy obligation.

## Explicit supersession

- **D5 is superseded** only in its physical accumulation rule, its committed whole-state snapshots, and its
  “untouched rows identical in the next file” consequence. Flatness, current-value semantics, and comparison
  at a genuine re-derivation are not decided here.
- **D7's frontier equation remains; its physical overwrite marker is superseded.** A later step's own
  `DERIVED` row marks retirement, and the generated fold performs the subtraction.
- In `S11_export_chain_decisions_v2.md`, **F4 is superseded insofar as it requires regeneration of an
  accumulated S10 export**; regeneration/evidence for S10-owned re-derivations remains. **F2, F3, F5, F6,
  and F7 remain. F1 was already superseded and is not revived.** In F3, “merged export” now means the
  generated whole-model view.
- For avoidance of a contradictory record, the accumulation clauses D6, D9, and D10 in the already-blocked
  `S11_sympy_export_decisions_codex.md` are also superseded. No blocked naming/identity design is revived.

The recorded cost of supersession is the loss of the ready one-file import, committed adjacent-state diff,
and write-time `KNOB → DERIVED` mutation; consumers needing the whole model require the generated view.

## Not settled here

The exact forward copy totals are unknowable until the remaining steps expose their own-row and
re-derivation counts; those exports would settle the projection. `DEFECT_REGISTER.md#c20`'s object-identity
problem remains common to both choices: neither a stored overwrite nor a generated fold can decide a true
meeting from the files that exist. It requires its own physics-bearing decision and two-leg review once a
live downstream re-derivation exists. This decision adds no fifth identity mechanism and solves no name.
