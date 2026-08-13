# S11 §Q8b — fix round 1. ⭐ One item: an obtained count is erased by an unrelated status.

## Authority and boundary

Edit `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` **only**. ⛔ Change no other
file. `CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` is the sole physics
authority. `research/pde_ledger_v3/directives/S11_q8b_build_directive.md` is **closed**.

⭐ **Baseline is `94e14e42`** (the reviewed §Q8b build). Verified there by two independent legs and not
to be revisited: the three-family pool and its 1–1 accounting; the promotion path; point-invariance of
every component-wide payload; equation-driven payloads; the count-record shape and companions; the
publish with the stratum rows aboard. Measurements:
`directives/_measurements/S11_q8b_fix_round1_brief.md`.

## ⛔⛔ 1 · BLOCKING — THE `VALUE` FIELD ANSWERS THE WRONG QUESTION

`emit_component_count` (`:1747-1758`) sets, at `:1754`:
`payload_value = generic_value if status in ("CONSTANT", "VARIES") else NOT_DEFINED_ON_COMPONENT` —
⛔ the `VALUE` field is conditioned on the **constancy status**, not on whether a component-wide count
was **obtained**. Measured independently by both review legs through the real promotion path: count
records whose computation bound an integer emit `VALUE: NOT_DEFINED_ON_COMPONENT` beside
`STATUS_TOKEN: UNDECIDED` — the engine's own obtained count, erased.
⛔ Against the physics authority: `S11_SHARED_PHYSICS.md:737-739` — `VALUE` is *the count on the
component where it is defined there, **whatever the status**; the single token
`NOT_DEFINED_ON_COMPONENT` where the engine **did not obtain** one*. `:754-756`'s
`UNDECIDED`+`NOT_DEFINED` case has *"cannot build a component-wide Q3/Q4 object"* as its antecedent —
built-but-constancy-open is not that case, and `VARIES` (`:742-745`) keeps the obtained `VALUE`.
`:710-711` makes the counts and statuses cross-engine comparison rows — an erased operand makes a
measured count uncomparable. ⚠ This is also the build's one value-conditional emission.
⇒ ⭐ **`VALUE` reports what was obtained: the obtained component count, whatever the status token;
`NOT_DEFINED_ON_COMPONENT` exactly where no component-wide count was obtained.** The status token is
where open constancy lives. ⛔ Do not specify the routine (rule 3). ⚠ If you read `:754-756` as
contradicting this, that is a §10 report — ⛔ not a licence to keep the erasure.

## Boundaries

- ⛔ No memory cap, no timeout, no handler that swallows a failure to make a run finish.
- ⛔ Do not change `PACKAGE_DIMS`, the `D` range, or any package.
- ⭐ Everything the two build legs verified (list above) must not change; ⛔ nothing outside the item's
  property may change what is emitted.
- ⛔ **No expected value and no acceptance criterion referencing one** (rule 5). A changed value is a §10
  finding, never something to tune away.
- ⭐⭐ Runtime (user, 2026-08-13): a cell may exceed 600 s provided its output is streaming; a long
  silent stretch is the failure.

## Acceptance

⚠ Checked against the item: at baseline, the demonstration below fails by both legs' measurement.
1. `/home/trevnorris/.s11_build/repro_d5.py` runs `run_cell('MAIN', 5)` to completion; literal stdout
   and timing.
2. **Through the real promotion path on a `/tmp` drive copy** (⛔ never in-tree; ⛔ never by hand-writing
   output tags): the demonstration is **invalid unless it produces at least one component-count record
   whose computation obtained a count and whose `_STATUS` is `UNDECIDED` at baseline** — an empty
   ordering, an unpromoted drive, or an all-`CONSTANT` promotion demonstrates nothing (all three were
   measured green-while-unfixed by the brief legs). ⭐ The **same driven component** is run at baseline
   and after the fix. Report, with literal stdout, every component-count record beside what that
   record's own computation bound. ⛔ **The criterion is equality, in both directions**: `VALUE` must be
   exactly the count the computation bound where one was obtained, and exactly
   `NOT_DEFINED_ON_COMPONENT` where none was — any inequality either way fails, a wrong non-sentinel
   value included. (Two run-produced operands; no expected value is introduced.) ⛔ Do not manufacture a
   value where none was obtained (spec `:744-745`).
3. Report every tag whose payload moved between baseline and the fix, and where each is reported. ⛔ A
   payload outside the item's property that moved is a finding to report under §10, not to keep silently.
4. ⛔ Do not run the full package loop on the in-tree engine. ⭐ Any engine-cell run happens with
   `/tmp/s11_watchdog.sh` armed by the builder (capture the pid; kill that pid when the run ends — ⛔
   never `pkill` by name).

## Deliverable

The edited script, and a note saying per item what changed, plus every value that moved anywhere and
where you reported it.
