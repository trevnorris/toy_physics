# S9 export re-key — fix round 1, decision list

**You author the directive prose and apply the change.** This is a decision list, not a directive.
Two independent review legs reported on your build. Four defects survived verification by the orchestrator.
The scope of this round is those four and nothing else.

Files in scope are unchanged from the previous round: `scripts/S9_light_requires_shear_sympy_audit.py` and
the generated `scripts/S9_exports.py`. ⛔ No other file. ⛔ No new file. ⛔ Not the `.wl`.

---

## F1 · The key's `D` is a typed literal. The key must be read out of the computation.

`sp.Integer(3)` is typed at six production sites, and the key formatter concatenates the literal string
`"_d3"` and raises on any other value. Measured by ablation: the spatial algebra was rebuilt at four
components with the production tags left exactly as you wrote them — **the key set came back
byte-identical** while a majority of the exported values moved, so a key positively misstated the `D` its
object was computed at.

**Decision: the component count in a key is read out of the same structure the computation uses**, at the
point of production, and the formatter accepts whatever count it is given. After this change there is no
literal component count anywhere on the path from the computation to the key.

⚠ The engine runs at one component count today. That is not a reason to type it — it is the reason the
defect is invisible without an ablation.

## F2 · The general/fixed cut is in the wrong place inside the dimension block.

The dimension block is tagged general-`D` as a unit. It is not a unit: **the solve's inputs are built from
the fixed-component action, and only the solve's outputs are symbolic in `D`.** `D` enters that solve
through the energy-density target, not through the term count — which is why the outputs can be general
while their inputs are not.

**Decision: put the cut between the solve's outputs and the solve's inputs.** Classify by which side of
that boundary an object falls on. ⛔ Do not hardcode a list of names, and ⛔ do not classify by inspecting
free symbols — that test was measured wrong in both directions before this round began.

⚠ Your directive already stated the correct principle for `lagrangian` and `equation_of_motion` — that a
`t, x, y, z` construction fixes an object to three spatial components — and then applied the opposite
classification to objects built the same way. Make the directive and the code agree.

## F3 · For four records the recorded classification is unmeasured.

The four entries carrying `_d3` before this round take that suffix from their hand-typed Python name, not
from their production tag. Measured by ablation: retagging that whole block general-`D` leaves the
generated export **byte-identical**. So for exactly the four entries your directive presents as the
precedent the convention adopts, the tag is connected to nothing.

**Decision: a key's suffix comes from the production tag alone.** An internal name must not be able to
contribute a component count to a key.

## F4 · The partition is computed and never printed.

The regenerated `.out` is byte-identical to the committed baseline, so **not one emitted line carries an
export key or the partition.** F1, F2 and F3 were therefore all invisible in the record, and no future run
could catch a recurrence.

**Decision: emit the D-partition as a computed object**, alongside the existing export tallies, following
the same print-then-guard shape those already use. ⛔ Not a boolean. ⛔ Not a conclusion. The operands and
the residual, then the guard.

---

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** Measured three times in this project: a half-fix
passed the new test, the whole suite and the full battery, and produced byte-identical output.

**Build the weaker implementations and show the test fails on each.** At minimum: the formatter fixed but
the production sites still typing a literal; the production sites reading the structure but the formatter
still hardcoding one suffix; F2 applied without F1; F3 applied without F1.

**Demonstrate F1 by ablation, not by argument.** Rebuild the computation at a different component count
with nothing else touched, and show what the keys do. Run a control at the unchanged component count first
and show it reproduces the current export, so the ablation harness itself contributes nothing.

⛔ Do not state, in the directive or in a comment, what any count, tally, term number or partition size
comes out as. Emit it and let it be read.

## Report, do not fix

If an object does not fall cleanly on one side of the F2 boundary, **name it and stop.** An ambiguous
object is a finding about the export boundary, not a puzzle to solve inside this round.

## Still out of scope

⛔ The three open S9 items (`wavevector_norm_dimension`'s name, the placeholder-naming class,
`q_dimension`). ⛔ Any change to the derivation, action, ansatz, assumptions or computed values. ⛔ The
cross-engine comparator. ⛔ Anything in S10.

## Deliverables

- The two changed files, and your updated directive at `directives/S9_export_dkey_directive.md`.
- The literal stdout of the re-run and the `.out` diff against the committed baseline.
- Every ablation script and its literal stdout, at named absolute paths, including the weaker-fix runs and
  the F1 control.
