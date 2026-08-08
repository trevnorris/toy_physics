# S9 export re-key — fix round 2, decision list

**You author the directive prose and apply the change.** This is a decision list, not a directive.

Two independent legs reviewed your previous round. Four items survived verification by the orchestrator,
and two of them are defects **this round introduced**. Scope is those four and nothing else.

⛔⛔ **This round REMOVES code in two places and adds no mechanism.** If you find yourself adding a check,
a guard, a helper or a collection, stop and say why — that is the signal this list is wrong, not that the
scope is bigger.

Files in scope: `scripts/S9_light_requires_shear_sympy_audit.py`, the generated `scripts/S9_exports.py`,
and — newly in scope, correcting an error in the previous list — `scripts/out/S9_light_requires_shear_sympy_audit.out`.

---

## F1 · The supplied dimension premises are classified three different ways. Classify them by one construction.

Four supplied reference dimensions are declared together as premises. They are **not** products of the
action: no spatial coordinate and no action term enters any of them. Today one is swept into the
solve-inputs group and carries a component count; one is held general by a hand-typed line naming it; two
reach the export through a separate path that bypasses classification entirely.

Measured consequence: one of them now exports under a key asserting a component count **while its value
still contains the symbol `D`** — and that value is byte-for-byte identical to another entry exported
under an unsuffixed key.

**Decision: the supplied references are one group, constructed as one group, and tagged as one group.**
Their classification follows from how they are built, like every other group.

## F2 · Delete the post-hoc classification override.

A single line assigns one output's classification by name, after its group was tagged. Your own directive
forbids exactly this: the boundary is expressed by constructing the collections separately, not by a
post-hoc list of names. Measured: delete that line and the run stays green while the object changes sides.

**Decision: the line goes.** F1 done properly makes it unnecessary.

## F3 · Delete the D-partition residual and its assert. Keep both operands emitted.

The "key-side" operand is parsed back out of the keys that the "production-side" operand just generated.
The two are related by an exact bijection, so the residual is structurally zero. Measured: a mis-tagged
object exits 0 with residual 0; **inverting the entire dimension-solve boundary exits 0 with residual 0**
and twelve keys swapped sides, leaving a symbolic value under a fixed key and a fixed value under a
general key. The residual polices the formatter, never the classification — and it is what let F1 ship.

⛔ **Do not repair it.** No check inside this engine can audit the classification, because the
classification is that check's own input. Auditing it requires running the computation at a second
component count, which is the ablation battery for the export writer that the plan forbids.

**Decision: emit both partition operands and stop there.** The partition becomes a printed computed object
a reader and a review leg can check. ⛔ No residual, ⛔ no assert, ⛔ no replacement guard.

## F4 · Regenerate the committed `.out`.

The repository's record of the S9 run predates the partition emissions and no longer reproduces from
source. The previous list's scope fence forbade touching it while requiring a change that alters it; that
was an error in the list, not in your build.

**Decision: re-run the engine and commit the regenerated stdout.** Report the complete diff against the
committed baseline.

---

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** Build the weaker implementations and show what
each does:

- F1 applied to some of the supplied references but not all of them;
- F2's line deleted without F1;
- the partition operands emitted but still derived from a single source.

**Demonstrate F1 by mutation.** Show, for each supplied reference, what its key does when the component
count the computation runs at is changed, and what it does when its group's tag is changed. Run a control
at the unchanged count first and show your harness reproduces the current export.

⛔ Do not state, in the directive or in a comment, what any count, tally, partition size or term number
comes out as. Emit it and let it be read.

## Report, do not fix

If an object does not fall cleanly on one side of the construction boundary, **name it and stop.**

## Out of scope

⛔ The three open S9 items (`wavevector_norm_dimension`'s name, the placeholder-naming class,
`q_dimension`). ⛔ The derivation, action, ansatz, assumptions and every computed value. ⛔ The `.wl`.
⛔ The cross-engine comparator. ⛔ Anything in S10. ⛔ Any change to what a key *means* — a key records the
construction that produced the object, not the object's sensitivity to `D`.

## Deliverables

- The three changed files and your updated directive at `directives/S9_export_dkey_directive.md`.
- The literal stdout of the re-run and its complete diff against the committed baseline.
- Every ablation script and its literal stdout at named absolute paths, including the weaker-fix runs and
  the control.
- A statement of how many lines this round removed and how many it added.
