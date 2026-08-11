# Independent review — the repaired `S11_SHARED_PHYSICS.md`

## Artifact under review

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`

⭐ It was just repaired to close `C17` and `C18`. The **baseline it was repaired from** is
`/tmp/claude-1000/-var-projects-toy-physics/36f37d88-e717-46ce-a790-6f9d1ef3d7bc/scratchpad/S11_SHARED_PHYSICS.baseline.md`
— ⭐ diff against it to see exactly what moved.

⚠ This file is the **shared spec both engines read**. An error in it makes **both** engines agree on the
same wrong thing, or **disagree for a reason that is not physics**. That is why it is gated twice.

You are one of two independent legs; the other is not visible to you.

## ⛔⛔ DO NOT READ

- `directives/S11_C17_C18_spec_repair_build_directive.md` — ⚠ **the directive the repair was built from.**
  ⭐ An artifact can satisfy its directive completely and still fail its source. ⛔ That case is exactly
  what this leg exists to catch; reading the directive blinds you to it.
- Any file under `_scratch/`.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not open the repaired file first

Reading the artifact first anchors you to its framing, which is the thing under test.

1. `research/pde_ledger_v3/DEFECT_REGISTER.md`, entries `C17` and `C18` **in full** — ⭐ including the
   measured divergences they record.
2. The committed outputs those entries cite:
   `mathematica/out/S11_stray_longitudinal_mathematica_audit.out` and
   `scripts/out/S11_stray_longitudinal_sympy_audit.out`. ⭐ **Open them.** ⛔ Do not take the register's
   quotation on trust — it has carried stale line numbers before.
3. `directives/S11_C17_C18_spec_repair_decisions_v2.md` — ⭐ the seven decisions `T1`–`T7` the repair owes.
4. **Write down, before step 5:** what a spec must say so that two independently built engines produce
   **comparable** stratum information, given one CAS can decide real-domain questions the other cannot.
5. **Only now** the repaired file, and the diff against the baseline.

## What to check

### 1. ⛔⛔ DOES THE REPAIRED FILE STATE, OR IMPLY, AN ANSWER?

⛔⛔ **This is the highest-severity check and it fails the whole repair.** The register entries contain
**measured results of this step's own system** — root counts, a named coincidence locus, which engine
emitted which stratum.
⭐ **The test is DERIVABILITY, ⛔ not literal presence.** A justification, an example, or a
"for instance" from which a builder could **derive** a count, a value, a sign or a named locus **is a
leak**, even though no number appears.
⇒ ⭐ Go clause by clause through every added sentence and ask: **could an engine builder reading only this
file work out what the answer is?**

### 2. ⭐⭐ IS `C17` CLOSED, OR RELOCATED AGAIN?

⚠ The previous attempt at this repair **relocated** the defect: it scoped the rerun to the component, and
the unsupported count simply moved from an arbitrary point to a generic free-parameter reading, with the
sub-locus where the count changes getting **no tag at all**.
⭐ **Test the repaired text against that.** Take a component from the committed outputs, follow the
repaired instructions literally, and report what an engine would emit. ⛔ Is there still a place where the
mode count moves and **nothing says so**?

### 3. ⭐⭐ IS `C18` CLOSED FOR THE MEASURED CASE?

⭐ Apply the repaired text to the `XFORM_EXTRA, D = 2` divergence in the committed outputs. ⛔ Under the new
wording, what does **each** engine emit, and what could a comparator then conclude?
⚠ *"Not compared"* is **not** a pass for a region where the mode count may move — ⭐ say whether the repair
produces a **computation**, a **finding**, or **silence**.

### 4. ⭐⭐ IS THE RESULT BUILDABLE BY TWO ENGINES THAT NEVER SPEAK?

⭐ For every new obligation: can **each** CAS discharge it independently? ⛔ Does anything require a
capability only one of them has? ⚠ `§5` at the top of the locus protocol already records a **measured**
asymmetry between the two CASes — ⭐ check the repair against it.
⛔ Does any new instruction ask an engine to **import**, **agree with**, or **adjust toward** the other?
⚠ That would destroy the only blind build in the project.

### 5. ⭐ DID THE REPAIR DO MORE THAN IT SHOULD?

⭐ Diff against the baseline. ⛔ Is anything changed that `T1`–`T6` do not require? ⭐ Name what was
damaged, and check in particular that **none of `§5`'s five original locus objects was deleted or
weakened** — they were written against a measured CAS limitation.

### 6. ⚠ STALE PROSE ABOVE THE FIX

⭐ For every edited section, read the section's **introduction** and any summary above it. ⛔ Is there a
sentence the edit made **false**? ⚠ This has happened in this file before and was caught a gate later.

### 7. ⭐ WHAT IS MISSING?

⭐ Compare your step-4 list. ⛔ Which of `T1`–`T7` is **not** actually discharged by the text? ⭐ Quote the
decision and the sentence that was supposed to satisfy it.

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, a way a **claim would be
unsupported**, or a way two engines built from this file would produce **wrong or incomparable** results.
⛔ Not style, ⛔ not wording taste, ⛔ not marker conventions.

## Method

- ⭐⭐ **Quote both sides of every finding**: the repaired file's text, and the source it fails against.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
  ⛔ A path that resolves is not a source. ⛔ A line number in the register is not a source — open it.
- ⭐⭐ Checks 2 and 3 are settled by **following the repaired instructions on a real case from the committed
  outputs** and reporting what an engine would emit. ⛔ A prose argument that it "should work" is discarded.
- ⭐ If you run a CAS: wrap **every** kernel run in `timeout 600`; a 600 s hit is a **failed** check —
  report it and move on. ⛔ Never raise the timeout, ⛔ never more than one kernel at a time (two seats).
  ⛔ Copy anything you execute to `/tmp` and run the copy.
- ⛔ Do **not** edit anything in the working tree. Read-only.
- ⭐ End with: which of `T1`–`T7` you judge discharged and which not; whether the file leaks; and what you
  checked that could have failed and did not.
