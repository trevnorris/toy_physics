# Independent physics review — the S11 PY engine decision list

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_py_decisions_v2.md`

It is **orchestrator-written**. It is the list a builder will be handed to rewrite
`research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`, which on running writes
`research/pde_ledger_v3/scripts/S11_exports.py`.

⚠⚠ **It is the one artifact in this pipeline that gets checked once and then trusted.** Everything
downstream of it is reviewed twice; it is reviewed here or nowhere. A defect in it becomes a build round
plus two more legs, and the last time this gate was skipped it cost six spec defects — one of them an
acceptance test that would have passed with the defect still in place.

## Read in this order — ⛔ the order is load-bearing

1. `/var/projects/toy_physics/CLAUDE.md` — the sixteen rules the list must satisfy.
2. `git show cf4a21a4:research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — the shared physics spec,
   1149 lines. This is the source of truth for **what** gets computed.
3. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_exports.py` — the committed artifact the
   list's measurements are about. **Open it and measure it yourself.**
4. `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — the
   engine that wrote it, i.e. the existing precedent for every mechanism the list decides.
5. `git log -1 --format=%B 11bf8e05` — the five findings that blocked v1 of this list, and
   `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_py_rebuild_decisions.md`, the blocked
   v1 itself. ⛔ v1 is a record, ⛔ not a directive.
6. `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` — `F1`–`F7`,
   already settled and two-legged. The list under review must be consistent with these or must say
   explicitly where it is not and why.
7. `/var/projects/toy_physics/research/pde_ledger_v3/S9_REWRITE_PLAN.md` — `D1`–`D12`.
8. ⭐ **Only now** open the artifact.

⚠ Do not read the artifact first. Form your own view of what the chain requires, then read what the list
decided. Reading it first anchors you to its framing, and its framing is the thing under test.

## Do not read

- The other review leg's output. Two legs are running on this identical prompt; you get no hint of the
  other's findings.

## What to check

### ⭐⭐ 1 · RE-RUN THE MEASUREMENTS. ⛔ Do not read them off the list.

The list stands or falls on a table of six measurements, `M1`–`M6`, each claimed to have been taken from a
committed artifact. **Take each one yourself and report your own numbers**, whether or not they match.

⛔⛔ **A prose re-derivation is worth nothing here.** Write a script, run it, and save **both the script and
its literal stdout** to named absolute paths under `/tmp`. Report those paths. Without them your
measurement claims will be discarded — a claim with no computation behind it is the exact defect class this
whole rebuild exists to remove, relocated into the review.

⭐ `M4` is a **whole-namespace census**: every `<QUANTITY>` name the spec defines, against every key in the
committed ledger. Re-run it your own way. ⚠ If the list's census missed a name, or counted one that is not
a quantity, say which and what it collides with. A single missed collision changes `P2`.

### ⭐⭐ 2 · `P2` IS THE LOAD-BEARING DECISION. ATTACK IT.

`P2` puts the run scope into the storage key and so departs from `F1`, which says keys stay flat because
**the collision is the measurement**. The list's defence is `M4`. Test it:

- **Find a counterexample.** Is there any object S11 derives that is genuinely **the same mathematical
  object** as one an earlier step already exported — so that a flat key would make them meet and compare,
  and `P2` silently prevents that meeting? If one exists, `P2` deletes a measurement, which is rule 6.
- Does `P2` break any live read of the imported ledger — `Q6r`'s `LEDGER[name]['dimension_key']` lookup, or
  any other consumer? Show the lookup and what it resolves to under `P2`.
- Is the "two disjoint key populations" claim true of the committed artifact, or is there a key that is
  both? A key that is in both populations under `P2` is a silent overwrite.
- ⭐ Is the list's distinction between "the engine that computed it" (`D11` forbids it in the name) and "the
  run scope the spec's tag grammar assigns" real, or is it a rationalisation for a producer prefix that
  `F1` blocked and both legs rejected once already?
- ⛔ Is there a **third option** neither `F1` nor `P2` considered?

### ⭐⭐ 3 · CAN EACH GUARD ACTUALLY FAIL? ⛔ Construct the failure; do not reason about it.

For `P4` and `P5`, build the bad case and show the literal output. A guard that cannot fail is worse than
no guard, and this list already carries one demonstrated case of that class.

- **`P5`, chain integrity.** Write the accidental write it is meant to catch — mutate an imported row the
  step never declared it would touch — and show that `P5`'s three operands make it visible. Then try to
  defeat it: is there a mutation that changes an imported value and still leaves the symmetric difference
  empty? What about a value that is `==`-equal but not the same object, or vice versa?
- **`P4`, same-object.** `P4` says compare **as objects**: `Symbol` identity including assumptions,
  otherwise a residual reducing to zero. Construct two objects that this predicate calls the same and that
  a physicist would not, and two it calls different that are the same. ⚠ `Symbol("D")` versus
  `Symbol("D", integer=True, positive=True)` is a known trap in this repo — does `P4` catch it, and does
  it survive the reconstruction round trip of `P3`?
- **`P4`'s vacuity operand.** The list argues that `P2` makes a tag-derived collision impossible, so the
  guard must publish its own vacuity. Is that argument right? If a tag-derived collision **is** still
  possible under `P2`, the guard is live and the vacuity claim is wrong — which is the more interesting
  finding.
- **`P1`.** Is the export set decidable from the tag alone, with no judgement left to the builder? Walk the
  spec's tag inventory and find a tag the rule cannot classify. ⛔ The list claims excluding every `_LOCAL_`
  tag loses nothing — check that claim against the spec, naming any object that would be lost.

### ⭐ 4 · THE RULES THE LIST MUST NOT BREAK

- **Rule 5 — the list says what to compute, ⛔ never what anything equals, is expected, or was measured
  about the physics.** Does any line hand a builder a value, a count, a membership, or a sign it is
  supposed to compute? ⚠ v1 was caught stating which coefficients resolve. Check `P7` and `F7` especially.
  ⚠ A **prohibition** leaks as surely as an assertion.
- **Rule 3 — name the object, ⛔ do not specify the recipe.** Which items specify *how* rather than *what
  must become true*? Nine rounds on the spec were spent removing exactly that, and thirteen of sixteen bred
  defects lived in machinery invented to stop two engines describing something differently.
- **Rule 2 — a script prints computed objects and never states conclusions.** Does every guard emit both
  operands and the residual?
- **Physics filter (standing user instruction).** Report an item only if it catches a way the **physics**
  could be wrong. ⛔ Do not report "the script would be wrong on a different input", and ⛔ do not report
  process or tooling preferences. If an item in the list is **process rather than physics**, say so — that
  is a finding, because it is cost with no correctness behind it.
- **Rule 11 — cost is never a reason** to drop a control, narrow a check, or skip a leg. If the list drops
  something to save work, that is a finding.

### ⭐ 5 · WHAT IS MISSING

⛔ The most expensive defects in this project have been **absent computations**, not wrong ones. Ask of the
build this list describes: what could go wrong in the export chain that **no** item here would detect?
⚠ Name it concretely, with the mechanism.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way a wrong result could pass
undetected. ⛔ Do not report style, naming taste, or "a builder might misread this" without showing the
second reading and what it would produce.

## ⛔ Constraints on how you run

- ⛔ **Read-only on the working tree.** Copy anything you want to modify to `/tmp` and modify the copy.
- ⛔ Wrap every CAS run in `timeout 600`. A 600s hit is a failed check — report it and move on. ⛔ Never
  raise the timeout.
- ⛔ Do not run Mathematica for this review; the artifact is a document about a SymPy engine and the licence
  has two seats.
- ⭐ Save every script and its literal stdout to named absolute paths under `/tmp`, and report those paths.

## Report

For each finding: **what is wrong**, **the mechanism by which a wrong result would survive**, **the literal
output that shows it**, and **what must become true instead** — ⛔ not a rewrite of the list.

⭐ Separate findings that block the build from findings that do not. ⭐ State plainly which of `M1`–`M6` you
reproduced and which you could not. ⛔ If you find nothing on an item, say which computation you ran that
could have found something and did not — "nothing found" with no ablation behind it is not a result.
