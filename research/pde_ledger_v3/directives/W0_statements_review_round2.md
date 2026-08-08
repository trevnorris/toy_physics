# Independent review — W0 builder and integrator statements, ROUND 2

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Artifacts — Codex-written, second round

- `/var/projects/toy_physics/research/pde_ledger_v3/directives/W0_builder_statement.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/W0_integrator_statement.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/directives/W0_restatement_verification_report.md`

⭐ The **builder statement** will be read by **four independent engine builds** (S9-WL, S9-Python, S10-WL,
S10-Python). ⛔⛔ **An error in it makes all four engines wrong the same way**, defeating dual-engine
verification by construction.

## ⛔⛔ ATTACK THE REPAIRS

⚠ **In this repository, repairs have repeatedly introduced NEW defects into the material just changed.**
Round 1 of these statements was reviewed by two legs and **nine** defects were fixed. ⭐ Your first job is
to check whether each repair actually closes what it claims **and** whether it broke something round 1 had
right.

⭐ The round-1 defects, as fixed. ⛔ Verify each independently; ⛔ do not accept the artifact's account:

1. **The defining residual did not pin the quotient.** ⭐ A residual now ties `v_T²` to the selected root
   and `K`. ⚠ **Note carefully:** in every in-scope cell that quotient is algebraically identical to a
   recombination from material coefficients, ⇒ ⛔ **the residual alone is zero for a recombined build too.**
   ⭐ So: does the surrounding machinery actually make recombination detectable, or has a zero-by-
   construction object simply been added next to another one? ⭐ **Exhibit a recombined build and follow it
   through every required artifact.**
2. **Selection was not bound to the predicate.** ⭐ It now claims membership is the truth set of the emitted
   predicate. ⭐ Can a build satisfy the words while selecting by something else?
3. **A required cardinality leaked an expected value.** ⭐ It is now removed. ⭐ Check nothing else in either
   document supplies a count, a target, or an expected verdict a builder could iterate toward.
4. **The sentinel mutation distinguished nothing.** ⭐ It now enters upstream of the quotient. ⭐ Work out
   whether an honest build and a recombined build genuinely diverge under it — ⛔ and whether a builder
   could satisfy it while still recombining.
5. **The control leaked the target.** ⚠ A non-negative scalar of dimension `[1,-1,0]` built from the two
   material coefficients is **uniquely determined**. ⭐ Check no wording reintroduces an engine-side
   construction of it, ⛔ and that no coefficient is named in the builder statement.
6. **Boundaries contradicted the regression bar.** ⭐ Check the exception is now exact — ⛔ neither too
   narrow (impossible build) nor too broad (a real change hides in it).
7. **The manifest was a self-consistency check.** ⭐ It now requires name **and** value and is declared a
   review input. ⭐ Does anything still let a build *pass* by agreeing with itself?
8. **Placement was unpinned across four builders.** ⭐ Check the grain is now unambiguous, including S9's
   direction specialisations, and that a later cross-engine row could bind the pairs.
9. **`L_T = 0` was offered to a reviewer as evidence.** ⭐ Check the integrator now says plainly what
   carries information and what does not.

## ⭐ Then attack what round 2 introduced

- ⭐ **The launch hold.** The builder statement is held pending a separate repair, because one of the four
  engines does not currently run. ⭐ Verify that engine's state yourself. ⭐ Is the hold's release condition
  precise enough to be checkable, or is it a sentence someone will wave through?
- ⭐ **The mutation-only run.** Two runs now exist — production and mutation-only — with different rules
  about what may change. ⭐ Can a builder blur them, or report a mutation artifact that never ran?
- ⭐ **Per-member speed bundles.** A bundle is now required for **every** predicate-true root. ⭐ Check
  against the committed outputs that this is well-defined in every in-scope cell.
- ⭐ **Can a compliant build still be worthless?** ⭐ Construct one, or show each known mode is excluded by
  the words as written.
- ⭐ **Four-builder divergence.** ⭐ Is anything still ambiguous enough that two builders reasonably diverge,
  manufacturing a specification artifact rather than a physics disagreement?

## Method

⭐ Read the sources of truth **first** — `reduction/relations.yaml`, `reduction/quantities.yaml`, both
`checks_S*.yaml`, the four engines, the committed outputs — form your own view, ⭐ **then** read the
artifacts. ⛔ Not in the other order.
⚠ `W0_decision_list_round2.md` is the orchestrator's input to this build. ⭐ Read it **after** forming your
own view; ⛔ it is under test, not a source of truth.

- ⛔ Do not modify the repository. ⛔ Do not write or run an engine change.
- ⛔ **Use absolute paths for every file you touch outside the repository.** ⚠ A `cd` into a temp directory
  has already failed silently in this session and edited live files.
- ⭐ Where you claim a build could comply and still be wrong, **exhibit the complying implementation** and
  follow it through every artifact the statement requires.
- ⭐ Give the literal command and output for every per-cell claim.
- ⛔ Wrap any long run in `timeout 900`; ⛔ never raise it; ⛔ never more than one Mathematica kernel at a
  time — the licence has **two** seats.

## Physics filter

> Report a finding only if it changes what the project computes, what it may claim, or what method it
> adopts.

⭐ Worth most: *"a recombined build still survives every required artifact"* · *"this repair broke
something round 1 had right"* · *"this still leaks"* · *"two builders would reasonably diverge here."*

⚠ A leg returning *"nothing survives the filter"* is weak evidence on its own. State what you checked, what
you could not, and what would have had to be true for you to find something. ⛔ Do not manufacture findings.

Rank most-severe first: **claim · evidence (quotation with file:line, or a literal command and its output)
· what must change.**
