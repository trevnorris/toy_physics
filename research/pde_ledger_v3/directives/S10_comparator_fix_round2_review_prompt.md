# Independent physics review — cross-engine comparator, fix round 2

## Artifact

Untracked, in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_cross_engine_comparator_directive.md
research/pde_ledger_v3/directives/S10_comparator_fix_round2_directive.md
```

**You are one of two independent legs with this identical prompt.**

## What this is

Two engines compute S10 independently — SymPy and Wolfram — and each emits objects under the standard name
of the object, so comparison is a **join on the name**. This script performs that join.

⚠ **This round repaired five things that are INVISIBLE on the current input.** They matter because the
instrument is meant to be reused on other steps, where the inputs differ. ⇒ ⭐⭐ **you cannot review this
round by running it on the committed transcripts alone — you must CONSTRUCT the inputs that exhibit each
defect.**

⚠ The two transcripts are committed inputs and must **not** be regenerated.

## ⛔⛔ The method this round demands

For each repaired property: ⭐ **build the input that exhibits the defect, run the repaired instrument on
it and show it is caught, then run the UNREPAIRED instrument on the same input and show it is not.**
⛔ Without both halves you have shown only that the current input is unaffected, which was already known.

⭐ Recover the unrepaired behaviour however you like — a `/tmp` copy with the change reverted is fine.

## What to check

1. **A derivative that carries a different dependence set.** Construct two payloads naming the same field
   but depending on different variables, and show whether they compare equal. ⚠ Then check the repair did
   not overreach: the two engines **order their function arguments differently**, and that must still not
   matter; and two genuinely different partial derivatives must still not become interchangeable.
2. **A payload symbol that collides with a CAS builtin.** Construct payloads using such names and show
   whether the symbol survives as itself. ⚠ Check the one genuine constant in the current input still
   denotes what both engines meant by it.
3. **An ordering normalisation was deleted.** Confirm nothing now sorts a sequence, and that no object
   that relied on that branch changed verdict.
4. **The accounting property.** Construct the input that would previously have broken it — a shared name
   emitted more than once — and show the categories still account for every shared name.
5. **A typed disagreement kind.** One failing row's kind was a literal. Confirm it is now read off the
   comparison, and that it reports correctly for a row you construct.
6. **⛔ Did this round break the previous one?** ⭐ **The nullspace residual must still pass BOTH halves** —
   fires on a wrong subspace, a partial rescale, a rank loss; silent under whole-vector rescale,
   permutation, and a general invertible change of basis. ⚠ A round that quietly breaks that is the
   specific failure this question exists to find.
7. **Does it report agreement it did not compute?** Account for every shared name; the categories must sum.
   Derive the join with your own parser and compare — ⚠ **a disagreement in either direction is a finding.**
8. **Account for every changed line** against the previous run.
9. **Conclusions.** Report any hardcoded expected value, count or tally in the script or the directive; any
   typed verdict standing in for a computation; any `assert` preceding the value it guards.

## Required method

⛔ **Write your own comparison script BEFORE opening the artifact**, run it on the two committed
transcripts, and save it and its literal stdout to named absolute paths. **Report your own numbers first,
then the artifact's.** Without this your claims will be discarded.

⛔ Do not state anywhere what you expect any count to be before you measure it.

## Operational constraints

- ⛔ Do not start `wolframscript`. ⛔ Do not modify or re-run either engine. ⛔ Never modify the working
  tree — copy to `/tmp` and ablate the copy.
- Wrap long runs in `timeout 1800`.
- Save every script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

The two engines' naming and representation divergences are the **output** of this instrument, not defects
in it. The py↔wl naming pass is separate work. Export-chain carried limits are out of scope. The baseline
process guard already fails and the exit code therefore carries no ablation signal — recorded
deliberately; only per-row residuals carry signal. The membership residual and the span residual do
different jobs and both are intended.
