# Independent physics review — cross-engine comparator, fix round 3

## Artifact

Untracked, in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_cross_engine_comparator_directive.md
research/pde_ledger_v3/directives/S10_comparator_fix_round3_directive.md
```

**You are one of two independent legs with this identical prompt.**

## What this is

Two engines compute S10 independently — SymPy and Wolfram — and each emits objects under the standard name
of the object, so comparison is a **join on the name**. This script performs that join.

⚠ **This round repairs a REGRESSION that the previous round introduced**, plus two properties that are
invisible on the committed input. ⇒ ⭐⭐ **you must CONSTRUCT the inputs that exhibit each one** — running on
the committed transcripts alone shows only that the current record is unaffected, which is already known.

⚠ The two transcripts are committed inputs and must **not** be regenerated.

## ⛔⛔ The method this round demands

For each repaired property: ⭐ **build the input that exhibits it, run the repaired instrument and show it
is right, then run the UNREPAIRED instrument on the same input and show it is not.** ⛔ Without both halves
you have shown nothing. ⭐ Recover the unrepaired behaviour however you like.

## What to check

1. **A name that denotes a mathematical constant must keep denoting it, in both engines.** The previous
   round captured such names as ordinary symbols; the damage was that two engines writing *the same*
   object were recorded as diverging, with a printed residual that reads as zero, and in one case a
   proposal that the imaginary unit and a physical symbol are one symbol under two spellings. ⭐ Verify
   constants survive **and** that a payload symbol which would otherwise be shadowed by a library callable
   or namespace still survives as itself. ⚠ Report anything in either population that lands on the wrong
   side, and **how the two populations are told apart** — ⛔ if the answer is a hardcoded list of names,
   that is a finding, because a list cannot be the mechanism.
2. **A proposed symbol pair must not bind two objects the transcripts already distinguish.** Construct a
   pair the transcripts refute (both names appearing elsewhere as distinct objects) and one they do not,
   and report which are proposed, which are marked refuted, and whether anything is silently dropped.
   ⚠ **A repointed physical symbol — a stiffness swapped for a density, one wavevector component for
   another — must not be proposable as a rename.** That is the recorded failure mode this exists to stop.
3. **Is the content decomposition honest?** The run now splits its content divergences into
   sub-populations. ⭐ Reproduce that split independently and report any row you would place differently.
   ⚠ A misclassification that shrinks the population carrying a genuine residual is the finding to look
   for.
4. **⛔ Did this round break an earlier one?** ⭐ **The nullspace residual must still pass BOTH halves** —
   fires on wrong subspace, partial rescale, rank loss; silent under whole-vector rescale, permutation and
   a general invertible change of basis. ⭐ **And round 2's repairs must still hold** — derivative
   dependence sets, the accounting property, the computed disagreement kind. ⚠ **A round that quietly
   breaks an earlier one has now happened once on this artifact.**
5. **Does it report agreement it did not compute?** Account for every shared name; the categories must sum.
   Derive the join with your own parser — ⚠ **a disagreement in either direction is a finding.**
6. **Account for every changed line** against the previous run.
7. **Conclusions.** Report any hardcoded expected value, count or tally in the script or directive; any
   typed verdict standing in for a computation; any `assert` preceding the value it guards.

## Required method

⛔ **Write your own comparison script BEFORE opening the artifact**, run it on the two committed
transcripts, and save it and its literal stdout to named absolute paths. **Report your own numbers first,
then the artifact's.** Without this your claims will be discarded.

⛔ Do not state anywhere what you expect any count to be before you measure it.

⚠ **If your instrument and the artifact disagree, adjudicate with a third route before reporting it** —
a previous leg found its own convention was the shortcut, not the artifact's.

## Operational constraints

- ⛔ Do not start `wolframscript`. ⛔ Do not modify or re-run either engine. ⛔ Never modify the working
  tree — copy to `/tmp` and ablate the copy.
- Wrap long runs in `timeout 1800`.
- Save every script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

The two engines' naming and representation divergences are the **output** of this instrument. The py↔wl
naming pass is separate work. Export-chain carried limits are out of scope. The baseline process guard
already fails and the exit code carries no ablation signal — recorded deliberately. The membership and
span residuals do different jobs and both are intended. The unparsed rows are one CAS's parser refusing a
Boolean containing a domain-membership assertion.
