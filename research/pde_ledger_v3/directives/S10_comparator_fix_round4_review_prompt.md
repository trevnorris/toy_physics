# Independent physics review — cross-engine comparator, fix round 4 (FINAL)

## Artifact

Untracked, in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_cross_engine_comparator_directive.md
research/pde_ledger_v3/directives/S10_comparator_fix_round4_directive.md
```

**You are one of two independent legs with this identical prompt.**

## What this is

Two engines compute S10 independently — SymPy and Wolfram — and each emits objects under the standard name
of the object, so comparison is a **join on the name**. This script performs that join.

⚠ **This is the LAST round on this artifact.** It fixes one classification defect and **deletes** two guard
terms that could not be nonzero. It adds nothing. ⇒ ⭐ **an addition is itself a finding**, and so is any
change to a verdict on the committed input.

⚠ The two transcripts are committed inputs and must **not** be regenerated.

## What to check

1. **⛔⛔ The classification defect, on the REAL transcripts.** Repoint one physical symbol to another that
   both engines distinguish elsewhere — a stiffness to a density, one wavevector component to another —
   and report which population the row lands in. ⚠ Previously such a row was recorded as a **naming**
   difference while the content population stayed unmoved, in the same run that printed the pair as
   transcript-refuted. ⭐ Show both halves: the repaired instrument puts it where it belongs, the
   unrepaired one does not.
   ⭐ **And show a GENUINE spelling difference still classifies as naming-only** — ⚠ a fix that pushes
   everything into content is no better than the defect.
2. **The two deleted guard terms.** Reconstruct each and try to make it fire. ⭐ Show by running that
   neither could have been nonzero. ⚠ Then confirm the property they named is still **readable** — the
   populations and totals must still be emitted so a reader can do the arithmetic.
3. **⛔ Did this round break an earlier one?** ⭐ The nullspace residual must still pass **both** halves —
   fires on wrong subspace, partial rescale, rank loss; silent under whole-vector rescale, permutation and
   a general invertible change of basis. ⭐ The constant/shadow split, derivative dependence sets, the
   accounting behaviour under duplicated names, and the computed disagreement kind must all still hold.
   ⚠ **A round quietly breaking an earlier one has happened once on this artifact.**
4. **Does it report agreement it did not compute?** Account for every shared name. Derive the join with
   your own parser — ⚠ **a disagreement in either direction is a finding.**
5. **Account for every changed line** against the previous run, and confirm **no verdict moved** on the
   committed input.
6. **Conclusions.** Any hardcoded expected value, count or tally in script or directive; any typed verdict
   standing in for a computation; any `assert` preceding the value it guards.

## Required method

⛔ **Write your own comparison script BEFORE opening the artifact**, run it on the two committed
transcripts, and save it and its literal stdout to named absolute paths. **Report your own numbers first.**
Without this your claims will be discarded.

⛔ Do not state anywhere what you expect any count to be before you measure it.

⚠ **If your instrument and the artifact disagree, adjudicate with a third route before reporting it** — a
previous leg found its own convention was the shortcut, not the artifact's.

## Operational constraints

- ⛔ Do not start `wolframscript`. ⛔ Do not modify or re-run either engine. ⛔ Never modify the working
  tree — copy to `/tmp` and ablate the copy.
- Wrap long runs in `timeout 1800`.
- Save every script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings — ⭐ these are RECORDED limits and are deliberately unrepaired

⛔ Do not report these; ⭐ **do** report a new instance or any change to one:
the constant/shadow rule is applied to one engine's payload parser only, so a library-callable spelling in
the other engine's transcript would fail **loudly** as not-computed rather than agreeing falsely; two
sub-population selectors test the emitted name's suffix rather than the residual; the process guard already
fails at baseline so the exit code carries no ablation signal — only per-row residuals do; the membership
and span residuals do different jobs and both are intended; the unparsed rows are one CAS's parser refusing
a Boolean containing a domain-membership assertion. Engine naming and representation divergences are the
**output** of this instrument. The py↔wl naming pass is separate work.
