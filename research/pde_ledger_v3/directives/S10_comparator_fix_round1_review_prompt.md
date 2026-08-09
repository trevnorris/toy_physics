# Independent physics review — cross-engine comparator, fix round 1

## Artifact

Untracked, in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_cross_engine_comparator_directive.md
research/pde_ledger_v3/directives/S10_comparator_fix_round1_directive.md
```

**You are one of two independent legs with this identical prompt.**

## What this is

Two engines compute S10 independently — SymPy and Wolfram — and each emits objects under the standard name
of the object, so comparison is a **join on the name**. This script performs that join. A previous review
found that a class of object was being **excluded** from comparison rather than compared; this round
replaces the exclusion with a residual.

⚠ The two transcripts are committed inputs and must **not** be regenerated:
`scripts/out/S10_brane_mode_spectrum_sympy_audit.out` and
`mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`.

## ⛔⛔ The acceptance, and it has TWO halves that must both hold

A residual for an object with a genuine representational freedom must be **invariant under exactly that
freedom, and sensitive to everything else.** ⚠ A check that fires on both is no better than one that fires
on neither.

⭐ **Half one — it must fire.** For the object class this round changed, construct mutations that are *not*
legitimate representational changes and show the residual moves. ⭐ Include at least one where the
substituted value does not satisfy the defining property of the object at all.

⭐ **Half two — it must not fire.** Construct a *legitimate* representational change of the same object and
show the residual does **not** move. ⛔ A residual that rejects a valid alternative representation is a
false alarm and is a finding.

⭐ **Verify against the defining property yourself**, from objects in the transcripts, rather than trusting
the artifact's notion of what the object is.

## What else to check — by mutation, not by reading

1. **Does it report agreement it did not compute?** A row it could not parse, could not normalise, or
   silently skipped **reads exactly like a row that agreed.** ⭐ Account for every shared name: the
   categories must sum to the number of names both engines emit.
2. **Derive the join yourself** — your own parser, your own counts — and compare with what the artifact
   reports. ⛔ Do not use the artifact's parser to check the artifact's counts. ⚠ **A disagreement in
   either direction is a finding.**
3. **Are the reported agreements real?** Verify a sample by your own route, including at least one pair
   written in visibly different form (expanded vs factored).
4. **Is the new decomposition honest?** The run now classifies its agreements and its disagreements into
   sub-populations. ⭐ Reproduce those classifications independently and report any row you would classify
   differently. ⚠ A misclassification that flatters the result is the finding to look for.
5. **Correspondences.** Two are sanctioned: a mechanical symbol-spelling rule and container/power syntax,
   plus a documented derivative-syntax parse. ⛔ Any other absorbed difference — a special case, a
   hardcoded pair, a normalisation that discards structure — is a finding. ⚠ Check whether sequence
   **order** is preserved.
6. **Did this round break anything it was not meant to touch?** Compare the run against the previous one
   and account for every changed line.
7. **Conclusions.** Per rule 2 a script prints computed objects, then guards. Report any hardcoded expected
   value, any typed verdict standing in for a computation, any `assert` preceding the value it guards.

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
in it — report them as inventory. The py↔wl naming pass is separate work. Export-chain carried limits are
out of scope. The baseline process guard is already failing and the exit code therefore carries no
ablation signal — that is recorded deliberately; only per-row residuals carry signal.
