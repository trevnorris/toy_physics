# Independent physics review — the cross-engine comparator (F-2)

## Artifact

Untracked, in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py
research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
research/pde_ledger_v3/directives/S10_cross_engine_comparator_directive.md
```

**You are one of two independent legs with this identical prompt.**

## What this is

Two engines compute S10 independently — a SymPy one and a Wolfram one — and each emits objects under **the
standard name of the object**, so cross-engine comparison should be a **join on the name**. Until now
nothing in the repository performed that join. A review leg once re-pointed a single standard name so the
light cone became an `ω²` instead of a speed squared, and **every check in the repository passed**. This
script is the instrument that is supposed to catch that.

⚠ The two transcripts are committed and are **not** to be regenerated:
`scripts/out/S10_brane_mode_spectrum_sympy_audit.out` and
`mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`.

## ⛔⛔ The first question, and it is the whole point

⭐⭐ **Does the comparator report agreement it did not compute?** A row it could not parse, could not
normalise, or silently skipped **reads exactly like a row that agreed**. That is the failure mode this
instrument exists to remove, and it is the one it is most likely to have.

⭐ **Account for every shared name.** Sum the reported categories and check the total equals the number of
names both engines emit. ⛔ A name that falls out of the join without being counted is a finding.

## What to check — by mutation, not by reading

1. **Derive the join yourself.** Parse both transcripts with your own code, count the names each engine
   emits and how many are shared, and compare against what the artifact reports. ⛔ Do not use the
   artifact's own parser to check the artifact's own counts.
2. **Is agreement real?** For a sample of rows the artifact calls agreeing, verify the two payloads are
   the same mathematical object by your own route. ⚠ Include at least one where the two engines write the
   object in visibly different form — one expanded and one factored, say — because a comparator that only
   ever matches text is not comparing mathematics.
3. **Can it fail?** Re-point a standard name in one transcript to a neighbouring object and show the
   residual move. Then do the harder one: change a payload's **value** slightly and confirm it is caught.
   ⛔ A comparator that cannot fail is decoration, and a comparator that never compares has already passed
   this project's acceptance criterion once.
4. **What correspondences does it apply, and are any of them lying?** Two are sanctioned: a mechanical
   symbol-spelling rule, and container/power syntax. ⛔ Any other absorbed difference — a special case, a
   hardcoded pair, a normalisation that discards structure — is a finding. ⚠ Check specifically whether
   list **order** is normalised, and whether that is safe for every object it is applied to.
5. **Non-canonical objects.** A nullspace basis is not canonical; comparing one as though it were produces
   disagreements on representation alone. Report how the artifact handles such objects and whether its
   choice is stated.
6. **One-sided names.** Names only one engine emits are inventory, not failure. Report whether the counts
   are emitted and whether anything is silently dropped.
7. **Does it state conclusions?** Per this project's rule 2, a script prints computed objects and never
   states conclusions. Report any typed verdict, any hardcoded expected value, and any `assert` that
   precedes the value it guards.

## Required method

⛔ **Write your own comparison script BEFORE opening the artifact**, run it on the two committed
transcripts, and save it and its literal stdout to named absolute paths. **Report your own numbers first,
then the artifact's.** Without this your claims will be discarded.

⚠ **A disagreement between your tally and the artifact's is a finding, in either direction** — that is why
two independent implementations exist.

⛔ Do not state anywhere what you expect any count to be before you measure it.

## Operational constraints

- ⛔ Do not start `wolframscript`. ⛔ Do not modify or re-run either engine. ⛔ Never modify the working
  tree — copy to `/tmp` and ablate the copy.
- Wrap long runs in `timeout 1800`.
- Save every script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

The two engines' underlying disagreements about naming are the **output** of this instrument, not defects
in it — report them as inventory. The py↔wl naming pass itself is a separate piece of work. Carried limits
of the export chain (authored-text carriers, assumption-channel coverage, parse-failure staleness,
`assert`-based guards, digests omitting the CAS version) are out of scope.
