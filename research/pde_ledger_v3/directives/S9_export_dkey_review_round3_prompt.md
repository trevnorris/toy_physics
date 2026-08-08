# Independent physics review — S9 export re-key, round 3

## Artifact

The uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
git diff -- research/pde_ledger_v3/scripts/S9_exports.py
git diff -- research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out
```

plus `research/pde_ledger_v3/directives/S9_export_dkey_directive.md`.

Written by Codex across three passes. You are one of two independent legs; the other has this identical
prompt. **Previous review rounds exist. Do not look for them. You are not being asked to confirm anyone.**

## What the artifact is

`scripts/S9_exports.py` is a generated flat `LEDGER` — one entry per object S9's `MAIN` package emits,
each `{display, value, dim?, class, step}`. The next step imports it, adds entries, overwrites what it
re-derives, and re-exports the merged dict. That step computes the same physics at several values of the
brane spatial dimension `D` in one run.

**A ledger key records the construction that produced the object** — an object built from the
fixed-component action carries that component count; an object built from unit algebra or from the
dimension solve's symbolic output does not. ⚠ The key deliberately does **not** record whether the
object's *value* is sensitive to `D`: that is not computable from a run at a single component count, and
claiming it would be an assertion the engine cannot back.

## What to check

Establish each of these for yourself, **by mutation and not by reading**:

1. **Is every object's classification a consequence of how it is constructed?** Anything classified by a
   hand-typed name, a lookup, a special case or a bypass is the defect this review exists to catch.
   Mutate a group's tag and see exactly which objects follow it; mutate the component count the
   computation runs at and see which keys move.
2. **Is the partition correct, object by object?** An object carrying a component count that no
   construction gave it, or an object missing one that its construction implies. Derive the partition
   yourself from the structure. ⛔ Free-symbol inspection gives the wrong answer in both directions here —
   if you use it, show what it gets wrong.
3. **The `.out` carries emitted lines describing the partition.** Are they computed objects that tell a
   reader something, or do they create an appearance of a check that is not happening? Read them as a
   sceptical consumer of the record would. If a guard exists, show it failing; if none exists, say what a
   reader could and could not conclude from what is printed.
4. **Some objects are newly emitted per package this round.** Confirm whether their values are what the
   engine computes, and whether any previously-emitted value moved.
5. **Did anything else move?** Payloads, class tallies, entry count, derivation, action, ansatz,
   assumptions, and every previously-emitted `.out` line. It is claimed nothing did.
6. **Did this round introduce anything new?** It is the second repair on the same files. A repair that
   breeds a fresh defect is the specific failure this question exists to find, and it has already happened
   once on this artifact.

## Required method

Derive independently. **Ablate every load-bearing check and report its literal output.**

⛔⛔ **A FORM ABLATION IS MANDATORY.** A coefficient rescale tests arithmetic; only a change of
**structure** tests the claim. Run a control at the unchanged structure first and show your harness
reproduces the current export — an ablation whose harness moves things is not evidence.

⛔⛔ **A TEST THAT PASSES ON A WEAKER FIX IS NOT A TEST.** For each guard you find, construct the weakest
change it should reject and show whether it does.

Probe for, and report with line numbers: a key, suffix, count or tally assembled from a literal beside the
computation rather than read out of it; a check whose expected value lives inside the artifact it checks;
a value verified using the predicate that produced it; an `assert` preceding the value it guards; a
conclusion emitted as an unconditional literal; **an emitted object that is redundant with another
emitted object**.

Ask of every claim in the directive: **which line computed this?** Give the line number or report it
uncomputed.

### Write your own script first

⛔ Write your own derivation/verification script **before** opening the artifact, and save both the script
and its literal stdout to named absolute paths. Report those paths. **Without them your derivation claims
will be discarded.**

### Reading order for the directive

Read the **engine and the generated export first**, form your own view of what the partition should be,
and only then read `S9_export_dkey_directive.md`. Quote both sides for every finding.

## Do not read

- `research/pde_ledger_v3/directives/S9_export_dkey_decisions.md`
- `research/pde_ledger_v3/directives/S9_export_dkey_fix_round1_decisions.md`
- `research/pde_ledger_v3/directives/S9_export_dkey_fix_round2_decisions.md`
- `research/pde_ledger_v3/directives/S9_export_dkey_review_prompt.md`
- `research/pde_ledger_v3/directives/S9_export_dkey_review_round2_prompt.md`

These are the orchestrator's decision lists and earlier review briefs. They state the partition, the
counts and the defects already found. Reading them turns this leg into a check that two copies of an
assumption agree. No Codex or review transcript is in the repository.

## Operational constraints

- Copy anything you ablate to `/tmp` and ablate the **copy**. ⛔ Never modify the working tree — it holds
  the uncommitted change under review and there is no committed baseline to restore it from.
- Wrap every script run in `timeout 600`. A 600s hit is a failed ablation: report it and move on. ⛔ Never
  raise the timeout.
- This artifact spawns no CAS kernel. If you find yourself starting `wolframscript`, stop.
- Save every ablation script AND its literal stdout to named absolute paths, and report those paths.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong. Do
not report "the script would be wrong on a different input."

## Not findings

Out of scope by decision and unchanged: `wavevector_norm_dimension` names `dim(|k|)` but holds `dim(k·k)`;
the placeholder-naming class has eight members; `q_dimension` is unpinned inside SymPy. Do not re-report
them unless this round changed one.
