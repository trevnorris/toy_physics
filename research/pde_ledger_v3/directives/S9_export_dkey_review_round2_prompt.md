# Independent physics review — S9 export re-key, round 2

## Artifact

The uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
git diff -- research/pde_ledger_v3/scripts/S9_exports.py
git diff -- research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out
```

plus `research/pde_ledger_v3/directives/S9_export_dkey_directive.md`.

All were written by Codex across two passes. You are one of two independent legs; the other has this
identical prompt and you will not see its output. **A previous review round exists. You are not to look
for it, and you are not being asked to confirm it.**

## What the artifact is

`scripts/S9_exports.py` is a generated flat `LEDGER` — one entry per object S9's `MAIN` package emits,
each `{display, value, dim?, class, step}`. The next step imports it, adds entries, overwrites what it
re-derives, and re-exports the merged dict. That next step computes the same physics at several values of
the brane spatial dimension `D` in one run, so **a ledger key must record the `D` at which its object was
computed**, and an object that is genuinely symbolic in `D` must not carry one.

## What to check

**The load-bearing claim: a key's component count is read out of the computation that produced the
object, and the general/fixed partition follows the structure of the computation rather than a judgement
typed beside it.**

Establish for yourself, by mutation and not by reading:

1. **Does the key follow the computation?** If the spatial algebra ran at a different component count,
   with the classification machinery untouched, what do the keys do? Run a control at the unchanged count
   first and show your harness contributes nothing — an ablation whose harness moves things is not
   evidence.
2. **Is the general/fixed partition correct, object by object?** An entry keyed as fixed-`D` that is
   actually general, or an unsuffixed entry that is actually a value at one `D`, is the defect this review
   exists to catch. Derive the partition yourself from the structure of the computation. ⛔ Free-symbol
   inspection has been measured to give the wrong answer in both directions — if you use it, say so and
   show what it gets wrong.
3. **Can a key's suffix still come from anywhere other than the recorded classification?** An internal
   variable name, a hand-typed string, a lookup table.
4. **The `.out` now carries new emitted lines.** Are they computed objects — both operands and a residual,
   emitted before any guard — or a conclusion? Can the guard fail? Show it failing.
5. **Did anything else move?** Payloads, class tallies, entry count, derivation, action, ansatz,
   assumptions, and every previously-emitted line of the `.out`. It is claimed nothing did.
6. **Did this round introduce anything new?** It is a repair of a previous round on the same files. A
   repair that breeds a fresh defect is the specific failure this question exists to find.

## Required method

Derive independently. **Ablate every load-bearing check and report its literal output.** Code-reading
alone has repeatedly missed real defects here — in the immediately preceding round, mutation found three
times what reading found.

⛔⛔ **A FORM ABLATION IS MANDATORY.** A coefficient rescale tests arithmetic; only a change of
**structure** tests the claim.

⛔⛔ **A TEST THAT PASSES ON A WEAKER FIX IS NOT A TEST.** Measured three times in this project: a half-fix
passed the new test, the whole suite and the full battery, and produced byte-identical output. For each
guard you find, construct the weakest change it should reject and show whether it does.

Probe for, and report with line numbers: a key, suffix, count or tally assembled from a literal beside the
computation rather than read out of it; a check whose expected value lives inside the artifact it checks;
a value verified using the predicate that produced it; an `assert` preceding the value it guards; a
conclusion emitted as an unconditional literal.

Ask of every claim in the directive: **which line computed this?** Give the line number or report it
uncomputed.

### Write your own script first

⛔ Write your own derivation/verification script **before** opening the artifact, and save both the script
and its literal stdout to named absolute paths. Report those paths. **Without them your derivation claims
will be discarded** — a prose re-derivation is the defect class this rebuild exists to remove, relocated
into the review.

### Reading order for the directive

Read the **engine and the generated export first**, form your own view of what the partition should be,
and only then read `S9_export_dkey_directive.md`. Reading it first anchors you to its framing, which is
under test. Quote both sides for every finding.

## Do not read

- `research/pde_ledger_v3/directives/S9_export_dkey_decisions.md`
- `research/pde_ledger_v3/directives/S9_export_dkey_fix_round1_decisions.md`
- `research/pde_ledger_v3/directives/S9_export_dkey_review_prompt.md`

These are the orchestrator's decision lists and the previous round's review brief. They state the
partition, the counts and the defects already found. Reading them turns this leg into a check that two
copies of an assumption agree. No Codex or review transcript is in the repository.

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

Three items are open on this step by decision and are out of scope: `wavevector_norm_dimension` names
`dim(|k|)` but holds `dim(k·k)`; the placeholder-naming class has eight members; `q_dimension` is unpinned
inside SymPy. Do not re-report them — unless this round changed one, which is in scope.
