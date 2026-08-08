# Independent physics review — S9 export re-key

## Artifact

The uncommitted working-tree change in `/var/projects/toy_physics`, on branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
git diff -- research/pde_ledger_v3/scripts/S9_exports.py
```

plus the new file `research/pde_ledger_v3/directives/S9_export_dkey_directive.md`.

All three were written by Codex in one pass. You are one of two independent legs; the other has been given
this identical prompt and you will not see its output.

## What the change is for

`scripts/S9_exports.py` is a generated flat `LEDGER` — one entry per object S9's `MAIN` package emits,
each `{display, value, dim?, class, step}`. The next step (S10) will import it, add its own entries,
overwrite what it re-derives, and re-export the merged dict. S10 computes the same physics at several
values of the brane spatial dimension `D` in one run, so a ledger key must record **the `D` at which the
object was computed**, or S10's results have nowhere to land.

This change re-keys S9's export accordingly. It is supposed to move **keys only**.

## What to check

**The load-bearing claim is that a ledger key's `D` is READ OUT OF the computation that produced the
object, and is not a literal typed beside it.** S9's spectrum algebra runs at a fixed number of spatial
components; its dimension solve is symbolic in `D`. A key that says `d3` because the string `"d3"` was
concatenated somewhere is indistinguishable, on this run, from a key that says `d3` because the
computation had three components — and the second is the only one that survives contact with S10 and with
the 4D bulk sector.

Also check, independently:

- that the re-key moved **no value, dimension, class, step or display payload**;
- that the partition between fixed-`D` and general-`D` entries is correct object by object — an entry
  keyed `_d3` that is actually general in `D`, or an unsuffixed entry that is actually a `D = 3` value, is
  the defect this review exists to catch;
- that the existing exact-`srepr` round-trip check still covers every record under its new key, and would
  fail if it did not;
- whether anything in the derivation, the action, the ansatz, the assumptions or the emitted `.out` moved.
  It is claimed nothing did.

## Required method

Derive independently. **Ablate every load-bearing check and report its literal output.** Code-reading
alone has repeatedly missed real defects in this repository.

⛔⛔ **A FORM ABLATION IS MANDATORY, NOT OPTIONAL.** A coefficient rescale tests arithmetic; only a change
of **structure** tests the claim. For this artifact the decisive structural question is: *does the key
follow the computation?* Construct the ablation that answers it and report the literal diff of the
generated export. If you conclude the ablation cannot be built, say so and say why — that is itself a
finding about the claim.

Probe for, and report with line numbers:

- a key, suffix or count assembled from a literal beside the computation rather than read out of it;
- a check whose expected value lives inside the artifact it checks;
- a value verified using the predicate or definition that produced it;
- an `assert` that precedes the value it guards — a perturbation strong enough to flip it kills the
  process, so you see only PASS-or-crash;
- a conclusion emitted as an unconditional literal.

Ask of every claim in the authored directive: **which line computed this?** Give the line number, or
report it as uncomputed.

### Write your own script first

⛔ **Write your own derivation/verification script BEFORE opening the artifact**, and save both the script
and its literal stdout to named absolute paths. Report those paths. **Without them your derivation claims
will be discarded** — a prose re-derivation is the same defect class this rebuild exists to remove,
relocated into the review.

### Reading order for the authored directive

Read the **engine and the committed export first** and form your own view of what the key policy should
be. Only then read `S9_export_dkey_directive.md`. Reading it first anchors you to its framing, which is
part of what is under test. Quote both sides for every finding.

## Do not read

- `research/pde_ledger_v3/directives/S9_export_dkey_decisions.md` — the orchestrator's decision list. It
  states the partition and the counts. Reading it makes this leg a check that two copies of an assumption
  agree.
- Any Codex or review transcript. The build transcript has been moved outside the repository.

## Operational constraints

- Copy anything you ablate to `/tmp` and ablate the **copy**. ⛔ Never modify the working tree — it holds
  the uncommitted change under review and there is no committed baseline to restore it from.
- Wrap every script run in `timeout 600`. A 600s hit is a failed ablation: report it and move on. ⛔ Never
  raise the timeout.
- This artifact spawns no CAS kernel; if you find yourself starting `wolframscript`, stop — the Wolfram
  engine is explicitly not part of this change.
- Save every ablation script AND its literal stdout to named absolute paths, and report those paths.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong. Do
not report "the script would be wrong on a different input."

## One thing that is NOT a finding

Three items are open on this step by decision and are explicitly out of scope: `wavevector_norm_dimension`
names `dim(|k|)` but holds `dim(k·k)`; the placeholder-naming class has eight members; `q_dimension` is
unpinned inside SymPy. Do not re-report them. If the re-key **changed** any of them, that is in scope.
