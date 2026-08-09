# Independent physics review — S9/S10 export chain, fix round 3

## Artifact

Uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/
git status --short research/pde_ledger_v3/
```

Both SymPy engines, `extract_knob_inventory.py`, both **generated** modules `scripts/S9_exports.py` and
`scripts/S10_exports.py`, both regenerated transcripts under `scripts/out/`, and a new directive under
`directives/`.

Codex, one build pass. **You are one of two independent legs with this identical prompt. Earlier rounds
exist; do not look for them.**

## What this is

Each step's SymPy engine imports the previous step's flat `LEDGER`, binds its objects rather than
re-declaring them, adds its own derivations across several component counts, overwrites what it
re-derives, and writes the merged module the next step imports. This round repairs seven properties of
that export. The physics was independently reproduced by an earlier leg and is **not** what this round
touches.

## What to check — by mutation, not by reading

1. **Does one quantity have exactly one name across the chain?** The central claim of this round is that a
   consumer binding a quantity from the merged ledger reaches **every** record that uses it. Test it as
   that consumer: bind each exported symbol and substitute, and report any record that does not move but
   uses the same quantity under another spelling. ⚠ Two spellings of one quantity produce a binding check
   that reads zero **by construction** — a name→object check can never see a quantity→name split. Look
   for the split directly.
2. **Were assumptions changed to make a rename work?** A rename that quietly widens or narrows a symbol's
   assumptions is a physics change wearing a spelling change's clothes. Compare the assumption set of
   every renamed symbol before and after, and report any that moved and what it moves.
3. **Is an authored word distinguishable from a bindable quantity?** Emit a new authored word and see
   whether it becomes a ledger record. Report any record whose whole value is a word, and whether a
   difference across component counts over such a record can be anything but zero.
4. **Can a consumer recover what units an exported object carries?** It is claimed this is discoverable.
   Take the consumer's position and try, without being told the spelling of any sibling key.
5. **Does the traversability guard reach inside a value?** Construct a value that is a CAS object at the
   top and carries a foreign container beneath it, and report whether the run rejects it.
6. **Does a failed run leave a generated module behind?** Make a guard fail on **each** engine and look on
   disk. Also check a stale module from a previous successful run cannot survive a later failure.
7. **Can a consumer mutate the ledger in place?** Try to assign into an exported value.
8. **Does the overwrite row still carry what two steps independently agreed on?** Report what the merged
   entry records about its own provenance, and whether that survives the merge or only the transcript.
9. **Is the generated module a deterministic function of the run?** Run under several hash seeds and
   compare bytes. ⚠ Do this for **both** engines — a previous round proved it for one and assumed the other.
10. **Did the derivation move?** Report what changed in each transcript beyond names, containers and
    export bookkeeping, and say whether any computed physics value moved.
11. **Did this round introduce anything new?** It is the third repair pass on these files. A round that
    breeds a fresh defect is the specific failure this question exists to find, and it has happened here.

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY** — change the *structure* of a load-bearing object, not a coefficient.
A coefficient rescale tests arithmetic; only a form change tests physics. Run a control at the unchanged
structure first and show your harness reproduces the artifact.

⛔⛔ **A test that passes on a weaker fix is not a test.** For each of the seven repairs, construct the
weakest change it should reject and show whether it does.

⚠ **Watch for the recurring defect:** a residual whose two operands descend from one source is zero by
construction. **Four** such checks have been built and deleted on this project; one passed while the
change it policed was reverted entirely.

⚠ **A mutation can erase itself.** A previous leg injected a duplicate symbol as `0 * Symbol(...)` and the
term simplified away, so the residual stayed zero and the injection was gone — a false pass it caught
only by re-running with a non-simplifying carrier. Verify your mutation is still present in the object
you think you mutated.

Report with line numbers: anything typed where it could be read; a check whose expected value lives inside
the artifact it checks; an `assert` preceding the value it guards; a redundant emitted object.

### Write your own script first

⛔ Write your own verification script **before** opening the artifact; save it and its literal stdout to
named absolute paths and report them. **Without them your claims will be discarded.** A prose
re-derivation is the same defect class this rebuild exists to remove, relocated into the review.

## Do not read

Anything under `directives/` beginning `S10_fix_round3_decisions`, `S10_fix_round2_`,
`S10_export_chain_`, `S10_comparator_decisions`, or `S9_export_` / `S9_structural_`.
Anything under `/tmp/claude-1000/` or `/tmp/s9s10_export_r2_review/` — that is other people's evidence.
⭐ Produce your own.

## Operational constraints

- Ablate **copies** under `/tmp`. ⛔ Never modify the working tree.
- ⚠ **The S10 engine is slow** (~2-3 minutes per run). Wrap every run in `timeout 1800`; a hit is a failed
  ablation — report it and move on. ⛔ Never raise it further.
- ⛔ Do not start `wolframscript`. You may **read** the `.wl` engine as evidence about naming; ⛔ never run
  or modify it.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

Carried limits, out of scope unless this round altered one: the field-dimension premise is unfalsifiable
within this engine now that the registry comparison is gone; the third overwritten record is the
difference of the other two and adds no falsification power; `field_dimension` is an alias of
`length_dimension`; `overwrites_upstream` is an author-set flag and nothing can check that a value was
derived rather than copied; the relational round-trip repair is inert because no allowed stratum occurs on
the export path; records holding equal values for genuinely different questions are results, not
duplication; cross-engine (py↔wl) tag naming and the comparator are a separate pass.
