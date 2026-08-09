# Independent physics review — S9/S10 export chain, fix round 4

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
re-derives, and writes the merged module the next step imports. This round repairs five properties of that
export. The physics was independently reproduced by an earlier leg and is **not** what this round touches.

⚠ **The specific failure this round exists to correct**: previous rounds were satisfied *literally* — a
named object was made true while the level above it stayed broken. **Check each repair at the level above
the obvious one.**

## What to check — by mutation, not by reading

1. **Can any failure leave a module on disk that a consumer cannot tell is stale?** Enumerate the ways a
   run can fail — during the run, at import of its own upstream, at parse, at write — and for each, plant
   a good module from an earlier run first, then fail that way, and report what is on disk afterwards.
   ⚠ The chain's own failure mode is: upstream engine fails and removes its module, so the downstream
   engine cannot import it. Test that one specifically.
2. **Can a consumer detect a stale or mismatched module?** It is claimed a module says what it was built
   from. Try to detect a module built from an upstream that has since changed.
3. **Can an authored conclusion reach the ledger through any carrier?** Try `str`, the sympy string
   object, a symbol whose name is a sentence, and anything else that carries text. Report every carrier
   that lands a record indistinguishable from a computed quantity. ⚠ Legitimate authored words exist —
   field names and branch outcomes. Report whether the honest ones survived and whether the population is
   emitted so a reader can bound it.
4. **Do the two assumption channels agree, and would a disagreement be reported?** An exported `Symbol`
   carries sympy assumptions; the engine separately refines under `Q.*` predicates. Change one channel
   without the other — in **both** directions — and report whether the run says anything.
   ⚠ Weakest change: widen a symbol's domain at its single declaration and see how much moves silently.
5. **Is the mapping the ledger hands a consumer editable?** Try key assignment, nested field assignment,
   insertion and deletion — not only assignment into a value.
6. **Does the unit-coverage measure vary?** A count that cannot be anything but one value carries no
   information. Change the coverage and show the emitted measure moves.
7. **Is the generated module a deterministic function of the run?** Several hash seeds, both engines, bytes.
8. **Did the derivation move?** Report what changed in each transcript beyond names, containers and export
   bookkeeping, and say whether any computed physics value moved.
9. **Did this round introduce anything new?** It is the fourth repair pass on these files. A round that
   breeds a fresh defect is the specific failure this question exists to find, and it has happened here
   three times.

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY** — change the *structure* of a load-bearing object, not a coefficient.
A coefficient rescale tests arithmetic; only a form change tests physics. Run a control at the unchanged
structure first and show your harness reproduces the artifact.

⛔⛔ **A test that passes on a weaker fix is not a test.** For each repair, construct the weakest change it
should reject and show whether it does.

⚠ **Watch for the recurring defect:** a residual whose two operands descend from one source is zero by
construction. **Four** such checks have been built and deleted on this project; one passed while the
change it policed was reverted entirely.

⚠ **A mutation can erase itself.** A previous leg injected a duplicate symbol as `0 * Symbol(...)`; the
term simplified away, the residual stayed zero and the injection was gone. **Verify your mutation is still
present in the object you think you mutated.**

Report with line numbers: anything typed where it could be read; a check whose expected value lives inside
the artifact it checks; an `assert` preceding the value it guards; a redundant emitted object.

### Write your own script first

⛔ Write your own verification script **before** opening the artifact; save it and its literal stdout to
named absolute paths and report them. **Without them your claims will be discarded.** A prose
re-derivation is the same defect class this rebuild exists to remove, relocated into the review.

## Do not read

Anything under `directives/` beginning `S10_fix_round4_decisions`, `S10_fix_round3_`, `S10_fix_round2_`,
`S10_export_chain_`, `S10_comparator_decisions`, or `S9_export_` / `S9_structural_`.
Anything under `/tmp/claude-1000/` or `/tmp/s9s10_export_r2_review/` — that is other people's evidence.
⭐ Produce your own.

## Operational constraints

- Ablate **copies** under `/tmp`. ⛔ Never modify the working tree.
- ⚠ **The S10 engine is slow** (~2-3 minutes per run). Wrap every run in `timeout 1800`; a hit is a failed
  ablation — report it and move on. ⛔ Never raise it further.
- ⛔ Do not start `wolframscript`. You may **read** the `.wl` engines; ⛔ never run or modify them.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

Carried limits, out of scope unless this round altered one: the field-dimension premise is unfalsifiable
within this engine; the third overwritten record is the difference of the other two; `field_dimension` is
an alias of `length_dimension`; `overwrites_upstream` is an author-set flag and nothing can check that a
value was derived rather than copied; the relational round-trip repair is inert; records holding equal
values for genuinely different questions are results, not duplication; cross-engine (py↔wl) tag naming and
the comparator are a separate pass.

⭐ **Also not a finding, and deliberately so:** the symbol-binding residual cannot see one quantity bound
under two names — it counts variants per name. Building a check for it would require a hand-written
name→quantity table, which this design exists to make permanent the deletion of. It is recorded as a
limit. ⛔ Do not report it; ⭐ **do** report any *new* instance of a split it would have to catch.
