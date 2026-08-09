# Independent physics review — S9/S10 export chain, round 5 (a REMOVAL round)

## Artifact

Uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/
git status --short research/pde_ledger_v3/
```

⚠ The committed baseline `HEAD` is **two rounds back**, so this diff contains an earlier round's work as
well. ⭐ Review **the working tree as it stands**; you are not asked to separate the rounds.

**You are one of two independent legs with this identical prompt. Earlier rounds exist; do not look for
them.**

## What this is

Each step's SymPy engine imports the previous step's flat `LEDGER`, binds its objects rather than
re-declaring them, adds its own derivations, overwrites what it re-derives, and writes the merged module
the next step imports.

⭐⭐ **This round was supposed to ADD NOTHING.** It deletes two residuals that could not be nonzero, fixes
one regression, and records limits in prose. ⇒ ⭐ **an addition is itself a finding.**

## ⛔⛔ THE CONTROL IS THE FIRST THING, AND IT IS NOT OPTIONAL

⚠ **A previous leg's entire measurement was invalidated because the engines changed underneath it.** Its
control run, its FORM ablation and every weakest-change test were established against bytes that no longer
exist on disk.

⭐ **Before anything else: run both engines from a clean sandbox and report whether they reproduce, byte
for byte, the two generated modules and the two committed transcripts in the working tree.** ⛔ If they do
not, that is the finding and you should report it before going further.

## What to check — by mutation, not by reading

1. **Did anything get added?** Enumerate every emitted tag, every guard and every residual, before and
   after. ⭐ Report additions, ⛔ not just removals.
2. **Could the two deleted residuals have been nonzero?** Reconstruct each and try to make it fire.
   ⭐ Show by running, ⛔ not by argument.
3. **Was the evidence kept when the residual went?** The digest record and the value-kind field were meant
   to survive. Report whether they did and whether they still carry what they claim.
4. **Two properties must hold AT ONCE:** importing either engine changes no file on disk, **and** the
   chain's own failure mode — upstream module missing — leaves no stale downstream module. Test both, and
   test them **together**. ⚠ A fix for either one alone has already broken the other on this project.
5. **Did the fix for (4) change anything else?** One engine's export write was also moved. Report what that
   changes about running the engine as a script, importing it, and anything else you can measure.
6. **Enumerate failure modes again**: during the run, at import of its own upstream, at parse, at write,
   at a tripped guard. Plant a good module from an earlier run first, then fail each way, and report what
   is on disk afterwards.
7. **Did the derivation move?** ⭐ Expect the transcript diff to be **deletions and nothing else.** Report
   anything else that moved, and whether any computed physics value changed.
8. **Is the generated module still a deterministic function of the run?** Several hash seeds, both engines.
9. **Re-establish the cross-step measurement on the CURRENT bytes.** It is claimed the S9→S10 overwrite
   residual is a genuine independent measurement rather than zero by construction. ⛔ Do not take that on
   trust from any earlier round — **a FORM change to a load-bearing object must drive it nonzero, and a
   coefficient rescale must not.**

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY** — change the *structure* of a load-bearing object, not a coefficient.

⛔⛔ **A test that passes on a weaker fix is not a test.**

⚠ **A mutation can erase itself.** Verify your mutation is still present in the object you think you
mutated — a leg lost an injection to `simplify` and got a false pass.

⚠ **Watch for the recurring defect:** a residual whose two operands descend from one source is zero by
construction. **Six** have now been built on this project; two were deleted in this very round.

### Write your own script first

⛔ Write your own verification script **before** opening the artifact; save it and its literal stdout to
named absolute paths and report them. **Without them your claims will be discarded.**

## Do not read

Anything under `directives/` beginning `S10_fix_round5_decisions`, `S10_fix_round4_`, `S10_fix_round3_`,
`S10_fix_round2_`, `S10_export_chain_`, or `S10_comparator_decisions`.
Anything under `/tmp/claude-1000/` or `/tmp/s9s10_export_r2_review/` — other people's evidence.

## Operational constraints

- Ablate **copies** under `/tmp`. ⛔ Never modify the working tree.
- ⚠ **The S10 engine is slow** (~2-3 min). Wrap every run in `timeout 1800`. ⛔ Never raise it.
- ⛔ Do not start `wolframscript`. You may **read** the `.wl` engines; ⛔ never run or modify them.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings — ⭐ these are RECORDED limits and this round deliberately did not fix them

⛔ Do not report any of these; ⭐ **do** report a *new* instance or any change to one:
`value_kind` marks one carrier of five (`Symbol`/`Function`/nested `Str`/`Eq(Str,Str)`/`display`); the
authored-word count does not bound the population; the assumption channel is one-way, `MAIN`-only and
scores zero when emptied; a parse error leaves a stale module; guards are `assert` so `PYTHONOPTIMIZE`
strips them; the digests omit the sympy and Python versions; `symbol_binding_residual` cannot see one
quantity under two names; `overwrites_upstream` is an author-set flag; `field_dimension` aliases
`length_dimension`; the relational round-trip repair is inert; equal values for genuinely different
questions are results, not duplication; py↔wl naming and the comparator are a separate pass.
