# Independent review — the S10 step record

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/steps/S10_two_transverse_photons.md`

**You are one of two independent legs with this identical prompt.** ⛔ Read-only: do not edit the working
tree. Verify with `git status` before you report.

## What this is

The permanent record for S10 — the step establishing light's transverse mode count from the brane's own
dimension. It was just rewritten, because the instrument it was originally written against has been
deleted from the repository.

⚠ **Its job is to state what was MEASURED and what may NOT be claimed on the strength of it.** A defect
here is a false claim about the physics entering the ledger's permanent record.

## ⛔⛔ Required reading order — this is the whole method

**Read the sources of truth FIRST and form your own view of what S10 establishes and what it does not.
⛔ Only THEN open the record.** Reading the record first anchors you to its framing, which is the thing
under test. ⭐ For a document, blindness comes from reading order, and it is load-bearing.

Sources of truth:

- `research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out`
- `research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`
- `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out` — the cross-engine join
- `research/pde_ledger_v3/scripts/S9_exports.py`, `research/pde_ledger_v3/scripts/S10_exports.py`
- both engine sources; `research/pde_ledger_v3/directives/S10_SHARED_PHYSICS.md`
- `research/pde_ledger_v3/steps/S10_PREREGISTERED_PREDICTION.md`

Write down what you think the artifacts establish **before** you open the record. Report that view.

## What to check

1. ⛔⛔ **Can any sentence in this record be honestly quoted to OVERSTATE S10's result?** ⭐ **This is the
   review.** Go looking for the quotable sentence — a headline a reader could lift without its
   qualification, an agreement presented as stronger than the comparison behind it, a limit stated so
   softly it reads as closed.
2. ⛔ **Every claim: which artifact establishes it, at which locus?** ⚠ **A resolving path is not a
   source** — open the cited lines and check they contain an operand supporting the claim. Report any
   claim with no source, or a source that does not support it.
3. ⛔ **Cross-engine claims: does the join actually contain that family of names?** A claim of agreement
   for objects the two engines never joined is the defect class this rebuild exists to remove.
4. ⛔ **Export-chain claims: are they bounded to the objects actually compared**, or do they escalate to
   the imported ledger, the action, or the spectrum?
5. ⛔⛔ **WHAT WAS DROPPED.** The previous record is at `git show HEAD:research/pde_ledger_v3/steps/S10_two_transverse_photons.md`
   and was substantially longer. ⭐ **Diff them and report every limit, qualification or finding that
   disappeared and is still true of the committed artifacts.** ⚠ A shortened record is the easiest place
   to lose a limit a reader needs, and this is the check most likely to find something.
6. ⛔⛔ **DOES ANY SENTENCE NAME A CONTROL BY A SUBSTANCE INSTEAD OF BY ITS FUNCTIONAL?** ⚠ The controls
   differ from `MAIN` only in the **stiffness functional** the shared spec assigns them. ⛔ Calling one
   *"an ordinary elastic solid"* — or any other medium — imports a substance claim the computation does
   not make, and specifically grants the medium a **shear modulus**. ⚠ **This model's bulk is a
   superfluid: zero shear modulus by definition.** The brane's shear rigidity is an **OPEN requirement**
   (`SUBSTRATE_REQUIREMENTS.md` `R-S1-02`) resting on one external source with no executing script.
   ⭐ Report any sentence that treats it as settled, or that borrows intuition from a medium we do not
   claim to have. ⚠ The previous record did this; check it has not returned.
7. ⛔ **Does the record state what it CANNOT claim, specifically enough to bite?** Common-mode (both
   engines take the action from one spec), the generic-versus-stratum distinction, anything resting on one
   engine's own objects.
7. **Numbers.** Any figure you cannot reproduce from an artifact. ⚠ Report any that is right as well as
   any that is wrong — a number that happens to be correct but is unsourced is still a defect.

## Required method

⭐ **For every claim you make about what an artifact contains, run the command and quote its literal
output.** Save any script you write and its stdout to named absolute paths and report them.
⛔ **A prose assertion about what a file contains will be discarded.**

⛔ Do not state anywhere what you expect any count or value to be before you measure it.

⚠ If your measurement and the record disagree, adjudicate with a second route before reporting it.

## Operational constraints

- ⛔ **Do not start `wolframscript`.** The Wolfram transcript is committed; read it.
- ⛔ Do not run or modify either engine or the comparator in the tree. Copy to `/tmp` if you must run
  anything, and never modify the working tree.
- ⭐ A long run that is printing is fine; a long run that has gone **silent** is a failed measurement —
  kill it and report, ⛔ do not raise its budget.

## Physics filter

Report a finding only if it catches a way **the physics, or the ledger's record of it, could be wrong or
could be honestly quoted to overstate the result.** ⛔ Not style, not formatting, not length.

⚠ **"Nothing survives the filter" is a permitted answer and is weak evidence on its own** — if that is
your conclusion, say what you attacked and how, so the strength of the pass is visible.
