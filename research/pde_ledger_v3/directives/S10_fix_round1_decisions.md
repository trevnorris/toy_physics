# S10 export chain — fix round 1. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` first.

Files: `scripts/S10_brane_mode_spectrum_sympy_audit.py` · `scripts/S10_exports.py` (generated) ·
`scripts/out/S10_brane_mode_spectrum_sympy_audit.out` (regenerated) ·
`directives/S10_export_chain_directive.md` (update it). ⛔ Nothing else. ⛔ Not S9. ⛔ Not any `.wl`.

⭐ **Two independent legs reviewed your build. One wrote its own derivation from the shared physics spec
first and reproduced 26 of 26 committed tags — the physics is right.** Every item below is the export
chain.

---

## F1 · The generated module must be reproducible.

A `frozenset` fixes the order of three exported records; sixteen hash seeds produce **all six**
permutations of them. ⛔ This breaks the only chain-integrity control the design has — a diff between
committed generated modules cannot distinguish tampering from a different hash seed.

**Decision: the generated module is a deterministic function of the run.** ⛔ No iteration order anywhere on
the export path may depend on set or hash ordering. ⭐ Demonstrate it: same input, several hash seeds,
byte-identical output.

## F2 · `omegaSquared` must be one object.

The module binds `omegaSquared` with an assumption; two construction sites build a **bare** symbol of the
same name, and four exported records carry it. Measured: substituting the reported solve variable into the
determinant leaves it **unchanged and silent**; substituting the module's symbol gives zero. Both print
identically.

**Decision: every construction uses the module's object.** ⭐ Then check the whole exported namespace the
same way — **any two same-named symbols with different assumptions anywhere in the merged ledger is the
same defect**, and this is the defect class the import binding exists to remove. Report what you find.

## F3 · The module must not be written before it is validated.

A run that failed its own guard left a complete generated module on disk. A consumer importing it cannot
tell it came from a failed run.

**Decision: validate, then write.** ⛔ A run that does not pass its guards leaves no generated module
behind. ⭐ Demonstrate by making a guard fail and showing no module is written.

## F4 · Reconstruction must not silently alter what it stores.

`srepr`→`eval` collapses an unevaluated numeric relational (`1 > 0` becomes `True`), so a construction with
an allowed stratum fails its own round trip and the run exits 1. ⚠ The shared physics spec requires that
admissibility test to be emitted **with its operands**; the stored record keeps only the verdict.

⚠ This path is untested in the committed artifact because `MAIN` currently produces no allowed stratum —
⛔ but a later step will, and the spec says a physics finding must exit 0.

**Decision: the stored record must reconstruct to what was stored.** ⭐ Demonstrate on a construction that
actually produces an allowed stratum. ⛔ If you conclude the round-trip format cannot carry it, **name that
and stop** — that is a finding about the format, not something to work around by dropping the operands.

## F5 · Record what the deleted registry comparison was doing.

⭐ **This is a record item, ⛔ not a code change. Do not build a replacement.**

The registry block you deleted was, per the shared physics spec, the one place the field-dimension premise
could genuinely fail. Measured: change that premise **on the import alone** and the overwrite guard fires;
change it **consistently on both sides** and the run exits 0 with every guard green, writing coefficient
dimensions wrong by two powers of length into the ledger the next step imports.

**Decision: state that plainly in your directive** — what was removed, what it policed, that nothing
replaces it, and that the premise is therefore unfalsifiable within this engine.

⛔ Do not add a check. The mechanism it used no longer exists, and a self-comparison would be the
zero-by-construction defect this project has now built and deleted three times.

## F6 · Record the redundancy.

One of the three overwritten records is computed as the difference of the other two, and the premise above
cancels in that difference — measured: its residual read 0 while both its inputs read 1. **State in the
directive that it adds no falsification power over the two it is built from.**

## F7 · Record the movement inventory.

Your directive requires movement outside the export boundary to be reported and **names none**. Two control
readouts moved under the imported assumptions. **Name them, with their before and after values**, in the
directive.

---

## Constraints

- ⭐ Every S9 record the build did not overwrite stays byte-identical; the three overwrites stay.
- ⛔ The derivation, action, ansatz and every computed physics value stay untouched.
- ⭐ Report the complete `.out` diff and the complete generated-module diff against the current working tree.
- ⭐ Report how many lines this round removed and how many it added.

## Acceptance

⛔⛔ **A test that passes on a weaker fix is not a test.** For each item, construct the weakest change that
should be rejected and show whether it is.

⛔ Do not add any in-run check on the export writer beyond the overwrite residual that already exists.
Three such checks have been built and deleted on S9; each compared two operands descending from one source.

⛔ Do not state anywhere what a count, tally or partition size comes out as. Emit it and let it be read.

## Deliverables

The four files, the literal stdout, both diffs, and every ablation script with its stdout at named absolute
paths **outside the repository**.
