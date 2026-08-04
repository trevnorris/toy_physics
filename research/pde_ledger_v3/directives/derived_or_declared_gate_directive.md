# BUILD — the DERIVED-OR-DECLARED gate

## The defect this exists to catch

An engine script can `emit("TAG", "some prose")` where the prose states a physics conclusion that **no
computation in the script produced**. The tag looks like output. It is a constant. Cross-engine agreement on
such a tag is vacuous, because both engines carry the same author's typed sentence.

⚠ This is **not** the "tautological assertion" failure (asserting `c² − x = 0` after defining `c ≔ √x`).
There is no assertion at all. ⛔ Reading the script does not reveal it — eight review legs across two
sub-steps did not.

## The rule to enforce

> Every emitted tag is either **DERIVED** — its printed text changes under some perturbation of the
> script's symbolic inputs — or **CONSTANT**. A CONSTANT tag is legitimate only if it is a declared
> supplied premise. It is never legitimate as a derived result.

## Deliverable

`/var/projects/toy_physics/research/pde_ledger_v3/reduction/derived_or_declared.py`

⭐ **Model its output vocabulary and CLI shape on the existing sibling**
`/var/projects/toy_physics/research/pde_ledger_v3/reduction/able_to_fail.py`, which does the analogous job
for the registry. ⛔ Do not build a new framework, a new config format, or a plugin system. One file.

### Required behaviour

1. Take a path to an engine script (`.py`; Mathematica support is a separate later task — ⛔ do not attempt
   it here) and run it unperturbed, parsing its emitted `TAG -> text` pairs.
2. Re-run it under each of a small set of **symbol-collapse perturbations**: force two independent free
   symbols to be the same object, so any quantity that distinguishes them must change.
   ⛔⛔ **The perturbation MUST be applied from OUTSIDE the engine** — the engine must not read an
   environment variable, a flag, or any cooperating hook. A control whose trigger lives inside the artifact
   it polices is that artifact agreeing with itself. ⭐ Patching the symbolic layer in `sys.modules` before
   executing the script (e.g. via `runpy`) is the intended shape; choose the mechanism you can make robust.
   ⛔ **Do not perturb by editing the source text.** Spelling-based transforms are evadable and brittle.
3. ⭐ Perturbations must be **FORM** changes (collapsing a distinction), ⛔ **not COEFFICIENT rescalings.**
   Scaling never leaves the family, so it cannot test whether a distinction is being computed at all.
4. Classify each tag: **DERIVED** if its text differs under at least one perturbation; **CONSTANT** if it
   never differs.
5. Read the step's declared supplied-premise list from a plain sidecar file next to the engine script
   (you choose the name and a simple line-per-tag format; ⛔ YAML or plain text, never JSON).
   Exit **1** if any CONSTANT tag is absent from that list. Exit **0** otherwise.
6. Print every CONSTANT tag with its text, and a count line, whether or not the gate passes.

### ⚠ Honest limits you must implement, not paper over

- A tag may legitimately not depend on any collapsed pair. ⇒ ⭐ the CONSTANT list is **for adjudication**,
  and your output must say so. ⛔ Do not invent a confidence score, and ⛔ do not silently suppress a tag.
- ⛔ If a perturbation makes the script crash or exit nonzero, that perturbation yields **no information**
  about any tag. Report it as SKIPPED with the reason; ⛔ never count a crash as "the tag changed."
- Report how many perturbations actually ran. ⛔ A run where every perturbation was skipped must not
  print a pass.

## ⭐⭐ ACCEPTANCE — a known-bad ground truth, ⛔ not a synthetic fixture

Run the finished gate against the **real committed engine**
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bB_interface_assembly_sympy_audit.py`.

⭐ These four tags are **known** to be uncomputed prose in that file and the gate **MUST classify all four
as CONSTANT**:

```
S11BB_COEFFICIENT_ADMISSIBILITY
S11BB_TRANSVERSE_COUPLING
S11BB_TRANSVERSE_DISPERSION
S11BB_TRANSVERSE_DISSIPATION
```

⭐ And these are **known** to be genuinely computed (their text interpolates SymPy objects) and the gate
**MUST classify them as DERIVED**:

```
S11BB_LONGITUDINAL_DISPERSION
S11BB_ROOTS
S11BB_STABILITY_CONDITION
```

⛔⛔ **Both halves must hold. A gate that flags everything is as useless as one that flags nothing.**
⛔ Do not tune the perturbation set until only the four appear — other tags may legitimately be constant.
⭐ The requirement is exactly: **those four appear as CONSTANT, and those three do not.**

Paste the literal output of that run into your final report.

## ⛔ Constraints

- ⛔ **Do not modify the engine script.** This build only reads it.
- ⛔ Do not add any check beyond the rule above. This project has previously died from checks-on-checks.
- Iterate until the script exits cleanly and the acceptance holds. ⛔ Never wrap anything in a shell
  `timeout`; no single run may exceed 10 minutes.
- ⛔ Do not commit. Leave the working tree dirty for review.

## Report — under 30 lines

1. The path written.
2. The literal acceptance output.
3. The perturbation set you chose and **why each is a FORM change**.
4. ⭐ Anything you could not make robust, stated plainly rather than worked around.
