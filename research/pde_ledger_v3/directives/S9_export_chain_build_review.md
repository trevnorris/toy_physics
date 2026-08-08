# Independent physics review — the S9 export-chain build

## Artifact

Written by Codex, uncommitted, in the working tree at `/var/projects/toy_physics`:

```
research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py     (modified)
research/pde_ledger_v3/scripts/S9_exports.py                              (new, generated)
research/pde_ledger_v3/scripts/extract_knob_inventory.py                  (new)
research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl  (modified)
research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out       (regenerated)
research/pde_ledger_v3/mathematica/out/S9_light_requires_shear_mathematica_audit.out (regenerated)
```

The pre-edit baseline of every one of these is `git show fb29bba2:<path>`.

The directive it was built from is
`research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md`, and the governing plan is
`research/pde_ledger_v3/S9_REWRITE_PLAN.md`. **Read the directive last, not first** — an artifact can
satisfy its directive and still be wrong, and that case is exactly what this leg exists to catch.

## What the build claims

S9's SymPy engine now emits a flat 139-entry `LEDGER` for the next step to import; twelve objects are
emitted by both engines under one standard name; the Wolfram engine changed in emitted name strings only;
a mapped diff against the baseline exits 0 on both sides; and a cross-engine residual block reports
`CROSS_ENGINE_RESIDUALS: PASS (12/12)`.

## ⛔⛔ The first question, and it may be fatal

**The SymPy engine imports nothing from the Wolfram engine — by design, that silo is the only blindness
control this architecture has. So where did the Wolfram operands in its cross-engine residual block come
from?**

Three possibilities and they are not equally bad:

1. it parses the committed Wolfram `.out` at runtime — then the "independent" check is a file dependency,
   and it is **circular** if that `.out` was regenerated from the same build;
2. the Wolfram values are **typed into the Python source as literals** — then the residual is two copies of
   an assumption agreeing, and it is fabrication;
3. something else.

**Determine which, quote the lines, and say plainly whether `RESIDUAL=0` carries any information.** This
project has measured a script emitting a physics conclusion as a typed sentence with nothing computing it,
missed by eight review legs, so treat a clean 12/12 as the thing most in need of proof.

## Required method — this is a SCRIPT review

**Derive independently, and write your own script to do it.** A prose re-derivation is worth nothing here:
> *"trusting how they 'rederive' is not trustworthy. Unless it's in CAS and we can see the output from the
> inputs, it's not to be trusted."*

Write your derivation script **before** opening the artifact and save both the script and its literal
stdout to named absolute paths. Report those paths. Without them your derivation claims are discarded.

**Ablation is mandatory, not optional, and it is the only thing that has ever caught the worst defect
here.** Copy the artifact to a scratch directory and ablate the **copy**; never modify the working tree.

⛔⛔ **A FORM ablation is required.** Change the *structure* of a load-bearing object — replace the
curl-only stiffness with a different form, collapse two independent symbols into one, flip a sign **and**
an off-diagonal. A **coefficient** rescale only tests arithmetic; only a **form** change tests physics,
because scaling never leaves the family. Report the literal diff in the output.

Ask of every claim: **which line computed this?** Give the line number or report it as uncomputed.
⚠ An `assert` that precedes the value it guards hides this — a perturbation strong enough to flip the check
kills the process, so you see only PASS-or-crash. Report any such `assert`.

## What to attack

**① Are the twelve standard-name emissions read out of the live computation, or assembled beside it?**
This is the project's signature defect: a tag that *declares* what the run used, built from a literal next
to the computation instead of read out of it. Nine such defects were found in one recent engine and **eight
of the nine were invisible by reading and visible only by mutation.** For each of the twelve, mutate the
object it should come from and confirm the tag moves.

**② Is the mapped diff able to fail?** The build reports `PY_MAPPED_DIFF_EXIT=0` and `WL_MAPPED_DIFF_EXIT=0`
using a `sed` rename map and an `rg -v` exclusion of the additive tags. Construct and run the cases: a
changed payload on a preserved tag; a **dropped** baseline tag; a changed payload on one of the *excluded*
additive tags. Does the exclusion pattern hide anything it should not?

**③ `COMPUTED_DIMENSION_COUNT: 21` against the 17 the directive measured.** Four extra. The directive
forbids manufacturing a dimension and requires the field to be **absent** where the engine's solve produced
none. Identify every one of the 21, and say for each whether a dimension solve actually produced it.

**④ Is `S9_exports.py` generated from the emitted collection itself?** The directive requires it, because a
second name→object mapping permits an entry bound to the *wrong* live object while the `.out` stays
byte-identical and every check passes. Verify the generation path. Then try to construct that defect: can
you make an export entry disagree with the emission of the same name?

**⑤ The 139 entries.** `MAIN` emitted 132 objects. Reconcile 139 with 132, name the extra entries, and
check that no `X1`–`X8` control leaked into the export — controls are ablation evidence, not exports.

**⑥ Class assignments.** The inventory reports the displacement field and both amplitude triples as
`PREMISE`, and every `MAIN` entry as `DERIVED`. Are those defensible? A wrong class propagates into the
knob frontier — the record of what this project still has not determined. Consider in particular whether an
object built from a posited stiffness form is honestly `DERIVED`.

**⑦ The Wolfram engine: is it really a rename?** D12 permits changing emitted name strings and nothing
else. Diff it against `git show fb29bba2:` and confirm no change to the derivation, the action, the ansatz,
the premises, or any computed value.

**⑧ Eleven counterpart pairs beyond the twelve** were reported as found with zero residuals but not
renamed. Verify that they are genuine counterparts and that the zero residuals are real — and see ⛔ above
about where the Wolfram operand came from.

## Physics filter

> Report a finding only if it catches a way the physics could be wrong. Do not report "the script would be
> wrong on a different input."

**This project's failure mode is adding machinery, not omitting it.** A finding that a piece of this is
unnecessary is worth as much as one that it is insufficient. Do not propose a new check unless its absence
lets a specific wrong answer through — and then name that answer.

A leg returning "nothing survives the filter" is weak evidence on its own: state what you checked, what you
could not, and what would have had to be true for you to find something.

## Operational constraints

- ⛔ Wrap **every** Mathematica kernel run in `timeout 600`. A 600 s hit is a FAILED ablation — report it
  and move on. ⛔ Never raise the timeout, and ⛔ never run more than one kernel at a time: the licence has
  **two seats** and another leg may be running.
- ⛔ Copy any artifact to a scratch directory and ablate the **copy**. ⛔ Never modify the working tree, and
  ⛔ never `git checkout` anything — the baseline here is **uncommitted** and a checkout would destroy it.
- ⛔ Do not commit. ⛔ Do not restore `research/pde_ledger_v3/reduction/`, which was deleted at `fb29bba2`.
- ⭐ Save every ablation script **and** its literal stdout to named absolute paths, and report those paths.
- ⭐ Use absolute paths. A `cd` into the wrong directory has corrupted live files in this project.

## Output

Rank most-severe first. For each finding: **the claim · the evidence (a literal command and its output, or
a quotation with `file:line`) · what must change.**
