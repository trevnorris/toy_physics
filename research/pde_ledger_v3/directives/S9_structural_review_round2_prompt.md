# Independent physics review — S9 structural export, round 2

## Artifact

Uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
git diff -- research/pde_ledger_v3/scripts/S9_exports.py
git diff -- research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out
```

plus `research/pde_ledger_v3/directives/S9_export_structural_directive.md`.

Codex, two passes. You are one of two independent legs with this identical prompt. **Earlier rounds exist;
do not look for them, and you are not being asked to confirm anyone.**

## What the change is for

`scripts/S9_exports.py` is a generated flat `LEDGER` the next step imports. Its declared-input path
previously admitted one class of declaration, so the brane spatial dimension symbol reached the ledger only
buried inside other entries' expressions — and a consumer declaring its own would get a **different
object**: same printed name, different assumptions, a difference that does not simplify to zero. The change
widens that path by declaration class so the symbol and the unit-basis markers are carried directly.

**This round is a repair of the previous one, and it is almost entirely deletion.** The previous round
added a battery of in-run checks; two legs measured that none of them could fail, and they have been
removed rather than repaired. One entry was dropped from the ledger and one declaration's construction was
changed to read a count rather than type it.

## What to check

By mutation, not by reading:

1. **Was anything removed that was doing real work?** A deletion round's characteristic failure is removing
   a check that could in fact fail. For each removed emission, determine whether anything it would have
   caught is now uncaught, and say so concretely.
2. **What does the run still establish, and what does it not?** The directive makes a claim about this.
   Test it: construct changes the run should reject and show whether it does; construct changes only an
   external diff would catch and show that the run does not.
3. **Is the ledger's population correct?** One entry was removed. Was it the right one, and is every
   remaining entry one a consumer can use at a component count other than the engine's? An entry whose
   content is frozen to this engine's count, carried under a key that reads as general, is the defect this
   question exists to catch.
4. **Are the classifications honest?** Several declarations carry a class chosen so that they fall outside
   the export rather than because the class describes them. Judge each against the closed vocabulary
   `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`, and say where a tag is being used as a
   mechanism rather than a description.
5. **Is anything still typed where the computation has it?** Names, counts, dimensions, matrix shapes.
6. **Did every previously exported entry survive unchanged** — key, value, dim, class, step — apart from the
   one deliberately removed? Establish this against the **committed baseline**, not against the run.
7. **Did any derivation, action, ansatz, assumption or computed value move?** It is claimed none did.

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY.** Run a control at the unchanged structure first and show your harness
reproduces the current export.

⛔⛔ **A test that passes on a weaker fix is not a test.** For anything the artifact still guards, construct
the weakest change it should reject and show whether it does.

⚠ **Watch for the recurring defect:** a residual whose two operands are computed from one source by an
invertible pair of steps is zero by construction and polices nothing. If any survive, name them.

Report with line numbers: anything typed where it could be read; a check whose expected value lives inside
the artifact it checks; an `assert` preceding the value it guards; an emitted object redundant with
another.

### Write your own script first

⛔ Write your own verification script **before** opening the artifact; save it and its literal stdout to
named absolute paths and report them. **Without them your claims will be discarded.**

## Do not read

Any file under `directives/` beginning `S9_export_dkey_`, `S9_export_structural_decisions`,
`S9_structural_fix_round1_decisions`, `S9_export_structural_review_prompt`, or `S10_`. No Codex or review
transcript is in the repository.

## Operational constraints

- Ablate **copies** under `/tmp`. ⛔ Never modify the working tree — it holds the uncommitted change and
  there is no committed baseline to restore it from.
- Wrap every run in `timeout 600`. ⛔ Never raise it. ⛔ Do not start `wolframscript`.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

Out of scope and unchanged by decision: `wavevector_norm_dimension` names `dim(|k|)` but holds `dim(k·k)`;
the placeholder-naming class has eight members; `q_dimension` is unpinned inside SymPy; the `.out` record
name for the four D-specialisations is a pre-existing hand-typed literal; the supplied field-dimension
premise is an alias of a unit-basis marker and is therefore not independently checkable. Do not re-report
these unless this round altered one.
