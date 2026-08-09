# Independent physics review — S9 export carries the structural parameter

## Artifact

The uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py
git diff -- research/pde_ledger_v3/scripts/S9_exports.py
git diff -- research/pde_ledger_v3/scripts/out/S9_light_requires_shear_sympy_audit.out
```

plus `research/pde_ledger_v3/directives/S9_export_structural_directive.md`.

Written by Codex in one pass. You are one of two independent legs with this identical prompt; you will
not see the other's output. **Earlier review rounds on this file exist. Do not look for them.**

## What the change is for

`scripts/S9_exports.py` is a generated flat `LEDGER` that the next step imports. Its declared-input path
previously admitted one class of declaration, so the brane spatial dimension symbol reached the ledger
only buried inside other entries' expressions. A consumer that declares its own `Symbol("D")` gets an
object that is **not equal** to the one S9 used — same printed name, different assumptions, and a
difference that does not simplify to zero. The change widens the path so the symbol is carried directly.

## What to check

Establish each of these yourself, **by mutation, not by reading**:

1. **Is the widening driven by how declarations are classified, or by their names?** A name-based
   exclusion anywhere on the export path — a skip list, a special case, a string test — is the defect this
   review exists to catch. Mutate a declaration's class and see exactly what follows it.

2. **Several declarations were re-classified in this change.** For each: is the new class *correct for
   what the object is*, or was it chosen to get the object out of the export? Judge against the closed
   class vocabulary `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`. ⚠ If some object is
   not well described by **any** term in that vocabulary, **say so** — a gap in the vocabulary is a
   finding about the vocabulary, and this project has recorded one before.

3. **Are the newly exported entries correctly keyed under the D-key convention?** A key carries the
   component count of the construction that produced the object; objects not built at a fixed component
   count carry none. ⚠ Check each new entry against what actually built it. ⛔ Free-symbol inspection gives
   the wrong answer in both directions and is not the test.

4. **The run emits new residuals. Can each of them fail?** For every one, construct the change it is
   supposed to reject and show whether it does. ⚠ The measured failure mode here is a residual whose two
   operands are computed from one source by an invertible pair of steps — zero by construction, and it
   polices nothing. If you find one, say what it does and does not establish.

5. **Did every previously exported entry survive unchanged** — key, value, dim, class, step? It is claimed
   they did.

6. **Did this change move any derivation, action, ansatz, assumption or computed value?** It is claimed
   nothing did.

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY.** Run a control at the unchanged structure first and show your harness
reproduces the current export.

⛔⛔ **A test that passes on a weaker fix is not a test.** At minimum: the widening applied to some
declarations of the class and not others. Measured three times in this project — a half-fix passed the new
test, the whole suite and the full battery, and produced byte-identical output.

Probe for, and report with line numbers: anything typed where it could be read out of the computation; a
check whose expected value lives inside the artifact it checks; an `assert` preceding the value it guards;
an emitted object redundant with another emitted object.

Ask of every claim in the directive: **which line computed this?** Give the line number or report it
uncomputed.

### Write your own script first

⛔ Write your own verification script **before** opening the artifact and save it and its literal stdout to
named absolute paths. Report those paths. **Without them your claims will be discarded.**

## Do not read

- Any file under `directives/` beginning `S9_export_dkey_` or `S10_`.
- `research/pde_ledger_v3/directives/S9_export_structural_decisions.md` — the orchestrator's decision list.
- Any Codex or review transcript; none is in the repository.

## Operational constraints

- Copy anything you ablate to `/tmp` and ablate the **copy**. ⛔ Never modify the working tree — it holds
  the uncommitted change and there is no committed baseline to restore it from.
- Wrap every run in `timeout 600`. ⛔ Never raise it. ⛔ Do not start `wolframscript`.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

Out of scope and unchanged by decision: `wavevector_norm_dimension` names `dim(|k|)` but holds `dim(k·k)`;
the placeholder-naming class has eight members; `q_dimension` is unpinned inside SymPy; the `.out` record
name for the four D-specialisations is a pre-existing hand-typed literal. Do not re-report these unless
this change altered one.
