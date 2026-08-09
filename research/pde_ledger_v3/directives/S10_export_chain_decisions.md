# S10 onto the export chain — decision list

**You author the directive prose and apply the change.** Read `CLAUDE.md`, then
`research/pde_ledger_v3/S9_REWRITE_PLAN.md`, then `research/pde_ledger_v3/REBUILD_HANDOFF.md`.

⭐⭐ **S9 is the reference implementation.** `scripts/S9_light_requires_shear_sympy_audit.py` and its
generated `scripts/S9_exports.py` (closed at `b23c5334`) are the shape S10 takes. **Read them first.**
⛔ Where S9 has an answer, use S9's answer rather than inventing one.

Files: `scripts/S10_brane_mode_spectrum_sympy_audit.py` · `scripts/S10_exports.py` (**generated**) ·
`scripts/out/S10_brane_mode_spectrum_sympy_audit.out` (regenerated) · your directive at
`directives/S10_export_chain_directive.md`. ⛔ Nothing else. ⛔ Not `S9_exports.py`, ⛔ not the S9 engine,
⛔ not any `.wl`.

---

## Stage 1 · Make it run, and fix the baseline. ⛔ Report before going further.

S10 imports `registry_read` from a deleted directory, so it **cannot run**, and its committed `.out` came
from a version that could ⇒ **there is no live baseline.**

⭐ Delete the registry import and the block that reads the registry and emits a comparison. ⛔ Replace them
with nothing. ⭐ Run, and **report the complete diff against the committed `.out`.**

⚠ Expect one non-registry line to move: a tag that inventories the run's local tag names has the registry
locals in its payload today. ⭐ **That line moving is expected. Anything else is not** — if the diff shows
more, **name it and stop**; that would mean the committed record does not reproduce from its own source.

## Stage 2 · What S10 must produce

1. ⭐ **It imports S9's `LEDGER` and binds the objects**, rather than re-declaring what S9 exports.
   ⚠ Inheriting S9's assumptions is the point (user): *"by committing we learn where it breaks."*
   ⛔ No S10 declaration may shadow an imported S9 object.
2. ⭐ **Every declaration carries `# TAG · English description`**, one tag from
   `KNOB · STRUCTURAL · COORDINATE · CONTROL · PREMISE · DERIVED`. ⛔ An annotation never states what the
   computation will produce.
3. ⭐⭐ **Its exported objects live under fixed names.** ⛔ A key must never be minted from the answer.
   ⚠ S10 currently labels outputs by *which root* the solver returned, so the label depends on how many
   roots there were. S9 instead holds the varying content **inside a stable named collection**. ⇒ do that.
4. ⭐ **It generates `scripts/S10_exports.py`**: S9's entries, plus its own, overwriting what it re-derives,
   exported as one flat dict. Record shape and exact-`srepr` round trip as S9's.
5. ⭐ **Keys carry the component count of the construction that produced the object**, exactly as S9 does —
   ⛔ never decided by inspecting the value. ⚠ Read S9's cut: the dimension solve's **inputs** are built at
   a fixed count and carry it; its **outputs** are symbolic in `D` and do not.
6. ⭐⭐ **Where S10 writes a key that already exists, it emits S9's value, its own value, and the residual,
   then guards.** ⇒ a wrong name landing on an existing entry becomes a **nonzero residual**, not a silent
   corruption. ⭐ These two operands are genuinely independent — one read from a committed file S10 did not
   write, one freshly computed.
7. ⭐ **Its derived values still match stage 1's baseline**, compared by name, and it **shows** that.
8. ⭐ **Nothing is typed that the computation could produce** — no value, count, dimension, rank or
   multiplicity, and no variable portion of any key. ⚠ The authored *name* of an object is exempt; a CAS
   expression does not carry its own name. ⚠ Supplied **premises** are exempt; a premise is an input.

⛔ **Work the rest out against S9 and report what you chose.** If an object does not classify cleanly, or a
correspondence is ambiguous, **name it and stop** — that is a finding about the export boundary, not a
puzzle to solve inside this change.

## What NOT to build

⛔⛔ **No in-run check on the export writer beyond item 6.** Three have been built and deleted on S9; every
one compared two operands descending from one source, and one of them passed while the change it policed
was reverted. **The export writer cannot verify itself.** Chain integrity is a **diff between committed
files** — the plan says so.
⛔ Cross-engine (py↔wl) naming, and ⛔ the comparator. Separate pass, decision not taken. ⛔ Do not touch
the Wolfram engine. ⛔ No YAML, registry, runner, harness or committed test framework. ⛔ Nothing in S11.

## Acceptance

- ⛔⛔ **A test that passes on a weaker fix is not a test.** For each guard, construct the weakest change it
  should reject and show whether it does.
- **Run a control first** and show your harness reproduces the artifact.
- **Demonstrate item 5 by mutation**, and say what your ablation separates and what it does not.
- **Demonstrate item 6 by mutation**: make S10 write a wrong value to an existing key, show the residual.
- **Report a perturbation matrix**: perturb each declared input, including at least one **FORM** change,
  and name every `DERIVED` entry nothing moves. ⛔ Do not manufacture a perturbation or reclassify an entry
  to make it come out clean — an honest *"nothing I perturbed moves this"* is the useful answer.
- ⛔ Do not state anywhere what a count, tally or partition size comes out as. Emit it and let it be read.

## Deliverables

The four files; stage 1's diff reported before stage 2's results; the literal stdout; every ablation script
and its stdout at named absolute paths **outside the repository**; and a statement of what S10 imports,
what it overwrites, and what it now derives that it previously re-declared, each with line numbers.
