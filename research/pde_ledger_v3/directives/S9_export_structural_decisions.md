# S9 export — carry the structural parameter. Decision list.

**You author the directive prose and apply the change.** Read `CLAUDE.md` and
`research/pde_ledger_v3/S9_REWRITE_PLAN.md` first.

Files: `scripts/S9_light_requires_shear_sympy_audit.py` · `scripts/S9_exports.py` (generated) ·
`scripts/out/S9_light_requires_shear_sympy_audit.out` (regenerated) · your directive at
`directives/S9_export_structural_directive.md`. ⛔ Nothing else. ⛔ Not the `.wl`.

---

## What must become true

⭐ **A consumer of `S9_exports.LEDGER` must be able to bind the brane spatial dimension symbol from the
ledger, rather than declaring its own.**

⚠ **Why, measured:** S10 declares `Symbol("D")` while S9's exported dimensions contain
`Symbol("D", integer=True, positive=True)`. Two SymPy symbols with the same printed name and different
assumptions are **distinct objects** — their difference does not simplify to zero and the residual
**prints as something that looks like zero and is not.** Today the symbol reaches the ledger only buried
inside other entries' expressions, so a consumer has to fish it out of an expression tree.

⭐ The export's declared-input path currently admits one class of declaration. **Widen it so the structural
parameter is carried**, and let the mechanism follow from how the declarations are already classified —
⛔ do not special-case a name.

## What you will run into, and what to do about it

⚠ **Not every declaration in the widened class is a physics object.** One of them is a Python container of
the engine's action packages, including the control packages that the export boundary deliberately
excludes.

⭐ **Report it; do not force it through and do not quietly drop it.** If its class annotation is wrong,
say so and say what the right one is — that is a finding about the annotation, not about the export. If
you change an annotation, that is the change, and the export follows from it.

⛔ Do not add a name-based exclusion list to the export writer.

## Constraints

- ⭐ Every existing entry keeps its key, value, dim, class and step. **Emit the operands and the residual**
  that establish that, in the shape the export already uses.
- ⭐ The exact-`srepr` round trip must cover whatever is added.
- ⛔ The derivation, the action, the ansatz, the assumptions and every computed value stay untouched.
- ⭐ Report the complete `.out` diff against the committed baseline.

## The governing rule this change is subject to

⭐ **The script is input-driven: the only route to an output is a derivation.** ⛔ Nothing about the added
entry may be typed where it could be read — not its name, not its class, not its assumptions.

## Acceptance

- ⛔⛔ **A test that passes on a weaker fix is not a test.** Show what happens if the widening is applied to
  some declarations of the class and not others.
- **Demonstrate that a consumer binding the symbol from the ledger gets the same object the engine used** —
  operands and residual, ⛔ not a claim.
- ⛔ Do not state in the directive or a comment what any count or tally comes out as. Emit it.

## Deliverables

The changed files, your directive, the literal stdout, the complete `.out` diff, and every ablation script
with its stdout at named absolute paths outside the repository.
