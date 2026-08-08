# W3 — let the registry declare a dimension as a LAW, not only a number

⚠ **Warning: `steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

## Why

Three independent parties derived the same thing: **nobody's physics is wrong.** On a D-dimensional brane
an inertial density is mass per D-volume and a first-gradient modulus is energy per D-volume, so
`[rho_br] = [-D,0,1]`, `[mu_R] = [2-D,-2,1]`, `[B_comp] = [2-D,-2,1]`. Both engines **solve** for those
laws symbolically. `reduction/quantities.yaml` records only the `D = 3` evaluation, with nothing saying so.

⇒ ⛔ **Changing `Q.brane.D_brane.value` today silently invalidates three declarations and nothing catches
it.** That is the gap to close.

## The object to build

⭐ **A dimension declaration that can be a law in a bound structural quantity**, and a loader that can hand
consumers either the law or its evaluation.

⛔ **Name the object, ⛔ not the recipe.** The schema form, the loader surface, and how consumers obtain an
evaluation are yours to decide from the existing code. ⭐ Whatever you choose must make the three brane
quantities express their true `D`-dependence and must let a consumer ask for the value at a stated `D`.

⚠ **Known constraints in the existing code — verify, do not trust this list:**
`registry_schema.yaml` requires integer exponent components and forbids extra fields on the dimension
object; `registry_read.py` calls `int()` on every component; `dimensional_homogeneity_gate.py` raises
`UndeterminedFinding` unless it sees exactly three bare integers; `engine_output_checks.py` loads a
complete `Registry` before reaching its own machinery. ⭐ Advance schema/convention versions where the
semantics genuinely change.

## ⛔⛔ THE REGRESSION BAR — this is the acceptance criterion

⭐ **Every result that exists today at `D = 3` must be unchanged.** Concretely, and report each as literal
stdout, before and after:

- `python dimensional_homogeneity_gate.py` — ⭐ same per-relation classifications, same population counts,
  same exit code.
- `python registry_read.py`, `python acceptance_check.py`, `python able_to_fail.py` and its cases.
- `python engine_output_checks.py --config checks_S9.yaml` and `--config checks_S10.yaml` over the
  committed outputs — ⭐ **the cross-engine and registry-residual verdicts must not move.**
- The full `reduction/` test suite.

⛔ **If any D=3 result changes, STOP and report it.** ⚠ It is either a defect you introduced or a defect
this change exposed — ⭐ and which one it is, is the finding. ⛔ Do not adjust anything to restore the old
output.

## What must become possible, and what must become catchable

1. ⭐ The three brane quantities declare their `D`-dependence, bound to `Q.brane.D_brane`.
2. ⭐ `Q.brane.c_gamma` and `Q.brane.c_L` stay constant vectors — ⚠ **verify that is right**: `[mu_R] −
   [rho_br] = [2,-2,0]` with `D` cancelling, and likewise via `B_comp`, so both speeds are `D`-independent.
   ⛔ If your check disagrees, that is a finding — report it, don't reconcile it.
3. ⛔ **A wrong law must be catchable.** ⭐ Show it: perturb a declared law and show what fires; perturb a
   binding (point a law at the wrong structural quantity) and show what fires. ⚠ Both directions, with
   literal stdout.
4. ⭐ **An unbound or unresolvable symbol must be `UNDETERMINED`**, ⛔ never presumed equal and ⛔ never
   silently specialised.

## ⚠ A blindness that already exists — record it, ⛔ do not claim to fix it

A common-mode shift applied to **all three** brane constituents leaves **all five relations homogeneous**;
only shifting one makes R4 and R5 fail. ⇒ ⭐ **the relation gate catches differential errors only**, and
this change does not alter that. ⭐ State it plainly in `reduction/README.md`. ⛔ Do not let anyone read a
green gate as evidence about absolute brane dimensions.

## Scope

⛔ **Do not touch `reduction/registry_dimension_witness*`** — it has its own outstanding fix round and will
be rebuilt against whatever you land here.
⛔ **Do not modify any engine, any step record, any committed `.out`.** ⛔ Do not commit.

## Rules

- ⭐ Print computed objects — both operands and the residual — then guard. ⛔ No `assert` before a value it
  guards. ⛔ No status typed rather than computed.
- ⛔ **Do not tune anything to make a run green.** ⭐ A disagreement is a finding and is the deliverable.

## Deliverables

1. The schema/loader change and whatever consumer changes it forces.
2. ⭐ The before/after regression evidence above, as literal stdout.
3. ⭐ The able-to-fail runs from item 3, with scripts at named absolute paths.
4. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right.
