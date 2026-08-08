# W2 — compare EMITTED dimension vectors against the REGISTRY's declared vectors

⚠ **Warning: `steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Do not read it.** That is a warning,
not a control.

## The object to build

`reduction/` currently compares dimension vectors **engine to engine**, and
`dimensional_homogeneity_gate.py` walks **registry relations**. ⛔ Nothing compares an **engine's emitted
dimension vector** for a quantity against that quantity's **declared exponents** in
`reduction/quantities.yaml`.

⭐ **Build that comparison.** For every registry quantity an engine emits a dimension vector for, emit:
the **declared** exponent vector, the **emitted** exponent vector, and their **difference** — then guard.

⛔ **Name the object, ⛔ not the recipe.** Where it lives, how rows are selected, and what the config looks
like are yours to decide from the existing code. ⛔ Do not restructure the harness to accommodate it.

## ⛔⛔ Two conditions — skipping either makes the check worthless

1. ⛔ **Compare the axis CONVENTION before comparing any exponent.** `quantities.yaml` declares
   `LTM-exponent-vector-v1`; the engines emit bare triples. ⚠ A wrong axis order makes every row agree or
   disagree trivially, and this repository has already measured an `M,L,T` vs `(L,T,M)` disagreement
   between an engine pair (`STATUS.md`, stage037). ⇒ carry a convention tag on **both** sides and compare
   it first. ⛔ If a side does not declare one, that is `UNDETERMINED` — ⛔ never an assumed match.
2. ⛔ **Handle the brane-dimension `D` explicitly.** The engines emit dimensions **symbolic in `D`**; the
   registry declares a numeric specialisation. ⛔ A `D`-symbolic vector must **never** silently compare
   equal to a numeric one. ⇒ either specialise explicitly and **emit the specialisation you applied** as
   part of the record, or report `UNDETERMINED`. ⭐ Whichever you choose, the record must say which.

## Rules that govern the code

- ⭐ **Print computed objects; ⛔ never state conclusions.** Emit both operands and the residual, then
  guard. ⛔ A residual asserted zero always prints `0` and carries no information.
- ⛔ **No `assert` before the value it guards** — a perturbation strong enough to trip it kills the process
  and blinds the ablating reviewer.
- ⭐ Distinguish statuses: agreement · disagreement · convention mismatch · undetermined · not emitted.
  ⛔ Do not fold any of these into another.
- ⛔ **Do not fix anything the check finds.** Report it.

## Able-to-fail — ⛔ both directions, or it is a demo

⭐ Demonstrate, with the literal stdout of each run:
1. perturb **one declared exponent** in a copy of `quantities.yaml` ⇒ the row flips;
2. perturb **one emitted exponent** in a copy of an engine output ⇒ the row flips **independently**;
3. perturb the **axis convention** on one side ⇒ the convention status fires **before** any exponent
   comparison.

⚠ ⛔ A check that only fires in one direction is not calibrated.

## Then run it, and report literally

⭐ Run against the **existing committed** S9 and S10 outputs under `scripts/out/` and `mathematica/out/`.
⛔ **Change nothing** to make it green. ⭐ **Report exactly what it prints**, including every
`UNDETERMINED` and every not-emitted quantity, with counts.

⚠ ⛔ **Do not tune the check until the output looks reasonable.** If it reports a disagreement, that is the
deliverable — a disagreement is a finding, and this is the first instrument in the project able to make
one against an oracle outside a step's own specification.

## Deliverables

1. The check, in `reduction/`.
2. The three able-to-fail runs, each with its script path and literal stdout.
3. The literal output of the run against committed S9 and S10 outputs.
4. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right — that is the step record's
   job, not yours.

⛔ **Do not commit.** ⛔ Do not modify any engine, any step record, or any committed `.out`.
