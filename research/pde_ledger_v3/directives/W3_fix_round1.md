# W3 fix round 1 — the bound dimension law is unpoliced

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

⛔ **W3 does not land as built.** Two legs reviewed it; the orchestrator reproduced the decisive ablations
independently. ⭐ The physics in the three declared laws is **correct** and matches two independent
derivations. ⛔ What is missing is that **nothing constrains it.**

## ⭐⭐ The measurement that blocks the landing — reproduce it before you change anything

Three runs against `reduction/`, differing only in `quantities.yaml`:

```
laws intact                 dimension_law_check exit=0 PASS   gate exit=0 HOMOGENEOUS=5
every law: block removed    dimension_law_check exit=0 PASS   gate exit=0 HOMOGENEOUS=5
D-coefficient wrong         dimension_law_check exit=0 PASS   gate exit=0 HOMOGENEOUS=5
  (declarations correct at the registry's own D, wrong at every other D)
```

⇒ ⛔⛔ **Every check W3 ships is insensitive to the only property W3 exists to express.**
⭐ Computed cause: the relation gate constrains **differences** between declarations, and the `D`
coefficient cancels out of every difference in `relations.yaml`. ⇒ the gate **cannot** see it, and no
amount of work inside the gate will make it able to.

## The objects to build

⛔ **Name the object, ⛔ not the recipe.** ⭐ Decide the shape from the existing code.

1. ⭐⭐ **A check that cannot pass when the thing it checks is absent.** `dimension_law_check.py` currently
   returns success with **no laws present at all**; its predicate asks the registry to confirm what the
   registry declares. ⇒ ⛔ a check whose expected value lives inside the artifact it checks.
2. ⭐ **An evaluation that is `UNDETERMINED` or a hard error when a binding cannot resolve to an integer.**
   Today a missing or non-integer binding value leaves the gate passing and **prints a symbol where a
   dimension belongs**. ⚠ W3's own directive required this and it was not delivered.
3. ⭐ **A declared reference evaluation that cannot silently disagree with the quantity it is bound to.**
   `reference_values` is never compared against the bound quantity's own declared value.
4. ⭐ **A test population that binds under the runner that will actually be used**, and that covers the two
   guards currently shipping untested — the line binding the retained triple to the law, and the
   structural-target check. ⚠ Disabling either leaves the suite green today.
5. ⭐ **One command that reports the health of the able-to-fail cases.** Today each case returns its own
   code, there is no aggregate, and an **escape exits 0**.

## ⛔⛔ A test must FAIL on every WEAKER fix

⚠ **Measured three times in this repository:** a half-fix passed the new test, the whole suite, the full
battery **and** produced byte-identical output. ⭐ For each object above, **build the weaker
implementations yourself and show the test fails on each**, with literal stdout. ⛔ A test that
demonstrates the invariant on one example is a demo, not a pin.

⭐ **The two ablations at the top of this file are the acceptance surface for object 1.** ⛔ If either still
passes after your change, the change did not land.

## ⚠ What to RECORD, ⛔ not to fix

⛔ **Do not claim to close these.** ⭐ State each plainly where a reader will meet it.

- ⛔⛔ **The `D` coefficient stays unpinned by anything in `reduction/`.** ⭐ The pin requires comparing the
  registry's **law** against the engines' **symbolic** dimension emission — S9 and S10 already emit it —
  and that belongs to the witness rebuild, ⛔ not here. ⭐ Say so, and say the gate constrains differences
  only.
- ⚠ **`Q.medium.rho0` and `Q.medium.K` carry structural dependence and are declared as bare constants.**
  ⭐ `K`'s exponent additionally depends on `Q.medium.n_eos`, and ⛔ **the shipped grammar cannot express
  it** — `relations.yaml`'s `Pow` exponent is a bare integer the dimension algebra treats as
  dimensionless. ⭐ Record the gap and what extending to the medium sector would require. ⛔ Do not extend
  it in this round.
- ⚠ `README.md` and `W3_REGISTRY_D_LAWS_REPORT.md` describe the blindness as *"a common-mode shift"*.
  ⭐ The invisible set is larger than a shift. ⭐ Restate it at its real width; ⛔ do not narrow it to what
  was measured first.
- ⚠ The report's regression baseline cites a pytest population that already included staged, uncommitted
  work. ⭐ Correct it against a pristine `HEAD`.

## ⚠ Also in scope — the same defect class, one line to prove

`python -m unittest discover` over `reduction/` reports **OK** while collecting a small fraction of the
population `pytest` collects. ⇒ ⛔ a green run that hides most of the suite. ⭐ Either make collection
agree, or make the intended runner explicit and make a wrong runner's output impossible to mistake for
coverage. ⛔ A check that was not run and a check that passed must not print the same way.

## Scope

⛔ Do not touch `reduction/registry_dimension_witness*` — it is being rebuilt separately.
⛔ Do not modify any engine, any step record, any committed `.out`, or any `checks_S*.yaml`. ⛔ Do not commit.

## ⭐ The regression bar

⭐ Every `D = 3` result must be unchanged: the gate, the loader, `acceptance_check.py`, `able_to_fail.py`
and its cases, `engine_output_checks.py` on both committed configs, and the full test population — the
last measured against a pristine `HEAD` tree.
⛔ **If a `D = 3` result moves, STOP and report it.** ⚠ It is either a defect you introduced or one you
exposed — ⭐ which one it is, is the finding. ⛔ Do not adjust anything to restore the old output.

## Rules

- ⭐ Print computed objects — both operands and the residual — then guard. ⛔ No `assert` before a value it
  guards. ⛔ No status typed rather than computed.
- ⛔ **Do not tune anything to make a run green.** ⭐ A disagreement is a finding and is the deliverable.

## Deliverables

1. The change, and whatever it forces.
2. ⭐ The two top-of-file ablations re-run against your build, literal stdout.
3. ⭐ The weaker-fix runs for each object, with scripts at named absolute paths and their literal stdout.
4. ⭐ The regression evidence.
5. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right.
