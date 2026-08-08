# W3 fix round 2 — the registry change breaks an engine, and the new pin is a duplicate

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

⛔ **W3 still does not land.** Round 1 fixed real defects. ⭐ Two things are now measured that block it, and
⛔ **one of them is a live breakage in the working tree.**

---

## ⛔⛔ B1 · The registry change breaks S10-python — measured, ⛔ not reported

```
$ python3 scripts/S10_brane_mode_spectrum_sympy_audit.py ; echo EXIT=$?
  File ".../S10_brane_mode_spectrum_sympy_audit.py", line 1696, in emit_registry_comparison
    declared = tuple(sp.Integer(component) for component in declared_dimension)
TypeError: 'BoundDimensionLaw' object is not iterable
EXIT=1
lines produced 4224   vs   committed 4227
```

⭐ The engine reads the registry **directly** and iterates a quantity's declared dimension. Under the
uncommitted `quantities.yaml`, that value is now a `BoundDimensionLaw`. ⇒ ⛔ **the engine cannot run, and
its committed output cannot be regenerated.**

### ⛔⛔ Why no gate caught it — ⭐ this is the governing lesson

⚠ Round 1's regression bar covered the gate, the loader, `acceptance_check.py`, `able_to_fail.py`,
`engine_output_checks.py` on both configs, and the test population. **All green. All byte-identical.**
⛔ **`engine_output_checks.py` reads committed `.out` files. Nothing re-ran an engine.**

⇒ ⭐⭐ **The registry's real consumers were outside the fence.** ⛔ A regression bar for anything under
`reduction/` that does not **re-run all four engines** is measuring the readers and not the consumers.

⭐ **This round's regression bar must re-run all four engines** — `S9-py`, `S9-wl`, `S10-py`, `S10-wl` —
and compare against the committed outputs, excluding only records that are run-dependent by construction.
⚠ `S10-wl` emits `WL_S10_RUNTIME_SECONDS`, which moves every run; ⭐ name every exclusion and justify it.
⚠ Mathematica: ⛔ one kernel at a time, the licence has **two** seats; ⛔ wrap each run in `timeout 900`.

## ⛔⛔ B2 · The D-coefficient "pin" is a second copy of the answer, ⛔ not a measurement

`dimension_law_check.py` now hardcodes the expected laws beside the registry's:

```python
BRANE_BOUND_DIMENSION_CONTRACT = {
    "Q.brane.rho_br": BoundDimensionExpectation({"D": "Q.brane.D_brane"}, (-_D_BRANE, 0, 1)), ...
```

⭐ Measured, in an isolated copy:

```
A  baseline                                        exit=0  DIMENSION_LAW_CHECK: PASS
B  BOTH copies changed to the SAME wrong law       exit=0  DIMENSION_LAW_CHECK: PASS
   (rho_br = 3-2D, mu_R/B_comp = 5-2D; truth -D, 2-D)
```

⇒ ⛔ It detects **disagreement between two declarations**, ⛔ never a wrong law — whoever writes the law
wrong writes the expectation wrong too. ⇒ **a control inside the thing it polices.**
⚠ ⛔ **This is worse than the honest gap it replaced**: a reader now sees *"wrong D coefficient → exit 1"*
and concludes the coefficient is policed.

⚠ Round 1's directive said explicitly: ⛔ **do not attempt this pin here.** ⭐ It was attempted anyway.

### What to do with it

⭐ **Remove the duplicated expectation, or replace it with a measurement that is not a second declaration.**
⛔ Do not keep a check whose expected value is typed beside the value it checks. ⭐ Whatever remains must
state plainly, where a reader meets it, that **the `D` coefficient is not policed by anything in
`reduction/`.**

## ⭐ The object the real pin needs — ⛔ RECORD it, ⛔ do not build it here

⭐ S9 and S10 already emit their **derived dimension symbolically in `D`**. The registry now declares a law
symbolically in `D`. ⇒ ⭐⭐ **comparing those two symbolic objects is the only pin available that is not a
second copy of the answer**, and it is an oracle **outside** the registry.
⛔ It belongs to the witness rebuild, ⛔ not to this round. ⭐ Record it as owed, and say why the checks that
live inside `reduction/` cannot supply it.

⚠ Note the engine that breaks is doing a **specialised-to-`D=3`** version of exactly that comparison.
⭐ Whatever you do for **B1** determines whether that comparison can later become symbolic. ⛔ Do not close
that door; ⭐ say which way you left it and why.

## ⭐ What round 1 got right — ⛔ do not regress it

⭐ Verified independently: with laws absent, and with a wrong `D` coefficient in the registry alone,
`dimension_law_check.py` now exits 1 where it previously exited 0. ⭐ Keep that.

## Scope

⛔ Do not touch `reduction/registry_dimension_witness*`. ⛔ Do not modify any step record or any committed
`.out`. ⛔ Do not modify any `checks_S*.yaml`. ⛔ Do not commit.
⚠ **If repairing B1 requires an engine change, STOP and report it** — ⭐ an engine change is physics-bearing
and is not yours to make in a fix round. ⭐ Say exactly what would have to change and why.

## Rules

- ⭐ Print computed objects — both operands and the residual — then guard. ⛔ No `assert` before a value it
  guards. ⛔ No status typed rather than computed. ⛔ No expected value typed beside the value it checks.
- ⛔ **Do not tune anything to make a run green.** ⭐ A disagreement is a finding and is the deliverable.
- ⛔ Never `git checkout`, `git stash`, or `git restore` — ⚠ the working tree holds uncommitted work that a
  revert destroys. ⭐ Use absolute paths; ⛔ never rely on the shell's working directory.

## Deliverables

1. The change, and whatever it forces.
2. ⭐ **All four engines re-run**, with literal exit codes and diffs against the committed outputs.
3. ⭐ The B2 ablation re-run against your build, literal stdout.
4. ⭐ The round-1 ablations re-run, showing they still fail.
5. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right.
