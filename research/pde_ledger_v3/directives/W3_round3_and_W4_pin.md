# W3 round 3 + W4 — replace the typed duplicate with the engine operand

⚠ **Warning: `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md` holds S9 result values. Don't read it.**

⭐ Round 2's engine repair is **sound and independently confirmed** by two legs and by the orchestrator: all
four fenced engines regenerate their committed outputs. ⛔ Do not disturb it.

⭐ **This round closes W3 and builds the pin it deferred.** They are one change because ⛔ W3's last defect
**is** the missing pin.

---

## ⛔⛔ D1 · The typed duplicate moved instead of dying — third time in this material

`reduction/test_dimension_laws.py:52` types the law evaluated at three points, which determines a linear
law uniquely:

```python
((2, (-2, 0, 1), (0, -2, 1)), (3, (-3, 0, 1), (-1, -2, 1)), (4, (-4, 0, 1), (-2, -2, 1))),
```

⭐ **Measured:** it is the **only** thing under `reduction/` that detects a wrong `D` coefficient — and
changing it together with the registry leaves **everything** green: every check exit 0, `125 passed`, and
all four engines byte-identical.

⚠ **This is the same defect in its third location.** Round 1 put it in `dimension_law_check.py`; round 2
removed it there; it now lives in the test table. ⇒ ⛔ **Do not relocate it a fourth time.**

## ⭐⭐ D2 · The operand that is not a copy of the answer already exists

⭐ The engines emit their **derived dimension symbolically in `D`**, in committed outputs, ⛔ and three of
them never read the registry at all:

```
mathematica/S9_light_requires_shear_mathematica_audit.wl    registry references: 0
mathematica/S10_brane_mode_spectrum_mathematica_audit.wl    registry references: 0
scripts/S9_light_requires_shear_sympy_audit.py              registry references: 0
scripts/S10_brane_mode_spectrum_sympy_audit.py              registry references: 5   (numeric view only)
scripts/S11_stray_longitudinal_sympy_audit.py               registry references: 2   (numeric view only)
```

⇒ ⭐⭐ **Build the comparison of each engine's emitted symbolic dimension against the registry's declared
law**, and ⭐ **emit both operands and the residual, then guard.**
⛔⛔ **This directive does not tell you what that residual is. Do not assume it. Do not tune anything toward
a value. ⭐ Whatever it computes is the deliverable** — a disagreement is a finding, not a failure.

### ⚠ Scope constraints — measured, ⭐ verify them

- ⛔ **`Q.brane.B_comp` appears ZERO times in all four S9/S10 committed outputs and 736 times in S11-py's.**
  ⇒ a pin scoped to S9+S10 cannot cover one of the three required quantities. ⭐ Span **S9, S10 and S11**.
- ⚠ **S11-py already emits both operands** and subtracts a `D = 3` **specialisation** rather than the law,
  leaving a standing nonzero residual. ⭐ That residual is a specialisation artefact, ⛔ not a physics
  disagreement — ⭐ say so wherever it is recorded.
- ⭐ This pin is computable from **committed outputs**; ⛔ it does not require running an engine.

⇒ ⭐ **Delete the typed triple and derive the expectation from the engine operand.** ⛔ If some test genuinely
cannot be expressed that way, ⭐ say which and why, ⛔ rather than keeping a second copy quietly.

## ⛔ D3 · `README.md:116-118` states the opposite of what is measured

> *"A wrong coefficient is not detected by anything in `reduction/`."*

⭐ Measured: a **single-sided** wrong coefficient fails the reduction-local test population loudly.
⇒ ⛔ The record is false. ⭐ The accurate statement has **three** parts — what the production check does not
do, what the tests do catch and by what means, and what remains deferred. ⭐ Rewrite it, in
`README.md` and in `REBUILD_HANDOFF.md`, ⛔ and do not replace one absolute claim with another.

## ⛔ D4 · Three checks that cannot fail where they claim to

1. `w3_acceptance_ablations.py:87-93` — its `expected` is **three identical exit pairs**, so intact,
   laws-absent and wrong-coefficient are indistinguishable to its verdict. ⭐ It prints discriminating
   sub-verdicts and then compares none of them.
2. `dimension_law_able_to_fail.py:281-286` (`demonstrate_escape`) — the expected operand, observed operand,
   residual and status are **four typed literals**, formatted exactly like computed operand lines
   elsewhere. ⚠ Also `:225` types `reference = 3` and prints it against a registry that may say otherwise.
3. ⚠ Replacing a computed status print with a literal `PASS` leaves every check exit 0 and the suite green.
   ⇒ ⭐ the **print path** is unguarded while the exit code is permanently nonzero.

⭐ Each must be able to fail for the case it claims to police. ⛔ **Build the weaker implementation of each
and show the test fails on it**, with literal stdout.

## ⛔ D5 · The regression fence was drawn by step number, not by consumer

⭐ Of the four fenced engines, **one** reads the registry; **S11-py reads it and was outside the fence.**
⇒ ⭐ **Fence by `imports registry_read`**, computed — ⛔ not by a hand-maintained list of engines.
⭐ Whatever you build must make that set answerable by a command.

## ⚠ D6 · `generate_rows.py` cannot emit a schema-valid row

⭐ It emits `LTM-exponent-vector-v1` against a schema requiring the v2 convention, and its `build_row` has
**no law support** — so any brane quantity it generates would be a bare `D = 3` vector, recreating the
defect W3 removed. ⭐ It validates its own output, so it fails loudly. ⭐ Repair or explicitly retire it.

---

## ⭐ The regression bar

⭐ Every `D = 3` result unchanged: gate, loader, `acceptance_check.py`, `able_to_fail.py` and its cases,
`engine_output_checks.py` on both committed configs, the full test population against a **pristine `HEAD`**,
and ⭐ **every engine that imports `registry_read`, re-run and diffed against its committed output.**
⚠ `WL_S10_RUNTIME_SECONDS` moves every run. ⚠ `wolframscript` writes configuration warnings to **stderr** —
⭐ separate the streams or you will read them as a diff.
⛔ One Mathematica kernel at a time — the licence has **two** seats; ⛔ wrap each run in `timeout 900`.

⛔ **If a `D = 3` result moves, STOP and report it.** ⭐ Which of "introduced" or "exposed" it is, is the
finding. ⛔ Do not adjust anything to restore the old output.

## Scope

⛔ Do not modify any engine, any step record, any committed `.out`, or any `checks_S*.yaml`. ⛔ Do not commit.
⚠ **If closing the pin requires an engine change, STOP and report it** — an engine change is
physics-bearing and is not yours to make in a fix round. ⭐ Say exactly what would have to change.

## Rules

- ⭐ Print computed objects — both operands and the residual — then guard. ⛔ No `assert` before a value it
  guards. ⛔ No status typed rather than computed. ⛔ **No expected value typed beside the value it checks.**
- ⛔ **Do not tune anything to make a run green.** ⭐ A disagreement is a finding and is the deliverable.
- ⛔ Never `git checkout`, `git stash`, or `git restore` — ⚠ the tree holds uncommitted work a revert
  destroys. ⭐ Absolute paths only; ⛔ never rely on the shell's working directory.

## Deliverables

1. The change, and whatever it forces.
2. ⭐ The pin's literal stdout — both operands and the residual, per engine and quantity.
3. ⭐ The weaker-implementation runs for each object in **D4**, scripts at named absolute paths.
4. ⭐ The regression evidence, including every registry-importing engine.
5. ⛔ A ≤40-line report. ⛔ No conclusions about whether the physics is right.
