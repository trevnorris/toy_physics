# REVIEW — a repair that claims to close a demonstrated hole

**Read-only.** ⛔ Do not edit or create any file under `/var/projects/toy_physics`. Use `/tmp` for
scratch and for any ablation copy. Report findings as text.

## The hole this repair claims to close

`research/pde_ledger_v3/reduction/relations.yaml` records relation `R5` as
`c_L − √(B_comp/ρ_br) = 0`. A reviewer demonstrated that rewriting it to `c_L − √(μ_R/ρ_br)` — silently
declaring the longitudinal wave speed equal to the transverse one, **the precise claim this derivation
step exists to settle** — leaves **all five gates green**: acceptance `MATCH`, dimensional gate `PASS`,
able-to-fail `PASS`, 11 tests, script `VERDICT: PASS`.

Three guards interlock and all miss it: the two moduli share a dimension so homogeneity is blind; the
designated output stays a fresh variable so the constraint count is unchanged; and the independent
Mathematica engine cannot help because it must not read the registry.

**The repair** — see `_scratch/S11_sympy_repair_directive.md` — was to add an assertion to
`research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` that substitutes the root the
script **derives** into the registry residual and asserts it vanishes. Inspect the change with
`git diff` (nothing here is committed).

---

## What to check

**P1 · ⭐⭐ IS THE NEW CONTROL REAL, OR DOES THE REGISTRY CHECK ITSELF?**

This is the whole review. The control is only worth anything if it compares **two independently obtained
things**: the expression the script *derived* from the Lagrangian, and the expression the registry
*records*.

⛔⛔ **The failure mode to hunt:** if the assertion builds its "derived" side by reading quantities out
of the registry and reassembling them, or by reusing a symbol the registry relation itself supplied, then
it is **the registry agreeing with itself** and the hole is still open — while now looking closed, which
is worse than before. Trace the provenance of **every symbol and expression** on both sides of that
assertion back to where it came from.

**P2 · ⭐ Ablate it.** In a `/tmp` copy, rewrite `R5` to `c_L − √(μ_R/ρ_br)` — exactly the mutation that
previously passed everything — and run the script. **The new assertion must FAIL.** Report the actual
output. ⚠ If it passes, that is a blocking finding regardless of how the code reads.
Then try at least one *different* corruption of `R5` (e.g. an exponent or a factor) and report whether
that is caught too, or whether the control only catches the one case it was built against.

**P3 · Was `R4` covered?** The directive said to extend the same control to `R4` (`c_γ`) if cheap. If it
was covered, ablate `R4` the same way and confirm it fails. If it was not, say so — that is a scope
report, not a defect.

**P4 · Are the pre-existing outputs unchanged?** The repair was to add a control, ⛔ not to alter physics.
Every `S11_*` line that existed before must print the same thing. Any changed physics value is a
blocking finding — it would mean something was silently absorbed.

**P5 · Did the verdict stay earned?** The new assertion must appear in the printed enumerated list with
its own PASS/FAIL, and the verdict must still be the conjunction of that list.

**P6 · The acceptance fixture comment.** `reduction/acceptance_check.py`'s comment was updated to record
that the payload has three independent derivations. ⛔ Verify **no number changed** — comment only. Is
the claim it now makes accurate to what the comment can support?

**P7 · Tautology recurrence.** Two independent scripts in this project have already been caught with the
same shape: **a value verified using the predicate or definition that produced it** (`c ≔ √(x)` then
asserting `c² − x = 0`). Check the new code for that shape specifically.

---

## ⭐⭐ THE FILTER

Report a finding **only if it catches a way the PHYSICS or the CONTROL'S EFFECTIVENESS could be wrong.**

⛔ Not: style, naming, formatting, comments, performance, "wrong on a different input", error handling.

⚠ A short report with one real finding beats a long one. If the repair is sound, say so plainly — "the
control is real, here is the ablation output" is a complete result.

## Output format

```
VERDICT: CONTROL IS REAL | CONTROL IS NOT REAL | BLOCKING FINDINGS
ABLATION      (P2: the exact mutation, the exact output, PASS or FAIL)
SECOND ABLATION (P2: a different corruption — caught or missed?)
SELF-CHECK RISK (P1: provenance of each side of the assertion — independent, or circular?)
R4            (P3 — covered and ablated, or not covered)
UNCHANGED     (P4 — confirmed, or list what moved)
OTHER         (P5/P6/P7 — list, or "none")
```
