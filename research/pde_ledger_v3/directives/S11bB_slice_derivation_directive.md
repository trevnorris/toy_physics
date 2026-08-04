# S11b-B — MAKE THE B5 SLICE MACHINE-DERIVED IN BOTH ENGINES

⛔⛔ **ONE CHANGE PER SCRIPT. Nothing else.**

## Files

- `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bB_interface_assembly_mathematica_audit.wl`
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bB_interface_assembly_sympy_audit.py`

Run each after editing (`math -script <path>`, `python3 <path>`); each must complete well under **10
minutes** and the `.py` must exit 0.

## The defect

Both engines report B5's roots on an explicit slice (`k = 0`, impermeable `Λ_A⁰ = Λ_V⁰ = 0`, zero
reciprocal coupling `Λ_X⁰ = 0`), and both report the same quadratic and the same separator. ⭐ **But
neither engine actually derives that quadratic from its own assembled dispersion determinant.**

- The SymPy script solves a **hand-written factored polynomial**.
- The Mathematica script assigns slice quantities that are **never used in any emit** — its slice text is
  prose.

⇒ ⚠ **The agreement on B5 is therefore agreement between two hand transcriptions, not between two machine
derivations.** If either transcription were wrong, nothing in either script would catch it. B5 is this
step's headline deliverable, so that gap matters.

## The fix — the same in both engines

⭐ **Derive the slice quadratic FROM the script's own assembled determinant, by substitution and
cancellation, and emit the result of that derivation.**

Concretely:
1. Take the script's **own** full dispersion determinant — the one it already computes and already emits.
2. Apply the slice conditions to it symbolically: `k → 0`, `Λ_A⁰ → 0`, `Λ_V⁰ → 0`, `Λ_X⁰ → 0`, and the
   prescribed `q_out` on that slice.
3. **Cancel any overall kinematic factor** and report what factor was removed. ⚠ The slice determinant is
   known to carry an overall `ω`-power from the rank loss at `ω = 0`; ⛔ do not silently divide it away —
   **state it.**
4. ⭐ **Compare the surviving polynomial against the quadratic the script currently solves.** Emit the
   symbolic difference. ⛔ Do not adjust either side to make them match.
5. Solve **that** polynomial — the derived one — for the roots, and emit them.

⭐ **The comparison in step 4 is the point of this task.** It is a check that can fail: a wrong
transcription, a wrong slice condition, a dropped factor, or a sign error anywhere between the assembled
determinant and the reported quadratic all produce a non-zero difference.

⛔ **If the difference is NOT zero, report it and stop.** ⚠ Do **not** silently correct the quadratic to
match, and do **not** correct the determinant to match. A non-zero difference is a finding about the script
and must be reported as one.

## ⛔⛔ CONSTRAINTS

1. ⛔ **Do not change any other emitted value.** Every other tag must be byte-identical to its current
   output. ⚠ This will be verified by diffing the full tag set before and after.
2. ⛔ **Do not add checks, harnesses, or reporting machinery** beyond the single derivation-and-comparison
   described above. ⚠ A prior repair pass in this project added 183 lines of `x === x` checks.
3. ⛔ **Do not restructure, rename, or reformat anything.**
4. ⛔ **Do not touch the row order or the sign convention of the assembled determinant** in either script.
   ⚠ One engine's row order was recently changed and had to be reverted; the two engines use **different**
   row orders and that is fine — each must derive from **its own**.
5. ⛔ **Report whatever roots the derivation gives.** A growing root is a first-class outcome; ⛔ do not
   suppress, re-branch, or add a stability assumption.

## Output

The two edited scripts, each run to completion, plus a report **under 20 lines**:
- the overall factor cancelled in each engine, stated explicitly;
- the symbolic difference between the derived polynomial and the previously-solved quadratic, per engine;
- the roots from the derived polynomial, per engine;
- confirmation that no other emitted tag changed.

⛔ Then stop and exit.
