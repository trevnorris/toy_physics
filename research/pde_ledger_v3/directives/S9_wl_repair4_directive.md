# S9 `.wl` REPAIR 4 — final cleanup: retire the superseded tags, and stop losing the controls

⛔ **Repair in place. Do not rewrite.** Three fixes, nothing else.

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

Verify: `cd /var/projects/toy_physics && timeout 600 math -script research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`
⚠ Two-seat licence; one kernel at a time, retry on a licence error. Iterate to exit 0.

⛔⛔ **Do not read any `.py` engine.** ⛔ Do not commit.

**The standing clauses hold:** emit computed objects, never prose and never hand-typed algebra; emit
operands and residual **before** any guard; no verdict tag, no conclusion in a tag name; physical symbols
combined by hand **only** in the actions and the ansatz; ⛔ emission never conditional on a payload's
**value**.

---

## H1 — ⛔⛔ RETIRE the `q`-substitution homogeneity tags; the scaling block supersedes them

The `q`-based homogeneity family — `HOMOGENEITY_DEFECT*`, `SPEED_VARIATION*`, `NUMERATOR_Q_DEGREE*`,
`DENOMINATOR_Q_DEGREE*`, and the root-expressed-in-`q` tag — is computed **after** a textual replacement
of the wavevector norm by `q`. ⛔ **When that replacement does not consume all wavevector dependence, the
derivative with respect to `q` treats the leftover as a constant and the result is wrong.** It fails in
**both** directions, and both were demonstrated by ablation:

- a **gapped** root reports a defect equal to minus the whole root, i.e. spuriously "dispersive";
- a **genuinely dispersive** root (an action carrying both a curl and a biharmonic stiffness) reports a
  defect of `0` and a `q`-degree of `1` — i.e. **spuriously non-dispersive**, with exit 0 and every gate
  green.

⭐ **The scaling block added in the previous repair is the correct test and is already able-to-fail.**
It needs no substitution and cannot silently no-op.

### The fix
⭐ **Delete the entire `q`-substitution homogeneity family**, including the wavevector-occurrence guard
that was added only to police it. ⛔ Do not keep a known-unreliable value in the output because something
adjacent flags it — a flagged wrong value gets quoted without its flag eventually.

⭐ **Keep** the candidate propagation speed. But compute it **without** the `q` replacement: obtain it from
the root and the computed wavevector norm directly, so no pattern match is involved. ⛔ If a control's
speed is not a wavevector-independent constant, ⛔ do not force it to look like one — emit what the
computation returns.

⭐ **Keep** the scaling block exactly as it is: the scaled expression, `lambdaScale²·ω²`, and their
residual, per root per action.

## H2 — a failure in one package must not suppress the others

The guards live inside the per-package emission, and packages emit in order, so **a hard stop while
emitting the base package suppresses every control package that follows.** ⚠ Measured: an induced failure
reduced output from 1487 lines to 150, and another to 12 — the controls were computed, or computable, and
never reached the reader.

⛔ This is *"a guard converts informative values into a crash"* at package granularity, and it is the
pattern this rebuild exists to remove.

⭐ **Fix: emit every package first, then stop.** Collect each package's guard outcome as computed data,
emit **all** packages and **all** their diagnostics, and only after the final emission apply the hard stop
with the appropriate non-zero exit code. ⭐ **Emit first, stop second** — a reviewer ablating one package
must still see the other eight.

⛔ Do not change which conditions are fatal, and ⛔ do not change any exit code's meaning. Only the
**ordering** changes.

## H3 — derive the wavevector-norm dimension instead of typing it

The dimension map assigns the wavevector norm's dimension by a **hand-typed** entry. ⚠ It is the one
dimension identity linking that scalar back to the wavevector, and it sits **inside the block it
polices** — a check whose input is typed by the thing it checks.

⭐ **Fix:** obtain it by running the script's own dimension walk on the **computed** wavevector norm
expression, and emit the result. ⛔ Do not type it, and ⛔ do not assert it equals anything.

---

## Report back — under 20 lines

1. Exit status, runtime, tag count, all-unique and no-untagged confirmation.
2. **H1:** the tag names you deleted, and the literal candidate-speed tags for the **main**, the
   **flexural** and the **bare-field** actions. ⛔ Do not interpret them.
3. **H2:** the tag count emitted when the **base** package's guard trips. ⛔ Do not interpret it.
4. **H3:** the literal derived dimension of the wavevector norm. ⛔ Do not interpret it.
5. What you could not do, and ⭐ anything that surprised you.
