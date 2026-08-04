# S9 `.py` REPAIR 3 — a wrong homogeneity value, a sign-blind test, and a crash that discards physics

⛔ **Repair in place. Do not rewrite.** Four fixes, nothing else.

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.py`
**Sidecar:** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S9_light_requires_shear_sympy_audit.premises`

Verify: `cd /var/projects/toy_physics/research/pde_ledger_v3 && timeout 600 python3 scripts/S9_light_requires_shear_sympy_audit.py`

⛔⛔ **Do not read any `.wl` engine or the registry.** ⭐ The two engines must stay independently
constructed so they can disagree. ⛔ Do not commit.

**The standing clauses hold:** emit computed objects, never prose and never hand-typed algebra; emit
operands and residual **before** any guard; no verdict tag, no conclusion in a tag name; physical symbols
combined by hand **only** in the actions and the ansatz; ⛔ emission never conditional on a payload's
**value**.

---

## G1 — ⛔⛔ A WRONG VALUE: `q_form` can silently no-op

`q_form` re-expresses a root in the scalar `q` by substituting the pattern `kx**2+ky**2+kz**2 -> q`.
⛔ **When the root is not a bare multiple of that exact `Add`, the substitution matches nothing and
returns the root unchanged — and no guard notices.** `HOMOGENEITY_DEFECTS = q*diff(root,q) − root` then
degenerates to `−root`, reporting "dispersive" merely because the symbol `q` failed to appear.

⚠ **Measured, in the committed output.** For the bare-field control the emitted defect **retains the full
wavevector dependence of the root**, which an independent computation shows should have cancelled — and
the companion `TRANSVERSE_ROOTS_Q` tag still shows `kx, ky, kz` where the substitution was meant to have
put `q`. ⛔ **The emitted value is wrong.**

### The fix — ⭐ test homogeneity by SCALING, never by substituting a symbol

Introduce a positive scalar `lambda_scale` and compute, for each root **as a function of the wavevector
components**:
```python
scaled   = root.subs({kx: lambda_scale*kx, ky: lambda_scale*ky, kz: lambda_scale*kz}, simultaneous=True)
residual = sp.simplify(scaled - lambda_scale**2 * root)
```
⭐ This vanishes **exactly** when `ω²` is homogeneous of degree 2 in the wavevector, and needs no
substitution and no pattern match.

Emit, for **every** root of **every** action: `scaled`, `lambda_scale**2 * root`, and their **residual**,
as three tags. ⭐ Independent operands, so this residual is not tautological.

⭐ **Keep the `q`-based tags too**, but ⛔ **guard the substitution**: emit a computed object recording
whether the wavevector symbols still occur in the substituted expression, so a failed substitution is
**visible in output** rather than silently poisoning a downstream value. ⛔ Emit the computed occurrence
data, not a word.

## G2 — the polarisation tests are blind to the SIGN of `ω²`

`E1`, `E2` and `E4` at a root with `ω² < 0` are bit-identical to those at a root with `ω² > 0`. Under the
sign-flipped stiffness control, **two exponentially growing modes carry the same signature as two
propagating waves.** ⚠ Nothing false is emitted — the sign is visible in the root — but no tag's *value*
is the positivity of `ω²`, so a consumer reading only the polarisation tags cannot distinguish a wave
from an instability, and that control cannot fail through them.

⭐ **Fix:** for every root of every action, emit a computed object carrying the **sign of `ω²`** under the
script's assumption set, as its own tag. ⛔ Not a word, ⛔ no check that it is positive, ⛔ no comment.

## G3 — ⛔ a late failure discards every already-computed physics tag

The control loop `break`s on the first dimension-walk failure and the script then emits **only** the
diagnostic and exits 1, throwing away `all_outputs`. ⚠ **Measured: a failure induced in one control
reduced output from 590 tags to 1**, destroying the main action's fully-computed matrices, determinant,
roots and polarisation tests — none of which had anything to do with the failure.

⛔ This is *"a guard converts an informative value into a crash"* at package granularity, and it is the
pattern this rebuild exists to remove.

⭐ **Fix: emit everything already computed, THEN fail.** On a dimension-walk failure, finish emitting every
successfully computed package and the diagnostic for the failing one, and only then exit non-zero.
⭐ **Emit first, stop second** — a reviewer ablating the script must still see the informative values.

## G4 — the premises sidecar omits the bare-field control

The `.premises` sidecar declares the supplied inputs (`LAGRANGIAN`, `PLANE_WAVE_ANSATZ`, `ASSUMPTIONS`)
for the main action and controls X1–X7 but **omits all three for X8**, so that control's constructed
action would be counted as a derived result. ⭐ Add them.

---

## Report back — under 20 lines

1. Exit status, runtime, tag count, all-unique and no-untagged confirmation.
2. **G1:** the literal scaling residual for the **main** action and for the **bare-field** control, and
   the literal occurrence object for the bare-field control. ⛔ Do not interpret them.
3. **G2:** the literal sign tags for every root of the **main** action and the **sign-flipped** control.
4. **G3:** the tag count emitted when a control's dimension walk fails. ⛔ Do not interpret it.
5. **G4:** the sidecar lines added.
6. What you could not do, and ⭐ anything that surprised you.
