# S9 `.wl` REPAIR 3 — a wrong homogeneity value, and a sign-blind test

⛔ **Repair in place. Do not rewrite.** Two fixes, nothing else.

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`

Verify: `cd /var/projects/toy_physics && timeout 600 math -script research/pde_ledger_v3/mathematica/S9_light_requires_shear_mathematica_audit.wl`
⚠ Two-seat licence; one kernel at a time, retry on a licence error. Iterate to exit 0.

⛔⛔ **Do not read any `.py` engine.** ⭐ The two engines must stay independently constructed so they can
disagree. ⛔ Do not commit.

**The standing clauses hold:** emit computed objects, never prose and never hand-typed algebra; emit
operands and residual **before** any guard; no verdict tag, no conclusion in a tag name; physical symbols
combined by hand **only** in the actions and the ansatz.

---

## G1 — ⛔⛔ A WRONG VALUE: the homogeneity defect is computed after a substitution that can silently fail

The dispersion is re-expressed in a scalar `q` by **substituting the pattern `kx^2+ky^2+kz^2 -> q`**, and
the homogeneity defect `q·∂_q(ω²) − ω²` is then formed from the result.

⛔ **When the root is not a bare multiple of that exact pattern, the substitution matches nothing and
returns the root unchanged — and nothing notices.** The defect then degenerates to `−ω²`, which reports
"dispersive" merely because the symbol `q` failed to appear.

⚠ **Measured, in the committed output.** For the bare-field control the emitted defect **retains the full
wavevector dependence of the root**, which an independent computation shows should have cancelled.
⛔ **The emitted value is wrong.** The candidate-speed tag for the same control is a second symptom: it
carries `q` in its denominator while its numerator still shows `kx, ky, kz` — the substitution matched in
one place and not the other.

### The fix — ⭐ test homogeneity by SCALING, never by substituting a symbol

Do not rely on `q` appearing. Introduce a positive scalar `lambdaScale` and compute, for each root
`omegaSquared` **as a function of the wavevector components**:

```
scalingResidual  =  ( omegaSquared /. {kx -> lambdaScale kx, ky -> lambdaScale ky, kz -> lambdaScale kz} )
                    −  lambdaScale^2 · omegaSquared
```

⭐ This vanishes **exactly** when `ω²` is homogeneous of degree 2 in the wavevector — i.e. when the phase
speed is wavevector-independent in magnitude — and it needs no substitution and no pattern match.

Emit, for **every** root of **every** action: the scaled expression, `lambdaScale^2 · omegaSquared`, and
their **residual**, as three tags. ⭐ Genuinely independent operands, so this residual is not tautological.

⭐ **Keep the existing `q`-based tags as well**, but ⛔ **guard the substitution**: emit a computed object
recording whether the wavevector components still occur in the substituted expression, so a failed
substitution is **visible in output** instead of silently poisoning a downstream value. ⛔ Do not describe
it in words — emit the computed occurrence data.

## G2 — the polarisation tests are blind to the SIGN of `ω²`

`E1`, `E2` and `E4` at a root with `ω² < 0` are bit-identical to those at a root with `ω² > 0`. Under the
sign-flipped stiffness control, **two exponentially growing modes carry the same signature as two
propagating waves.** ⚠ Nothing false is emitted — the sign is visible in the root itself — but there is
**no tag whose value is the positivity of `ω²`**, so a consumer reading only the polarisation tags cannot
tell a wave from an instability, and that control cannot fail through them.

⭐ **Fix:** for every root of every action, emit a computed object carrying the **sign of `ω²`** under the
script's assumption set — e.g. `Simplify[Sign[root], assumptions]` — as its own tag.
⛔ Do not emit a word, ⛔ do not add a check that it is positive, and ⛔ do not comment on it. Emit the
computed sign.

---

## Report back — under 20 lines

1. Exit status, runtime, tag count, all-unique and no-untagged confirmation.
2. **G1:** the literal scaling residual for the **main** action and for the **bare-field** control, and
   the literal `q`-substitution occurrence object for the bare-field control. ⛔ Do not interpret them.
3. **G2:** the literal sign tags for every root of the **main** action and of the **sign-flipped**
   control. ⛔ Do not interpret them.
4. What you could not do, and what blocked it.
5. ⭐ Anything that surprised you.
