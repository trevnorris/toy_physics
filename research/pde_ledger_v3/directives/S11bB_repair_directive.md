# S11b-B — REPAIR PASS on the two audit scripts

⛔⛔ **SCOPED REPAIR. Fix exactly the defects listed below and NOTHING else.**

Both engines built, ran, and were compared. Most of the physics agreed. ⭐ **The comparison is what found
these** — that is the architecture working, not a sign the scripts are unsound.

## Files

- `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bB_interface_assembly_mathematica_audit.wl`
- `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bB_interface_assembly_sympy_audit.py`

The specification both implement:
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md`

⚠ Run each script after editing. `math -script <path>` and `python3 <path>`; each must complete under
**10 minutes** and the `.py` must exit 0.

## ⛔⛔⛔ WHAT THIS PASS MUST NOT DO — read before touching anything

⚠⚠ **A previous repair pass in this project MADE A SCRIPT WORSE.** It added 183 lines whose checks were of
the form `assoc["key"] === variable` where the association had been built as `"key" -> variable` — i.e.
`x === x` — because the repair directive demanded self-referential verification while forbidding
independent re-derivation.

⛔ **Do NOT add any check whose expected value is computed by the same expression it is checking.**
⛔ **Do NOT add checks, harnesses, wrappers, or reporting machinery of any kind.**
⛔ **Do NOT restructure, rename, reformat, or "clean up" anything.**
⛔ **Do NOT change any value that the two engines already agree on.**
⭐ **Every edit must be traceable to exactly one numbered defect below.** If an edit is not, revert it.

## THE DEFECTS

### ⭐⭐ R1 — SymPy: B5's roots were never solved for. **This is the important one.**

`S11BB_ROOTS`, `S11BB_IMAGINARY_PART`, `S11BB_ROOT_STABILITY_CLASS` and the two artifact-diagnostic tags are
emitted **definitionally** ("for every growing root, ..."), so the script never produces a root, a sign, or
a stability classification. ⇒ **B5's central mandate is undelivered**, and the step's growth verdict
currently rests on one engine.

⭐ **Fix:** actually solve the dispersion relation on a tractable slice and report **concrete** roots — each
root, the **sign of its imaginary part**, its **decay/growth classification**, and the condition on the
moduli and `C` separating the two.
⚠ You may choose the slice (a `k = 0`, impermeable, zero-reciprocal-coupling reduction is tractable in
closed form). ⭐ **State the slice explicitly** and ⛔ **state plainly that it is a slice**, not the general
dispersion.
⛔ **Report whatever you find.** A growing root is a first-class, reportable outcome — ⛔ do not suppress
one, do not re-branch to remove one, do not add a stability assumption. If you find one, also report the
diagnostics the specification requires alongside it (§0, B5).

### R2 — Mathematica: `S11BB_FACE_RESPONSE_MU_COEFF` is wrong by a factor of `ρ_br⁰/ρ_m`.

The emitted coefficient is dimensionless, which **contradicts this script's own** `WL_DIM_FACE_P_OVER_MUTHETA`
and its own `WL_DIM_ROUTE_KIND_MU_S`. The specification defines `μ_s ≡ μ_θ/ρ_br⁰`.
⭐ **Fix the emitted coefficient and its accompanying prose so they are consistent with the script's own
dimension tags.**
⚠ **The assembly is correct** — the rows and determinant are right. ⛔ Do not touch them; this is confined
to the reported line.

### R3 — SymPy: the control-A note asserts a channel its own expression does not contain.

`S11BB_CONTROL_NO_THICKNESS` states that a pressure-feedback `τ_V` dependence survives when the thickness
channel is removed. ⚠ Holding `δW = 0` sets the face velocity to zero, so the velocity-coupled term enters
nothing. ⭐ **Verify this against your own expression** — evaluate your control-A result at two different
`τ_V` and compare — then correct the note to match what your expression actually does.

### R4 — SymPy: `S11BB_STABILITY_CONDITION` overstates its own boundary.

It reports the positivity boundary of the **unconstrained** 3-field Hessian as the conservative stability
boundary. ⚠ On the constrained slice the in-plane row forces one field to zero, so the actual boundary is
the Hessian evaluated **on the constrained direction** — a strictly weaker condition. ⭐ **Report both, and
say which is which**; ⛔ do not present the stronger one as the boundary.

### R5 — Mathematica: printed determinant disagrees in sign with the script's own cofactor formula.

The emitted determinant equals `+W₀ ×` the other engine's; the cofactor expression printed alongside it
would give `−W₀ ×`. ⛔ Harmless where the determinant is set to zero, but it is an internal inconsistency
the script's own `VERDICT` did not detect. ⭐ **Make the printed formula and the printed expression agree.**

## Output

1. The two edited scripts, each run to completion.
2. A report **under 30 lines**: per defect, what you changed and the new emitted value. ⭐ For **R1**, give
   the roots, their imaginary-part signs, the classification, and the separating condition.
3. ⭐ Confirm explicitly that you added no checks and changed no value the engines already agreed on.

⛔ Then stop and exit. ⛔ Do not write a summary document, and ⛔ do not edit the specification.
