**CLEAR FOR THIS BOUNDED REFERENCE-GRAZING BUILD**

I found no mathematical or validation blocker. I read the source and JSON only and ran nothing. I re-derived the algebra by hand against the saved operands and the native c1/c2 text.

## What I checked against source

- **FIRST_SHAPE (`worker.py:174-184`).** After the on-shell dispersion substitution, the numerator reduces to `iμ(qo−qi)/qo·η h + μk/(qi qo)·σ j`. `eta_bg` and `sigma_W` stay independent, and the edge momenta 1/5 and 1/10 enter as 1/20. Both lower modes (grades (1,0) and (0,1)) are selected by `mode_value` from the saved list and join the same targets (`:189-194`). The flat delta coefficient `μ/qo` also joins.
- **Closure (`:195-216`).** The defining identity `(I+aZ)F=Z` is checked on all nine entries, and the first-row closed forms satisfy it by hand. The `omega=3` anchors are checked entry by entry for both faces. The saved coefficient `(100+30i)/109` equals `a(3)`, and `β(3)=(30+9i)/109`.
- **Trace (`:226-234`, `c2:521-524`).** The four native assignment fragments are `exec`'d on the saved `valueCoefficient=1` and `height_constant=0`. `T[0,2]=0`, so the row-0 residual `T·ref−P` closes. The `middle-trace` and `input-trace` zeros tie `T[0,1]` and `T[1,2]` to `i·eta·q·h` for both signs. `normal-binding` ties the saved ±`I q_o` to `I·face·qo`.
- **Mixed terms (`:241-265`).** `phs`, `psh` and `tracehs` match my hand derivation of (2) and (7). `phs+tracehs=C` holds because `β=aμ`. The three-leg reference matrix `ref` has the right ordering (`ref[1,2]` is updated before `ref[0,2]`).
- **PV and contact (`:308-332`).** I re-derived `qm²−qi²=−t(t+2k)` and the density `J` from `transformed_psh`. Both contacts vanish at the right points. The odd PV piece is removed, and the left-height contact `C·(W/4)j(Q)` survives.
- **Sheet and census (`:334-348`, `:358-367`).** The middle-sheet symbol census is exact (`k`, `increment_transfer`, `increment_effective_bulk_speed`). The `<`/`>` signs of the branch conditions are handled correctly. The gap certificates and the `|q|≤5` joins check by hand.
- **Controls (`:418-453`).** The point `cs=√(10/7)`, `k=3/2`, `l=2`, `m=30/13` gives `qi=2`, `qo=3/2` and `qm=25/26` from the native sheet. The nonzero checks are exact Gaussian-rational components, not thresholds. The wrong-lower coefficient comes from `trace_subtraction_coefficient(corruptLowerT, …)` and is consumed by `wrongLower`.
- **Bounds.** The `0.16>0.0901` bound for μ, the `|1−iΩτ|≥1` bound for `a`, the tail polynomials (`z²/2+6z`, `7/8z²+12z`) and the 200 constant (`159<200`) all check out. The Hölder constants and the tail inequalities follow the method.
- **Evidence order, helpers and launcher.** Operands and raw residuals are saved before each guard. The helper ASTs are exact copies, and the saved copies are re-hashed. The gate chain covers worker, manifest, helpers, launcher, method record, build record and authority. Containment enforces 4 GiB, no swap, one CPU and 32 tasks, and nothing sets a deadline. The launcher arms the hook first, and `RUN.mkdir(exist_ok=False)` allows only one launch.

## Non-blocking observations

None of these needs to hold up the run. Each is fail-closed or already disclosed.

1. **Final-slot check is circular (`:277`).** `action=T[0,1]/(i·face·qm)` makes the result equal `ref[0,col]` by construction. The plan says so, and the sign content is covered by `middle-trace`. A non-circular version would use `faces[label]['height'].xreplace({k:m})` and a separate `i·face·qm`.
2. **Lower-height control (`:432`).** It negates the actual lower `T[0,1]` instead of re-running the native `trace_three` fragment on `-labheight`. The two are identical by linearity, so this is optional.
3. **`qm²` rule (`:321`).** The substitution `qm²→qi²−t(t+2k)` is stipulated and tied to the radicand identity only by narrative. A symbolic square check, `native_middle.args[0].expr**2 − realrad`, would bind it.
4. **Possible spurious failure at `:317-320`.** If `sp.together` leaves an uncancelled `(qm+β)` factor, the numerator can reach `qm³`. The degree guard then fails. That guard is fail-closed, but it would spend the single authorized run on bookkeeping. A `sp.cancel` of the numerator, or multiplying through by the known denominator, avoids this.
5. **Gate hardening.** `verify_gate` does not require `g['sharedGuard']` and `g['supervisor']` to equal the executed command paths. It also does not require their hashes to equal the `sourcePins` entries. The launcher's literal command check plus the pins make this safe in practice. `codex_job_watch.py` is not pinned.
6. **Direct tag.** `directTag` is joined to the saved density by name and argument map only. The `R(qo)·D·R(qi)` structure is not asserted algebraically.
7. **Weak points that are declared.** The dimension check is coarse. The Hölder `|p+k|≤7` constant and the exponential, monotonicity and L1 limits are assessed mathematics, as the worker labels them.

## What a pass would establish

- Both faces' native FIRST_SHAPE height and slope coefficients, with `eta` and `sigma` independent.
- Both unrestricted closure identities with their `omega=3` anchors.
- The native-derived generic trace maps and reference kernels (4)–(6) and (8).
- The uncombined mixed assignments, the surviving contact, the PV factorization to `J`, and the middle-sheet and elementary-bound joins.
- Four sensitive controls.
- A once-tagged restored direct addend.

All of this is per unit native source on `cs∈[1,2]`, `k,l∈[-3,3]`.

It would not establish any of the following:

- A source-field contraction, or the full slab operator.
- A new Fourier map.
- A finite inverse, or the cutoff-4/6 domain.
- κ=0 or β=0.
- Any integral value, or a machine proof of the limits.
- Differentiability, a defect result, or a loss claim.

This is source-only clearance, not result acceptance. Restoration, exact residuals, control nonzeroness and resource success are still unproven until the run.