# Directive pathA_16 — Chunk B2c remediation #2: honest report framing (no false precision; lead with the rigorous verdict)

**Date:** 2026-06-19
**Owner:** Codex (codes/designs). Claude reviews only.
**Builds on:** `pathA_15` (remediation #1). The verdict LOGIC, the assembly/dual-engine math, the warm-start code,
the negative controls, and the no-sub-floor-extrapolation guard all passed fidelity + adversarial audit and are
CORRECT — **do NOT change them.** This is a **report/emission framing fix only.** The underlying scientific
conclusion (an unnatural knife-edge MISS) is right; the report over-claims precision it does not have and
under-states the result that IS rigorous.

## Why (the adversarial findings to fix — all presentation/over-claim, severity MAJOR)
1. **The rigorous verdict is buried.** The fit-independent, *measured* result is not the headline.
2. **False precision.** `τ_crit/τ*` are reported to 12 digits with a τ* "error bar" implying ~5 sig figs, but they
   are extrapolations ~3 orders below the data; three reasonable fits span ~3 orders of magnitude.
3. **A failed numerical push is dressed as success.** Warm-start converted ZERO new converged points (the deepest
   point converges COLD); the descent is ~3% (0.03→0.029), ~3× short of the required `1e-3–1e-2` band.
4. **"Numerical not Schur" is under-evidenced** — no Jacobian/D0 was measured at the failing τ.

## F-requirements (acceptance in §A)

### F1 — Lead the Verdict with the RIGOROUS, measured, fit-independent result
Compute and emit (these need only the already-converged points + algebra — NO new solves, NO extrapolation):
- **Measured deficit:** `P0` vs target (`54/5=10.8`) at EVERY converged τ — e.g. `P0=2.79e-9` (τ=1) … `1.22e-6`
  (τ=0.029) → a **~7-to-9 order-of-magnitude deficit on the natural side**, directly measured.
- **Algebraic knife-edge (fit-independent):** since `P0=N0/D0` exactly, any stable root has
  `D0(τ*)=N0(τ*)/10.8`. Because `N0≪K` at every measured point (`N0/10.8 ≈ 2.5e-8` vs `K=O(0.1–8)`), **any stable
  root is a knife-edge cancellation `D0(τ*)≪K`, regardless of where τ* lands.**
- **Conclusion to state plainly:** the τ=1 "structural win" (no fragile cancellation) does **NOT** survive
  calibration; reaching the GR anchor requires riding the `D0→0` stability boundary, i.e. it is not a legitimate /
  natural calibration → a MISS (→ κ_PV in a later chunk, never a knob here).
This is the headline of the `## Verdict` section.

### F2 — Remove the false precision; report τ_crit/τ* only as a bound + order-of-magnitude range
- Delete the 12-digit `τ_crit`/`τ*` presentation and the `local_tau_star…` row in the Error Bars table that implies
  sig figs.
- Replace with: (a) the rigorous **upper bound** `τ* < τ_floor (≈0.029)`; and (b) an explicit statement that a
  drift-aware estimate is **fit-dependent across ~3 orders of magnitude** — record the actual spread from at least
  the available fits (frozen-deepest-coeff ≈3e-5; wide-pair τ=1↔0.03 power-law ≈7e-4; close-pair fit degenerate /
  above floor) and state τ\* is **NOT pinned** by this run. Machine-record the spread in the bundle (a small
  `edge_estimate_fit_spread` block), not just prose.

### F3 — Soften "numerical not Schur" to match the evidence actually held
State only what is supported: no `D0`/Jacobian-smallest-mode was measured at the failing τ (the solve never
converges there); the deepest *converged* `D0=O(0.1)` and `K/(B0+Z0)=849`, and the failures are smooth line-search /
continuation stalls with residuals rising monotonically as τ drops. Conclude: **consistent with a numerical
(continuation/conditioning) floor, NOT a diagnosed physical marginal-stability edge** — do not assert a Schur/soft-mode
measurement that was not taken.

### F4 — State the numerics status honestly (no success-dressing)
The report must say: the directive's `1e-3–1e-2` target band was **NOT reached** (deepest converged τ≈0.029, ~3×
above the band); the warm-start, while functional, **converted no new converged point** (the deepest point converges
cold); the ~3% floor descent is **not a material improvement** over the original 0.03 floor. Frame as a known
limitation. (Pinning τ* exactly would need a small-τ preconditioner/linear-solver effort — note it as deferred and
**not required** for the verdict.)

### F5 — Minor: note the κ̂ provenance
The bundle stores `κ̂=7.7209353105` (analytic continuum) while the seed `K=7.7209443582` comes from the B2a
discretized eigen-solve (~1e-6 rel difference = discretization). Either use one consistently or add a one-line note
documenting the analytic-vs-discrete distinction. Immaterial to the verdict; just don't leave it silently
inconsistent.

### F6 — Do NOT touch the certified pieces
No change to: the verdict-logic in `calibrate_from_records` (frozen-prior-free), `assemble_primary/secondary` and the
dual-engine, the warm-start / wall-predictor code, the negative controls, the sub-floor extrapolation guard, the
spot-convergence check. This is report text + the small machine-recorded fit-spread block ONLY. Target-blind
firewall stays intact; no new DOF. All `test_patha*` tests still pass (update only if the report-string assertions
need it).

## A. Acceptance criteria
1. `## Verdict` leads with the measured deficit + the algebraic knife-edge (F1), stating the structural win does not
   survive.
2. No 12-digit τ_crit/τ*; replaced by the bound + fit-spread range, machine-recorded, with τ* declared not pinned
   (F2).
3. "Numerical floor" wording softened to the held evidence (F3).
4. Numerics status stated honestly (band not reached; warm-start converted no new point; ~3% descent) (F4).
5. κ̂ provenance noted/reconciled (F5).
6. Certified logic/math/warm-start untouched; tests pass (F6).

## B. Review stack (orchestrator, after Codex reports clean)
1 adversarial pass (distrust-all-clean): confirm the report is now honest (rigorous verdict foregrounded, no false
precision, no success-dressing, numerical-vs-physical not over-claimed) AND that the certified verdict logic /
assembly math / warm-start / negative controls are byte-for-byte unchanged in behavior. (Full fidelity re-audit not
needed — no math changed; the reviewer confirms the no-change.)
