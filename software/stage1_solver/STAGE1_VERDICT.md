# Stage-1 Verdict — first physical, target-blind falsification result (effective-closure branch)

**Date:** 2026-06-17
**Status:** §7 step-9 Stage-1 verdict DELIVERED for the pre-registered effective-closure branch.
**Freeze hash:** `834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8`
**Tracked artifact:** `software/stage1_solver/frozen/m1c/834835c9…/` (frozen packet, summary, diagnostics,
WP1 backgrounds, reproducibility report).

## The test
A target-blind, dual-engine realization of the moving-throat branch: torch solves the self-consistent WP1
background → Mathematica derives the BdG+wall + ℓ=2 vector-harmonic Maxwell outgoing transfer → the existing
Python V2 chain judges it → the GR-normalized residual `R_norm = mhat0²·S_port·N0/D0 − 54Gc_s⁵/(5a⁵c⁵)`
(target `R_norm = 0`). GATE A frozen target-blind (decision-07, structured branch): `a=1, L=37/20,
R_exit=3/2, R0(w)=1+½(3x²−2x³)` cubic smoothstep, wall functions `μ_η=T_w=T_Ω=K_η=1`, natural-unit constants,
`chi_Q` extracted-not-frozen. Freeze written + hashed BEFORE the solve; no value mutated after any residual.

## The result — a robust MISS
- **`R_norm = −10.7999993`.** The derived outgoing transfer `P0 = N0/D0 = 6.6e-7` is **~7 orders of magnitude
  below** the GR target `54/5 = 10.8`. (`R_pole = −0.83`, `P2 = −2.1e-7`, `P4 = 6.4e-8`, `chi_Q ≈ 1.6e7`
  = `10.8/P0`, near-singular.)
- **Magnitude miss, NOT a sign flip** (Claude+Codex sign-trace agreed): target is `R_norm=0` not `10.8`;
  `N0=P²/Δ²` is a square (≥0, no sign flips it negative or creates a 7-order gain); the spatial-gauge sign
  flip changes `N0` by exactly 0. Only `transfer=+10.8` gives `R_norm=0`; we got `6.6e-7`.
- **Robust, not a resolution artifact:** the 91% relative §J floor is on near-null `min_density`, NOT the
  `R_norm` path; `N0` is resolved to ~few % (transfer-mesh 3.4%, interp 1.6%); absolute `R_norm` budget `6e-7`
  ≪ the `10.8` miss. Byte-reproducible (two runs identical; arbiter reproduces from the committed packet).

## Review (all 4 gates clean)
fidelity CORRECT (frozen config = decision-07; freeze-before-solve; chain faithful; `N0` genuine,
mesh-converged, >0) · adversarial CLEAN (robust miss, budget honest, freeze discipline impeccable, no firewall
touched) · arbiter (reproduces from the committed frozen packet) · sign-trace (Claude+Codex: genuine
magnitude miss).

## Interpretation — what it means (and does NOT mean)
- **NOT "the theory is bad."** This is ONE pre-registered branch of a DELIBERATELY-INCOMPLETE version (the
  effective-closure approach, which omits the matter/gauge→wall *return* coupling `S_η^(A)` by construction).
  A single-branch miss falsifies that branch+posits, not the framework. The reduced algebra is known to be
  able to hit the target with the right coefficients (the calibrated control gives `R_norm=0`).
- **NOT "find the values that fit."** Searching the free_choice values for a fit is the forbidden tuning the
  target-blind discipline exists to prevent (a fitted pass is worthless); and a 7-order structural gap does
  not close by reparametrizing posited wall functions.
- **What it IS:** the quantified signature of the DEFERRED `S_η`/Path-A coupling. The forward-only transfer
  (effective closure) is ~0; whether the FULL (self-consistent, Path-A) model reproduces GR is **STILL OPEN** —
  this test could not answer it by design. The result CONDITIONAL on the frozen target-blind posits.

## What we gained
A working end-to-end pipeline that produces real, frozen, target-blind, dual-engine, reproducible
falsification numbers — and a first result that sharply localizes the open question to one place (the deferred
return coupling).

## Open after Stage-1 (deferred)
- **Path A** — derive `S_η^(ψ,A)` (the matter/gauge→wall return coupling); the decisive test of whether the
  FULL model reproduces GR. A subprogram + a NEW pre-registration (conceptual change → user-gated). Cannot
  promise it closes the gap.
- Secondary observables `R_pole, P2, P4, chi_Q`, the WP3 d-ln targets (response-gated); higher resolution
  (GPU/cluster); a branch sweep to confirm the effective-closure gap is universal.
