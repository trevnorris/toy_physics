# Directive pathA_15 — Chunk B2c remediation #1: locate the REAL root, fix the verdict logic, push the numerics

**Date:** 2026-06-18
**Owner:** Codex (codes/designs). Claude reviews only.
**Builds on:** `pathA_14` (the B2c build). The B2c machinery (assembly/observable/error math, dual-engine,
negative controls) passed fidelity audit and is CORRECT — **do not change it.** This remediation fixes the
**verdict**, which a fidelity-clean + adversarial review found **FATAL: not established.**

## Why the current verdict is FATAL (the review findings you must fix)
1. **The verdict is decided by a discredited number.** `calibrate_from_records` sets
   `status="root_below_tau_floor"` from `tau_frozen(5.06e-6) < tau_floor(0.03)` (lines ~887-892). `tau_frozen` is the
   **frozen-background prior** (τ=1 coeffs only). The report's OWN drift table shows the background is wildly NOT
   frozen (N0 ×22, B0 ×10 over τ=1→0.03). The real root-finder (`root_find_primary/secondary`, `root_agreement`)
   is **never invoked** (`root_agreement=null`). The locational claim rests on a quantity the same report calls a
   negative control.
2. **A stable-side root almost certainly EXISTS, ~140× above the frozen prior, as a knife-edge.** Because `B0+Z0`
   GROWS as the wall softens, the real Schur edge is `τ_crit ≈ 7e-4` (not 5e-6). As `τ↓τ_crit⁺`, `D0=K−B0−Z0→0⁺`
   while `N0>0`, so `P0=N0/D0→+∞` and `R_norm→+∞`; since `R_norm(0.03)=−10.8<0`, a sign change (root) is forced in
   `(τ_crit, 0.03)`. Drift-aware estimate: `τ*≈7e-4`, `D0(τ*)≈1.3e-6`, `K/(B0+Z0)≈1.0002`, `|ln τ*|≈7.3`. The root
   sits in the band the current run could NOT converge.
3. **The "floor" at 0.03 is NUMERICAL, not physical.** At the failed τ=0.02/0.01 probes `D0≈0.15/0.077` — far from
   any physical edge. Every failed probe stalls at the SAME continuation milestone `eos_K=0.08` regardless of target
   τ → a fixed continuation schedule choking, not a loss of equilibrium. The floor was never bracketed or diagnosed.
4. **The "structural win" (no fragile cancellation) does NOT survive calibration** — at the real root D0 is a ~4-digit
   knife-edge essentially on the wall-mode marginal-stability boundary. The report's naturalness section reports τ=0.03
   (comfortable), not the actual root.

## R-requirements (acceptance in §A)

### R1 — Verdict logic must derive from REAL re-solves, never the frozen prior
- Remove any code path where the FINAL status/location is decided by `tau_frozen`/`schur_critical_tau` computed from a
  single (e.g. τ=1) coefficient set. The frozen-background root stays ONLY as the R6 negative control.
- The verdict must come from the real, drift-aware `{B,Z,N,K}(τ)` sampled at genuinely converged τ points + the
  existing validated assembly model + the existing real root-finder. `root_agreement` must be actually computed (not
  null) whenever a bracket exists.

### R2 — Compute the real, drift-aware Schur edge and trace R_norm(τ) on the resolvable stable band
- Build `{B,Z,N}(τ)` from REAL converged re-solves at a ladder of τ descending from 0.03 into the `1e-3–1e-2` band
  (and lower if reachable). Each point: converged closed background (residual ≤ `1e-6`, target the `~1e-9` class as at
  τ=0.03/1), then B2a BdG (converged modal sweep) + B2b Maxwell on it.
- From these, locate the real `τ_crit` (where `D0=K−B0−Z0` crosses 0) and trace `D0(τ), P0(τ), R_norm(τ)`.
- **Spot convergence check at the deepest converged τ** (modal K-sweep for B_n; one grid refinement for Z/N): the
  small-τ geometry differs from τ=1, so do not assume the τ=1 grid/mode adequacy carries down — confirm it at the
  bottom, log the achieved bars.

### R3 — Push the continuation numerics (route is yours; the diagnosis points the way)
The observed failure is a **fixed eos_K continuation restarting from scratch and stalling at `eos_K=0.08`** at every
small τ. Make the closed-background solve robust enough to converge deeper. Standard options (YOUR design choice, not a
mandate): natural-parameter **τ-homotopy warm-starting** from the neighbouring converged background (step τ down,
reuse the previous state as the initial guess) instead of re-running the eos_K ramp cold; adaptive continuation step;
damped/line-searched Newton with backtracking; Anderson acceleration; or pseudo-arclength near the edge. Each script
stays ≤10 min under `timeout 600` (exit 124 = reformulate, never raise the cap) — warm-starting is also the natural way
to beat the τ=1e-4 timeout.

### R4 — Locate-and-confirm OR bound-and-characterize (be explicit which, and why)
- **Preferred:** bracket the `R_norm` sign change with two genuinely converged re-solves (one with `R_norm<0`, one with
  `R_norm>0` i.e. `P0>10.8`), run the existing Brent/bisection **dual-engine root agreement**, and **CONFIRM τ\*** with a
  fresh full re-solve (`|R_norm(τ*)|≤tol`, coeffs within §J bars). Getting `R_norm>0` requires a converged background at
  `D0≈N0/10.8` (very small) — attempt it.
- **If the knife-edge is genuinely unreachable** (the background goes marginal/singular before `R_norm>0` is
  achievable): bound the root from above by the deepest converged τ, estimate `τ_crit` and `τ*` from the validated
  `{B,Z,N}(τ)` model, and **determine whether the deep failure is PHYSICAL or NUMERICAL** — the decisive diagnostic:
  does the closed-background Jacobian's smallest mode / the wall-mode `D0(τ)` collapse toward 0 as τ approaches the
  failure (→ physical soft-mode / marginal-stability edge), or does the solve fail while `D0` is still `O(0.1)` with a
  healthy Jacobian (→ numerical)? Report the evidence either way. No extrapolation of any bundle below the deepest
  converged τ.

### R5 — Naturalness + held-out reported AT the real root (not the floor)
- Report `τ*` (or its rigorous bound), `D0(τ*)`, `|ln τ*|`, `K/(B0+Z0)` and digit-cancellation **at the located/bounded
  root**. State plainly whether the "no fragile cancellation" structural win survives calibration.
- If (and only if) a τ\* is CONFIRMED (R4 preferred branch), emit the held-out `R_pole/P2/P4` predictions at τ\* with §J
  error bars. If only bounded, say so and emit the floor diagnostics labelled as such.

### R6 — Discipline preserved (unchanged from pathA_14)
- Target-blind: τ calibrated ONLY on `R_norm`; held-out reported as predictions, never tuned to; **no new DOF / no
  rescue parameter.** A confirmed knife-edge / marginal-edge / unreachable-root is a legitimate scientific MISS →
  κ_PV handled in a LATER chunk (diagnose + DERIVE missing physics), never a knob here.
- Keep the genuine dual-engine + the 3 negative controls. ADD the physical-vs-numerical edge diagnostic (R4).
- Do NOT touch the assembly/observable/error math (fidelity-clean). Do NOT touch `research/pde_audit/simulation/` or
  `physical_export_permitted`. Outputs under `runs/patha_b2c_calibration/`. YAML/markdown for human/LLM files, JSON for
  machine bundles. CPU sparse-direct (GPU off). Mathematica ≤2 concurrent seats.

## A. Acceptance criteria
1. No final-status decision derives from `tau_frozen`; the verdict traces to real converged re-solves (R1).
2. A τ ladder of REAL converged backgrounds descends into `1e-3–1e-2` (or to a documented genuine edge), each with
   logged residual; deepest-τ spot convergence (modal + grid) confirmed (R2).
3. Real `τ_crit` and `R_norm(τ)` trace reported; the verdict is either (a) located+confirmed τ\* via the real
   dual-engine root-finder + confirmation gate, or (b) bounded τ\* + an explicit PHYSICAL-vs-NUMERICAL determination
   of the deep failure with evidence (R3, R4).
4. Naturalness + (if confirmed) held-out reported AT the root, with a clear statement on whether the structural win
   survives (R5).
5. Target-blind firewall intact; no new DOF (R6).
6. Tests updated (the verdict-logic change must be covered; the frozen-prior-as-verdict path must be gone — add a test
   that the status cannot be set to a "below floor" conclusion purely from frozen coeffs when real samples exist).
   All `test_patha*` pass.
7. Report `reports/patha_b2c_calibration_report.md` rewritten to state the real verdict honestly (no naturalness
   advertised at a non-root τ; the frozen prior shown only as a negative control).

## B. Review stack (orchestrator runs after Codex reports clean)
1. Transliteration-fidelity audit (1 clean agent) on any NEW/CHANGED math (continuation route, edge diagnostic,
   verdict logic). The assembly/observable math is already certified — confirm it is UNTOUCHED.
2. Adversarial audit (≥1 round, distrust-all-clean): hunt for a new can't-fail floor, a verdict still secretly leaning
   on the frozen prior, a "converged" deep point that is actually unconverged (residual gate honesty), extrapolation
   below the floor, a physical-vs-numerical claim asserted without evidence, and any held-out leakage.
