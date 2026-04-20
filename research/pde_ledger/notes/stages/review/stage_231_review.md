# Review: Stage 231 — Dynamic event-chain compiler

**Batch:** Batch 27 — Relaxed Dynamic Event Chain
**Status:** Hardened and verified (dual-CAS symbolic + dual-CAS numerical PASS, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md`
- **Script:** `scripts/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- **Numerical stress:** `scripts/numerical/stage231_event_chain_stress.py`
- **Numerical stress (Mathematica):** `mathematica/numerical/stage231_event_chain_stress.wl`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Mathematica script faithfully implements notes
- [ ] Scripts run without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template: -->

### Agent: Codex GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

The Stage 231 note keeps the dynamic claim boundary correct. It does not claim
to produce a new same-charge interaction law; it promotes the already-lowered
Stage-230 stationary front end into the reduced event-chain objects used by the
Session-II audit. The key outputs are the exact energy integral, the finite-
radius threshold-speed formulas, the Coulomb turning-point/WKB reference, the
near-top parabolic action law, and the dynamic turning-point diagnostics.

**Script Review:**

The SymPy audit covers the right split:

1. the symbolic event-chain conservation and threshold/WKB formulas,
2. the exact Coulomb antiderivative and near-top action normal form, and
3. the declared Session-II benchmark specialization used for the reported
   numerical readback.

I added a Mathematica mirror at
`mathematica/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
that keeps the symbolic theorem path separate from the benchmark-only numeric
layer and checks both in the second CAS.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

The checkpoint claim remains correctly scoped. Stage `231` does not promise a
new microscopic same-charge law; it compiles the reduced dynamic corridor,
subbarrier transport packet, Coulomb comparison object, and turning-point
diagnostics carried by the relaxed branch.

**Script Review:**

The symbolic SymPy and Mathematica audits are still the main theorem checks, but
they no longer stand alone numerically. The new shared-stress layer now probes:

1. three admissible event-chain families with explicit threshold-window and
   turning-point checks,
2. direct Coulomb-action quadrature against the closed-form formula in both CAS
   layers, and
3. three independent near-top parabolic probes that compare direct quadrature
   against the exact top-action law.

The JSON inputs are probe-only stress parameters, not hidden theorem constants.
All threshold speeds, Coulomb turning points, action values, transmission
ratios, and `\lambda_{\rm th}` readbacks are derived inside the harnesses from
those probe families.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

The stage claim is now checked more literally than before. The symbolic layer
still keeps the benchmark-only numeric banner separate from the theorem path,
but it now also solves the threshold-speed equations directly, verifies the
Coulomb endpoint evaluation from the antiderivative, and checks the implicit
turning-point transport law
\(dr_\pm/dE = 1 / V'(r_\pm(E))\).

**Script Review:**

The numerical stress harness is no longer replay-based.

1. It starts from raw lowered-barrier parameters for an explicit one-peak
   barrier family rather than from `I_new/I_coul` targets.
2. It numerically locates the peak and both turning points from the barrier
   itself.
3. It computes `I_new` by direct quadrature, compares the Coulomb quadrature to
   the exact closed form, and verifies transmission improvement forward from
   those actions.
4. It derives `\Xi_{\rm turn}` from a monotone probe profile and derives
   `\lambda_{\rm th}` from `E / |V'(r_+)|` using the same barrier.
5. It verifies the turning transport law numerically by a high-order
   finite-difference check.
6. It includes genuine expected-failure probes: one with no admissible peak and
   one with multiple stationary points, both of which are detected as failures
   by both CAS layers.

The Mathematica numerical mirror now executes the same raw-case classification
flow as the Python harness, so the stress layer is dual-CAS rather than
Python-only.

**Issues Found:**

None in the hardened Stage `231` path.

---
