# Review: Stage 222 — Rigid-mouth physical normal form

**Batch:** Batch 24 — Rigid-Mouth Normal Form
**Status:** Hardened and verified (dual-CAS PASS, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage222_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md`
- **Script:** `scripts/moving_throat_pde_stage222_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage222_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`

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

The Stage 222 note does the right structural job. It turns the rigid-mouth
packet from a triangular microscopic bookkeeping device into a physical
two-coordinate normal form in `(U,V)`, then proves that the dependent-plane
compiler and the orbit-lock condition are exactly Cartesian in those variables.
The claim boundary remains correct: this is an exact rigid-mouth reduction
theorem, not a proof of full nonlinear branch realization.

**Script Review:**

The SymPy audit checks the core theorem path:

1. the diagonal physical logarithmic chart,
2. the target-ratio reconstruction,
3. the complementary physical projectors and commuting finite legs,
4. the exact physical-to-microscopic dependent-plane compiler,
5. the split between transfer-shape and dressing images,
6. the static-only, post-static, and full orbit-restoring corrections,
7. support-blindness, and
8. the Cartesian orbit-lock theorem plus its first-order form.

I added a Mathematica mirror at
`mathematica/moving_throat_pde_stage222_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`
that checks the same symbolic chart/compiler logic in the second CAS.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

The Stage `222` claim boundary is now matched more honestly by the executable
audits. This stage is still a rigid-mouth physical normal-form and correction
compiler, not a new proof of nonlinear branch realization, but the hardened
audits now reconstruct the carried inputs from the right upstream sources:
Stage `221` for the physical branch chart and Stage `219` for the dependent
plane.

**Script Review:**

The hardened SymPy and Mathematica audits now verify:

1. the Stage-`221` transfer-shape identity
   `R_target T^2 = Lambda0 (1 - eps_eta)` on the carried branch formula;
2. the rigid-mouth reduction
   `q_nt = U`, `q_eta = V` from the Stage-`221` packet factorization rather
   than by free-symbol renaming;
3. the Stage-`219` dependent-plane compiler propagated into the physical chart,
   including a left inverse rebuilt from the dependent-plane equations
   `Delta_Keta = -q_eta`, `Delta_mu = q_nt - q_eta`;
4. the exact transfer/dressing microscopic images and the static / dressing /
   full orbit-restoring correction split;
5. propagation of Stage-`221` support-blindness through `U`, `V`, and the
   correction packet using explicit support-blind derivative hypotheses on the
   carried observables, rather than by differentiating log-of-free-symbols;
6. the Cartesian orbit-lock equivalence by solving the dependent-plane defect
   equations for `(U,V)` and by checking the physical equilibrium equations
   `T^2 = T_ref^2`, `eps_eta = eps_eta_ref`, rather than by `Exp[0] - 1`
   tautologies;
7. the first-order compiler and the static-blind equal-drift dressing ray as
   before.

This closes the main defect from the earlier version: the old script proved
only that already-independent symbols stayed independent, while the new script
proves the Stage-`222` statements as inherited propagation facts from the real
Stage-`219/221` formulas.

**Issues Found:**

None in the hardened Stage `222` path.

---
