# Review: Stage 225 — Actual twin-support placement and coherent orbit-lock compiler

**Batch:** Batch 25 — Selected Twin-Support Placement
**Status:** Hardened and verified (dual-CAS PASS, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage225_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.md`
- **Script:** `scripts/moving_throat_pde_stage225_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage225_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`

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

The Stage 225 note keeps the right claim boundary. It does not claim to solve
the full same-charge realization problem; it compiles the already-selected
twin-support branch into the actual coherent placement variables and cleanly
separates support placement from orbit lock. The key outputs are the realized
support coordinate `varrho_phys = 2 (1 - epsilon) / 3`, the threshold rewrites
in the physical `epsilon` variable, the support-blind orbit packet, and the
direct observable compiler for `(R_tr, R_target, epsilon_eta)`.

**Script Review:**

The SymPy audit checks the correct theorem path:

1. the actual selected-twin placement formulas,
2. the threshold rewrites and selected-window inequalities,
3. the reduced-state bridge and support-blind observables,
4. the finite orbit packet,
5. the infinitesimal observable compilers and orbit-lock implication, and
6. the exact support/orbit packet split.

I added a Mathematica mirror at
`mathematica/moving_throat_pde_stage225_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`
that checks the same symbolic placement/compiler identities in the second CAS.
The explicit rational sample point is kept probe-only and is not used as part
of the proof path.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

The Stage `225` note is now matched more precisely by the executable layer.
The support-placement algebra remains the same, but the orbit side is now
checked as a direct-observable compiler carried by the coherent placement map,
not as a bundle of `zeta`-free formulas plus a single algebraic reconstruction
identity.

**Script Review:**

The hardened SymPy and Mathematica audits now verify:

1. the imported selected-branch placement identities
   `varrho_phys = 2 (1 - epsilon) / 3`,
   `sigma_phys = 2 epsilon / (1 - epsilon)`,
   with the Stage-`223` loading ratio and Stage-`010` `2/11` coefficient called
   out explicitly as imported constants;
2. the threshold rewrites and selected-window inequalities as before;
3. the reduced-state bridge together with honest `zeta`-absence checks on the
   coherent placement map itself (`epsilon`, `R_tr`, `R_target`);
4. support-blind propagation from abstract direct observables
   `(R_tr, R_target, epsilon_eta)` into the finite orbit packet, rather than
   differentiating packet formulas that simply omit `zeta`;
5. the infinitesimal direct-observable compiler
   `Theta_1 = dln R_tr`,
   `Xi_1 = -dln R_target - (eps_eta/(1-eps_eta)) dln eps_eta`,
   `R_1 = dln R_target`,
   together with the exact matrix compiler, its determinant, and its inverse;
6. the orbit-lock implication as a full invertible packet recovery:
   the zero orbit packet forces the zero direct-drift triple, rather than being
   checked by unwinding only the `Xi_1` definition; and
7. the mixed-only product law
   `R_target M_mix = C_mix` and the exact support-packet sensitivity
   `partial_zeta M_tr`.

This closes the main weakness in the earlier version. The old script had real
placement algebra, but its support-blindness checks were by omission and its
orbit-lock step reduced to solving one definition for `dln eps_eta`. The new
version proves the stage as a carried compiler on the observable packet.

**Issues Found:**

None in the hardened Stage `225` path.

---
