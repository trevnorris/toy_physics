# Review: Stage 236 — Physical calibration and material thresholds

**Batch:** Batch 28 — Physical Calibration Companion
**Status:** Hardened (dual-CAS symbolic pass plus primitive-input numerical-stress coverage, 2026-04-20)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage236_physical_calibration_and_material_threshold_companion_from_the_stage235_export_and_cold_survival_compiler_sympy_audit.md`
- **Script:** `scripts/moving_throat_pde_stage236_physical_calibration_and_material_threshold_companion_from_the_stage235_export_and_cold_survival_compiler_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage236_physical_calibration_and_material_threshold_companion_from_the_stage235_export_and_cold_survival_compiler_mathematica_audit.wl`
- **Numerical stress:** `scripts/numerical/stage236_material_threshold_stress.py`
- **Numerical stress (Mathematica):** `mathematica/numerical/stage236_material_threshold_stress.wl`

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

Stage 236 is correctly scoped as a calibration companion rather than a new PDE
theorem stage. Its job is to translate the exact Stage-235 microscopic export
and cold-survival outputs into condensed-matter screening formulas while keeping
the physical unit map and the transport projection factor explicit. The
derivation cleanly separates:

1. the exact symbolic threshold compilers,
2. the legacy Session-V recovery slice, and
3. the benchmark-only numeric readback used to recover the previously reported
   Session-V material formulas.

**Script Review:**

The SymPy audit already covered the symbolic calibration identities, the legacy
slice recovery, the harmonic trigger/stiffness compiler, the Korringa ceiling,
and the material-screening ratios. I added a Mathematica mirror at
`mathematica/moving_throat_pde_stage236_physical_calibration_and_material_threshold_companion_from_the_stage235_export_and_cold_survival_compiler_mathematica_audit.wl`
that checks the same symbolic path in a second CAS and keeps the decimal
quantities confined to an explicit benchmark-only specialization section.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

Stage `236` is still correctly framed as a calibration/screening companion
rather than a new PDE theorem stage. The important point is that the threshold
stack is exact, while host-specific pass/fail outcomes depend on supplied
physical constants and calibration choices.

**Script Review:**

The symbolic SymPy and Mathematica audits still carry the theorem path. The new
shared-stress layer now adds explicit host-screen probes in both CAS layers:

1. a raw microscopic all-pass slice,
2. an exact legacy-boundary slice that replays the Session-V calibration
   recovery,
3. an isolated thermal-failure slice,
4. an isolated geometric-failure slice, and
5. an isolated stiffness-failure slice.

Those probes are driven by declared target margins `Pi_ep`, `Pi_chi`, `Pi_k`,
and `Pi_t`. The harnesses then derive `\lambda_{\rm ep}\omega_D`,
`r_{\rm turn}`, `k_{\rm eff}`, and `T` from the exact Stage-236 threshold
formulas rather than hard-coding those candidate values directly.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-20 (hardening pass)
**Verdict:** PASS after hardening

**Issue Closed:**

The prior numerical stress layer was circular: it reconstructed host inputs
from chosen `Pi_*` targets and then replayed those targets. That made the
targeted fail slices structurally weak even though the symbolic theorem path
was still useful.

**Fix Applied:**

The SymPy and Mathematica stress harnesses now read primitive host inputs
directly:

1. `\lambda_{\rm ep}\omega_D`,
2. `r_{\rm turn}`,
3. `r_{\rm turn}^{\rm phys}`,
4. `|V'_{\rm red}|(r_{\rm turn})`,
5. `k_{\rm eff}`,
6. and `T`.

From those raw inputs, the harnesses:

1. verify the physical dictionary (`r_{\rm turn}^{\rm phys}` and `K_{\rm turn}`
   consistency),
2. derive the Stage-236 thresholds forward,
3. compare the direct screening ratios against their threshold-ratio forms, and
4. check genuine pass/fail classifications for all four `\Pi_*` screens.

The new sample bank includes:

1. an all-pass microscopic host,
2. a legacy-boundary replay,
3. a dedicated electron-phonon failure,
4. a dedicated geometry failure,
5. a dedicated stiffness failure,
6. and a dedicated thermal failure.

**Issues Found:**

None after the harness rewrite.

---
