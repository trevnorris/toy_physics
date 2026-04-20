# Review: Stage 204 — Resonance/linewidth tradeoff

**Batch:** Batch 23 — Same-Charge Dynamic Window
**Status:** Verified (1× PASS, 2026-04-20)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage204_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.md`
- **Script:** `scripts/moving_throat_pde_stage204_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage204_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`

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

The Stage 204 note keeps the correct claim boundary. It does not claim that
resonance solves the same-charge problem; it proves the exact local pole
tradeoff inside the already-fixed one-port mixed bundle. The key outputs are the
Breit-Wigner normal form, the wall-like specialization, the exact
dispersive/absorptive no-free-lunch theorem, and the resulting residue-to-
linewidth survival inequality.

**Script Review:**

The SymPy audit covers the right theorem path:

1. the generic simple-pole normal form,
2. the carried Stage-203 derivative identity `dD_Pi/dPi = -N(omega)`,
3. the wall-like pole specialization,
4. the exact line-shape, optimum, and low-loss identities,
5. the barrier / absorbed-power ratio and quality-factor detuning formulas,
6. the linear survival-window inequality.

I added a Mathematica mirror at
`mathematica/moving_throat_pde_stage204_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
that checks the same symbolic theorem path in the second CAS, while keeping the
constructive numeric slice explicitly probe-only.

**Issues Found:**

None.

---
