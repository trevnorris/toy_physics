# Review: Stage 047 — Equilibrium alignment

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage047_equilibrium_alignment.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage047_equilibrium_alignment_sympy_audit.py`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Script runs without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template:

### Agent: [Model Name] — [Date]
**Verdict:** [PASS | MINOR | ISSUE | BLOCK]

**Notes Derivation Review:**

**Script Review:**

**Issues Found:**

**Questions:**

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** PASS

**Notes Derivation Review:**

1. **Equation-level correctness.** Equilibrium response chi_sigma = g_phi chi_phi/H(y) from static linearization. Overlaps: O_sp = g_phi I_1, N_ss = g_phi^2 I_2 correct. Coherence C^2 = I_1^2/(N_pp I_2) — g_phi^2 cancels. Cauchy-Schwarz I_1^2 <= N_pp I_2 verified (apply CS with measure chi_phi^2 d^3y). Equality iff 1/H constant on supp(chi_phi). Constant-compressibility: I_1 = N_pp/H_w, I_2 = N_pp/H_w^2 → C^2 = 1. Structural identity Theta_eq = Lambda_eq = g_phi^2 I_1 on equilibrium branch verified. G_eq = g_phi^2 I_1/K_X. Contact with Stage 45: matches G_micro at C^2=1, rho_*=rho_w.

2. **Logical flow.** Clean: equilibrium response → overlaps → coherence → Cauchy → constant-H → softening → contact.

3. **Physical interpretation.** Key insight: source profile is fixed by parent equilibrium, not a free choice.

**Script Review:** Overlap invariants, constant-H limit, discrete two-point Cauchy gap, eliminated-source softening. Minor: Theta_eq and Lambda_eq defined identically in script (lines 113-118), making ratio check partially tautological. All other checks genuine. All pass (exit code 0).

**Issues Found:** None (minor tautology noted but non-blocking).

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

The equilibrium-alignment derivation is correct: the source profile is fixed by `1/H(y)`, the overlap invariants reduce as claimed, and the constant-compressibility limit gives `C^2 = 1` and the matched-layer gain used to connect back to Stages 45-46. The discrete two-point Cauchy check is a good independent sanity test.

**Script Review:**

The script is mostly faithful and the saved output matches the notes, but the final eliminated-source softening check is partly tautological because `Theta_eq` and `Lambda_eq` are both assigned the same constant-layer expression before the `Lambda^2/Theta` comparison. Earlier checks do verify the real overlap structure, so this is a coverage issue rather than a mathematical blocker.

**Issues Found:**

Minor script-faithfulness gap in the final softening check.

---
