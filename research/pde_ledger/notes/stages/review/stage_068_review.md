# Review: Stage 068 — Softening depth normal form

**Batch:** 2 — Wall Profiles & Loading
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage034_softening_depth_normal_form.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`

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

1. **Equation-level correctness.** All verified: x := A - lambda_- substitution correct. Secular equation in x-form: alpha_0 = x(x+DK)/[kappa_0^2(x+DK) + kappa_1^2 x] verified. Overlap s_- = D^2/[kappa_0^2(x+DK)^2 + kappa_1^2 x^2] derived from 1/(dalpha/dx) with dalpha/dx = [kappa_0^2(x+DK)^2 + kappa_1^2 x^2]/D^2 verified by quotient rule. Normalization product in x-form correct. Monotonicity dalpha/dx manifestly positive (ratio of positive quantities). Loading decomposition direct from Stage 16.

2. **Logical flow.** Clean: variable substitution → differentiation → composition → Stage-16 rewrite.

3. **Assumptions.** All inherited from prior stages. Hellmann-Feynman from Stage 13.

4. **Completeness.** Edge cases x=0 (unloaded) and x→A (softening) addressed.

5. **Notation consistent** with Stage 16.

6. **Physical interpretation.** Correct: trades eigenvector data for single scalar deformation coordinate.

**Script Review:**

**B.1-B.7.** Faithful. No bugs. No hardcoded values (all symbolic). Checks: substitution equivalence (genuine), monotonicity derivative (sp.diff vs closed form — genuine), loading decomposition. All pass. Minor: check 17.3 (alpha_mix + gBreq - alpha_x = 0) is mildly tautological by construction, and no explicit s_- = 1/diff(alpha_x, x) check, but substitution check implicitly covers this.

**Issues Found:**

1. **(MINOR, script)** Check 17.3 is tautological by construction (gBreq := alpha_x - alpha_mix). Does not add verification power.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

- The softening-depth reparameterization is correct: substituting `x = A - lambda_-` turns the selected-branch secular equation into the rational form shown in the note, and the resulting `alpha_0(x)` and `N_-(x)` expressions match the saved audit output.
- I independently rechecked the monotonicity identity with SymPy. `d alpha_0 / dx` reduces to the manifestly positive rational form in the note, so the branch remains one-to-one in the stable interval.
- The support-loading rewrite is also algebraically consistent: once `alpha_mix` is separated out, the remaining required support term is exactly the difference `alpha_x - alpha_mix`.

**Issues Found:**

- **[MINOR] The final support-loading check is still self-referential.** The script defines `g_B,req^2/varpi^2 := alpha_x - alpha_mix` and then verifies `alpha_mix + g_B,req^2/varpi^2 - alpha_x = 0`, so the last assertion does not independently test the formula. That matches the existing audit gap and is acceptable only as a coverage limitation, not as a separate derivation check.

---
