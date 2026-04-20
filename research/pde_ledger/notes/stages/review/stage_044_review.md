# Review: Stage 044 — Microscopic gain thresholds

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage044_microscopic_gain_thresholds.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage044_microscopic_gain_thresholds_sympy_audit.py`

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

1. **Equation-level correctness.** Factorization Xi_micro = kappa G_micro with G_micro = chi_sigma Lambda_phi^2/K_X verified. Gain thresholds G_fail = Pe_req/(kappa Delta_inf), G_suff = Pe_req/(kappa Delta_0) correct. Threshold surfaces for chi_sigma and Lambda_phi^2 correct. Soft-support limits: Delta_0 → 1/2, Delta_inf → 1 verified by expanding hyperbolic functions around alpha=0. Compliant-mouth limits: Delta_0 → (1-sech alpha)/kappa, Delta_inf → tanh(alpha)/alpha verified. Stiff-support asymptotics within compliant-mouth: G_fail → Pe_req/sqrt(kappa), G_suff → Pe_req. Physical: compliant mouth saturates sufficiency threshold at O(Pe_req).

2. **Logical flow.** Clean: G_micro factorization → exact thresholds → soft-support limits → compliant-mouth limits → combined asymptotics.

**Script Review:** Factorization check, gain thresholds, threshold surface consistency, soft-support limits (sp.limit), compliant-mouth limits (sp.limit), combined asymptotics. All genuine computations. All pass (exit code 0). Most comprehensive script in the batch.

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The microscopic factorization `Xi_micro = kappa G_micro` is algebraically correct, and the derived thresholds `G_fail = Pe_req / [kappa Delta_inf]` and `G_suff = Pe_req / [kappa Delta_0]` follow directly from Stage 42.
2. The threshold surfaces for `chi_sigma` and `Lambda_phi^2` are consistent rewrites of the same inequalities, with the fail/succeed ordering preserved.
3. The limiting results are correct: `Delta_0 -> 1/2`, `Delta_inf -> 1` as `kappa -> 0+`, the compliant-mouth limits match the hyperbolic closed forms, and the stiff-support asymptotics reduce to the stated scalings.

**Script Review:**

The audit script performs genuine limit computations and threshold substitutions rather than symbolic restatements. The reported pass/fail surfaces, soft-support limits, compliant-mouth limits, and asymptotic checks all agree with the notes and the output exits cleanly.

**Issues Found:** None.

---
