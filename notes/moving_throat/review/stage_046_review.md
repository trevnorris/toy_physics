# Review: Stage 046 — Parent thresholds

**Batch:** 8 — Operator & Gain Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage046_parent_thresholds.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage046_parent_thresholds_sympy_audit.py`

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

1. **Equation-level correctness.** Parent amplitude thresholds: g_{phi,fail}^2 from inverting G_micro = G_fail. Coherence-factor form: N_ss cancels when O_sp^2 = C^2 N_ss N_pp substituted. Coherence thresholds C_fail^2, C_suff^2 correct. Cauchy no-go: C^2 <= 1 gives G_max = rho_* g_phi^2 N_pp/(m c_s^2 K_X), rigorous. Expanded thresholds with kappa: K_X cancels — amplitude threshold is tension/compressibility-controlled. Ordering C_fail^2 <= C_suff^2 from G_fail <= G_suff. Matched-profile special case (C^2=1, N_ss=N_pp=O_sp=1) collapses correctly.

2. **Logical flow.** Clean: amplitude thresholds → coherence thresholds → Cauchy no-go → expanded forms → special case.

**Script Review:** Overlap-form thresholds, coherence substitution, coherence thresholds and ratio, best-case gain, kappa insertion with K_X cancellation in 4 forms. All genuine, all pass (exit code 0).

**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The parent-threshold algebra is consistent throughout: the `g_phi^2` fail/succeed surfaces follow directly from the Stage-45 gain, the coherence-factor rewrite is exact, and the Cauchy no-go statement is the correct `C^2 <= 1` consequence. The `kappa`-substituted forms also cancel `K_X` in the same way the notes claim.

**Script Review:**

The script matches the notes and the saved output. It independently checks the overlap-threshold identities, the coherence-floor formulas, the best-case `G_max` factorization, and the `kappa = K_X L^2/T_X` substitutions without relying on hardcoded final answers.

**Issues Found:**

None.

---
