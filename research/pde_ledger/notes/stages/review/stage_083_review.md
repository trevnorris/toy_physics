# Review: Stage 083 — Outgoing normalization factorization

**Batch:** 13 — Outgoing DtN
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage083_outgoing_normalization_factorization.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage083_outgoing_normalization_factorization_sympy_audit.py`

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

**Notes Derivation Review:** One-pole ansatz with chi_Q outgoing factor. Kbar_2, Kbar_4 correct at omega^2, omega^4. Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega_Q^5) correct at omega^5 imaginary part. Even ratios = N_Q, odd ratio = chi_Q N_Q. Observable condition mhat_0^2 chi_Q N_Q = 1. All verified.

**Script Review:** Series expansion, coefficient extraction, all five identities checked. Genuine non-tautological. All pass (exit code 0).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The one-pole retarded ansatz with an outgoing renormalization factor `chi_Q` is used correctly. Its low-frequency expansion preserves the conservative `omega^2` and `omega^4` structure while multiplying the odd `omega^5` coefficient by exactly `chi_Q`.
2. The factorization is clean: the even ratios are still governed only by `N_Q = Kbar_0 / Kbar_0^target`, while the odd ratio becomes `chi_Q N_Q`.
3. The observable normalization law `mhat_0^2 chi_Q N_Q = 1` follows directly once the source-map factor is included, so the last reduced 2.5PN defect is properly separated into conservative and outgoing pieces.

**Script Review:**

The script faithfully expands the retarded branch, extracts the `omega^2`, `omega^4`, and imaginary `omega^5` coefficients, and checks the even/odd target ratios plus the normalized odd-condition solution for `N_Q`. I reran the audit locally and independently spot-checked the `9 chi_Q /(32 Omega_Q^5)` odd coefficient.

**Issues Found:** None.

---
