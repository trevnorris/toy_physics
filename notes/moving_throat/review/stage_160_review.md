# Review: Stage 160 — Weak axisymmetric outgoing slippage

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage160_weak_axisymmetric_outgoing_slippage.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage160_weak_axisymmetric_outgoing_slippage_sympy_audit.py`

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

### Agent: Claude Opus 4.6 — 2026-04-03
**Verdict:** PASS
**Notes Derivation Review:** Weak-axisymmetric transport: m_r, i_r, h_r inherit eps*lambda_A prefactor. Portwise sigma_r = 2m_r + 2I/(1+I) i_r + 2H/(1-H) h_r correct. Grouped collapse Xi_1 = sum rho_r sigma_r. Trace/anomaly: Xi_bar=0, b_Xi=3 a_Xi verified. P_1/P_0 = Xi_1 identification correct on conservative-shape branch.
**Script Review:** Series expansion of log(M_p/M) etc. verifies slope factorization. Full Sigma_exact from Lambda definition vs sigma_r formula. Trace/anomaly checked. Prefactor slope identity. 10 expect_zero. All genuine. All pass (exit code 0).
**Issues Found:** None.

---

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS
**Notes Derivation Review:** The weak-axisymmetric outgoing-slippage chain is consistent end to end: the microscopic slopes `m_r`, `i_r`, and `h_r` match the perturbation definitions, the collapsed `sigma_r` formula is correct, the grouped `Xi_bar`, `a_Xi`, and `b_Xi` relations check out, and the conservative-branch identification `Xi_1 = P_1 / P_0` is valid.
**Script Review:** The audit script recomputes the log-slopes of `M_r`, `I_r`, and `H_r`, verifies `Sigma_{A,r} = lambda_A sigma_r`, checks the grouped trace/anomaly identities, and confirms the outgoing-prefactor slope relation. The saved output matches the note exactly.
**Issues Found:** None.
