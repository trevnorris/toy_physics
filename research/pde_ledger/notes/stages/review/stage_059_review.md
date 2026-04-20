# Review: Stage 059 — N5 wall depth lock

**Batch:** 9 — Wall Branch & Family-1 Geometry
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage059_n5_wall_depth_lock.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage059_n5_wall_depth_lock_sympy_audit.py`

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

**Notes Derivation Review:** n=5 EOS: h = m c_s^2/4 (factor 4 correct for n=5). Substituting mu_* = lambda_mu h_w into Theta_w: lambda_mu^2 m^2 rho_w^2 c_sw^2/(4 hbar^2). Healing-width: Theta_w = lambda_mu^2 rho_w^2/(16 ell^2). With ell/a=1/20: Theta_w = 25 lambda_mu^2 rho_w^2. All verified.

**Script Review:** All genuine substitution checks. All pass (exit code 0).

**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The `n=5` enthalpy identity and the wall-depth reductions are correct end to end. The healing-lock step gives `Theta_w = lambda_mu^2 rho_w^2/(16 ell^2)`, and the normalized reference-branch substitution reproduces the stated `25 lambda_mu^2 rho_w^2` factor.

**Script Review:**

The audit script faithfully implements the same chain and the saved output matches the notes exactly. The `ell = hbar/(2 m c_sw)` substitution and the final normalized branch check both reduce to zero as expected.

**Issues Found:**

None.

---
