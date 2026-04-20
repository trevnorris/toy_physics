# Review: Stage 125 — Selfconsistent mouth branch

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage125_selfconsistent_mouth_branch.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage125_selfconsistent_mouth_branch_sympy_audit.py`

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
**Notes Derivation Review:** Identifies g_c = g_Pi from Stage 113. R_q(Pi) = (g_Pi-r_F1)^2/(1+r_F1^2) from Stage 121. Inverted fixed-point Sigma_0 = Pi/(1-R_q S_q) correct. Self-matched That_m = sqrt(9 Sigma_0/20) from Stage 123. At compensation point Pi_*: R_q = 1/4 (algebraic identity), Sigma_0 ~ 1.806, That ~ 0.901 — all consistent with Stages 113/117/122/123 to 13-15 sig figs.
**Script Review:** Symbolic R_q(g_minus)=1/4 check. Pi_* via nsolve. All branch quantities evaluated at 30-digit precision. Output matches notes. All pass (exit code 0).
**Issues Found:** None (minor: no discussion of 1-R_q S_q singularity locus, but safely positive at canonical point).

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The self-consistent mouth-branch closure is correct. Identifying `g_c` with the explicit exponential mouth-bias map gives the stated `R_q(Pi)`, and the fixed-point relation for `Sigma_0(Pi)` and the traction map `\widehat T_m(Pi)` follows directly. The canonical compensation point reproduces the earlier outlet-consistent branch values.

**Script Review:**

The audit script faithfully checks the `R_q(g_minus)=1/4` identity, solves for `Pi_*`, and evaluates the branch quantities at high precision. The saved output matches the notes, including the canonical values at `Pi_*`.

**Issues Found:**

None.

---
