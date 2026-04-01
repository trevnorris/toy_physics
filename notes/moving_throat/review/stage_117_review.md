# Review: Stage 117 — Family1 mouth fixedpoint

**Batch:** 15 — Positive Source & Mouth Dynamics
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage117_family1_mouth_fixedpoint.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage117_family1_mouth_fixedpoint_sympy_audit.py`

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
**Notes Derivation Review:** Two-channel specialization: kappa_s=0 gives S(Pi,0)=1 (static shell, verified by limit). kappa_q=pi/2 (D/N half-wave). Fixed-point Pi = M_s + M_q S_q(Pi) correct. Compensation line M_s = Pi_* - S_q(Pi_*) M_q correct. S_q(Pi_*) computed to 30 digits.
**Script Review:** S(Pi,kappa) from Stage 116. Static limit via sp.limit. S_q by substitution. Numerical evaluation at Pi_*. All genuine. All pass (exit code 0).
**Issues Found:**
1. **(MINOR)** S_q(Pi_*) displayed as ...429 in notes but 30-digit computation gives ...428494, so correct 15-sig-fig rounding is ...428. Cosmetic 1-ULP discrepancy, no downstream impact.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The Family-1 reduction is faithful to Stage 116: the static shell lane contributes the exact unit response, and the mixed lane is the first D/N half-wave at `kappa_q = pi/2`.
2. The fixed-point law `Pi = M_s + M_q S_q(Pi)` and the canonical gain line `M_s = Pi_* - S_q(Pi_*) M_q` are consistent with the note and the saved audit output.
3. The stage is an appropriate narrowing step. It does not claim a final PDE proof; it just turns the mouth-layer solve into a two-gain fixed-point problem.

**Script Review:**

The script faithfully computes the static limit, the mixed-lane response at `kappa = pi/2`, and the canonical gain line at `Pi_*`. The saved output matches the note, and the numerical evaluation is nontrivial.

**Issues Found:** None.

**Questions:** None.

---
