# Review: Stage 122 — Family1 actual mouth gains

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage122_family1_actual_mouth_gains.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage122_family1_actual_mouth_gains_sympy_audit.py`

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
**Verdict:** MINOR
**Notes Derivation Review:** Family-1 geometric inputs from Stage 104/117. Fixed-point M_s = Pi_*/(1-R_q S_q) correct. Natural branch: R_q^{nat} ~ 0.1455, M_s^{nat} ~ 1.669, M_q^{nat} ~ -0.243 all verified. Compensated branch: g_c from lower branch, R_q = 1/4 exactly (algebraic identity confirmed to 25 digits), M_s^{comp} ~ 1.806, M_q^{comp} ~ -0.451 all verified. Shell gain shift ~8.2% correct.
**Script Review:** All formulas from notes implemented. R_q = 1/4 assertion genuine. 30-digit precision. All pass (exit code 0).
**Issues Found:**
1. **(MINOR)** Notes Section 4: mixed gain magnitude ratio stated as ~1.86097 but script computes ~1.86028. Confirmed by manual check (0.4515/0.2427 ~ 1.860). Qualitative conclusion ("about factor 1.86") unaffected.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The branch formulas are correct: `M_s = Pi_*/(1 - R_q S_q(Pi_*))` and `M_q = -R_q M_s` reproduce the reported natural and compensated Family-1 gains.
2. The compensated branch reduction to `R_q = 1/4` is exact, and the self-consistency with the earlier compensation family is intact.
3. The only issue is the numeric ratio in the prose: the text says `|M_q^{comp}|/|M_q^{nat}| ≈ 1.86097385480`, while the script output gives `1.86028418097`, so the magnitude comparison should be corrected.

**Script Review:**

The script computes the reference branch values directly and the `R_q = 1/4` assertion passes. I independently rechecked the reported ratio from the printed gains, and it matches the script output rather than the prose value.

**Issues Found:**

1. The mixed-gain magnitude ratio in the notes is numerically inconsistent with the script/output by about `6.9e-4`.

---
