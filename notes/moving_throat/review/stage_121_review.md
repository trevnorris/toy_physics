# Review: Stage 121 — Normalized mouth gain family

**Batch:** 16 — Core-to-Mouth Gain
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage121_normalized_mouth_gain_family.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage121_normalized_mouth_gain_family_sympy_audit.py`

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
**Notes Derivation Review:** M_s = Sigma_0. K_s g_q - lambda g_s = g_s sqrt(K_s K_q)(g_c-r) verified by substitution. M_q = -Sigma_0(g_c-r)^2/(1+r^2) verified by algebra. R_q = (g_c-r)^2/(1+r^2). On compensation branch g_c = r ± sqrt(1+r^2)/2: R_q = 1/4 exactly. M_s = 4 Sigma_m, M_q = -Sigma_m recovers Stage 118 closure. Mouth law Pi = Sigma_0[1-R_q S_q(Pi)] correct.
**Script Review:** Stage 120 raw gains, normalized substitutions, R_q ratio, both compensation branches asserted = 1/4. All genuine. All pass (exit code 0).
**Issues Found:**
1. **(MINOR)** Notes line 45: LaTeX typo `=rac{` should be `=\frac{`.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The normalization map is correct: `M_s = Sigma_0`, and substituting `lambda = r sqrt(K_s K_q)` together with `g_q = g_c g_s sqrt(K_q / K_s)` gives `M_q = -Sigma_0 (g_c - r)^2 / (1 + r^2)`.
2. The mixed-to-shell ratio `R_q = (g_c - r)^2 / (1 + r^2)` is therefore exact, and on the compensated core-balance surface `g_c = r ± sqrt(1 + r^2)/2` it collapses to `1/4` on both branches.
3. That immediately recovers the Stage-118 outlet-consistent closure `M_s = 4 Sigma_m`, `M_q = -Sigma_m`, so the note’s main conceptual claim is sound.

**Script Review:**

The script correctly normalizes the Stage-120 gain pair, forms `R_q`, and checks both compensation branches symbolically. I reran the audit locally and it matches the saved output.

**Issues Found:**

1. **(MINOR)** The note has a small LaTeX typo in the displayed `R_q` definition (`=rac{` instead of `=\frac{`), but the intended formula and downstream derivation are clear.

---
