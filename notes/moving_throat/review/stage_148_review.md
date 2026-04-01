# Review: Stage 148 — Exact branch drifts

**Batch:** 18 — Linear Defect Transport & Final
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage148_exact_branch_drifts.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage148_exact_branch_drifts_sympy_audit.py`

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
**Notes Derivation Review:** D/N co-transport d ln LW = d ln a. Fixed-r channel solves for dv: (1/2) d ln(Z_q c_s^2 c_sw^3/(rho_w a^5)). Fixed-g channel solves for dT. Exact amplitude laws with numerical prefactors. Product/ratio factorization. n=5 wall-EOS reduction to 4 irreducible drifts (Z_q, rho_w, c_s, a). Stage-147 delta_perp vanishes by construction. All hand-verified.
**Script Review:** Two constraint equations via sp.solve for dv, dT. Product/ratio, n=5 reduction, back-substitution all checked. All pass (exit code 0). Minor: D/N law check is dead code; numerical prefactors not verified.
**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The lower-branch drift laws are correct. The D/N co-transport law, the fixed-r and fixed-g solved drifts, the product/ratio factorization, and the `n=5` reduction all line up with the note. The collapse of the Stage-147 off-family scalar to zero on the exact lower branch is the right conclusion.

**Script Review:**

The script solves the two-channel drift system, verifies the product and ratio laws, and checks the `n=5` reduction and the vanishing of the off-family channels. The saved output matches the note.

**Issues Found:**

None.

---
