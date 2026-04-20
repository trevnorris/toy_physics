# Review: Stage 029 — Tracking branch bounds

**Batch:** 4 — Kernel Continuation
**Status:** Verified (1× MINOR, 1× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage029_tracking_branch_bounds.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage029_tracking_branch_bounds_sympy_audit.py`

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

**Notes Derivation Review:**

1. **Equation-level correctness.** G_tr and F_tr correctly from Stage 028. Flat recovery at R=1 correct (9+2=11). Loading monotonicity dG_tr/dR = -36R xi^2(delta+xi)/[...]^2 verified by quotient rule — negative for R>0. Normalization monotonicity dF_tr/dR factored with polynomial P_R having all positive coefficients — positive for physical range. Exact comparison G_tr-G_flat and F_flat-F_tr with factored forms including (1-R^2) and (1-R) factors correct. Endpoint formulas: G_tr(R=0)=xi, F_tr(R=0)=1/(1-xi) verified. Bounds G_flat <= G_tr <= xi and 1/(1-xi) <= F_tr <= F_flat consistent with monotonicity directions. Residual comparison E_tr - E_flat = F_flat - F_tr > 0 correct.

2. **Logical flow.** Clean: tracking functions → monotonicity → comparison → endpoints → bounds → residual theorem.

3. **Physical interpretation.** Sound: R < 1 tracking branch requires less loading but also produces smaller normalization — a double-edged effect.

**Script Review:**

**B.1-B.7.** Faithful: endpoint limits, exact derivatives (sp.diff vs closed forms), flat-branch comparison formulas all verified. Polynomial coefficients P_R, P1, P2 transcribed and verified against SymPy computation. No tautologies. All pass (exit code 0).

**Issues Found:**

1. **(MINOR)** Script does not programmatically verify that all coefficients in P_R, P1, P2 are positive — it verifies the formula correctness but not the sign claim. Positivity is manifest by inspection of coefficient lists, but a programmatic check would strengthen the audit.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

- The tracking-branch formulas are correct. I independently rechecked `G_tr` and `F_tr` and the endpoint limits `R=0` and `R=1`; the algebra matches the saved output exactly.
- The derivative identities are also right. `dG_tr/dR` is strictly negative for physical `R>0`, while `dF_tr/dR` factors through a polynomial with positive coefficients, so the monotonicity claims are justified.
- The flat-branch comparison formulas are exact, and the sign conclusions follow correctly from the factorization structure. The residual ordering `E_tr > E_flat` on `0 < R < 1` is therefore sound.

**Script Review:**

- The script faithfully checks the endpoint formulas, derivative identities, flat-branch comparison formulas, and the strong-split bounds.
- I did not find a coding bug or a tautological check. The output matches the note, and the algebra is nontrivial enough to support the branch ordering theorem.

**Issues Found:**

None.

**Questions:**

None.

---
