# Review: Stage 033 — Zeta threshold comparison

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage033_zeta_threshold_comparison.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage033_zeta_threshold_comparison_sympy_audit.py`

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

1. **Equation-level correctness.** Doubling theorem S(1;eps)=2 verified (trivial cancellation). Equivalence zeta_req <= 1 ⟺ S_req <= 2 via factored form of zeta_req-1 — sign analysis correct. Impossibility tower: zeta_req > 1/(2n+1)^2 rules out nth harmonic — correct from Stage 32 suppression. Stiffness threshold x_max = [1/((2n+1)^2 zeta_req)-1]/(n(n+1)) verified by solving rational inequality. Enhancement ceiling S_n^{max} = 1+(1-eps)/((2n+1)^2-eps) correct at x=0 limit. Explicit S_1^{max} = (10-2eps)/(9-eps) verified.

2. **Logical flow.** Clean: doubling theorem → equivalence → impossibility tower → stiffness threshold → enhancement ceiling.

3. **Physical interpretation.** Sound: lowest twin is universal doubling branch (S=2), higher harmonics provide only modest enhancement. Key question (whether real PDE lands on lowest twin with zeta_req <= 1) appropriately narrow.

**Script Review:**

**B.1-B.7.** Checks doubling theorem, equivalence factored form, x_max formula (sp.solve), suppression factor, enhancement ceiling including n=1,2 values. No tautologies (each involves nontrivial cancellation/solving). All pass (exit code 0).

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The threshold comparison is algebraically exact. I independently rechecked the key boundary identities: `zeta_req = 1` occurs exactly at `S_req = 2`, and the Stage-32 same-operator twin family gives the stated `zeta_n^(twin) = 1/[(2n+1)^2(1+x n(n+1))]`.
2. The doubling theorem is correct. `S(1; eps) = 2` follows by direct substitution, and the higher-harmonic bound `zeta_req > (2n+1)^(-2)` is the right immediate no-go criterion for `n >= 1`.
3. The stiffness threshold `x <= x_max(n; zeta_req)` is the correct exact solution of `zeta_n^(twin) >= zeta_req`, and the enhancement ceiling `S_n^(max)` matches the `x=0` saturation limit.
4. Physically, the stage correctly separates a universal lowest-twin doubling branch from the strongly suppressed higher-harmonic family. No missing branch or sign issue showed up in the comparison.

**Script Review:**

1. The audit script checks the doubling theorem, the `zeta_req` boundary identities, the exact `x_max` solve, the suppression factor, and the enhancement ceiling at `n=1,2`.
2. The saved output matches the note, and the checks are nontrivial rather than tautological.

**Issues Found:**

None.

**Questions:**

None.

---
