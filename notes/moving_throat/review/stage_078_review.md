# Review: Stage 078 — Second order geometry contamination

**Batch:** 12 — Geometry Lane
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage078_second_order_geometry_contamination.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage078_second_order_geometry_contamination_sympy_audit.py`

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

**Notes Derivation Review:** Schur complement D_eff = D_q - chi^2 M0^2/D_g correct. Series expansion: K0corr ~ -chi^2 M0^2/G0, K2corr ~ +chi^2 M0^2 G2/G0^2, K4corr ~ chi^2 M0^2(G0 G4-G2^2)/G0^3 — signs verified (double negative from Schur). eps_2, eps_4 both O(chi^2), so c_pole perturbation is O(chi^2). Linear expansion (1/4)(1+eps_4-2eps_2) correct. Perturbative stability established.

**Script Review:** Series of -chi^2 M0^2/D_g computed. chi^2 factoring verified. Unperturbed limit c_pole(chi=0)=1/4 checked. First derivative at chi=0 vanishes (O(chi^2) confirmed). All pass (exit code 0).

**Issues Found:** None.

### Agent: Codex GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The mixed-mode Schur-complement derivation is correct. The contamination coefficients enter with the right signs, both `eps_2` and `eps_4` start at `O(chi^2)`, and the pole-fraction deformation is therefore quadratic in the anisotropy mixing as claimed.

**Script Review:**

The audit script computes the low-frequency expansion directly and checks the first derivative of `c_pole` at `chi=0`, which is the right nontrivial test here. The saved output matches the notes and the symbolic checks pass.

**Issues Found:**

None.

---
