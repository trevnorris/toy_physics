# Review: Stage 037 — Robin softening support lane

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage037_robin_softening_support_lane.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage037_robin_softening_sympy_audit.py`

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

1. **Equation-level correctness.** Robin BVP psi''+k^2 psi=0 with psi'(0)=h psi(0), psi'(L)=0 gives characteristic equation k tan(kL)=h, dimensionless y tan(y)=eta. Softening factor A_K = 1/(1-x/4+xy^2/pi^2) with x = pi^2 T_X/(L^2 K_W^{eff}), range 0 < x < 4. Endpoints: A_K(y=pi/2)=1, A_K(y→0)=4/(4-x). Pure-softening threshold zeta_req <= 4/(4-x). Robin eigenvalue threshold y^2 <= (pi^2/x)(1/zeta_req - 1 + x/4) verified by rearranging A_K >= zeta_req.

2. **Logical flow.** Clean: Robin BVP → stiffness formula → softening factor → endpoints → thresholds.

3. **Physical interpretation.** Weaker mouth pinning = stronger support softening. Soft-mouth limit as closure endpoint.

**Script Review:** Robin BVP from mode solution. Characteristic equation verified. A_K in (x,y) form verified. Endpoints by substitution and limit. Pure-softening floor by sp.solve. All genuine, all pass (exit code 0).

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

- The Robin-compliance derivation is correct. Solving the boundary-value problem yields the characteristic equation `k tan(kL) = h`, and the dimensionless form `y tan y = eta` follows exactly.
- The softening factor `A_K = K_W^(eff) / K_(phi,0)^(eff)` is also correct. Rewriting it in terms of `x = pi^2 T_X/(L^2 K_W^(eff))` gives the stated `1 / [1 - x/4 + x y^2/pi^2]`, with the D/N and soft-mouth limits landing at `1` and `4/(4-x)` respectively.
- I checked the threshold algebra independently. The pure-softening rescue criterion `zeta_req <= 4/(4-x)` and the corresponding `y_req^2` / `eta_req` expressions both reduce cleanly from `A_K >= zeta_req`.
- Physically, the stage does what it claims: weaker mouth pinning increases support softening, but only within a finite operator-controlled window.

**Script Review:**

- The script correctly derives the Robin characteristic equation, the dimensionless reduction, the softening factor in both `x,y` and endpoint forms, and the pure-softening threshold.
- The saved output matches the note and I did not find a bug or a trivial self-check.

**Issues Found:**

None.

**Questions:**

None.

---
