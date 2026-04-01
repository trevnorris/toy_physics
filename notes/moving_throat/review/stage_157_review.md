# Review: Stage 157 — Static self similarity

**Batch:** Batch 19 — Grouped Outlet & Similarity Closure
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage157_static_self_similarity.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage157_static_self_similarity_sympy_audit.py`

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
**Notes Derivation Review:** Weighted decomposition delta_D = omega_K delta_K - omega_B delta_B - omega_Z delta_Z verified. Weight identity omega_K-omega_B-omega_Z=1. Wall-referenced Xi_load form correct. Microscopic log slopes for BdG (2 delta ln(c/varpi)), Z (quotient rule), N (P^2/Delta^2) sectors verified. Self-similarity theorem: all Sigma=0 → Xi_load=0. Clean.
**Script Review:** Weight identity, delta_D decomposition, wall-referenced form, all three microscopic weighted averages (BdG 2-mode, Z 2-port, N 2-port), self-similarity, mismatch-field. All genuine. All pass (exit code 0).
**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

The self-similarity reduction is coherent. The weight identity closes exactly, the `delta_D` decomposition is consistent with the grouped `K/B/Z` lanes, and the wall-referenced `Xi_load` form matches the staged microscopic averages. The theorem that common static self-similarity forces `Xi_load = 0` is the correct algebraic endpoint.

**Script Review:**

The audit script checks the weight identity, the weighted decomposition, the wall-referenced form, the microscopic averages, and the self-similarity theorem. The saved output shows each check returning zero, so the implementation is faithful.

**Issues Found:**

None.

---
