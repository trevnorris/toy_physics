# Review: Stage 002 — Breathing reduction

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage002_breathing_reduction.md`
- **Script:** None (status/consolidation stage)

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

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** PASS

**Notes Derivation Review:**

**A.1 Equation-level correctness**

- **Y_00 normalization (Section 1):** `Y_00 = 1/(2 sqrt(pi))`, `int_{S^2} Y_00^2 dOmega = 1`, and `(1/4pi) int_{S^2} Y_00 dOmega = 1/(2 sqrt(pi))` — all correct. Confirmed by SymPy audit Section I.1.

- **Normalization bridge (Section 1):** `q_00(0,t) = 2 sqrt(pi) delta_a(t)`. Consistent with Stage 001 Section 4 and SymPy audit Section I.2. Correct.

- **Two-mode ansatz (Section 1):** `eta_0(w,t) = 2 sqrt(pi) [alpha_a(w) delta_a(t) + alpha_L(w) delta_L(t)]`. The factor `2 sqrt(pi)` matches the `q_00 = 2 sqrt(pi) delta_a` normalization. Verified that mouth average of `eta` gives the correct identification.

- **Reduced Lagrangian (Section 2):** The angular integral uses `int_{S^2} Y_00^2 dOmega = 1`. Traced the mass matrix: `M_{AB} = 4 pi int dw mu_eta alpha_A alpha_B` — the factor `4 pi` arises from `(2 sqrt(pi))^2 * 1 = 4 pi`. Stiffness matrix similarly correct. Confirmed.

- **Euler-Lagrange equations (Section 3):** `M_{AB} ddot{Q}^B + K_{AB} Q^B = 0` is the standard result for a quadratic Lagrangian. Confirmed by SymPy audit Section VI.

- **Grouped P2 reduction (Section 4):** For one-mode truncation `eta_{2m} = beta_2(w) q_{2m}(t) Y_{2m}(Omega)`, effective mass `M_2 = int dw mu_eta beta_2^2` and stiffness `K_2 = int dw [T_w beta_2'^2 + (K_eta + 6 T_Omega) beta_2^2]`. The factor `l(l+1) = 6` for `l = 2` is correct. No extra `4 pi` factor because `int (Y_{2m}^{real})^2 dOmega = 1` directly. Confirmed by SymPy audit Section VII.

**A.2 Logical flow** — Tightly structured: fix normalization → insert ansatz into action → derive E-L equations → match to old closure → repeat for P2 → summarize. Each step follows from the preceding one.

**A.3 Assumptions** — Two-mode truncation explicitly stated as "minimal truncation that can recover the old (a,L) closure." Higher axisymmetric modes deferred. One-mode P2 truncation explicit. P2 degeneracy qualified with "before any symmetry breaking or matter/gauge coupling."

**A.4 Completeness** — Correctly identifies what the next coupled calculation must do (Section 5): renormalize matrices, split or preserve P2 degeneracy, shift pole data, produce passive/outgoing odd parts.

**A.5 Notation consistency** — `K_0(w) = K_eta(w)` matches Stage 001 `l = 0` case. Metric factor `sqrt(gamma_0)` absorbed into integration measure, standard practice. `Q^A = (delta a, delta L)` notation clean.

**A.6 Physical interpretation** — Sound. Distributed wall field reduces to old `(a,L)` closure in lowest-mode truncation (key continuity check). P2 degeneracy on isotropic background is a consequence of geometry action alone. Correct identification of what breaks degeneracy.

**Script-backed status:** References `moving_throat_pde_master_sympy_audit.py` which passed (exit code 0) and covers Sections I through VII with substantive `expect_zero` assertions using `euler_equations`, `simplify`, `diff`, `integrate`.

**Issues Found:**

None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

- The normalization bridge is correct. I independently rechecked `Y_00 = 1/(2 sqrt(pi))`, the mouth average `1/(2 sqrt(pi))`, and the resulting `(2 sqrt(pi))^2 = 4 pi` factor that appears in the reduced `l=0` mass and stiffness matrices.
- With that normalization, the two-mode axisymmetric truncation does reduce to the expected conservative matrix form `M_AB ddot Q^B + K_AB Q^B = 0`, so the main continuity claim back to the old `(a,L)` closure is sound.
- The grouped-real `P2` one-mode reduction is also structurally correct: the common uncoupled stiffness shift is `l(l+1) T_Omega = 6 T_Omega`, so the isotropic reference branch gives a degenerate `P2` bundle before any symmetry breaking or matter/gauge dressing.
- The note also correctly identifies what the next coupled stages must change: matrix renormalization, pole shifts, and possible splitting of the grouped `P2` block.

**Issues Found:**

- **[MINOR] The reduction silently drops the surface-weight factor inherited from Stage 001.** Stage 001 wrote the wall action with an explicit `sqrt(gamma_0)`, but the reduced matrices here are written as `4 pi int dw ...` and `int dw ...` without either carrying `sqrt(gamma_0)` into the overlap formulas or stating that the coefficient functions have been densitized. That does not break the qualitative conclusion, but the note should explicitly choose one convention; otherwise the reduced matrices and the Stage 001 modal operator are only correct up to that unstated redefinition.

---
