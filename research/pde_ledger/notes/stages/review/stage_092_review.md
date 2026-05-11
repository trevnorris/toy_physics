# Review: Stage 092 — Overlap isotropy

**Batch:** 1 — Geometry Lift & Coupling
**Status:** Verified (current PASS after isotropic overlap closure, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage024_overlap_isotropy.md`
- **SymPy:** `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Mathematica script faithfully implements notes
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

**1. Equation-level correctness.**

- **STF harmonics (Sec. 1):** Normalization Y_A = sqrt(15/(8 pi)) E_A^{ij} n_i n_j with Tr(E_A E_B) = delta_AB. Fourth-moment formula int dOmega n_i n_j n_k n_l = (4 pi/15)(delta_ij delta_kl + perms) is standard isotropic result. Substituting and using trace-free property gives delta_AB. Correct.

- **Source-map identity (Sec. 2):** mhat_ang = 1 follows from orthonormality. Correct.

- **Overlap-integral factorization (Sec. 3):** O(3) invariance forces angular overlaps diagonal in Y_A basis via Schur orthogonality / Wigner-Eckart theorem. All lanes become identical. Isotropic coefficients re-quoted from Stage 023 match exactly (with A-lane index dropped).

- **Axisymmetric splitting (Sec. 4):** Triple-overlap matrix M^(20)_AB = int dOmega Y_A Y_20 Y_B verified diagonal with pattern diag(1, 1/2, 1/2, -1, -1) * kappa_* where kappa_* = sqrt(5)/(7 sqrt(pi)). Grouped defect algebra: with splitting (1, 1/2, -1), formulas a_x = eps*x1/4 and b_x = 3*eps*x1/4 follow from standard definitions. b = 3a law confirmed.

- **Normalization transport (Sec. 5):** P_A = P_0 + eps lambda_A P_1 with P_1 = (N_1 D_0 - N_0 D_1)/D_0^2 is standard quotient rule. Grouped defects inherit b = 3a law. Correct.

**2. Logical flow.** Clean from Stage 023. Angular closure is the natural first consequence of O(3) structure.

**3. Assumptions.** All explicit: O(3) invariant reference throat, separable ansatz, weak axisymmetric perturbation as single Y_20 background.

**4. Completeness.** Isotropy theorem closes angular part only; radial/axial problem correctly identified as remaining. Non-axisymmetric perturbations noted as producing different patterns.

**5. Notation consistency.** H_r introduced for g_U^2 + g_W^2, replacing Stage 023's G_{A,r} (avoids Newton's G clash). Change motivated but not explicitly flagged in notes.

**6. Physical interpretation.** Clear: orthonormality closes angular normalization gap, O(3) forces lane collapse, first symmetry-breaking produces diagnostic fingerprint.

**Script Review:**

**B.1 Faithful.** Covers STF Gram matrix (I), angular source-map (I), isotropic collapse (II), triple-overlap matrix from sixth moment (III), splitting weights and b=3a law (III), first-order transport (IV).

**B.2 No bugs.** STF basis matrices correct. `pairings` function generates all 15 complete pairings. I4 and I6 implement correct isotropic moment prefactors. `triple_overlap` sums correctly.

**B.3 Hardcoded values.** 4pi/15 and 4pi/122 are standard sphere moment prefactors. kappa_* derived from computation, not assumed.

**B.4 Tautological checks.** Section II (isotropic collapse: setting x20=x21=x22=x0 and checking a=b=0) is trivially true by definition — weakest check. Section III (25 triple integrals via sixth moment) is genuinely non-trivial. Section IV (transport via epsilon expansion) is non-trivial.

**B.5 Symbol assumptions.** D0, D1, N0, N1 nonzero+real (appropriate). STF basis uses exact SymPy rationals.

**B.6 Output agreement.** All pass. Gram = 5×5 identity. Triple-overlap = kappa_* * diag(1, 1/2, 1/2, -1, -1). Transport law verified.

**B.7 Coverage.** O(3) factorization theorem itself (Schur lemma argument) not amenable to SymPy verification — consequences tested instead. Inherited Z_n, N_n not re-verified (acceptable).

**Issues Found:**

1. **(MINOR, notation)** H_r replacing G_{A,r} is not explicitly flagged in the notes. A brief note about the symbol change would help cross-stage readers.

2. **(MINOR, tautological)** Section II of script is trivially true by definition. Would be stronger to verify non-equal lanes produce non-zero defects.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The angular part is solid: the normalized real STF basis is orthonormal, the angular source-map identity is exact, and the `O(3)` isotropy result correctly collapses the grouped lanes to one common value.
2. The axisymmetric `Y_20` triple-overlap matrix and the resulting grouped splitting law are correct, and the first-order normalization transport formula preserves the `b = 3 a` defect line exactly.

**Script Review:**

1. The script reproduces the harmonic Gram matrix, the `Y_20` overlap matrix, the grouped splitting weights, and the first-order transport law. The cached output matches the note claims.
2. I also ran a small independent SymPy check for the `b - 3 a` relation and it returned zero, matching the stage result.

**Issues Found:**

1. **(MINOR)** The `H_r` symbol substitution from Stage 6 is not called out explicitly in the note, so cross-stage readers have to infer that the rename is cosmetic.
2. **(MINOR)** The isotropy theorem is validated through its consequences, but the radial/axial overlap integrals themselves are still inherited from Stage 6 rather than recomputed here.

**Questions:**

None.

---

### Agent: GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

1. The checkpoint claim is the angular closure layer of the grouped
   `20/21/22` system: normalized real-STF harmonics, exact angular source-map
   identity, isotropic lane collapse under an `O(3)`-invariant kernel, the
   exact `Y_20` triple-overlap matrix, the grouped splitting signature
   `(1,1/2,-1)`, and the inherited `b = 3 a` transport law. I did not find a
   mathematical mismatch between the current note and the live CAS outputs.
2. The note now explicitly says that the radial/axial overlap objects are the
   carried Stage-6 reduced-overlap layer, and it explicitly flags `H_r` as the
   renamed combined gauge/mixed strength chosen to avoid confusion with
   Newton's `G`. That removes the earlier notation caveat.

**Script Review:**

1. The earlier “tautological Section II” caveat is stale. The current SymPy
   audit still checks equal-lane collapse, but it also includes unequal-lane
   witness tests showing that the grouped defects become nonzero away from the
   isotropic locus.
2. The Mathematica mirror independently replays the same audit surface:
   harmonic Gram matrix, exact angular source map, unequal-lane grouped
   witnesses, exact `Y_20` triple-overlap matrix, grouped splitting weights,
   and first-order transport of `P_A = N_A / D_A`.
3. I did not find a remaining checkpoint-level trust defect. The stage is
   intentionally reduced in the radial/axial direction, but that scope limit is
   explicit and matches the theorem claim.

**Issues Found:**

- None. The older notation and tautology caveats are both closed in the live
  note/script surface.

**Questions:**

- None.

---

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

The remaining stage-level gap was the radial/axial overlap output rather than
the angular closure. That gap is now closed at the correct reduced level:
finite witness modes reconstruct the boxed isotropic overlap coefficients and
the low-frequency conservative compiler without pretending Stage `024` is
re-deriving the full Stage `023` continuum from scratch.

**Script Review:**

1. Both CAS layers now verify representative instances of the boxed overlap
   formulas
   `C_\alpha = \lambda_{B,\alpha} I_{\eta\alpha}`,
   `B_n = \sum_\alpha C_\alpha^2 / \varpi_\alpha^{2n}`,
   and the one-pair reduced
   `Z_n`, `N_n`, `D(\omega) = D_0 + D_2 \omega^2 + D_4 \omega^4`
   decomposition.
2. The audit remains intentionally finite-mode at this stage, but that is
   enough to close the referee-visible gap that the boxed overlap outputs had
   no executable verification at all.

**Issues Found:**

- None.

---
