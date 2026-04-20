# Review: Stage 072 — Family1 minimal isotropic verdict

**Batch:** 11 — Loading Ratio & Isotropic Verdict
**Status:** Verified (current PASS after Family-1 replay, 2026-04-20)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage072_family1_minimal_isotropic_verdict.md`
- **SymPy:** `scripts/moving_throat_pde_stage072_family1_minimal_isotropic_verdict_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage072_family1_minimal_isotropic_verdict_mathematica_audit.wl`

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

**Notes Derivation Review:** Compares rho_alpha = 4/3 against Family-1 window. Margins all ~2.13 (huge). Zero-transport-bias result: zeta_req = 1/3 < A_F1 ~ 1.0001, so demand met before any transport bias. Strong conclusion. Remaining gap: whether actual geometry realizes minimal isotropic module. Ordering 1.333 < 3.466 < 3.468 trivially correct.

**Script Review:** Exact rho_min = 4/3, zeta_min = 1/3 via sp.Rational. Five margins computed. Four ordering assertions verified. All pass (exit code 0).

**Issues Found:** None.

### Agent: GPT-5 — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The comparison against the frozen Family-1 windows is correct: `4/3` is far below the loading-ratio thresholds and `1/3` is below both `A_F1` and the hard support ceiling.
2. The regime classification is consistent with Stage 71 and the earlier support-ratio maps. The minimal isotropic branch stays in the symmetric-lowest-twin window and does not require transport bias.
3. The conclusion is appropriately bounded. Stage 72 shows the explicit Family-1 branch is comfortably compatible with the minimal isotropic precursor, but it does not claim to solve the deeper branch-realization problem.

**Script Review:**

The script faithfully performs the intended arithmetic checks: it compares the minimal isotropic ratio against the explicit Family-1 thresholds, verifies the strict ordering, and confirms `zeta_min < A_F1`, which implies zero transport bias on the explicit Family-1 map. The saved output matches the notes.

**Issues Found:** None.

**Questions:** None.

---

### Agent: GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

1. The checkpoint claim is narrow and exact: once the minimal isotropic
   precursor is accepted, the explicit Family-1 support/source branch succeeds
   already at zero transport bias. The note does not overclaim beyond that
   arithmetic verdict.
2. The theorem packet is clean: `rho_min = 4/3`, `zeta_min = 1/3`, strict
   ordering against the carried Family-1 loading window, and the decisive
   `zeta_min < A_F1` zero-bias check.

**Script Review:**

1. The SymPy audit passes and checks exactly the intended verdict surface:
   loading-ratio ordering, symmetric-lowest-twin regime membership, the
   zero-bias `A_F1` comparison, and the hard ceiling check.
2. The Mathematica mirror passes and independently replays the same explicit
   comparisons with the carried Family-1 thresholds.
3. This is not pretending to be a fresh support-window derivation. It is a
   closed arithmetic theorem conditional on the upstream minimal-isotropic
   module, and the current dual-CAS surface is strong enough for that narrower
   claim.

**Issues Found:** None.

**Questions:** None.

---

### Agent: GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

1. The previous provenance weakness is now closed. The stage no longer relies
   on naked Family-1 decimals for `A_F1`, `zeta_max`, or the carried `rho_*`
   thresholds.
2. The script now rebuilds `A_F1` from the Stage-62 root/closed form, derives
   `zeta_max^(F1) = A_F1 pi^2/4`, and reconstructs the `rho_*` values from the
   Stage-63/69 `zeta_F1(Pe)` and `Q(zeta;0)` formulas at `lambda_mu = 1`.

**Script Review:**

1. SymPy now contains explicit anchor checks:
   `zeta_max - A_F1 pi^2/4 = 0`,
   `Q(zeta;0) - (1+zeta) = 0`,
   and the resulting `rho_suff`, `rho_fail`, `rho_max` bridge checks.
2. The Mathematica mirror now follows the same derived route with `FindRoot`
   for `y_F1` and explicit `Q(zeta;0)` anchor checks.

**Issues Found:** None.

**Questions:** None.

---
