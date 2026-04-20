# Review: Stage 035 — Nontwin asymmetry threshold

**Batch:** 6 — Support & Threshold Analysis
**Status:** Verified (2× PASS, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage035_nontwin_asymmetry_threshold.md`
- **Script:** `scripts/moving_throat/moving_throat_pde_stage035_nontwin_asymmetry_threshold_sympy_audit.py`

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

1. **Equation-level correctness.** Three-regime classification from Pi_tr vs C_mix and 2C_mix correct. zeta_req = (Pi_tr-C_mix)/[C_mix-eps(2C_mix-Pi_tr)] verified by substituting S_req = Pi_tr/C_mix. Derivative dzeta_req/dPi = C_mix(1-eps)/[...]^2 verified by quotient rule — positive for 0 < eps < 1. Excess Delta_zeta = (1-eps)(Pi_tr-2C_mix)/[...] verified. General lowest-lane zeta_0^{phys} = (K_W/K_{phi,0}) Omega_0^2 consistent with Stage-32 framework. Softening fraction 1-1/zeta_req = (1-eps)(Pi_tr-2C_mix)/(Pi_tr-C_mix) verified.

2. **Logical flow.** Clean: regime classification → monotonicity → excess → rescue thresholds (overlap boost, stiffness softening) → twin diagnostics.

**Script Review:** Anchor values, derivative (sp.diff vs closed form), excess formula, both threshold round-trips, twin diagnostic, softening fraction cross-check. All genuine, all pass (exit code 0).

**Issues Found:** None.

---

### Agent: GPT-5 Codex — 2026-04-03
**Verdict:** PASS

**Notes Derivation Review:**

1. The regime split is algebraically exact. Rechecking the closed form `zeta_req = (Pi_tr - C_mix) / [C_mix - eps(2 C_mix - Pi_tr)]` gives the correct anchors `zeta_req(C_mix)=0` and `zeta_req(2 C_mix)=1`, and the derivative with respect to `Pi_tr` is manifestly positive on the stated support-needed branch.
2. The excess beyond the symmetric twin and the rescue thresholds are also correct. `Delta_zeta = zeta_req - 1` reduces to the stated rational form, and the overlap/stiffness thresholds for the non-twin lowest lane follow directly.
3. The symmetric twin diagnostics are physically sensible: `Omega_0 = 1` and `K_(phi,0)^(eff) = K_W^(eff)` indeed give `zeta_0^(twin)=1`, so any `zeta_req > 1` forces a genuine non-twin deformation.
4. I did not find a missing branch or sign mismatch in the three-regime classification.

**Script Review:**

1. The audit script computes the regime formula, derivative, excess beyond twin, fixed-stiffness and fixed-overlap thresholds, and the symmetric-twin diagnostics directly from SymPy.
2. The saved output matches the note and the checks are nontrivial.

**Issues Found:**

None.

**Questions:**

None.

---
