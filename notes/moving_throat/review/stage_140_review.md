# Review: Stage 140 — Core mouth coevolution status

**Batch:** 17 — Rigidity & Corrections
**Status:** Verified (1× PASS, 1× MINOR, 2026-04-03)

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage140_core_mouth_coevolution_status.md`
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

<!-- Agents: append your review below this line using the template:

### Agent: [Model Name] — [Date]
**Verdict:** [PASS | MINOR | ISSUE | BLOCK]

**Notes Derivation Review:**

**Script Review:**

**Issues Found:**

**Questions:**

### Agent: Claude Opus 4.6 — 2026-04-02
**Verdict:** MINOR
**Notes Derivation Review:** Qualitative content correct: consolidates stages 137-139. Branch structure claims (upper impossible, equal-normalized singular, lower uniquely regular, renormalized fixed point) all traceable. Formulas reproduced correctly. Forward framing appropriate.
**Issues Found:**
1. **(MINOR)** Frozen-traction fixed-point values (g_fp, R_fp, Pi_fp) differ from Stage 138 at 5th-7th significant figure. Status stage should quote prior values exactly.
2. **(MINOR)** Renormalized canonical branch values (Sigma_0^can, T_hat, Pi_can) differ from Stage 139 at 5th significant figure. Possibly from different iteration parameters — should be reconciled or noted.
3. **(MINOR)** No explicit stage citations for quoted numerical values.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The consolidation narrative is correct: the stage properly states that full co-evolution yields a renormalized finite branch rather than preserving the original canonical point.
2. The quoted values are internally consistent as a standalone summary, but they do not exactly match the more precise values reported in Stages 138-139, and the note does not explain that discrepancy.
3. The open question is still framed properly. The remaining issue is numerical consolidation and citation, not a broken derivation.

**Script Review:**

No script is attached for this consolidation note, which is acceptable. The issue is not missing computation in this file but the fact that its quoted numbers shift slightly relative to the earlier stage outputs without explicit explanation.

**Issues Found:**

1. **MINOR:** The summary numbers for the frozen and renormalized branches differ slightly from Stages 138-139 and are not cited or reconciled.

**Questions:** None.

---
