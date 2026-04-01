# Review: Stage 067 — Full reduced pde writeup

**Batch:** 10 — Full Reduced PDE Write-Up [CP]
**Status:** Verified (2× MINOR, 2026-04-03)

**This is a CHECKPOINT stage.** Also verify cross-stage consistency (Protocol C).

## Files Under Review

- **Notes:** `notes/moving_throat/moving_throat_pde_stage067_full_reduced_pde_writeup.md`
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

**Notes Derivation Review:**

**Cross-stage consistency of formulas:** All reproduced formulas exactly match their source stages: coupled operator (Stage 041/065), zeta_phys (Stages 039-040/065), zeta_req (Stages 034-035/065), master residual R_quad (Stage 065), fixed-point law and bracket (Stage 041), Family-1 data (kappa_F1=12321/5, eta_F1=37, Xi_F1=136900 Theta_w) consistent with Stages 056-058/065-066.

**Cross-stage consistency of numerical values:** Three zeta values at lambda_mu=1 match Stage 066 digit-for-digit.

**Structural faithfulness:** Claim that support/source side is "no longer the serious unresolved issue" well-supported by Stages 061/066. Identification of outgoing quadrupole normalization as remaining gap consistent with Stages 061/065/066. No unsupported claims. No silently papered-over issues.

**Issues Found:**

1. **(MINOR)** Omega_Pe used but not defined in Section 1.3. The explicit formula from Stage 039/065 should be included or cited.

2. **(MINOR)** Conservative-floor (J) branch values mentioned in prose but only shell-weighted (chi) values actually presented. Either narrow prose or include (J) window.

3. **(MINOR)** No specific stage citations — only "Stages 1-66" collectively. Individual provenance for each major result would strengthen traceability.

### Agent: GPT-5 — 2026-04-03
**Verdict:** MINOR

**Notes Derivation Review:**

1. The consolidated formulas are consistent with the reviewed source stages. I spot-checked the Stage-41 fixed-point law and brackets, the Stage-65 residual definition and `Q` inversion, and the Stage-66 Family-1 values and `zeta` window quoted here; they line up.
2. The structural summary is also accurate: by this point the support/source side has been reduced to an operator-selected scalar supply, and the remaining unresolved object is the outgoing quadrupole-normalization branch that fixes `Pi_tr`, `C_mix`, and `eps_blk`.
3. As a checkpoint write-up skeleton, though, this note is still a summary layer rather than a self-contained theorem statement, so the main weaknesses are presentation and provenance rather than algebra.

**Script Review:**

No script is expected for this consolidation stage.

**Issues Found:**

1. **(MINOR)** `Omega_Pe` is used in Section 1.3 without giving its explicit formula or an immediate stage citation, even though the checkpoint is supposed to function as a write-up skeleton.
2. **(MINOR)** Section 3 says both the shell-weighted and conservative-floor Family-1 branches produce direct operator windows, but the note only prints the shell-weighted numeric window. The prose should be narrowed or the missing numbers added.
3. **(MINOR)** The checkpoint still lacks fine-grained provenance. Reproduced formulas are correct, but without per-result stage citations the review trail is weaker than it should be at a checkpoint.

---
