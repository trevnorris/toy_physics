# Review: Stage 201 — Final local mixed-ray closure theorem

**Batch:** Batch 22 — Realization Compiler
**Status:** Hardened and verified (dual-CAS PASS, 2026-04-21)

## Files Under Review

- **Notes:** `notes/moving_throat_pde_stage201_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure.md`
- **Script:** `scripts/moving_throat_pde_stage201_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
- **Mathematica:** `mathematica/moving_throat_pde_stage201_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`

## Review Checklist

- [ ] Equation-level correctness (signs, factors, indices, limits)
- [ ] Logical flow from prior stage(s)
- [ ] Assumptions stated and justified
- [ ] Notation consistent with prior stages
- [ ] Physical interpretation sensible
- [ ] SymPy script faithfully implements notes
- [ ] Mathematica script faithfully implements notes
- [ ] Scripts run without error
- [ ] Script output matches notes claims
- [ ] No missing edge cases or branches

## Agent Reviews

<!-- Agents: append your review below this line using the template: -->

### Agent: Codex GPT-5 — 2026-04-20
**Verdict:** PASS

**Notes Derivation Review:**

The Stage 201 note stays inside the right claim boundary. It does not pretend to
solve the nonlinear PDE; it closes the declared local mixed-ray search class by
splicing the carried Stage-198 support-`<=4` ledger to the carried Stage-200
unique support-5 interior packet. The boundary-identification theorem, the
support-cardinality ceiling, and the final splice theorem are all stated as
ledger-closure results rather than as realized-branch existence claims.

**Script Review:**

The existing SymPy audit checks the right closure structure:

1. the unique five-simplex face stratification and boundary coverage counts,
2. the support ceiling at five primitive free axes,
3. the splice formulas for the imported support-`<=4` and support-5 intervals,
4. exhaustive integer-sample verification of the splice / improvement /
   no-improvement theorems,
5. the carried Stage-198 and Stage-200 budget arithmetic.

I added a Mathematica mirror at
`mathematica/moving_throat_pde_stage201_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
that checks the same combinatorial, splice, and budget structure in the second
CAS.

**Issues Found:**

None.

### Agent: Codex GPT-5 — 2026-04-21
**Verdict:** PASS

**Notes Derivation Review:**

The Stage `201` claim is now checked against the actual imported ledger objects
rather than against generic interval placeholders. The audit still treats the
result as a symbolic/combinatorial closure theorem, but it now reconstructs the
free-quintuple boundary from the five imported Stage-198 face packets, rebuilds
the support-`<=4` and support-`<=5` ledgers from those packets plus the unique
Stage-200 interior packet, and checks the three theorem classes on that finite
packet family.

**Script Review:**

The hardened SymPy and Mathematica audits now verify:

1. the exact Stage-198 boundary packet structure:
   one support-`<=3` packet plus five primitive quadruple face packets on the
   unique quintuple;
2. the exact boundary-identification incidence counts for all 30 proper faces;
3. the imported support-`<=4` ledger
   `min(tau_le3, tau_face_lambda, ..., tau_face_W)`;
4. the imported support-`<=5` splice
   `min(boundary ledger, support-5 interior packet)` together with the flattened
   packet form;
5. three exact finite family checks on the actual packet map:
   strict support-5 improvement, strict no-improvement, and finite overlap;
6. the Stage-200 support-five candidate filters via the real compiler degree
   patterns `(3,3,3,3,2)` and `(5,5,5,6)` with canonical screens
   `{gradient-optimal, equal-mix}`;
7. the final budgets rebuilt from the imported packet counts rather than from
   naked totals.

This closes the main weakness in the previous version: the old “exhaustive
splice theorem” was mostly a universal `min` lemma, while the new audit is tied
to the actual Stage-198/200 ledger topology.

**Issues Found:**

None in the hardened Stage `201` path.

---
