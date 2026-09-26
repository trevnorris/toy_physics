PLAN: NEEDS REVISION

The route is physically sound and fits comfortably in one guarded job. It needs two small, specific corrections before production. Otherwise its most likely result, a zero, would come with no evidence that the final assembly step could have shown anything else, and one reuse premise is misdescribed.

## My own reading of the physics

- **The centre is a separate degree of freedom.** a `dof_fields` and `build_face_source` give face displacement ζ_s = ζ_c + s·W0·e_W/2. In the graph representation, face work is −Σ_s s·P_s·δζ_s. Here P_s is the physical-face normal load −(δp_s + Λ_X𝒜_s), evaluated at the shifted face, so it includes the shift term.
  - The thickness load comes out as a sum over faces, Σ P_s.
  - The centre load comes out as an odd combination, Σ s·P_s.
  - The supplied U and T carry no ζ_c. The material constraint Σ = ρ_4D·W does not involve ζ_c either. So the only centre force is face work. That matches b `face_generalized_force_rows`, where the centre row comes only from `virtual_work_shape_deriv` along the (DELTA_W, ZETA_C) axis.
- **Normalization.** b multiplies only the thickness and in-plane face rows by `ACTION_TO_STORED_ROW_MULTIPLIER` (b lines 634–637). The centre row is left in action orientation.
- **The saved literals agree with this.**
  - The saved centre row (`saved-operands.json`, hash `bc1592…`, the same in all four cases) contains only the four face slots. It has no μ_θ, field or σ_W terms.
  - The normalized E_W face terms (`normalized-face-literal-terms.json`) have the same within-face structure.
  - Reading them by eye, not with CAS: for each face, the centre coefficient divided by the normalized E_W coefficient looks like −2s/W0 for both the pressure slot and the jet slot. That is a constant. So the plan's step-4 requirement of a field-free scalar will probably be met.
- **Sum/difference bookkeeping.** In c2 `main`, FACE_DIFFERENCE is defined as (f₊ − f₋)/2. With an odd multiplier, the closed centre load becomes a multiple of FACE_DIFFERENCE(E_W).
- **Which sector this tests.** c2 `build_face` drives the closure only with face velocity along the DELTA_W direction. So the diagnostic tests the centre source coming from the ζ_c = 0 sector. That is the correct necessary condition to check.

## Substantive findings

**1. Blocking: a zero result would not be tested by the controls as designed.**
- In all four cases the saved `FACE_DIFFERENCE/E_W` is the literal `Integer(0)` (`reuse-inputs.json` lines 326, 565, 804, 1043).
- The plan's step 5 builds the result from the saved parity blocks. With a face-odd multiplier, the output is therefore the multiplier times a saved zero.
- Both proposed controls (plan lines 106–113) break the relation inside one face:
  - reversing one pressure-slot term makes that face's pressure ratio and jet ratio disagree;
  - swapping one face's jet slot with the other face's removes one face's jet term.
- Both therefore stop at step 3, "no safe relation". Neither ever exercises the step-5 assembly. Amendment §3.3 does not accept an insensitive zero comparison.
- **Required correction:**
  - Add one control that keeps the within-face relation valid but flips the parity across faces. For example, reverse both slots of one face together in the imported centre row. The multipliers then come out equal on the two faces, and step 5 must return a load proportional to the saved FACE_SUM, which is non-zero. One case is enough.
  - State in the report that any zero centre load is derived from c2's saved FACE_DIFFERENCE. It carries c2's debts, including the face-force sign convention that c2 §1a says has not been checked across engines. It is not an independent recomputation of the closed centre force.

**2. Mandatory: the reuse premise names the wrong post-closure operations.**
- The plan says the per-face and parity returns come "after the native weak restriction" (line 52) and that the multiplier must pass through the "weak restriction" (line 87).
- In c2 `build_case` (lines 294–305) and `main` (lines 445–451), `CLOSED_SLAB_OPERATOR_TERM_ORIGINS` and `PARITY_BLOCKS` are the full closed rows after three steps: substitution, `retained_shape` (truncation to η ≤ 1, σ ≤ 1) and `physical_fields`. The weak restriction `extract` is applied only to `CLOSED_COUPLING_KERNEL*`.
- **Required correction:**
  - List the actual pass-through steps: closure substitution, profile `xreplace`, `retained_shape` and `physical_fields`. A constant multiplier commutes with all of them.
  - Describe the centre load as an unrestricted row object.
  - Make the step-4 join literal. c2 imported b's `slab_operator` E_W EXPANDED row. In b `build_operator` (lines 682–695) that row equals the reduced row plus the kinetic term plus `face_e`. So:
    - add that row, or its export hash, as an input;
    - lexically confirm that its slot-bearing terms equal the slot terms of `FACE_VIRTUAL_WORK/ROWS/E_W`;
    - confirm that the b `.out` record is the same object as the `af560257` export c2 consumed.

## Nonblocking observations

- **Units.** The saved FACE_DIFFERENCE carries dimension `[0,0,0]` only because its value is zero. Take units from the multiplier (per length) times the face-row units `[-1,-2,1]`, and check the result against the centre row's `[-2,-2,1]`.
- **Orientation.** Record the action-to-stored orientation factor separately from the geometric factor 2s/W0, using the saved `ROW_NORMALIZATION`.
- **Terms outside the slots.** The instruction to account for any part outside the slots (step 2) will probably find nothing; report that explicitly.
- **Advection wording.** The shared background-advection term is already inside c2's V_s(DELTA_W). It does not enter the centre/E_W slot relation. In the a source, the Eulerian-route LAB_HELD branch carries u·∇W_bg in its virtual displacement and velocity, while MATERIAL carries it in the height perturbation dh. The phrase "material direction records" is therefore imprecise.
- **Reduction claim.** If the centre load is identically zero on arbitrary trial fields, it is also zero on the solved sector, so (x_e, ζ_c = 0) is *a* solution. Uniqueness then depends on:
  - the centre's own operator: the odd-parity face response to ζ_c-driven motion, with no centre inertia or stiffness in the supplied T and U;
  - whether that operator is nondegenerate on the domain.
  
  Coupling from ζ_c back into the five fields is not needed for consistency. Name this missing premise concretely. Reflection symmetry w → −w explains why a zero is expected, but c1 §2a/§3a says parity is a computed result, so the plan is right not to assume it.

## Source coverage

I read INDEX, the plan and the accepted scope in full, plus:
- b physics §§1–3 and b source excerpts (`face_generalized_force_rows`, `mechanical_work_row_normalization`, `kinetic_balance_from_energy`, `constraint_fold_from_source`, `build_operator`);
- c2 physics §§0–3 and c2 source (`build_face`, `build_case`, `retained_shape`, `main`);
- the relevant c1 sections §§1c–3a and a source (`dof_fields`, `build_material_face_source`, `build_face_source`, `supplied_face_maps`);
- saved-operands (centre row and KINETIC literals), normalized-face-literal-terms, reuse-inputs (all eight parity fields) and the start of closure-inputs.

I did not open navigation-manifest.json or the full `.out` values. I ran no CAS and did no numerical check. The −2s/W0 ratio is by eye only.

## Adequacy of joins and stopping conditions

- The b/c2 slot joins are adequate once finding 2 is fixed.
- The stopping rules are honest and allow a legitimate stop.
- The resource cap is ample: the new algebra works on roughly 2 KB and 200 KB expressions, and the added control needs one FACE_SUM record of about 1.2 MB.

## Minimum safe next action

The author should make the two text edits (the parity-flip control plus the statement of where the zero comes from, and the corrected pass-through list plus the literal join). The coordinator can then clear the plan without another full review cycle, and the single guarded job can run.
