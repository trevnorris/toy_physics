# Measurements twin — S11c_a_wl_engine_repair_directive.md

## Legs (both quoted verbatim; orchestrator verified each at source)
- Agent leg af9664b5ad52f8672: scripts /tmp/s11ca_wl/{derive_sympy,probe*,parse_test2,parse_scenarios}.* — found Defect 1 (T-h parse), proved from exact bytes.
- Grok leg bf6eufsrh: log ~/.s11_build/S11c_a_wl_engine_grok2.log; scripts /tmp/s11ca_wl/{derive_s11ca,probe_routes,probe_dof_and_routes,engine_stripped,engine_form_ablate}.* — found Defects 2,3,4,5.

## Source verifications (CMD: sed -n on the committed baseline 277f3fe7)
- 481-484 sigmaEulerianSource: multi-line implicit product; first line complete ⇒ newline drops (1+θ)+thickness. T-0 645 Together[rho4·widthAnsatz] correct. ⇒ Defect 1 CONFIRMED.
- 346 kinematic = n·v_bulk - V - flux/rhoM (same Module) ⇒ ≡0; emit 797-808 TEST_OBJECT relationalObject[shapeDerivative[exact],0]. ⇒ Defect 2 CONFIRMED (bare Equal[0,0]).
- 204-213 flatteningCoordinateSource uses raw zetaCenter/deltaWidth; 129-131 zetaSource uses dofCenter/dofWidth; 179-188 heightSource→zetaSource. ⇒ route-2 not DOF-keyed. Defect 3 CONFIRMED.
- 490-524 virtualConstraintDirectSource vs 526-549 virtualConstraintMaterialSource: both rhoFour*(1+θ)*W*jacobian, thickness==denominator, no Solve[w']/graphThicknessSource ⇒ byte-identical. Defect 4 CONFIRMED.
- 1106-1109 CONTROL_INDEPENDENCE: corrupted=applyProfileJets[applyPhysicalDof[directShape[...],dof],0,sign==1]; 419 directShape=shapeDerivative[...] (RESULT); 164-165 applyProfileJets reversal flips w1Jet[1,0,0]. ⇒ mutates a RESULT. Defect 5 CONFIRMED.

## Spec authority for the fixes (route-2 = material w' flattening; one-sided at the SOURCE)
CMD: sed -n '452,481p' research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md

## BLINDNESS: the directive quotes ONLY the WL engine's own defects (legs derived from the SPEC) + the spec §5a. NO sibling (SymPy) construction, value, or result referenced. Blind-safe.
