# Measurements behind S11b_wl_engine_repair_directive.md (rule 2)
# The F-WL-1/2/3 defects are measured in steps/S11b_wl_engine_review_disposition.md (committed) and
# its _measurements twin; the repair directive targets the same, at these engine lines:
```
367:equationZeroForm[equal_Equal] := equal[[1]] - equal[[2]];
368:equationZeroForm[expression_] := expression;
374:  zeroForms = Together /@ (equationZeroForm /@ equations);
565:faceEquationZeroForms = equationZeroForm /@ generalFace["EQUATIONS"];
813:  (equationZeroForm /@ compressionEliminationEquations);
859:  equationZeroForm[rawPressureSliceClosure[[3]]], facePressure];
865:emitShared["ZPERM_SLICE_MAP", <|
939:allFluxCleared = Factor[equationZeroForm[allFluxCrossEquation]];
940:allForceCleared = Factor[equationZeroForm[allForceCrossEquation]];
1102:emitShared["KERNEL_ORIENTATION_IDENTITIES", kernelOrientationIdentities];
```
F-WL-1: rawPressureCoefficient = Coefficient[equationZeroForm[flux], facePressure] (zero-form negates
  the flux coefficient; equationZeroForm[a==b]=a-b at L367) + emits the dynamic Lambda_A(omega) form.
F-WL-2: PRESSURE_WORK_SIGN_CHECK residual = -1/2Re[X]+1/2Re[X] identically 0; energy/two-port checks
  consume hand-typed energyEquationRules, not fullSystem (byte-identical under a traction FORM ablation).
F-WL-3a: KERNEL_ORIENTATION_IDENTITIES from closureLawForExtraction shadow; 3b CAUSALITY_CHECK Cases->{}
  And@@{}=True; 3c GRAZING Limit[q*fullSystem[MATRIX],q->0]=zero matrix rank 0.
