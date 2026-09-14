(* Source-pinned pressure-slot audit on the constant-height, diagonal-momentum
   locus. An exact counterexample on this locus can establish a defect, but a
   zero result cannot establish the full ordered eta/sigma rectangle. *)
$HistoryLength = 0;
ClearAll["Global`*"];
$Messages = {OutputStream["stderr", 2]};
auditRoot = Environment["S11C_WOLFRAM_AUDIT_ROOT"];
auditHeldC2 = ToExpression[Import[FileNameJoin[{auditRoot, "sources", "mathematica", "S11c_c2_N6_mathematica_audit.wl"}], "Text"], InputForm, HoldComplete];
auditHeldC1 = ToExpression[Import[FileNameJoin[{auditRoot, "sources", "mathematica", "S11c_c1_bulk_closure_mathematica_audit.wl"}], "Text"], InputForm, HoldComplete];
auditRecords = <||>;
auditEmit[key_String, value_, unit_, epsilon_:0] := Module[{record},
  record = <|"value" -> value, "dimensionLTM" -> unit, "epsilonOrder" -> epsilon,
    "etaSigmaDomain" -> "eta retained through 1; sigma=0 exact constant-height locus",
    "lambda" -> "eta homotopy"|>;
  WriteString[$Output, key, " = ", ToString[record, InputForm, PageWidth -> Infinity], "\n"];
  If[KeyExistsQ[auditRecords, key], Quit[91]]; AssociateTo[auditRecords, key -> record]];
auditTake[held_, pattern_, label_, count_:1, levels_:{2}] := Module[{found},
  found = Cases[held, x:pattern :> HoldComplete[x], levels];
  auditEmit["wlTraceSource" <> label, <|"count" -> Length[found], "sha256" -> Hash[found, "SHA256", "HexString"], "definitions" -> found|>, {0, 0, 0}];
  If[count =!= Automatic && Length[found] =!= count, Quit[92]];
  If[found === {}, Quit[93]]; found];
auditLoad[held_, pattern_, label_, count_:1] := Scan[ReleaseHold, auditTake[held, pattern, label, count]];
auditTop = Cases[auditHeldC2, x_ :> HoldComplete[x], {2}];
auditFirst = First[FirstPosition[auditTop, HoldPattern[HoldComplete[Set[jetRegistry, _]]], Missing["Absent"], {1}]];
auditLast = First[FirstPosition[auditTop, HoldPattern[HoldComplete[SetDelayed[el[___], _]]], Missing["Absent"], {1}]];
auditSetup = Take[auditTop, {auditFirst, auditLast}];
auditEmit["wlTraceSourceSetup", <|"range" -> {auditFirst, auditLast}, "sha256" -> Hash[auditSetup, "SHA256", "HexString"], "definitions" -> auditSetup|>, {0, 0, 0}];
If[!IntegerQ[auditFirst] || !IntegerQ[auditLast], Quit[94]];
Scan[ReleaseHold, auditSetup];
auditLoad[auditHeldC2, HoldPattern[SetDelayed[relationalObject[___], _]], "Relation"];
auditLoad[auditHeldC2, HoldPattern[Set[momentum | qLeg | virtualU, _]], "LegsVirtual", 3];
auditLoad[auditHeldC2, HoldPattern[SetDelayed[(linearShape | shapeBackground | firstBackgroundJet | geometryExpansion | graphGeometry | faceLaws | sourceConstruction | constructKernel | normalContinuation | eulerianSlabFace | materialFaceFold)[___], _]], "GeometryAndResponse", 12];
auditLoad[auditHeldC2, HoldPattern[SetDelayed[(add | neg | sub | mul | power | circuitExpression | gAdd | gSub | gMul | gScale | graded | responseFamilies)[___], _]], "GradedResponse", 21];
auditLoad[auditHeldC1, HoldPattern[SetDelayed[flatFaceResponseSolve[___], _]], "C1FlatSolve"];
auditCaseDefinition = First[auditTake[auditHeldC2, HoldPattern[SetDelayed[buildCase[___], _]], "BuildCase"]];
auditImageAssignments = auditTake[auditCaseDefinition,
  HoldPattern[Set[pressureImage, _] | AssociateTo[imageRules, _]], "NativePressureImages", 3, Infinity];
auditDensityAssignments = auditTake[auditCaseDefinition, HoldPattern[Set[density4 | density3, _]], "NativeDensity", 2, Infinity];
auditFaceDefinition = First[auditTake[auditHeldC2, HoldPattern[SetDelayed[faceLaws[___], _]], "FaceLaw"]];
auditTraceAssignments = auditTake[auditFaceDefinition, HoldPattern[Set[p | trace, _]], "NativeTrace", 2, Infinity];
(* Explicit coefficient extraction also discards a rational residual whose first
   nonzero power lies beyond the requested order. Keep that raw remainder. *)
auditFirstOrder[value_List] := auditFirstOrder /@ value;
auditFirstOrder[value_] := Expand[Sum[SeriesCoefficient[value, {etaBg, 0, n}] etaBg^n, {n, 0, 1}]];
auditPhysical[key_, value_, unit_, epsilon_:1] := Module[{raw, retained},
  raw = Together[value]; retained = auditFirstOrder[raw];
  auditEmit[key, <|"raw" -> raw, "retained" -> retained,
    "discarded" -> Together[raw - retained],
    "etaSigmaSupport" -> Select[{{0, 0}, {1, 0}},
      AnyTrue[Flatten[{Coefficient[retained, etaBg, #[[1]]]}], !TrueQ[# === 0] &] &]|>, unit, epsilon]];
auditDecode[value_] := Activate[circuitExpression[value]];
auditDiagonalRules = Join[Thread[momentum["kOut"] -> momentum["kIn"]],
  Thread[momentum["kMiddle"] -> momentum["kIn"]], {qLeg["Out"] -> qLeg["In"], qLeg["Middle"] -> qLeg["In"]}];
auditEmit["wlTraceDomain", <|"anchorings" -> {"LAB_HELD", "MATERIAL_ADVECTED"},
  "densities" -> {"RHO4_CONSTANT", "RHOBR_CONSTANT"}, "faces" -> faces,
  "locus" -> "stationary constant nonzero height contrast; kOut=kMiddle=kIn; all width jets zero",
  "exclusions" -> {"W0=0", "rhoBr=0", "rhoM=0", "memory denominators=0", "response matching denominators=0"},
  "unexecuted" -> "nonconstant-profile ordered two/three-leg eta/sigma trace inverse; no global or spectral clearance"|>, {0, 0, 0}];
auditSourceSolve = sourceConstruction[];
auditCaseIndex = 0;
Do[
  auditCaseIndex++;
  auditNativeDensity = Block[{densityKind = auditDensityKind, density4, density3},
    Scan[ReleaseHold, auditDensityAssignments]; density3] /. WBg -> W0 (1 + etaBg auditContrast);
  auditKernel = constructKernel[auditSign];
  auditResponse = responseFamilies[auditKernel];
  auditPhysical["wlTraceDiagonalFirstKernel" <> ToString[auditCaseIndex],
    Together[(auditDecode /@ auditResponse["FOURIER_KOUT_KIN_Y"]) /. auditDiagonalRules], {0, -1, 1}, 0];
  auditPhysical["wlTraceDiagonalMiddleKernel" <> ToString[auditCaseIndex],
    Together[(auditDecode /@ auditResponse["FOURIER_KOUT_KIN_Y_MIDDLE"]) /. auditDiagonalRules], {3, -1, 1}, 0];
  auditDrive = auditSourceSolve["SOURCE"] /. {rhoFace -> auditNativeDensity};
  auditImages = Block[{face = auditSign, response, sourceForGuard = retain[auditDrive],
      familyName = "FOURIER_KOUT_Y", pressureImage, imageRules = <||>},
    response = <|auditSign -> auditResponse|>;
    Scan[ReleaseHold, auditImageAssignments]; imageRules];
  auditPressureSlot = pressureSlots[[If[auditSign === 1, 1, 3]]];
  auditJetSlot = pressureSlots[[If[auditSign === 1, 2, 4]]];
  auditImageExpressions = Map[Total[MapThread[Times, {auditDecode /@ #, (etaBg^#[[1]] sigmaW^#[[2]] & /@ gradeIndices)}]] &, auditImages];
  auditImageExpressions = auditImageExpressions /. auditDiagonalRules;
  auditTraceExpression = Block[{s = auditSign, p, trace}, Scan[ReleaseHold, auditTraceAssignments]; trace] /.
    {shapeParameter -> 1, WBg -> W0 (1 + etaBg auditContrast)};
  auditNativeFold = auditTraceExpression /. Normal[auditImageExpressions];

  (* An independent exact wave and boundary closure. h is an input face map;
     all pressure/reference values below are obtained by differentiation/Solve. *)
  auditHeight = etaBg W0 auditContrast/2;
  auditWave = auditAmplitude Exp[I (momentum["kIn"].{auditX1, auditX2, auditX3} + auditSign qLeg["In"] auditW - omega auditTime)];
  auditPhase = Exp[I (momentum["kIn"].{auditX1, auditX2, auditX3} - omega auditTime)];
  auditPressureField = -rhoM D[auditWave, auditTime]/auditPhase;
  auditVelocityField = auditSign D[auditWave, auditW]/auditPhase;
  auditPressureAtFace = Simplify[auditPressureField /. auditW -> auditSign auditHeight];
  auditVelocityAtFace = Simplify[auditVelocityField /. auditW -> auditSign auditHeight];
  auditFluxEquation = auditFlux == lambdaA (muTheta/auditNativeDensity - auditPressureAtFace/rhoM) + lambdaV vFace;
  auditVelocityEquation = auditVelocityAtFace == vFace + auditFlux/rhoM;
  auditSolution = First[Solve[{auditFluxEquation, auditVelocityEquation}, {auditAmplitude, auditFlux}]];
  auditExactFacePressure = Simplify[auditPressureAtFace /. auditSolution];
  auditExactReferencePressure = Simplify[(auditPressureField /. auditW -> 0) /. auditSolution];
  auditReferenceJet = Simplify[(D[auditPressureField, auditW] /. auditW -> 0) /. auditSolution];
  auditIndependentFold = auditTraceExpression /. {auditPressureSlot -> auditExactReferencePressure, auditJetSlot -> auditReferenceJet};
  auditC1 = flatFaceResponseSolve[auditKernel["FLAT"]["In"], auditNativeDensity];
  auditC1Pressure = auditC1["PRESSURE"] /. {epsilonShape -> 1, faceVelocityInput -> vFace};
  auditPhysical["wlTraceC1PressureJoin" <> ToString[auditCaseIndex],
    {auditC1Pressure, auditExactFacePressure, auditC1Pressure - auditExactFacePressure}, {-2, -2, 1}];
  auditPhysical["wlTraceNativeResponseJoin" <> ToString[auditCaseIndex],
    {auditImageExpressions[auditPressureSlot], auditExactFacePressure,
     auditImageExpressions[auditPressureSlot] - auditExactFacePressure}, {-2, -2, 1}];
  auditPhysical["wlTraceReferencePressure" <> ToString[auditCaseIndex], auditExactReferencePressure, {-2, -2, 1}];
  auditPhysical["wlTraceReferenceNormalJet" <> ToString[auditCaseIndex], auditReferenceJet, {-3, -2, 1}];
  auditPhysical["wlTraceIndependentAffineResidual" <> ToString[auditCaseIndex], auditIndependentFold - auditExactFacePressure, {-2, -2, 1}];
  auditPhysical["wlTraceNativeFoldResidual" <> ToString[auditCaseIndex], auditNativeFold - auditExactFacePressure, {-2, -2, 1}];
  (* A concrete dispersion-compatible lossless witness is evaluated only after
     the symbolic operands. It has nonzero height, q, frequency and source. *)
  auditWitness = {W0 -> 2, auditContrast -> 1, rhoBr -> 3, rhoM -> 5,
    omega -> 13, qLeg["In"] -> 5, cS0 -> 1,
    muTheta -> 7, vFace -> 11, LambdaA0 -> 0, LambdaV0 -> 0, LambdaX0 -> 0,
    tauA -> 0, tauV -> 0, tauX -> 0, etaBg -> 1/100};
  auditWitness = Join[auditWitness, Thread[momentum["kIn"] -> {12, 0, 0}]];
  auditEmit["wlTraceWitness" <> ToString[auditCaseIndex], <|
    "case" -> {auditAnchor, auditDensityKind, auditSign},
    "bindingsInRestoredUnitFrame" -> auditWitness,
    "dispersionResidual" -> Together[auditKernel["BULK_EQUATION_OPERAND"] /. auditWitness],
    "matchingDenominators" -> Together[Values[auditResponse["DENOMINATORS"]] /. auditDiagonalRules /. auditWitness],
    "nativeFoldResidual" -> Together[auditFirstOrder[auditNativeFold - auditExactFacePressure] /. auditWitness]|>,
    <|"dispersionResidual" -> {0, -2, 0}, "matchingDenominators" -> {-1, 0, 0}, "nativeFoldResidual" -> {-2, -2, 1}|>, 1],
  {auditAnchor, {"LAB_HELD", "MATERIAL_ADVECTED"}},
  {auditDensityKind, {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}, {auditSign, faces}];
auditEmit["wlTraceCompletion", <|"cases" -> auditCaseIndex, "emittedRecords" -> Length[auditRecords]|>, {0, 0, 0}];
auditSummary = <|"cases" -> auditCaseIndex, "records" -> Length[auditRecords],
  "c1JoinResiduals" -> Table[ToString[auditRecords["wlTraceC1PressureJoin" <> ToString[i]]["value"]["retained"][[3]], InputForm], {i, auditCaseIndex}],
  "nativeResponseJoinResiduals" -> Table[ToString[auditRecords["wlTraceNativeResponseJoin" <> ToString[i]]["value"]["retained"][[3]], InputForm], {i, auditCaseIndex}],
  "independentAffineResiduals" -> Table[ToString[auditRecords["wlTraceIndependentAffineResidual" <> ToString[i]]["value"]["retained"], InputForm], {i, auditCaseIndex}],
  "nativeFoldResiduals" -> Table[ToString[auditRecords["wlTraceNativeFoldResidual" <> ToString[i]]["value"]["retained"], InputForm], {i, auditCaseIndex}],
  "witnessResiduals" -> Table[ToString[auditRecords["wlTraceWitness" <> ToString[i]]["value"]["nativeFoldResidual"], InputForm], {i, auditCaseIndex}]|>;
Export[FileNameJoin[{auditRoot, "summary.json"}], auditSummary, "RawJSON"];
