(* Source-pinned repaired pressure-slot audit. Generic ordered three-leg matrix
   solves retain every eta/sigma coefficient; a separate exact shifted-wave
   regression and native face-row coefficient checks cover the consumer boundary. *)
$HistoryLength = 0;
ClearAll["Global`*"];
$Messages = {OutputStream["stderr", 2]};
auditRoot = Environment["S11C_WOLFRAM_AUDIT_ROOT"];
auditHeldC2 = ToExpression[Import[FileNameJoin[{auditRoot, "sources", "mathematica", "S11c_c2_N6_mathematica_audit.wl"}], "Text"], InputForm, HoldComplete];
auditHeldC1 = ToExpression[Import[FileNameJoin[{auditRoot, "sources", "mathematica", "S11c_c1_bulk_closure_mathematica_audit.wl"}], "Text"], InputForm, HoldComplete];
auditRecords = <||>;
auditEmit[key_String, value_, unit_, epsilon_:0] := Module[{record},
  record = <|"value" -> value, "dimensionLTM" -> unit, "epsilonOrder" -> epsilon,
    "etaSigmaDomain" -> "independent eta/sigma rectangle through (1,1); lambda=eta physical homotopy",
    "lambda" -> "eta homotopy"|>;
  WriteString[$Output, key, " = ", ToString[record, InputForm, PageWidth -> Infinity], "\n"];
  If[KeyExistsQ[auditRecords, key], Quit[91]]; AssociateTo[auditRecords, key -> record]];
auditTake[held_, pattern_, label_, count_:1, levels_:{2}] := Module[{found},
  found = Cases[held, x:pattern :> HoldComplete[x], levels];
  auditEmit["wlRepairedTraceSource" <> label, <|"count" -> Length[found], "sha256" -> Hash[found, "SHA256", "HexString"], "definitions" -> found|>, {0, 0, 0}];
  If[count =!= Automatic && Length[found] =!= count, Quit[92]];
  If[found === {}, Quit[93]]; found];
auditLoad[held_, pattern_, label_, count_:1] := Scan[ReleaseHold, auditTake[held, pattern, label, count]];
auditTop = Cases[auditHeldC2, x_ :> HoldComplete[x], {2}];
auditFirst = First[FirstPosition[auditTop, HoldPattern[HoldComplete[Set[jetRegistry, _]]], Missing["Absent"], {1}]];
auditLast = First[FirstPosition[auditTop, HoldPattern[HoldComplete[SetDelayed[el[___], _]]], Missing["Absent"], {1}]];
auditSetup = Take[auditTop, {auditFirst, auditLast}];
auditEmit["wlRepairedTraceSourceSetup", <|"range" -> {auditFirst, auditLast}, "sha256" -> Hash[auditSetup, "SHA256", "HexString"], "definitions" -> auditSetup|>, {0, 0, 0}];
If[!IntegerQ[auditFirst] || !IntegerQ[auditLast], Quit[94]];
Scan[ReleaseHold, auditSetup];
auditLoad[auditHeldC2, HoldPattern[SetDelayed[relationalObject[___], _]], "Relation"];
auditLoad[auditHeldC2, HoldPattern[Set[momentum | qLeg | virtualU, _]], "LegsVirtual", 3];
auditLoad[auditHeldC2, HoldPattern[SetDelayed[(linearShape | shapeBackground | firstBackgroundJet | geometryExpansion | graphGeometry | faceLaws | sourceConstruction | constructKernel | normalContinuation | eulerianSlabFace | materialFaceFold | facePressureTrace | materialGeometry | carrier)[___], _]], "GeometryAndResponse", 15];
auditLoad[auditHeldC2, HoldPattern[SetDelayed[(add | neg | sub | mul | power | circuitExpression | gAdd | gSub | gMul | gScale | graded | responseFamilies | referencePressureFamilies | gInverse)[___], _]], "GradedResponse", 23];
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
auditFirstOrder[value_] := Module[{rational = Together[value], numerator, denominator},
  numerator = Numerator[rational]; denominator = Denominator[rational];
  If[FreeQ[denominator, etaBg | sigmaW],
    Total[Table[etaBg^g[[1]] sigmaW^g[[2]] Coefficient[Coefficient[
      Expand[numerator], etaBg, g[[1]]], sigmaW, g[[2]]]/denominator, {g, gradeIndices}]],
    Total[Table[etaBg^g[[1]] sigmaW^g[[2]] SeriesCoefficient[SeriesCoefficient[
      rational, {etaBg, 0, g[[1]]}], {sigmaW, 0, g[[2]]}], {g, gradeIndices}]]]];
auditPhysical[key_, value_, unit_, epsilon_:1] := Module[{raw, retained},
  raw = Together[value]; retained = Together[auditFirstOrder[raw]];
  auditEmit[key, <|"raw" -> raw, "retained" -> retained,
    "discarded" -> Together[raw - retained],
    "etaSigmaSupport" -> Select[gradeIndices,
      AnyTrue[Flatten[{Coefficient[Coefficient[retained, etaBg, #[[1]]], sigmaW, #[[2]]]}], !TrueQ[# === 0] &] &]|>, unit, epsilon]];
auditDecode[value_] := Activate[circuitExpression[value]];
auditGradeExpression[values_List] := Total[MapThread[Times,
  {auditDecode /@ values, (etaBg^#[[1]] sigmaW^#[[2]] & /@ gradeIndices)}]];
auditFamilyExpression[response_, family_] := auditGradeExpression[response[family]];
auditFamilies = {"FOURIER_KOUT_Y", "FOURIER_KOUT_KIN_Y", "FOURIER_KOUT_KIN_Y_MIDDLE"};
auditCases = If[Environment["S11C_WOLFRAM_REPAIR_CASE"] === "LAB_HELD_RHO4_CONSTANT",
  {{"LAB_HELD", "RHO4_CONSTANT"}}, Tuples[{{"LAB_HELD", "MATERIAL_ADVECTED"}, {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}]];
auditLegNames = {"Out", "Middle", "In"};
auditMeasures = {auditMeasureOut, auditMeasureMiddle, auditMeasureIn};
Scan[declare[#, {-3, 0, 0}] &, auditMeasures];
auditLegRule[i_, j_, profile_] := Join[
  Thread[momentum["kOut"] -> momentum["k" <> auditLegNames[[i]]]],
  Thread[momentum["kIn"] -> momentum["k" <> auditLegNames[[j]]]],
  {qLeg["Out"] -> qLeg[auditLegNames[[i]]], qLeg["In"] -> qLeg[auditLegNames[[j]]], profileTransform -> profile}];
auditEdgeProfile[i_, j_] := Switch[{i, j}, {1, 2}, profileTransformLeft,
  {2, 3}, profileTransformRight, {1, 3}, profileTransform];
auditResponseMatrix[response_] := Table[
  Which[i === j, auditFamilyExpression[response, First[auditFamilies]] /. qLeg["Out"] -> qLeg[auditLegNames[[i]]],
    i < j, (auditFamilyExpression[response, auditFamilies[[2]]] /. auditLegRule[i, j, auditEdgeProfile[i, j]]) auditMeasures[[j]] +
      If[{i, j} === {1, 3}, auditFamilyExpression[response, auditFamilies[[3]]] auditMeasures[[2]] auditMeasures[[3]], 0],
    True, 0], {i, 3}, {j, 3}];

auditEmit["wlRepairedTraceDomain", <|"orderedLegs" -> auditLegNames,
  "measures" -> "independent quadrature weights of dimension L^-3 include Fourier normalization",
  "matrixMeaning" -> "upper-triangular three-site coefficient test retaining direct and ordered intermediate paths; not a spectral discretization claim",
  "coverage" -> "both faces, generic rational coefficients, all four eta/sigma grades; enumerated anchoring/density cases in both coordinate routes",
  "rowCases" -> auditCases,
  "excluded" -> "matching denominators, rhoBr, rhoM, W0 and trace inverse denominators zero; no threshold or global-sheet clearance"|>, {0, 0, 0}];
auditSourceSolve = sourceConstruction[];
auditResidualKeys = {};
auditCheck[key_, operands_, unit_, epsilon_:0] := (
  auditPhysical[key, operands, unit, epsilon]; AppendTo[auditResidualKeys, key]);
auditMatrixCase = 0;
Do[
  auditMatrixCase++;
  auditKernel = constructKernel[auditSign];
  auditPhysicalResponse = responseFamilies[auditKernel];
  auditReference = referencePressureFamilies[auditPhysicalResponse, auditSign];
  auditP0 = DiagonalMatrix[Table[First[auditKernel["TRACE_PRESSURE"]], {i, 3}]];
  auditP = Table[Which[i === j, auditP0[[i, j]], i < j,
    (auditGradeExpression[auditPhysicalResponse["TRACE_PRESSURE_GRADES"]] /.
      auditLegRule[i, j, auditEdgeProfile[i, j]]) auditMeasures[[j]], True, 0], {i, 3}, {j, 3}];
  auditN = Table[Which[i === j, First[auditKernel["TRACE_CONORMAL"]] /. qLeg["In"] -> qLeg[auditLegNames[[i]]],
    i < j, (auditGradeExpression[auditPhysicalResponse["TRACE_CONORMAL_GRADES"]] /.
        auditLegRule[i, j, auditEdgeProfile[i, j]]) auditMeasures[[j]] +
      If[{i, j} === {1, 3}, auditGradeExpression[auditPhysicalResponse["TRACE_MIXED_CONORMAL_GRADES"]] auditMeasures[[2]] auditMeasures[[3]], 0],
    True, 0], {i, 3}, {j, 3}];
  (* Independent linear solve of the boundary-amplitude equation. No native
     reference-pressure recurrence is called to obtain either matrix solution. *)
  auditD = auditN + lambdaA auditP/rhoM^2;
  auditEmit["wlRepairedTraceMatrixStage" <> ToString[auditMatrixCase] <> "BoundarySolve", "constructed inputs", {0, 0, 0}];
  auditAmplitude = LinearSolve[auditD, IdentityMatrix[3]];
  auditPhysicalIndependent = auditFirstOrder[auditP.auditAmplitude];
  auditReferenceIndependent = auditFirstOrder[auditP0.auditAmplitude];
  auditNativePhysical = auditFirstOrder[auditResponseMatrix[auditPhysicalResponse]];
  auditNativeReference = auditFirstOrder[auditResponseMatrix[auditReference]];
  auditA = Table[Which[i === j, auditGradeExpression[auditReference["TRACE_MAP_FLAT"]],
    i < j, (auditGradeExpression[auditReference["TRACE_MAP_FIRST"]] /. auditLegRule[i, j, auditEdgeProfile[i, j]]) auditMeasures[[j]],
    True, 0], {i, 3}, {j, 3}];
  auditExactTrace = auditP.Inverse[auditP0];
  auditIndependentInverse = auditFirstOrder[LinearSolve[auditA, auditPhysicalIndependent]];
  auditEmit["wlRepairedTraceMatrixStage" <> ToString[auditMatrixCase] <> "Joins", "constructed independent solutions", {0, 0, 0}];
  auditMatrixPrefix = "wlRepairedTraceMatrix" <> ToString[auditMatrixCase];
  auditPhysical[auditMatrixPrefix <> "PressureTrace", auditP, {-4, -1, 1}, 0];
  auditPhysical[auditMatrixPrefix <> "ConormalTrace", auditN, {-1, 0, 0}, 0];
  auditPhysical[auditMatrixPrefix <> "MatchingDeterminant", Det[auditD], {-3, 0, 0}, 0];
  auditCheck[auditMatrixPrefix <> "PhysicalJoin", {auditNativePhysical, auditPhysicalIndependent,
    auditNativePhysical - auditPhysicalIndependent}, {-3, -1, 1}];
  auditCheck[auditMatrixPrefix <> "ReferenceJoin", {auditNativeReference, auditReferenceIndependent,
    auditNativeReference - auditReferenceIndependent}, {-3, -1, 1}];
  auditCheck[auditMatrixPrefix <> "TraceMapJoin", {auditA, auditExactTrace, auditA - auditExactTrace}, {0, 0, 0}];
  auditCheck[auditMatrixPrefix <> "TraceInverseJoin", {auditNativeReference, auditIndependentInverse,
    auditNativeReference - auditIndependentInverse}, {-3, -1, 1}];
  auditCheck[auditMatrixPrefix <> "Reconstruction", {auditFirstOrder[auditA.auditNativeReference], auditPhysicalIndependent,
    auditFirstOrder[auditA.auditNativeReference - auditPhysicalIndependent]}, {-3, -1, 1}];
  (* At zero frequency solve before any division by the pressure-amplitude
     factor. The pressure-coordinate trace ratio is not used on this locus. *)
  auditZeroAmplitude = LinearSolve[auditD /. omega -> 0, IdentityMatrix[3]];
  auditZeroPhysical = (auditP /. omega -> 0).auditZeroAmplitude;
  auditZeroReference = (auditP0 /. omega -> 0).auditZeroAmplitude;
  auditCheck[auditMatrixPrefix <> "ZeroFrequencyPhysical", {auditNativePhysical /. omega -> 0,
    auditZeroPhysical, (auditNativePhysical /. omega -> 0) - auditZeroPhysical}, {-3, -1, 1}];
  auditCheck[auditMatrixPrefix <> "ZeroFrequencyReference", {auditNativeReference /. omega -> 0,
    auditZeroReference, (auditNativeReference /. omega -> 0) - auditZeroReference}, {-3, -1, 1}];
  (* Deliberate input mutations of the face evaluation, then a new independent
     solve; the unchanged physical face map is used to observe the error. *)
  Do[
    auditMutatedA = IdentityMatrix[3] + auditMutation (auditA - IdentityMatrix[3]);
    auditMutatedReference = auditFirstOrder[LinearSolve[auditMutatedA, auditPhysicalIndependent]];
    auditPhysical[auditMatrixPrefix <> "Mutation" <> ToString[auditMutation + 2],
      {auditMutatedReference, auditFirstOrder[auditA.auditMutatedReference - auditPhysicalIndependent]}, {-3, -1, 1}, 0],
    {auditMutation, {0, 2, -1}}];
  (* Distinct-leg mutations of the native conversion input, including both
     ordered factors. Keep their nonzero eta*sigma response separate. *)
  auditPhysical[auditMatrixPrefix <> "DistinctMomentumDifference", auditNativeReference -
    (auditNativeReference /. qLeg["Middle"] -> qLeg["In"]), {-3, -1, 1}, 0];
  auditPhysical[auditMatrixPrefix <> "SourceResiduals", {
    auditReference["TRACE_SOURCE_AFFINE_RESIDUAL"], Sequence @@ auditReference["TRACE_COEFFICIENT_EXTRACTION_RESIDUAL"]},
    {{-2, -2, 1}, {0, 0, 0}, {1, 0, 0}}, 0], {auditSign, faces}];

(* Native geometry/row joins on arbitrary background jets. The pressure response
   is an independent prescribed input; no truncated N15 energy is substituted. *)
auditRowCase = 0;
Do[
  auditRowCase++;
  auditDensity = Block[{densityKind = auditDensityKind, density4, density3},
    Scan[ReleaseHold, auditDensityAssignments]; density3];
  auditGeometry = If[auditRoute === "EULERIAN", graphGeometry[auditAnchor, auditSign], materialGeometry[auditAnchor, auditSign, 0]];
  auditLaw = faceLaws[auditGeometry, auditSign, muTheta, auditDensity];
  auditRows = If[auditRoute === "EULERIAN", eulerianSlabFace[auditLaw, auditMassInput], materialFaceFold[auditLaw, auditMassInput]];
  auditSlots = If[auditSign === 1, Take[pressureSlots, 2], Take[pressureSlots, -2]];
  auditActualTrace = linearShape[facePressureTrace[auditSign]];
  auditTraceC = D[auditActualTrace, #] & /@ auditSlots;
  auditRowC = Table[D[auditRows[[i]], slot], {i, 5}, {slot, auditSlots}];
  auditEmit["wlRepairedTraceNativeRowCase" <> ToString[auditRowCase],
    {auditAnchor, auditDensityKind, auditRoute, auditSign}, {0, 0, 0}, 0];
  Do[auditCheck["wlRepairedTraceNativeRow" <> ToString[auditRowCase] <> "Component" <> ToString[i],
    finish /@ {auditRowC[[i, 1]] auditTraceC[[2]], auditRowC[[i, 2]] auditTraceC[[1]],
      auditRowC[[i, 1]] auditTraceC[[2]] - auditRowC[[i, 2]] auditTraceC[[1]]},
    If[i <= 3, {1, 0, 0}, If[i === 4, {0, 1, 0}, {2, 0, 0}]], 0], {i, 5}],
  {auditCase, auditCases}, {auditAnchor, {auditCase[[1]]}}, {auditDensityKind, {auditCase[[2]]}},
  {auditRoute, {"EULERIAN", "MATERIAL"}}, {auditSign, faces}];

(* Uniform profile collapse: the coefficient of each Fourier profile is its
   constant-profile delta amplitude after the matching convolution. This is
   solely the constant-height regression; the generic ordered test above keeps
   all Fourier profiles, measures and momenta independent. *)
auditUniformRules = Join[Thread[momentum["kOut"] -> momentum["kIn"]],
  Thread[momentum["kMiddle"] -> momentum["kIn"]],
  {qLeg["Out"] -> qLeg["In"], qLeg["Middle"] -> qLeg["In"], sigmaW -> 0,
   profileTransform -> auditContrast, profileTransformLeft -> auditContrast, profileTransformRight -> auditContrast}];
auditUniformCase = 0;
Do[
  auditUniformCase++;
  auditDensity = Block[{densityKind = auditDensityKind, density4, density3},
    Scan[ReleaseHold, auditDensityAssignments]; density3] /. WBg -> W0 (1 + etaBg auditContrast);
  auditKernel = constructKernel[auditSign]; auditPhysicalResponse = responseFamilies[auditKernel];
  auditReference = referencePressureFamilies[auditPhysicalResponse, auditSign];
  auditDrive = auditSourceSolve["SOURCE"] /. rhoFace -> auditDensity;
  auditImagesByFamily = Table[Block[{face = auditSign, response = <|auditSign -> auditReference|>,
      sourceForGuard = retain[auditDrive], familyName = auditFamily, pressureImage, imageRules = <||>},
    Scan[ReleaseHold, auditImageAssignments]; Map[auditGradeExpression, imageRules]], {auditFamily, auditFamilies}];
  auditSlot = pressureSlots[[If[auditSign === 1, 1, 3]]]; auditJetSlot = pressureSlots[[If[auditSign === 1, 2, 4]]];
  auditReferenceImage = auditFirstOrder[Total[#[auditSlot] & /@ auditImagesByFamily] /. auditUniformRules];
  auditJetImage = auditFirstOrder[Total[#[auditJetSlot] & /@ auditImagesByFamily] /. auditUniformRules];
  auditTrace = linearShape[facePressureTrace[auditSign]] /. WBg -> W0 (1 + etaBg auditContrast);
  auditNativeFold = auditTrace /. {auditSlot -> auditReferenceImage, auditJetSlot -> auditJetImage};
  auditHeight = etaBg W0 auditContrast/2;
  auditWave = auditWaveAmplitude Exp[I (momentum["kIn"].{auditX1, auditX2, auditX3} + auditSign qLeg["In"] auditW - omega auditTime)];
  auditPhase = Exp[I (momentum["kIn"].{auditX1, auditX2, auditX3} - omega auditTime)];
  auditPField = -rhoM D[auditWave, auditTime]/auditPhase;
  auditVField = auditSign D[auditWave, auditW]/auditPhase;
  auditPAtFace = Simplify[auditPField /. auditW -> auditSign auditHeight];
  auditVAtFace = Simplify[auditVField /. auditW -> auditSign auditHeight];
  auditSolution = First[Solve[{auditFlux == lambdaA (muTheta/auditDensity - auditPAtFace/rhoM) + lambdaV vFace,
    auditVAtFace == vFace + auditFlux/rhoM}, {auditWaveAmplitude, auditFlux}]];
  auditFacePressure = Simplify[auditPAtFace /. auditSolution];
  auditReferencePressure = Simplify[(auditPField /. auditW -> 0) /. auditSolution];
  auditExactJet = Simplify[(D[auditPField, auditW] /. auditW -> 0) /. auditSolution];
  auditPrefix = "wlRepairedTraceUniform" <> ToString[auditUniformCase];
  auditCheck[auditPrefix <> "Face", {auditNativeFold, auditFacePressure, auditNativeFold - auditFacePressure}, {-2, -2, 1}, 1];
  auditCheck[auditPrefix <> "Reference", {auditReferenceImage, auditReferencePressure, auditReferenceImage - auditReferencePressure}, {-2, -2, 1}, 1];
  auditCheck[auditPrefix <> "Jet", {auditJetImage, auditExactJet, auditJetImage - auditExactJet}, {-3, -2, 1}, 1],
  {auditDensityKind, {"RHO4_CONSTANT", "RHOBR_CONSTANT"}}, {auditSign, faces}];

auditSummary = <|"records" -> Length[auditRecords], "matrixFaceCases" -> auditMatrixCase,
  "nativeRowCases" -> auditRowCase, "uniformCases" -> auditUniformCase,
  "residuals" -> Association[Map[# -> Map[ToString[#, InputForm] &, Flatten[{auditRecords[#]["value"]["retained"][[3]]}]] &, auditResidualKeys]],
  "mutationNonzeroComponents" -> Association[Table[
    key -> Count[Flatten[auditRecords[key]["value"]["retained"][[2]]], Except[0]],
    {key, Select[Keys[auditRecords], StringContainsQ[#, "Mutation"] &]}]],
  "orderedMomentumNonzeroComponents" -> Table[Count[Flatten[auditRecords["wlRepairedTraceMatrix" <> ToString[i] <> "DistinctMomentumDifference"]["value"]["retained"]], Except[0]], {i, auditMatrixCase}]|>;
auditEmit["wlRepairedTraceCompletion", KeyDrop[auditSummary, "residuals"], {0, 0, 0}];
Export[FileNameJoin[{auditRoot, "summary.json"}], auditSummary, "RawJSON"];
