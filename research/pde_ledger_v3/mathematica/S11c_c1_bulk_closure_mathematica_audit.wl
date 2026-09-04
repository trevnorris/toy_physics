$HistoryLength = 0;
ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};

(* ---------------------------------------------------------------------- *)
(* Flushed tag-stream boundary.                                          *)
(* ---------------------------------------------------------------------- *)

sharedObject[name_String] := HoldComplete[SharedS11cc1Object, name];
localObject[name_String] := HoldComplete[LocalS11cc1Object, name];

standardEmissionName[HoldComplete[SharedS11cc1Object, name_String]] :=
  "WL_S11CC1_" <> name;
standardEmissionName[HoldComplete[LocalS11cc1Object, name_String]] :=
  "WL_S11CC1_LOCAL_" <> name;

stripConditional[object_] := Replace[object,
  ConditionalExpression[value_, condition_] :>
    <|"CONDITIONAL_VALUE" -> value,
      "CONDITION_OPERAND" -> HoldForm[condition]|>, {0, Infinity}];

generatedSymbolNameQ[name_String] := StringMatchQ[name,
  RegularExpression[".*\\$[0-9]+$"]];
generatedSymbolBaseName[name_String] := StringReplace[name,
  RegularExpression["\\$[0-9]+$"] -> ""];
canonicalizeGeneratedSymbols[object_, scope_String] := Module[
  {bindings, sequence, scopeCode},
  bindings = <||>;
  sequence = 0;
  scopeCode = IntegerString[Abs[Hash[scope, "CRC32"]]];
  Replace[object, symbol_Symbol /;
      generatedSymbolNameQ[SymbolName[Unevaluated[symbol]]] :>
    RuleCondition[Module[{identity, baseName},
      identity = Context[Unevaluated[symbol]] <>
        SymbolName[Unevaluated[symbol]];
      baseName = generatedSymbolBaseName[
        SymbolName[Unevaluated[symbol]]];
      If[!KeyExistsQ[bindings, identity],
        sequence++;
        AssociateTo[bindings, identity -> Symbol[
          Context[Unevaluated[symbol]] <> baseName <> "Local" <>
            scopeCode <> "N" <> IntegerString[sequence]]]];
      bindings[identity]]], {0, Infinity}]
];

emittedNames = <||>;
localNames = {};
emit[key_, payload_] := Module[{name, rendered, stream, cleanPayload},
  name = standardEmissionName[key];
  If[!StringQ[name], Quit[90]];
  If[KeyExistsQ[emittedNames, name], Quit[91]];
  AssociateTo[emittedNames, name -> 1];
  If[StringStartsQ[name, "WL_S11CC1_LOCAL_"], AppendTo[localNames, name]];
  cleanPayload = canonicalizeGeneratedSymbols[
    stripConditional[payload], name];
  rendered = ToString[cleanPayload, InputForm, PageWidth -> Infinity];
  stream = First[$Output];
  WriteString[stream, name <> ": " <> rendered <> "\n"];
  Flush[stream];
];
emitShared[name_String, payload_] := emit[sharedObject[name], payload];
emitLocal[name_String, payload_] := emit[localObject[name], payload];

beginAssociationEmission[key_] := Module[{name, stream},
  name = standardEmissionName[key];
  If[!StringQ[name], Quit[90]];
  If[KeyExistsQ[emittedNames, name], Quit[91]];
  AssociateTo[emittedNames, name -> 1];
  If[StringStartsQ[name, "WL_S11CC1_LOCAL_"], AppendTo[localNames, name]];
  currentAssociationEmissionName = name;
  stream = First[$Output];
  WriteString[stream, name <> ": <|"];
  Flush[stream];
];

appendAssociationEmission[key_String, value_, first_:False] := Module[
  {stream, rendered},
  stream = First[$Output];
  rendered = ToString[canonicalizeGeneratedSymbols[
    stripConditional[value], currentAssociationEmissionName <> "|" <> key],
    InputForm, PageWidth -> Infinity];
  WriteString[stream, If[TrueQ[first], "", ", "] <>
    ToString[key, InputForm, PageWidth -> Infinity] <> " -> " <>
    rendered];
  Flush[stream];
];

endAssociationEmission[] := Module[{stream},
  stream = First[$Output];
  WriteString[stream, "|>\n"];
  Flush[stream];
  Clear[currentAssociationEmissionName];
];

relationalObject[left_, right_] := Inactive[Equal][left, right];

(* ---------------------------------------------------------------------- *)
(* Shared names, coordinates, branch data, and supplied laws.            *)
(* ---------------------------------------------------------------------- *)

braneDimension = 3;
spatialCoordinates = {xOne, xTwo, xThree};
materialCoordinates = {capitalXOne, capitalXTwo, capitalXThree};
momentumOutput = {kOne, kTwo, kThree};
momentumInput = {kPrimeOne, kPrimeTwo, kPrimeThree};
momentumTransfer = momentumOutput - momentumInput;
directions = Range[braneDimension];
faces = {1, -1};
anchorings = {"LAB_HELD", "MATERIAL_ADVECTED"};
densityRepresentatives = {"RHO4_CONSTANT", "RHOBR_CONSTANT"};
regimes = {"PROPAGATING", "EVANESCENT", "GRAZING"};

(* Reserved objects retain their mechanical lower-camel spellings. *)
reservedObjects = <|
  "W_bg" -> WBg, "w1_profile" -> w1Profile, "L_W" -> LW,
  "sigma_W" -> sigmaW, "eta_bg" -> etaBg,
  "epsilon_shape" -> epsilonShape,
  "rho_br" -> rhoBr,
  "rho_br_bg_rho4_constant" -> rhoBrBgRho4Constant,
  "Lambda_A_0" -> LambdaA0, "Lambda_V_0" -> LambdaV0,
  "Lambda_X_0" -> LambdaX0,
  "tau_A" -> tauA, "tau_V" -> tauV, "tau_X" -> tauX,
  "q_out" -> qOut, "omega" -> omega, "rho_m" -> rhoM,
  "c_s0" -> cS0, "W_0" -> W0, "e_W" -> eW,
  "v_bulk_normal_0" -> vBulkNormal0, "mu_theta" -> muTheta|>;

kappa = omega/cS0;
lambdaA = LambdaA0/(1 - I omega tauA);
lambdaV = LambdaV0/(1 - I omega tauV);
lambdaX = LambdaX0/(1 - I omega tauX);

bulkDispersionEquation[momentum_, normalWaveNumber_] :=
  relationalObject[normalWaveNumber^2 + momentum.momentum, kappa^2];
bulkRootCandidates[momentum_] := normalWaveNumber /. Solve[
  normalWaveNumber^2 + momentum.momentum == kappa^2,
  normalWaveNumber];

radiationBranchRecord[momentum_] := <|
  "DISPERSION_EQUATION" -> bulkDispersionEquation[momentum,
    qOut[omega, momentum]],
  "ROOT_CANDIDATES" -> bulkRootCandidates[momentum],
  "SELECTED_ROOT" -> qOut[omega, momentum],
  "REAL_AXIS_PROPAGATING_SELECTION" ->
    relationalObject[qOut[omega, momentum]/omega,
      Abs[qOut[omega, momentum]/omega]],
  "REAL_AXIS_EVANESCENT_SELECTION" ->
    relationalObject[Im[qOut[omega, momentum]],
      Abs[Im[qOut[omega, momentum]]]],
  "COMPLEX_CONTINUATION" ->
    AnalyticContinuation[
      qOut[omega + I upperRimInfinitesimal, momentum],
      FixedBranchPoints -> {
        -cS0 Sqrt[momentum.momentum],
        cS0 Sqrt[momentum.momentum]},
      ContinuationPath -> FixedRealPartRay]
|>;

chiVector = Through[{chiOne, chiTwo, chiThree}[
  xOne, xTwo, xThree, time]];
anchorArguments["LAB_HELD"] := spatialCoordinates;
anchorArguments["MATERIAL_ADVECTED"] := chiVector;
scaledAnchorArguments[anchoring_] := anchorArguments[anchoring]/LW;

widthAnsatz[anchoring_] := W0 (1 + etaBg
  w1Profile @@ scaledAnchorArguments[anchoring]);
widthBinding[anchoring_] := relationalObject[
  WBg @@ anchorArguments[anchoring], widthAnsatz[anchoring]];

profileTransform[anchoring_] := Inactive[FourierTransform][
  w1Profile @@ scaledAnchorArguments[anchoring],
  spatialCoordinates, momentumTransfer];

profileSource[anchoring_, slopeWeights_List] := Module[
  {profileHat, heightHat, slopeHat},
  profileHat = profileTransform[anchoring];
  heightHat = etaBg W0 profileHat/2;
  slopeHat = I sigmaW profileHat
    (slopeWeights (LW momentumTransfer))/2;
  <|
    "WIDTH_BINDING" -> widthBinding[anchoring],
    "SIGMA_BINDING" -> relationalObject[sigmaW, etaBg W0/LW],
    "HEIGHT_FOURIER" -> heightHat,
    "SLOPE_FOURIER" -> slopeHat,
    "SLOPE_WEIGHTS" -> slopeWeights|>
];

densityAtAnchor["RHO4_CONSTANT", anchoring_] :=
  rhoBrBgRho4Constant @@ anchorArguments[anchoring];
densityAtAnchor["RHOBR_CONSTANT", anchoring_] := rhoBr;
densityBinding["RHO4_CONSTANT", anchoring_] := relationalObject[
  rhoBrBgRho4Constant @@ anchorArguments[anchoring],
  (rhoBr/W0) (WBg @@ anchorArguments[anchoring])];
densityBinding["RHOBR_CONSTANT", anchoring_] :=
  HeldFixedConstant[rhoBr];

(* ---------------------------------------------------------------------- *)
(* Re-derived S11c-a face-shape substrate.                               *)
(* ---------------------------------------------------------------------- *)

gradient[scalar_] := D[scalar, #] & /@ spatialCoordinates;
divergence[vector_] := Sum[D[vector[[index]],
  spatialCoordinates[[index]]], {index, directions}];
shapeCoefficient[expression_, parameter_] :=
  Simplify[D[expression, parameter] /. parameter -> 0];
firstBackgroundShape[expression_] := Module[{scaled},
  scaled = expression /. {etaBg -> backgroundScale etaBg,
    sigmaW -> backgroundScale sigmaW};
  Simplify[(scaled /. backgroundScale -> 0) +
    shapeCoefficient[scaled, backgroundScale]]
];

zetaForFace[face_] := zetaCenter @@ Append[spatialCoordinates, time] +
  face deltaW @@ Append[spatialCoordinates, time]/2;

backgroundHeight[anchoring_, face_] :=
  face widthAnsatz[anchoring]/2;
waveHeight[face_] := zetaForFace[face];

orientedNormalFromLevelSet[anchoring_, face_, waveScale_] := Module[
  {height, levelGradient, oriented},
  height = backgroundHeight[anchoring, face] +
    waveScale waveHeight[face];
  levelGradient = Join[-gradient[height], {1}];
  oriented = face levelGradient/Sqrt[levelGradient.levelGradient];
  oriented
];

faceMeasureFromMap[anchoring_, face_, waveScale_] := Module[{height},
  height = backgroundHeight[anchoring, face] +
    waveScale waveHeight[face];
  Sqrt[1 + gradient[height].gradient[height]]
];

faceMapFromSource["LAB_HELD", face_, waveScale_] := Join[
  spatialCoordinates + waveScale
    Through[{uOne, uTwo, uThree} @@ Append[spatialCoordinates, time]],
  {face widthAnsatz["LAB_HELD"]/2 + waveScale waveHeight[face]}];
faceMapFromSource["MATERIAL_ADVECTED", face_, waveScale_] := Join[
  spatialCoordinates + waveScale
    Through[{uOne, uTwo, uThree} @@ Append[spatialCoordinates, time]],
  {face W0 (1 + etaBg
      w1Profile @@ (materialCoordinates/LW))/2 +
    waveScale waveHeight[face]}];

faceVelocityFromMap[anchoring_, face_, waveScale_] :=
  D[faceMapFromSource[anchoring, face, waveScale], time];

shiftedTraceFromMap[anchoring_, face_, waveScale_] := Module[
  {height, tracedField},
  height = backgroundHeight[anchoring, face] +
    waveScale waveHeight[face];
  tracedField = bulkFieldBackground @@ Join[spatialCoordinates,
      {height, time}] + waveScale
    bulkFieldPerturbation @@ Join[spatialCoordinates, {height, time}];
  shapeCoefficient[tracedField, waveScale]
];

backgroundWidthJet[anchoring_] := sigmaW Through[
  {w1JetOne, w1JetTwo, w1JetThree} @@
    scaledAnchorArguments[anchoring]];
waveFaceJet[face_] := Through[
  {zetaJetOne, zetaJetTwo, zetaJetThree}[
    face, xOne, xTwo, xThree, time]];
waveFaceVelocityVector[anchoring_, face_] := Module[
  {inPlaneVelocity, normalVelocityComponent, transportTerm},
  inPlaneVelocity = Through[{uVelocityOne, uVelocityTwo,
    uVelocityThree}[xOne, xTwo, xThree, time]];
  transportTerm = If[anchoring === "LAB_HELD",
    face backgroundWidthJet[anchoring].inPlaneVelocity/2, 0];
  normalVelocityComponent = transportTerm +
    zetaCenterVelocity[xOne, xTwo, xThree, time] +
    face deltaWVelocity[xOne, xTwo, xThree, time]/2;
  Join[inPlaneVelocity, {normalVelocityComponent}]
];
bulkVelocityBackgroundVector = ConstantArray[0, 4];
bulkVelocityPerturbationVector = Through[
  {bulkVelocityOne, bulkVelocityTwo, bulkVelocityThree,
    bulkVelocityNormal}[xOne, xTwo, xThree, time]];

shapeSubstrate[anchoring_, face_, density_] := Module[
  {backgroundSlope, waveSlope, levelGradient, normalExact,
   normalBackground, normalVariation, measureExact,
   measureBackground, measureVariation, conormalExact,
   faceVelocityExact, faceNormalVelocity, bulkVelocityExact,
   relativeFluxExact, relativeFluxVariation, kinematicExact,
   affinityExact, tractionExact, closureExact, faceShiftSource,
   faceShift, backgroundFaceLocation, backgroundTrace,
   perturbationTrace, backgroundTraceDerivative},
  backgroundSlope = face backgroundWidthJet[anchoring]/2;
  waveSlope = waveFaceJet[face];
  levelGradient = Join[-(backgroundSlope + waveScale waveSlope), {1}];
  normalExact = face levelGradient/Sqrt[levelGradient.levelGradient];
  normalBackground = firstBackgroundShape[
    normalExact /. waveScale -> 0];
  normalVariation = firstBackgroundShape[
    shapeCoefficient[normalExact, waveScale]];
  measureExact = Sqrt[1 +
    (backgroundSlope + waveScale waveSlope).
      (backgroundSlope + waveScale waveSlope)];
  measureBackground = firstBackgroundShape[
    measureExact /. waveScale -> 0];
  measureVariation = firstBackgroundShape[
    shapeCoefficient[measureExact, waveScale]];
  conormalExact = Total[MapThread[Times, {normalExact,
    {Inactive[D][conormalTestField, xOne],
      Inactive[D][conormalTestField, xTwo],
      Inactive[D][conormalTestField, xThree],
      Inactive[D][conormalTestField, w]}}]];
  faceVelocityExact = waveScale
    waveFaceVelocityVector[anchoring, face];
  faceNormalVelocity = Together[normalExact.faceVelocityExact];
  bulkVelocityExact = bulkVelocityBackgroundVector +
    waveScale bulkVelocityPerturbationVector;
  relativeFluxExact = rhoM Together[
    (bulkVelocityExact - faceVelocityExact).normalExact];
  relativeFluxVariation = shapeCoefficient[relativeFluxExact, waveScale];
  kinematicExact = Together[normalExact.bulkVelocityExact -
    faceNormalVelocity - relativeFluxExact/rhoM];
  affinityExact = waveScale muTheta/densityAtAnchor[
      density, anchoring] -
    waveScale pressureTrace[anchoring, face]/rhoM;
  tractionExact = -(waveScale pressureTrace[anchoring, face] +
      lambdaX affinityExact) normalExact;
  closureExact = relativeFluxExact - lambdaA affinityExact -
    lambdaV faceNormalVelocity;
  backgroundFaceLocation = backgroundHeight[anchoring, face];
  backgroundTrace = bulkFieldBackground @@ Join[spatialCoordinates,
    {backgroundFaceLocation, time}];
  perturbationTrace = bulkFieldPerturbation @@ Join[spatialCoordinates,
    {backgroundFaceLocation, time}];
  backgroundTraceDerivative = Derivative[0, 0, 0, 1, 0][
    bulkFieldBackground] @@ Join[spatialCoordinates,
      {backgroundFaceLocation, time}];
  faceShiftSource = backgroundTrace + waveScale
    (perturbationTrace + waveHeight[face]
      backgroundTraceDerivative);
  faceShift = firstBackgroundShape[
    shapeCoefficient[faceShiftSource, waveScale]];
  <|
    "T_A_FACE_NORMAL" -> <|"BACKGROUND" -> normalBackground,
      "SHAPE_DERIVATIVE" -> normalVariation|>,
    "T_A_PRIME_CONORMAL" -> <|"EXACT" -> conormalExact,
      "SHAPE_DERIVATIVE" -> firstBackgroundShape[
        shapeCoefficient[conormalExact, waveScale]]|>,
    "T_A_DOUBLEPRIME_FACE_MEASURE" -> <|
      "BACKGROUND" -> measureBackground,
      "SHAPE_DERIVATIVE" -> measureVariation|>,
    "T_B_FACE_VELOCITY" -> <|"EXACT" -> faceNormalVelocity,
      "SHAPE_DERIVATIVE" -> firstBackgroundShape[
        shapeCoefficient[faceNormalVelocity, waveScale]]|>,
    "T_C_RELATIVE_FLUX" ->
      firstBackgroundShape[relativeFluxVariation],
    "T_C_PRIME_KINEMATIC_BALANCE" ->
      firstBackgroundShape[shapeCoefficient[kinematicExact, waveScale]],
    "T_D_TRACTION" -> firstBackgroundShape[
      shapeCoefficient[tractionExact, waveScale]],
    "T_E_FACE_SHIFT" -> faceShift,
    "T_I_CLOSURE_SHAPE_DERIVATIVE" ->
      firstBackgroundShape[shapeCoefficient[closureExact, waveScale]]|>
];

(* The background true-area factor is obtained from the map before the
   first-shape truncation. *)
backgroundMeasureFirstShape[anchoring_, face_] := Module[
  {measure, scaled, slope},
  slope = face backgroundWidthJet[anchoring]/2;
  measure = Sqrt[1 + slope.slope];
  scaled = measure /. {etaBg -> geometryScale etaBg,
    sigmaW -> geometryScale sigmaW};
  Simplify[(scaled /. geometryScale -> 0) +
    shapeCoefficient[scaled, geometryScale]]
];

(* ---------------------------------------------------------------------- *)
(* Direct outgoing Neumann boundary matching.                            *)
(* ---------------------------------------------------------------------- *)

directBoundaryDerivation[source_Association, qOutput_, qInput_] := Module[
  {heightHat, slopeHat, amplitudeZero, amplitudeOne, inputAmplitude,
   flatEquation, flatRule, boundaryShift, conormalTilt,
   firstEquation, firstRule, shiftedTrace, pressureTraceCorrection,
   kernel, flatTrace, flatPressure, flatImpedance,
   originKernel},
  heightHat = source["HEIGHT_FOURIER"];
  slopeHat = source["SLOPE_FOURIER"];
  flatEquation = I qInput amplitudeZero == inputAmplitude;
  flatRule = First[Solve[flatEquation, amplitudeZero]];
  boundaryShift = -qInput^2 heightHat amplitudeZero;
  conormalTilt = -slopeHat.(I momentumInput) amplitudeZero;
  firstEquation = I qOutput amplitudeOne + boundaryShift +
    conormalTilt == 0;
  firstRule = First[Solve[firstEquation /. flatRule, amplitudeOne]];
  shiftedTrace = heightHat I qInput amplitudeZero;
  pressureTraceCorrection = Simplify[
    I rhoM omega (amplitudeOne + shiftedTrace) / inputAmplitude /.
      flatRule /. firstRule];
  kernel = Together[pressureTraceCorrection];
  flatTrace = amplitudeZero /. flatRule;
  flatPressure = I rhoM omega flatTrace;
  flatImpedance = Simplify[flatPressure/inputAmplitude];
  originKernel = <|
    "BOUNDARY_SHIFT" -> Simplify[I rhoM omega
      (amplitudeOne /. First[Solve[
        I qOutput amplitudeOne + boundaryShift == 0,
        amplitudeOne]] /. flatRule)/inputAmplitude],
    "CONORMAL_TILT" -> Simplify[I rhoM omega
      (amplitudeOne /. First[Solve[
        I qOutput amplitudeOne + conormalTilt == 0,
        amplitudeOne]] /. flatRule)/inputAmplitude],
    "SHIFTED_TRACE" -> Simplify[I rhoM omega
      shiftedTrace /. flatRule]/inputAmplitude|>;
  <|
    "FLAT_BOUNDARY_EQUATION" -> HoldForm[flatEquation],
    "FIRST_BOUNDARY_EQUATION" -> HoldForm[firstEquation],
    "FLAT_IMPEDANCE" -> flatImpedance,
    "KERNEL" -> kernel,
    "TERM_ORIGINS" -> originKernel,
    "ORIGIN_SUM_RESIDUAL" ->
      Together[Total[Values[originKernel]] - kernel],
    "SHIFTED_TRACE" -> shiftedTrace|>
];

(* ---------------------------------------------------------------------- *)
(* Radiation-selected layer-potential/Dirichlet route and inverse.       *)
(* ---------------------------------------------------------------------- *)

layerPotentialDerivation[source_Association, qOutput_, qInput_] := Module[
  {heightHat, slopeHat, dirichletInput, flatPotentialAmplitude,
   correctionAmplitude, flatDirichletEquation, flatRule,
   boundaryValueEquation, correctionRule, conormalTrace,
   dtnFlat, dtnFirst, inverseZeroOutput, inverseZeroInput,
   inverseIdentityEquation, inverseFirstRule, ntdFirst,
   impedanceFirst},
  heightHat = source["HEIGHT_FOURIER"];
  slopeHat = source["SLOPE_FOURIER"];
  flatDirichletEquation = flatPotentialAmplitude == dirichletInput;
  flatRule = First[Solve[flatDirichletEquation,
    flatPotentialAmplitude]];
  boundaryValueEquation = correctionAmplitude +
    heightHat I qInput flatPotentialAmplitude == 0;
  correctionRule = First[Solve[
    boundaryValueEquation /. flatRule, correctionAmplitude]];
  conormalTrace = I qOutput correctionAmplitude -
    qInput^2 heightHat flatPotentialAmplitude -
    slopeHat.(I momentumInput) flatPotentialAmplitude;
  dtnFlat = I qInput;
  dtnFirst = Simplify[conormalTrace/dirichletInput /.
    flatRule /. correctionRule];
  inverseZeroOutput = 1/(I qOutput);
  inverseZeroInput = 1/(I qInput);
  inverseIdentityEquation = inverseFirstCoefficient +
    inverseZeroOutput dtnFirst inverseZeroInput == 0;
  inverseFirstRule = First[Solve[inverseIdentityEquation,
    inverseFirstCoefficient]];
  ntdFirst = Simplify[inverseFirstCoefficient /. inverseFirstRule];
  impedanceFirst = Together[I rhoM omega ntdFirst];
  <|
    "OUTGOING_GREEN_SYMBOL_OUTPUT" -> inverseZeroOutput,
    "OUTGOING_GREEN_SYMBOL_INPUT" -> inverseZeroInput,
    "DIRICHLET_BOUNDARY_EQUATION" -> HoldForm[boundaryValueEquation],
    "DTN_FLAT" -> dtnFlat,
    "DTN_FIRST_SHAPE" -> dtnFirst,
    "INVERSE_IDENTITY_EQUATION" -> HoldForm[inverseIdentityEquation],
    "KERNEL" -> impedanceFirst|>
];

sourceForCase[anchoring_, slopeWeights_List] :=
  profileSource[anchoring, slopeWeights];

directCase[anchoring_, slopeWeights_List, qOutput_, qInput_] :=
  directBoundaryDerivation[
    sourceForCase[anchoring, slopeWeights], qOutput, qInput];
layerCase[anchoring_, slopeWeights_List, qOutput_, qInput_] :=
  layerPotentialDerivation[
    sourceForCase[anchoring, slopeWeights], qOutput, qInput];

qOutputLive = qOut[omega, momentumOutput];
qInputLive = qOut[omega, momentumInput];
unitSlopeWeights = ConstantArray[1, braneDimension];

directMain = Association@Table[anchoring -> directCase[anchoring,
    unitSlopeWeights, qOutputLive, qInputLive],
  {anchoring, anchorings}];
layerMain = Association@Table[anchoring -> layerCase[anchoring,
    unitSlopeWeights, qOutputLive, qInputLive],
  {anchoring, anchorings}];

(* A second, independent flat half-space solve is retained for the S11b
   regression operand. *)
deriveFlatS11bImpedance[momentum_] := Module[
  {bulkAmplitudeB, outwardVelocityB, equationB, ruleB,
   pressureB},
  equationB = I qOut[omega, momentum] bulkAmplitudeB ==
    outwardVelocityB;
  ruleB = First[Solve[equationB, bulkAmplitudeB]];
  pressureB = I rhoM omega bulkAmplitudeB /. ruleB;
  Simplify[pressureB/outwardVelocityB]
];

flatSymbolC1 = directBoundaryDerivation[
  <|"HEIGHT_FOURIER" -> 0,
    "SLOPE_FOURIER" -> ConstantArray[0, braneDimension]|>,
  qOutputLive, qOutputLive]["FLAT_IMPEDANCE"];
flatSymbolS11b = deriveFlatS11bImpedance[momentumOutput];

(* Emit the first completed object immediately; later work is deliberately
   streamed behind it.  Its units are derived here from the supplied bulk
   laws so the payload is complete at this point in evaluation. *)
earlyDimensionL = {1, 0, 0};
earlyDimensionT = {0, 1, 0};
earlyDimensionM = {0, 0, 1};
earlyDimensionVelocity = earlyDimensionL - earlyDimensionT;
earlyDimensionRhoM = earlyDimensionM - 4 earlyDimensionL;
earlyDimensionPhi = earlyDimensionVelocity + earlyDimensionL;
earlyDimensionPressure = earlyDimensionRhoM + earlyDimensionPhi -
  earlyDimensionT;
earlyDimensionImpedance = earlyDimensionPressure -
  earlyDimensionVelocity;
earlyFlatGrades = DeleteDuplicates[({Exponent[#, epsilonShape],
      Exponent[#, etaBg], Exponent[#, sigmaW]} & /@
    If[Head[Expand[flatSymbolC1]] === Plus,
      List @@ Expand[flatSymbolC1], {Expand[flatSymbolC1]}])];
emitShared["DTN_FLAT_SYMBOL", <|
  "C1_INDEPENDENT_DERIVATION" -> <|
    "EXPRESSION" -> flatSymbolC1,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> earlyFlatGrades,
    "DIMENSION_L_T_M" -> earlyDimensionImpedance|>,
  "RADIATION_BRANCH" -> radiationBranchRecord[momentumOutput],
  "PRESSURE_PER_OUTWARD_NORMAL_VELOCITY" ->
    relationalObject[flatPressureRatio, flatSymbolC1]|>];
Clear[earlyDimensionL, earlyDimensionT, earlyDimensionM,
  earlyDimensionVelocity, earlyDimensionRhoM, earlyDimensionPhi,
  earlyDimensionPressure, earlyDimensionImpedance, earlyFlatGrades];

substrateByAnchoringAndFace = Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face] <> "|" <> density) ->
    shapeSubstrate[anchoring, face, density],
  {anchoring, anchorings}, {face, faces},
  {density, densityRepresentatives}];

operatorCompositionFromDerivation[anchoring_] := Module[
  {source, nZero, gZero, multiplication, dtnVariation,
   inverseVariation, flatOperator, firstOperator},
  source = sourceForCase[anchoring, unitSlopeWeights];
  nZero = FourierMultiplier[1/(I qOut[omega, momentumOutput])];
  gZero = FourierMultiplier[I qOut[omega, momentumOutput]];
  multiplication = MultiplicationOperator[
    (WBg @@ anchorArguments[anchoring] - W0)/2];
  dtnVariation = OperatorSum[
    -OperatorComposition[gZero, multiplication, gZero],
    -DivergenceOperator[multiplication, GradientOperator],
    -ScalarMultiplier[kappa^2] multiplication];
  inverseVariation = -OperatorComposition[
    nZero, dtnVariation, nZero];
  flatOperator = FourierMultiplier[flatSymbolC1];
  firstOperator = OperatorComposition[
    ScalarMultiplier[I rhoM omega], inverseVariation];
  OperatorSum[flatOperator, firstOperator]
];

(* ---------------------------------------------------------------------- *)
(* Regime, adjoint, parity, rigid-shift, and branch controls.             *)
(* ---------------------------------------------------------------------- *)

qRegimeValue["PROPAGATING", momentum_] :=
  Sign[omega] Sqrt[kappa^2 - momentum.momentum];
qRegimeValue["EVANESCENT", momentum_] :=
  I Sqrt[momentum.momentum - kappa^2];
qRegimeValue["GRAZING", momentum_] := 0;

regimeCondition["PROPAGATING", momentum_] :=
  Inactive[Greater][kappa^2 - momentum.momentum, 0];
regimeCondition["EVANESCENT", momentum_] :=
  Inactive[Less][kappa^2 - momentum.momentum, 0];
regimeCondition["GRAZING", momentum_] :=
  relationalObject[kappa^2 - momentum.momentum, 0];

grazingExpression[expression_, outputRegime_, inputRegime_] := Module[
  {generic, outputApplied, inputApplied},
  generic = expression /. {qOutputLive -> qOutputProbe,
    qInputLive -> qInputProbe};
  outputApplied = If[outputRegime === "GRAZING",
    Quiet[Limit[generic, qOutputProbe -> 0]],
    generic /. qOutputProbe -> qRegimeValue[outputRegime,
      momentumOutput]];
  inputApplied = If[inputRegime === "GRAZING",
    Quiet[Limit[outputApplied, qInputProbe -> 0]],
    outputApplied /. qInputProbe -> qRegimeValue[inputRegime,
      momentumInput]];
  inputApplied
];

regimePairRecord[expression_, outputRegime_, inputRegime_] := <|
  "OUTPUT_CONDITION" -> regimeCondition[outputRegime, momentumOutput],
  "INPUT_CONDITION" -> regimeCondition[inputRegime, momentumInput],
  "OUTPUT_BRANCH_VALUE" -> qRegimeValue[outputRegime, momentumOutput],
  "INPUT_BRANCH_VALUE" -> qRegimeValue[inputRegime, momentumInput],
  "KERNEL_OBJECT" ->
    objectRecord[grazingExpression[expression, outputRegime,
      inputRegime], dimensionImpedance],
  "OUTPUT_GRAZING_POLE_TEST" ->
    relationalObject[1/qOutputProbe, 0],
  "INPUT_GRAZING_POLE_TEST" ->
    relationalObject[1/qInputProbe, 0]|>;

swapMomentumLegs[expression_] := expression /. Join[
  Thread[momentumOutput -> momentumSwapOutput],
  Thread[momentumInput -> momentumOutput],
  Thread[momentumSwapOutput -> momentumInput]];

weightedAdjointKernel[expression_, anchoring_, face_] := Module[
  {measure, swapped},
  measure = backgroundMeasureFirstShape[anchoring, face];
  swapped = Conjugate[swapMomentumLegs[expression]];
  Simplify[InverseBoundaryMeasure[measure] swapped
    BoundaryMeasure[measure]]
];

hermitianKernel[expression_, anchoring_, face_] := Module[{adjoint},
  adjoint = weightedAdjointKernel[expression, anchoring, face];
  Simplify[(expression + adjoint)/2]
];
reactiveKernel[expression_, anchoring_, face_] := Module[{adjoint},
  adjoint = weightedAdjointKernel[expression, anchoring, face];
  Simplify[(expression - adjoint)/(2 I)]
];

faceToParityMatrix = {
  {1/2, 1},
  {1/2, -1}};
parityOperatorMatrix[zPlus_, zMinus_] := Simplify[
  Inverse[faceToParityMatrix].DiagonalMatrix[{zPlus, zMinus}].
    faceToParityMatrix];

rigidShiftSource[anchoring_] := <|
  "HEIGHT_FOURIER" -> etaBg W0 rigidProfileAmplitude/2,
  "SLOPE_FOURIER" -> ConstantArray[0, braneDimension],
  "WIDTH_BINDING" -> widthBinding[anchoring]|>;
rigidShiftKernel[anchoring_] := Module[{derived},
  derived = directBoundaryDerivation[rigidShiftSource[anchoring],
    qRigid, qRigid];
  Simplify[derived["KERNEL"] /. {
    momentumInput -> momentumOutput,
    momentumTransfer -> ConstantArray[0, braneDimension]}]
];
deriveTranslatedFlatImpedance[shift_] := Module[
  {translatedAmplitude, translatedVelocity, translatedPressure,
   translatedInput, translatedRule},
  translatedVelocity = D[translatedAmplitude Exp[
      I qOut[omega, momentumOutput]
        (outwardCoordinate - shift)], outwardCoordinate] /.
    outwardCoordinate -> shift;
  translatedPressure = I rhoM omega translatedAmplitude Exp[
      I qOut[omega, momentumOutput]
        (outwardCoordinate - shift)] /.
    outwardCoordinate -> shift;
  translatedRule = First[Solve[
    translatedVelocity == translatedInput, translatedAmplitude]];
  Simplify[(translatedPressure /. translatedRule)/translatedInput]
];
flatRigidShiftReference = shapeCoefficient[
  deriveTranslatedFlatImpedance[rigidShiftParameter],
  rigidShiftParameter];

branchCorruptionRules = <|
  "SIGNFLIP_INPUT" -> {qOutputLive -> qOutputLive,
    qInputLive -> -qInputLive},
  "SIGNFLIP_OUTPUT" -> {qOutputLive -> -qOutputLive,
    qInputLive -> qInputLive},
  "MOMENTUMFREEZE_OUTPUT" -> {qOutputLive -> qInputLive,
    qInputLive -> qInputLive},
  "MOMENTUMFREEZE_INPUT" -> {qOutputLive -> qOutputLive,
    qInputLive -> qOutputLive}|>;

recomputeBranchKernel[anchoring_, corruption_] := Module[
  {rules, qOutCorrupted, qInCorrupted},
  rules = branchCorruptionRules[corruption];
  qOutCorrupted = qOutputLive /. rules;
  qInCorrupted = qInputLive /. rules;
  directCase[anchoring, unitSlopeWeights,
    qOutCorrupted, qInCorrupted]["KERNEL"]
];

(* ---------------------------------------------------------------------- *)
(* Face closure solved from the supplied laws.                           *)
(* ---------------------------------------------------------------------- *)

flatFaceResponseSolve[zSymbol_, rhoBackground_] := Module[
  {pressureUnknown, fluxUnknown, faceVelocityDrive, chemicalDrive,
   bulkRelation, fluxRelation, solution, pressureSolution,
   fluxSolution, affinitySolution, tractionScalarSolution},
  faceVelocityDrive = epsilonShape faceVelocityInput;
  chemicalDrive = epsilonShape muTheta/rhoBackground;
  bulkRelation = pressureUnknown == zSymbol
    (faceVelocityDrive + fluxUnknown/rhoM);
  fluxRelation = fluxUnknown == lambdaA
    (chemicalDrive - pressureUnknown/rhoM) +
    lambdaV faceVelocityDrive;
  solution = First[Solve[{bulkRelation, fluxRelation},
    {pressureUnknown, fluxUnknown}]];
  pressureSolution = Together[pressureUnknown /. solution];
  fluxSolution = Together[fluxUnknown /. solution];
  affinitySolution = Together[chemicalDrive - pressureSolution/rhoM];
  tractionScalarSolution = Together[-(pressureSolution +
    lambdaX affinitySolution)];
  <|
    "SOURCE_EQUATIONS" -> HoldForm[{bulkRelation, fluxRelation}],
    "PRESSURE" -> pressureSolution,
    "RELATIVE_FLUX" -> fluxSolution,
    "AFFINITY" -> affinitySolution,
    "TRACTION_SCALAR" -> tractionScalarSolution,
    "COEFFICIENT_MATRIX" -> CoefficientArrays[
      Subtract @@@ ({bulkRelation, fluxRelation} /. Equal -> List),
      {pressureUnknown, fluxUnknown}][[2]]|>
];

curvedFaceResponseKernelSolve[zOutput_, zInput_, zFirst_,
    rhoBackground_, anchoring_, face_, memoryKernels_:Automatic] := Module[
  {pressureZero, fluxZero, pressureOne, fluxOne,
   faceVelocityDrive, chemicalDrive, flatBulk, flatClosure,
   curvedBulk, curvedClosure, solution, pZero, jZero, pOne, jOne,
   affinityZero, affinityOne, tractionScalarZero,
   tractionScalarOne, normalZero, normalOne, tractionVectorOne,
   source, slope, rawNormal, normalDerivative, lambdaASolve,
   lambdaVSolve, lambdaXSolve},
  {lambdaASolve, lambdaVSolve, lambdaXSolve} = If[
    memoryKernels === Automatic, {lambdaA, lambdaV, lambdaX},
    memoryKernels];
  faceVelocityDrive = epsilonShape faceVelocityInput;
  chemicalDrive = epsilonShape muTheta/rhoBackground;
  flatBulk = pressureZero == zInput
    (faceVelocityDrive + fluxZero/rhoM);
  flatClosure = fluxZero == lambdaASolve
    (chemicalDrive - pressureZero/rhoM) +
    lambdaVSolve faceVelocityDrive;
  curvedBulk = pressureOne == zOutput fluxOne/rhoM +
    zFirst (faceVelocityDrive + fluxZero/rhoM);
  curvedClosure = fluxOne == -lambdaASolve pressureOne/rhoM;
  solution = First[Solve[
    {flatBulk, flatClosure, curvedBulk, curvedClosure},
    {pressureZero, fluxZero, pressureOne, fluxOne}]];
  {pZero, jZero, pOne, jOne} = Together[
    {pressureZero, fluxZero, pressureOne, fluxOne} /. solution];
  affinityZero = Together[chemicalDrive - pZero/rhoM];
  affinityOne = Together[-pOne/rhoM];
  tractionScalarZero = Together[-(pZero + lambdaXSolve affinityZero)];
  tractionScalarOne = Together[-(pOne + lambdaXSolve affinityOne)];
  source = sourceForCase[anchoring, unitSlopeWeights];
  slope = source["SLOPE_FOURIER"];
  rawNormal = Join[-normalScale slope, {face}]/
    Sqrt[1 + normalScale^2 slope.slope];
  normalZero = rawNormal /. normalScale -> 0;
  normalDerivative = shapeCoefficient[rawNormal, normalScale];
  normalOne = normalDerivative;
  tractionVectorOne = Simplify[
    tractionScalarOne normalZero + tractionScalarZero normalOne];
  <|
    "SOURCE_EQUATIONS" -> HoldForm[
      {flatBulk, flatClosure, curvedBulk, curvedClosure}],
    "PRESSURE_FLAT" -> pZero,
    "RELATIVE_FLUX_FLAT" -> jZero,
    "TRACTION_FLAT_VECTOR" -> tractionScalarZero normalZero,
    "PRESSURE_FIRST_SHAPE_KERNEL" -> pOne,
    "RELATIVE_FLUX_FIRST_SHAPE_KERNEL" -> jOne,
    "TRACTION_FIRST_SHAPE_KERNEL_VECTOR" -> tractionVectorOne,
    "AFFINITY_FLAT" -> affinityZero,
    "AFFINITY_FIRST_SHAPE_KERNEL" -> affinityOne|>
];

faceResponseCase[anchoring_, face_, density_, kernel_:Automatic,
    qOutput_:qOutputLive, qInput_:qInputLive,
    memoryKernels_:Automatic] :=
 Module[{zFirst, rhoBackground},
  zFirst = If[kernel === Automatic,
    directMain[anchoring]["KERNEL"], kernel];
  rhoBackground = densityAtAnchor[density, anchoring];
  Append[curvedFaceResponseKernelSolve[
    flatSymbolC1 /. qOutputLive -> qOutput,
    flatSymbolC1 /. qOutputLive -> qInput,
    zFirst, rhoBackground, anchoring, face, memoryKernels],
    "DENSITY_BINDING" -> densityBinding[density, anchoring]]
];

mainFaceResponses = Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face] <> "|" <> density) ->
    faceResponseCase[anchoring, face, density],
  {anchoring, anchorings}, {face, faces},
  {density, densityRepresentatives}];

responseOutputDimension[name_] := If[
  StringContainsQ[name, "FLUX"], dimensionFlux,
  If[StringContainsQ[name, "TRACTION"], dimensionTraction,
    dimensionPressure]];
responseCoefficients[response_Association] := AssociationMap[
  Function[name, <|
    "FACE_VELOCITY_COEFFICIENT" -> objectRecord[
      Together[Coefficient[response[name], faceVelocityInput]],
      responseOutputDimension[name] - dimensionVelocity],
    "MU_THETA_COEFFICIENT" -> objectRecord[
      Together[Coefficient[response[name], muTheta]],
      responseOutputDimension[name] - dimensionMuTheta]|>],
  {"PRESSURE_FLAT", "RELATIVE_FLUX_FLAT",
   "TRACTION_FLAT_VECTOR", "PRESSURE_FIRST_SHAPE_KERNEL",
   "RELATIVE_FLUX_FIRST_SHAPE_KERNEL",
   "TRACTION_FIRST_SHAPE_KERNEL_VECTOR"}];

(* An independently solved flat B0c object is used only by regressions. *)
deriveFlatS11bFaceResponse[momentum_, rhoBackground_] := Module[
  {zReference, pressureReference, fluxReference, velocityReference,
   chemicalReference, equationsReference, solutionReference},
  zReference = deriveFlatS11bImpedance[momentum];
  velocityReference = epsilonShape faceVelocityInput;
  chemicalReference = epsilonShape muTheta/rhoBackground;
  equationsReference = {
    pressureReference == zReference
      (velocityReference + fluxReference/rhoM),
    fluxReference == lambdaA
      (chemicalReference - pressureReference/rhoM) +
      lambdaV velocityReference};
  solutionReference = First[Solve[equationsReference,
    {pressureReference, fluxReference}]];
  <|
    "PRESSURE" -> Together[pressureReference /. solutionReference],
    "RELATIVE_FLUX" -> Together[fluxReference /. solutionReference],
    "TRACTION_SCALAR" -> Together[-((pressureReference /.
      solutionReference) + lambdaX (chemicalReference -
      (pressureReference /. solutionReference)/rhoM))],
    "SOURCE_EQUATIONS" -> HoldForm[equationsReference]|>
];

(* ---------------------------------------------------------------------- *)
(* Algebraic loci and typed tests.                                       *)
(* ---------------------------------------------------------------------- *)

flatResponseForLocus = flatFaceResponseSolve[flatSymbolC1, rhoBr];
flatCoefficientMatrix = flatResponseForLocus["COEFFICIENT_MATRIX"];
flatDegeneracyPolynomial = Together[Det[flatCoefficientMatrix]];
flatDegeneracyEquationActive = flatDegeneracyPolynomial == 0;
flatDegeneracyEquation = relationalObject[flatDegeneracyPolynomial, 0];
flatDegeneracySolutions = Solve[flatDegeneracyEquationActive,
  LambdaA0];

identityNumerator = Numerator[Together[flatDegeneracyPolynomial]];
identityComputed = TrueQ[Simplify[identityNumerator === 0]];
identityStatus = If[identityComputed, "PROVED_TRUE",
  If[PolynomialQ[identityNumerator, LambdaA0], "PROVED_FALSE",
    "UNDECIDED"]];
inconsistentComputed = TrueQ[flatDegeneracySolutions === {}];
inconsistentStatus = If[inconsistentComputed, "PROVED_TRUE",
  "PROVED_FALSE"];

realAdmissibilityForBranch[branchRule_] := Module[
  {branchEquations, realTest, reduced, status},
  branchEquations = Thread[Equal @@@ (List @@@ branchRule)];
  realTest = Inactive[Reduce][
    And @@ Join[branchEquations,
      {Element[{LambdaA0, tauA, omega, rhoM, cS0}, Reals],
       tauA >= 0}],
    {LambdaA0, tauA, omega, rhoM, cS0}, Reals];
  reduced = Quiet[Reduce[
    And @@ Join[branchEquations,
      {Element[{LambdaA0, tauA, omega, rhoM, cS0}, Reals],
       tauA >= 0}],
    {LambdaA0, tauA, omega, rhoM, cS0}, Reals]];
  status = Which[
    TrueQ[reduced === False], "EXCLUDED",
    TrueQ[reduced === True], "ADMISSIBLE",
    Head[reduced] =!= Reduce, "ADMISSIBLE",
    True, "UNDECIDED"];
  <|"BRANCH" -> branchRule, "STATUS_TOKEN" -> status,
    "TEST_OBJECT" -> realTest,
    "OPERANDS" -> <|"BRANCH_EQUATIONS" ->
      (relationalObject[#[[1]], #[[2]]] & /@ branchEquations),
      "REAL_PREMISES" -> HoldForm[
        {LambdaA0, tauA, omega, rhoM, cS0} ∈ Reals && tauA >= 0]|>,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
    "DIMENSION_L_T_M" -> {0, 0, 0}|>
];

realAdmissibilityPayload = If[flatDegeneracySolutions === {},
  {<|"BRANCH" -> EmptySolutionSet,
    "STATUS_TOKEN" -> "EXCLUDED",
    "TEST_OBJECT" -> Inactive[Reduce][flatDegeneracyEquationActive,
      LambdaA0, Reals],
    "OPERANDS" -> {flatDegeneracyEquation},
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
    "DIMENSION_L_T_M" -> {0, 0, 0}|>},
  realAdmissibilityForBranch /@ flatDegeneracySolutions];

(* ---------------------------------------------------------------------- *)
(* Port Hermitian form, memory limits, and energy construction.          *)
(* ---------------------------------------------------------------------- *)

portMapFromResponse[response_Association, rhoBackground_,
    lambdaXForPort_:Automatic] := Module[
  {pressure, flux, affinity, totalTractionPressure, outputVector,
   inputVector, matrix, muSubstitution, lambdaXPort},
  lambdaXPort = If[lambdaXForPort === Automatic,
    lambdaX, lambdaXForPort];
  muSubstitution = muTheta -> muSPort rhoBackground;
  pressure = response["PRESSURE"] /. muSubstitution;
  flux = response["RELATIVE_FLUX"] /. muSubstitution;
  affinity = response["AFFINITY"] /. muSubstitution;
  totalTractionPressure = Together[pressure + lambdaXPort affinity];
  outputVector = {totalTractionPressure, flux};
  inputVector = {faceVelocityInput, muSPort};
  matrix = Table[Together[Coefficient[outputVector[[row]],
      inputVector[[column]]]], {row, 2}, {column, 2}];
  <|"PORT_MATRIX" -> matrix,
    "HERMITIAN_FORM" -> Simplify[(matrix + ConjugateTranspose[matrix])/2],
    "INPUT_VECTOR" -> inputVector,
    "OUTPUT_VECTOR" -> outputVector|>
];

curvedPortResponse[anchoring_, face_, density_,
    memoryKernels_:Automatic] := Module[
  {rhoBackground, response, lambdaXForPort},
  rhoBackground = densityAtAnchor[density, anchoring];
  response = If[memoryKernels === Automatic,
    mainFaceResponses[
      anchoring <> "|FACE_" <> ToString[face] <> "|" <> density],
    faceResponseCase[anchoring, face, density, Automatic,
      qOutputLive, qInputLive, memoryKernels]];
  lambdaXForPort = If[memoryKernels === Automatic,
    lambdaX, memoryKernels[[3]]];
  portMapFromResponse[<|
    "PRESSURE" -> response["PRESSURE_FLAT"] +
      response["PRESSURE_FIRST_SHAPE_KERNEL"],
    "RELATIVE_FLUX" -> response["RELATIVE_FLUX_FLAT"] +
      response["RELATIVE_FLUX_FIRST_SHAPE_KERNEL"],
    "AFFINITY" -> response["AFFINITY_FLAT"] +
      response["AFFINITY_FIRST_SHAPE_KERNEL"]|>, rhoBackground,
      lambdaXForPort]
];

portHermitianBaseCases = Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face] <> "|" <> density) ->
    curvedPortResponse[anchoring, face, density],
  {anchoring, anchorings}, {face, faces},
  {density, densityRepresentatives}];

principalMinorObjects[matrix_] := Module[{minors},
  minors = {matrix[[1, 1]], Det[matrix]};
  <|"PRINCIPAL_MINORS" -> minors,
    "NONNEGATIVITY_RELATIONS" ->
      (Inactive[GreaterEqual][#, 0] & /@ minors)|>
];

realMemoryParameterRules = {
  Conjugate[omega] -> omega, Conjugate[rhoM] -> rhoM,
  Conjugate[cS0] -> cS0, Conjugate[rhoBr] -> rhoBr,
  Conjugate[LambdaA0] -> LambdaA0,
  Conjugate[LambdaV0] -> LambdaV0,
  Conjugate[LambdaX0] -> LambdaX0,
  Conjugate[tauA] -> tauA, Conjugate[tauV] -> tauV,
  Conjugate[tauX] -> tauX};
memoryLawEndpoint[timeSymbol_, endpoint_] :=
  memoryLawEndpoint[timeSymbol, endpoint] = Map[
    Quiet[FullSimplify[Limit[#, timeSymbol -> endpoint],
      Assumptions -> Element[omega, Reals] && omega != 0]] &,
    {lambdaA, lambdaV, lambdaX}];
memoryKernelHandles = {memoryKernelA, memoryKernelV, memoryKernelX};
memoryPortBase[anchoring_, face_, density_] :=
  memoryPortBase[anchoring, face, density] =
    curvedPortResponse[anchoring, face, density, memoryKernelHandles];
memoryRegimeHermitian[anchoring_, face_, density_, outputRegime_,
    inputRegime_] := memoryRegimeHermitian[anchoring, face, density,
      outputRegime, inputRegime] = Map[
    grazingExpression[#, outputRegime, inputRegime] &,
    memoryPortBase[anchoring, face, density]["HERMITIAN_FORM"], {2}];
memoryMatrixEndpoint[anchoring_, face_, density_, outputRegime_,
    inputRegime_, timeSymbol_, endpoint_] :=
  memoryMatrixEndpoint[anchoring, face, density, outputRegime,
    inputRegime, timeSymbol, endpoint] = Module[
  {hermitianAtEndpoint, endpointRules},
  endpointRules = Thread[memoryKernelHandles ->
    memoryLawEndpoint[timeSymbol, endpoint]];
  hermitianAtEndpoint = memoryRegimeHermitian[anchoring, face,
    density, outputRegime, inputRegime] /. endpointRules;
  hermitianAtEndpoint /. realMemoryParameterRules
];
memoryLimitPackage[anchoring_, face_, density_, outputRegime_,
    inputRegime_] := memoryLimitPackage[anchoring, face, density,
      outputRegime, inputRegime] = Association@Table[
  ToString[timeSymbol, InputForm] -> <|
    "OMEGA_TAU_OPERAND" -> omega timeSymbol,
    "ZERO_LIMIT" -> memoryMatrixEndpoint[anchoring, face, density,
      outputRegime, inputRegime, timeSymbol, 0],
    "INFINITE_LIMIT" -> memoryMatrixEndpoint[anchoring, face, density,
      outputRegime, inputRegime, timeSymbol, Infinity]|>,
  {timeSymbol, {tauA, tauV, tauX}}];

deriveEnergyOperands[anchoring_, face_, tractionOrientation_:-1] := Module[
  {measure, normal, faceVelocityVector, bulkAmplitude,
   boundaryEquation, amplitudeRule, pressureAtFace,
   tractionVector, facePairing, farFieldPressure,
   farFieldVelocity, outgoingFlux, comparisonBulkOperand},
  measure = backgroundMeasureFirstShape[anchoring, face];
  normal = substrateByAnchoringAndFace[
    anchoring <> "|FACE_" <> ToString[face] <>
      "|RHOBR_CONSTANT"][
      "T_A_FACE_NORMAL"]["BACKGROUND"];
  normal = Simplify[(normal /. sigmaW -> geometryScale sigmaW) /.
      geometryScale -> 0] + shapeCoefficient[
    normal /. sigmaW -> geometryScale sigmaW, geometryScale];
  faceVelocityVector = energyVelocity normal;
  boundaryEquation = I qOut[omega, momentumOutput] bulkAmplitude ==
    energyVelocity;
  amplitudeRule = First[Solve[boundaryEquation, bulkAmplitude]];
  pressureAtFace = I rhoM omega bulkAmplitude /. amplitudeRule;
  tractionVector = tractionOrientation pressureAtFace normal;
  facePairing = Simplify[measure Re[
    tractionVector.Conjugate[faceVelocityVector]]/2];
  farFieldPressure = I rhoM omega bulkAmplitude farFieldPhase /.
    amplitudeRule;
  farFieldVelocity = I qOut[omega, momentumOutput]
    bulkAmplitude farFieldPhase /. amplitudeRule;
  outgoingFlux = Simplify[measure Re[
    farFieldPressure Conjugate[farFieldVelocity]]/2];
  comparisonBulkOperand = -outgoingFlux;
  <|"FACE_TRACTION_PAIRING" -> facePairing,
    "OUTGOING_FARFIELD_POYNTING_FLUX" -> outgoingFlux,
    "COMPARISON_ORIENTED_BULK_OPERAND" -> comparisonBulkOperand,
    "RESIDUAL_A_MINUS_B" ->
      Simplify[facePairing - comparisonBulkOperand],
    "SOURCE_BOUNDARY_EQUATION" -> HoldForm[boundaryEquation]|>
];

energyCases = Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face]) ->
    deriveEnergyOperands[anchoring, face],
  {anchoring, anchorings}, {face, faces}];
energySignCorruptedCases = Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face]) ->
    deriveEnergyOperands[anchoring, face, 1],
  {anchoring, anchorings}, {face, faces}];

(* ---------------------------------------------------------------------- *)
(* Dimensions and multigrades derived from the supplied equations.       *)
(* ---------------------------------------------------------------------- *)

dimensionL = {1, 0, 0};
dimensionT = {0, 1, 0};
dimensionM = {0, 0, 1};
dimensionZero = {0, 0, 0};
dimensionGradient = -dimensionL;
dimensionTimeDerivative = -dimensionT;
dimensionVelocity = dimensionL - dimensionT;
dimensionRhoM = dimensionM - 4 dimensionL;
dimensionRhoBr = dimensionM - 3 dimensionL;
dimensionPhi = dimensionVelocity - dimensionGradient;
dimensionPressure = dimensionRhoM + dimensionTimeDerivative +
  dimensionPhi;
dimensionImpedance = dimensionPressure - dimensionVelocity;
dimensionMuTheta = dimensionRhoBr + 2 dimensionVelocity;
dimensionMuS = dimensionMuTheta - dimensionRhoBr;
dimensionFlux = dimensionRhoM + dimensionVelocity;
dimensionLambdaA = dimensionFlux - dimensionMuS;
dimensionLambdaV = dimensionFlux - dimensionVelocity;
dimensionLambdaX = dimensionPressure - dimensionMuS;
dimensionAddedMass = dimensionPressure -
  (dimensionLengthAcceleration = dimensionL - 2 dimensionT);
dimensionTraction = dimensionPressure;

dimensionDerivationEquations = <|
  "VELOCITY_FROM_POTENTIAL" -> relationalObject[
    dimensionVelocity, dimensionPhi + dimensionGradient],
  "PRESSURE_FROM_POTENTIAL" -> relationalObject[
    dimensionPressure,
    dimensionRhoM + dimensionPhi + dimensionTimeDerivative],
  "IMPEDANCE_PRESSURE_PER_VELOCITY" -> relationalObject[
    dimensionImpedance, dimensionPressure - dimensionVelocity],
  "MU_THETA_FROM_ENERGY_VARIATION" -> relationalObject[
    dimensionMuTheta, dimensionRhoBr + 2 dimensionVelocity],
  "MU_S_FROM_MU_THETA_PER_RHO_BR" -> relationalObject[
    dimensionMuS, dimensionMuTheta - dimensionRhoBr],
  "FLUX_FROM_RHO_M_VELOCITY" -> relationalObject[
    dimensionFlux, dimensionRhoM + dimensionVelocity],
  "LAMBDA_A_FROM_FLUX_PER_AFFINITY" -> relationalObject[
    dimensionLambdaA, dimensionFlux - dimensionMuS],
  "LAMBDA_V_FROM_FLUX_PER_VELOCITY" -> relationalObject[
    dimensionLambdaV, dimensionFlux - dimensionVelocity],
  "LAMBDA_X_FROM_TRACTION_PER_AFFINITY" -> relationalObject[
    dimensionLambdaX, dimensionPressure - dimensionMuS],
  "ADDED_MASS_PRESSURE_PER_ACCELERATION" -> relationalObject[
    dimensionAddedMass, dimensionPressure - dimensionLengthAcceleration]|>;

dimensionByObject = <|
  "DTN_FLAT_SYMBOL" -> dimensionImpedance,
  "DTN_OPERATOR" -> dimensionImpedance,
  "DTN_KERNEL" -> dimensionImpedance,
  "FACE_PRESSURE" -> dimensionPressure,
  "RELATIVE_FLUX" -> dimensionFlux,
  "TRACTION" -> dimensionTraction,
  "MU_THETA" -> dimensionMuTheta,
  "MU_S" -> dimensionMuS,
  "LAMBDA_A_0" -> dimensionLambdaA,
  "LAMBDA_V_0" -> dimensionLambdaV,
  "LAMBDA_X_0" -> dimensionLambdaX,
  "TAU_A" -> dimensionT, "TAU_V" -> dimensionT,
  "TAU_X" -> dimensionT, "ADDED_MASS" -> dimensionAddedMass|>;

gradeTerms[expression_] := Module[{expanded, terms, grades},
  expanded = Expand[expression];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  grades = ({Exponent[#, epsilonShape], Exponent[#, etaBg],
      Exponent[#, sigmaW]} & /@ terms);
  DeleteDuplicates[grades]
];

objectRecord[expression_, dimension_] := <|
  "EXPRESSION" -> expression,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> gradeTerms[expression],
  "DIMENSION_L_T_M" -> dimension|>;

mapDifference[left_Association, right_Association] := AssociationMap[
  Function[key, mapDifference[left[key], right[key]]],
  Intersection[Keys[left], Keys[right]]];
mapDifference[left_List, right_List] /; Length[left] == Length[right] :=
  MapThread[mapDifference, {left, right}];
mapDifference[left_, right_] := Together[left - right];

(* ---------------------------------------------------------------------- *)
(* Primary DtN emissions.                                                *)
(* ---------------------------------------------------------------------- *)

dtnOperatorPayload = Association@Table[anchoring -> <|
    "COMPOSITION" -> objectRecord[
      operatorCompositionFromDerivation[anchoring], dimensionImpedance],
    "TWO_DISCONNECTED_FACE_BLOCK" -> DiagonalMatrix[
      {operatorCompositionFromDerivation[anchoring],
       operatorCompositionFromDerivation[anchoring]}],
    "KERNEL_EQUIVALENCE_RESIDUAL" -> Together[
      directMain[anchoring]["KERNEL"] -
        layerMain[anchoring]["KERNEL"]],
    "PROFILE_SOURCE" -> sourceForCase[anchoring, unitSlopeWeights]
    |>, {anchoring, anchorings}];
emitShared["DTN_OPERATOR", dtnOperatorPayload];
Clear[dtnOperatorPayload];

beginAssociationEmission[sharedObject["DTN_KERNEL"]];
firstKernelEmission = True;
Do[
  Do[
    kernelKey = anchoring <> "|FACE_" <> ToString[face];
    appendAssociationEmission[kernelKey, <|
      "OUTPUT_LEG" -> momentumOutput,
      "INPUT_LEG" -> momentumInput,
      "MOMENTUM_TRANSFER" -> momentumTransfer,
      "KERNEL" -> objectRecord[
        directMain[anchoring]["KERNEL"], dimensionImpedance],
      "SOURCE" -> sourceForCase[anchoring, unitSlopeWeights]|>,
      firstKernelEmission];
    firstKernelEmission = False,
    {face, faces}],
  {anchoring, anchorings}];
endAssociationEmission[];
Clear[firstKernelEmission, kernelKey];

rigidOperandPayload = Association@Table[anchoring -> <|
    "CURVED_DIAGONAL_OPERAND_A" ->
      objectRecord[rigidShiftKernel[anchoring], dimensionImpedance],
    "FLAT_TRANSLATED_OPERAND_B" ->
      objectRecord[flatRigidShiftReference, dimensionImpedance],
    "SOURCE" -> rigidShiftSource[anchoring]|>,
  {anchoring, anchorings}];
emitShared["DTN_RIGID_SHIFT_OPERAND", rigidOperandPayload];
emitShared["DTN_RIGID_SHIFT_RESIDUAL", Association@Table[
  anchoring -> objectRecord[
    rigidShiftKernel[anchoring] - flatRigidShiftReference,
    dimensionImpedance], {anchoring, anchorings}]];
Clear[rigidOperandPayload];

beginAssociationEmission[sharedObject["DTN_BY_REGIME_PAIR"]];
firstRegimeEmission = True;
Do[
  regimeKey = anchoring <> "|FACE_" <> ToString[face] <>
    "|OUTPUT_" <> outputRegime <> "|INPUT_" <> inputRegime;
  appendAssociationEmission[regimeKey,
    regimePairRecord[directMain[anchoring]["KERNEL"],
      outputRegime, inputRegime], firstRegimeEmission];
  firstRegimeEmission = False,
  {anchoring, anchorings}, {face, faces}, {outputRegime, regimes},
  {inputRegime, regimes}];
endAssociationEmission[];
Clear[firstRegimeEmission, regimeKey];

emitShared["DTN_BY_PARITY", Association@Table[
  anchoring -> objectRecord[parityOperatorMatrix[
    flatSymbolC1 + directMain[anchoring]["KERNEL"],
    flatSymbolC1 + directMain[anchoring]["KERNEL"]],
    dimensionImpedance], {anchoring, anchorings}]];

emitShared["DTN_HERMITIAN_PART", Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face]) -> objectRecord[
    hermitianKernel[flatSymbolC1 + directMain[anchoring]["KERNEL"],
      anchoring, face], dimensionImpedance],
  {anchoring, anchorings}, {face, faces}]];

emitShared["DTN_REACTIVE_PART", Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face]) -> objectRecord[
    reactiveKernel[flatSymbolC1 + directMain[anchoring]["KERNEL"],
      anchoring, face], dimensionImpedance],
  {anchoring, anchorings}, {face, faces}]];

inertialLoading = Module[{zEvanescent, pressureFromDisplacement,
    outwardAcceleration},
  zEvanescent = flatSymbolC1 /. qOutputLive ->
    qRegimeValue["EVANESCENT", momentumOutput];
  pressureFromDisplacement = Simplify[
    zEvanescent (-I omega outwardDisplacement)];
  outwardAcceleration = -omega^2 outwardDisplacement;
  Simplify[pressureFromDisplacement/outwardAcceleration]
];
emitShared["DTN_INERTIAL_LOADING", Association@Table[
  (anchoring <> "|FACE_" <> ToString[face] <>
    "|PURELY_REACTIVE_BLOCK") -> <|
      "DEFINING_RELATION" -> relationalObject[
        inertialPressure, inertialLoading outwardAcceleration],
      "PER_FACE_OBJECT" -> objectRecord[
        inertialLoading, dimensionAddedMass],
      "OUTWARD_ACCELERATION" -> outwardAcceleration|>,
  {anchoring, anchorings}, {face, faces}]];

emitShared["DTN_GRAZING_BEHAVIOUR", Association@Flatten@Table[
  (anchoring <> "|FACE_" <> ToString[face] <>
    "|OUTPUT_" <> outputRegime <>
    "|INPUT_" <> inputRegime) -> <|
      "BEHAVIOUR_OBJECT" -> grazingExpression[
        directMain[anchoring]["KERNEL"], outputRegime, inputRegime],
      "DIMENSION_L_T_M" -> dimensionImpedance,
      "MULTIGRADE_EPSILON_ETA_SIGMAW" -> gradeTerms[
        grazingExpression[directMain[anchoring]["KERNEL"],
          outputRegime, inputRegime]],
      "COALESCENCE_RELATIONS" -> {
        relationalObject[1/qOutputProbe, 0],
        relationalObject[1/qInputProbe, 0]},
      "STRICT_REST_FRAME_PARAMETER" ->
        relationalObject[vBulkNormal0, 0]|>,
  {anchoring, anchorings}, {face, faces}, {outputRegime, regimes},
  {inputRegime, regimes}]];

emitShared["DTN_TERM_ORIGINS", Association@Table[
  anchoring -> <|
    "ORIGINS" -> Map[objectRecord[#, dimensionImpedance] &,
      directMain[anchoring]["TERM_ORIGINS"]],
    "ORIGIN_SUM_RESIDUAL" -> objectRecord[
      directMain[anchoring]["ORIGIN_SUM_RESIDUAL"],
      dimensionImpedance],
    "SHAPE_SUBSTRATE" -> Association@Flatten@Table[
      "FACE_" <> ToString[face] <> "|" <> density ->
        substrateByAnchoringAndFace[
          anchoring <> "|FACE_" <> ToString[face] <> "|" <> density],
      {face, faces}, {density, densityRepresentatives}]|>,
  {anchoring, anchorings}]];

(* ---------------------------------------------------------------------- *)
(* Face response, Fredholm object, finite loci, and dissipation.          *)
(* ---------------------------------------------------------------------- *)

faceResponsePayload = Map[
  Function[response, <|
    "PRESSURE" -> objectRecord[response["PRESSURE_FLAT"] +
      response["PRESSURE_FIRST_SHAPE_KERNEL"], dimensionPressure],
    "RELATIVE_FLUX" -> objectRecord[
      response["RELATIVE_FLUX_FLAT"] +
      response["RELATIVE_FLUX_FIRST_SHAPE_KERNEL"], dimensionFlux],
    "TRACTION" -> objectRecord[response["TRACTION_FLAT_VECTOR"] +
      response["TRACTION_FIRST_SHAPE_KERNEL_VECTOR"],
      dimensionTraction],
    "SOURCE_EQUATIONS" -> response["SOURCE_EQUATIONS"]|>],
  mainFaceResponses];

faceResponseTotal[response_, "PRESSURE"] :=
  response["PRESSURE_FLAT"] + response["PRESSURE_FIRST_SHAPE_KERNEL"];
faceResponseTotal[response_, "RELATIVE_FLUX"] :=
  response["RELATIVE_FLUX_FLAT"] +
    response["RELATIVE_FLUX_FIRST_SHAPE_KERNEL"];
faceResponseTotal[response_, "TRACTION"] :=
  response["TRACTION_FLAT_VECTOR"] +
    response["TRACTION_FIRST_SHAPE_KERNEL_VECTOR"];
parityFromFaceCoefficients[plus_List, minus_List] := MapThread[
  parityOperatorMatrix, {plus, minus}];
parityFromFaceCoefficients[plus_, minus_] :=
  parityOperatorMatrix[plus, minus];

parityResponseView[anchoring_, density_] := Module[
  {plusResponse, minusResponse},
  plusResponse = mainFaceResponses[
    anchoring <> "|FACE_1|" <> density];
  minusResponse = mainFaceResponses[
    anchoring <> "|FACE_-1|" <> density];
  AssociationMap[Function[fieldName, <|
    "VELOCITY_PARITY_OPERATOR" -> parityFromFaceCoefficients[
      Coefficient[faceResponseTotal[plusResponse, fieldName],
        faceVelocityInput],
      Coefficient[faceResponseTotal[minusResponse, fieldName],
        faceVelocityInput]],
    "MU_THETA_FACE_VECTOR" -> {
      Coefficient[faceResponseTotal[plusResponse, fieldName], muTheta],
      Coefficient[faceResponseTotal[minusResponse, fieldName], muTheta]}
    |>], {"PRESSURE", "RELATIVE_FLUX", "TRACTION"}]
];

faceResponsePayload = Join[faceResponsePayload,
  Association@Flatten@Table[
    ("PARITY|" <> anchoring <> "|" <> density) ->
      parityResponseView[anchoring, density],
    {anchoring, anchorings}, {density, densityRepresentatives}]];
emitShared["FACE_RESPONSE", faceResponsePayload];
Clear[faceResponsePayload];

emitShared["FACE_RESPONSE_COEFFS", Map[
  Function[response, responseCoefficients[response]],
  mainFaceResponses]];

fredholmPayload = Association@Table[anchoring -> <|
  "OPERATOR" -> OperatorSum[IdentityOperator,
    OperatorComposition[ScalarMultiplier[lambdaA/rhoM^2],
      operatorCompositionFromDerivation[anchoring]]],
  "NONINVERTIBILITY_RELATION" -> relationalObject[
    Inactive[Det][OperatorSum[IdentityOperator,
      OperatorComposition[ScalarMultiplier[lambdaA/rhoM^2],
        operatorCompositionFromDerivation[anchoring]]]], 0],
  "FLAT_DIAGONAL_SYMBOL_RELATION" -> flatDegeneracyEquation,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>,
  {anchoring, anchorings}];
emitShared["NONINVERTIBILITY_CONDITION", fredholmPayload];
Clear[fredholmPayload];

emitShared["DEGENERATE_LOCI_EQUATIONS", <|
  "VARIABLES" -> {LambdaA0},
  "EQUATIONS" -> {flatDegeneracyEquation},
  "OPERANDS" -> <|"COEFFICIENT_MATRIX" -> flatCoefficientMatrix,
    "DETERMINANT" -> flatDegeneracyPolynomial|>,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>];

emitShared["DEGENERATE_LOCI_SOLUTION", <|
  "VARIABLES" -> {LambdaA0},
  "DOMAIN" -> UnrestrictedDomain,
  "SOLUTION" -> flatDegeneracySolutions,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>];

emitShared["DEGENERATE_LOCI_IDENTICALLY_SATISFIED", <|
  "STATUS_TOKEN" -> identityStatus,
  "TEST_OBJECT" -> relationalObject[identityNumerator, 0],
  "OPERANDS" -> <|"NUMERATOR" -> identityNumerator,
    "VARIABLES" -> {LambdaA0}|>,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>];

emitShared["DEGENERATE_LOCI_INCONSISTENT", <|
  "STATUS_TOKEN" -> inconsistentStatus,
  "TEST_OBJECT" -> Inactive[Equivalent][
    Inactive[Solve][{flatDegeneracyEquation}, {LambdaA0}],
    EmptySolutionSet],
  "OPERANDS" -> <|"EQUATIONS" -> {flatDegeneracyEquation},
    "SOLUTION" -> flatDegeneracySolutions|>,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>];

emitShared["DEGENERATE_LOCI_REAL_ADMISSIBLE",
  realAdmissibilityPayload];

portCaseForAxes[anchoring_, face_, density_, outputRegime_,
    inputRegime_] := Module[{baseKey, basePort},
  baseKey = anchoring <> "|FACE_" <> ToString[face] <> "|" <> density;
  basePort = portHermitianBaseCases[baseKey];
  Map[grazingExpression[#, outputRegime, inputRegime] &, basePort]
];
portCaseKey[anchoring_, face_, density_, outputRegime_, inputRegime_,
    parity_] := anchoring <> "|FACE_" <> ToString[face] <> "|" <>
  density <> "|OUTPUT_" <> outputRegime <> "|INPUT_" <> inputRegime <>
  "|PARITY_" <> parity;
portCaseGrade[anchoring_, face_, density_] :=
  portCaseGrade[anchoring, face, density] = gradeTerms[
    epsilonShape (flatSymbolC1 + directMain[anchoring]["KERNEL"])
  ];

beginAssociationEmission[sharedObject["PERMEABLE_PORT_HERMITIAN"]];
firstPortEmission = True;
Do[
  portCase = portCaseForAxes[anchoring, face, density,
    outputRegime, inputRegime];
  appendAssociationEmission[portCaseKey[anchoring, face, density,
    outputRegime, inputRegime, parity], <|
      "PORT_MATRIX" -> portCase["PORT_MATRIX"],
      "HERMITIAN_FORM" -> portCase["HERMITIAN_FORM"],
      "POWER_CONJUGATE_INPUT" -> portCase["INPUT_VECTOR"],
      "POWER_CONJUGATE_OUTPUT" -> portCase["OUTPUT_VECTOR"],
      "INPUT_DIMENSIONS_L_T_M" -> {dimensionVelocity, dimensionMuS},
      "OUTPUT_DIMENSIONS_L_T_M" -> {dimensionPressure, dimensionFlux},
      "POWER_DIMENSION_L_T_M" -> dimensionPressure + dimensionVelocity,
      "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
        portCaseGrade[anchoring, face, density],
      "NONDEGENERATE_SUBSPACE_TEST" ->
        principalMinorObjects[portCase["HERMITIAN_FORM"]],
      "ZEROTH_ORDER_NULLSPACE_STATUS" ->
        "NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER"|>, firstPortEmission];
  firstPortEmission = False;
  Clear[portCase];
  ClearSystemCache[],
  {anchoring, anchorings}, {face, faces},
  {density, densityRepresentatives}, {outputRegime, regimes},
  {inputRegime, regimes}, {parity, {"DELTA_W", "ZETA_C"}}];
endAssociationEmission[];

beginAssociationEmission[
  sharedObject["PERMEABLE_DISSIPATION_VS_OMEGA_TAU"]];
firstMemoryEmission = True;
Do[
  portCase = portCaseForAxes[anchoring, face, density,
    outputRegime, inputRegime];
  appendAssociationEmission[portCaseKey[anchoring, face, density,
    outputRegime, inputRegime, parity], <|
      "HERMITIAN_FORM" -> portCase["HERMITIAN_FORM"],
      "INPUT_DIMENSIONS_L_T_M" -> {dimensionVelocity, dimensionMuS},
      "OUTPUT_DIMENSIONS_L_T_M" -> {dimensionPressure, dimensionFlux},
      "POWER_DIMENSION_L_T_M" -> dimensionPressure + dimensionVelocity,
      "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
        portCaseGrade[anchoring, face, density],
      "INDEPENDENT_MEMORY_LIMITS" ->
        memoryLimitPackage[anchoring, face, density,
          outputRegime, inputRegime]|>,
    firstMemoryEmission];
  firstMemoryEmission = False;
  Clear[portCase];
  ClearSystemCache[],
  {anchoring, anchorings}, {face, faces},
  {density, densityRepresentatives}, {outputRegime, regimes},
  {inputRegime, regimes}, {parity, {"DELTA_W", "ZETA_C"}}];
endAssociationEmission[];
Clear[firstPortEmission, firstMemoryEmission, portCase];

emitShared["ENERGY_FACE_TRACTION_OPERAND", Map[
  Function[energyCase, objectRecord[
    energyCase["FACE_TRACTION_PAIRING"],
    dimensionPressure + dimensionVelocity]], energyCases]];

emitShared["ENERGY_BULK_FARFIELD_FLUX_OPERAND", Map[
  Function[energyCase, <|
    "OUTGOING_FLUX" -> objectRecord[
      energyCase["OUTGOING_FARFIELD_POYNTING_FLUX"],
      dimensionPressure + dimensionVelocity],
    "COMPARISON_ORIENTED_OPERAND" -> objectRecord[
      energyCase["COMPARISON_ORIENTED_BULK_OPERAND"],
      dimensionPressure + dimensionVelocity],
    "SOURCE_BOUNDARY_EQUATION" ->
      energyCase["SOURCE_BOUNDARY_EQUATION"]|>], energyCases]];

emitShared["ENERGY_RESIDUAL", Association@KeyValueMap[
  Function[{energyKey, energyCase}, energyKey -> <|
    "OPERAND_A" -> energyCase["FACE_TRACTION_PAIRING"],
    "OPERAND_B" -> energyCase["COMPARISON_ORIENTED_BULK_OPERAND"],
    "RESIDUAL_A_MINUS_B" -> objectRecord[
      energyCase["RESIDUAL_A_MINUS_B"],
      dimensionPressure + dimensionVelocity],
    "TEST_OBJECT" -> relationalObject[
      energyCase["FACE_TRACTION_PAIRING"],
      energyCase["COMPARISON_ORIENTED_BULK_OPERAND"]],
    "TRACTION_SIGN_SOURCE_CORRUPTION" -> 1,
    "CORRUPTED_OPERAND_A" ->
      energySignCorruptedCases[energyKey]["FACE_TRACTION_PAIRING"],
    "UNCHANGED_OPERAND_B" ->
      energySignCorruptedCases[energyKey][
        "COMPARISON_ORIENTED_BULK_OPERAND"],
    "CORRUPTED_RESIDUAL_A_MINUS_B" ->
      energySignCorruptedCases[energyKey]["RESIDUAL_A_MINUS_B"]|>],
  energyCases]];

(* ---------------------------------------------------------------------- *)
(* Representation route comparison.                                     *)
(* ---------------------------------------------------------------------- *)

representationEulerianPayload = Association@Table[
  anchoring -> <|
    "DTN" -> objectRecord[directMain[anchoring]["KERNEL"],
      dimensionImpedance],
    "FACE_RESPONSE" -> Association@Flatten@Table[
      ("FACE_" <> ToString[face] <> "|" <> density) ->
        faceResponseCase[anchoring, face, density,
        directMain[anchoring]["KERNEL"]],
      {face, faces}, {density, densityRepresentatives}]|>,
  {anchoring, anchorings}];

representationHanzawaPayload = Association@Table[
  anchoring -> <|
    "DTN" -> objectRecord[layerMain[anchoring]["KERNEL"],
      dimensionImpedance],
    "FACE_RESPONSE" -> Association@Flatten@Table[
      ("FACE_" <> ToString[face] <> "|" <> density) ->
        faceResponseCase[anchoring, face, density,
        layerMain[anchoring]["KERNEL"]],
      {face, faces}, {density, densityRepresentatives}],
    "RADIATION_PRESERVING_SECOND_ROUTE" ->
      RadiationPreservingLayerPotential[outgoingGreenKernel,
        IdentityMapAtInfinity]|>,
  {anchoring, anchorings}];

emitShared["REP_INVARIANCE_EULERIAN_OPERAND",
  representationEulerianPayload];
emitShared["REP_INVARIANCE_HANZAWA_OPERAND",
  representationHanzawaPayload];
emitShared["REP_INVARIANCE_RESIDUAL", Association@Table[
  anchoring -> <|
    "DTN_OPERAND_A" -> directMain[anchoring]["KERNEL"],
    "DTN_OPERAND_B" -> layerMain[anchoring]["KERNEL"],
    "DTN_RESIDUAL_A_MINUS_B" -> Together[
      directMain[anchoring]["KERNEL"] -
      layerMain[anchoring]["KERNEL"]],
    "FACE_RESPONSE_RESIDUAL_A_MINUS_B" -> Association@Flatten@Table[
      ("FACE_" <> ToString[face] <> "|" <> density) -> mapDifference[
        representationEulerianPayload[anchoring]["FACE_RESPONSE"][
          "FACE_" <> ToString[face] <> "|" <> density],
        representationHanzawaPayload[anchoring]["FACE_RESPONSE"][
          "FACE_" <> ToString[face] <> "|" <> density]],
      {face, faces}, {density, densityRepresentatives}]|>,
  {anchoring, anchorings}]];

Clear[representationEulerianPayload, representationHanzawaPayload];

(* ---------------------------------------------------------------------- *)
(* One-sided upper-face first-jet reversal.                              *)
(* ---------------------------------------------------------------------- *)

normalFourierFromProfileSource[source_Association, face_] := Module[
  {rawNormal},
  rawNormal = Join[-normalScale source["SLOPE_FOURIER"], {face}]/
    Sqrt[1 + normalScale^2
      source["SLOPE_FOURIER"].source["SLOPE_FOURIER"]];
  shapeCoefficient[rawNormal, normalScale]
];

upperReversedSlopeWeights = {-1, 1, 1};
independenceBase = Association@Table[
  anchoring -> <|
    "DTN_UPPER_FACE" -> directCase[anchoring, unitSlopeWeights,
      qOutputLive, qInputLive]["KERNEL"],
    "T_A_UPPER_FACE_NORMAL" -> normalFourierFromProfileSource[
      sourceForCase[anchoring, unitSlopeWeights], 1],
    "FACE_RESPONSE" -> Association@Table[density ->
      faceResponseCase[anchoring, 1, density,
        directCase[anchoring, unitSlopeWeights,
          qOutputLive, qInputLive]["KERNEL"]],
      {density, densityRepresentatives}],
    "SOURCE" -> sourceForCase[anchoring, unitSlopeWeights]|>,
  {anchoring, anchorings}];

independenceCorrupted = Association@Table[
  anchoring -> <|
    "DTN_UPPER_FACE" -> directCase[anchoring,
      upperReversedSlopeWeights, qOutputLive, qInputLive]["KERNEL"],
    "T_A_UPPER_FACE_NORMAL" -> normalFourierFromProfileSource[
      sourceForCase[anchoring, upperReversedSlopeWeights], 1],
    "FACE_RESPONSE" -> Association@Table[density ->
      faceResponseCase[anchoring, 1, density,
        directCase[anchoring, upperReversedSlopeWeights,
          qOutputLive, qInputLive]["KERNEL"]],
      {density, densityRepresentatives}],
    "MUTATED_SOURCE" -> sourceForCase[anchoring,
      upperReversedSlopeWeights]|>,
  {anchoring, anchorings}];

emitShared["CONTROL_INDEPENDENCE_BASE_OPERAND", independenceBase];
emitShared["CONTROL_INDEPENDENCE_CORRUPTED_OPERAND",
  independenceCorrupted];
emitShared["CONTROL_INDEPENDENCE_RESIDUAL",
  mapDifference[independenceBase, independenceCorrupted]];
Clear[independenceBase, independenceCorrupted];

(* ---------------------------------------------------------------------- *)
(* Per-direction first-jet form ablations, rebuilt from the source.       *)
(* ---------------------------------------------------------------------- *)

formWeights[direction_, value_] := ReplacePart[
  unitSlopeWeights, direction -> value];
formResponseFromKernel[anchoring_, density_, kernel_] :=
  Association@Table["FACE_" <> ToString[face] ->
    faceResponseCase[anchoring, face, density, kernel],
    {face, faces}];
formObjectPackage[objectName_, objectValue_, source_] := <|
  "OBJECT" -> objectValue,
  "T_A_FACE_NORMAL_FROM_SOURCE" -> Association@Table[
    "FACE_" <> ToString[face] ->
      normalFourierFromProfileSource[source, face], {face, faces}],
  "SOURCE" -> source|>;

beginAssociationEmission[sharedObject["CONTROL_FORM_BASE_OPERAND"]];
firstFormEmission = True;
Do[
  Do[
    formKey = objectName <> "|" <> anchoring <> "|" <> density <>
      "|DIRECTION_" <> ToString[direction];
    baseKernel = directCase[anchoring, formWeights[direction, 1],
      qOutputLive, qInputLive]["KERNEL"];
    baseObject = If[objectName === "DTN",
      baseKernel,
      formResponseFromKernel[anchoring, density, baseKernel]];
    appendAssociationEmission[formKey,
      formObjectPackage[objectName, baseObject,
        sourceForCase[anchoring, formWeights[direction, 1]]],
      firstFormEmission];
    firstFormEmission = False,
    {objectName, {"DTN", "FACE_RESPONSE"}}],
  {anchoring, anchorings}, {density, densityRepresentatives},
  {direction, directions}];
endAssociationEmission[];

beginAssociationEmission[sharedObject["CONTROL_FORM_ABLATED_OPERAND"]];
firstFormEmission = True;
Do[
  Do[
    formKey = objectName <> "|" <> anchoring <> "|" <> density <>
      "|DIRECTION_" <> ToString[direction];
    ablatedKernel = directCase[anchoring, formWeights[direction, 0],
      qOutputLive, qInputLive]["KERNEL"];
    ablatedObject = If[objectName === "DTN",
      ablatedKernel,
      formResponseFromKernel[anchoring, density, ablatedKernel]];
    appendAssociationEmission[formKey,
      formObjectPackage[objectName, ablatedObject,
        sourceForCase[anchoring, formWeights[direction, 0]]],
      firstFormEmission];
    firstFormEmission = False,
    {objectName, {"DTN", "FACE_RESPONSE"}}],
  {anchoring, anchorings}, {density, densityRepresentatives},
  {direction, directions}];
endAssociationEmission[];

beginAssociationEmission[sharedObject["CONTROL_FORM_RESIDUAL"]];
firstFormEmission = True;
Do[
  Do[
    formKey = objectName <> "|" <> anchoring <> "|" <> density <>
      "|DIRECTION_" <> ToString[direction];
    baseKernel = directCase[anchoring, formWeights[direction, 1],
      qOutputLive, qInputLive]["KERNEL"];
    ablatedKernel = directCase[anchoring, formWeights[direction, 0],
      qOutputLive, qInputLive]["KERNEL"];
    baseObject = If[objectName === "DTN", baseKernel,
      formResponseFromKernel[anchoring, density, baseKernel]];
    ablatedObject = If[objectName === "DTN", ablatedKernel,
      formResponseFromKernel[anchoring, density, ablatedKernel]];
    basePackage = formObjectPackage[objectName, baseObject,
      sourceForCase[anchoring, formWeights[direction, 1]]];
    ablatedPackage = formObjectPackage[objectName, ablatedObject,
      sourceForCase[anchoring, formWeights[direction, 0]]];
    appendAssociationEmission[formKey, <|
      "OPERAND_A" -> basePackage,
      "OPERAND_B" -> ablatedPackage,
      "RESIDUAL_A_MINUS_B" ->
        mapDifference[basePackage, ablatedPackage]|>,
      firstFormEmission];
    firstFormEmission = False,
    {objectName, {"DTN", "FACE_RESPONSE"}}],
  {anchoring, anchorings}, {density, densityRepresentatives},
  {direction, directions}];
endAssociationEmission[];
Clear[firstFormEmission, formKey, baseKernel, ablatedKernel,
  baseObject, ablatedObject, basePackage, ablatedPackage];

(* ---------------------------------------------------------------------- *)
(* Uniform and zero-jet independently reconstructed regressions.          *)
(* ---------------------------------------------------------------------- *)

regimeSpecializeObject[object_Association, outputRegime_, inputRegime_] :=
  Map[regimeSpecializeObject[#, outputRegime, inputRegime] &, object];
regimeSpecializeObject[object_List, outputRegime_, inputRegime_] :=
  (regimeSpecializeObject[#, outputRegime, inputRegime] & /@ object);
regimeSpecializeObject[object_, outputRegime_, inputRegime_] :=
  grazingExpression[object, outputRegime, inputRegime];

uniformC1 = Association@Flatten@Table[
  (objectName <> "|" <> outputRegime <> "|" <> inputRegime <>
    "|PARITY_" <> parity <> "|" <> density) ->
    regimeSpecializeObject[If[objectName === "DTN",
      flatSymbolC1,
      flatFaceResponseSolve[flatSymbolC1, rhoBr]],
      outputRegime, inputRegime],
  {objectName, {"DTN", "FACE_RESPONSE"}},
  {outputRegime, regimes}, {inputRegime, regimes},
  {parity, {"DELTA_W", "ZETA_C"}},
  {density, densityRepresentatives}];

uniformS11b = Association@Flatten@Table[
  (objectName <> "|" <> outputRegime <> "|" <> inputRegime <>
    "|PARITY_" <> parity <> "|" <> density) ->
    regimeSpecializeObject[If[objectName === "DTN",
      flatSymbolS11b,
      deriveFlatS11bFaceResponse[momentumOutput, rhoBr]],
      outputRegime, inputRegime],
  {objectName, {"DTN", "FACE_RESPONSE"}},
  {outputRegime, regimes}, {inputRegime, regimes},
  {parity, {"DELTA_W", "ZETA_C"}},
  {density, densityRepresentatives}];

emitShared["UNIFORM_LIMIT_S11CC1_OPERAND", uniformC1];
emitShared["UNIFORM_LIMIT_S11B_OPERAND", uniformS11b];
emitShared["UNIFORM_LIMIT_RESIDUAL",
  mapDifference[uniformC1, uniformS11b]];
Clear[uniformC1, uniformS11b];

zeroJetSource[anchoring_] := <|
  "HEIGHT_FOURIER" -> etaBg W0 w1Profile/2,
  "SLOPE_FOURIER" -> ConstantArray[0, braneDimension],
  "WIDTH_BINDING" -> relationalObject[
    WBg @@ anchorArguments[anchoring],
    W0 (1 + etaBg w1Profile)]|>;

zeroJetC1 = Association@Table[anchoring ->
  directBoundaryDerivation[zeroJetSource[anchoring],
    qOutputLive, qOutputLive]["FLAT_IMPEDANCE"] +
  directBoundaryDerivation[zeroJetSource[anchoring],
    qOutputLive, qOutputLive]["KERNEL"],
  {anchoring, anchorings}];
zeroJetS11b = Association@Table[anchoring -> flatSymbolS11b,
  {anchoring, anchorings}];

emitShared["ZERO_JET_S11CC1_OPERAND", zeroJetC1];
emitShared["ZERO_JET_S11B_OPERAND", zeroJetS11b];
emitShared["ZERO_JET_RESIDUAL",
  mapDifference[zeroJetC1, zeroJetS11b]];
Clear[zeroJetC1, zeroJetS11b];

(* ---------------------------------------------------------------------- *)
(* Branch-sign and momentum-leg liveness controls.                       *)
(* ---------------------------------------------------------------------- *)

branchObjectFromKernel[anchoring_, kernel_, qOutput_:qOutputLive,
    qInput_:qInputLive] := <|
  "DTN" -> kernel,
  "FACE_RESPONSE" -> Association@Flatten@Table[
    ("FACE_" <> ToString[face] <> "|" <> density) ->
      faceResponseCase[anchoring, face, density, kernel,
        qOutput, qInput],
    {face, faces}, {density, densityRepresentatives}]|>;

branchBaseline = Association@Table[anchoring ->
  Append[branchObjectFromKernel[anchoring,
    directCase[anchoring, unitSlopeWeights,
      qOutputLive, qInputLive]["KERNEL"]],
    "RADIATION_BRANCH_CHECKS" -> {
      radiationBranchRecord[momentumOutput],
      radiationBranchRecord[momentumInput]}],
  {anchoring, anchorings}];

branchSignFlipInput = Association@Table[anchoring ->
  branchObjectFromKernel[anchoring,
    recomputeBranchKernel[anchoring, "SIGNFLIP_INPUT"],
    qOutputLive, -qInputLive],
  {anchoring, anchorings}];
branchSignFlipOutput = Association@Table[anchoring ->
  branchObjectFromKernel[anchoring,
    recomputeBranchKernel[anchoring, "SIGNFLIP_OUTPUT"],
    -qOutputLive, qInputLive],
  {anchoring, anchorings}];
branchMomentumFreezeOutput = Association@Table[anchoring ->
  branchObjectFromKernel[anchoring,
    recomputeBranchKernel[anchoring, "MOMENTUMFREEZE_OUTPUT"],
    qInputLive, qInputLive],
  {anchoring, anchorings}];
branchMomentumFreezeInput = Association@Table[anchoring ->
  branchObjectFromKernel[anchoring,
    recomputeBranchKernel[anchoring, "MOMENTUMFREEZE_INPUT"],
    qOutputLive, qOutputLive],
  {anchoring, anchorings}];

emitShared["CONTROL_BRANCH_BASE_OPERAND", branchBaseline];
emitShared["CONTROL_BRANCH_SIGNFLIP_INPUT_OPERAND",
  branchSignFlipInput];
emitShared["CONTROL_BRANCH_SIGNFLIP_OUTPUT_OPERAND",
  branchSignFlipOutput];
emitShared["CONTROL_BRANCH_MOMENTUMFREEZE_OUTPUT_OPERAND",
  branchMomentumFreezeOutput];
emitShared["CONTROL_BRANCH_MOMENTUMFREEZE_INPUT_OPERAND",
  branchMomentumFreezeInput];
emitShared["CONTROL_BRANCH_RESIDUAL", <|
  "SIGNFLIP_INPUT" ->
    mapDifference[branchBaseline, branchSignFlipInput],
  "SIGNFLIP_OUTPUT" ->
    mapDifference[branchBaseline, branchSignFlipOutput],
  "MOMENTUMFREEZE_OUTPUT" ->
    mapDifference[branchBaseline, branchMomentumFreezeOutput],
  "MOMENTUMFREEZE_INPUT" ->
    mapDifference[branchBaseline, branchMomentumFreezeInput]|>];
Clear[branchBaseline, branchSignFlipInput, branchSignFlipOutput,
  branchMomentumFreezeOutput, branchMomentumFreezeInput];

(* ---------------------------------------------------------------------- *)
(* Dimension inventory and a source-level homogeneity corruption.         *)
(* ---------------------------------------------------------------------- *)

emitShared["DIMENSIONS", <|
  "DERIVATION_EQUATIONS" -> dimensionDerivationEquations,
  "OBJECT_DIMENSIONS" -> dimensionByObject,
  "RATIO_DEFINITIONS" -> <|
    "DTN" -> PressurePerturbation/OutwardBulkNormalVelocity,
    "ADDED_MASS" -> PressurePerturbation/OutwardFaceAcceleration,
    "FACE_RESPONSE_PRESSURE" ->
      PressurePerturbation/{FaceVelocity, muTheta},
    "FACE_RESPONSE_FLUX" -> RelativeMassFlux/{FaceVelocity, muTheta}|>|>];

homogeneityBaseDimensions = <|
  "BULK_RELATION" -> {
    dimensionPressure,
    dimensionImpedance + dimensionVelocity,
    dimensionImpedance + dimensionFlux - dimensionRhoM},
  "FLUX_CLOSURE" -> {
    dimensionFlux,
    dimensionLambdaA + dimensionMuS,
    dimensionLambdaV + dimensionVelocity},
  "TRACTION" -> {
    dimensionTraction, dimensionPressure,
    dimensionLambdaX + dimensionMuS}|>;

homogeneityControlDimensions = ReplacePart[
  homogeneityBaseDimensions,
  {"FLUX_CLOSURE", 2} ->
    (dimensionLambdaA + dimensionL + dimensionMuS)];

homogeneityRelations[dimensionMap_Association] := AssociationMap[
  Function[key, relationalObject[#, First[dimensionMap[key]]] & /@
    dimensionMap[key]], Keys[dimensionMap]];

homogeneityBasePayload = <|
  "SOURCE_DIMENSIONS" -> dimensionByObject,
  "TERM_DIMENSIONS" -> homogeneityBaseDimensions,
  "RELATIONAL_OBJECTS" ->
    homogeneityRelations[homogeneityBaseDimensions]|>;
homogeneityControlPayload = <|
  "SOURCE_DIMENSIONS" -> ReplacePart[dimensionByObject,
    "LAMBDA_A_0" -> (dimensionLambdaA + dimensionL)],
  "CORRUPTED_SOURCE_OBJECT" -> LambdaA0,
  "TERM_DIMENSIONS" -> homogeneityControlDimensions,
  "RELATIONAL_OBJECTS" ->
    homogeneityRelations[homogeneityControlDimensions]|>;

emitShared["HOMOGENEITY_BASE_OPERAND", homogeneityBasePayload];
emitShared["HOMOGENEITY_CONTROL_OPERAND", homogeneityControlPayload];
emitShared["HOMOGENEITY_RESIDUAL", <|
  "OPERAND_A" -> homogeneityBasePayload,
  "OPERAND_B" -> homogeneityControlPayload,
  "RESIDUAL_A_MINUS_B" ->
    mapDifference[homogeneityBaseDimensions,
      homogeneityControlDimensions],
  "TEST_OBJECT" -> relationalObject[
    homogeneityBaseDimensions, homogeneityControlDimensions]|>];

(* One local inventory records the implementation-only heads and its own tag. *)
reservedGlobalNames = ("Global`" <> SymbolName[#] &) /@
  Values[reservedObjects];
localSymbolInventory = Sort[Select[
  Complement[Names["Global`*"], reservedGlobalNames],
  !generatedSymbolNameQ[#] &]];
emitLocal["TAG_NAMES", <|
  "LOCAL_TAG_NAMES" -> Append[localNames,
    "WL_S11CC1_LOCAL_TAG_NAMES"],
  "ENGINE_LOCAL_SYMBOLS" -> localSymbolInventory|>];
