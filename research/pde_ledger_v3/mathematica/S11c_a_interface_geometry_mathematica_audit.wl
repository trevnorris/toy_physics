$HistoryLength = 0;
ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};

(* Public names exist only at this boundary.  All physics below is rebuilt
   from the supplied S11c-a equations. *)
sharedObject[name_String] := HoldComplete[SharedS11caObject, name];
localObject[name_String] := HoldComplete[LocalS11caObject, name];

standardEmissionName[HoldComplete[SharedS11caObject, name_String]] :=
  "WL_S11CA_" <> name;
standardEmissionName[HoldComplete[LocalS11caObject, name_String]] :=
  "WL_S11CA_LOCAL_" <> name;

stripConditional[object_] := Replace[object,
  ConditionalExpression[value_, condition_] :>
    <|"CONDITIONAL_VALUE" -> value,
      "CONDITION_OPERAND" -> HoldForm[condition]|>, {0, Infinity}];

emittedNames = <||>;
localNames = {};
emit[key_, payload_] := Module[{name, rendered, stream, cleanPayload},
  name = standardEmissionName[key];
  If[!StringQ[name], Quit[90]];
  If[KeyExistsQ[emittedNames, name], Quit[91]];
  AssociateTo[emittedNames, name -> 1];
  If[StringStartsQ[name, "WL_S11CA_LOCAL_"], AppendTo[localNames, name]];
  cleanPayload = stripConditional[payload];
  rendered = ToString[cleanPayload, InputForm, PageWidth -> Infinity];
  stream = First[$Output];
  WriteString[stream, name <> ": " <> rendered <> "\n"];
  Flush[stream];
];
emitShared[name_String, payload_] := emit[sharedObject[name], payload];
emitLocal[name_String, payload_] := emit[localObject[name], payload];

relationalObject[left_, right_] := Inactive[Equal][left, right];
greaterObject[left_, right_] := Inactive[Greater][left, right];
zeroForm[equal_Equal] := equal[[1]] - equal[[2]];
zeroForm[expression_] := expression;
shapeDerivative[expression_List] := shapeDerivative /@ expression;
shapeDerivative[expression_Association] :=
  AssociationMap[shapeDerivative, expression];
shapeDerivative[expression_] := Together[Expand[
  D[expression, waveOrder] /. waveOrder -> 0]];
backgroundValue[expression_] := expression /. waveOrder -> 0;
backgroundTruncate[expression_List] := backgroundTruncate /@ expression;
backgroundTruncate[expression_Association] :=
  AssociationMap[backgroundTruncate, expression];
backgroundTruncate[expression_] := Expand[Normal[Series[
  Normal[Series[expression, {etaBg, 0, 1}]], {sigmaW, 0, 1}]]];
computedExpression[expression_] := backgroundTruncate[expression];

dimensionL = {1, 0, 0};
dimensionT = {0, 1, 0};
dimensionM = {0, 0, 1};
dimensionZero = {0, 0, 0};
dimensionVelocity = dimensionL - dimensionT;
dimensionBulkDensity = dimensionM - 4 dimensionL;
dimensionBraneDensity = dimensionM - 3 dimensionL;
dimensionPressure = dimensionM - 2 dimensionL - 2 dimensionT;
dimensionEnergyBrane = dimensionM - dimensionL - 2 dimensionT;
dimensionFlux = dimensionBulkDensity + dimensionVelocity;
dimensionAffinity = dimensionPressure - dimensionBulkDensity;
dimensionEvolution = dimensionBraneDensity - dimensionT;
dimensionVirtualWork = dimensionPressure + dimensionL;

gradeTerms[expression_List] := DeleteDuplicates[Sort[
  Flatten[gradeTerms /@ expression, 1]]];
gradeTerms[expression_Association] := DeleteDuplicates[Sort[
  Flatten[gradeTerms /@ Values[expression], 1]]];
gradeTerms[expression_] := Module[{expanded, terms},
  expanded = Expand[expression];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  DeleteDuplicates[Sort[({1, Exponent[#, etaBg],
          Exponent[#, sigmaW]} & /@ terms)]]
];
objectRecord[expression_, dimension_, epsilonOrder_:1] := Module[
  {computed, grades},
  computed = computedExpression[expression];
  grades = gradeTerms[computed] /. {1, a_, b_} :> {epsilonOrder, a, b};
  <|"EXPRESSION" -> computed, "MULTIGRADE_EPSILON_ETA_SIGMAW" -> grades,
    "DIMENSION_L_T_M" -> dimension|>
];
suppliedRecord[expression_, dimension_] := <|
  "EXPRESSION" -> expression,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> gradeTerms[expression] /.
    {1, a_, b_} :> {Exponent[expression, waveOrder], a, b},
  "DIMENSION_L_T_M" -> dimension|>;
comparisonPayload[expression_, dimension_] := Join[
  objectRecord[expression, dimension],
  <|"TEST_OBJECT" -> relationalObject[computedExpression[expression], 0]|>];

dimensionalSum[dimensions_List] := Module[{first},
  first = First[dimensions];
  <|"TERM_DIMENSIONS" -> dimensions,
    "COMMON_DIMENSION" -> first,
    "RELATIONS" -> (relationalObject[#, first] & /@ Rest[dimensions])|>
];

dimBrane = 3;
spatialCoordinates = {x1, x2, x3};
materialCoordinates = {capitalX1, capitalX2, capitalX3};
spatialDirections = Range[dimBrane];
faceSigns = {1, -1};
branchNames = {"LAB_HELD", "MATERIAL_ADVECTED"};
densityNames = {"RHO4_CONSTANT", "RHOBR_CONSTANT"};
dofNames = {"DELTA_W", "ZETA_C"};

uVector[coordinates_List] := Through[
  {uOne, uTwo, uThree}[Sequence @@ Append[coordinates, time]]];
virtualUVector[coordinates_List] := Through[
  {virtualUOne, virtualUTwo, virtualUThree}[
    Sequence @@ Append[coordinates, time]]];
bulkVelocityZero[coordinates_List, normalCoordinate_] := Through[
  {zeroVelocityOne, zeroVelocityTwo, zeroVelocityThree, zeroVelocityW}[
      Sequence @@ Append[coordinates, {normalCoordinate, time}]]] /.
    {zeroVelocityOne[___] -> 0, zeroVelocityTwo[___] -> 0,
      zeroVelocityThree[___] -> 0, zeroVelocityW[___] -> 0};
bulkVelocityWave[coordinates_List, normalCoordinate_] := Through[
  {bulkVelocityWaveOne, bulkVelocityWaveTwo, bulkVelocityWaveThree,
    bulkVelocityWaveW}[
      Sequence @@ Append[coordinates, {normalCoordinate, time}]]];
bulkCurrentX[coordinates_List, normalCoordinate_] := Through[
  {bulkCurrentOne, bulkCurrentTwo, bulkCurrentThree}[
    Sequence @@ Append[coordinates, {normalCoordinate, time}]]];

zetaSource[coordinates_List, sign_] :=
  dofCenter (zetaCenter @@ Append[coordinates, time]) +
    sign dofWidth (deltaWidth @@ Append[coordinates, time])/2;
virtualZetaSource[coordinates_List, sign_] :=
  virtualDofCenter (virtualZetaCenter @@ Append[coordinates, time]) +
    sign virtualDofWidth (virtualDeltaWidth @@
      Append[coordinates, time])/2;

widthValue = W0 (1 + etaBg w1Jet[0, 0, 0]);
muRValue = muR (1 + etaBg m1Jet[0, 0, 0]);

widthJetRule[orders_List] := With[{total = Total[orders]},
  Apply[Derivative, orders][widthProfile][x1, x2, x3] ->
    (sigmaW/LW^(total - 1)) (w1Jet @@ orders)];
reversedWidthJetRule[orders_List] := With[{total = Total[orders]},
  Apply[Derivative, orders][widthProfileReversedX1][x1, x2, x3] ->
    (sigmaW/LW^(total - 1)) If[orders === {1, 0, 0},
      -(w1Jet @@ orders), w1Jet @@ orders]];
muRJetRule[orders_List] := With[{total = Total[orders]},
  Apply[Derivative, orders][muRProfile][x1, x2, x3] ->
    (muR/W0) (sigmaW/LW^(total - 1)) (m1Jet @@ orders)];

profileJetRules = Join[
  {widthProfile[x1, x2, x3] -> widthValue,
    muRProfile[x1, x2, x3] -> muRValue},
  Flatten[Table[If[1 <= i + j + k <= 3,
      {widthJetRule[{i, j, k}], muRJetRule[{i, j, k}]}, Nothing],
    {i, 0, 3}, {j, 0, 3}, {k, 0, 3}]]];
reversedWidthProfileJetRules = Join[
  {widthProfileReversedX1[x1, x2, x3] -> widthValue},
  Flatten[Table[If[1 <= i + j + k <= 3,
      {reversedWidthJetRule[{i, j, k}]}, Nothing],
    {i, 0, 3}, {j, 0, 3}, {k, 0, 3}]]];
corruptedProfileJetRules = Join[profileJetRules,
  reversedWidthProfileJetRules];

materialToEulerianRules = Thread[materialCoordinates -> spatialCoordinates];

applySelectedProfileJets[expression_, rules_List,
    ablatedDirection_:0] := Module[{result, target},
  result = expression /. materialToEulerianRules /. rules;
  If[MemberQ[spatialDirections, ablatedDirection],
    target = UnitVector[dimBrane, ablatedDirection];
    result = result /. {
      w1Jet @@ target -> 0,
      m1Jet @@ target -> m1Jet @@ target}];
  computedExpression[result]
];
applyProfileJets[expression_, ablatedDirection_:0] :=
  applySelectedProfileJets[expression, profileJetRules, ablatedDirection];
applyCorruptedProfileJets[expression_, ablatedDirection_:0] :=
  applySelectedProfileJets[expression, corruptedProfileJetRules,
    ablatedDirection];

directWidthProfileSource[coordinates_List, reversedFirstJet_:False] :=
  If[TrueQ[reversedFirstJet],
    widthProfileReversedX1 @@ coordinates, widthProfile @@ coordinates];

rho4Profile["RHO4_CONSTANT", coordinates_List] := rhoBr/W0;
rhoBrProfile["RHO4_CONSTANT", coordinates_List] :=
  (rhoBr/W0) (widthProfile @@ coordinates);
rhoBrProfile["RHOBR_CONSTANT", coordinates_List] := rhoBr;
rho4Profile["RHOBR_CONSTANT", coordinates_List] :=
  rhoBr/(widthProfile @@ coordinates);

inverseCoordinates[] := spatialCoordinates - waveOrder uVector[spatialCoordinates];
materialMap[] := materialCoordinates + waveOrder uVector[materialCoordinates];

heightSource[branch_String, sign_, reversedFirstJet_:False] := Module[
  {inverse},
  inverse = inverseCoordinates[];
  Switch[branch,
    "LAB_HELD",
      sign directWidthProfileSource[spatialCoordinates,
        reversedFirstJet]/2 +
        waveOrder zetaSource[inverse, sign],
    "MATERIAL_ADVECTED",
      sign directWidthProfileSource[inverse, reversedFirstJet]/2 +
        waveOrder zetaSource[inverse, sign]]
];

levelSetSource[branch_String, sign_, reversedFirstJet_:False] :=
  normalCoordinate - heightSource[branch, sign, reversedFirstJet];

parametricFaceMapSource[branch_String, sign_, reversedFirstJet_:False] := Module[
  {mapped, vertical},
  mapped = materialMap[];
  vertical = Switch[branch,
    "LAB_HELD", sign directWidthProfileSource[mapped,
      reversedFirstJet]/2 +
      waveOrder zetaSource[materialCoordinates, sign],
    "MATERIAL_ADVECTED", sign directWidthProfileSource[
      materialCoordinates, reversedFirstJet]/2 +
      waveOrder zetaSource[materialCoordinates, sign]];
  Append[mapped, vertical]
];

flatteningCoordinateSource[branch_String] := Module[{mapped, denominator},
  mapped = materialMap[];
  denominator = Switch[branch,
    "LAB_HELD", (widthProfile @@ mapped) +
      waveOrder dofWidth (deltaWidth @@ Append[materialCoordinates, time]),
    "MATERIAL_ADVECTED", (widthProfile @@ materialCoordinates) +
      waveOrder dofWidth (deltaWidth @@ Append[materialCoordinates, time])];
  (normalCoordinate - waveOrder dofCenter (zetaCenter @@
      Append[materialCoordinates, time]))/denominator
];

graphNormalSource[branch_String, sign_, reversedFirstJet_:False] := Module[
  {height, gradient},
  height = heightSource[branch, sign, reversedFirstJet];
  gradient = D[height, #] & /@ spatialCoordinates;
  sign Join[-gradient, {1}]/Sqrt[1 + gradient.gradient]
];

graphMeasureSource[branch_String, sign_, reversedFirstJet_:False] := Module[
  {height, gradient},
  height = heightSource[branch, sign, reversedFirstJet];
  gradient = D[height, #] & /@ spatialCoordinates;
  Sqrt[1 + gradient.gradient]
];

parametricTangentSource[branch_String, sign_] := Module[{faceMap},
  faceMap = Append[materialMap[], faceHeightFromFlattening[branch, sign]];
  Table[D[faceMap, materialCoordinates[[index]]],
    {index, spatialDirections}]
];

parametricNormalSource[branch_String, sign_] := Module[
  {tangents, spatialJacobian, verticalGradient, conormal},
  tangents = parametricTangentSource[branch, sign];
  spatialJacobian = tangents[[All, 1 ;; dimBrane]];
  verticalGradient = tangents[[All, 4]];
  conormal = Join[-Inverse[spatialJacobian].verticalGradient,
    {1}];
  sign conormal/Sqrt[conormal.conormal]
];

parametricMeasureSource[branch_String, sign_] := Module[
  {tangents, spatialJacobian},
  tangents = parametricTangentSource[branch, sign];
  spatialJacobian = tangents[[All, 1 ;; dimBrane]];
  Sqrt[Det[tangents.Transpose[tangents]]]/Det[spatialJacobian]
];

faceVelocitySource[branch_String, sign_, reversedFirstJet_:False] :=
  D[parametricFaceMapSource[branch, sign, reversedFirstJet], time] /.
    materialToEulerianRules;

virtualFaceDisplacementSource[branch_String, sign_,
    reversedFirstJet_:False] := Module[
  {virtualMap, mapped, vertical},
  mapped = materialCoordinates + waveOrder uVector[materialCoordinates] +
    virtualOrder virtualUVector[materialCoordinates];
  vertical = Switch[branch,
    "LAB_HELD", sign directWidthProfileSource[mapped,
      reversedFirstJet]/2 +
      waveOrder zetaSource[materialCoordinates, sign] +
      virtualOrder virtualZetaSource[materialCoordinates, sign],
    "MATERIAL_ADVECTED", sign directWidthProfileSource[
      materialCoordinates, reversedFirstJet]/2 +
      waveOrder zetaSource[materialCoordinates, sign] +
      virtualOrder virtualZetaSource[materialCoordinates, sign]];
  virtualMap = Append[mapped, vertical];
  D[virtualMap, virtualOrder] /. virtualOrder -> 0 /.
    materialToEulerianRules
];

faceHeightFromFlattening[branch_String, sign_] := Module[
  {solutions},
  solutions = Solve[flatteningCoordinateSource[branch] == sign/2,
    normalCoordinate];
  normalCoordinate /. First[solutions]
];

traceSource[fieldZero_, fieldWave_, branch_String, sign_,
    reversedFirstJet_:False] := Module[
  {height},
  height = heightSource[branch, sign, reversedFirstJet];
  fieldZero[spatialCoordinates, height] +
    waveOrder fieldWave[spatialCoordinates, height]
];

pressureZero[coordinates_List, normal_] :=
  0;
pressureWave[coordinates_List, normal_] :=
  pressurePerturbation @@ Append[coordinates, {normal, time}];
rhoBulkZero[coordinates_List, normal_] :=
  rhoBulkBackground @@ Append[coordinates, {normal, time}];
rhoBulkWave[coordinates_List, normal_] :=
  rhoBulkPerturbation @@ Append[coordinates, {normal, time}];
currentWZero[coordinates_List, normal_] :=
  currentWBackground @@ Append[coordinates, {normal, time}];
currentWWave[coordinates_List, normal_] :=
  currentWPerturbation @@ Append[coordinates, {normal, time}];
currentXZero[index_][coordinates_List, normal_] :=
  Symbol["currentXBackground" <> ToString[index]] @@
    Append[coordinates, {normal, time}];
currentXWave[index_][coordinates_List, normal_] :=
  Symbol["currentXPerturbation" <> ToString[index]] @@
    Append[coordinates, {normal, time}];

bulkVelocityTraceSource[branch_String, sign_,
    reversedFirstJet_:False] := Module[{height},
  height = heightSource[branch, sign, reversedFirstJet];
  bulkVelocityZero[spatialCoordinates, height] +
    waveOrder bulkVelocityWave[spatialCoordinates, height]
];

pressureTraceSource[branch_String, sign_, reversedFirstJet_:False] :=
  traceSource[pressureZero, pressureWave, branch, sign, reversedFirstJet];

branchedProfileCoordinates[branch_String] := Switch[branch,
  "LAB_HELD", spatialCoordinates,
  "MATERIAL_ADVECTED", inverseCoordinates[]];

muThetaSource[branch_String] := Switch[branch,
  "LAB_HELD", waveOrder (muThetaOperand @@
    Append[spatialCoordinates, time]),
  "MATERIAL_ADVECTED", waveOrder (muThetaOperand @@
    Append[inverseCoordinates[], time])];

affinitySource[branch_String, sign_, density_String,
    reversedFirstJet_:False] := Module[
  {coordinates, densityProfile},
  coordinates = branchedProfileCoordinates[branch];
  densityProfile = rhoBrProfile[density, coordinates];
  muThetaSource[branch]/densityProfile -
    pressureTraceSource[branch, sign, reversedFirstJet]/rhoM
];

lambdaA = lambdaAZero/(1 - I omega tauA);
lambdaV = lambdaVZero/(1 - I omega tauV);
lambdaX = lambdaXZero/(1 - I omega tauX);

directFaceObjectsSource[branch_String, sign_, density_String,
    reversedFirstJet_:False] := Module[
  {normal, measure, velocity, normalVelocity, bulkVelocity, flux,
    affinity, pressure, traction, kinematicOperandA, kinematicOperandB,
    kinematic, closure},
  normal = graphNormalSource[branch, sign, reversedFirstJet];
  measure = graphMeasureSource[branch, sign, reversedFirstJet];
  velocity = faceVelocitySource[branch, sign, reversedFirstJet];
  normalVelocity = normal.velocity;
  bulkVelocity = bulkVelocityTraceSource[branch, sign, reversedFirstJet];
  flux = rhoM (bulkVelocity - velocity).normal;
  affinity = affinitySource[branch, sign, density, reversedFirstJet];
  pressure = pressureTraceSource[branch, sign, reversedFirstJet];
  traction = -(pressure + lambdaX affinity) normal;
  kinematicOperandA = normal.bulkVelocity;
  kinematicOperandB = normalVelocity + flux/rhoM;
  kinematic = kinematicOperandA - kinematicOperandB;
  closure = flux - lambdaA affinity - lambdaV normalVelocity;
  <|"NORMAL" -> normal, "MEASURE" -> measure,
    "FACE_VELOCITY_VECTOR" -> velocity, "NORMAL_VELOCITY" -> normalVelocity,
    "BULK_VELOCITY_TRACE" -> bulkVelocity, "FLUX" -> flux,
    "AFFINITY" -> affinity, "PRESSURE_TRACE" -> pressure,
    "TRACTION" -> traction, "KINEMATIC" -> kinematic,
    "KINEMATIC_OPERAND_A" -> kinematicOperandA,
    "KINEMATIC_OPERAND_B" -> kinematicOperandB,
    "CLOSURE" -> closure|>
];

materialFaceObjectsSource[branch_String, sign_, density_String] := Module[
  {normal, measure, velocity, normalVelocity, faceMap, height,
    bulkVelocity, pressure, coordinates, densityProfile, muTheta,
    affinity, flux, traction, kinematic, closure},
  normal = parametricNormalSource[branch, sign];
  measure = parametricMeasureSource[branch, sign];
  velocity = faceVelocitySource[branch, sign];
  height = faceHeightFromFlattening[branch, sign];
  faceMap = Append[materialMap[], height];
  bulkVelocity = bulkVelocityZero[Most[faceMap], height] +
    waveOrder bulkVelocityWave[Most[faceMap], height];
  pressure = pressureZero[Most[faceMap], height] +
    waveOrder pressureWave[Most[faceMap], height];
  coordinates = Switch[branch,
    "LAB_HELD", Most[faceMap],
    "MATERIAL_ADVECTED", materialCoordinates];
  densityProfile = rhoBrProfile[density, coordinates];
  muTheta = waveOrder (muThetaOperand @@ Append[coordinates, time]);
  affinity = muTheta/densityProfile - pressure/rhoM;
  normalVelocity = normal.velocity;
  flux = rhoM (bulkVelocity - velocity).normal;
  traction = -(pressure + lambdaX affinity) normal;
  kinematic = normal.bulkVelocity - normalVelocity - flux/rhoM;
  closure = flux - lambdaA affinity - lambdaV normalVelocity;
  (<|"NORMAL" -> normal, "MEASURE" -> measure,
    "FACE_VELOCITY_VECTOR" -> velocity,
    "NORMAL_VELOCITY" -> normalVelocity,
    "BULK_VELOCITY_TRACE" -> bulkVelocity, "FLUX" -> flux,
    "AFFINITY" -> affinity, "PRESSURE_TRACE" -> pressure,
    "TRACTION" -> traction, "KINEMATIC" -> kinematic,
    "CLOSURE" -> closure|>) /. materialToEulerianRules
];

applyPhysicalDof[expression_, dof_String] := Switch[dof,
  "DELTA_W", expression /. {dofWidth -> 1, dofCenter -> 0},
  "ZETA_C", expression /. {dofWidth -> 0, dofCenter -> 1}];
applyVirtualDof[expression_, dof_String] := Switch[dof,
  "DELTA_W", expression /. {virtualDofWidth -> 1,
    virtualDofCenter -> 0},
  "ZETA_C", expression /. {virtualDofWidth -> 0,
    virtualDofCenter -> 1}];
applyDof[expression_, dof_String] := applyVirtualDof[
  applyPhysicalDof[expression, dof], dof];

derivedFaceCase[source_, dof_String, dimension_, ablatedDirection_:0] :=
  objectRecord[
  applyProfileJets[applyDof[shapeDerivative[source], dof],
    ablatedDirection], dimension];

(* The constant bindings and the reserved operand are deliberately inert. *)
inheritedConstants = {cs0, muR, rhoBr, W0, eW, rhoM,
  vBulkNormalZero, muThetaOperand};

faceLabel[1] = "PLUS";
faceLabel[-1] = "MINUS";
faceCaseKey[branch_, sign_, dof_] :=
  branch <> "|FACE_" <> faceLabel[sign] <> "|DOF_" <> dof;
densityFaceCaseKey[branch_, density_, sign_, dof_] :=
  branch <> "|" <> density <> "|FACE_" <> faceLabel[sign] <>
    "|DOF_" <> dof;
densityCaseKey[branch_, density_, dof_] :=
  branch <> "|" <> density <> "|DOF_" <> dof;

directShape[branch_String, sign_, density_String, name_String] :=
  directShape[branch, sign, density, name] =
    shapeDerivative[directFaceObjectsSource[branch, sign, density][name]];
corruptedUpperDirectShape[branch_String, density_String, name_String] :=
  shapeDerivative[directFaceObjectsSource[branch, 1, density, True][name]];
materialShape[branch_String, sign_, density_String, name_String] :=
  materialShape[branch, sign, density, name] =
    shapeDerivative[materialFaceObjectsSource[branch, sign, density][name]];

conormalSource[branch_String, sign_] := Module[
  {height, normal, scalarField, gradient},
  height = heightSource[branch, sign];
  normal = graphNormalSource[branch, sign];
  scalarField = (conormalBackground @@
      Append[spatialCoordinates, {normalCoordinate, time}]) +
    waveOrder (conormalPerturbation @@
      Append[spatialCoordinates, {normalCoordinate, time}]);
  gradient = Join[D[scalarField, #] & /@ spatialCoordinates,
    {D[scalarField, normalCoordinate]}];
  (normal.gradient) /. normalCoordinate -> height
];

traceFieldInventory = <|
  "PRESSURE" -> {pressureZero, pressureWave, dimensionPressure},
  "BULK_DENSITY" -> {rhoBulkZero, rhoBulkWave, dimensionBulkDensity},
  "NORMAL_CURRENT" -> {currentWZero, currentWWave, dimensionFlux},
  "CURRENT_X1" -> {currentXZero[1], currentXWave[1], dimensionFlux},
  "CURRENT_X2" -> {currentXZero[2], currentXWave[2], dimensionFlux},
  "CURRENT_X3" -> {currentXZero[3], currentXWave[3], dimensionFlux},
  "BULK_VELOCITY_X1" -> {
    Function[{coordinates, normal},
      First[bulkVelocityZero[coordinates, normal]]],
    Function[{coordinates, normal},
      First[bulkVelocityWave[coordinates, normal]]], dimensionVelocity},
  "BULK_VELOCITY_X2" -> {
    Function[{coordinates, normal},
      bulkVelocityZero[coordinates, normal][[2]]],
    Function[{coordinates, normal},
      bulkVelocityWave[coordinates, normal][[2]]], dimensionVelocity},
  "BULK_VELOCITY_X3" -> {
    Function[{coordinates, normal},
      bulkVelocityZero[coordinates, normal][[3]]],
    Function[{coordinates, normal},
      bulkVelocityWave[coordinates, normal][[3]]], dimensionVelocity},
  "BULK_VELOCITY_W" -> {
    Function[{coordinates, normal},
      Last[bulkVelocityZero[coordinates, normal]]],
    Function[{coordinates, normal},
      Last[bulkVelocityWave[coordinates, normal]]], dimensionVelocity}
|>;

virtualWorkSource[route_String, branch_String, density_String] := Total[
  Table[With[{objects = If[route === "EULERIAN",
        directFaceObjectsSource[branch, sign, density],
        materialFaceObjectsSource[branch, sign, density]]},
    objects["MEASURE"] objects["TRACTION"].
      virtualFaceDisplacementSource[branch, sign]], {sign, faceSigns}]];

graphThicknessSource[branch_String] :=
  heightSource[branch, 1] - heightSource[branch, -1];

branchedRho4Source[branch_String, density_String] :=
  rho4Profile[density, branchedProfileCoordinates[branch]];

sigmaEulerianSource[branch_String, density_String] :=
  branchedRho4Source[branch, density] *
    (1 + waveOrder (thetaWave @@ Append[spatialCoordinates, time])) *
    graphThicknessSource[branch];

materialJacobianSource[] := Det[Table[
  D[materialMap[][[row]], materialCoordinates[[column]]],
  {row, spatialDirections}, {column, spatialDirections}]];

virtualConstraintDirectSource[branch_String, density_String,
    corrupted_:False] := Module[
  {virtualMap, inverseVirtual, thickness, rhoFour, sigmaEulerian,
    jacobian, evaluationCoordinates},
  virtualMap = materialCoordinates + virtualOrder *
    virtualUVector[materialCoordinates];
  inverseVirtual = spatialCoordinates - virtualOrder *
    virtualUVector[spatialCoordinates];
  evaluationCoordinates = If[TrueQ[corrupted], materialCoordinates,
    virtualMap];
  thickness = Switch[branch,
    "LAB_HELD", (widthProfile @@ evaluationCoordinates) +
      virtualOrder virtualDofWidth (virtualDeltaWidth @@
        Append[materialCoordinates, time]),
    "MATERIAL_ADVECTED", (widthProfile @@
        If[TrueQ[corrupted], inverseVirtual /. Thread[
          spatialCoordinates -> materialCoordinates], materialCoordinates]) +
      virtualOrder virtualDofWidth (virtualDeltaWidth @@
        Append[materialCoordinates, time])];
  rhoFour = Switch[branch,
    "LAB_HELD", rho4Profile[density, evaluationCoordinates],
    "MATERIAL_ADVECTED", rho4Profile[density,
      If[TrueQ[corrupted], inverseVirtual /. Thread[
        spatialCoordinates -> materialCoordinates], materialCoordinates]]];
  sigmaEulerian = rhoFour *
    (1 + virtualOrder (virtualTheta @@ Append[materialCoordinates, time])) *
    thickness;
  jacobian = Det[IdentityMatrix[dimBrane] + virtualOrder Table[
    D[virtualUVector[materialCoordinates][[row]],
      materialCoordinates[[column]]], {row, spatialDirections},
      {column, spatialDirections}]];
  Together[(D[sigmaEulerian jacobian, virtualOrder] /.
      virtualOrder -> 0)/(sigmaEulerian /. virtualOrder -> 0)] /.
    materialToEulerianRules
];

virtualFlatteningCoordinateSource[branch_String] := Module[
  {virtualMap, center, denominator},
  virtualMap = materialCoordinates + virtualOrder *
    virtualUVector[materialCoordinates];
  center = virtualOrder * virtualDofCenter *
    (virtualZetaCenter @@ Append[materialCoordinates, time]);
  denominator = Switch[branch,
    "LAB_HELD", (widthProfile @@ virtualMap) +
      virtualOrder * virtualDofWidth *
        (virtualDeltaWidth @@ Append[materialCoordinates, time]),
    "MATERIAL_ADVECTED", (widthProfile @@ materialCoordinates) +
      virtualOrder * virtualDofWidth *
        (virtualDeltaWidth @@ Append[materialCoordinates, time])];
  (normalCoordinate - center)/denominator
];

virtualFlatteningMapSource[branch_String] := Module[
  {solutions, physicalNormal},
  solutions = Solve[
    virtualFlatteningCoordinateSource[branch] == flattenedNormalCoordinate,
    normalCoordinate];
  physicalNormal = normalCoordinate /. First[solutions];
  Append[materialCoordinates + virtualOrder *
    virtualUVector[materialCoordinates], physicalNormal]
];

virtualConstraintMaterialSource[branch_String, density_String] := Module[
  {flatteningMap, sourceCoordinates, fullJacobian, profileCoordinates,
    rhoFour, pulledFourDensity, pulledMass},
  flatteningMap = virtualFlatteningMapSource[branch];
  sourceCoordinates = Append[materialCoordinates, flattenedNormalCoordinate];
  fullJacobian = Det[Table[
    D[flatteningMap[[row]], sourceCoordinates[[column]]],
    {row, Range[dimBrane + 1]}, {column, Range[dimBrane + 1]}]];
  profileCoordinates = Switch[branch,
    "LAB_HELD", Most[flatteningMap],
    "MATERIAL_ADVECTED", materialCoordinates];
  rhoFour = Switch[branch,
    "LAB_HELD", rho4Profile[density, profileCoordinates],
    "MATERIAL_ADVECTED", rho4Profile[density, materialCoordinates]];
  pulledFourDensity = rhoFour *
    (1 + virtualOrder *
      (virtualTheta @@ Append[materialCoordinates, time])) *
    fullJacobian;
  pulledMass = Integrate[pulledFourDensity,
    {flattenedNormalCoordinate, -1/2, 1/2}];
  Together[(D[pulledMass, virtualOrder] /. virtualOrder -> 0)/
    (pulledMass /. virtualOrder -> 0)] /. materialToEulerianRules
];

evolutionTermsDirectSource[branch_String, density_String,
    corrupted_:False] := Module[{sigma, velocity, sources, transport},
  sigma = sigmaEulerianSource[branch, density];
  velocity = waveOrder D[uVector[spatialCoordinates], time];
  transport = If[TrueQ[corrupted],
    sigma Total[Table[D[velocity[[index]], spatialCoordinates[[index]]],
      {index, spatialDirections}]],
    Total[Table[D[sigma velocity[[index]], spatialCoordinates[[index]]],
      {index, spatialDirections}]]];
  sources = Total[Table[
    directFaceObjectsSource[branch, sign, density]["MEASURE"]
      directFaceObjectsSource[branch, sign, density]["FLUX"],
    {sign, faceSigns}]];
  <|"TIME_DERIVATIVE" -> D[sigma, time],
    "INPLANE_TRANSPORT" -> transport,
    "TRUE_AREA_FACE_FLUX" -> sources|>
];

evolutionTermsMaterialSource[branch_String, density_String] := Module[
  {mapped, sigmaAtMap, jacobian, sources},
  mapped = materialMap[];
  sigmaAtMap = sigmaEulerianSource[branch, density] /.
    Thread[spatialCoordinates -> mapped];
  jacobian = materialJacobianSource[];
  sources = Total[Table[
    materialFaceObjectsSource[branch, sign, density]["MEASURE"]
      materialFaceObjectsSource[branch, sign, density]["FLUX"],
    {sign, faceSigns}]];
  (<|"PULLED_TIME_DERIVATIVE" -> D[sigmaAtMap jacobian, time]/jacobian,
    "PULLED_TRUE_AREA_FACE_FLUX" -> sources|>) /.
      materialToEulerianRules
];

windowSource[branch_String] := Module[{gplus, gminus},
  gplus = levelSetSource[branch, 1];
  gminus = -levelSetSource[branch, -1];
  windowFunction[gplus, gminus]
];

staticFlatWindow = windowFunction[
  normalCoordinate - W0/2, -normalCoordinate - W0/2];

projectionTermsSource[branch_String, dynamic_:True] := Module[
  {window, rho, currentX, currentW, leftTime, leftSpace,
    windowTime, windowSpace, windowNormal},
  window = If[TrueQ[dynamic], windowSource[branch], staticFlatWindow];
  rho = rhoBulkZero[spatialCoordinates, normalCoordinate] +
    waveOrder rhoBulkWave[spatialCoordinates, normalCoordinate];
  currentX = Table[
    currentXZero[index][spatialCoordinates, normalCoordinate] +
      waveOrder currentXWave[index][spatialCoordinates, normalCoordinate],
    {index, spatialDirections}];
  currentW = currentWZero[spatialCoordinates, normalCoordinate] +
    waveOrder currentWWave[spatialCoordinates, normalCoordinate];
  leftTime = D[Inactive[Integrate][window rho,
    {normalCoordinate, -Infinity, Infinity}], time];
  leftSpace = Total[Table[D[Inactive[Integrate][
    window currentX[[index]], {normalCoordinate, -Infinity, Infinity}],
    spatialCoordinates[[index]]], {index, spatialDirections}]];
  windowTime = Inactive[Integrate][D[window, time] rho,
    {normalCoordinate, -Infinity, Infinity}];
  windowSpace = Inactive[Integrate][Total[Table[
    D[window, spatialCoordinates[[index]]] currentX[[index]],
    {index, spatialDirections}]], {normalCoordinate, -Infinity, Infinity}];
  windowNormal = Inactive[Integrate][D[window, normalCoordinate] currentW,
    {normalCoordinate, -Infinity, Infinity}];
  <|"TIME_DERIVATIVE_OF_PROJECTED_DENSITY" -> leftTime,
    "DIVERGENCE_OF_PROJECTED_CURRENT" -> leftSpace,
    "WINDOW_TIME_ORIGIN" -> -windowTime,
    "WINDOW_SPACE_ORIGIN" -> -windowSpace,
    "WINDOW_NORMAL_ORIGIN" -> -windowNormal|>
];

projectionOperand[terms_Association] := Total[Values[terms]];

(* ---------------------------------------------------------------------- *)
(* Supplied background maps, face maps, state, and support premise.       *)
(* ---------------------------------------------------------------------- *)

widthAnsatz = W0 (1 + etaBg w1Profile[y1/LW, y2/LW, y3/LW]);
muRAnsatz = muR (1 + etaBg m1Profile[y1/LW, y2/LW, y3/LW]);
widthGradientAnsatz = (D[widthAnsatz, #] & /@ {y1, y2, y3}) /.
  etaBg W0/LW -> sigmaW;
muRGradientAnsatz = (D[muRAnsatz, #] & /@ {y1, y2, y3}) /.
  etaBg muR/LW -> (muR/W0) sigmaW;

backgroundDensityMap = Association[Table[
  density -> Module[{rhoFour, rhoIntegrated, sigmaZero, gradient},
    rhoFour = Switch[density,
      "RHO4_CONSTANT", rhoBr/W0,
      "RHOBR_CONSTANT", rhoBr/widthAnsatz];
    rhoIntegrated = Switch[density,
      "RHO4_CONSTANT", (rhoBr/W0) widthAnsatz,
      "RHOBR_CONSTANT", rhoBr];
    sigmaZero = Together[rhoFour widthAnsatz];
    gradient = (D[sigmaZero, #] & /@ {y1, y2, y3}) /.
      etaBg W0/LW -> sigmaW;
    <|
      "RHO4_BG" -> objectRecord[rhoFour, dimensionBulkDensity, 0],
      "RHOBR_BG" -> objectRecord[rhoIntegrated,
        dimensionBraneDensity, 0],
      "SIGMA_E_ZERO" -> objectRecord[sigmaZero,
        dimensionBraneDensity, 0],
      "GRADIENT_SIGMA_E_ZERO" -> objectRecord[gradient,
        dimensionBraneDensity - dimensionL, 0]|>
  ], {density, densityNames}]];
emitShared["BACKGROUND_DENSITY_MAP", backgroundDensityMap];

backgroundStatePayload = Association[Flatten[Table[
  With[{key = branch <> "|" <> density}, key -> <|
    "EXPRESSION" -> <|
      "W_BG" -> widthAnsatz,
      "MU_R_BG" -> muRAnsatz,
      "RHO4_BG" -> Switch[density, "RHO4_CONSTANT", rhoBr/W0,
        "RHOBR_CONSTANT", rhoBr/widthAnsatz],
      "RHOBR_BG" -> Switch[density,
        "RHO4_CONSTANT", (rhoBr/W0) widthAnsatz,
        "RHOBR_CONSTANT", rhoBr],
      "THETA_ZERO" -> 0,
      "FACE_VELOCITY_ZERO" -> AssociationThread[faceSigns, {0, 0}],
      "FACE_FLUX_ZERO" -> AssociationThread[faceSigns, {0, 0}],
      "AFFINITY_ZERO" -> AssociationThread[faceSigns, {0, 0}],
      "BOUNDARY_LOADS" -> AssociationThread[faceSigns,
        {tractionHoldPlusZero, tractionHoldMinusZero}]|>,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
      (gradeTerms[widthAnsatz + muRAnsatz] /.
        {1, a_, b_} :> {0, a, b}),
    "DIMENSION_L_T_M" -> <|"W_BG" -> dimensionL,
      "MU_R_BG" -> dimensionEnergyBrane,
      "RHO4_BG" -> dimensionBulkDensity,
      "RHOBR_BG" -> dimensionBraneDensity,
      "THETA_ZERO" -> dimensionZero,
      "FACE_VELOCITY_ZERO" -> dimensionVelocity,
      "FACE_FLUX_ZERO" -> dimensionFlux,
      "AFFINITY_ZERO" -> dimensionAffinity,
      "BOUNDARY_LOADS" -> dimensionPressure|>|>],
  {branch, branchNames}, {density, densityNames}], 1]];
emitShared["BACKGROUND_STATE", backgroundStatePayload];

admissibilityPremisePayload = Association[Table[
  branch -> <|
    "EXPRESSION" -> <|"SUPPORT_BUNDLE" ->
      {forceHoldZero[x1, x2, x3],
        AssociationThread[faceSigns,
          {tractionHoldPlusZero, tractionHoldMinusZero}]},
      "SUPPORT_STABILISED_BACKGROUND" -> HoldForm[
        backgroundState[branch] supportedBy supportHoldZero]|>,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
    "DIMENSION_L_T_M" -> <|"BODY_SUPPORT" ->
      dimensionEnergyBrane - dimensionL,
      "FACE_SUPPORT" -> dimensionPressure|>|>,
  {branch, branchNames}]];
emitShared["ADMISSIBILITY_PREMISE", admissibilityPremisePayload];

faceMapLabPayload = Association[Table[
  faceLabel[sign] -> <|
    "EXPRESSION" -> Inactive[Equal][
      faceMapLab[sign][materialCoordinates, time],
      parametricFaceMapSource["LAB_HELD", sign] /.
        {dofWidth -> 1, dofCenter -> 1}],
    "ZETA_DECOMPOSITION" -> Inactive[Equal][zetaFace[sign],
      zetaCenter + sign deltaWidth/2],
    "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
      gradeTerms[applyProfileJets[shapeDerivative[
        parametricFaceMapSource["LAB_HELD", sign]]]],
    "DIMENSION_L_T_M" -> dimensionL|>, {sign, faceSigns}]];
emitShared["FACE_MAP_LAB_HELD", faceMapLabPayload];

faceMapMaterialPayload = Association[Table[
  faceLabel[sign] -> <|
    "EXPRESSION" -> Inactive[Equal][
      faceMapMaterialAdvected[sign][materialCoordinates, time],
      parametricFaceMapSource["MATERIAL_ADVECTED", sign] /.
        {dofWidth -> 1, dofCenter -> 1}],
    "ZETA_DECOMPOSITION" -> Inactive[Equal][zetaFace[sign],
      zetaCenter + sign deltaWidth/2],
    "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
      gradeTerms[applyProfileJets[shapeDerivative[
        parametricFaceMapSource["MATERIAL_ADVECTED", sign]]]],
    "DIMENSION_L_T_M" -> dimensionL|>, {sign, faceSigns}]];
emitShared["FACE_MAP_MATERIAL_ADVECTED", faceMapMaterialPayload];

(* ---------------------------------------------------------------------- *)
(* T-a through T-e: graph geometry and the bound face laws.               *)
(* ---------------------------------------------------------------------- *)

normalPayload = Association[Flatten[Table[
  With[{key = faceCaseKey[branch, sign, dof],
      exact = graphNormalSource[branch, sign]}, key -> <|
    "EXACT_SOURCE" -> exact,
    "BACKGROUND" -> objectRecord[applyProfileJets[
      backgroundValue[exact]], dimensionZero, 0],
    "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof, dimensionZero],
    "ORIENTATION_OBJECT" -> greaterObject[
      sign exact[[4]], 0],
    "DIMENSION_L_T_M" -> dimensionZero|>],
  {branch, branchNames}, {sign, faceSigns}, {dof, dofNames}]]];
emitShared["FACE_NORMAL", normalPayload];

conormalPayload = Association[Flatten[Table[
  With[{key = faceCaseKey[branch, sign, dof],
      exact = conormalSource[branch, sign]}, key -> <|
    "GRAPH_EVALUATED_SOURCE" -> exact,
    "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof, -dimensionL],
    "OPERATOR_DIMENSION_L_T_M" -> -dimensionL|>],
  {branch, branchNames}, {sign, faceSigns}, {dof, dofNames}]]];
emitShared["CONORMAL_DERIV", conormalPayload];

measurePayload = Association[Flatten[Table[
  With[{key = faceCaseKey[branch, sign, dof],
      exact = graphMeasureSource[branch, sign]}, key -> <|
    "EXACT_SOURCE" -> exact,
    "BACKGROUND" -> objectRecord[applyProfileJets[
      backgroundValue[exact]], dimensionZero, 0],
    "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof, dimensionZero],
    "DIMENSION_L_T_M" -> dimensionZero|>],
  {branch, branchNames}, {sign, faceSigns}, {dof, dofNames}]]];
emitShared["FACE_MEASURE_SHAPE_DERIV", measurePayload];

velocityPayload = Association[Flatten[Table[
  With[{key = faceCaseKey[branch, sign, dof],
      exact = directFaceObjectsSource[branch, sign,
        "RHO4_CONSTANT"]["NORMAL_VELOCITY"]}, key -> <|
    "EXACT_SOURCE" -> exact,
    "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof,
      dimensionVelocity],
    "FACE_VECTOR_SOURCE" -> faceVelocitySource[branch, sign],
    "DIMENSION_L_T_M" -> dimensionVelocity|>],
  {branch, branchNames}, {sign, faceSigns}, {dof, dofNames}]]];
emitShared["FACE_VELOCITY", velocityPayload];

fluxPayload = Association[Flatten[Table[
  With[{key = densityFaceCaseKey[branch, density, sign, dof],
      exact = directFaceObjectsSource[branch, sign, density]["FLUX"]},
    key -> <|"EXACT_RELATIVE_FLUX_SOURCE" -> exact,
      "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof,
        dimensionFlux],
      "TRUE_AREA_DENSITY" -> objectRecord[
        applyProfileJets[applyPhysicalDof[
          directFaceObjectsSource[branch, sign, density]["MEASURE"], dof]],
        dimensionZero, 0],
      "DIMENSION_L_T_M" -> dimensionFlux|>],
  {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {dof, dofNames}]]];
emitShared["RELATIVE_FLUX", fluxPayload];

kinematicPayload = Association[Flatten[Table[
  With[{key = densityFaceCaseKey[branch, density, sign, dof],
      operandASource = directFaceObjectsSource[branch, sign, density][
        "KINEMATIC_OPERAND_A"],
      operandBSource = directFaceObjectsSource[branch, sign, density][
        "KINEMATIC_OPERAND_B"]}, Module[{operandA, operandB, residual},
    operandA = applyProfileJets[applyPhysicalDof[
      shapeDerivative[operandASource], dof]];
    operandB = applyProfileJets[applyPhysicalDof[
      shapeDerivative[operandBSource], dof]];
    residual = Together[operandA - operandB];
    key -> <|"BOUND_SOURCE_LAW" ->
        relationalObject[operandASource, operandBSource],
      "OPERAND_A_SHAPE_DERIVATIVE" -> objectRecord[operandA,
        dimensionVelocity],
      "OPERAND_B_SHAPE_DERIVATIVE" -> objectRecord[operandB,
        dimensionVelocity],
      "RESIDUAL_A_MINUS_B" -> objectRecord[residual, dimensionVelocity],
      "DIMENSION_L_T_M" -> dimensionVelocity|>]],
  {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {dof, dofNames}]]];
emitShared["KINEMATIC_BALANCE", kinematicPayload];

tractionPayload = Association[Flatten[Table[
  With[{key = densityFaceCaseKey[branch, density, sign, dof],
      exact = directFaceObjectsSource[branch, sign, density]["TRACTION"]},
    key -> <|"EXACT_SOURCE" -> exact,
      "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof,
        dimensionPressure],
      "DIMENSION_L_T_M" -> dimensionPressure|>],
  {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {dof, dofNames}]]];
emitShared["TRACTION", tractionPayload];

virtualWorkPayload = Association[Flatten[Table[
  With[{key = densityCaseKey[branch, density, physicalDof] <>
      "|VIRTUAL_DOF_" <> virtualDof,
      exact = virtualWorkSource["EULERIAN", branch, density]}, key -> <|
    "EXACT_TRUE_AREA_SOURCE" -> exact,
    "SHAPE_DERIVATIVE" -> objectRecord[applyProfileJets[
      applyVirtualDof[applyPhysicalDof[shapeDerivative[exact],
        physicalDof], virtualDof]], dimensionVirtualWork],
    "FACE_DISPLACEMENTS_FROM_MAP" -> Association[Table[
      faceLabel[sign] -> applyVirtualDof[
        virtualFaceDisplacementSource[branch, sign], virtualDof],
      {sign, faceSigns}]],
    "DIMENSION_L_T_M" -> dimensionVirtualWork|>],
  {branch, branchNames}, {density, densityNames},
  {physicalDof, dofNames}, {virtualDof, dofNames}]]];
emitShared["VIRTUAL_WORK_SHAPE_DERIV", virtualWorkPayload];

faceShiftPayload = Association[Flatten[Table[
  With[{inventory = traceFieldInventory[fieldName],
      key = branch <> "|FACE_" <> faceLabel[sign] <> "|FIELD_" <>
        fieldName <> "|DOF_" <> dof},
    key -> <|"EXACT_TRACE_SOURCE" -> traceSource[inventory[[1]],
        inventory[[2]], branch, sign],
      "SHAPE_DERIVATIVE" -> derivedFaceCase[
        traceSource[inventory[[1]], inventory[[2]], branch, sign], dof,
        inventory[[3]]],
      "SHIFT_LAW_OPERANDS" -> HoldForm[
        deltaTrace[fieldName, sign] == deltaField[fieldName, sign] +
          deltaHeight[branch, sign]
            normalDerivativeBackground[fieldName, sign]],
      "DIMENSION_L_T_M" -> inventory[[3]]|>],
  {branch, branchNames}, {sign, faceSigns},
  {fieldName, Keys[traceFieldInventory]}, {dof, dofNames}]]];
emitShared["FACE_SHIFT", faceShiftPayload];

(* ---------------------------------------------------------------------- *)
(* T-f: dynamic-window projection and its static-flat route.              *)
(* ---------------------------------------------------------------------- *)

projectionDynamicTerms = AssociationMap[
  projectionTermsSource[#, True] &, branchNames];
projectionStaticTerms = AssociationMap[
  projectionTermsSource[#, False] &, branchNames];

projectionShapePayload = Association[Flatten[Table[
  With[{key = branch <> "|DOF_" <> dof,
      exact = projectionOperand[projectionDynamicTerms[branch]]}, key -> <|
    "DYNAMIC_WINDOW" -> windowSource[branch],
    "EXACT_PROJECTED_IDENTITY_ZERO_FORM" -> exact,
    "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof,
      dimensionEvolution],
    "DIMENSION_L_T_M" -> dimensionEvolution|>],
  {branch, branchNames}, {dof, dofNames}]]];
emitShared["PROJECTION_SHAPE_DERIV", projectionShapePayload];

projectionDynamicPayload = Association[Flatten[Table[
  With[{key = branch <> "|DOF_" <> dof,
      operand = applyProfileJets[applyPhysicalDof[
        shapeDerivative[projectionOperand[projectionDynamicTerms[branch]]],
        dof]]}, key -> objectRecord[operand, dimensionEvolution]],
  {branch, branchNames}, {dof, dofNames}]]];
emitShared["PROJECTION_DYNAMIC_OPERAND", projectionDynamicPayload];

projectionStaticPayload = Association[Flatten[Table[
  With[{key = branch <> "|DOF_" <> dof,
      operand = applyProfileJets[applyPhysicalDof[
        shapeDerivative[projectionOperand[projectionStaticTerms[branch]]],
        dof]]}, key -> objectRecord[operand, dimensionEvolution]],
  {branch, branchNames}, {dof, dofNames}]]];
emitShared["PROJECTION_STATIC_OPERAND", projectionStaticPayload];

projectionResidualPayload = AssociationMap[
  Function[key, Module[{dynamic, static, residual},
    dynamic = projectionDynamicPayload[key]["EXPRESSION"];
    static = projectionStaticPayload[key]["EXPRESSION"];
    residual = Together[dynamic - static];
    comparisonPayload[residual, dimensionEvolution]
  ]], Keys[projectionDynamicPayload]];
emitShared["PROJECTION_RESIDUAL", projectionResidualPayload];

projectionOriginsPayload = Association[Flatten[Table[
  With[{key = branch <> "|" <> origin <> "|DOF_" <> dof,
      dynamicTerm = projectionDynamicTerms[branch][origin],
      staticTerm = projectionStaticTerms[branch][origin]}, key -> <|
    "DYNAMIC_SOURCE_TERM" -> dynamicTerm,
    "DYNAMIC_SHAPE_DERIVATIVE" -> derivedFaceCase[dynamicTerm, dof,
      dimensionEvolution],
    "STATIC_SOURCE_TERM" -> staticTerm,
    "STATIC_SHAPE_DERIVATIVE" -> derivedFaceCase[staticTerm, dof,
      dimensionEvolution],
    "DIMENSION_L_T_M" -> dimensionEvolution|>],
  {branch, branchNames}, {origin, Keys[projectionDynamicTerms[branch]]},
  {dof, dofNames}]]];
emitShared["PROJECTION_TERM_ORIGINS", projectionOriginsPayload];

(* ---------------------------------------------------------------------- *)
(* T-g through T-i: material constraint, evolution, and face closure.     *)
(* ---------------------------------------------------------------------- *)

virtualConstraintPayload = Association[Flatten[Table[
  With[{key = densityCaseKey[branch, density, dof],
      expression = virtualConstraintDirectSource[branch, density]},
    key -> <|
      "NORMALIZED_VIRTUAL_MASS_VARIATION" -> objectRecord[
        applyProfileJets[applyVirtualDof[expression, dof]],
        dimensionZero],
      "EXACT_EW_BG_MAP" -> suppliedRecord[
        W0 eW/widthProfile[x1, x2, x3], dimensionZero],
      "VIRTUAL_DISPLACEMENT_OPERAND" ->
        virtualUVector[spatialCoordinates],
      "DIMENSION_L_T_M" -> dimensionZero|>],
  {branch, branchNames}, {density, densityNames},
  {dof, dofNames}]]];
emitShared["VIRTUAL_CONSTRAINT", virtualConstraintPayload];

evolutionPayload = Association[Flatten[Table[
  With[{key = densityCaseKey[branch, density, dof],
      terms = evolutionTermsDirectSource[branch, density]},
    key -> <|
      "EXACT_TRUE_AREA_BALANCE_ZERO_FORM" -> Total[Values[terms]],
      "SHAPE_DERIVATIVE" -> objectRecord[applyProfileJets[
        applyPhysicalDof[shapeDerivative[Total[Values[terms]]], dof]],
        dimensionEvolution],
      "SLAB_VELOCITY" -> waveOrder D[uVector[spatialCoordinates], time],
      "DIMENSION_L_T_M" -> dimensionEvolution|>],
  {branch, branchNames}, {density, densityNames},
  {dof, dofNames}]]];
emitShared["EVOLUTION_MASS_BALANCE", evolutionPayload];

evolutionOriginsPayload = Association[Flatten[Table[
  With[{terms = evolutionTermsDirectSource[branch, density],
      key = densityCaseKey[branch, density, dof] <> "|ORIGIN_" <> origin},
    key -> <|"EXACT_SOURCE_TERM" -> terms[origin],
      "SHAPE_DERIVATIVE" -> objectRecord[applyProfileJets[
        applyPhysicalDof[shapeDerivative[terms[origin]], dof]],
        dimensionEvolution],
      "DIMENSION_L_T_M" -> dimensionEvolution|>],
  {branch, branchNames}, {density, densityNames},
  {dof, dofNames},
  {origin, Keys[evolutionTermsDirectSource[branch, density]]}]]];
emitShared["EVOLUTION_TERM_ORIGINS", evolutionOriginsPayload];

closurePayload = Association[Flatten[Table[
  With[{key = densityFaceCaseKey[branch, density, sign, dof],
      exact = directFaceObjectsSource[branch, sign, density]["CLOSURE"]},
    key -> <|"BOUND_CLOSURE_ZERO_FORM" -> exact,
      "MU_THETA_RESERVED_OPERAND" -> muThetaSource[branch],
      "MU_S_BRANCHED_MAP" ->
        muThetaSource[branch]/rhoBrProfile[density,
          branchedProfileCoordinates[branch]],
      "SHAPE_DERIVATIVE" -> derivedFaceCase[exact, dof,
        dimensionFlux],
      "TEST_OBJECT" -> relationalObject[
        applyProfileJets[applyPhysicalDof[shapeDerivative[exact], dof]], 0],
      "DIMENSION_L_T_M" -> dimensionFlux|>],
  {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {dof, dofNames}]]];
emitShared["CLOSURE_SHAPE_DERIV", closurePayload];

(* ---------------------------------------------------------------------- *)
(* Representation routes and their one-sided source mutations.            *)
(* ---------------------------------------------------------------------- *)

processedVirtualWorkShape[route_String, branch_String, density_String,
    physicalDof_String, virtualDof_String, ablatedDirection_:0,
    reverseUpper_:False] := Total[Table[Module[
  {objects, faceTerm, processed, reverseDirectSource},
  reverseDirectSource = route === "EULERIAN" && TrueQ[reverseUpper] &&
    sign == 1;
  objects = If[route === "EULERIAN",
    directFaceObjectsSource[branch, sign, density, reverseDirectSource],
    materialFaceObjectsSource[branch, sign, density]];
  faceTerm = objects["MEASURE"] objects["TRACTION"].
    virtualFaceDisplacementSource[branch, sign, reverseDirectSource];
  processed = applyVirtualDof[applyPhysicalDof[
    shapeDerivative[faceTerm], physicalDof], virtualDof];
  If[TrueQ[reverseDirectSource],
    applyCorruptedProfileJets[processed, ablatedDirection],
    applyProfileJets[processed, ablatedDirection]]
], {sign, faceSigns}]];

repEulerianOperand = <||>;
repMaterialOperand = <||>;

Do[Module[{key, eulerian, material},
  key = "VIRTUAL_CONSTRAINT|" <> densityCaseKey[branch, density, dof];
  eulerian = applyProfileJets[applyVirtualDof[
    virtualConstraintDirectSource[branch, density], dof]];
  material = applyProfileJets[applyVirtualDof[
    virtualConstraintMaterialSource[branch, density], dof]];
  AssociateTo[repEulerianOperand, key -> objectRecord[eulerian,
    dimensionZero]];
  AssociateTo[repMaterialOperand, key -> objectRecord[material,
    dimensionZero]];
], {branch, branchNames}, {density, densityNames}, {dof, dofNames}];

Do[Module[{key, eulerian, material},
  key = "EVOLUTION_MASS_BALANCE|" <>
    densityCaseKey[branch, density, dof];
  eulerian = applyProfileJets[applyPhysicalDof[shapeDerivative[
    Total[Values[evolutionTermsDirectSource[branch, density]]]], dof]];
  material = applyProfileJets[applyPhysicalDof[shapeDerivative[
    Total[Values[evolutionTermsMaterialSource[branch, density]]]], dof]];
  AssociateTo[repEulerianOperand, key -> objectRecord[eulerian,
    dimensionEvolution]];
  AssociateTo[repMaterialOperand, key -> objectRecord[material,
    dimensionEvolution]];
], {branch, branchNames}, {density, densityNames}, {dof, dofNames}];

Do[Module[{key, eulerian, material, dimension},
  dimension = Switch[object,
    "RELATIVE_FLUX", dimensionFlux,
    "TRACTION", dimensionPressure,
    "CLOSURE_SHAPE_DERIV", dimensionFlux];
  key = object <> "|" <> densityFaceCaseKey[branch, density, sign, dof];
  eulerian = applyProfileJets[applyPhysicalDof[
    directShape[branch, sign, density, Switch[object,
      "RELATIVE_FLUX", "FLUX", "TRACTION", "TRACTION",
      "CLOSURE_SHAPE_DERIV", "CLOSURE"]], dof]];
  material = applyProfileJets[applyPhysicalDof[
    materialShape[branch, sign, density, Switch[object,
      "RELATIVE_FLUX", "FLUX", "TRACTION", "TRACTION",
      "CLOSURE_SHAPE_DERIV", "CLOSURE"]], dof]];
  AssociateTo[repEulerianOperand, key -> objectRecord[eulerian, dimension]];
  AssociateTo[repMaterialOperand, key -> objectRecord[material, dimension]];
], {object, {"RELATIVE_FLUX", "TRACTION", "CLOSURE_SHAPE_DERIV"}},
  {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {dof, dofNames}];

Do[Module[{key, eulerian, material},
  key = "VIRTUAL_WORK_SHAPE_DERIV|" <>
    densityCaseKey[branch, density, physicalDof] <>
    "|VIRTUAL_DOF_" <> virtualDof;
  eulerian = processedVirtualWorkShape["EULERIAN", branch, density,
    physicalDof, virtualDof];
  material = processedVirtualWorkShape["MATERIAL", branch, density,
    physicalDof, virtualDof];
  AssociateTo[repEulerianOperand, key -> objectRecord[eulerian,
    dimensionVirtualWork]];
  AssociateTo[repMaterialOperand, key -> objectRecord[material,
    dimensionVirtualWork]];
], {branch, branchNames}, {density, densityNames},
  {physicalDof, dofNames}, {virtualDof, dofNames}];

emitShared["REP_INVARIANCE_EULERIAN_OPERAND", repEulerianOperand];
emitShared["REP_INVARIANCE_MATERIAL_OPERAND", repMaterialOperand];

repResidual = AssociationMap[Function[key, Module[{left, right},
  left = repEulerianOperand[key]["EXPRESSION"];
  right = repMaterialOperand[key]["EXPRESSION"];
  comparisonPayload[Together[left - right],
    repEulerianOperand[key]["DIMENSION_L_T_M"]]
]], Keys[repEulerianOperand]];
emitShared["REP_INVARIANCE_RESIDUAL", repResidual];

independenceBaseOperand = <||>;
independenceCorruptedOperand = <||>;

Do[Module[{key, base, corrupted},
  key = "VIRTUAL_CONSTRAINT|" <> densityCaseKey[branch, density, dof];
  base = applyProfileJets[applyVirtualDof[
    virtualConstraintDirectSource[branch, density], dof]];
  corrupted = applyProfileJets[applyVirtualDof[
    virtualConstraintDirectSource[branch, density, True], dof]];
  AssociateTo[independenceBaseOperand, key -> objectRecord[base,
    dimensionZero]];
  AssociateTo[independenceCorruptedOperand, key -> objectRecord[corrupted,
    dimensionZero]];
], {branch, branchNames}, {density, densityNames}, {dof, dofNames}];

Do[Module[{key, base, corrupted},
  key = "EVOLUTION_MASS_BALANCE|" <>
    densityCaseKey[branch, density, dof];
  base = applyProfileJets[applyPhysicalDof[shapeDerivative[
    Total[Values[evolutionTermsDirectSource[branch, density]]]], dof]];
  corrupted = applyProfileJets[applyPhysicalDof[shapeDerivative[
    Total[Values[evolutionTermsDirectSource[branch, density, True]]]],
    dof]];
  AssociateTo[independenceBaseOperand, key -> objectRecord[base,
    dimensionEvolution]];
  AssociateTo[independenceCorruptedOperand, key -> objectRecord[corrupted,
    dimensionEvolution]];
], {branch, branchNames}, {density, densityNames}, {dof, dofNames}];

Do[Module[{key, sourceName, dimension, base, corrupted},
  sourceName = Switch[object, "RELATIVE_FLUX", "FLUX",
    "TRACTION", "TRACTION", "CLOSURE_SHAPE_DERIV", "CLOSURE"];
  dimension = Switch[object, "TRACTION", dimensionPressure, _,
    dimensionFlux];
  key = object <> "|" <> densityFaceCaseKey[branch, density, sign, dof];
  base = applyProfileJets[applyPhysicalDof[
    directShape[branch, sign, density, sourceName], dof]];
  corrupted = If[sign == 1,
    applyCorruptedProfileJets[applyPhysicalDof[
      corruptedUpperDirectShape[branch, density, sourceName], dof]],
    applyProfileJets[applyPhysicalDof[
      directShape[branch, sign, density, sourceName], dof]]];
  AssociateTo[independenceBaseOperand, key -> objectRecord[base, dimension]];
  AssociateTo[independenceCorruptedOperand, key -> objectRecord[corrupted,
    dimension]];
], {object, {"RELATIVE_FLUX", "TRACTION", "CLOSURE_SHAPE_DERIV"}},
  {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {dof, dofNames}];

Do[Module[{key, base, corrupted},
  key = "VIRTUAL_WORK_SHAPE_DERIV|" <>
    densityCaseKey[branch, density, physicalDof] <>
    "|VIRTUAL_DOF_" <> virtualDof;
  base = processedVirtualWorkShape["EULERIAN", branch, density,
    physicalDof, virtualDof];
  corrupted = processedVirtualWorkShape["EULERIAN", branch, density,
    physicalDof, virtualDof, 0, True];
  AssociateTo[independenceBaseOperand, key -> objectRecord[base,
    dimensionVirtualWork]];
  AssociateTo[independenceCorruptedOperand, key -> objectRecord[corrupted,
    dimensionVirtualWork]];
], {branch, branchNames}, {density, densityNames},
  {physicalDof, dofNames}, {virtualDof, dofNames}];

emitShared["CONTROL_INDEPENDENCE_BASE_OPERAND", independenceBaseOperand];
emitShared["CONTROL_INDEPENDENCE_CORRUPTED_OPERAND",
  independenceCorruptedOperand];

independenceResidual = AssociationMap[Function[key, Module[{base, corrupt},
  base = independenceBaseOperand[key]["EXPRESSION"];
  corrupt = independenceCorruptedOperand[key]["EXPRESSION"];
  comparisonPayload[Together[base - corrupt],
    independenceBaseOperand[key]["DIMENSION_L_T_M"]]
]], Keys[independenceBaseOperand]];
emitShared["CONTROL_INDEPENDENCE_RESIDUAL", independenceResidual];

(* ---------------------------------------------------------------------- *)
(* Per-direction C-1 source-form ablations.                               *)
(* ---------------------------------------------------------------------- *)

formBaseOperand = <||>;
formAblatedOperand = <||>;

Do[Module[{source, dimension, key, base, ablated},
  source = Switch[object,
    "FACE_NORMAL", shapeDerivative[graphNormalSource[branch, sign]],
    "CONORMAL_DERIV", shapeDerivative[conormalSource[branch, sign]],
    "FACE_MEASURE_SHAPE_DERIV",
      shapeDerivative[graphMeasureSource[branch, sign]],
    "FACE_VELOCITY", directShape[branch, sign,
      "RHO4_CONSTANT", "NORMAL_VELOCITY"]];
  dimension = Switch[object,
    "FACE_NORMAL", dimensionZero,
    "CONORMAL_DERIV", -dimensionL,
    "FACE_MEASURE_SHAPE_DERIV", dimensionZero,
    "FACE_VELOCITY", dimensionVelocity];
  key = object <> "|" <> faceCaseKey[branch, sign, dof] <>
    "|DIRECTION_" <> ToString[direction];
  base = applyProfileJets[applyPhysicalDof[source, dof]];
  ablated = applyProfileJets[applyPhysicalDof[source, dof], direction];
  AssociateTo[formBaseOperand, key -> objectRecord[base, dimension]];
  AssociateTo[formAblatedOperand, key -> objectRecord[ablated, dimension]];
], {object, {"FACE_NORMAL", "CONORMAL_DERIV",
    "FACE_MEASURE_SHAPE_DERIV", "FACE_VELOCITY"}},
  {branch, branchNames}, {sign, faceSigns}, {dof, dofNames},
  {direction, spatialDirections}];

Do[Module[{sourceName, source, dimension, key, base, ablated},
  sourceName = Switch[object,
    "RELATIVE_FLUX", "FLUX", "KINEMATIC_BALANCE", "KINEMATIC",
    "TRACTION", "TRACTION", "CLOSURE_SHAPE_DERIV", "CLOSURE"];
  source = directShape[branch, sign, density, sourceName];
  dimension = Switch[object,
    "KINEMATIC_BALANCE", dimensionVelocity,
    "TRACTION", dimensionPressure, _, dimensionFlux];
  key = object <> "|" <> densityFaceCaseKey[branch, density, sign, dof] <>
    "|DIRECTION_" <> ToString[direction];
  base = applyProfileJets[applyPhysicalDof[source, dof]];
  ablated = applyProfileJets[applyPhysicalDof[source, dof], direction];
  AssociateTo[formBaseOperand, key -> objectRecord[base, dimension]];
  AssociateTo[formAblatedOperand, key -> objectRecord[ablated, dimension]];
], {object, {"RELATIVE_FLUX", "KINEMATIC_BALANCE", "TRACTION",
    "CLOSURE_SHAPE_DERIV"}}, {branch, branchNames},
  {density, densityNames}, {sign, faceSigns}, {dof, dofNames},
  {direction, spatialDirections}];

Do[Module[{inventory, source, dimension, key, base, ablated},
  inventory = traceFieldInventory[fieldName];
  source = shapeDerivative[traceSource[inventory[[1]], inventory[[2]],
    branch, sign]];
  dimension = inventory[[3]];
  key = "FACE_SHIFT|" <> faceCaseKey[branch, sign, dof] <>
    "|FIELD_" <> fieldName <> "|DIRECTION_" <> ToString[direction];
  base = applyProfileJets[applyPhysicalDof[source, dof]];
  ablated = applyProfileJets[applyPhysicalDof[source, dof], direction];
  AssociateTo[formBaseOperand, key -> objectRecord[base, dimension]];
  AssociateTo[formAblatedOperand, key -> objectRecord[ablated, dimension]];
], {branch, branchNames}, {sign, faceSigns},
  {fieldName, Keys[traceFieldInventory]}, {dof, dofNames},
  {direction, spatialDirections}];

Do[Module[{objects, faceTerm, source, key, base, ablated},
  objects = directFaceObjectsSource[branch, sign, density];
  faceTerm = objects["MEASURE"] objects["TRACTION"].
    virtualFaceDisplacementSource[branch, sign];
  source = applyVirtualDof[applyPhysicalDof[
    shapeDerivative[faceTerm], physicalDof], virtualDof];
  key = "VIRTUAL_WORK_SHAPE_DERIV|" <>
    densityFaceCaseKey[branch, density, sign, physicalDof] <>
    "|VIRTUAL_DOF_" <> virtualDof <> "|DIRECTION_" <>
    ToString[direction];
  base = applyProfileJets[source];
  ablated = applyProfileJets[source, direction];
  AssociateTo[formBaseOperand, key -> objectRecord[base,
    dimensionVirtualWork]];
  AssociateTo[formAblatedOperand, key -> objectRecord[ablated,
    dimensionVirtualWork]];
], {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {physicalDof, dofNames}, {virtualDof, dofNames},
  {direction, spatialDirections}];

Do[Module[{source, key, base, ablated},
  source = shapeDerivative[projectionOperand[projectionDynamicTerms[branch]]];
  key = "PROJECTION_SHAPE_DERIV|" <> faceCaseKey[branch, sign, dof] <>
    "|DIRECTION_" <> ToString[direction];
  base = applyProfileJets[applyPhysicalDof[source, dof]];
  ablated = applyProfileJets[applyPhysicalDof[source, dof], direction];
  AssociateTo[formBaseOperand, key -> objectRecord[base,
    dimensionEvolution]];
  AssociateTo[formAblatedOperand, key -> objectRecord[ablated,
    dimensionEvolution]];
], {branch, branchNames}, {sign, faceSigns}, {dof, dofNames},
  {direction, spatialDirections}];

Do[Module[{virtualSource, evolutionSource, virtualKey, evolutionKey,
    virtualBase, virtualAblated, evolutionBase, evolutionAblated},
  virtualSource = applyVirtualDof[
    virtualConstraintDirectSource[branch, density], dof];
  evolutionSource = applyPhysicalDof[shapeDerivative[
    Total[Values[evolutionTermsDirectSource[branch, density]]]], dof];
  virtualKey = "VIRTUAL_CONSTRAINT|" <>
    densityFaceCaseKey[branch, density, sign, dof] <>
    "|DIRECTION_" <> ToString[direction];
  evolutionKey = "EVOLUTION_MASS_BALANCE|" <>
    densityFaceCaseKey[branch, density, sign, dof] <>
    "|DIRECTION_" <> ToString[direction];
  virtualBase = applyProfileJets[virtualSource];
  virtualAblated = applyProfileJets[virtualSource, direction];
  evolutionBase = applyProfileJets[evolutionSource];
  evolutionAblated = applyProfileJets[evolutionSource, direction];
  AssociateTo[formBaseOperand, {
    virtualKey -> objectRecord[virtualBase, dimensionZero],
    evolutionKey -> objectRecord[evolutionBase, dimensionEvolution]}];
  AssociateTo[formAblatedOperand, {
    virtualKey -> objectRecord[virtualAblated, dimensionZero],
    evolutionKey -> objectRecord[evolutionAblated, dimensionEvolution]}];
], {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {dof, dofNames}, {direction, spatialDirections}];

emitShared["CONTROL_FORM_BASE_OPERAND", formBaseOperand];
emitShared["CONTROL_FORM_ABLATED_OPERAND", formAblatedOperand];

formResidual = AssociationMap[Function[key, Module[{base, ablated},
  base = formBaseOperand[key]["EXPRESSION"];
  ablated = formAblatedOperand[key]["EXPRESSION"];
  comparisonPayload[Together[base - ablated],
    formBaseOperand[key]["DIMENSION_L_T_M"]]
]], Keys[formBaseOperand]];
emitShared["CONTROL_FORM_RESIDUAL", formResidual];

(* ---------------------------------------------------------------------- *)
(* Uniform-limit regression: independent profile cut and flat sources.    *)
(* ---------------------------------------------------------------------- *)

uniformProfileRules = {
  HoldPattern[Derivative[orders__][widthProfile][arguments__]] :> 0,
  HoldPattern[Derivative[orders__][muRProfile][arguments__]] :> 0,
  HoldPattern[widthProfile[arguments__]] :> W0,
  HoldPattern[muRProfile[arguments__]] :> muR
};
uniformS11ca[expression_] := computedExpression[
  applyProfileJets[expression] /. {etaBg -> 0, sigmaW -> 0}];
uniformS11b[expression_] := computedExpression[
  expression /. materialToEulerianRules /. uniformProfileRules];

uniformS11caOperand = <||>;
uniformS11bOperand = <||>;

Do[Module[{key, sigmaSource, gradientSource},
  sigmaSource = Switch[density,
    "RHO4_CONSTANT", (rhoBr/W0) widthAnsatz,
    "RHOBR_CONSTANT", rhoBr];
  gradientSource = D[sigmaSource, #] & /@ {y1, y2, y3};
  key = "BACKGROUND_DENSITY_MAP|" <> density <> "|SIGMA_E_ZERO";
  AssociateTo[uniformS11caOperand, key -> objectRecord[
    computedExpression[sigmaSource /. etaBg -> 0],
    dimensionBraneDensity, 0]];
  AssociateTo[uniformS11bOperand, key -> objectRecord[rhoBr,
    dimensionBraneDensity, 0]];
  key = "BACKGROUND_DENSITY_MAP|" <> density <>
    "|GRADIENT_SIGMA_E_ZERO";
  AssociateTo[uniformS11caOperand, key -> objectRecord[
    computedExpression[gradientSource /. etaBg -> 0],
    dimensionBraneDensity - dimensionL, 0]];
  AssociateTo[uniformS11bOperand, key -> objectRecord[{0, 0, 0},
    dimensionBraneDensity - dimensionL, 0]];
], {density, densityNames}];

Do[Module[{source, dimension, key, s11ca, s11b},
  source = Switch[object,
    "FACE_NORMAL", shapeDerivative[graphNormalSource[branch, sign]],
    "CONORMAL_DERIV", shapeDerivative[conormalSource[branch, sign]],
    "FACE_MEASURE_SHAPE_DERIV",
      shapeDerivative[graphMeasureSource[branch, sign]],
    "FACE_VELOCITY", directShape[branch, sign,
      "RHO4_CONSTANT", "NORMAL_VELOCITY"]];
  dimension = Switch[object,
    "FACE_NORMAL", dimensionZero, "CONORMAL_DERIV", -dimensionL,
    "FACE_MEASURE_SHAPE_DERIV", dimensionZero,
    "FACE_VELOCITY", dimensionVelocity];
  key = object <> "|" <> faceCaseKey[branch, sign, dof];
  source = applyPhysicalDof[source, dof];
  s11ca = uniformS11ca[source];
  s11b = uniformS11b[source];
  AssociateTo[uniformS11caOperand, key -> objectRecord[s11ca, dimension]];
  AssociateTo[uniformS11bOperand, key -> objectRecord[s11b, dimension]];
], {object, {"FACE_NORMAL", "CONORMAL_DERIV",
    "FACE_MEASURE_SHAPE_DERIV", "FACE_VELOCITY"}},
  {branch, branchNames}, {sign, faceSigns}, {dof, dofNames}];

Do[Module[{sourceName, source, dimension, key},
  sourceName = Switch[object,
    "RELATIVE_FLUX", "FLUX", "KINEMATIC_BALANCE", "KINEMATIC",
    "TRACTION", "TRACTION", "CLOSURE_SHAPE_DERIV", "CLOSURE"];
  source = applyPhysicalDof[directShape[branch, sign, density, sourceName],
    dof];
  dimension = Switch[object, "KINEMATIC_BALANCE", dimensionVelocity,
    "TRACTION", dimensionPressure, _, dimensionFlux];
  key = object <> "|" <> densityFaceCaseKey[branch, density, sign, dof];
  AssociateTo[uniformS11caOperand, key -> objectRecord[
    uniformS11ca[source], dimension]];
  AssociateTo[uniformS11bOperand, key -> objectRecord[
    uniformS11b[source], dimension]];
], {object, {"RELATIVE_FLUX", "KINEMATIC_BALANCE", "TRACTION",
    "CLOSURE_SHAPE_DERIV"}}, {branch, branchNames},
  {density, densityNames}, {sign, faceSigns}, {dof, dofNames}];

Do[Module[{inventory, source, key, dimension},
  inventory = traceFieldInventory[fieldName];
  source = applyPhysicalDof[shapeDerivative[traceSource[inventory[[1]],
    inventory[[2]], branch, sign]], dof];
  dimension = inventory[[3]];
  key = "FACE_SHIFT|" <> faceCaseKey[branch, sign, dof] <>
    "|FIELD_" <> fieldName;
  AssociateTo[uniformS11caOperand, key -> objectRecord[
    uniformS11ca[source], dimension]];
  AssociateTo[uniformS11bOperand, key -> objectRecord[
    uniformS11b[source], dimension]];
], {branch, branchNames}, {sign, faceSigns},
  {fieldName, Keys[traceFieldInventory]}, {dof, dofNames}];

Do[Module[{source, key},
  source = applyPhysicalDof[shapeDerivative[
    projectionOperand[projectionDynamicTerms[branch]]], dof];
  key = "PROJECTION_SHAPE_DERIV|" <> branch <> "|DOF_" <> dof;
  AssociateTo[uniformS11caOperand, key -> objectRecord[
    uniformS11ca[source], dimensionEvolution]];
  AssociateTo[uniformS11bOperand, key -> objectRecord[
    uniformS11b[source], dimensionEvolution]];
], {branch, branchNames}, {dof, dofNames}];

Do[Module[{virtualSource, evolutionSource, virtualKey, evolutionKey},
  virtualSource = applyVirtualDof[
    virtualConstraintDirectSource[branch, density], dof];
  evolutionSource = applyPhysicalDof[shapeDerivative[
    Total[Values[evolutionTermsDirectSource[branch, density]]]], dof];
  virtualKey = "VIRTUAL_CONSTRAINT|" <>
    densityCaseKey[branch, density, dof];
  evolutionKey = "EVOLUTION_MASS_BALANCE|" <>
    densityCaseKey[branch, density, dof];
  AssociateTo[uniformS11caOperand, {
    virtualKey -> objectRecord[uniformS11ca[virtualSource], dimensionZero],
    evolutionKey -> objectRecord[uniformS11ca[evolutionSource],
      dimensionEvolution]}];
  AssociateTo[uniformS11bOperand, {
    virtualKey -> objectRecord[uniformS11b[virtualSource], dimensionZero],
    evolutionKey -> objectRecord[uniformS11b[evolutionSource],
      dimensionEvolution]}];
], {branch, branchNames}, {density, densityNames}, {dof, dofNames}];

Do[Module[{objects, faceTerm, source, key},
  objects = directFaceObjectsSource[branch, sign, density];
  faceTerm = objects["MEASURE"] objects["TRACTION"].
    virtualFaceDisplacementSource[branch, sign];
  source = applyVirtualDof[applyPhysicalDof[shapeDerivative[faceTerm],
    physicalDof], virtualDof];
  key = "VIRTUAL_WORK_SHAPE_DERIV|" <>
    densityFaceCaseKey[branch, density, sign, physicalDof] <>
    "|VIRTUAL_DOF_" <> virtualDof;
  AssociateTo[uniformS11caOperand, key -> objectRecord[
    uniformS11ca[source], dimensionVirtualWork]];
  AssociateTo[uniformS11bOperand, key -> objectRecord[
    uniformS11b[source], dimensionVirtualWork]];
], {branch, branchNames}, {density, densityNames}, {sign, faceSigns},
  {physicalDof, dofNames}, {virtualDof, dofNames}];

emitShared["UNIFORM_LIMIT_S11CA_OPERAND", uniformS11caOperand];
emitShared["UNIFORM_LIMIT_S11B_OPERAND", uniformS11bOperand];

uniformResidual = AssociationMap[Function[key, Module[{left, right},
  left = uniformS11caOperand[key]["EXPRESSION"];
  right = uniformS11bOperand[key]["EXPRESSION"];
  comparisonPayload[Together[left - right],
    uniformS11caOperand[key]["DIMENSION_L_T_M"]]
]], Keys[uniformS11caOperand]];
emitShared["UNIFORM_LIMIT_RESIDUAL", uniformResidual];

(* ---------------------------------------------------------------------- *)
(* Restored dimensions and source-level homogeneity corruption.           *)
(* ---------------------------------------------------------------------- *)

dimensionLambdaA = dimensionFlux - dimensionAffinity;
dimensionLambdaV = dimensionFlux - dimensionVelocity;
dimensionLambdaX = dimensionPressure - dimensionAffinity;

dimensionInventory = <|
  "BACKGROUND_DENSITY_MAP" -> <|"RHO4_BG" -> dimensionBulkDensity,
    "RHOBR_BG" -> dimensionBraneDensity,
    "SIGMA_E_ZERO" -> dimensionBraneDensity,
    "GRADIENT_SIGMA_E_ZERO" -> dimensionBraneDensity - dimensionL|>,
  "BACKGROUND_STATE" -> <|"W_BG" -> dimensionL,
    "MU_R_BG" -> dimensionEnergyBrane,
    "RHO4_BG" -> dimensionBulkDensity,
    "RHOBR_BG" -> dimensionBraneDensity,
    "FACE_VELOCITY" -> dimensionVelocity, "FACE_FLUX" -> dimensionFlux,
    "AFFINITY" -> dimensionAffinity, "BOUNDARY_LOAD" -> dimensionPressure|>,
  "ADMISSIBILITY_PREMISE" -> <|
    "BODY_SUPPORT" -> dimensionEnergyBrane - dimensionL,
    "FACE_SUPPORT" -> dimensionPressure|>,
  "FACE_MAP_LAB_HELD" -> dimensionL,
  "FACE_MAP_MATERIAL_ADVECTED" -> dimensionL,
  "FACE_NORMAL" -> dimensionZero,
  "CONORMAL_DERIV" -> -dimensionL,
  "FACE_MEASURE_SHAPE_DERIV" -> dimensionZero,
  "FACE_VELOCITY" -> dimensionVelocity,
  "RELATIVE_FLUX" -> dimensionFlux,
  "KINEMATIC_BALANCE" -> dimensionVelocity,
  "TRACTION" -> dimensionPressure,
  "VIRTUAL_WORK_SHAPE_DERIV" -> dimensionVirtualWork,
  "FACE_SHIFT" -> <|"PRESSURE" -> dimensionPressure,
    "BULK_DENSITY" -> dimensionBulkDensity,
    "CURRENT" -> dimensionFlux, "BULK_VELOCITY" -> dimensionVelocity|>,
  "PROJECTION_SHAPE_DERIV" -> dimensionEvolution,
  "PROJECTION_STATIC_OPERAND" -> dimensionEvolution,
  "PROJECTION_DYNAMIC_OPERAND" -> dimensionEvolution,
  "PROJECTION_RESIDUAL" -> dimensionEvolution,
  "PROJECTION_TERM_ORIGINS" -> dimensionEvolution,
  "VIRTUAL_CONSTRAINT" -> dimensionZero,
  "EVOLUTION_MASS_BALANCE" -> dimensionEvolution,
  "EVOLUTION_TERM_ORIGINS" -> dimensionEvolution,
  "CLOSURE_SHAPE_DERIV" -> dimensionFlux,
  "REP_INVARIANCE_EULERIAN_OPERAND" ->
    {dimensionZero, dimensionEvolution, dimensionFlux,
      dimensionPressure, dimensionVirtualWork},
  "REP_INVARIANCE_MATERIAL_OPERAND" ->
    {dimensionZero, dimensionEvolution, dimensionFlux,
      dimensionPressure, dimensionVirtualWork},
  "REP_INVARIANCE_RESIDUAL" ->
    {dimensionZero, dimensionEvolution, dimensionFlux,
      dimensionPressure, dimensionVirtualWork},
  "CONTROL_INDEPENDENCE_BASE_OPERAND" ->
    {dimensionZero, dimensionEvolution, dimensionFlux,
      dimensionPressure, dimensionVirtualWork},
  "CONTROL_INDEPENDENCE_CORRUPTED_OPERAND" ->
    {dimensionZero, dimensionEvolution, dimensionFlux,
      dimensionPressure, dimensionVirtualWork},
  "CONTROL_INDEPENDENCE_RESIDUAL" ->
    {dimensionZero, dimensionEvolution, dimensionFlux,
      dimensionPressure, dimensionVirtualWork},
  "CONTROL_FORM_BASE_OPERAND" -> Values[dimensionInventoryForm = <|
    "NORMAL" -> dimensionZero, "CONORMAL" -> -dimensionL,
    "MEASURE" -> dimensionZero, "VELOCITY" -> dimensionVelocity,
    "FLUX" -> dimensionFlux, "TRACTION" -> dimensionPressure,
    "WORK" -> dimensionVirtualWork, "PROJECTION" -> dimensionEvolution,
    "CONSTRAINT" -> dimensionZero, "EVOLUTION" -> dimensionEvolution|>],
  "CONTROL_FORM_ABLATED_OPERAND" -> Values[dimensionInventoryForm],
  "CONTROL_FORM_RESIDUAL" -> Values[dimensionInventoryForm],
  "UNIFORM_LIMIT_S11CA_OPERAND" -> Values[dimensionInventoryForm],
  "UNIFORM_LIMIT_S11B_OPERAND" -> Values[dimensionInventoryForm],
  "UNIFORM_LIMIT_RESIDUAL" -> Values[dimensionInventoryForm],
  "HOMOGENEITY_BASE_OPERAND" -> dimensionZero,
  "HOMOGENEITY_CONTROL_OPERAND" -> dimensionZero,
  "HOMOGENEITY_RESIDUAL" -> dimensionZero,
  "DIMENSIONS" -> dimensionZero|>;

dimensionsPayload = <|"EXPRESSION" -> dimensionInventory,
  "BASE_UNITS" -> {lengthUnit, timeUnit, massUnit},
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>;
emitShared["DIMENSIONS", dimensionsPayload];

homogeneityBaseDimensions = <|
  "RELATIVE_FLUX" -> {dimensionFlux,
    dimensionBulkDensity + dimensionVelocity},
  "KINEMATIC_BALANCE" -> {dimensionVelocity, dimensionVelocity,
    dimensionFlux - dimensionBulkDensity},
  "TRACTION" -> {dimensionPressure,
    dimensionLambdaX + dimensionAffinity},
  "VIRTUAL_WORK" -> {dimensionZero + dimensionPressure + dimensionL,
    dimensionVirtualWork},
  "PROJECTION" -> ConstantArray[dimensionEvolution, 5],
  "VIRTUAL_CONSTRAINT" -> ConstantArray[dimensionZero, 3],
  "EVOLUTION" -> ConstantArray[dimensionEvolution, 3],
  "CLOSURE" -> {dimensionFlux,
    dimensionLambdaA + dimensionAffinity,
    dimensionLambdaV + dimensionVelocity}|>;

homogeneityControlDimensions = AssociationMap[
  Function[name, ReplacePart[homogeneityBaseDimensions[name],
    1 -> (First[homogeneityBaseDimensions[name]] + dimensionL)]],
  Keys[homogeneityBaseDimensions]];

homogeneityBasePayload = AssociationMap[Function[name, <|
  "SOURCE_TERM_DIMENSIONS" -> homogeneityBaseDimensions[name],
  "HOMOGENEITY_OBJECT" -> dimensionalSum[
    homogeneityBaseDimensions[name]],
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>],
  Keys[homogeneityBaseDimensions]];
emitShared["HOMOGENEITY_BASE_OPERAND", homogeneityBasePayload];

homogeneityControlPayload = AssociationMap[Function[name, <|
  "SOURCE_LEVEL_CORRUPTION" -> HoldForm[
    firstSourceTerm[name] -> lengthCorruptor firstSourceTerm[name]],
  "SOURCE_TERM_DIMENSIONS" -> homogeneityControlDimensions[name],
  "HOMOGENEITY_OBJECT" -> dimensionalSum[
    homogeneityControlDimensions[name]],
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>],
  Keys[homogeneityControlDimensions]];
emitShared["HOMOGENEITY_CONTROL_OPERAND", homogeneityControlPayload];

homogeneityResidualPayload = AssociationMap[Function[name, Module[
  {base, control, residual},
  base = homogeneityBaseDimensions[name];
  control = homogeneityControlDimensions[name];
  residual = MapThread[#1 - #2 &, {base, control}];
  <|"BASE_OPERAND" -> base, "CONTROL_OPERAND" -> control,
    "RESIDUAL" -> residual,
    "TEST_OBJECTS" -> MapThread[relationalObject, {base, control}],
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
    "DIMENSION_L_T_M" -> dimensionZero|>
]], Keys[homogeneityBaseDimensions]];
emitShared["HOMOGENEITY_RESIDUAL", homogeneityResidualPayload];

localInventoryName = standardEmissionName[localObject["TAG_NAMES"]];
emitLocal["TAG_NAMES", <|"EXPRESSION" -> Append[localNames,
    localInventoryName],
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>];
