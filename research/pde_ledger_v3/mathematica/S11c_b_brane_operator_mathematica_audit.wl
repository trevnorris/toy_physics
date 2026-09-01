$HistoryLength = 0;
ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};

(* The public strings are kept at the output boundary.  Every symbolic object
   below is constructed from the three supplied shared-physics specifications. *)
sharedObject[name_String] := HoldComplete[SharedS11cbObject, name];
localObject[name_String] := HoldComplete[LocalS11cbObject, name];

standardEmissionName[HoldComplete[SharedS11cbObject, name_String]] :=
  "WL_S11CB_" <> name;
standardEmissionName[HoldComplete[LocalS11cbObject, name_String]] :=
  "WL_S11CB_LOCAL_" <> name;

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
  If[StringStartsQ[name, "WL_S11CB_LOCAL_"], AppendTo[localNames, name]];
  cleanPayload = stripConditional[payload];
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
  If[StringStartsQ[name, "WL_S11CB_LOCAL_"], AppendTo[localNames, name]];
  stream = First[$Output];
  WriteString[stream, name <> ": <|"];
  Flush[stream];
];

appendAssociationEmission[key_String, value_, first_:False] := Module[
  {stream, rendered},
  stream = First[$Output];
  rendered = ToString[stripConditional[value], InputForm,
    PageWidth -> Infinity];
  WriteString[stream, If[TrueQ[first], "", ", "] <>
    ToString[key, InputForm, PageWidth -> Infinity] <> " -> " <>
    rendered];
  Flush[stream];
];

endAssociationEmission[] := Module[{stream},
  stream = First[$Output];
  WriteString[stream, "|>\n"];
  Flush[stream];
];

relationalObject[left_, right_] := Inactive[Equal][left, right];

(* ---------------------------------------------------------------------- *)
(* Dimensions and multigrades.                                           *)
(* ---------------------------------------------------------------------- *)

dimensionL = {1, 0, 0};
dimensionT = {0, 1, 0};
dimensionM = {0, 0, 1};
dimensionZero = {0, 0, 0};

dimensionX = dimensionL;
dimensionTime = dimensionT;
dimensionU = dimensionL;
dimensionTheta = dimensionZero;
dimensionEw = dimensionZero;
dimensionGradient = -dimensionL;
dimensionRhoBr = dimensionM - 3 dimensionL;
dimensionRhoFour = dimensionRhoBr - dimensionL;
dimensionVelocity = dimensionL - dimensionT;
dimensionEnergy = dimensionRhoBr + 2 dimensionVelocity;
dimensionBulkForce = dimensionEnergy - dimensionU;
dimensionFaceTraction = dimensionBulkForce;
dimensionMassEvolution = dimensionRhoBr - dimensionT;
dimensionFlux = dimensionRhoFour + dimensionVelocity;
dimensionPressure = dimensionRhoFour + 2 dimensionVelocity;
dimensionAffinity = dimensionPressure - dimensionRhoFour;
dimensionLambdaA = dimensionFlux - dimensionAffinity;
dimensionLambdaV = dimensionFlux - dimensionVelocity;
dimensionLambdaX = dimensionPressure - dimensionAffinity;
dimensionMuTheta = dimensionEnergy;
dimensionThicknessRow = dimensionEnergy;

truncateScalar[expression_] := Module[
  {operators, tokens, protected, first, second, restored},
  operators = DeleteDuplicates[Cases[expression,
    operator : Inactive[_][___] :> operator, {0, Infinity}]];
  If[operators === {},
    first = Quiet[Check[Normal[Series[expression, {etaBg, 0, 1}]],
        expression]];
    Return[Quiet[Check[Normal[Series[first, {sigmaW, 0, 1}]],
      first]]]
  ];
  tokens = Table[inactiveOperatorToken[index],
    {index, Length[operators]}];
  protected = expression /. Thread[operators -> tokens];
  first = Quiet[Check[Normal[Series[protected, {etaBg, 0, 1}]],
      protected]];
  second = Quiet[Check[Normal[Series[first, {sigmaW, 0, 1}]], first]];
  restored = second /. Thread[tokens -> (truncateBackground /@ operators)];
  restored
];
truncateBackground[expression_Association] :=
  Map[truncateBackground, expression];
truncateBackground[expression_List] := truncateBackground /@ expression;
truncateBackground[expression_Rule] :=
  Rule[expression[[1]], truncateBackground[expression[[2]]]];
truncateBackground[Inactive[head_][arguments___]] :=
  Inactive[head] @@ (truncateBackground /@ {arguments});
truncateBackground[expression_] := truncateScalar[expression];

gradeTerms[expression_Association] := DeleteDuplicates[Sort[
  Flatten[gradeTerms /@ Values[expression], 1]]];
gradeTerms[expression_List] := DeleteDuplicates[Sort[
  Flatten[gradeTerms /@ expression, 1]]];
gradeTerms[expression_] := Module[{coefficient, grades, operators,
    protectedGrades, operatorGrades},
  If[!FreeQ[expression, _String], Return[{{0, 0}}]];
  operators = DeleteDuplicates[Cases[expression,
    operator : Inactive[_][___] :> operator, {0, Infinity}]];
  If[operators =!= {},
    protectedGrades = gradeTerms[expression /. Thread[operators -> 0]];
    operatorGrades = Flatten[gradeTerms /@ (List @@@ operators), 1];
    Return[DeleteDuplicates[Sort[Join[protectedGrades,
      operatorGrades]]]]
  ];
  grades = Reap[Do[
    coefficient = Coefficient[Coefficient[expression, etaBg, etaOrder],
      sigmaW, sigmaOrder];
    If[!TrueQ[PossibleZeroQ[coefficient]],
      Sow[{etaOrder, sigmaOrder}]],
    {etaOrder, 0, 1}, {sigmaOrder, 0, 1}]][[2]];
  If[grades === {}, {}, DeleteDuplicates[Sort[First[grades]]]]
];

withEpsilonGrade[expression_, epsilonOrder_Integer] :=
  ({epsilonOrder, #[[1]], #[[2]]} & /@ gradeTerms[expression]);

objectRecord[expression_, dimension_, epsilonOrder_:1] := Module[
  {computed, grades},
  computed = truncateBackground[expression];
  grades = gradeTerms[computed] /. {a_, b_} :> {epsilonOrder, a, b};
  <|"EXPRESSION" -> computed,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> grades,
    "DIMENSION_L_T_M" -> dimension|>
];

preparedObjectRecord[expression_, dimension_, epsilonOrder_:1] := <|
  "EXPRESSION" -> expression,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
    (gradeTerms[expression] /. {a_, b_} :> {epsilonOrder, a, b}),
  "DIMENSION_L_T_M" -> dimension|>;

comparisonRecord[left_, right_, dimension_, epsilonOrder_:1] := Module[
  {leftComputed, rightComputed, residual},
  leftComputed = truncateBackground[left];
  rightComputed = truncateBackground[right];
  residual = truncateBackground[leftComputed - rightComputed];
  <|"OPERAND_A" -> leftComputed, "OPERAND_B" -> rightComputed,
    "RESIDUAL_A_MINUS_B" -> residual,
    "TEST_OBJECT" -> relationalObject[leftComputed, rightComputed],
    "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
      (gradeTerms[residual] /. {a_, b_} :> {epsilonOrder, a, b}),
    "DIMENSION_L_T_M" -> dimension|>
];

mapDifference[left_Association, right_Association] := AssociationMap[
  Function[key, mapDifference[left[key], right[key]]], Keys[left]];
mapDifference[left_List, right_List] /; Length[left] == Length[right] :=
  MapThread[mapDifference, {left, right}];
(* Every operand passed here has already been truncated at its construction
   source.  Re-expanding the difference would change no retained order and
   needlessly duplicates the largest operator trees. *)
mapDifference[left_, right_] := left - right;

(* ---------------------------------------------------------------------- *)
(* Coordinates, the inherited constants, and fresh background profiles.  *)
(* ---------------------------------------------------------------------- *)

braneDimension = 3;
spatialCoordinates = {xOne, xTwo, xThree};
materialCoordinates = {capitalXOne, capitalXTwo, capitalXThree};
directions = Range[braneDimension];
faces = {1, -1};
branches = {"LAB_HELD", "MATERIAL_ADVECTED"};
densities = {"RHO4_CONSTANT", "RHOBR_CONSTANT"};
profileSources = {"W_BG", "MU_R_BG"};

(* Inherited constants and the reserved operand have their own blind-engine
   symbols.  Fresh profile symbols below are not aliases of these constants. *)
inheritedConstants = <|
  "c_s0" -> csZero,
  "mu_R" -> muR,
  "rho_br" -> rhoBr,
  "W_0" -> WZero,
  "e_W" -> eWField,
  "rho_m" -> rhoM,
  "v_bulk_normal_0" -> vBulkNormalZero,
  "mu_theta" -> muThetaReserved|>;

(* The three supplied response laws retain independent relaxation times. *)
lambdaAResponse = lambdaAZero/(1 - I frequency tauA);
lambdaVResponse = lambdaVZero/(1 - I frequency tauV);
lambdaXResponse = lambdaXZero/(1 - I frequency tauX);

uVector = Through[{uOne, uTwo, uThree}[
    xOne, xTwo, xThree, time]];
transverseTrialPotentialVector = Through[{transverseTrialPotentialOne,
    transverseTrialPotentialTwo, transverseTrialPotentialThree}[
      xOne, xTwo, xThree, time]];
transverseTestPotentialVector = Through[{transverseTestPotentialOne,
    transverseTestPotentialTwo, transverseTestPotentialThree}[
      xOne, xTwo, xThree, time]];
longitudinalTrialPotentialField = longitudinalTrialPotential[
  xOne, xTwo, xThree, time];
longitudinalTestPotentialField = longitudinalTestPotential[
  xOne, xTwo, xThree, time];
thetaTrialField = thetaTrial[xOne, xTwo, xThree, time];
thetaTestField = thetaTest[xOne, xTwo, xThree, time];
ewTrialField = ewTrial[xOne, xTwo, xThree, time];
ewTestField = ewTest[xOne, xTwo, xThree, time];
virtualUVector = Through[{virtualUOne, virtualUTwo, virtualUThree}[
    xOne, xTwo, xThree, time]];
thetaField = thetaWave[xOne, xTwo, xThree, time];
eWField = eWave[xOne, xTwo, xThree, time];
zetaCenterField = zetaCenter[xOne, xTwo, xThree, time];
virtualThetaField = virtualTheta[xOne, xTwo, xThree, time];
virtualEwField = virtualEw[xOne, xTwo, xThree, time];
virtualZetaCenterField = virtualZetaCenter[xOne, xTwo, xThree, time];

gradient[scalar_] := D[scalar, #] & /@ spatialCoordinates;
divergence[vector_] := Sum[D[vector[[index]],
    spatialCoordinates[[index]]], {index, directions}];
curl[vector_] := {
  D[vector[[3]], xTwo] - D[vector[[2]], xThree],
  D[vector[[1]], xThree] - D[vector[[3]], xOne],
  D[vector[[2]], xOne] - D[vector[[1]], xTwo]};

functionAtCoordinates[expression_] := Function[Evaluate[
  expression /. Thread[Append[spatialCoordinates, time] ->
    {#1, #2, #3, #4}]]];

waveFieldHeadRules[uReplacement_List, thetaReplacement_, ewReplacement_,
    centerReplacement_:0, upperPressureReplacement_:0,
    lowerPressureReplacement_:0] := Join[
  MapThread[#1 -> functionAtCoordinates[#2] &,
    {{uOne, uTwo, uThree}, uReplacement}],
  {thetaWave -> functionAtCoordinates[thetaReplacement],
    eWave -> functionAtCoordinates[ewReplacement],
    zetaCenter -> functionAtCoordinates[centerReplacement],
    pressureUpper -> functionAtCoordinates[upperPressureReplacement],
    pressureLower -> functionAtCoordinates[lowerPressureReplacement]}];

eulerScalar[density_, field_] := D[density, field] -
  Sum[D[D[density, D[field, spatialCoordinates[[index]]]],
    spatialCoordinates[[index]]], {index, directions}];
eulerVector[density_, fields_List] := eulerScalar[density, #] & /@ fields;

inactiveDivergence[vector_] := Inactive[Div][vector, spatialCoordinates];
variationalSource[density_, field_] :=
  D[density, field] - inactiveDivergence[
    D[density, D[field, #]] & /@ spatialCoordinates];

linearVirtualVariation[density_, field_] :=
  Coefficient[density, field] - inactiveDivergence[
    (Coefficient[density, D[field, #]] & /@ spatialCoordinates)];

widthHeads = {widthBase, widthCutOne, widthCutTwo, widthCutThree};
modulusHeads = {modulusBase, modulusCutOne, modulusCutTwo,
  modulusCutThree};

profileHeads[source_String, direction_Integer] := Switch[source,
  "BASE", {widthBase, modulusBase},
  "W_BG", {widthHeads[[direction + 1]], modulusBase},
  "MU_R_BG", {widthBase, modulusHeads[[direction + 1]]}];

widthJetSymbols = {w1JetOne, w1JetTwo, w1JetThree};
modulusJetSymbols = {m1JetOne, m1JetTwo, m1JetThree};
widthValue = WZero (1 + etaBg w1ProfileZero);
modulusValue = muR (1 + etaBg m1ProfileZero);

profileDerivativeRules[head_Symbol, value_, jets_List, scale_,
    cutDirection_Integer] := Module[{firstRules, higherRules},
  firstRules = Table[
    With[{orders = UnitVector[braneDimension, index]},
      Apply[Derivative, orders][head][Sequence @@ spatialCoordinates] ->
        If[index == cutDirection, 0, scale jets[[index]]]],
    {index, directions}];
  higherRules = Flatten[Table[
    If[2 <= i + j + k <= 3,
      {Derivative[i, j, k][head][Sequence @@ spatialCoordinates] -> 0},
      Nothing], {i, 0, 3}, {j, 0, 3}, {k, 0, 3}]];
  Join[firstRules, higherRules,
    {head[Sequence @@ spatialCoordinates] -> value}]
];

profileRules[source_String, direction_Integer] := Module[
  {heads, widthCut, modulusCut, widthJets, modulusJets},
  heads = {widthBase, modulusBase};
  widthCut = If[source === "W_BG", direction, 0];
  modulusCut = If[source === "MU_R_BG", direction, 0];
  widthJets = If[MemberQ[{"PARAMETRIC", "PARAMETRIC_W"}, source],
    MapThread[Times, {{widthFormOne, widthFormTwo, widthFormThree},
      widthJetSymbols}], widthJetSymbols];
  modulusJets = If[MemberQ[{"PARAMETRIC", "PARAMETRIC_MU"}, source],
    MapThread[Times, {{modulusFormOne, modulusFormTwo,
      modulusFormThree}, modulusJetSymbols}], modulusJetSymbols];
  Join[
    profileDerivativeRules[heads[[1]], widthValue, widthJets,
      sigmaW, widthCut],
    profileDerivativeRules[heads[[2]], modulusValue, modulusJets,
      (muR/WZero) sigmaW, modulusCut]]
];

profileJetSymbol["W_BG", {1, 0, 0}] := w1JetOne;
profileJetSymbol["W_BG", {0, 1, 0}] := w1JetTwo;
profileJetSymbol["W_BG", {0, 0, 1}] := w1JetThree;
profileJetSymbol["MU_R_BG", {1, 0, 0}] := m1JetOne;
profileJetSymbol["MU_R_BG", {0, 1, 0}] := m1JetTwo;
profileJetSymbol["MU_R_BG", {0, 0, 1}] := m1JetThree;
profileJetSymbol["W_BG", orders_List] := widthProfileJet @@ orders;
profileJetSymbol["MU_R_BG", orders_List] := modulusProfileJet @@ orders;

profileRulesRetainingGeneratedJets[] := Module[{orders},
  Flatten[Table[
    orders = {i, j, k};
    Which[
      i + j + k == 0,
        {widthBase[Sequence @@ spatialCoordinates] -> widthValue,
          modulusBase[Sequence @@ spatialCoordinates] -> modulusValue},
      i + j + k >= 1,
        {Derivative[i, j, k][widthBase][
            Sequence @@ spatialCoordinates] ->
            sigmaW profileJetSymbol["W_BG", orders]/
              LWidth^(i + j + k - 1),
          Derivative[i, j, k][modulusBase][
            Sequence @@ spatialCoordinates] ->
            (muR/WZero) sigmaW
              profileJetSymbol["MU_R_BG", orders]/
              LWidth^(i + j + k - 1)},
      True, Nothing],
    {i, 0, 3}, {j, 0, 3}, {k, 0, 3}]]
];

applyBackgroundProfileWithGeneratedJets[expression_] :=
  truncateBackground[expression /. profileRulesRetainingGeneratedJets[]];

applyProfile[expression_, source_String:"BASE", direction_Integer:0] :=
  truncateBackground[expression /. profileRules[source, direction]];

backgroundAnchor["LAB_HELD", head_Symbol, wave_] :=
  head[Sequence @@ spatialCoordinates];
backgroundAnchor["MATERIAL_ADVECTED", head_Symbol, wave_] :=
  head[Sequence @@ (spatialCoordinates - wave uVector)];

rhoFourBackground["RHO4_CONSTANT", width_] := rhoBr/WZero;
rhoFourBackground["RHOBR_CONSTANT", width_] := rhoBr/width;
rhoBrBackground["RHO4_CONSTANT", width_] := (rhoBr/WZero) width;
rhoBrBackground["RHOBR_CONSTANT", width_] := rhoBr;
sigmaBackground[density_String, width_] :=
  rhoFourBackground[density, width] width;

(* ---------------------------------------------------------------------- *)
(* Generated variable-coefficient energy quotient representative.        *)
(* ---------------------------------------------------------------------- *)

uniformLabels = {
  "CURL_U_SQUARED", "DIV_U_SQUARED", "THETA_SQUARED",
  "ELOCAL_SQUARED", "THETA_ELOCAL", "GRAD_THETA_SQUARED",
  "GRAD_ELOCAL_SQUARED", "GRAD_THETA_DOT_GRAD_ELOCAL",
  "THETA_DIV_U", "ELOCAL_DIV_U"};
carriedUniformLabels = {"CURL_U_SQUARED", "THETA_SQUARED",
  "THETA_ELOCAL", "ELOCAL_SQUARED", "GRAD_ELOCAL_SQUARED"};

(* The O(3) family is generated from tensor ranks.  Every scalar is one
   vector spurion times two DOF factors, with all slots paired by delta.
   Levi-Civita contractions are generated only for the parity diagnostic. *)
dofFactorSpecifications = {
  <|"NAME" -> "U", "RANK" -> 1, "DIMENSION" -> dimensionU|>,
  <|"NAME" -> "GRADU", "RANK" -> 2,
    "DIMENSION" -> dimensionGradient + dimensionU|>,
  <|"NAME" -> "THETA", "RANK" -> 0,
    "DIMENSION" -> dimensionTheta|>,
  <|"NAME" -> "GRADTHETA", "RANK" -> 1,
    "DIMENSION" -> dimensionGradient + dimensionTheta|>,
  <|"NAME" -> "EW", "RANK" -> 0, "DIMENSION" -> dimensionEw|>,
  <|"NAME" -> "GRADEW", "RANK" -> 1,
    "DIMENSION" -> dimensionGradient + dimensionEw|>};

perfectMatchings[{}] := {{}};
perfectMatchings[slots_List] /; EvenQ[Length[slots]] := Module[
  {first, rest},
  first = First[slots];
  rest = Rest[slots];
  Flatten[Table[
    (Prepend[#, {first, rest[[position]]}] & /@
      perfectMatchings[Delete[rest, position]]),
    {position, Length[rest]}], 1]
];

partnerOfSlot[matching_List, slot_List] := Module[{pair},
  pair = First[Select[matching, MemberQ[#, slot] &]];
  First[DeleteCases[pair, slot]]
];

pairingCode[matching_List] := Module[{slotCode, pairCodes},
  slotCode[{1, slot_}] := "S" <> ToString[slot];
  slotCode[{2, slot_}] := "A" <> ToString[slot];
  slotCode[{3, slot_}] := "B" <> ToString[slot];
  pairCodes = Sort[StringRiffle[Sort[slotCode /@ #], ""] & /@
    matching];
  "K" <> StringRiffle[pairCodes, "_"]
];

invariantLabelSuffix[factors_List, matching_List] := Module[
  {ranks, names, scalarPosition, vectorPosition, tensorPosition,
    scalarName, vectorName, vectorBlock, tensorBlock, spurionPartner},
  ranks = factors[[All, "RANK"]];
  names = factors[[All, "NAME"]];
  Which[
    Sort[ranks] === {0, 1},
      scalarPosition = First[FirstPosition[ranks, 0]];
      vectorPosition = First[FirstPosition[ranks, 1]];
      scalarName = names[[scalarPosition]];
      vectorName = names[[vectorPosition]];
      If[vectorName === "U", vectorName <> "_" <> scalarName,
        scalarName <> "_" <> vectorName],
    Sort[ranks] === {1, 2},
      vectorPosition = First[FirstPosition[ranks, 1]];
      tensorPosition = First[FirstPosition[ranks, 2]];
      vectorName = names[[vectorPosition]];
      vectorBlock = vectorPosition + 1;
      tensorBlock = tensorPosition + 1;
      spurionPartner = partnerOfSlot[matching, {1, 1}];
      Which[
        spurionPartner === {vectorBlock, 1},
          If[vectorName === "U", "U_DIVU",
            "DIVU_" <> vectorName],
        spurionPartner === {tensorBlock, 1},
          "GRADU_" <> vectorName,
        spurionPartner === {tensorBlock, 2},
          "GRADUT_" <> vectorName,
        True, StringRiffle[names, "_"] <> "_" <> pairingCode[matching]],
    True, StringRiffle[names, "_"] <> "_" <> pairingCode[matching]]
];

enumerateO3KroneckerBlueprints[] := Module[
  {factorPairs, candidates},
  factorPairs = Join[Subsets[Range[Length[dofFactorSpecifications]],
    {2}], ({#, #} & /@ Range[Length[dofFactorSpecifications]])];
  factorPairs = Select[factorPairs, EvenQ[1 + Total[
    dofFactorSpecifications[[#, "RANK"]]]] &];
  candidates = Flatten[Table[Module[
      {factors, ranks, slots, matchings},
      factors = dofFactorSpecifications[[factorPair]];
      ranks = Join[{1}, factors[[All, "RANK"]]];
      slots = Flatten[Table[Table[{block, slot},
        {slot, ranks[[block]]}], {block, Length[ranks]}], 1];
      matchings = perfectMatchings[slots];
      Map[Function[matching, <|
        "FACTOR_NAMES" -> factors[[All, "NAME"]],
        "FACTOR_RANKS" -> factors[[All, "RANK"]],
        "FACTOR_DIMENSIONS" -> Join[{dimensionGradient},
          factors[[All, "DIMENSION"]]],
        "PAIRING" -> matching,
        "LABEL_SUFFIX" -> invariantLabelSuffix[factors, matching],
        "STORED_INVARIANT_DIMENSION" ->
          Total[Join[{dimensionGradient},
            factors[[All, "DIMENSION"]]]]|>], matchings]
    ], {factorPair, factorPairs}], 1];
  DeleteDuplicatesBy[candidates,
    {#1["FACTOR_NAMES"], #1["PAIRING"]} &]
];

o3KroneckerBlueprints = enumerateO3KroneckerBlueprints[];

(* This list names only the committed form-ablation slice.  It never sources
   the completed family or its rank. *)
legacyRestrictedLabelSuffixes = {"U_THETA", "U_EW", "U_DIVU",
  "DIVU_GRADTHETA", "GRADU_GRADTHETA", "DIVU_GRADEW",
  "GRADU_GRADEW", "THETA_GRADEW"};
restrictNewInvariantEnumerationToLegacyEight = False;

selectedO3Blueprints[
    restricted_:restrictNewInvariantEnumerationToLegacyEight] :=
  If[TrueQ[restricted],
  Select[o3KroneckerBlueprints,
    MemberQ[legacyRestrictedLabelSuffixes, #1["LABEL_SUFFIX"]] &],
  o3KroneckerBlueprints];

tensorComponent[tensor_, {}] := tensor;
tensorComponent[tensor_, indices_List] := Extract[tensor, indices];
pairPositionForSlot[matching_List, slot_List] :=
  First[FirstPosition[matching, pair_ /; MemberQ[pair, slot]]];

kroneckerFullContraction[tensors_List, ranks_List,
    matching_List] := Total[Map[Function[assignment,
  Times @@ Table[Module[{indices},
    indices = Table[assignment[[pairPositionForSlot[matching,
      {block, slot}]]], {slot, ranks[[block]]}];
    tensorComponent[tensors[[block]], indices]],
    {block, Length[tensors]}]],
  Tuples[directions, Length[matching]]]];

instantiateO3KroneckerFamily[spurion_List, u_List, theta_, ew_,
    restricted_:restrictNewInvariantEnumerationToLegacyEight] := Module[
  {tensorValues, blueprints, records},
  tensorValues = <|
    "U" -> u,
    "GRADU" -> Table[D[u[[i]], spatialCoordinates[[j]]],
      {i, directions}, {j, directions}],
    "THETA" -> theta,
    "GRADTHETA" -> gradient[theta],
    "EW" -> ew,
    "GRADEW" -> gradient[ew]|>;
  blueprints = selectedO3Blueprints[restricted];
  records = Map[Function[blueprint, Append[blueprint,
    "EXPRESSION" -> kroneckerFullContraction[
      Join[{spurion}, Lookup[tensorValues,
        blueprint["FACTOR_NAMES"]]],
      Join[{1}, blueprint["FACTOR_RANKS"]],
      blueprint["PAIRING"]]]], blueprints];
  DeleteDuplicatesBy[records, Expand[#1["EXPRESSION"]] &]
];

newInvariantLabels[source_String,
    restricted_:restrictNewInvariantEnumerationToLegacyEight] :=
  Module[{prefix},
  prefix = If[source === "W_BG", "WJET", "MUJET"];
  prefix <> "_" <> #1["LABEL_SUFFIX"] & /@
    selectedO3Blueprints[restricted]
];

newCoefficientSymbol[source_String, label_String] := Symbol[
  "gamma" <> StringReplace[label, "_" -> ""]];
newCoefficientStandardName[source_String, label_String] := Module[
  {sourceCode, suffix},
  sourceCode = If[source === "W_BG", "W", "mu"];
  suffix = StringDrop[label, StringLength[
    If[source === "W_BG", "WJET_", "MUJET_"]]];
  "gamma_" <> sourceCode <> "_" <> ToLowerCase[suffix]
];

newCoefficientSymbols = Association[Table[source ->
  (newCoefficientSymbol[source, #] & /@ newInvariantLabels[source, False]),
  {source, profileSources}]];
newCoefficientNames = Association[Table[source ->
  (newCoefficientStandardName[source, #] & /@
    newInvariantLabels[source, False]),
  {source, profileSources}]];
newCoefficientSymbolMaps = Association[Table[source ->
  AssociationThread[newInvariantLabels[source, False],
    newCoefficientSymbols[source]], {source, profileSources}]];
newCoefficientNameMaps = Association[Table[source ->
  AssociationThread[newInvariantLabels[source, False],
    newCoefficientNames[source]], {source, profileSources}]];

newInvariantExpressions[spurion_List, u_List, theta_, ew_,
    restricted_:restrictNewInvariantEnumerationToLegacyEight] :=
  instantiateO3KroneckerFamily[
  spurion, u, theta, ew, restricted][[All, "EXPRESSION"]];

quadraticCoefficient[expression_] :=
  D[expression, {waveOrder, 2}]/2 /. waveOrder -> 0;

constructEnergyData[branch_String, widthHead_Symbol,
    modulusHead_Symbol,
    restricted_:restrictNewInvariantEnumerationToLegacyEight,
    useThicknessMap_:True] := Module[
  {wave, anchoredWidth, anchoredModulus, uWave, thetaWave, ewWave,
    localEw, newEw, gradTheta, gradLocalEw, divU, curlU,
    uniformInvariants,
    uniformCoefficients, uniformFactors, uniformTerms,
    widthSpurion, modulusSpurion, newRecords, allRecords},
  wave = waveOrder;
  anchoredWidth = backgroundAnchor[branch, widthHead, wave];
  anchoredModulus = backgroundAnchor[branch, modulusHead, wave];
  uWave = wave uVector;
  thetaWave = wave thetaField;
  ewWave = wave eWField;
  localEw = (WZero/anchoredWidth) ewWave;
  newEw = If[TrueQ[useThicknessMap], localEw, ewWave];
  gradTheta = gradient[thetaWave];
  gradLocalEw = gradient[localEw];
  divU = divergence[uWave];
  curlU = curl[uWave];

  uniformInvariants = {
    Dot[curlU, curlU],
    divU^2,
    thetaWave^2,
    localEw^2,
    thetaWave localEw,
    Dot[gradTheta, gradTheta],
    Dot[gradLocalEw, gradLocalEw],
    Dot[gradTheta, gradLocalEw],
    thetaWave divU,
    localEw divU};

  uniformCoefficients = {
    anchoredModulus, modulusDivU, bRho anchoredWidth,
    kW anchoredWidth^2, cCross anchoredWidth,
    gradientThetaCoefficient,
    kappaW anchoredWidth^4,
    gradientThetaEwCoefficient,
    thetaDivUCoefficient,
    ewDivUCoefficient anchoredWidth/WZero};
  uniformFactors = {1/2, 1/2, 1/2, 1/2, 1, 1/2, 1/2, 1, 1, 1};

  (* The carried gradient term is constructed with the physical thickness
     anchoredWidth localEw before the exact map simplifies it. *)
  uniformInvariants[[7]] = anchoredWidth^(-2) Dot[
    gradient[anchoredWidth localEw],
    gradient[anchoredWidth localEw]];

  uniformTerms = MapThread[
    Function[{label, invariant, coefficient, factor}, <|
      "LABEL" -> label,
      "INVARIANT" -> quadraticCoefficient[invariant],
      "COEFFICIENT" -> coefficient /. waveOrder -> 0,
      "ENERGY_TERM" -> quadraticCoefficient[
        factor coefficient invariant],
      "CARRIED_UNIFORM_INPUT" -> MemberQ[carriedUniformLabels, label],
      "ORIGIN" -> UniformQuotientRepresentative|>],
    {uniformLabels, uniformInvariants, uniformCoefficients,
      uniformFactors}];

  widthSpurion = gradient[anchoredWidth]/WZero;
  modulusSpurion = gradient[anchoredModulus]/muR;
  newRecords = Flatten[Table[Module[
      {source, spurion, expressions, labels, coefficients, names},
      source = sourceName;
      spurion = If[source === "W_BG", widthSpurion, modulusSpurion];
      expressions = newInvariantExpressions[spurion, uWave,
        thetaWave, newEw, restricted];
      labels = newInvariantLabels[source, restricted];
      coefficients = Lookup[newCoefficientSymbolMaps[source], labels];
      names = Lookup[newCoefficientNameMaps[source], labels];
      MapThread[Function[{label, invariant, coefficient, standardName},
        <|"LABEL" -> label,
          "INVARIANT" -> quadraticCoefficient[invariant],
          "COEFFICIENT" -> coefficient,
          "COEFFICIENT_STANDARD_NAME" -> standardName,
          "ENERGY_TERM" -> quadraticCoefficient[
            coefficient invariant],
          "CARRIED_UNIFORM_INPUT" -> False,
          "ORIGIN" -> FirstBackgroundJetRepresentative|>],
        {labels, expressions, coefficients, names}]
    ], {sourceName, profileSources}], 1];
  allRecords = Join[uniformTerms, newRecords];
  <|"RECORDS" -> allRecords,
    "ENERGY" -> Total[allRecords[[All, "ENERGY_TERM"]]],
    "ANCHORED_WIDTH" -> anchoredWidth /. waveOrder -> 0,
    "ANCHORED_MODULUS" -> anchoredModulus /. waveOrder -> 0|>
];

constructFullFieldBackgroundEnergy[branch_String, widthHead_Symbol,
    modulusHead_Symbol] := Module[
  {order, anchoredWidth, anchoredModulus, uVariation, thetaVariation,
    ewVariation, localEw, fullWidth, gradTheta, gradFullEw, divU,
    curlU, uniformInvariants, uniformCoefficients, uniformFactors,
    uniformTerms, widthSpurion, modulusSpurion, newTerms},
  order = backgroundOrder;
  anchoredWidth = backgroundAnchor[branch, widthHead, order];
  anchoredModulus = backgroundAnchor[branch, modulusHead, order];
  uVariation = order uVector;
  thetaVariation = order thetaField;
  ewVariation = order eWField;
  localEw = (WZero/anchoredWidth) ewVariation;
  fullWidth = anchoredWidth + WZero ewVariation;
  gradTheta = gradient[thetaVariation];
  gradFullEw = anchoredWidth^(-1) gradient[fullWidth];
  divU = divergence[uVariation];
  curlU = curl[uVariation];

  uniformInvariants = {
    Dot[curlU, curlU],
    divU^2,
    thetaVariation^2,
    localEw^2,
    thetaVariation localEw,
    Dot[gradTheta, gradTheta],
    anchoredWidth^(-2) Dot[gradient[fullWidth], gradient[fullWidth]],
    Dot[gradTheta, gradFullEw],
    thetaVariation divU,
    localEw divU};
  uniformCoefficients = {
    anchoredModulus, modulusDivU, bRho anchoredWidth,
    kW anchoredWidth^2, cCross anchoredWidth,
    gradientThetaCoefficient, kappaW anchoredWidth^4,
    gradientThetaEwCoefficient, thetaDivUCoefficient,
    ewDivUCoefficient anchoredWidth/WZero};
  uniformFactors = {1/2, 1/2, 1/2, 1/2, 1, 1/2, 1/2, 1, 1, 1};
  uniformTerms = MapThread[#1 #2 #3 &,
    {uniformFactors, uniformCoefficients, uniformInvariants}];

  widthSpurion = gradient[anchoredWidth]/WZero;
  modulusSpurion = gradient[anchoredModulus]/muR;
  newTerms = Flatten[Table[MapThread[Times,
      {Lookup[newCoefficientSymbolMaps[sourceName],
          newInvariantLabels[sourceName]],
        newInvariantExpressions[
          If[sourceName === "W_BG", widthSpurion, modulusSpurion],
          uVariation, thetaVariation, localEw]}],
    {sourceName, profileSources}], 1];
  Total[Join[uniformTerms, newTerms][[basisRepresentativeIndices]]]
];

fieldAlgebraRules = Join[
  Thread[uVector -> {uAlgebraOne, uAlgebraTwo, uAlgebraThree}],
  Flatten[Table[
    D[uVector[[i]], spatialCoordinates[[j]]] ->
      uAlgebraGradient[i, j], {i, directions}, {j, directions}]],
  {thetaField -> thetaAlgebra, eWField -> ewAlgebra},
  Table[D[thetaField, spatialCoordinates[[i]]] ->
    thetaAlgebraGradient[i], {i, directions}],
  Table[D[eWField, spatialCoordinates[[i]]] ->
    ewAlgebraGradient[i], {i, directions}]];

rankVariables = Join[
  {uAlgebraOne, uAlgebraTwo, uAlgebraThree, thetaAlgebra, ewAlgebra,
    etaBg, sigmaW, w1ProfileZero, m1ProfileZero},
  Flatten[Table[uAlgebraGradient[i, j], {i, directions},
    {j, directions}]],
  Table[thetaAlgebraGradient[i], {i, directions}],
  Table[ewAlgebraGradient[i], {i, directions}],
  widthJetSymbols, modulusJetSymbols];

rankMatrix[expressions_List] := Module[
  {polynomials, ruleLists, monomials},
  polynomials = Expand[(applyProfile[#, "BASE", 0] /.
      fieldAlgebraRules)] & /@ expressions;
  ruleLists = CoefficientRules[#, rankVariables] & /@ polynomials;
  monomials = DeleteDuplicates[Flatten[Keys /@ ruleLists, 1]];
  Table[Lookup[Association[ruleLists[[row]]], Key[monomial], 0],
    {row, Length[ruleLists]}, {monomial, monomials}]
];

polynomialFeatureMatrix[expressionVectors_List, variables_List] := Module[
  {ruleLists, featureKeys},
  ruleLists = Map[Function[expressionVector, Flatten[MapIndexed[
    Function[{expression, position},
      ({First[position], First[#]} -> Last[#] & /@
        CoefficientRules[Expand[expression], variables])],
    Flatten[{expressionVector}]], 1]], expressionVectors];
  featureKeys = DeleteDuplicates[Flatten[Keys /@ ruleLists, 1]];
  Table[Lookup[Association[ruleLists[[row]]], Key[feature], 0],
    {row, Length[ruleLists]}, {feature, featureKeys}]
];

independentRepresentativeIndices[expressions_List] := Module[
  {matrix, selected, current, oldRank, newRank},
  matrix = rankMatrix[expressions];
  selected = {};
  current = {};
  oldRank = 0;
  Do[
    newRank = MatrixRank[Append[current, matrix[[index]]]];
    If[newRank > oldRank,
      AppendTo[selected, index];
      AppendTo[current, matrix[[index]]];
      oldRank = newRank],
    {index, Length[matrix]}];
  selected
];

baseEnergyDataByBranch = Association[Table[branch ->
  constructEnergyData[branch, widthBase, modulusBase],
  {branch, branches}]];
basisRepresentativeIndicesByBranch = AssociationMap[
  independentRepresentativeIndices[
    baseEnergyDataByBranch[#]["RECORDS"][[All, "INVARIANT"]]] &,
  branches];

anchoringExpressions = AssociationMap[
  baseEnergyDataByBranch[#]["RECORDS"][[All, "INVARIANT"]] &,
  branches];
anchoringRankMatrices = AssociationMap[rankMatrix[
  anchoringExpressions[#]] &, branches];
anchoringCommonMatrix = rankMatrix[Join[
  anchoringExpressions["LAB_HELD"],
  anchoringExpressions["MATERIAL_ADVECTED"]]];
anchoringRecordCount = Length[anchoringExpressions["LAB_HELD"]];
anchoringLabSelectedRows = anchoringCommonMatrix[[
  basisRepresentativeIndicesByBranch["LAB_HELD"]]];
anchoringMaterialSelectedRows = anchoringCommonMatrix[[
  anchoringRecordCount +
    basisRepresentativeIndicesByBranch["MATERIAL_ADVECTED"]]];
anchoringSpanCombinedRank = MatrixRank[Join[
  anchoringLabSelectedRows, anchoringMaterialSelectedRows]];
anchoringRankPayload = <|
  "PER_ANCHORING" -> AssociationMap[<|
    "RANK_OPERAND" -> MatrixRank[anchoringRankMatrices[#]],
    "RANK_SELECTED_INDICES" -> basisRepresentativeIndicesByBranch[#]|> &,
    branches],
  "COMBINED_SELECTED_SPAN_RANK" -> anchoringSpanCombinedRank,
  "SPAN_EQUIVALENCE_RESIDUAL" -> <|
    "COMBINED_MINUS_LAB" -> anchoringSpanCombinedRank -
      MatrixRank[anchoringLabSelectedRows],
    "COMBINED_MINUS_MATERIAL" -> anchoringSpanCombinedRank -
      MatrixRank[anchoringMaterialSelectedRows]|>,
  "DIMENSION_L_T_M" -> dimensionZero,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>;
emitShared["ENERGY_BASIS_ANCHORING_RANK", anchoringRankPayload];

(* Common indices enter downstream constructors only after the independent
   anchoring ranks and span residual have been emitted. *)
basisRepresentativeIndices =
  basisRepresentativeIndicesByBranch["LAB_HELD"];
baseEnergyData = baseEnergyDataByBranch["LAB_HELD"];
basisLabels = baseEnergyData["RECORDS"][[basisRepresentativeIndices,
  "LABEL"]];
newBasisLabels = Select[basisLabels,
  !MemberQ[uniformLabels, #] &];
Clear[anchoringExpressions, anchoringRankMatrices, anchoringCommonMatrix,
  anchoringLabSelectedRows, anchoringMaterialSelectedRows,
  anchoringRankPayload];

epsilonPseudoscalarCandidates[branch_String, source_String] := Module[
  {wave, anchoredWidth, anchoredModulus, uWave, thetaWave, ewWave,
    localEw, spurion, vectorNames, vectorValues, vectorPairs},
  wave = waveOrder;
  anchoredWidth = backgroundAnchor[branch, widthBase, wave];
  anchoredModulus = backgroundAnchor[branch, modulusBase, wave];
  uWave = wave uVector;
  thetaWave = wave thetaField;
  ewWave = wave eWField;
  localEw = (WZero/anchoredWidth) ewWave;
  spurion = If[source === "W_BG", gradient[anchoredWidth]/WZero,
    gradient[anchoredModulus]/muR];
  vectorNames = Cases[dofFactorSpecifications,
    specification_ /; specification["RANK"] == 1 :>
      specification["NAME"]];
  vectorValues = <|"U" -> uWave,
    "GRADTHETA" -> gradient[thetaWave],
    "GRADEW" -> gradient[localEw]|>;
  vectorPairs = Subsets[vectorNames, {2}];
  Map[Function[pair, <|"FACTOR_NAMES" -> pair,
    "EXPRESSION" -> quadraticCoefficient[Sum[
      Signature[{i, j, k}] spurion[[i]]
        vectorValues[pair[[1]]][[j]] vectorValues[pair[[2]]][[k]],
      {i, directions}, {j, directions}, {k, directions}]]|>],
    vectorPairs]
];

(* ---------------------------------------------------------------------- *)
(* Independent virtual-constraint, evolution, and tilted-face routes.    *)
(* ---------------------------------------------------------------------- *)

jacobianLinear[vector_List] := 1 + waveOrder divergence[vector];

virtualConstraintSource[route_String, branch_String, density_String,
    widthHead_Symbol, corrupted_:False] := Module[
  {wave, virtualMap, widthAtMaterial, sigmaAtMaterial, alphaAtMaterial,
    exactMaterialMass, directSigma, directMappedSigma},
  wave = waveOrder;
  virtualMap = spatialCoordinates + wave virtualUVector;
  widthAtMaterial = widthHead[Sequence @@ materialCoordinates];

  If[route === "MATERIAL",
    sigmaAtMaterial = Switch[branch,
      "LAB_HELD", sigmaBackground[density,
        widthHead[Sequence @@ virtualMap]],
      "MATERIAL_ADVECTED", sigmaBackground[density,
        widthHead[Sequence @@ spatialCoordinates]]];
    alphaAtMaterial = WZero/Switch[branch,
      "LAB_HELD", widthHead[Sequence @@ virtualMap],
      "MATERIAL_ADVECTED", widthHead[Sequence @@ spatialCoordinates]];
    exactMaterialMass = sigmaAtMaterial
      (1 + wave virtualThetaField + wave alphaAtMaterial virtualEwField)
      jacobianLinear[virtualUVector];
    Return[(D[exactMaterialMass, wave] /. wave -> 0)/
      sigmaBackground[density,
        widthHead[Sequence @@ spatialCoordinates]]]
  ];

  directSigma = Switch[branch,
    "LAB_HELD", sigmaBackground[density,
      widthHead[Sequence @@ spatialCoordinates]],
    "MATERIAL_ADVECTED", sigmaBackground[density,
      widthHead[Sequence @@ (spatialCoordinates - wave virtualUVector)]]];
  directSigma = directSigma (1 + wave virtualThetaField +
    wave (WZero/widthHead[Sequence @@ spatialCoordinates])
      virtualEwField);
  directMappedSigma = If[TrueQ[corrupted], directSigma,
    directSigma /. Thread[spatialCoordinates -> virtualMap]];
  (D[directMappedSigma jacobianLinear[virtualUVector], wave] /.
      wave -> 0)/sigmaBackground[density,
      widthHead[Sequence @@ spatialCoordinates]]
];

evolutionSource[route_String, branch_String, density_String,
    widthHead_Symbol] := Module[
  {wave, sigmaZero, alpha, sigmaEuler, accumulation, advection},
  wave = waveOrder;
  sigmaZero = sigmaBackground[density,
    widthHead[Sequence @@ spatialCoordinates]];
  alpha = WZero/widthHead[Sequence @@ spatialCoordinates];
  sigmaEuler = Switch[branch,
    "LAB_HELD", sigmaZero,
    "MATERIAL_ADVECTED", sigmaBackground[density,
      widthHead[Sequence @@ (spatialCoordinates - wave uVector)]]];
  sigmaEuler = sigmaEuler (1 + wave thetaField + wave alpha eWField);
  If[route === "MATERIAL",
    sigmaEuler = Switch[branch,
      "LAB_HELD", sigmaBackground[density,
        widthHead[Sequence @@ (spatialCoordinates + wave uVector)]],
      "MATERIAL_ADVECTED", sigmaZero] (1 + wave thetaField +
        wave alpha eWField)];
  accumulation = D[sigmaEuler, time];
  advection = divergence[sigmaEuler wave D[uVector, time]];
  <|"ACCUMULATION" -> (D[accumulation, wave] /. wave -> 0),
    "ADVECTIVE" -> (D[advection, wave] /. wave -> 0)|>
];

faceLabel[1] := "UPPER";
faceLabel[-1] := "LOWER";

pressureField[1] := pressureUpper[xOne, xTwo, xThree, time];
pressureField[-1] := pressureLower[xOne, xTwo, xThree, time];

faceSources[route_String, branch_String, sign_Integer,
    widthHead_Symbol, muTheta_, rhoBrValue_] := Module[
  {wave, virtualWave, zeta, virtualZeta, graphHeight, graphNormal,
    graphMeasure, graphVelocity, materialMap, flatteningDenominator,
    parametricHeight,
    parametricVelocity, normalVelocity, virtualGraphHeight,
    virtualNormalDisplacement, affinity, flux, tractionPressure,
    virtualWork},
  wave = waveOrder;
  virtualWave = virtualOrder;
  zeta = zetaCenterField + sign WZero eWField/2;
  virtualZeta = virtualZetaCenterField +
    sign WZero virtualEwField/2;
  graphHeight = Switch[branch,
    "LAB_HELD", sign widthHead[Sequence @@ spatialCoordinates]/2 +
      wave zeta,
    "MATERIAL_ADVECTED", sign widthHead[Sequence @@
      (spatialCoordinates - wave uVector)]/2 + wave zeta];
  graphNormal = sign Join[-gradient[graphHeight], {1}]/
    Sqrt[1 + Dot[gradient[graphHeight], gradient[graphHeight]]];
  graphMeasure = Sqrt[1 + Dot[gradient[graphHeight],
    gradient[graphHeight]]];

  If[route === "EULERIAN",
    normalVelocity = D[sign D[graphHeight, time]/graphMeasure,
      wave] /. wave -> 0;
    virtualGraphHeight = Switch[branch,
      "LAB_HELD", sign widthHead[Sequence @@ spatialCoordinates]/2 +
        virtualWave virtualZeta,
      "MATERIAL_ADVECTED", sign widthHead[Sequence @@
        (spatialCoordinates - virtualWave virtualUVector)]/2 +
        virtualWave virtualZeta];
    virtualNormalDisplacement = D[
      sign (virtualGraphHeight - sign
        widthHead[Sequence @@ spatialCoordinates]/2)/
        (graphMeasure /. wave -> 0), virtualWave] /. virtualWave -> 0,

    materialMap = spatialCoordinates + wave uVector;
    (* This denominator is the exact material-coordinate flattening
       denominator W_bg + delta W; the numerator is w-zeta_c. *)
    flatteningDenominator = Switch[branch,
      "LAB_HELD", widthHead[Sequence @@ materialMap],
      "MATERIAL_ADVECTED", widthHead[Sequence @@ spatialCoordinates]] +
      wave WZero eWField;
    parametricHeight = wave zetaCenterField +
      sign flatteningDenominator/2;
    parametricVelocity = Join[D[materialMap, time],
      {D[parametricHeight, time]}];
    normalVelocity = D[(graphNormal . parametricVelocity), wave] /.
      wave -> 0;
    materialMap = spatialCoordinates + virtualWave virtualUVector;
    flatteningDenominator = Switch[branch,
      "LAB_HELD", widthHead[Sequence @@ materialMap],
      "MATERIAL_ADVECTED", widthHead[Sequence @@ spatialCoordinates]] +
      virtualWave WZero virtualEwField;
    parametricHeight = virtualWave virtualZetaCenterField +
      sign flatteningDenominator/2;
    virtualNormalDisplacement = D[(graphNormal /. wave -> 0) .
      Join[materialMap, {parametricHeight}], virtualWave] /.
        virtualWave -> 0
  ];

  affinity = muTheta/rhoBrValue - pressureField[sign]/rhoM;
  flux = lambdaAResponse affinity + lambdaVResponse normalVelocity;
  tractionPressure = pressureField[sign] + lambdaXResponse affinity;
  virtualWork = -graphMeasure /. wave -> 0;
  virtualWork = virtualWork tractionPressure virtualNormalDisplacement;
  <|"NORMAL" -> (graphNormal /. wave -> 0),
    "MEASURE" -> (graphMeasure /. wave -> 0),
    "NORMAL_VELOCITY" -> normalVelocity,
    "AFFINITY" -> affinity,
    "RELATIVE_FLUX" -> flux,
    "TRACTION_PRESSURE" -> tractionPressure,
    "VIRTUAL_NORMAL_DISPLACEMENT" -> virtualNormalDisplacement,
    "VIRTUAL_WORK" -> virtualWork|>
];

projectedFaceFlux[faceAssociation_Association] := Total[Values[
  Map[Function[face, face["MEASURE"] face["RELATIVE_FLUX"]],
    faceAssociation]]];

(* ---------------------------------------------------------------------- *)
(* Slab operator constructor.                                            *)
(* ---------------------------------------------------------------------- *)

constrainedRows[energy_, constraint_] := Module[
  {muTheta, ewDerivative, uDerivative, thetaRule, variedDensity,
    uRows, ewRow},
  muTheta = variationalSource[energy, thetaField];
  ewDerivative = variationalSource[energy, eWField];
  uDerivative = variationalSource[energy, #] & /@ uVector;
  thetaRule = -(constraint -
      Coefficient[constraint, virtualThetaField] virtualThetaField)/
      Coefficient[constraint, virtualThetaField];
  variedDensity = muTheta thetaRule + ewDerivative virtualEwField +
    Dot[uDerivative, virtualUVector];
  uRows = linearVirtualVariation[variedDensity, #] & /@ virtualUVector;
  ewRow = linearVirtualVariation[variedDensity, virtualEwField];
  <|"MU_THETA" -> muTheta, "U_INTERNAL" -> uRows,
    "EW_INTERNAL" -> ewRow, "THETA_VARIATION_RULE" -> thetaRule|>
];

faceGeneralizedRows[faceAssociation_Association] := Module[
  {work, uRows, ewRow, centerRow},
  work = Total[faceAssociation[[All, "VIRTUAL_WORK"]]];
  uRows = -(linearVirtualVariation[work, #] & /@ virtualUVector);
  ewRow = -linearVirtualVariation[work, virtualEwField];
  centerRow = -linearVirtualVariation[work, virtualZetaCenterField];
  <|"U_FACE" -> uRows, "EW_FACE" -> ewRow,
    "CENTER_FACE" -> centerRow|>
];

rawModel[route_String, branch_String, density_String,
    source_String:"BASE", direction_Integer:0, corrupted_:False,
    formOnly_:False] := Module[
  {heads, energyData, selectedRecords, energy, constraint, rows,
    width, rhoBrValue, facesData, faceRows, evolutionTerms,
    kineticU, kineticEw, operator, origins, divergenceSource},
  heads = profileHeads[source, direction];
  energyData = constructEnergyData[branch, heads[[1]], heads[[2]]];
  selectedRecords = energyData["RECORDS"][[basisRepresentativeIndices]];
  energy = Total[selectedRecords[[All, "ENERGY_TERM"]]];
  constraint = virtualConstraintSource[route, branch, density,
    heads[[1]], corrupted];
  rows = constrainedRows[energy, constraint];
  width = heads[[1]][Sequence @@ spatialCoordinates];
  rhoBrValue = rhoBrBackground[density, width];
  facesData = Association[Table[faceLabel[sign] ->
    faceSources[route, branch, sign, heads[[1]], rows["MU_THETA"],
      rhoBrValue], {sign, faces}]];
  faceRows = faceGeneralizedRows[facesData];
  evolutionTerms = evolutionSource[route, branch, density, heads[[1]]];
  kineticU = rhoBrValue D[uVector, {time, 2}];
  kineticEw = muW WZero^2 D[eWField, {time, 2}];
  operator = <|
    "U_MOMENTUM_ROWS" -> (kineticU + rows["U_INTERNAL"] +
      faceRows["U_FACE"]),
    "THICKNESS_ROW" -> (kineticEw + rows["EW_INTERNAL"] +
      faceRows["EW_FACE"]),
    "MASS_EVOLUTION_ROW" -> (
      evolutionTerms["ACCUMULATION"] + evolutionTerms["ADVECTIVE"] +
      projectedFaceFlux[facesData]),
    "CENTER_FACE_GENERALIZED_ROW" -> faceRows["CENTER_FACE"]|>;
  If[TrueQ[formOnly], Return[<|
    "SELECTED_RECORDS" -> selectedRecords,
    "ENERGY" -> energy, "OPERATOR" -> operator|>]];
  origins = <|
    "KINETIC" -> <|"U_MOMENTUM_ROWS" -> kineticU,
      "THICKNESS_ROW" -> kineticEw|>,
    "BULK_ENERGY" -> <|"U_MOMENTUM_ROWS" -> rows["U_INTERNAL"],
      "THICKNESS_ROW" -> rows["EW_INTERNAL"],
      "ENERGY_BASIS_LABELS" -> selectedRecords[[All, "LABEL"]]|>,
    "FACE" -> <|"U_MOMENTUM_ROWS" -> faceRows["U_FACE"],
      "THICKNESS_ROW" -> faceRows["EW_FACE"],
      "CENTER_FACE_GENERALIZED_ROW" -> faceRows["CENTER_FACE"]|>,
    "FLUX" -> projectedFaceFlux[facesData],
    "ADVECTIVE" -> evolutionTerms["ADVECTIVE"],
    "ACCUMULATION" -> evolutionTerms["ACCUMULATION"]|>;
  divergenceSource = <|
    "MU_THETA" -> variationalSource[energy, thetaField],
    "U_ENERGY_VARIATIONS" ->
      (variationalSource[energy, #] & /@ uVector),
    "EW_ENERGY_VARIATION" -> variationalSource[energy, eWField]|>;
  <|"SELECTED_RECORDS" -> selectedRecords,
    "ENERGY" -> energy, "CONSTRAINT" -> constraint,
    "THETA_VARIATION_RULE" -> rows["THETA_VARIATION_RULE"],
    "MU_THETA" -> rows["MU_THETA"], "OPERATOR" -> operator,
    "ORIGINS" -> origins, "FACE_SUBSTRATE" -> facesData,
    "DIVERGENCE_FORM_SOURCE" -> divergenceSource|>
];

processModel[model_Association, source_String:"BASE",
    direction_Integer:0] := <|
  "SELECTED_RECORDS" -> applyProfile[model["SELECTED_RECORDS"],
    source, direction],
  "ENERGY" -> applyProfile[model["ENERGY"], source, direction],
  "CONSTRAINT" -> applyProfile[model["CONSTRAINT"], source, direction],
  "THETA_VARIATION_RULE" -> applyProfile[
    model["THETA_VARIATION_RULE"], source, direction],
  "MU_THETA" -> applyProfile[model["MU_THETA"], source, direction],
  "OPERATOR" -> applyProfile[model["OPERATOR"], source, direction],
  "ORIGINS" -> applyProfile[model["ORIGINS"], source, direction],
  "FACE_SUBSTRATE" -> applyProfile[model["FACE_SUBSTRATE"],
    source, direction],
  "DIVERGENCE_FORM_SOURCE" -> model["DIVERGENCE_FORM_SOURCE"]|>;

processFormModel[model_Association, source_String:"BASE",
    direction_Integer:0] := <|
  "SELECTED_RECORDS" -> applyProfile[model["SELECTED_RECORDS"],
    source, direction],
  "ENERGY" -> applyProfile[model["ENERGY"], source, direction],
  "OPERATOR" -> applyProfile[model["OPERATOR"], source, direction]|>;

evaluatedModel[route_String, branch_String, density_String,
    source_String:"BASE", direction_Integer:0, corrupted_:False,
    formOnly_:False] := Module[
  {energyData, selectedRecords, energy, constraint, rows, width,
    rhoBrValue, facesData, faceRows, evolutionTerms, kineticU,
    kineticEw, operator, origins, divergenceSource},
  energyData = constructEnergyData[branch, widthBase, modulusBase];
  selectedRecords = applyProfile[
    energyData["RECORDS"][[basisRepresentativeIndices]],
    source, direction];
  energy = Total[selectedRecords[[All, "ENERGY_TERM"]]];
  constraint = applyProfile[virtualConstraintSource[route, branch,
    density, widthBase, corrupted], source, direction];
  rows = constrainedRows[energy, constraint];
  width = applyProfile[widthBase[Sequence @@ spatialCoordinates],
    source, direction];
  rhoBrValue = applyProfile[rhoBrBackground[density,
    widthBase[Sequence @@ spatialCoordinates]], source, direction];
  facesData = Association[Table[faceLabel[sign] -> applyProfile[
    faceSources[route, branch, sign, widthBase, rows["MU_THETA"],
      rhoBrValue], source, direction], {sign, faces}]];
  faceRows = faceGeneralizedRows[facesData];
  evolutionTerms = applyProfile[evolutionSource[route, branch, density,
    widthBase], source, direction];
  kineticU = rhoBrValue D[uVector, {time, 2}];
  kineticEw = muW WZero^2 D[eWField, {time, 2}];
  operator = <|
    "U_MOMENTUM_ROWS" -> (kineticU + rows["U_INTERNAL"] +
      faceRows["U_FACE"]),
    "THICKNESS_ROW" -> (kineticEw + rows["EW_INTERNAL"] +
      faceRows["EW_FACE"]),
    "MASS_EVOLUTION_ROW" -> (evolutionTerms["ACCUMULATION"] +
      evolutionTerms["ADVECTIVE"] +
      projectedFaceFlux[facesData]),
    "CENTER_FACE_GENERALIZED_ROW" -> faceRows["CENTER_FACE"]|>;
  If[TrueQ[formOnly], Return[<|
    "SELECTED_RECORDS" -> selectedRecords,
    "ENERGY" -> energy, "OPERATOR" -> operator|>]];
  origins = <|
    "KINETIC" -> <|"U_MOMENTUM_ROWS" -> kineticU,
      "THICKNESS_ROW" -> kineticEw|>,
    "BULK_ENERGY" -> <|"U_MOMENTUM_ROWS" -> rows["U_INTERNAL"],
      "THICKNESS_ROW" -> rows["EW_INTERNAL"],
      "ENERGY_BASIS_LABELS" -> selectedRecords[[All, "LABEL"]]|>,
    "FACE" -> <|"U_MOMENTUM_ROWS" -> faceRows["U_FACE"],
      "THICKNESS_ROW" -> faceRows["EW_FACE"],
      "CENTER_FACE_GENERALIZED_ROW" -> faceRows["CENTER_FACE"]|>,
    "FLUX" -> projectedFaceFlux[facesData],
    "ADVECTIVE" -> evolutionTerms["ADVECTIVE"],
    "ACCUMULATION" -> evolutionTerms["ACCUMULATION"]|>;
  divergenceSource = <|
    "MU_THETA" -> variationalSource[energy, thetaField],
    "U_ENERGY_VARIATIONS" ->
      (variationalSource[energy, #] & /@ uVector),
    "EW_ENERGY_VARIATION" -> variationalSource[energy, eWField]|>;
  <|"SELECTED_RECORDS" -> selectedRecords, "ENERGY" -> energy,
    "CONSTRAINT" -> constraint,
    "THETA_VARIATION_RULE" -> rows["THETA_VARIATION_RULE"],
    "MU_THETA" -> rows["MU_THETA"], "OPERATOR" -> operator,
    "ORIGINS" -> origins, "FACE_SUBSTRATE" -> facesData,
    "DIVERGENCE_FORM_SOURCE" -> divergenceSource,
    "CONSTRUCTION_METADATA" -> <|"ROUTE" -> route,
      "BRANCH" -> branch, "DENSITY" -> density,
      "SOURCE" -> source, "DIRECTION" -> direction,
      "CORRUPTED" -> corrupted|>|>
];

activateSpatialDivergences[expression_] := FixedPoint[ReplaceAll[#,
  divergenceObject : Inactive[Div][vector_List, coordinates_List] :>
    If[SameQ[coordinates, spatialCoordinates],
      Sum[D[vector[[index]], coordinates[[index]]],
        {index, Length[coordinates]}], divergenceObject]] &, expression];

elRowVector[energy_, constraint_] := Module[{rows},
  rows = constrainedRows[energy, constraint];
  activateSpatialDivergences[Join[{rows["MU_THETA"]},
    rows["U_INTERNAL"], {rows["EW_INTERNAL"]}]]
];

operatorFieldAtoms[expression_] := DeleteDuplicates[Join[
  Cases[expression,
    derivative : HoldPattern[Derivative[orders___][head_][arguments___]] /;
      MemberQ[{uOne, uTwo, uThree, thetaWave, eWave}, head] :>
        derivative, {0, Infinity}],
  Cases[expression,
    field : HoldPattern[head_[arguments___]] /;
      MemberQ[{uOne, uTwo, uThree, thetaWave, eWave}, head] :>
        field, {0, Infinity}]]];

operatorFreezeRankDiagnostic[branch_String, density_String] := Module[
  {records, coefficientRules, constraint, frozenRows, liveRows, allAtoms,
    atomTokens, atomRules, variables, frozenMatrix, liveMatrix,
    frozenRank, liveRank},
  records = baseEnergyDataByBranch[branch]["RECORDS"][[
    basisRepresentativeIndicesByBranch[branch]]];
  coefficientRules = Thread[Flatten[Values[newCoefficientSymbols]] -> 1];
  constraint = virtualConstraintSource["EULERIAN", branch, density,
    widthBase];
  frozenRows = Map[Function[record, elRowVector[
    applyProfile[record["ENERGY_TERM"] /. coefficientRules,
      "BASE", 0], applyProfile[constraint, "BASE", 0]]], records];
  liveRows = Map[Function[record,
    applyBackgroundProfileWithGeneratedJets[elRowVector[
      record["ENERGY_TERM"] /. coefficientRules, constraint]]], records];
  allAtoms = operatorFieldAtoms[Join[frozenRows, liveRows]];
  atomTokens = Table[Symbol["elFieldAtom" <> ToString[index]],
    {index, Length[allAtoms]}];
  atomRules = Thread[allAtoms -> atomTokens];
  variables = DeleteDuplicates[Join[atomTokens, rankVariables,
    generatedBackgroundJetAtoms]];
  frozenMatrix = polynomialFeatureMatrix[frozenRows /. atomRules,
    variables];
  liveMatrix = polynomialFeatureMatrix[liveRows /. atomRules, variables];
  frozenRank = MatrixRank[frozenMatrix];
  liveRank = MatrixRank[liveMatrix];
  <|"EL_ROW_COMPONENT_ORDER" -> {"MU_THETA", "U_ONE", "U_TWO",
      "U_THREE", "E_W"},
    "ENERGY_TERM_CONTRIBUTION_COUNT" -> Length[records],
    "PRODUCTION_FROZEN_EL_RANK" -> frozenRank,
    "LIVE_BACKGROUND_EL_RANK" -> liveRank,
    "LIVE_MINUS_FROZEN_RANK" -> liveRank - frozenRank,
    "PRODUCTION_FROZEN_ZERO_CONTRIBUTION_INDICES" ->
      Select[Range[Length[frozenMatrix]],
        MatrixRank[{frozenMatrix[[#]]}] == 0 &],
    "LIVE_BACKGROUND_ZERO_CONTRIBUTION_INDICES" ->
      Select[Range[Length[liveMatrix]],
        MatrixRank[{liveMatrix[[#]]}] == 0 &],
    "DIMENSION_L_T_M" -> dimensionZero,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>
];

modelDimensions = <|
  "U_MOMENTUM_ROWS" -> dimensionBulkForce,
  "THICKNESS_ROW" -> dimensionThicknessRow,
  "MASS_EVOLUTION_ROW" -> dimensionMassEvolution,
  "CENTER_FACE_GENERALIZED_ROW" -> dimensionThicknessRow|>;

divergenceSourceDimensions = <|"MU_THETA" -> dimensionMuTheta,
  "U_ENERGY_VARIATIONS" -> dimensionBulkForce,
  "EW_ENERGY_VARIATION" -> dimensionEnergy|>;
backgroundFirstJetDimensions = <|"W_BG" -> dimensionZero,
  "MU_R_BG" -> dimensionEnergy - dimensionL|>;
faceSubstrateDimensions = <|"NORMAL" -> dimensionZero,
  "MEASURE" -> dimensionZero, "NORMAL_VELOCITY" -> dimensionVelocity,
  "AFFINITY" -> dimensionAffinity, "RELATIVE_FLUX" -> dimensionFlux,
  "TRACTION_PRESSURE" -> dimensionPressure,
  "VIRTUAL_NORMAL_DISPLACEMENT" -> dimensionL,
  "VIRTUAL_WORK" -> dimensionEnergy|>;
faceSubstrateEpsilonOrders = <|"NORMAL" -> 0, "MEASURE" -> 0,
  "NORMAL_VELOCITY" -> 1, "AFFINITY" -> 1,
  "RELATIVE_FLUX" -> 1, "TRACTION_PRESSURE" -> 1,
  "VIRTUAL_NORMAL_DISPLACEMENT" -> 0, "VIRTUAL_WORK" -> 1|>;
faceSubstrateRecord[faces_Association] := Map[Function[face,
  AssociationMap[Function[field, preparedObjectRecord[face[field],
    faceSubstrateDimensions[field], faceSubstrateEpsilonOrders[field]]],
    Keys[face]]], faces];

caseKey[branch_String, density_String] :=
  branch <> "|" <> density;

modelRecord[model_Association] := <|
  "DIVERGENCE_FORM_SOURCE" -> preparedObjectRecord[
    model["DIVERGENCE_FORM_SOURCE"], divergenceSourceDimensions, 1],
  "BACKGROUND_FIRST_JETS" -> preparedObjectRecord[<|
    "W_BG" -> sigmaW widthJetSymbols,
    "MU_R_BG" -> (muR/WZero) sigmaW modulusJetSymbols|>,
    backgroundFirstJetDimensions, 0],
  "ROWS" -> AssociationMap[Function[row,
    objectRecord[model["OPERATOR"][row], modelDimensions[row], 1]],
    Keys[model["OPERATOR"]]],
  "VIRTUAL_CONSTRAINT" -> objectRecord[model["CONSTRAINT"],
    dimensionZero, 0],
  "FACE_SHAPE_SUBSTRATE" ->
    faceSubstrateRecord[model["FACE_SUBSTRATE"]]|>;

(* ---------------------------------------------------------------------- *)
(* Weak local differential-sector restrictions and adjoint operands.     *)
(* ---------------------------------------------------------------------- *)

transverseTrialVector = curl[transverseTrialPotentialVector];
transverseTestVector = curl[transverseTestPotentialVector];
longitudinalTrialVector = gradient[longitudinalTrialPotentialField];
longitudinalTestVector = gradient[longitudinalTestPotentialField];

linearRestrictedOperator[operator_Association, trialU_List, trialTheta_,
    trialEw_] := Module[{rules, scaled},
  rules = waveFieldHeadRules[sectorOrder trialU,
    sectorOrder trialTheta, sectorOrder trialEw];
  scaled = operator /. rules;
  Map[(D[#, sectorOrder] /. sectorOrder -> 0) &, scaled]
];

weakPairTerm[test_,
    coefficient_. Inactive[Div][flux_, coordinates_List]] /;
      SameQ[coordinates, spatialCoordinates] :=
  -Dot[gradient[coefficient test], flux];
weakPairTerm[test_, term_] := test term;

topLevelTerms[row_] := If[Head[row] === Plus, List @@ row, {row}];
divergenceTermQ[term_] := MatchQ[term,
  HoldPattern[coefficient_. Inactive[Div][flux_, coordinates_List]] /;
    SameQ[coordinates, spatialCoordinates]];

weakPairScalar[test_, row_] := Module[{terms},
  terms = topLevelTerms[row];
  Total[weakPairTerm[test, #] & /@ terms]
];
weakPairVector[test_List, rows_List] := Total[
  MapThread[weakPairScalar, {test, rows}]];

splitDivergenceRows[rows_List] := Module[{termRows, divergenceRows,
    directRows},
  termRows = topLevelTerms /@ rows;
  divergenceRows = Select[#, divergenceTermQ] & /@ termRows;
  directRows = Total[Select[#, !divergenceTermQ[#] &]] & /@ termRows;
  {divergenceRows, directRows}
];

weakLongitudinalPair[potential_, testVector_List, rows_List] := Module[
  {split, divergenceRows, directRows, divergencePairing},
  split = splitDivergenceRows[rows];
  divergenceRows = split[[1]];
  directRows = split[[2]];
  divergencePairing = Total[MapThread[
    Function[{test, terms}, Total[weakPairTerm[test, #] & /@ terms]],
    {testVector, divergenceRows}]];
  divergencePairing - potential divergence[directRows]
];

weakTransversePair[potentialVector_List, testVector_List,
    rows_List] := Module[
  {split, divergenceRows, directRows, divergencePairing},
  split = splitDivergenceRows[rows];
  divergenceRows = split[[1]];
  directRows = split[[2]];
  divergencePairing = Total[MapThread[
    Function[{test, terms}, Total[weakPairTerm[test, #] & /@ terms]],
    {testVector, divergenceRows}]];
  divergencePairing + Dot[potentialVector, curl[directRows]]
];

weakPairingRecord[density_] := <|
  "PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP" -> density,
  "IN_PLANE_BOUNDARY_TERM" -> 0,
  "SUPPORT" -> CompactInPlaneSupport|>;

extractCoupling[model_Association] := Module[
  {operator, transverseRows, thicknessRows, transverseToThickness,
    thicknessToTransverse, forwardDensity, reverseDensity,
    reverseRelabelRules,
    reverseRelabeledDensity, adjointnessResidual},
  operator = model["OPERATOR"];
  transverseRows = linearRestrictedOperator[operator,
    transverseTrialVector, 0, 0];
  thicknessRows = linearRestrictedOperator[operator,
    longitudinalTrialVector, thetaTrialField, ewTrialField];

  forwardDensity =
    weakPairScalar[thetaTestField,
      transverseRows["MASS_EVOLUTION_ROW"]] +
    weakPairScalar[ewTestField, transverseRows["THICKNESS_ROW"]] +
    weakLongitudinalPair[longitudinalTestPotentialField,
      longitudinalTestVector,
      transverseRows["U_MOMENTUM_ROWS"]];
  reverseDensity =
    weakTransversePair[transverseTestPotentialVector,
      transverseTestVector,
      thicknessRows["U_MOMENTUM_ROWS"]];
  reverseRelabelRules = {
    transverseTestPotentialOne -> functionAtCoordinates[
      transverseTrialPotentialVector[[1]]],
    transverseTestPotentialTwo -> functionAtCoordinates[
      transverseTrialPotentialVector[[2]]],
    transverseTestPotentialThree -> functionAtCoordinates[
      transverseTrialPotentialVector[[3]]],
    longitudinalTrialPotential -> functionAtCoordinates[
      longitudinalTestPotentialField],
    thetaTrial -> functionAtCoordinates[thetaTestField],
    ewTrial -> functionAtCoordinates[ewTestField]};
  reverseRelabeledDensity = reverseDensity /. reverseRelabelRules;
  adjointnessResidual = forwardDensity - reverseRelabeledDensity;

  transverseToThickness = <|
    "WEAK_PAIRING" -> weakPairingRecord[forwardDensity],
    "TEST_FIELDS" -> <|"THETA" -> thetaTestField,
      "E_W" -> ewTestField, "U_L" -> longitudinalTestVector|>|>;
  thicknessToTransverse = <|
    "WEAK_PAIRING" -> weakPairingRecord[reverseRelabeledDensity],
    "TEST_FIELD" -> transverseTestVector,
    "FORMAL_ADJOINT_RELABELING" -> reverseRelabelRules,
    "RELABELLED_WEAK_PAIRING_DENSITY" -> reverseRelabeledDensity|>;
  <|"TRANSVERSE_TO_THETA_EW_UL" -> transverseToThickness,
    "THETA_EW_UL_TO_TRANSVERSE" -> thicknessToTransverse,
    "ADJOINTNESS_OPERAND_FORWARD" -> weakPairingRecord[forwardDensity],
    "ADJOINTNESS_OPERAND_REVERSE" ->
      weakPairingRecord[reverseRelabeledDensity],
    "ADJOINTNESS_RELATION" -> relationalObject[forwardDensity,
      reverseRelabeledDensity],
    "ADJOINTNESS_RESIDUAL" -> weakPairingRecord[adjointnessResidual],
    "SECTOR_LABELS" -> <|
      "TRANSVERSE_TRIAL" -> transverseTrialVector,
      "TRANSVERSE_TEST" -> transverseTestVector,
      "LONGITUDINAL_TRIAL" -> longitudinalTrialVector,
      "LONGITUDINAL_TEST" -> longitudinalTestVector,
      "TRANSVERSE_CONSTRAINTS" -> {
        relationalObject[divergence[transverseTrialVector], 0],
        relationalObject[divergence[transverseTestVector], 0]},
      "LONGITUDINAL_CONSTRAINTS" -> {
        relationalObject[curl[longitudinalTrialVector], {0, 0, 0}],
        relationalObject[curl[longitudinalTestVector], {0, 0, 0}]},
      "PAIRING_DOMAIN" -> CompactInPlaneSupport|>|>
];

kernelDimensions = <|
  "TRANSVERSE_TO_THETA_EW_UL" -> <|
    "WEAK_PAIRING" -> dimensionEnergy|>,
  "THETA_EW_UL_TO_TRANSVERSE" -> <|
    "WEAK_PAIRING" -> dimensionEnergy,
    "RELABELLED_WEAK_PAIRING_DENSITY" -> dimensionEnergy|>,
  "ADJOINTNESS_OPERANDS" -> dimensionEnergy|>;

kernelRecord[kernel_Association] := <|
  "EXPRESSION" -> kernel,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
    (gradeTerms[kernel] /. {a_, b_} :> {1, a, b}),
  "DIMENSION_L_T_M" -> kernelDimensions|>;

removeZeroInactiveDivergences[expression_] := FixedPoint[
  ReplaceAll[#, divergenceObject : Inactive[Div][vector_List,
      coordinates_List] :> If[
        SameQ[coordinates, spatialCoordinates] &&
          AllTrue[vector, PossibleZeroQ], 0, divergenceObject]] &,
  expression];

simplifyWeakKernel[kernel_Association] := Module[
  {forward, reverse, residual, forwardBlock, reverseBlock},
  forward = FullSimplify[
    removeZeroInactiveDivergences[
      kernel["TRANSVERSE_TO_THETA_EW_UL"]["WEAK_PAIRING"][
        "PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP"]]];
  reverse = FullSimplify[
    removeZeroInactiveDivergences[
      kernel["THETA_EW_UL_TO_TRANSVERSE"]["WEAK_PAIRING"][
        "PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP"]]];
  residual = FullSimplify[forward - reverse];
  forwardBlock = Join[kernel["TRANSVERSE_TO_THETA_EW_UL"], <|
    "WEAK_PAIRING" -> weakPairingRecord[forward]|>];
  reverseBlock = Join[kernel["THETA_EW_UL_TO_TRANSVERSE"], <|
    "WEAK_PAIRING" -> weakPairingRecord[reverse],
    "RELABELLED_WEAK_PAIRING_DENSITY" -> reverse|>];
  Join[kernel, <|
    "TRANSVERSE_TO_THETA_EW_UL" -> forwardBlock,
    "THETA_EW_UL_TO_TRANSVERSE" -> reverseBlock,
    "ADJOINTNESS_OPERAND_FORWARD" -> weakPairingRecord[forward],
    "ADJOINTNESS_OPERAND_REVERSE" -> weakPairingRecord[reverse],
    "ADJOINTNESS_RELATION" -> relationalObject[forward, reverse],
    "ADJOINTNESS_RESIDUAL" -> weakPairingRecord[residual]|>]
];

kernelFromOrigin[originOperator_Association] := Module[{fakeModel},
  fakeModel = <|"OPERATOR" -> <|
    "U_MOMENTUM_ROWS" -> Lookup[originOperator,
      "U_MOMENTUM_ROWS", {0, 0, 0}],
    "THICKNESS_ROW" -> Lookup[originOperator, "THICKNESS_ROW", 0],
    "MASS_EVOLUTION_ROW" -> Lookup[originOperator,
      "MASS_EVOLUTION_ROW", 0]|>, "ENERGY" -> 0|>;
  KeyTake[extractCoupling[fakeModel],
    {"TRANSVERSE_TO_THETA_EW_UL", "THETA_EW_UL_TO_TRANSVERSE"}]
];

(* ---------------------------------------------------------------------- *)
(* Main variable-coefficient objects.                                    *)
(* ---------------------------------------------------------------------- *)

energyRecordsByBranch = Association[Table[branch -> applyProfile[
    If[branch === "LAB_HELD", baseEnergyData,
      constructEnergyData[branch, widthBase, modulusBase]]["RECORDS"][[
        basisRepresentativeIndices]]], {branch, branches}]];

energyBasisPayload = Association[Table[branch -> Module[
    {records, variableLabels},
    records = energyRecordsByBranch[branch];
    variableLabels = Select[records,
      !FreeQ[#["ENERGY_TERM"], etaBg | sigmaW] &][[All, "LABEL"]];
    <|"BACKGROUND_ANCHOR" -> If[branch === "LAB_HELD",
        widthBackground[xOne, xTwo, xThree],
        widthBackground[chiOne[xOne, xTwo, xThree, time],
          chiTwo[xOne, xTwo, xThree, time],
          chiThree[xOne, xTwo, xThree, time]]],
      "REPRESENTATIVE" -> Association[Map[
        #["LABEL"] -> <|"INVARIANT" -> #["INVARIANT"],
          "COEFFICIENT" -> #["COEFFICIENT"],
          "ENERGY_TERM" -> #["ENERGY_TERM"]|> &, records]],
      "UNIFORM_INVARIANTS_WITH_VARIABLE_COEFFICIENT" ->
        variableLabels,
      "EXACT_THICKNESS_MAP" -> relationalObject[eWBackground,
        (WZero/widthBackground[xOne, xTwo, xThree]) eWField],
      "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
        (gradeTerms[records[[All, "ENERGY_TERM"]]] /.
          {a_, b_} :> {2, a, b}),
      "DIMENSION_L_T_M" -> dimensionEnergy|>
  ], {branch, branches}]];
emitShared["ENERGY_BASIS_VARIABLE", energyBasisPayload];

energyBasisCountPayload = Association[Table[branch -> <|
    "COUNT_OPERAND" -> Length[energyRecordsByBranch[branch]],
    "RANK_SELECTED_INDICES" -> basisRepresentativeIndices,
    "DIMENSION_L_T_M" -> dimensionZero,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>,
  {branch, branches}]];
emitShared["ENERGY_BASIS_COUNT", energyBasisCountPayload];

(* Dedicated form ablation to the committed eight-form slice. *)
restrictedEnergyDataByBranch = Association[Table[branch ->
  constructEnergyData[branch, widthBase, modulusBase, True, True],
  {branch, branches}]];
restrictedIndicesByBranch = AssociationMap[
  independentRepresentativeIndices[
    restrictedEnergyDataByBranch[#]["RECORDS"][[All, "INVARIANT"]]] &,
  branches];
restrictedCountPayload = AssociationMap[Function[branch, Module[
  {count, labels},
  count = Length[restrictedIndicesByBranch[branch]];
  labels = restrictedEnergyDataByBranch[branch]["RECORDS"][[
    restrictedIndicesByBranch[branch], "LABEL"]];
  <|"RESTRICT_TO_COMMITTED_FORMS_SWITCH" -> True,
    "COUNT_OPERAND" -> count,
    "REFERENCE_OPERAND" -> 26,
    "RESIDUAL_COUNT_MINUS_REFERENCE" -> count - 26,
    "RANK_SELECTED_INDICES" -> restrictedIndicesByBranch[branch],
    "NEW_INVARIANT_LABELS" -> Select[labels,
      !MemberQ[uniformLabels, #] &],
    "DIMENSION_L_T_M" -> dimensionZero,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>]], branches];
emitShared["ENERGY_BASIS_RESTRICTED_FORM_OPERAND",
  restrictedCountPayload];

formAblationCountPayload = AssociationMap[Function[branch, <|
  "COMPLETED_COUNT_OPERAND" ->
    Length[basisRepresentativeIndicesByBranch[branch]],
  "RESTRICTED_COUNT_OPERAND" -> Length[restrictedIndicesByBranch[branch]],
  "COMPLETED_MINUS_RESTRICTED" ->
    Length[basisRepresentativeIndicesByBranch[branch]] -
      Length[restrictedIndicesByBranch[branch]],
  "COMPLETED_SELECTED_INDICES" ->
    basisRepresentativeIndicesByBranch[branch],
  "RESTRICTED_SELECTED_INDICES" -> restrictedIndicesByBranch[branch],
  "DIMENSION_L_T_M" -> dimensionZero,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>], branches];
emitShared["ENERGY_BASIS_FORM_ABLATION_COUNT", formAblationCountPayload];

o3KroneckerEnumerationPayload = <|
  "SPURION_TENSOR_RANK" -> 1,
  "DOF_BUILDING_BLOCKS" -> dofFactorSpecifications,
  "CONTRACTION_TENSOR" -> HoldForm[KroneckerDelta],
  "GENERATED_BLUEPRINT_COUNT" -> Length[o3KroneckerBlueprints],
  "GENERATED_BLUEPRINTS" -> Map[KeyTake[#,
    {"FACTOR_NAMES", "FACTOR_RANKS", "PAIRING", "LABEL_SUFFIX"}] &,
    o3KroneckerBlueprints],
  "PAIRING_SLOT_COVERAGE_RESIDUAL" -> Map[Function[blueprint, Module[
    {ranks, expectedSlots, pairedSlots},
    ranks = Join[{1}, blueprint["FACTOR_RANKS"]];
    expectedSlots = Flatten[Table[Table[{block, slot},
      {slot, ranks[[block]]}], {block, Length[ranks]}], 1];
    pairedSlots = Flatten[blueprint["PAIRING"], 1];
    <|"MISSING_SLOTS" -> Complement[expectedSlots, pairedSlots],
      "DUPLICATED_SLOT_COUNT" -> Length[pairedSlots] -
        Length[DeleteDuplicates[pairedSlots]]|>
  ]], o3KroneckerBlueprints],
  "DIMENSION_L_T_M" -> dimensionZero,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>;
emitShared["ENERGY_BASIS_O3_KRONECKER_ENUMERATION",
  o3KroneckerEnumerationPayload];

recordsForProfileSource[data_Association, source_String] := Module[
  {prefix},
  prefix = If[source === "W_BG", "WJET_", "MUJET_"];
  Select[data["RECORDS"], StringStartsQ[#1["LABEL"], prefix] &]
];

completenessPayload = Association[Flatten[Table[
  (branch <> "|" <> source) -> Module[
    {fullRecords, legacyRecords, fullExpressions, legacyExpressions,
      fullRank, combinedRank, epsilonRecord, epsilonRank,
      selectedIndices, selectedMatrix, zeroRows, dependentPairs},
    fullRecords = recordsForProfileSource[
      baseEnergyDataByBranch[branch], source];
    legacyRecords = recordsForProfileSource[
      restrictedEnergyDataByBranch[branch], source];
    fullExpressions = fullRecords[[All, "INVARIANT"]];
    legacyExpressions = legacyRecords[[All, "INVARIANT"]];
    fullRank = MatrixRank[rankMatrix[fullExpressions]];
    combinedRank = MatrixRank[rankMatrix[
      Join[fullExpressions, legacyExpressions]]];
    epsilonRecord = First[epsilonPseudoscalarCandidates[branch, source]];
    epsilonRank = MatrixRank[rankMatrix[Append[fullExpressions,
      epsilonRecord["EXPRESSION"]]]];
    selectedIndices = independentRepresentativeIndices[fullExpressions];
    selectedMatrix = rankMatrix[fullExpressions][[selectedIndices]];
    zeroRows = Select[Range[Length[selectedMatrix]],
      MatrixRank[{selectedMatrix[[#]]}] == 0 &];
    dependentPairs = Select[Subsets[Range[Length[selectedMatrix]], {2}],
      MatrixRank[selectedMatrix[[#]]] < 2 &];
    <|"GENERATED_FAMILY_RANK" -> fullRank,
      "COMMITTED_FORM_RANK" -> MatrixRank[rankMatrix[legacyExpressions]],
      "GENERATED_PLUS_COMMITTED_RANK" -> combinedRank,
      "COMMITTED_RANK_INCREMENT_ON_GENERATED" -> combinedRank - fullRank,
      "EPSILON_TRIAL_FACTOR_CONTENT" ->
        epsilonRecord["FACTOR_NAMES"],
      "EPSILON_TRIAL_EXPRESSION" -> applyProfile[
        epsilonRecord["EXPRESSION"], "BASE", 0],
      "EPSILON_TRIAL_RANK_OPERAND" -> epsilonRank,
      "EPSILON_TRIAL_RANK_INCREMENT" -> epsilonRank - fullRank,
      "GENERATED_CONTRACTION_TENSOR" -> HoldForm[KroneckerDelta],
      "EPSILON_TRIAL_CONTRACTION_TENSOR" -> HoldForm[Signature],
      "SELECTED_COUNT_OPERAND" -> Length[selectedIndices],
      "SELECTED_RANK_OPERAND" -> MatrixRank[selectedMatrix],
      "ZERO_ROW_INDICES" -> zeroRows,
      "DEPENDENT_DIRECTION_PAIRS" -> dependentPairs,
      "DIMENSION_L_T_M" -> dimensionZero,
      "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>
  ], {branch, branches}, {source, profileSources}], 1]];
emitShared["ENERGY_BASIS_KRONECKER_COMPLETENESS_OPERAND",
  completenessPayload];

coefficientRescalingPayload = AssociationMap[Function[branch, Module[
  {expressions, labels, firstNewPosition, rescaledExpressions,
    baseSelected, rescaledSelected},
  expressions = baseEnergyDataByBranch[branch]["RECORDS"][[All,
    "INVARIANT"]];
  labels = baseEnergyDataByBranch[branch]["RECORDS"][[All, "LABEL"]];
  firstNewPosition = First[FirstPosition[labels,
    label_String /; !MemberQ[uniformLabels, label]]];
  rescaledExpressions = ReplacePart[expressions,
    firstNewPosition -> 7 expressions[[firstNewPosition]]];
  baseSelected = independentRepresentativeIndices[expressions];
  rescaledSelected = independentRepresentativeIndices[rescaledExpressions];
  <|"RESCALED_LABEL" -> labels[[firstNewPosition]],
    "NONZERO_RESCALING_OPERAND" -> 7,
    "BASE_COUNT_OPERAND" -> Length[baseSelected],
    "RESCALED_COUNT_OPERAND" -> Length[rescaledSelected],
    "COUNT_RESIDUAL" -> Length[rescaledSelected] - Length[baseSelected],
    "BASE_SELECTED_INDICES" -> baseSelected,
    "RESCALED_SELECTED_INDICES" -> rescaledSelected,
    "SELECTED_INDEX_RESIDUAL" ->
      Complement[Join[baseSelected, rescaledSelected],
        Intersection[baseSelected, rescaledSelected]],
    "DIMENSION_L_T_M" -> dimensionZero,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>]], branches];
emitShared["ENERGY_BASIS_COEFFICIENT_RESCALING_CONTROL",
  coefficientRescalingPayload];

rawThicknessEnergyDataByBranch = Association[Table[branch ->
  constructEnergyData[branch, widthBase, modulusBase, False, False],
  {branch, branches}]];
thicknessMapResidualPayload = Association[Flatten[Table[
  (branch <> "|" <> source) -> Module[
    {mappedRecords, rawRecords},
    mappedRecords = recordsForProfileSource[
      baseEnergyDataByBranch[branch], source];
    rawRecords = recordsForProfileSource[
      rawThicknessEnergyDataByBranch[branch], source];
    Association[MapThread[Function[{mappedRecord, rawRecord},
      mappedRecord["LABEL"] -> <|
        "MAPPED_LOCAL_FIELD_OPERAND" -> applyProfile[
          mappedRecord["INVARIANT"], "BASE", 0],
        "RAW_FIELD_OPERAND" -> applyProfile[
          rawRecord["INVARIANT"], "BASE", 0],
        "RESIDUAL_MAPPED_MINUS_RAW" -> applyProfile[
          mappedRecord["INVARIANT"] - rawRecord["INVARIANT"],
          "BASE", 0],
        "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
          withEpsilonGrade[applyProfile[mappedRecord["INVARIANT"] -
            rawRecord["INVARIANT"], "BASE", 0], 2],
        "DIMENSION_L_T_M" ->
          dimensionGradient + Total[First[Select[
            o3KroneckerBlueprints,
            StringEndsQ[mappedRecord["LABEL"],
              #1["LABEL_SUFFIX"]] &]]["FACTOR_DIMENSIONS"][[2 ;;]]]|>
    ], {mappedRecords, rawRecords}]]
  ], {branch, branches}, {source, profileSources}], 1]];
emitShared["ENERGY_BASIS_THICKNESS_MAP_RESIDUAL",
  thicknessMapResidualPayload];

generatedBackgroundJetAtoms = DeleteDuplicates[Flatten[Table[
  If[1 <= i + j + k <= 3,
    Table[profileJetSymbol[source, {i, j, k}],
      {source, profileSources}], Nothing],
  {i, 0, 3}, {j, 0, 3}, {k, 0, 3}]]];
backgroundSecondJetAtoms = DeleteDuplicates[Flatten[Table[
  If[i + j + k == 2,
    Table[profileJetSymbol[source, {i, j, k}],
      {source, profileSources}], Nothing],
  {i, 0, 2}, {j, 0, 2}, {k, 0, 2}]]];
energyHessianGuardPayload = AssociationMap[Function[branch, Module[
  {expressions, liveExpressions, zeroHessianExpressions, variables,
    liveRank, zeroHessianRank},
  expressions = baseEnergyDataByBranch[branch]["RECORDS"][[
    basisRepresentativeIndicesByBranch[branch], "INVARIANT"]];
  liveExpressions = (Expand[
    applyBackgroundProfileWithGeneratedJets[#] /. fieldAlgebraRules] & /@
    expressions);
  zeroHessianExpressions = liveExpressions /.
    Thread[backgroundSecondJetAtoms -> 0];
  variables = DeleteDuplicates[Join[rankVariables,
    generatedBackgroundJetAtoms]];
  liveRank = MatrixRank[polynomialFeatureMatrix[
    liveExpressions, variables]];
  zeroHessianRank = MatrixRank[polynomialFeatureMatrix[
    zeroHessianExpressions, variables]];
  <|"WITH_BACKGROUND_SECOND_JET_ATOMS_RANK" -> liveRank,
    "SECOND_JET_ATOMS_ZERO_RANK" -> zeroHessianRank,
    "RANK_DIFFERENCE" -> liveRank - zeroHessianRank,
    "SECOND_JET_ATOMS" -> backgroundSecondJetAtoms,
    "DIMENSION_L_T_M" -> dimensionZero,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>]], branches];
emitShared["ENERGY_BASIS_HESSIAN_IN_ENERGY_GUARD",
  energyHessianGuardPayload];

operatorFreezeDiagnosticPayload = Association[Flatten[Table[
  caseKey[branch, density] ->
    operatorFreezeRankDiagnostic[branch, density],
  {branch, branches}, {density, densities}], 1]];
emitShared["OPERATOR_BACKGROUND_FREEZE_DIAGNOSTIC",
  operatorFreezeDiagnosticPayload];

newInvariantPayload = Association[Table[branch -> Module[
    {records},
    records = Select[
      energyRecordsByBranch[branch],
      !MemberQ[uniformLabels, #["LABEL"]] &];
    Association[Map[#["LABEL"] -> <|
      "INVARIANT" -> #["INVARIANT"],
      "COEFFICIENT" -> #["COEFFICIENT"],
      "COEFFICIENT_STANDARD_NAME" -> #["COEFFICIENT_STANDARD_NAME"],
      "ENERGY_TERM" -> #["ENERGY_TERM"],
      "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
        (gradeTerms[#["ENERGY_TERM"]] /. {a_, b_} :> {2, a, b}),
      "DIMENSION_L_T_M" -> dimensionEnergy|> &, records]]
  ], {branch, branches}]];
emitShared["ENERGY_BASIS_NEW_INVARIANTS", newInvariantPayload];

energyBasisOmissionsPayload = Association[Table[branch -> Module[
    {records, selectedUniform, omitted},
    records = energyRecordsByBranch[branch];
    selectedUniform = Select[records,
      MemberQ[uniformLabels, #["LABEL"]] &][[All, "LABEL"]];
    omitted = Complement[selectedUniform, carriedUniformLabels];
    <|"CARRIED_INPUT_LABELS" -> carriedUniformLabels,
      "CONSTRUCTED_UNIFORM_LABELS" -> selectedUniform,
      "OMITTED_FROM_CARRIED_INPUT" -> omitted,
      "OMISSION_COUNT_OPERAND" -> Length[omitted],
      "DIMENSION_L_T_M" -> dimensionZero,
      "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>
  ], {branch, branches}]];
emitShared["ENERGY_BASIS_OMISSIONS", energyBasisOmissionsPayload];

Clear[energyBasisPayload, energyBasisCountPayload, newInvariantPayload,
  energyBasisOmissionsPayload, energyRecordsByBranch, baseEnergyData,
  baseEnergyDataByBranch, restrictedEnergyDataByBranch,
  restrictedIndicesByBranch, restrictedCountPayload,
  formAblationCountPayload, completenessPayload,
  coefficientRescalingPayload, rawThicknessEnergyDataByBranch,
  thicknessMapResidualPayload, energyHessianGuardPayload,
  operatorFreezeDiagnosticPayload, generatedBackgroundJetAtoms,
  backgroundSecondJetAtoms];
ClearSystemCache[];

mainModels = <||>;
Do[Module[{processed},
  processed = evaluatedModel["EULERIAN", branch, density];
  AssociateTo[mainModels, caseKey[branch, density] -> processed];
  Clear[processed];
  ClearSystemCache[]],
  {branch, branches}, {density, densities}];

slabOperatorPayload = Map[modelRecord, mainModels];
emitShared["SLAB_OPERATOR", slabOperatorPayload];
Clear[slabOperatorPayload];

operatorOriginsPayload = Map[Function[model,
  preparedObjectRecord[model["ORIGINS"], modelDimensions, 1]],
  mainModels];
emitShared["SLAB_OPERATOR_TERM_ORIGINS", operatorOriginsPayload];
Clear[operatorOriginsPayload];

muThetaPayload = Association[KeyValueMap[Function[{key, model},
  key -> <|
    "VARIABLE_COEFFICIENT_VARIATIONAL_OPERAND" ->
      objectRecord[model["MU_THETA"], dimensionMuTheta, 1],
    "HELD_FIXED_FIELDS" -> {uHeldFixed, eWHeldFixed,
      allOtherFieldsHeldFixed},
    "BACKGROUND_ANCHOR" -> If[
      StringStartsQ[key, "MATERIAL_ADVECTED"],
      branchAnchorMappedFromChi, branchAnchorAtEulerianPosition],
    "RESERVED_OPERAND" -> muThetaReserved|>], mainModels]];
emitShared["MU_THETA_OPERATOR", muThetaPayload];
Clear[muThetaPayload];

mainKernels = Map[extractCoupling, mainModels];
couplingKernelPayload = Map[kernelRecord, mainKernels];
emitShared["COUPLING_KERNEL", couplingKernelPayload];
Clear[couplingKernelPayload];

kernelOriginsPayload = Map[Function[model, preparedObjectRecord[
  Map[Function[origin,
    kernelFromOrigin[If[AssociationQ[origin], origin, <||>]]],
    model["ORIGINS"]], kernelDimensions, 1]], mainModels];
emitShared["COUPLING_KERNEL_TERM_ORIGINS", kernelOriginsPayload];
Clear[kernelOriginsPayload];

(* ---------------------------------------------------------------------- *)
(* Background-order admissibility pairing.                               *)
(* ---------------------------------------------------------------------- *)

backgroundBalanceFromModel[model_Association] := Module[
  {branch, fullFieldEnergy, firstVariation, bodyRows, ewBodyForce,
    faceTractions},
  branch = model["CONSTRUCTION_METADATA"]["BRANCH"];
  fullFieldEnergy = constructFullFieldBackgroundEnergy[branch,
    widthBase, modulusBase];
  firstVariation = D[fullFieldEnergy, backgroundOrder] /.
    backgroundOrder -> 0;
  bodyRows = <|
    "U" -> applyBackgroundProfileWithGeneratedJets[
      eulerVector[firstVariation, uVector]],
    "THETA" -> applyBackgroundProfileWithGeneratedJets[
      eulerScalar[firstVariation, thetaField]],
    "E_W" -> applyBackgroundProfileWithGeneratedJets[
      eulerScalar[firstVariation, eWField]]|>;
  ewBodyForce = bodyRows["E_W"];
  faceTractions = AssociationMap[Function[face,
    With[{normal = model["FACE_SUBSTRATE"][face]["NORMAL"]},
      truncateBackground[(ewBodyForce/WZero) normal]]],
    {"UPPER", "LOWER"}];
  <|"BULK_DOF_BODY_FORCE" -> bodyRows,
    "PER_FACE_TRACTION" -> faceTractions|>
];

admissibilityOperatorPayload = Map[Function[model, <|
  "BACKGROUND_ORDER_BALANCE" -> backgroundBalanceFromModel[model],
  "ORDER_SOURCE" -> fullFieldBackgroundFirstVariationAtBZero,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
    withEpsilonGrade[backgroundBalanceFromModel[model], 0],
  "DIMENSION_L_T_M" -> <|
    "BULK_DOF_BODY_FORCE" -> <|"U" -> dimensionBulkForce,
      "THETA" -> dimensionEnergy, "E_W" -> dimensionEnergy|>,
    "PER_FACE_TRACTION" -> dimensionFaceTraction|>|>], mainModels];
emitShared["ADMISSIBILITY_OPERATOR_OPERAND",
  admissibilityOperatorPayload];

supportForCase[branch_String, density_String] := <|
  "BULK_DOF_BODY_FORCE" -> <|
    "U" -> Through[{forceHoldOne, forceHoldTwo, forceHoldThree}[
      xOne, xTwo, xThree]],
    "THETA" -> forceHoldTheta[xOne, xTwo, xThree],
    "E_W" -> forceHoldEw[xOne, xTwo, xThree]|>,
  "PER_FACE_TRACTION" -> <|
    "UPPER" -> Through[{tractionHoldUpperOne, tractionHoldUpperTwo,
      tractionHoldUpperThree, tractionHoldUpperW}[
        xOne, xTwo, xThree]],
    "LOWER" -> Through[{tractionHoldLowerOne, tractionHoldLowerTwo,
      tractionHoldLowerThree, tractionHoldLowerW}[
        xOne, xTwo, xThree]]|>|>;

admissibilitySupportPayload = Association[Flatten[Table[
  caseKey[branch, density] -> <|
    "DECLARED_SUPPORT_BUNDLE" -> supportForCase[branch, density],
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
    "DIMENSION_L_T_M" -> <|
      "BULK_DOF_BODY_FORCE" -> <|"U" -> dimensionBulkForce,
        "THETA" -> dimensionEnergy, "E_W" -> dimensionEnergy|>,
      "PER_FACE_TRACTION" -> dimensionFaceTraction|>|>,
  {branch, branches}, {density, densities}], 1]];
emitShared["ADMISSIBILITY_SUPPORT_OPERAND", admissibilitySupportPayload];

admissibilityResidualPayload = AssociationMap[Function[key, Module[
  {operator, support, residual},
  operator = admissibilityOperatorPayload[key]["BACKGROUND_ORDER_BALANCE"];
  support = admissibilitySupportPayload[key]["DECLARED_SUPPORT_BUNDLE"];
  residual = mapDifference[operator, support];
  <|"RESIDUAL_OPERATOR_MINUS_SUPPORT" -> residual,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" ->
      withEpsilonGrade[residual, 0],
    "DIMENSION_L_T_M" -> admissibilitySupportPayload[key][
      "DIMENSION_L_T_M"]|>
]], Keys[admissibilityOperatorPayload]];
emitShared["ADMISSIBILITY_RESIDUAL", admissibilityResidualPayload];
Clear[admissibilityOperatorPayload, admissibilitySupportPayload,
  admissibilityResidualPayload];
ClearSystemCache[];

admissibilityPairingDimensions = <|
  "BULK_DOF_BODY_FORCE" -> <|"U" -> dimensionBulkForce,
    "THETA" -> dimensionEnergy, "E_W" -> dimensionEnergy|>,
  "PER_FACE_TRACTION" -> dimensionFaceTraction|>;
comparisonPackageDimensions = <|
  "SLAB_OPERATOR" -> modelDimensions,
  "COUPLING_KERNEL" -> kernelDimensions,
  "ADMISSIBILITY_OPERATOR" -> admissibilityPairingDimensions|>;
comparisonPackageEpsilonOrders = <|
  "SLAB_OPERATOR" -> 1, "COUPLING_KERNEL" -> 1,
  "ADMISSIBILITY_OPERATOR" -> 0|>;
comparisonPackageRecord[package_Association] := AssociationMap[
  Function[object, preparedObjectRecord[package[object],
    comparisonPackageDimensions[object],
    comparisonPackageEpsilonOrders[object]]], Keys[package]];
comparisonPackageResidualRecord[left_Association,
    right_Association] := AssociationMap[Function[object,
  preparedObjectRecord[mapDifference[left[object], right[object]],
    comparisonPackageDimensions[object],
    comparisonPackageEpsilonOrders[object]]], Keys[left]];

(* ---------------------------------------------------------------------- *)
(* Representation routes and one-sided source corruption.                *)
(* ---------------------------------------------------------------------- *)

representationEulerian = <||>;
representationMaterial = <||>;
Do[Module[{key, eulerian, material},
  key = caseKey[branch, density];
  eulerian = mainModels[key];
  material = evaluatedModel["MATERIAL", branch, density];
  AssociateTo[representationEulerian, key -> <|
    "SLAB_OPERATOR" -> eulerian["OPERATOR"],
    "COUPLING_KERNEL" -> mainKernels[key],
    "ADMISSIBILITY_OPERATOR" -> backgroundBalanceFromModel[eulerian]|>];
  AssociateTo[representationMaterial, key -> <|
    "SLAB_OPERATOR" -> material["OPERATOR"],
    "COUPLING_KERNEL" -> extractCoupling[material],
    "ADMISSIBILITY_OPERATOR" -> backgroundBalanceFromModel[material]|>]],
  {branch, branches}, {density, densities}];

emitShared["REP_INVARIANCE_EULERIAN_OPERAND",
  Map[comparisonPackageRecord, representationEulerian]];
emitShared["REP_INVARIANCE_MATERIAL_OPERAND",
  Map[comparisonPackageRecord, representationMaterial]];

representationResidual = AssociationMap[Function[key, <|
  "RESIDUAL_A_MINUS_B" -> comparisonPackageResidualRecord[
    representationEulerian[key], representationMaterial[key]]|>],
  Keys[representationEulerian]];
emitShared["REP_INVARIANCE_RESIDUAL", representationResidual];

independenceBase = <||>;
independenceCorrupted = <||>;
Do[Module[{key, base, corrupted},
  key = caseKey[branch, density];
  base = mainModels[key];
  corrupted = evaluatedModel["EULERIAN", branch, density,
    "BASE", 0, True];
  AssociateTo[independenceBase, key -> <|
    "SLAB_OPERATOR" -> base["OPERATOR"],
    "COUPLING_KERNEL" -> mainKernels[key],
    "ADMISSIBILITY_OPERATOR" -> backgroundBalanceFromModel[base]|>];
  AssociateTo[independenceCorrupted, key -> <|
    "SLAB_OPERATOR" -> corrupted["OPERATOR"],
    "COUPLING_KERNEL" -> extractCoupling[corrupted],
    "ADMISSIBILITY_OPERATOR" -> backgroundBalanceFromModel[corrupted]|>]],
  {branch, branches}, {density, densities}];
emitShared["CONTROL_INDEPENDENCE_BASE_OPERAND",
  Map[comparisonPackageRecord, independenceBase]];
emitShared["CONTROL_INDEPENDENCE_CORRUPTED_OPERAND",
  Map[comparisonPackageRecord, independenceCorrupted]];
emitShared["CONTROL_INDEPENDENCE_RESIDUAL", AssociationMap[
  Function[key, <|"RESIDUAL_A_MINUS_B" ->
    comparisonPackageResidualRecord[independenceBase[key],
      independenceCorrupted[key]]|>],
  Keys[independenceBase]]];

Clear[representationEulerian, representationMaterial,
  representationResidual, independenceBase, independenceCorrupted];
ClearSystemCache[];

(* ---------------------------------------------------------------------- *)
(* Per-source, per-direction form ablations.                              *)
(* ---------------------------------------------------------------------- *)

formObjectValue["ENERGY_BASIS_VARIABLE", model_Association,
    kernel_] := Association[Map[
  #["LABEL"] -> #["ENERGY_TERM"] &, model["SELECTED_RECORDS"]]];
formObjectValue["SLAB_OPERATOR", model_Association, kernel_] :=
  model["OPERATOR"];
formObjectValue["COUPLING_KERNEL", model_Association, kernel_] := kernel;

formObjectNames = {"ENERGY_BASIS_VARIABLE", "SLAB_OPERATOR",
  "COUPLING_KERNEL"};
formObjectDimension["ENERGY_BASIS_VARIABLE"] := dimensionEnergy;
formObjectDimension["SLAB_OPERATOR"] := modelDimensions;
formObjectDimension["COUPLING_KERNEL"] := kernelDimensions;
formObjectEpsilonOrder["ENERGY_BASIS_VARIABLE"] := 2;
formObjectEpsilonOrder["SLAB_OPERATOR"] := 1;
formObjectEpsilonOrder["COUPLING_KERNEL"] := 1;
formObjectRecord[object_String, value_] := preparedObjectRecord[value,
  formObjectDimension[object], formObjectEpsilonOrder[object]];

beginAssociationEmission[sharedObject["CONTROL_FORM_BASE_OPERAND"]];
firstFormEmission = True;
Do[Module[{baseline, baselineKernel, value, key},
  baseline = mainModels[caseKey[branch, density]];
  baselineKernel = mainKernels[caseKey[branch, density]];
  Do[
    key = object <> "|" <> branch <> "|" <> density <> "|" <>
      source <> "|DIRECTION_" <> ToString[direction];
    value = formObjectValue[object, baseline, baselineKernel];
    appendAssociationEmission[key, formObjectRecord[object, value],
      firstFormEmission];
    firstFormEmission = False,
    {source, profileSources}, {direction, directions},
    {object, formObjectNames}]],
  {branch, branches}, {density, densities}];
endAssociationEmission[];

Clear[mainModels, mainKernels];
ClearSystemCache[];

formParameters = {widthFormOne, widthFormTwo, widthFormThree,
  modulusFormOne, modulusFormTwo, modulusFormThree};
formBaseRules = Thread[formParameters -> 1];
formAblationRules[source_String, direction_Integer] := Module[
  {values, position},
  values = ConstantArray[1, Length[formParameters]];
  position = If[source === "W_BG", direction,
    braneDimension + direction];
  values[[position]] = 0;
  Thread[formParameters -> values]
];
parametricSourceName["W_BG"] := "PARAMETRIC_W";
parametricSourceName["MU_R_BG"] := "PARAMETRIC_MU";

beginAssociationEmission[sharedObject["CONTROL_FORM_ABLATED_OPERAND"]];
firstFormEmission = True;
Do[Module[{parametric, parametricKernel, ablated, ablatedKernel,
    value, key, rules},
  parametric = evaluatedModel["EULERIAN", branch, density,
    parametricSourceName[source], 0, False, True];
  parametricKernel = extractCoupling[parametric];
  Do[
    rules = formAblationRules[source, direction];
    ablated = parametric /. rules;
    ablatedKernel = parametricKernel /. rules;
    Do[
      key = object <> "|" <> branch <> "|" <> density <> "|" <>
      source <> "|DIRECTION_" <> ToString[direction];
      value = formObjectValue[object, ablated, ablatedKernel];
      appendAssociationEmission[key, formObjectRecord[object, value],
        firstFormEmission];
      firstFormEmission = False,
      {object, formObjectNames}];
    Clear[ablated, ablatedKernel, value, rules];
    ClearSystemCache[];
    Share[],
    {direction, directions}];
  Clear[parametric, parametricKernel];
  ClearSystemCache[];
  Share[]],
  {branch, branches}, {density, densities}, {source, profileSources}];
endAssociationEmission[];

beginAssociationEmission[sharedObject["CONTROL_FORM_RESIDUAL"]];
firstFormEmission = True;
Do[Module[{parametric, parametricKernel, baseline, baselineKernel,
    ablated, ablatedKernel, baseValue, ablatedValue, residualValue,
    key, rules},
  parametric = evaluatedModel["EULERIAN", branch, density,
    parametricSourceName[source], 0, False, True];
  parametricKernel = extractCoupling[parametric];
  baseline = parametric /. formBaseRules;
  baselineKernel = parametricKernel /. formBaseRules;
  Do[
    rules = formAblationRules[source, direction];
    ablated = parametric /. rules;
    ablatedKernel = parametricKernel /. rules;
    Do[
      key = object <> "|" <> branch <> "|" <> density <> "|" <>
        source <> "|DIRECTION_" <> ToString[direction];
      baseValue = formObjectValue[object, baseline, baselineKernel];
      ablatedValue = formObjectValue[object, ablated, ablatedKernel];
      residualValue = mapDifference[baseValue, ablatedValue];
      appendAssociationEmission[key,
        formObjectRecord[object, residualValue], firstFormEmission];
      firstFormEmission = False,
      {object, formObjectNames}];
    Clear[ablated, ablatedKernel, ablatedValue, residualValue, rules];
    ClearSystemCache[];
    Share[],
    {direction, directions}];
  Clear[parametric, parametricKernel, baseline, baselineKernel,
    baseValue];
  ClearSystemCache[];
  Share[]],
  {branch, branches}, {density, densities}, {source, profileSources}];
endAssociationEmission[];

Clear[firstFormEmission, formObjectNames, formParameters, formBaseRules];
ClearSystemCache[];
Share[];

(* ---------------------------------------------------------------------- *)
(* Independently reconstructed uniform S11b operands.                    *)
(* ---------------------------------------------------------------------- *)

constructUniformS11bEnergy[] := Module[
  {curlU, divU, gradTheta, gradEw, invariants, coefficients, factors},
  curlU = curl[uVector];
  divU = divergence[uVector];
  gradTheta = gradient[thetaField];
  gradEw = gradient[eWField];
  invariants = {Dot[curlU, curlU], divU^2, thetaField^2,
    eWField^2, thetaField eWField, Dot[gradTheta, gradTheta],
    Dot[gradEw, gradEw], Dot[gradTheta, gradEw],
    thetaField divU, eWField divU};
  coefficients = {muR, modulusDivU, bRho WZero, kW WZero^2,
    cCross WZero, gradientThetaCoefficient, kappaW WZero^4,
    gradientThetaEwCoefficient, thetaDivUCoefficient,
    ewDivUCoefficient};
  factors = {1/2, 1/2, 1/2, 1/2, 1, 1/2, 1/2, 1, 1, 1};
  Total[MapThread[#1 #2 #3 &, {factors, coefficients, invariants}]]
];

uniformS11bModel[] := Module[
  {energy, constraint, rows, muTheta, faceData, faceRows, operator,
    massRow},
  energy = constructUniformS11bEnergy[];
  constraint = virtualThetaField + virtualEwField +
    divergence[virtualUVector];
  rows = constrainedRows[energy, constraint];
  muTheta = rows["MU_THETA"];
  faceData = Association[Table[faceLabel[sign] -> Module[
    {velocity, affinity, flux, pressure, virtualDisplacement},
    velocity = WZero D[eWField, time]/2 +
      sign D[zetaCenterField, time];
    affinity = muTheta/rhoBr - pressureField[sign]/rhoM;
    flux = lambdaAResponse affinity + lambdaVResponse velocity;
    pressure = pressureField[sign] + lambdaXResponse affinity;
    virtualDisplacement = WZero virtualEwField/2 +
      sign virtualZetaCenterField;
    <|"NORMAL_VELOCITY" -> velocity, "AFFINITY" -> affinity,
      "RELATIVE_FLUX" -> flux, "TRACTION_PRESSURE" -> pressure,
      "VIRTUAL_WORK" -> -pressure virtualDisplacement|>],
    {sign, faces}]];
  faceRows = faceGeneralizedRows[faceData];
  massRow = rhoBr D[thetaField + eWField + divergence[uVector], time] +
    Total[faceData[[All, "RELATIVE_FLUX"]]];
  operator = <|
    "U_MOMENTUM_ROWS" -> rhoBr D[uVector, {time, 2}] +
      rows["U_INTERNAL"] + faceRows["U_FACE"],
    "THICKNESS_ROW" -> muW WZero^2 D[eWField, {time, 2}] +
      rows["EW_INTERNAL"] + faceRows["EW_FACE"],
    "MASS_EVOLUTION_ROW" -> massRow,
    "CENTER_FACE_GENERALIZED_ROW" -> faceRows["CENTER_FACE"]|>;
  <|"ENERGY" -> energy, "MU_THETA" -> muTheta,
    "OPERATOR" -> operator|>
];

uniformLimit[expression_] := truncateBackground[
  expression /. {etaBg -> 0, sigmaW -> 0}];

planeWaveCoefficient[expression_] := Module[{rules, reduced},
  rules = waveFieldHeadRules[{0, transverseAmplitude
      Exp[I (waveNumber xOne - frequency time)], 0}, 0, 0];
  reduced = expression /. rules;
  Together[reduced[[2]]/
    (transverseAmplitude Exp[I (waveNumber xOne - frequency time)])]
];

uniformS11b = uniformS11bModel[];
uniformS11bKernel = simplifyWeakKernel[extractCoupling[uniformS11b]];

uniformLimitS11cbOperand = <||>;
uniformLimitS11bOperand = <||>;
Do[Module[{key, model, operatorLimit, kernelLimit},
  key = caseKey[branch, density];
  model = evaluatedModel["EULERIAN", branch, density,
    "BASE", 0, False, True];
  operatorLimit = uniformLimit[model["OPERATOR"]];
  kernelLimit = simplifyWeakKernel[
    uniformLimit[extractCoupling[model]]];
  AssociateTo[uniformLimitS11cbOperand, key -> <|
    "SLAB_OPERATOR" -> operatorLimit,
    "COUPLING_KERNEL" -> kernelLimit,
    "TRANSVERSE_DISPERSION" -> planeWaveCoefficient[
      operatorLimit["U_MOMENTUM_ROWS"]]|>];
  AssociateTo[uniformLimitS11bOperand, key -> <|
    "SLAB_OPERATOR" -> uniformS11b["OPERATOR"],
    "COUPLING_KERNEL" -> uniformS11bKernel,
    "TRANSVERSE_DISPERSION" -> planeWaveCoefficient[
      uniformS11b["OPERATOR"]["U_MOMENTUM_ROWS"]]|>];
  Clear[model, operatorLimit, kernelLimit];
  ClearSystemCache[];
  Share[]],
  {branch, branches}, {density, densities}];

uniformPackageDimensions = <|"SLAB_OPERATOR" -> modelDimensions,
  "COUPLING_KERNEL" -> kernelDimensions,
  "TRANSVERSE_DISPERSION" -> dimensionBulkForce - dimensionU|>;
uniformPackageEpsilonOrders = <|"SLAB_OPERATOR" -> 1,
  "COUPLING_KERNEL" -> 1, "TRANSVERSE_DISPERSION" -> 0|>;
uniformPackageRecord[package_Association] := AssociationMap[
  Function[object, preparedObjectRecord[package[object],
    uniformPackageDimensions[object],
    uniformPackageEpsilonOrders[object]]], Keys[package]];
uniformPackageResidualRecord[left_Association,
    right_Association] := AssociationMap[Function[object,
  preparedObjectRecord[mapDifference[left[object], right[object]],
    uniformPackageDimensions[object],
    uniformPackageEpsilonOrders[object]]], Keys[left]];

emitShared["UNIFORM_LIMIT_S11CB_OPERAND",
  Map[uniformPackageRecord, uniformLimitS11cbOperand]];
emitShared["UNIFORM_LIMIT_S11B_OPERAND",
  Map[uniformPackageRecord, uniformLimitS11bOperand]];
emitShared["UNIFORM_LIMIT_RESIDUAL", AssociationMap[Function[key, <|
  "RESIDUAL_A_MINUS_B" -> uniformPackageResidualRecord[
    uniformLimitS11cbOperand[key], uniformLimitS11bOperand[key]]|>],
  Keys[uniformLimitS11cbOperand]]];
Clear[uniformLimitS11cbOperand, uniformLimitS11bOperand, uniformS11b,
  uniformS11bKernel, mainModels, mainKernels];
ClearSystemCache[];

(* ---------------------------------------------------------------------- *)
(* Dimensions, homogeneity, and a source-level dimension corruption.     *)
(* ---------------------------------------------------------------------- *)

uniformInvariantDimensions = <|
  "CURL_U_SQUARED" -> 2 (dimensionGradient + dimensionU),
  "DIV_U_SQUARED" -> 2 (dimensionGradient + dimensionU),
  "THETA_SQUARED" -> 2 dimensionTheta,
  "ELOCAL_SQUARED" -> 2 dimensionEw,
  "THETA_ELOCAL" -> dimensionTheta + dimensionEw,
  "GRAD_THETA_SQUARED" -> 2 (dimensionGradient + dimensionTheta),
  "GRAD_ELOCAL_SQUARED" -> 2 (dimensionGradient + dimensionEw),
  "GRAD_THETA_DOT_GRAD_ELOCAL" -> 2 dimensionGradient +
    dimensionTheta + dimensionEw,
  "THETA_DIV_U" -> dimensionTheta + dimensionGradient + dimensionU,
  "ELOCAL_DIV_U" -> dimensionEw + dimensionGradient + dimensionU|>;

factorDimensionByName = Association[Map[
  #1["NAME"] -> #1["DIMENSION"] &, dofFactorSpecifications]];
deriveInvariantDimensionFromFactorContent[blueprint_Association] :=
  dimensionGradient + Total[Lookup[factorDimensionByName,
    blueprint["FACTOR_NAMES"]]];

newInvariantDimensions = Association[Flatten[Table[
  MapThread[#1 -> deriveInvariantDimensionFromFactorContent[#2] &,
    {newInvariantLabels[source, False], o3KroneckerBlueprints}],
  {source, profileSources}], 1]];

newInvariantDimensionPayload = Association[Flatten[Table[
  MapThread[Function[{label, blueprint}, label -> <|
    "FACTOR_CONTENT" -> Join[{"BACKGROUND_FIRST_JET"},
      blueprint["FACTOR_NAMES"]],
    "FACTOR_DIMENSIONS" -> Join[{dimensionGradient},
      Lookup[factorDimensionByName, blueprint["FACTOR_NAMES"]]],
    "DERIVED_DIMENSION_OPERAND" ->
      deriveInvariantDimensionFromFactorContent[blueprint],
    "STORED_DIMENSION_OPERAND" ->
      blueprint["STORED_INVARIANT_DIMENSION"],
    "RESIDUAL_DERIVED_MINUS_STORED" ->
      deriveInvariantDimensionFromFactorContent[blueprint] -
        blueprint["STORED_INVARIANT_DIMENSION"],
    "DIMENSION_L_T_M" -> dimensionZero,
    "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>],
    {newInvariantLabels[source, False], o3KroneckerBlueprints}],
  {source, profileSources}], 1]];
emitShared["ENERGY_BASIS_NEW_INVARIANT_DIMENSION_DERIVATION",
  newInvariantDimensionPayload];

coefficientDimensions = Join[
  Map[dimensionEnergy - # &, uniformInvariantDimensions],
  Map[dimensionEnergy - # &, newInvariantDimensions]];

dimensionsPayload = <|
  "BASE_DERIVATION" -> <|
    "ENERGY_FROM_KINETIC" -> relationalObject[dimensionEnergy,
      dimensionRhoBr + 2 (dimensionU - dimensionTime)],
    "RHO4_FROM_RHOBR_AND_WIDTH" -> relationalObject[dimensionRhoFour,
      dimensionRhoBr - dimensionL],
    "MU_THETA_FROM_ENERGY_VARIATION" -> relationalObject[
      dimensionMuTheta, dimensionEnergy - dimensionTheta]|>,
  "INHERITED_AND_BACKGROUND_SYMBOLS" -> <|
    "c_s0" -> dimensionVelocity, "mu_R" -> dimensionEnergy,
    "rho_br" -> dimensionRhoBr, "W_0" -> dimensionL,
    "e_W" -> dimensionZero, "rho_m" -> dimensionRhoFour,
    "v_bulk_normal_0" -> dimensionVelocity,
    "mu_theta" -> dimensionMuTheta, "W_bg" -> dimensionL,
    "w1_profile" -> dimensionZero, "L_W" -> dimensionL,
    "sigma_W" -> dimensionZero, "mu_R_bg" -> dimensionEnergy,
    "m1_profile" -> dimensionZero,
    "rho_4D_bg_rho4_constant" -> dimensionRhoFour,
    "rho_br_bg_rho4_constant" -> dimensionRhoBr,
    "rho_4D_bg_rhobr_constant" -> dimensionRhoFour,
    "rho_br_bg_rhobr_constant" -> dimensionRhoBr,
    "e_W_bg" -> dimensionZero, "eta_bg" -> dimensionZero|>,
  "RESPONSE_KERNELS" -> <|
    "Lambda_A" -> dimensionLambdaA,
    "Lambda_V" -> dimensionLambdaV,
    "Lambda_X" -> dimensionLambdaX,
    "tau_A" -> dimensionT, "tau_V" -> dimensionT,
    "tau_X" -> dimensionT,
    "frequency" -> -dimensionT|>,
  "ENERGY_INVARIANT_COEFFICIENTS" -> coefficientDimensions,
  "NAMED_OBJECTS" -> <|
    "ENERGY_BASIS_VARIABLE" -> dimensionEnergy,
    "ENERGY_BASIS_COUNT" -> dimensionZero,
    "ENERGY_BASIS_NEW_INVARIANTS" -> dimensionEnergy,
    "ENERGY_BASIS_OMISSIONS" -> dimensionZero,
    "SLAB_OPERATOR" -> modelDimensions,
    "SLAB_OPERATOR_TERM_ORIGINS" -> modelDimensions,
    "MU_THETA_OPERATOR" -> dimensionMuTheta,
    "COUPLING_KERNEL" -> kernelDimensions,
    "COUPLING_KERNEL_TERM_ORIGINS" -> kernelDimensions,
    "ADMISSIBILITY_OPERATOR_OPERAND" -> <|
      "BODY" -> dimensionBulkForce, "FACE" -> dimensionFaceTraction|>,
    "ADMISSIBILITY_SUPPORT_OPERAND" -> <|
      "BODY" -> dimensionBulkForce, "FACE" -> dimensionFaceTraction|>,
    "ADMISSIBILITY_RESIDUAL" -> <|
      "BODY" -> dimensionBulkForce, "FACE" -> dimensionFaceTraction|>,
    "REPRESENTATION_PACKAGE" -> <|"OPERATOR" -> modelDimensions,
      "KERNEL" -> kernelDimensions|>,
    "CONTROL_PACKAGES" -> <|"OPERATOR" -> modelDimensions,
      "KERNEL" -> kernelDimensions, "ENERGY" -> dimensionEnergy|>,
    "UNIFORM_LIMIT_PACKAGE" -> <|"OPERATOR" -> modelDimensions,
      "KERNEL" -> kernelDimensions,
      "TRANSVERSE_DISPERSION" -> dimensionBulkForce - dimensionU|>,
    "HOMOGENEITY_PACKAGE" -> dimensionZero|>,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>;
emitShared["DIMENSIONS", dimensionsPayload];

homogeneityTermDimensions[coefficientDimensionMap_Association] := <|
  "ENERGY_TERMS" -> AssociationMap[
    coefficientDimensionMap[#] + Join[uniformInvariantDimensions,
      newInvariantDimensions][#] &, Keys[coefficientDimensionMap]],
  "U_MOMENTUM_ROWS" -> {dimensionRhoBr + dimensionU -
      2 dimensionTime, dimensionBulkForce, dimensionFaceTraction},
  "THICKNESS_ROW" -> {dimensionMuW + 2 dimensionL + dimensionEw -
      2 dimensionTime, dimensionThicknessRow, dimensionThicknessRow},
  "MASS_EVOLUTION_ROW" -> {dimensionRhoBr - dimensionTime,
      dimensionMassEvolution, dimensionFlux},
  "MU_THETA" -> {dimensionEnergy - dimensionTheta,
      dimensionMuTheta},
  "ADMISSIBILITY_PAIRING" -> {dimensionBulkForce,
      dimensionFaceTraction}|>;

dimensionMuW = dimensionEnergy - 2 dimensionVelocity;
baseHomogeneityDimensions = homogeneityTermDimensions[coefficientDimensions];
corruptedCoefficientDimensions = Association[coefficientDimensions];
firstNewDimensionKey = First[Keys[newInvariantDimensions]];
AssociateTo[corruptedCoefficientDimensions,
  firstNewDimensionKey ->
    (corruptedCoefficientDimensions[firstNewDimensionKey] + dimensionL)];
controlHomogeneityDimensions =
  homogeneityTermDimensions[corruptedCoefficientDimensions];

dimensionVectorQ[value_] := ListQ[value] && Length[value] == 3 &&
  VectorQ[value, IntegerQ];
homogeneityRelations[dimensionObject_Association] := Module[{values},
  values = Values[dimensionObject];
  If[values =!= {} && AllTrue[values, dimensionVectorQ],
    AssociationThread[Keys[dimensionObject],
      (relationalObject[#, First[values]] & /@ values)],
    Map[Function[value, If[AssociationQ[value],
      homogeneityRelations[value],
      relationalObject[#, First[value]] & /@ value]], dimensionObject]]
];

homogeneityBasePayload = <|
  "SOURCE_COEFFICIENT_DIMENSIONS" -> coefficientDimensions,
  "TERM_DIMENSIONS" -> baseHomogeneityDimensions,
  "RELATIONAL_OBJECTS" -> homogeneityRelations[baseHomogeneityDimensions],
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>;
emitShared["HOMOGENEITY_BASE_OPERAND", homogeneityBasePayload];

homogeneityControlPayload = <|
  "SOURCE_COEFFICIENT_DIMENSIONS" -> corruptedCoefficientDimensions,
  "CORRUPTED_SOURCE_KEY" -> firstNewDimensionKey,
  "TERM_DIMENSIONS" -> controlHomogeneityDimensions,
  "RELATIONAL_OBJECTS" ->
    homogeneityRelations[controlHomogeneityDimensions],
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>;
emitShared["HOMOGENEITY_CONTROL_OPERAND", homogeneityControlPayload];

emitShared["HOMOGENEITY_RESIDUAL", <|
  "OPERAND_A" -> homogeneityBasePayload,
  "OPERAND_B" -> homogeneityControlPayload,
  "RESIDUAL_A_MINUS_B" -> mapDifference[
    baseHomogeneityDimensions, controlHomogeneityDimensions],
  "TEST_OBJECT" -> relationalObject[baseHomogeneityDimensions,
    controlHomogeneityDimensions],
  "DIMENSION_L_T_M" -> dimensionZero,
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}}|>];

emitLocal["TAG_NAMES", <|"EXPRESSION" -> Append[localNames,
    "WL_S11CB_LOCAL_TAG_NAMES"],
  "MULTIGRADE_EPSILON_ETA_SIGMAW" -> {{0, 0, 0}},
  "DIMENSION_L_T_M" -> dimensionZero|>];
