(* S11b-A: blind, standalone Mathematica audit. *)
$HistoryLength = 0;
ClearAll["Global`*"];

(* Output and normalization helpers. *)
stripConditionalZero[expr_] := expr //. ConditionalExpression[0, _] :> 0;
inactiveIntegralRules = {
  HoldPattern[Inactive[Integrate][0, {_, _, _}]] :> 0,
  HoldPattern[Inactive[Integrate][_, {_, same_, same_}]] :> 0
};
clean[expr_] := stripConditionalZero[FullSimplify[expr /. inactiveIntegralRules]];
emit[tag_String, value_] :=
  Print["WL_" <> tag <> ": " <> ToString[value, InputForm, PageWidth -> Infinity]];
int[expr_, {var_, same_, same_}] := 0;
int[expr_, iterator_] := Inactive[Integrate][expr, iterator];

(* All dimensional vectors below are ordered {L,T,M}. *)
dimL = {1, 0, 0};
dimT = {0, 1, 0};
dimM = {0, 0, 1};
dimless = {0, 0, 0};

assumptions = Element[{rhoM, cs0, omega, k, qR, alpha, a, bigL, W0}, Reals] &&
  rhoM > 0 && cs0 > 0 && omega > 0 && k >= 0 && qR > 0 && alpha > 0 &&
  a > 0 && bigL > 0 && W0 > 0;
closureAssumptions = assumptions && Element[{lambdaP0, lambdaV0, x}, Reals] && x >= 0;

(* A1: obtain the finite source by the product rule/integration-by-parts identity. *)
productRuleResidual = Expand[D[om[w] jw[w], w] - om[w] D[jw[w], w] - D[om[w], w] jw[w]];
boundaryProduct[lo_, hi_] := om[lo] jw[lo] - om[hi] jw[hi];
weightedNormalCurrent[lo_, hi_] := int[D[om[w], w] jw[w], {w, lo, hi}];
sourceFinite = boundaryProduct[w1, w2] + weightedNormalCurrent[w1, w2];
sourceInfinite = int[D[om[w], w] jw[w], {w, -Infinity, Infinity}];
infiniteEndpointCondition = HoldForm[
  Limit[om[w] jw[w], w -> -Infinity] == 0 && Limit[om[w] jw[w], w -> Infinity] == 0];
projectionFinite = <|
  "projectedBalance" -> (dt[rho3] + divX[j3] == sourceFinite),
  "rho3" -> int[om[w] rho[w], {w, w1, w2}],
  "j3" -> int[om[w] jX[w], {w, w1, w2}],
  "source" -> sourceFinite
|>;
projectionInfinite = <|
  "endpointCondition" -> infiniteEndpointCondition,
  "source" -> sourceInfinite
|>;

(* A2: exact parity reduction on the explicitly symmetric finite interval [-L,L]. *)
evenBoundaryReduced = clean[om[bigL] jEven[bigL] - om[bigL] jEven[bigL]];
evenWeightedCurrentReduced = clean[
  int[om'[u] jEven[u], {u, 0, bigL}] - int[om'[u] jEven[u], {u, 0, bigL}]];
sourceEvenSymmetric = clean[evenBoundaryReduced + evenWeightedCurrentReduced];
sourceOddSymmetric = -2 om[bigL] jOdd[bigL] + 2 int[om'[u] jOdd[u], {u, 0, bigL}];
sourceEvenInfinite = clean[
  int[om'[u] jEven[u], {u, 0, Infinity}] - int[om'[u] jEven[u], {u, 0, Infinity}]];
sourceOddInfinite = 2 int[om'[u] jOdd[u], {u, 0, Infinity}];

(* Polynomial witnesses verify the parity reductions without assigning either result. *)
omWitness[w_] := 1 + w^2;
jEvenWitness[w_] := 2 + w^4;
jOddWitness[w_] := w + w^3;
finiteSourceConcrete[window_, current_, lo_, hi_] :=
  clean[(window /. w -> lo) (current /. w -> lo) -
    (window /. w -> hi) (current /. w -> hi) +
    Integrate[D[window, w] current, {w, lo, hi}]];
evenParityCheck = TrueQ[finiteSourceConcrete[omWitness[w], jEvenWitness[w], -bigL, bigL] == 0];
oddParityDirect = finiteSourceConcrete[omWitness[w], jOddWitness[w], -bigL, bigL];
oddParityReduced = clean[-2 omWitness[bigL] jOddWitness[bigL] +
  2 Integrate[omWitness'[u] jOddWitness[u], {u, 0, bigL}]];
oddParityCheck = TrueQ[clean[oddParityDirect - oddParityReduced] == 0];

parityInterval = <|
  "finiteInterval" -> {-bigL, bigL},
  "symmetricAboutW0" -> True,
  "finiteStatementsExactOn" -> HoldForm[{-bigL, bigL}],
  "infiniteStatementsExactOn" -> HoldForm[{-Infinity, Infinity}],
  "infiniteEndpointCondition" -> HoldForm[om[w] jw[w] -> 0]
|>;
parityEven = <|
  "finiteSymmetricSource" -> sourceEvenSymmetric,
  "infiniteSymmetricSource" -> sourceEvenInfinite,
  "exact" -> evenParityCheck
|>;
parityOdd = <|
  "finiteSymmetricSource" -> sourceOddSymmetric,
  "infiniteSymmetricSource" -> sourceOddInfinite,
  "exact" -> oddParityCheck,
  "notForcedToVanish" -> TrueQ[oddParityReduced =!= 0]
|>;

(* A3: product rules in t and the three tangential coordinates. *)
dynamicCoordinates = {x1, x2, x3};
dynamicArguments = {w, x1, x2, x3, t};
dynamicExtraTime = int[(rhoDyn @@ dynamicArguments) D[omDyn @@ dynamicArguments, t], {w, w1, w2}];
dynamicExtraTangential = int[
  Sum[jDyn[i, w, x1, x2, x3, t] D[omDyn @@ dynamicArguments, dynamicCoordinates[[i]]], {i, 1, 3}],
  {w, w1, w2}];
dynamicTimeProductTerm = Expand[D[omTime[t] rhoTime[t], t] - omTime[t] D[rhoTime[t], t]];
dynamicTangentialProductTerm = Expand[
  D[omTangential[x1] jTangential[x1], x1] - omTangential[x1] D[jTangential[x1], x1]];
dynamicExtraSigns = clean[{
  dynamicTimeProductTerm/(rhoTime[t] omTime'[t]),
  dynamicTangentialProductTerm/(jTangential[x1] omTangential'[x1])
}];
dynamicWindowExtras = <|
  "timeProductRuleTerm" -> dynamicExtraTime,
  "tangentialProductRuleTerm" -> dynamicExtraTangential,
  "numberOfTermsAbsentFromStaticA1" -> Length[{dynamicExtraTime, dynamicExtraTangential}],
  "signsInProjectedBalance" -> dynamicExtraSigns
|>;

(* A4: harmonic acoustic half-space solution and radiation/decay branch selection. *)
qSquared = omega^2/cs0^2 - k^2;
upperPotential = ampU Exp[I q (w - W0/2)];
lowerPotential = ampD Exp[-I q (w + W0/2)];
upperNormalWaveNumber = q;
lowerNormalWaveNumber = -q;

(* Time-averaged normal acoustic energy-flux signs for a real normal wave number. *)
energyFluxSign[qw_] := FullSimplify[Sign[rhoM omega qw/2], assumptions];
propagatingBranchCheck = {
  TrueQ[energyFluxSign[qR] == 1],
  TrueQ[energyFluxSign[-qR] == -1]
};
evanescentBranchCheck = {
  TrueQ[FullSimplify[Im[I alpha] > 0, assumptions]],
  TrueQ[FullSimplify[Im[-I alpha] < 0, assumptions]]
};
radiationBranches = <|
  "harmonicFactor" -> HoldForm[Exp[I (kVec . xVec - omega t)]],
  "qSquared" -> qSquared,
  "propagating_qSquaredPositive" -> <|"upperNormalWaveNumber" -> qR, "lowerNormalWaveNumber" -> -qR,
    "upperEnergyDirection" -> energyFluxSign[qR], "lowerEnergyDirection" -> energyFluxSign[-qR]|>,
  "evanescent_qSquaredNegative" -> <|"upperNormalWaveNumber" -> I alpha, "lowerNormalWaveNumber" -> -I alpha,
    "upperDecayDirection" -> FullSimplify[Sign[Im[I alpha]], assumptions],
    "lowerDecayDirection" -> FullSimplify[Sign[Im[-I alpha]], assumptions]|>,
  "grazing_qSquaredZero" -> <|"upperNormalWaveNumber" -> 0, "lowerNormalWaveNumber" -> 0|>,
  "noIncomingWaves" -> And @@ Join[propagatingBranchCheck, evanescentBranchCheck]
|>;

(* At either face, outward fluid velocity is i q A and pressure is i rho omega A. *)
impedanceSolve = First[Solve[I q amp == faceVelocity, amp]];
pressureAmplitude = I rhoM omega amp;
zImpermeable = clean[(pressureAmplitude/faceVelocity) /. impedanceSolve];
zPropagating = clean[zImpermeable /. q -> qR];
zEvanescent = clean[zImpermeable /. q -> I alpha];
zPropParts = clean[ComplexExpand[{Re[zPropagating], Im[zPropagating]}]];
zEvanParts = clean[ComplexExpand[{Re[zEvanescent], Im[zEvanescent]}]];

(* Radiation leaves one amplitude in each half-space; the two face conditions fix both. *)
zetaThickness = {deltaW/2, -deltaW/2};
zetaCenter = {zetaC, zetaC};
outwardSigns = {1, -1};
outwardVelocity[zetas_List] := clean[outwardSigns (-I omega zetas)];
pressureFor[zetas_List, zLocal_] := clean[zLocal outwardVelocity[zetas]];
radiatingAmplitudes = {ampU, ampD};
bulkPotentials = {upperPotential, lowerPotential};
outwardBulkFaceVelocities = clean[{
  D[upperPotential, w] /. w -> W0/2,
  -D[lowerPotential, w] /. w -> -W0/2
}];
bulkFacePressures = clean[{
  I rhoM omega upperPotential /. w -> W0/2,
  I rhoM omega lowerPotential /. w -> -W0/2
}];
faceConditions[zetas_List] := Thread[outwardBulkFaceVelocities == outwardVelocity[zetas]];
amplitudeSolution[zetas_List] := First[Solve[faceConditions[zetas], radiatingAmplitudes]];
solvedFacePressures[zetas_List] := clean[bulkFacePressures /. amplitudeSolution[zetas]];
solvedPerFaceResponse[zetas_List] := clean[
  solvedFacePressures[zetas]/outwardVelocity[zetas]];
solvedModeResponse[zetas_List] := <|
  "motionGlobalW" -> zetas,
  "outwardFaceVelocities" -> outwardVelocity[zetas],
  "amplitudeSolution" -> amplitudeSolution[zetas],
  "facePressures" -> solvedFacePressures[zetas],
  "ZPerFace" -> solvedPerFaceResponse[zetas]
|>;
thicknessModeResponse = solvedModeResponse[zetaThickness];
centerModeResponse = solvedModeResponse[zetaCenter];

(* This probe is neither a pure thickness motion nor a pure centre shift. *)
nonCombinationMotion = {zetaProbe, 2 zetaProbe};
thicknessProbeMotion = {zetaProbe, -zetaProbe};
centerProbeMotion = {zetaProbe, zetaProbe};
nonCombinationModeResponse = solvedModeResponse[nonCombinationMotion];
nonCombinationChangesSolvedQuantities = And[
  amplitudeSolution[nonCombinationMotion] =!= amplitudeSolution[thicknessProbeMotion],
  amplitudeSolution[nonCombinationMotion] =!= amplitudeSolution[centerProbeMotion],
  solvedFacePressures[nonCombinationMotion] =!= solvedFacePressures[thicknessProbeMotion],
  solvedFacePressures[nonCombinationMotion] =!= solvedFacePressures[centerProbeMotion]
];
perFaceResponseDependsOnCombination = Not[TrueQ[
  clean[solvedPerFaceResponse[zetaThickness] - solvedPerFaceResponse[zetaCenter]] == {0, 0} &&
  clean[solvedPerFaceResponse[zetaThickness] - solvedPerFaceResponse[nonCombinationMotion]] == {0, 0}
]];
radiationAmplitudeCount = Length[radiatingAmplitudes];
faceConditionCount = Length[faceConditions[{zetaU, zetaD}]];
allRadiatingAmplitudesFixed = TrueQ[
  Sort[First /@ amplitudeSolution[{zetaU, zetaD}]] === Sort[radiatingAmplitudes]];

(* Added mass is deliberately computed against global-w face acceleration, per face. *)
globalFaceAcceleration[zetas_List] := clean[-omega^2 zetas];
addedMassByFace = clean[pressureFor[{zetaU, zetaD}, zEvanescent]/globalFaceAcceleration[{zetaU, zetaD}]];
addedMassExpected = {rhoM/alpha, -rhoM/alpha};
addedMassCheck = TrueQ[clean[addedMassByFace - addedMassExpected] == {0, 0}];
addedMassSigns = FullSimplify[Sign /@ addedMassByFace, assumptions];
addedMassOutput = <|
  "purelyImaginaryRegime" -> HoldForm[q^2 < 0],
  "perFaceUpperLowerGlobalWDefinition" -> addedMassByFace,
  "deltaWPerFace" -> addedMassByFace,
  "zetaCPerFace" -> addedMassByFace,
  "signsUpperLower" -> addedMassSigns,
  "check" -> addedMassCheck
|>;

(* The prescribed pole test uses 1/expression == 0. *)
grazingPoleTest = TrueQ[FullSimplify[(1/zImpermeable == 0) /. q -> 0, assumptions]];
grazingAddedMassPoleTest = And @@ (TrueQ[FullSimplify[(1/# == 0) /. alpha -> 0, assumptions]] & /@ addedMassExpected);
grazingFaceResidual = clean[I q amp - faceVelocity];
grazingDrivenSolutionSet = Reduce[
  ((grazingFaceResidual /. q -> 0) == 0) && faceVelocity != 0, amp, Complexes];
grazingUndrivenResidual = clean[grazingFaceResidual /. {q -> 0, faceVelocity -> 0}];
grazingUndrivenAmplitudeCoefficient = clean[Coefficient[grazingUndrivenResidual, amp]];
grazingDrivenStatus = If[grazingPoleTest && TrueQ[grazingDrivenSolutionSet === False],
  "NO_SOLUTION_WITH_FINITE_BULK_AMPLITUDE", "NOT_ESTABLISHED"];
grazingUndrivenStatus = If[
  TrueQ[grazingUndrivenAmplitudeCoefficient == 0 && (grazingUndrivenResidual /. amp -> 0) == 0],
  "CONSTANT_BULK_AMPLITUDE_FREE", "NOT_ESTABLISHED"];
propagatingRePoleTest = TrueQ[FullSimplify[(1/zPropParts[[1]] == 0) /. qR -> 0, assumptions]];
evanescentImPoleTest = TrueQ[FullSimplify[(1/zEvanParts[[2]] == 0) /. alpha -> 0, assumptions]];
grazingBehavior = <|
  "q" -> 0,
  "finiteAmplitudeOutwardVelocity" -> clean[(I q amp) /. q -> 0],
  "drivenNonzeroFaceVelocity" -> grazingDrivenStatus,
  "undrivenFace" -> grazingUndrivenStatus,
  "ZPoleByReciprocalTest" -> grazingPoleTest,
  "propagatingApproach" -> <|"ReZPoleByReciprocalTest" -> propagatingRePoleTest,
    "ImZ" -> zPropParts[[2]]|>,
  "evanescentApproach" -> <|"ReZ" -> zEvanParts[[1]],
    "ImZPoleByReciprocalTest" -> evanescentImPoleTest|>,
  "addedMassPerFaceReciprocalsAtGrazing" -> clean[(1/addedMassExpected) /. alpha -> 0],
  "addedMassSignsUpperLower" -> addedMassSigns,
  "addedMassPoleByReciprocalTest" -> grazingAddedMassPoleTest
|>;

zByRegime = <|
  "qSquaredPositive" -> <|"q" -> qR, "Z" -> zPropagating, "ReIm" -> zPropParts|>,
  "qSquaredNegative" -> <|"q" -> I alpha, "Z" -> zEvanescent, "ReIm" -> zEvanParts|>,
  "qSquaredZero" -> grazingBehavior
|>;
zImpermeableOutput = <|"branches" -> radiationBranches, "genericZ" -> zImpermeable|>;
genericSolvedPerFaceResponse = solvedPerFaceResponse[{zetaU, zetaD}];
propagatingSolvedPerFaceResponse = clean[genericSolvedPerFaceResponse /. q -> qR];
evanescentSolvedPerFaceResponse = clean[genericSolvedPerFaceResponse /. q -> I alpha];
perFaceRealImaginaryParts[responses_List] :=
  (clean[ComplexExpand[{Re[#], Im[#]}]] & /@ responses);
parityRegimeResponse = <|
  "qSquaredPositive" -> <|"ZPerFace" -> propagatingSolvedPerFaceResponse,
    "ReImPerFace" -> perFaceRealImaginaryParts[propagatingSolvedPerFaceResponse]|>,
  "qSquaredNegative" -> <|"ZPerFace" -> evanescentSolvedPerFaceResponse,
    "ReImPerFace" -> perFaceRealImaginaryParts[evanescentSolvedPerFaceResponse]|>,
  "qSquaredZero" -> grazingBehavior
|>;
zByParity = <|
  "amplitudesAfterRadiationCondition" -> radiatingAmplitudes,
  "amplitudeCountAfterRadiationCondition" -> radiationAmplitudeCount,
  "faceConditionCount" -> faceConditionCount,
  "faceConditionsFixAllAmplitudes" -> allRadiatingAmplitudesFixed,
  "deltaW" -> Append[thicknessModeResponse, "byRegime" -> parityRegimeResponse],
  "zetaC" -> Append[centerModeResponse, "byRegime" -> parityRegimeResponse],
  "perFaceResponseDependsOnCombination" -> perFaceResponseDependsOnCombination,
  "nonCombinationAblation" -> <|
    "motionGlobalW" -> nonCombinationMotion,
    "solvedResponse" -> nonCombinationModeResponse,
    "changesSolvedAmplitudesAndPressures" -> nonCombinationChangesSolvedQuantities
  |>
|>;

(* Relative-flux channels follow from outward face orientation, not an invertible relabelling. *)
faceFluxes = {jPlus, jMinus};
netAccretion = clean[-Total[faceFluxes]];
throughFlowLowerToUpper = clean[(outwardSigns . faceFluxes)/2];
fluxChannelVector = {netAccretion, throughFlowLowerToUpper};
fluxChannelMatrix = Table[Coefficient[fluxChannelVector[[i]], faceFluxes[[j]]],
  {i, Length[fluxChannelVector]}, {j, Length[faceFluxes]}];
fluxOrientationMatrix = {
  -ConstantArray[1, Length[faceFluxes]],
  outwardSigns/Length[faceFluxes]
};
fluxInverse = Thread[faceFluxes -> clean[Inverse[fluxChannelMatrix] . {acc, thru}]];
jPlusDefinition = HoldForm[rhoM (vWUpper - zetaUDot)];
jMinusDefinition = HoldForm[-rhoM (vWLower - zetaDDot)];
fluxChannelCheck = TrueQ[fluxChannelMatrix === fluxOrientationMatrix &&
  Det[fluxChannelMatrix] =!= 0];
relativeFluxChannels = <|
  "JplusDefinition" -> jPlusDefinition,
  "JminusDefinition" -> jMinusDefinition,
  "netAccretionBySlab" -> netAccretion,
  "throughFlowOrientedLowerToUpper" -> throughFlowLowerToUpper,
  "inverse" -> fluxInverse
|>;

(* A5: one shared, finite relaxation time; x denotes omega tau. *)
relaxationDenominator = 1 - I x;
lambdaP = lambdaP0/relaxationDenominator;
lambdaV = lambdaV0/relaxationDenominator;
acousticOutwardVelocityFromPressure = q pressure/(rhoM omega);
closureEquation = rhoM (acousticOutwardVelocityFromPressure - faceVelocity) ==
  lambdaP pressure + lambdaV faceVelocity;
permeablePressureSolve = First[Solve[closureEquation, pressure]];
zPermeable = clean[(pressure/faceVelocity) /. permeablePressureSolve];
zPermeableExpected = clean[(rhoM relaxationDenominator + lambdaV0)/
  ((q/omega) relaxationDenominator - lambdaP0)];
zPermeableCheck = TrueQ[clean[zPermeable - zPermeableExpected] == 0];

zPermProp = clean[zPermeable /. q -> omega r];
zPermEvan = clean[zPermeable /. q -> I omega s];
zPermGrazing = clean[zPermeable /. q -> 0];
reProp = clean[ComplexExpand[Re[zPermProp]]];
imProp = clean[ComplexExpand[Im[zPermProp]]];
reEvan = clean[ComplexExpand[Re[zPermEvan]]];
imEvan = clean[ComplexExpand[Im[zPermEvan]]];
reGrazing = clean[ComplexExpand[Re[zPermGrazing]]];
imGrazing = clean[ComplexExpand[Im[zPermGrazing]]];
signIndefiniteQ[expr_, positivePoint_List, negativePoint_List] :=
  TrueQ[N[expr /. positivePoint] > 0 && N[expr /. negativePoint] < 0];
propSignIndefinite = signIndefiniteQ[reProp,
  {rhoM -> 1, r -> 1, x -> 0, lambdaV0 -> 0, lambdaP0 -> 0},
  {rhoM -> 1, r -> 1, x -> 0, lambdaV0 -> 0, lambdaP0 -> 2}];
evanSignIndefinite = signIndefiniteQ[reEvan,
  {rhoM -> 1, s -> 1, x -> 0, lambdaV0 -> 0, lambdaP0 -> -1},
  {rhoM -> 1, s -> 1, x -> 0, lambdaV0 -> 0, lambdaP0 -> 1}];
grazingSignIndefinite = signIndefiniteQ[reGrazing,
  {rhoM -> 1, x -> 0, lambdaV0 -> 0, lambdaP0 -> -1},
  {rhoM -> 1, x -> 0, lambdaV0 -> 0, lambdaP0 -> 1}];

smallX[expr_] := clean[Limit[expr, x -> 0, Direction -> "FromAbove"]];
largeX[expr_] := clean[Limit[expr, x -> Infinity]];
largeXScaled[expr_] := clean[Limit[expr/x, x -> Infinity]];

inPhaseByRegime = <|
  "definitions" -> <|"r" -> HoldForm[qR/omega], "s" -> HoldForm[alpha/omega],
    "x" -> HoldForm[omega tau]|>,
  "qSquaredPositive" -> <|
    "deltaW" -> reProp, "zetaC" -> reProp,
    "genericPresence" -> TrueQ[reProp =!= 0],
    "radiativeBaselineAtZeroCoefficients" -> clean[reProp /. {lambdaP0 -> 0, lambdaV0 -> 0}],
    "dependsOnLambdaP0" -> TrueQ[D[reProp, lambdaP0] =!= 0],
    "dependsOnLambdaV0" -> TrueQ[D[reProp, lambdaV0] =!= 0],
    "zeroLocus" -> (Numerator[Together[reProp]] == 0),
    "signDefiniteForFreeRealCoefficients" -> Not[propSignIndefinite]|>,
  "qSquaredNegative" -> <|
    "deltaW" -> reEvan, "zetaC" -> reEvan,
    "genericPresence" -> TrueQ[reEvan =!= 0],
    "baselineAtZeroCoefficients" -> clean[reEvan /. {lambdaP0 -> 0, lambdaV0 -> 0}],
    "dependsOnLambdaP0" -> TrueQ[D[reEvan, lambdaP0] =!= 0],
    "dependsOnLambdaV0" -> TrueQ[D[reEvan, lambdaV0] =!= 0],
    "zeroLocus" -> (Numerator[Together[reEvan]] == 0),
    "signDefiniteForFreeRealCoefficients" -> Not[evanSignIndefinite]|>,
  "qSquaredZero_lambdaP0Nonzero" -> <|
    "deltaW" -> reGrazing, "zetaC" -> reGrazing,
    "genericPresence" -> TrueQ[reGrazing =!= 0],
    "dependsOnLambdaP0" -> TrueQ[D[reGrazing, lambdaP0] =!= 0],
    "dependsOnLambdaV0" -> TrueQ[D[reGrazing, lambdaV0] =!= 0],
    "zeroLocus" -> (Numerator[Together[reGrazing]] == 0),
    "signDefiniteForFreeRealCoefficients" -> Not[grazingSignIndefinite]|>
|>;

dissipationVsX = <|
  "definitions" -> <|"r" -> HoldForm[qR/omega], "s" -> HoldForm[alpha/omega],
    "x" -> HoldForm[omega tau]|>,
  "qSquaredPositive" -> <|"smallX_Z" -> smallX[zPermProp], "orderUnity_Z" -> zPermProp,
    "largeX_Z" -> largeX[zPermProp], "smallX_ReZ" -> smallX[reProp],
    "orderUnity_ReZ" -> reProp, "largeX_ReZ" -> largeX[reProp]|>,
  "qSquaredNegative" -> <|"smallX_Z" -> smallX[zPermEvan], "orderUnity_Z" -> zPermEvan,
    "largeX_Z" -> largeX[zPermEvan], "smallX_ReZ" -> smallX[reEvan],
    "orderUnity_ReZ" -> reEvan, "largeX_ReZ" -> largeX[reEvan]|>,
  "qSquaredZero_lambdaP0Nonzero" -> <|"smallX_Z" -> smallX[zPermGrazing], "orderUnity_Z" -> zPermGrazing,
    "largeXLeadingCoefficientOfX" -> largeXScaled[zPermGrazing], "smallX_ReZ" -> smallX[reGrazing],
    "orderUnity_ReZ" -> reGrazing, "largeX_ReZ" -> largeX[reGrazing]|>
|>;

tauZeroLimit = <|
  "qSquaredPositive" -> smallX[zPermProp],
  "qSquaredNegative" -> smallX[zPermEvan],
  "qSquaredZero_lambdaP0Nonzero" -> smallX[zPermGrazing],
  "status" -> "SPECIAL_MEMORYLESS_LIMIT_NOT_PREMISE"
|>;

(* Degeneracy is the vanishing coefficient of the bulk pressure/amplitude. *)
bulkAmplitudeCoefficient = clean[(q/omega) relaxationDenominator - lambdaP0];
driveCoefficient = clean[rhoM relaxationDenominator + lambdaV0];
propCoefficientParts = clean[ComplexExpand[{Re[bulkAmplitudeCoefficient /. q -> omega r],
  Im[bulkAmplitudeCoefficient /. q -> omega r]}]];
evanCoefficientParts = clean[ComplexExpand[{Re[bulkAmplitudeCoefficient /. q -> I omega s],
  Im[bulkAmplitudeCoefficient /. q -> I omega s]}]];
propDegenerate = Reduce[And @@ Thread[propCoefficientParts == {0, 0}] && r > 0 && x >= 0,
  {lambdaP0, x}, Reals];
evanDegenerate = Reduce[And @@ Thread[evanCoefficientParts == {0, 0}] && s > 0 && x >= 0,
  {lambdaP0, x}, Reals];
grazingDegenerate = Reduce[(bulkAmplitudeCoefficient /. q -> 0) == 0, lambdaP0, Reals];
driveCoefficientZero = Reduce[And @@ Thread[ComplexExpand[{Re[driveCoefficient], Im[driveCoefficient]}] == {0, 0}] &&
  rhoM > 0 && x >= 0, {lambdaV0, x}, Reals];
propDrivenNoSolution = Reduce[propDegenerate && lambdaV0 != -rhoM, {lambdaP0, lambdaV0, x}, Reals];
propDrivenFree = Reduce[propDegenerate && lambdaV0 == -rhoM, {lambdaP0, lambdaV0, x}, Reals];
grazingDrivenNoSolution = Reduce[grazingDegenerate && rhoM > 0 && x >= 0 && Not[driveCoefficientZero],
  {lambdaP0, lambdaV0, x}, Reals];
grazingDrivenFree = Reduce[grazingDegenerate && rhoM > 0 && x >= 0 && driveCoefficientZero,
  {lambdaP0, lambdaV0, x}, Reals];
regularizedClosureResidual[coefficient_, drive_, prescribedVelocity_] :=
  clean[coefficient pressure - drive prescribedVelocity];
undrivenFreeAmplitudeStatus[regimeCoefficient_, locus_, regimeAssumptions_] := Module[
  {undrivenResidual, coefficientVanishing, residualVanishing},
  undrivenResidual = regularizedClosureResidual[regimeCoefficient, driveCoefficient, 0];
  coefficientVanishing = TrueQ[FullSimplify[
    Coefficient[undrivenResidual, pressure] == 0, regimeAssumptions && locus]];
  residualVanishing = TrueQ[FullSimplify[
    (undrivenResidual /. pressure -> 0) == 0, regimeAssumptions && locus]];
  If[coefficientVanishing && residualVanishing, "FREE_AMPLITUDE", "NOT_ESTABLISHED"]
];
propUndrivenStatus = undrivenFreeAmplitudeStatus[
  bulkAmplitudeCoefficient /. q -> omega r, propDegenerate,
  Element[{r, lambdaP0, x}, Reals] && r > 0 && x >= 0];
grazingUndrivenClosureStatus = undrivenFreeAmplitudeStatus[
  bulkAmplitudeCoefficient /. q -> 0, grazingDegenerate,
  Element[{lambdaP0, x}, Reals] && x >= 0];
degenerateLoci = <|
  "qSquaredPositive" -> <|"locus" -> propDegenerate,
    "drivenNoSolutionLocus" -> propDrivenNoSolution,
    "drivenFreeAmplitudeLocus" -> propDrivenFree,
    "undriven" -> propUndrivenStatus|>,
  "qSquaredNegative" -> <|"locus" -> evanDegenerate|>,
  "qSquaredZero" -> <|"locus" -> grazingDegenerate,
    "drivenNoSolutionLocus" -> grazingDrivenNoSolution,
    "drivenFreeAmplitudeLocus" -> grazingDrivenFree,
    "undriven" -> grazingUndrivenClosureStatus|>
|>;

closureScopeStructuralCheck = TrueQ[Together[lambdaP/lambdaP0 - lambdaV/lambdaV0] == 0];
closureScopeStatus = If[closureScopeStructuralCheck && ! FreeQ[zPermeable, x],
  "ONE_SHARED_RELAXATION_TIME; SEPARATE_TIMES_AND_FULL_KERNEL_NOT_ATTEMPTED",
  "NOT_ESTABLISHED"];
closureScope = <|
  "sharedRelaxationTime" -> closureScopeStructuralCheck,
  "x" -> HoldForm[omega tau],
  "tauDomain" -> HoldForm[tau >= 0],
  "freeRealConstantsDomain" -> HoldForm[Element[{lambdaP0, lambdaV0, tau}, Reals] && tau >= 0],
  "scopeStatus" -> closureScopeStatus,
  "memorylessLaw" -> If[TrueQ[smallX[zPermProp] =!= zPermProp],
    "TAU_TO_ZERO_ONLY_AS_SPECIAL_LIMIT", "NOT_ESTABLISHED"]
|>;

(* A6: dimensional derivations from the stated equations and four-space geometry. *)
gradDim = -dimL;
dtDim = -dimT;
velocityDim = dimL - dimT;
phiDim = clean[velocityDim - gradDim];
rhoDim = clean[dimM - 4 dimL];
pressureDim = clean[rhoDim + dtDim + phiDim];

csUnknown = {csL, csT, csM};
csSolution = First[Solve[Thread[2 csUnknown + phiDim + 2 gradDim == phiDim + 2 dtDim], csUnknown]];
csDim = clean[csUnknown /. csSolution];
v0Dim = velocityDim;
jwDim = clean[rhoDim + velocityDim];
zDim = clean[pressureDim - velocityDim];
accelerationDim = dimL - 2 dimT;
addedMassDim = clean[pressureDim - accelerationDim];
lambdaPDim = clean[jwDim - pressureDim];
lambdaVDim = clean[jwDim - velocityDim];
omegaDim = -dimT;
tauUnknown = {tauL, tauT, tauM};
tauSolution = First[Solve[Thread[omegaDim + tauUnknown == dimless], tauUnknown]];
tauDim = clean[tauUnknown /. tauSolution];
sourceDimFromBoundary = jwDim;
sourceDimFromProjectedContinuity = clean[rhoDim + dimL + dtDim];
sourceDim = sourceDimFromBoundary;

dimensionValues = <|
  "Z" -> zDim, "m_add" -> addedMassDim, "rho_m" -> rhoDim, "c_s0" -> csDim,
  "v0" -> v0Dim, "Lambda_p0" -> lambdaPDim, "Lambda_V0" -> lambdaVDim,
  "tau" -> tauDim, "A1_source" -> sourceDim
|>;
dimensionRoutes = <|
  "Z" -> HoldForm[dim[zImpedance] == dim[pressure/velocity]],
  "m_add" -> HoldForm[dim[mAdd] == dim[pressure/acceleration]],
  "rho_m" -> HoldForm[dim[rhoM] == dim[mass/(length^4)]],
  "c_s0" -> HoldForm[dim[dt^2 phi] == dim[cs0^2 laplacian phi]],
  "v0" -> HoldForm[dim[v0] == dim[length/time]],
  "Lambda_p0" -> HoldForm[dim[lambdaP0] == dim[massFlux/pressure]],
  "Lambda_V0" -> HoldForm[dim[lambdaV0] == dim[massFlux/velocity]],
  "tau" -> HoldForm[dim[tau] == dim[1/omega]],
  "A1_source" -> HoldForm[dim[dt rho3 + divX[j3] - a1Source] == dim[0]]
|>;
dimensionRouteTargets = <|
  "Z" -> zImpedance, "m_add" -> mAdd, "rho_m" -> rhoM, "c_s0" -> cs0,
  "v0" -> v0, "Lambda_p0" -> lambdaP0, "Lambda_V0" -> lambdaV0,
  "tau" -> tau, "A1_source" -> a1Source
|>;
routeKindFromAssertion[route_, target_] := Module[
  {heldRoute, equationSides, heldTarget},
  heldRoute = route /. HoldForm[body_] :> HoldComplete[body];
  equationSides = Cases[heldRoute,
    HoldComplete[Equal[left_, right_]] :> {HoldComplete[left], HoldComplete[right]}, {0}];
  heldTarget = HoldComplete[dim[routeTargetPlaceholder]] /. routeTargetPlaceholder -> target;
  Which[
    Length[equationSides] == 1 && MemberQ[First[equationSides], heldTarget], "DEFINITIONAL",
    ! FreeQ[heldRoute, target], "INDEPENDENT",
    True, "NOT_ESTABLISHED"
  ]
];
routeKinds = AssociationMap[
  routeKindFromAssertion[dimensionRoutes[#], dimensionRouteTargets[#]] &,
  Keys[dimensionValues]];
dimensionChecks = {
  TrueQ[2 csDim + phiDim + 2 gradDim == phiDim + 2 dtDim],
  TrueQ[jwDim == lambdaPDim + pressureDim],
  TrueQ[jwDim == lambdaVDim + velocityDim],
  TrueQ[omegaDim + tauDim == dimless],
  TrueQ[sourceDimFromBoundary == sourceDimFromProjectedContinuity]
};

(* A7: the two prescribed form controls, retaining general even/odd currents. *)
omegaB[w_] := Sech[(w - b)/a]^2;
omegaBPrime[z_] := clean[-(2/a) Sech[(z - b)/a]^2 Tanh[(z - b)/a]];
omegaZero[w_] := clean[omegaB[w] /. b -> 0];
omegaZeroPrime[z_] := clean[-(2/a) Sech[z/a]^2 Tanh[z/a]];

windowControlEven = clean[(omegaB[-bigL] - omegaB[bigL]) jEven[bigL]] +
  int[(omegaBPrime[u] + omegaBPrime[-u]) jEven[u], {u, 0, bigL}];
windowControlOdd = clean[-(omegaB[-bigL] + omegaB[bigL]) jOdd[bigL]] +
  int[(omegaBPrime[u] - omegaBPrime[-u]) jOdd[u], {u, 0, bigL}];

(* Concrete parity-respecting witnesses establish non-identity in b. *)
windowWitnessSource[current_, bValue_] := NIntegrate[
  Evaluate[-omegaB[w] D[current, w] /. {a -> 1, bigL -> 2, b -> bValue}], {w, -2, 2}];
evenBChanged = Abs[windowWitnessSource[w^2, 1/2] - windowWitnessSource[w^2, 0]] > 10^-10;
oddBChanged = Abs[windowWitnessSource[w, 1/2] - windowWitnessSource[w, 0]] > 10^-10;

windowFamily = HoldForm[Sech[(w - b)/a]^2];
windowControlInterval = {-bigL, bigL};
windowBControlDomain = HoldForm[b != 0];
windowParityControl = <|
  "family" -> windowFamily,
  "interval" -> windowControlInterval,
  "evenCurrentSourceAsFunctionOfB" -> windowControlEven,
  "oddCurrentSourceAsFunctionOfB" -> windowControlOdd,
  "evenCurrentIdenticallyIndependentOfB" -> Not[evenBChanged],
  "oddCurrentIdenticallyIndependentOfB" -> Not[oddBChanged],
  "b0EvenSource" -> clean[windowControlEven /. b -> 0],
  "bControlDomain" -> windowBControlDomain
|>;

intervalControlEven = clean[omegaZero[bigL] jEven[bigL] - omegaZero[bigL + c] jEven[bigL + c]] +
  int[omegaZeroPrime[u] jEven[u], {u, bigL, bigL + c}];
intervalControlOdd = clean[-omegaZero[bigL] jOdd[bigL] - omegaZero[bigL + c] jOdd[bigL + c]] +
  2 int[omegaZeroPrime[u] jOdd[u], {u, 0, bigL}] +
  int[omegaZeroPrime[u] jOdd[u], {u, bigL, bigL + c}];
evenCDerivativeWitness = clean[-omegaZero[bigL + c] 2 (bigL + c) /. c -> 0];
oddCDerivativeWitness = clean[-omegaZero[bigL + c] /. c -> 0];
evenCChanged = TrueQ[FullSimplify[evenCDerivativeWitness != 0, assumptions]];
oddCChanged = TrueQ[FullSimplify[oddCDerivativeWitness != 0, assumptions]];
intervalControlFamily = {-bigL, bigL + c};
intervalControlWindow = HoldForm[Sech[w/a]^2];
intervalControlDomain = HoldForm[c != 0 && c > -2 bigL];
intervalSymmetryControl = <|
  "family" -> intervalControlFamily,
  "window" -> intervalControlWindow,
  "evenCurrentSourceAsFunctionOfC" -> intervalControlEven,
  "oddCurrentSourceAsFunctionOfC" -> intervalControlOdd,
  "evenCurrentIdenticallyIndependentOfC" -> Not[evenCChanged],
  "oddCurrentIdenticallyIndependentOfC" -> Not[oddCChanged],
  "c0EvenSource" -> clean[intervalControlEven /. c -> 0],
  "c0OddSource" -> clean[intervalControlOdd /. c -> 0],
  "controlDomain" -> intervalControlDomain
|>;

(* A8: independent dimensionless small parameters and convective expansion. *)
backgroundQuantityDim = rhoDim;
backgroundRateDim = clean[backgroundQuantityDim + dtDim];
tBgDim = clean[backgroundQuantityDim - backgroundRateDim];
timescaleParameter = clean[1/(omega tBg)];
flowSpeedParameter = clean[Abs[v0]/cs0];
timescaleDimensionCheck = TrueQ[clean[-omegaDim - tBgDim] == dimless];
flowSpeedDimensionCheck = TrueQ[clean[v0Dim - csDim] == dimless];
validityTimescale = <|
  "backgroundEvolutionTime" -> tBg,
  "definition" -> HoldForm[tBg == Abs[backgroundQuantity/backgroundRate]],
  "waveTime" -> HoldForm[1/omega],
  "smallParameter" -> timescaleParameter,
  "condition" -> HoldForm[LessLess[1/(omega tBg), 1]],
  "dimensionless" -> timescaleDimensionCheck,
  "notInferredFromV0" -> TrueQ[FreeQ[timescaleParameter, v0]]
|>;
validityFlowSpeed = <|
  "MachParameter" -> flowSpeedParameter,
  "condition" -> HoldForm[LessLess[Abs[v0]/cs0, 1]],
  "dimensionless" -> flowSpeedDimensionCheck,
  "independentOfTimescaleCondition" -> TrueQ[FreeQ[flowSpeedParameter, omega | tBg]]
|>;
convectedOperator = Expand[(dtOp + v0 dwOp)^2];
leadingConvectivePower = Exponent[convectedOperator - dtOp^2, v0, Min];
leadingConvectiveTerm = clean[
  Coefficient[convectedOperator - dtOp^2, v0, leadingConvectivePower] v0^leadingConvectivePower];
spatialToTemporalScalingRule = dwOp -> dtOp/cs0;
leadingRelativeCorrection = clean[
  (leadingConvectiveTerm/dtOp^2) /. spatialToTemporalScalingRule];
relativeOrderFromCorrection[relativeCorrection_] := If[
  TrueQ[Exponent[Numerator[Together[relativeCorrection]], v0, Min] == 1] &&
    TrueQ[clean[relativeCorrection/(v0/cs0)] =!= 0],
  HoldForm[O[v0/cs0]], "NOT_ESTABLISHED"];
discardedRelativeOrder = relativeOrderFromCorrection[leadingRelativeCorrection];
discardedConvectiveOrder = <|
  "expandedConvectedTimeOperator" -> convectedOperator,
  "leadingPowerOfV0" -> leadingConvectivePower,
  "spatialToTemporalScaling" -> HoldForm[dwOp == dtOp/cs0],
  "leadingRelativeCorrection" -> leadingRelativeCorrection,
  "relativeOrder" -> discardedRelativeOrder,
  "alsoDiscardedInMassFlux" -> HoldForm[deltaRho v0]
|>;

(* Each conclusion-bearing emitted object is itself on the checked data path. *)
projectionFiniteCheck = And[
  TrueQ[projectionFinite["source"] === sourceFinite],
  TrueQ[projectionFinite["rho3"] === int[om[w] rho[w], {w, w1, w2}]],
  TrueQ[projectionFinite["j3"] === int[om[w] jX[w], {w, w1, w2}]],
  TrueQ[projectionFinite["projectedBalance"] ===
    (dt[rho3] + divX[j3] == projectionFinite["source"])],
  TrueQ[productRuleResidual == 0]
];
projectionInfiniteCheck = And[
  TrueQ[projectionInfinite["endpointCondition"] === infiniteEndpointCondition],
  TrueQ[projectionInfinite["source"] === sourceInfinite]
];
parityEvenOutputCheck = And[
  TrueQ[parityEven["finiteSymmetricSource"] === sourceEvenSymmetric],
  TrueQ[parityEven["infiniteSymmetricSource"] === sourceEvenInfinite],
  TrueQ[parityEven["exact"] === evenParityCheck],
  TrueQ[clean[parityEven["finiteSymmetricSource"]] == 0],
  TrueQ[clean[parityEven["infiniteSymmetricSource"]] == 0]
];
parityOddOutputCheck = And[
  TrueQ[parityOdd["finiteSymmetricSource"] === sourceOddSymmetric],
  TrueQ[parityOdd["infiniteSymmetricSource"] === sourceOddInfinite],
  TrueQ[parityOdd["exact"] === oddParityCheck],
  TrueQ[parityOdd["notForcedToVanish"] === TrueQ[oddParityReduced =!= 0]]
];
dynamicWindowExtrasCheck = And[
  TrueQ[dynamicWindowExtras["timeProductRuleTerm"] === dynamicExtraTime],
  TrueQ[dynamicWindowExtras["tangentialProductRuleTerm"] === dynamicExtraTangential],
  TrueQ[dynamicWindowExtras["numberOfTermsAbsentFromStaticA1"] ===
    Length[{dynamicWindowExtras["timeProductRuleTerm"],
      dynamicWindowExtras["tangentialProductRuleTerm"]}]],
  TrueQ[dynamicWindowExtras["signsInProjectedBalance"] === dynamicExtraSigns],
  TrueQ[dynamicExtraSigns === {1, 1}]
];
zByRegimeCheck = And[
  TrueQ[zByRegime["qSquaredPositive"]["q"] === qR],
  TrueQ[clean[zByRegime["qSquaredPositive"]["Z"] -
    (zImpermeable /. q -> zByRegime["qSquaredPositive"]["q"])] == 0],
  TrueQ[zByRegime["qSquaredPositive"]["ReIm"] === clean[ComplexExpand[{
    Re[zByRegime["qSquaredPositive"]["Z"]], Im[zByRegime["qSquaredPositive"]["Z"]]}]]],
  TrueQ[zByRegime["qSquaredNegative"]["q"] === I alpha],
  TrueQ[clean[zByRegime["qSquaredNegative"]["Z"] -
    (zImpermeable /. q -> zByRegime["qSquaredNegative"]["q"])] == 0],
  TrueQ[zByRegime["qSquaredNegative"]["ReIm"] === clean[ComplexExpand[{
    Re[zByRegime["qSquaredNegative"]["Z"]], Im[zByRegime["qSquaredNegative"]["Z"]]}]]],
  TrueQ[zByRegime["qSquaredZero"] === grazingBehavior]
];
zByParityCheck = And[
  TrueQ[zByParity["amplitudesAfterRadiationCondition"] === radiatingAmplitudes],
  TrueQ[zByParity["amplitudeCountAfterRadiationCondition"] ===
    Length[zByParity["amplitudesAfterRadiationCondition"]]],
  TrueQ[zByParity["faceConditionCount"] === faceConditionCount],
  TrueQ[zByParity["faceConditionsFixAllAmplitudes"] === allRadiatingAmplitudesFixed],
  TrueQ[zByParity["deltaW"]["motionGlobalW"] === zetaThickness],
  TrueQ[zByParity["deltaW"]["outwardFaceVelocities"] === outwardVelocity[zetaThickness]],
  TrueQ[zByParity["deltaW"]["amplitudeSolution"] === amplitudeSolution[zetaThickness]],
  TrueQ[zByParity["deltaW"]["facePressures"] === solvedFacePressures[zetaThickness]],
  TrueQ[zByParity["deltaW"]["ZPerFace"] === solvedPerFaceResponse[zetaThickness]],
  TrueQ[zByParity["deltaW"]["byRegime"] === parityRegimeResponse],
  TrueQ[zByParity["zetaC"]["motionGlobalW"] === zetaCenter],
  TrueQ[zByParity["zetaC"]["outwardFaceVelocities"] === outwardVelocity[zetaCenter]],
  TrueQ[zByParity["zetaC"]["amplitudeSolution"] === amplitudeSolution[zetaCenter]],
  TrueQ[zByParity["zetaC"]["facePressures"] === solvedFacePressures[zetaCenter]],
  TrueQ[zByParity["zetaC"]["ZPerFace"] === solvedPerFaceResponse[zetaCenter]],
  TrueQ[zByParity["zetaC"]["byRegime"] === parityRegimeResponse],
  TrueQ[zByParity["perFaceResponseDependsOnCombination"] ===
    perFaceResponseDependsOnCombination],
  TrueQ[zByParity["nonCombinationAblation"]["motionGlobalW"] === nonCombinationMotion],
  TrueQ[zByParity["nonCombinationAblation"]["solvedResponse"] ===
    nonCombinationModeResponse],
  TrueQ[zByParity["nonCombinationAblation"]["changesSolvedAmplitudesAndPressures"] ===
    nonCombinationChangesSolvedQuantities],
  TrueQ[clean[nonCombinationMotion[[1]] - nonCombinationMotion[[2]]] =!= 0],
  TrueQ[clean[nonCombinationMotion[[1]] + nonCombinationMotion[[2]]] =!= 0]
];
grazingBehaviorCheck = And[
  TrueQ[grazingBehavior["q"] === 0],
  TrueQ[grazingBehavior["finiteAmplitudeOutwardVelocity"] === clean[(I q amp) /. q -> 0]],
  TrueQ[grazingBehavior["drivenNonzeroFaceVelocity"] === grazingDrivenStatus],
  TrueQ[grazingBehavior["undrivenFace"] === grazingUndrivenStatus],
  TrueQ[grazingBehavior["ZPoleByReciprocalTest"] === grazingPoleTest],
  TrueQ[grazingBehavior["propagatingApproach"]["ReZPoleByReciprocalTest"] ===
    propagatingRePoleTest],
  TrueQ[grazingBehavior["propagatingApproach"]["ImZ"] === zPropParts[[2]]],
  TrueQ[grazingBehavior["evanescentApproach"]["ReZ"] === zEvanParts[[1]]],
  TrueQ[grazingBehavior["evanescentApproach"]["ImZPoleByReciprocalTest"] ===
    evanescentImPoleTest],
  TrueQ[grazingBehavior["addedMassPerFaceReciprocalsAtGrazing"] ===
    clean[(1/addedMassExpected) /. alpha -> 0]],
  TrueQ[grazingBehavior["addedMassSignsUpperLower"] === addedMassSigns],
  TrueQ[grazingBehavior["addedMassPoleByReciprocalTest"] === grazingAddedMassPoleTest]
];
relativeFluxChannelsCheck = Module[{emittedChannels, emittedMatrix},
  emittedChannels = {relativeFluxChannels["netAccretionBySlab"],
    relativeFluxChannels["throughFlowOrientedLowerToUpper"]};
  emittedMatrix = Table[Coefficient[emittedChannels[[i]], faceFluxes[[j]]],
    {i, Length[emittedChannels]}, {j, Length[faceFluxes]}];
  And[
    TrueQ[relativeFluxChannels["JplusDefinition"] === jPlusDefinition],
    TrueQ[relativeFluxChannels["JminusDefinition"] === jMinusDefinition],
    TrueQ[emittedMatrix === fluxOrientationMatrix],
    TrueQ[relativeFluxChannels["inverse"] ===
      Thread[faceFluxes -> clean[Inverse[emittedMatrix] . {acc, thru}]]]
  ]
];
windowParityControlCheck = And[
  TrueQ[windowParityControl["family"] === windowFamily],
  TrueQ[windowParityControl["interval"] === windowControlInterval],
  TrueQ[windowParityControl["evenCurrentSourceAsFunctionOfB"] === windowControlEven],
  TrueQ[windowParityControl["oddCurrentSourceAsFunctionOfB"] === windowControlOdd],
  TrueQ[windowParityControl["evenCurrentIdenticallyIndependentOfB"] === Not[evenBChanged]],
  TrueQ[windowParityControl["oddCurrentIdenticallyIndependentOfB"] === Not[oddBChanged]],
  TrueQ[windowParityControl["b0EvenSource"] === clean[
    windowParityControl["evenCurrentSourceAsFunctionOfB"] /. b -> 0]],
  TrueQ[windowParityControl["bControlDomain"] === windowBControlDomain]
];
intervalSymmetryControlCheck = And[
  TrueQ[intervalSymmetryControl["family"] === intervalControlFamily],
  TrueQ[intervalSymmetryControl["window"] === intervalControlWindow],
  TrueQ[intervalSymmetryControl["evenCurrentSourceAsFunctionOfC"] === intervalControlEven],
  TrueQ[intervalSymmetryControl["oddCurrentSourceAsFunctionOfC"] === intervalControlOdd],
  TrueQ[intervalSymmetryControl["evenCurrentIdenticallyIndependentOfC"] === Not[evenCChanged]],
  TrueQ[intervalSymmetryControl["oddCurrentIdenticallyIndependentOfC"] === Not[oddCChanged]],
  TrueQ[intervalSymmetryControl["c0EvenSource"] === clean[
    intervalSymmetryControl["evenCurrentSourceAsFunctionOfC"] /. c -> 0]],
  TrueQ[intervalSymmetryControl["c0OddSource"] === clean[
    intervalSymmetryControl["oddCurrentSourceAsFunctionOfC"] /. c -> 0]],
  TrueQ[intervalSymmetryControl["controlDomain"] === intervalControlDomain]
];
discardedConvectiveOrderCheck = Module[
  {emittedOperator, emittedPower, emittedLeadingTerm, emittedRelativeCorrection},
  emittedOperator = discardedConvectiveOrder["expandedConvectedTimeOperator"];
  emittedPower = Exponent[emittedOperator - dtOp^2, v0, Min];
  emittedLeadingTerm = clean[
    Coefficient[emittedOperator - dtOp^2, v0, emittedPower] v0^emittedPower];
  emittedRelativeCorrection = clean[
    (emittedLeadingTerm/dtOp^2) /. spatialToTemporalScalingRule];
  And[
    TrueQ[discardedConvectiveOrder["leadingPowerOfV0"] === emittedPower],
    TrueQ[discardedConvectiveOrder["spatialToTemporalScaling"] ===
      HoldForm[dwOp == dtOp/cs0]],
    TrueQ[discardedConvectiveOrder["leadingRelativeCorrection"] ===
      emittedRelativeCorrection],
    TrueQ[discardedConvectiveOrder["relativeOrder"] ===
      relativeOrderFromCorrection[discardedConvectiveOrder["leadingRelativeCorrection"]]],
    TrueQ[discardedConvectiveOrder["alsoDiscardedInMassFlux"] === HoldForm[deltaRho v0]]
  ]
];
addedMassOutputCheck = And[
  TrueQ[addedMassOutput["purelyImaginaryRegime"] === HoldForm[q^2 < 0]],
  TrueQ[addedMassOutput["perFaceUpperLowerGlobalWDefinition"] === addedMassByFace],
  TrueQ[addedMassOutput["deltaWPerFace"] ===
    addedMassOutput["perFaceUpperLowerGlobalWDefinition"]],
  TrueQ[addedMassOutput["zetaCPerFace"] ===
    addedMassOutput["perFaceUpperLowerGlobalWDefinition"]],
  TrueQ[addedMassOutput["signsUpperLower"] === FullSimplify[
    Sign /@ addedMassOutput["perFaceUpperLowerGlobalWDefinition"], assumptions]],
  TrueQ[addedMassOutput["check"] === TrueQ[clean[
    addedMassOutput["perFaceUpperLowerGlobalWDefinition"] - addedMassExpected] == {0, 0}]]
];
degenerateLociCheck = And[
  TrueQ[degenerateLoci["qSquaredPositive"]["locus"] === propDegenerate],
  TrueQ[degenerateLoci["qSquaredPositive"]["drivenNoSolutionLocus"] ===
    propDrivenNoSolution],
  TrueQ[degenerateLoci["qSquaredPositive"]["drivenFreeAmplitudeLocus"] === propDrivenFree],
  TrueQ[degenerateLoci["qSquaredPositive"]["undriven"] ===
    undrivenFreeAmplitudeStatus[bulkAmplitudeCoefficient /. q -> omega r,
      degenerateLoci["qSquaredPositive"]["locus"],
      Element[{r, lambdaP0, x}, Reals] && r > 0 && x >= 0]],
  TrueQ[degenerateLoci["qSquaredNegative"]["locus"] === evanDegenerate],
  TrueQ[degenerateLoci["qSquaredZero"]["locus"] === grazingDegenerate],
  TrueQ[degenerateLoci["qSquaredZero"]["drivenNoSolutionLocus"] ===
    grazingDrivenNoSolution],
  TrueQ[degenerateLoci["qSquaredZero"]["drivenFreeAmplitudeLocus"] === grazingDrivenFree],
  TrueQ[degenerateLoci["qSquaredZero"]["undriven"] ===
    undrivenFreeAmplitudeStatus[bulkAmplitudeCoefficient /. q -> 0,
      degenerateLoci["qSquaredZero"]["locus"],
      Element[{lambdaP0, x}, Reals] && x >= 0]]
];
dimensionRouteKindChecks = KeyValueMap[
  TrueQ[#2 === routeKindFromAssertion[dimensionRoutes[#1], dimensionRouteTargets[#1]]] &,
  routeKinds];

(* Cross-check exact real/imaginary reconstructions and limiting relations. *)
responseChecks = {
  TrueQ[clean[zImpermeable - rhoM omega/q] == 0],
  TrueQ[solvedPerFaceResponse[zetaThickness] == {zImpermeable, zImpermeable}],
  TrueQ[solvedPerFaceResponse[zetaCenter] == {zImpermeable, zImpermeable}],
  addedMassCheck,
  grazingPoleTest,
  zPermeableCheck,
  TrueQ[clean[reProp + I imProp - zPermProp] == 0],
  TrueQ[clean[reEvan + I imEvan - zPermEvan] == 0],
  TrueQ[clean[reGrazing + I imGrazing - zPermGrazing] == 0],
  TrueQ[clean[smallX[zPermProp] - ((rhoM + lambdaV0)/(r - lambdaP0))] == 0],
  TrueQ[clean[largeX[zPermProp] - rhoM/r] == 0],
  TrueQ[clean[largeX[zPermEvan] + I rhoM/s] == 0]
};
controlChecks = {
  TrueQ[clean[windowControlEven /. b -> 0] == 0],
  evenBChanged, oddBChanged, evenCChanged, oddCChanged
};
gatedOutputChecks = {
  projectionFiniteCheck, projectionInfiniteCheck, parityEvenOutputCheck,
  parityOddOutputCheck, dynamicWindowExtrasCheck, zByRegimeCheck, zByParityCheck,
  grazingBehaviorCheck, relativeFluxChannelsCheck, windowParityControlCheck,
  intervalSymmetryControlCheck, discardedConvectiveOrderCheck,
  addedMassOutputCheck, degenerateLociCheck
};
allChecks = Join[
  {TrueQ[productRuleResidual == 0], evenParityCheck, oddParityCheck,
    And @@ propagatingBranchCheck, And @@ evanescentBranchCheck, fluxChannelCheck,
    closureScopeStructuralCheck, timescaleDimensionCheck, flowSpeedDimensionCheck},
  responseChecks, dimensionChecks, controlChecks, gatedOutputChecks,
  dimensionRouteKindChecks
];
internalConsistency = And @@ allChecks;

(* Required tagged output. *)
emit["S11BA_PROJECTION_FINITE", projectionFinite];
emit["S11BA_PROJECTION_INFINITE", projectionInfinite];
emit["S11BA_PARITY_EVEN_JW", parityEven];
emit["S11BA_PARITY_ODD_JW", parityOdd];
emit["S11BA_PARITY_INTERVAL", parityInterval];
emit["S11BA_DYNAMIC_WINDOW_EXTRA_TERMS", dynamicWindowExtras];
emit["S11BA_Z_IMPERMEABLE", zImpermeableOutput];
emit["S11BA_Z_BY_REGIME", zByRegime];
emit["S11BA_Z_BY_PARITY", zByParity];
emit["S11BA_ADDED_MASS", addedMassOutput];
emit["S11BA_GRAZING_BEHAVIOUR", grazingBehavior];
emit["S11BA_RELATIVE_FLUX_CHANNELS", relativeFluxChannels];
emit["S11BA_PERMEABLE_COEFF_DIMS", <|"Lambda_p0" -> lambdaPDim, "Lambda_V0" -> lambdaVDim,
  "tau" -> tauDim, "basisOrder" -> {"L", "T", "M"}|>];
emit["S11BA_Z_PERMEABLE", <|"x" -> HoldForm[omega tau], "genericZ" -> zPermeable,
  "deltaW_ZPerFace" -> {zPermeable, zPermeable}, "zetaC_ZPerFace" -> {zPermeable, zPermeable}|>];
emit["S11BA_DISSIPATIVE_BY_REGIME_AND_PARITY", inPhaseByRegime];
emit["S11BA_DISSIPATION_VS_OMEGA_TAU", dissipationVsX];
emit["S11BA_TAU_ZERO_LIMIT", tauZeroLimit];
emit["S11BA_CLOSURE_SCOPE_LIMIT", closureScope];
emit["S11BA_DEGENERATE_LOCI", degenerateLoci];

KeyValueMap[(emit["S11BA_DIM_" <> #1, <|"exponentsLTM" -> #2, "route" -> dimensionRoutes[#1]|>]) &,
  dimensionValues];
KeyValueMap[(emit["S11BA_DIM_ROUTE_KIND_" <> #1, routeKinds[#1]]) &, dimensionValues];

emit["S11BA_CONTROL_WINDOW_PARITY", windowParityControl];
emit["S11BA_CONTROL_INTERVAL_SYMMETRY", intervalSymmetryControl];
emit["S11BA_VALIDITY_TIMESCALE", validityTimescale];
emit["S11BA_VALIDITY_FLOW_SPEED", validityFlowSpeed];
emit["S11BA_DISCARDED_CONVECTIVE_ORDER", discardedConvectiveOrder];

ablationCoverageTags = {
  "S11BA_PROJECTION_FINITE", "S11BA_PROJECTION_INFINITE",
  "S11BA_PARITY_EVEN_JW", "S11BA_PARITY_ODD_JW",
  "S11BA_DYNAMIC_WINDOW_EXTRA_TERMS", "S11BA_Z_BY_REGIME",
  "S11BA_Z_BY_PARITY", "S11BA_GRAZING_BEHAVIOUR",
  "S11BA_RELATIVE_FLUX_CHANNELS", "S11BA_CONTROL_WINDOW_PARITY",
  "S11BA_CONTROL_INTERVAL_SYMMETRY", "S11BA_DISCARDED_CONVECTIVE_ORDER",
  "S11BA_ADDED_MASS", "S11BA_DEGENERATE_LOCI",
  "S11BA_DIM_ROUTE_KIND_Z", "S11BA_DIM_ROUTE_KIND_m_add",
  "S11BA_DIM_ROUTE_KIND_rho_m", "S11BA_DIM_ROUTE_KIND_c_s0",
  "S11BA_DIM_ROUTE_KIND_v0", "S11BA_DIM_ROUTE_KIND_Lambda_p0",
  "S11BA_DIM_ROUTE_KIND_Lambda_V0", "S11BA_DIM_ROUTE_KIND_tau",
  "S11BA_DIM_ROUTE_KIND_A1_source"
};
emit["S11BA_ABLATION_COVERAGE", ablationCoverageTags];

Print["VERDICT: " <> If[TrueQ[internalConsistency], "PASS", "INTERNAL_CONTRADICTION"]];
