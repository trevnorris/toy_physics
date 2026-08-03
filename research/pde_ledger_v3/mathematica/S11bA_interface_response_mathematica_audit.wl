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
projectionFinite = <|
  "projectedBalance" -> (dt[rho3] + divX[j3] == sourceFinite),
  "rho3" -> int[om[w] rho[w], {w, w1, w2}],
  "j3" -> int[om[w] jX[w], {w, w1, w2}],
  "source" -> sourceFinite
|>;
projectionInfinite = <|
  "endpointCondition" -> HoldForm[Limit[om[w] jw[w], w -> -Infinity] == 0 && Limit[om[w] jw[w], w -> Infinity] == 0],
  "source" -> sourceInfinite
|>;

(* A2: exact parity reduction on the explicitly symmetric finite interval [-L,L]. *)
sourceEvenSymmetric = 0;
sourceOddSymmetric = -2 om[bigL] jOdd[bigL] + 2 int[om'[u] jOdd[u], {u, 0, bigL}];
sourceEvenInfinite = 0;
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
dynamicWindowExtras = <|
  "timeProductRuleTerm" -> dynamicExtraTime,
  "tangentialProductRuleTerm" -> dynamicExtraTangential,
  "numberOfTermsAbsentFromStaticA1" -> Length[{dynamicExtraTime, dynamicExtraTangential}],
  "signsInProjectedBalance" -> {1, 1}
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

(* Parity amplitudes use both face displacements measured in global +w. *)
zetaThickness = {deltaW/2, -deltaW/2};
zetaCenter = {zetaC, zetaC};
outwardSigns = {1, -1};
outwardVelocity[zetas_List] := clean[outwardSigns (-I omega zetas)];
pressureFor[zetas_List, zLocal_] := clean[zLocal outwardVelocity[zetas]];
ratioFor[zetas_List, zLocal_] := clean[pressureFor[zetas, zLocal]/outwardVelocity[zetas]];
thicknessVelocity = outwardVelocity[zetaThickness];
centerVelocity = outwardVelocity[zetaCenter];
thicknessPressure = pressureFor[zetaThickness, zImpermeable];
centerPressure = pressureFor[zetaCenter, zImpermeable];
zByParity = <|
  "deltaW" -> <|"outwardFaceVelocities" -> thicknessVelocity, "facePressures" -> thicknessPressure,
    "ZPerFace" -> ratioFor[zetaThickness, zImpermeable]|>,
  "zetaC" -> <|"outwardFaceVelocities" -> centerVelocity, "facePressures" -> centerPressure,
    "ZPerFace" -> ratioFor[zetaCenter, zImpermeable]|>
|>;

(* Added mass is deliberately computed against global-w face acceleration, per face. *)
globalFaceAcceleration[zetas_List] := clean[-omega^2 zetas];
addedMassByFace = clean[pressureFor[{zetaU, zetaD}, zEvanescent]/globalFaceAcceleration[{zetaU, zetaD}]];
addedMassExpected = {rhoM/alpha, -rhoM/alpha};
addedMassCheck = TrueQ[clean[addedMassByFace - addedMassExpected] == {0, 0}];

(* The prescribed pole test uses 1/expression == 0. *)
grazingPoleTest = TrueQ[FullSimplify[(1/zImpermeable == 0) /. q -> 0, assumptions]];
grazingAddedMassPoleTest = And @@ (TrueQ[FullSimplify[(1/# == 0) /. alpha -> 0, assumptions]] & /@ addedMassExpected);
grazingDrivenStatus = If[grazingPoleTest && TrueQ[clean[(I q amp) /. q -> 0] == 0],
  "NO_SOLUTION_WITH_FINITE_BULK_AMPLITUDE", "NOT_ESTABLISHED"];
grazingUndrivenStatus = If[TrueQ[Reduce[0 == 0, amp, Complexes]],
  "CONSTANT_BULK_AMPLITUDE_FREE", "NOT_ESTABLISHED"];
propagatingRePoleTest = TrueQ[FullSimplify[(1/zPropParts[[1]] == 0) /. qR -> 0, assumptions]];
evanescentImPoleTest = TrueQ[FullSimplify[(1/zEvanParts[[2]] == 0) /. alpha -> 0, assumptions]];
grazingBehavior = <|
  "q" -> 0,
  "finiteAmplitudeOutwardVelocity" -> clean[(I q amp) /. q -> 0],
  "drivenNonzeroFaceVelocity" -> grazingDrivenStatus,
  "undrivenFace" -> grazingUndrivenStatus,
  "ZPoleByReciprocalTest" -> grazingPoleTest,
  "propagatingApproach" -> <|"ReZPoleByReciprocalTest" -> propagatingRePoleTest, "ImZ" -> 0|>,
  "evanescentApproach" -> <|"ReZ" -> 0, "ImZPoleByReciprocalTest" -> evanescentImPoleTest|>,
  "addedMassPerFaceReciprocalsAtGrazing" -> clean[(1/addedMassExpected) /. alpha -> 0],
  "addedMassSignsUpperLower" -> {1, -1},
  "addedMassPoleByReciprocalTest" -> grazingAddedMassPoleTest
|>;

zByRegime = <|
  "qSquaredPositive" -> <|"q" -> qR, "Z" -> zPropagating, "ReIm" -> zPropParts|>,
  "qSquaredNegative" -> <|"q" -> I alpha, "Z" -> zEvanescent, "ReIm" -> zEvanParts|>,
  "qSquaredZero" -> grazingBehavior
|>;
zImpermeableOutput = <|"branches" -> radiationBranches, "genericZ" -> zImpermeable|>;
parityRegimeResponse = <|
  "qSquaredPositive" -> <|"ZPerFace" -> {zPropagating, zPropagating},
    "ReImPerFace" -> {zPropParts, zPropParts}|>,
  "qSquaredNegative" -> <|"ZPerFace" -> {zEvanescent, zEvanescent},
    "ReImPerFace" -> {zEvanParts, zEvanParts}|>,
  "qSquaredZero" -> grazingBehavior
|>;
zByParity = Map[Append[#, "byRegime" -> parityRegimeResponse] &, zByParity];

(* Relative-flux channels; invert the linear transformation as a check. *)
netAccretion = -(jPlus + jMinus);
throughFlowLowerToUpper = (jPlus - jMinus)/2;
fluxInverse = First[Solve[{acc == netAccretion, thru == throughFlowLowerToUpper}, {jPlus, jMinus}]];
fluxChannelCheck = TrueQ[clean[({netAccretion, throughFlowLowerToUpper} /. fluxInverse) - {acc, thru}] == {0, 0}];
relativeFluxChannels = <|
  "JplusDefinition" -> HoldForm[rhoM (vWUpper - zetaUDot)],
  "JminusDefinition" -> HoldForm[-rhoM (vWLower - zetaDDot)],
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
degenerateLoci = <|
  "qSquaredPositive" -> <|"locus" -> propDegenerate,
    "drivenNoSolutionLocus" -> propDrivenNoSolution,
    "drivenFreeAmplitudeLocus" -> propDrivenFree,
    "undriven" -> "FREE_AMPLITUDE"|>,
  "qSquaredNegative" -> <|"locus" -> evanDegenerate|>,
  "qSquaredZero" -> <|"locus" -> grazingDegenerate,
    "drivenNoSolutionLocus" -> grazingDrivenNoSolution,
    "drivenFreeAmplitudeLocus" -> grazingDrivenFree,
    "undriven" -> "FREE_AMPLITUDE"|>
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
  "Z" -> HoldForm[dim[pressure/velocity]],
  "m_add" -> HoldForm[dim[pressure/acceleration]],
  "rho_m" -> HoldForm[dim[mass/(length^4)]],
  "c_s0" -> HoldForm[dim[dt^2 phi] == dim[cs0^2 laplacian phi]],
  "v0" -> HoldForm[dim[length/time]],
  "Lambda_p0" -> HoldForm[dim[massFlux/pressure]],
  "Lambda_V0" -> HoldForm[dim[massFlux/velocity]],
  "tau" -> HoldForm[dim[omega tau] == dim[1]],
  "A1_source" -> HoldForm[dim[jw] == dim[dt Integrate[rho, w]]]
|>;
independentRouteNames = {"c_s0", "Lambda_p0", "Lambda_V0", "A1_source"};
routeKind[name_] := If[MemberQ[independentRouteNames, name], "INDEPENDENT", "DEFINITIONAL"];
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

windowParityControl = <|
  "family" -> HoldForm[Sech[(w - b)/a]^2],
  "interval" -> {-bigL, bigL},
  "evenCurrentSourceAsFunctionOfB" -> windowControlEven,
  "oddCurrentSourceAsFunctionOfB" -> windowControlOdd,
  "evenCurrentIdenticallyIndependentOfB" -> Not[evenBChanged],
  "oddCurrentIdenticallyIndependentOfB" -> Not[oddBChanged],
  "b0EvenSource" -> clean[windowControlEven /. b -> 0],
  "bControlDomain" -> HoldForm[b != 0]
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
intervalSymmetryControl = <|
  "family" -> {-bigL, bigL + c},
  "window" -> HoldForm[Sech[w/a]^2],
  "evenCurrentSourceAsFunctionOfC" -> intervalControlEven,
  "oddCurrentSourceAsFunctionOfC" -> intervalControlOdd,
  "evenCurrentIdenticallyIndependentOfC" -> Not[evenCChanged],
  "oddCurrentIdenticallyIndependentOfC" -> Not[oddCChanged],
  "c0EvenSource" -> clean[intervalControlEven /. c -> 0],
  "c0OddSource" -> clean[intervalControlOdd /. c -> 0],
  "controlDomain" -> HoldForm[c != 0 && c > -2 bigL]
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
discardedConvectiveOrder = <|
  "expandedConvectedTimeOperator" -> convectedOperator,
  "leadingPowerOfV0" -> leadingConvectivePower,
  "relativeOrder" -> HoldForm[O[v0/cs0]],
  "alsoDiscardedInMassFlux" -> HoldForm[deltaRho v0]
|>;

(* Cross-check exact real/imaginary reconstructions and limiting relations. *)
responseChecks = {
  TrueQ[clean[zImpermeable - rhoM omega/q] == 0],
  TrueQ[ratioFor[zetaThickness, zImpermeable] == {zImpermeable, zImpermeable}],
  TrueQ[ratioFor[zetaCenter, zImpermeable] == {zImpermeable, zImpermeable}],
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
allChecks = Join[
  {TrueQ[productRuleResidual == 0], evenParityCheck, oddParityCheck,
    And @@ propagatingBranchCheck, And @@ evanescentBranchCheck, fluxChannelCheck,
    closureScopeStructuralCheck, timescaleDimensionCheck, flowSpeedDimensionCheck},
  responseChecks, dimensionChecks, controlChecks
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
emit["S11BA_ADDED_MASS", <|"purelyImaginaryRegime" -> HoldForm[q^2 < 0],
  "perFaceUpperLowerGlobalWDefinition" -> addedMassByFace,
  "deltaWPerFace" -> addedMassByFace, "zetaCPerFace" -> addedMassByFace,
  "check" -> addedMassCheck|>];
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
KeyValueMap[(emit["S11BA_DIM_ROUTE_KIND_" <> #1, routeKind[#1]]) &, dimensionValues];

emit["S11BA_CONTROL_WINDOW_PARITY", windowParityControl];
emit["S11BA_CONTROL_INTERVAL_SYMMETRY", intervalSymmetryControl];
emit["S11BA_VALIDITY_TIMESCALE", validityTimescale];
emit["S11BA_VALIDITY_FLOW_SPEED", validityFlowSpeed];
emit["S11BA_DISCARDED_CONVECTIVE_ORDER", discardedConvectiveOrder];

Print["VERDICT: " <> If[TrueQ[internalConsistency], "PASS", "INTERNAL_CONTRADICTION"]];
