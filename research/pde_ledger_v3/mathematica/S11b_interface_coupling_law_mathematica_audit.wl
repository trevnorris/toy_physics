$HistoryLength = 0;
ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};

(* The emitted name is constructed only here.  Internal keys are deliberately
   not strings in the public tag grammar. *)
sharedObject[name_String] := HoldComplete[SharedS11bObject, name];
localObject[name_String] := HoldComplete[LocalS11bObject, name];

standardEmissionName[HoldComplete[SharedS11bObject, name_String]] :=
  "WL_S11B_" <> name;
standardEmissionName[HoldComplete[LocalS11bObject, name_String]] :=
  "WL_S11B_LOCAL_" <> name;

emittedNames = <||>;
localNames = {};
emit[key_, payload_] := Module[{name, rendered, stream},
  name = standardEmissionName[key];
  If[!StringQ[name], Quit[90]];
  If[KeyExistsQ[emittedNames, name], Quit[91]];
  AssociateTo[emittedNames, name -> True];
  If[StringStartsQ[name, "WL_S11B_LOCAL_"], AppendTo[localNames, name]];
  rendered = ToString[payload, InputForm, PageWidth -> Infinity];
  stream = First[$Output];
  WriteString[stream, name <> ": " <> rendered <> "\n"];
  Flush[stream];
];
emitShared[name_String, payload_] := emit[sharedObject[name], payload];
emitLocal[name_String, payload_] := emit[localObject[name], payload];

casTest[object_] := Simplify[object];
comparisonRecord[a_, b_] := Module[{residual},
  residual = Together[a - b];
  <|"OPERAND_A" -> a, "OPERAND_B" -> b, "RESIDUAL" -> residual,
    "TEST_OBJECT" -> casTest[residual == 0]|>
];
coefficientVector[expression_, variables_List] :=
  (Coefficient[expression, #] & /@ variables);
normalizedPolynomial[polynomial_, variables_List] := Module[
  {expanded, rules, pivot},
  expanded = Expand[polynomial];
  rules = CoefficientRules[expanded, variables];
  If[rules === {}, Return[expanded]];
  pivot = Last[First[rules]];
  Cancel[expanded/pivot]
];
rationalIdentityRecord[a_, b_, variable_] := Module[
  {difference, cancelled, numerator},
  difference = Together[a - b];
  cancelled = Cancel[difference];
  numerator = Numerator[Together[cancelled]];
  <|"OPERAND_A" -> a, "OPERAND_B" -> b, "RESIDUAL" -> cancelled,
    "NUMERATOR" -> numerator,
    "TEST_OBJECT" -> casTest[PolynomialQ[numerator, variable] &&
       numerator === 0]|>
];

tokenFromBoolean[object_] := Which[
  TrueQ[object], "PROVED_TRUE",
  TrueQ[object === False], "PROVED_FALSE",
  True, "UNDECIDED"
];

(* ---------------------------------------------------------------------- *)
(* Supplied geometry, ansatz, acoustic equations, energy terms and laws.   *)
(* ---------------------------------------------------------------------- *)

dimBrane = 3;
spatialIndices = Range[dimBrane];

$Assumptions = Element[{w0, rho4, rhoBr, rhoM, cs0, muR, br3, bRho,
      muW, kW, kappaW, crossC, lambdaA0, lambdaV0, lambdaX0,
      tauA, tauV, tauX, k, vDr}, Reals] &&
    w0 > 0 && rho4 > 0 && rhoBr > 0 && rhoM > 0 && cs0 > 0 &&
    muW > 0 && k >= 0 && tauA >= 0 && tauV >= 0 && tauX >= 0;

harmonicFactor = Exp[I (k x - omega t)];
qSquared = Expand[omega^2/cs0^2 - k^2];
qAlgebraicEquation = q^2 == qSquared;
qPrincipalProduct[z_] := Sqrt[z - cs0 k] Sqrt[z + cs0 k]/cs0;
qUpperRim[z_] := qPrincipalProduct[z];
qContinued[z_] := Piecewise[{
   {-qPrincipalProduct[z], Im[z] < 0 && -cs0 k < Re[z] < cs0 k}},
  qPrincipalProduct[z]];

lambdaA = lambdaA0/(1 - I omega tauA);
lambdaV = lambdaV0/(1 - I omega tauV);
lambdaX = lambdaX0/(1 - I omega tauX);

rhoBrDefinition = rhoBr == rho4 w0;
eWDefinition = dW == w0 eW;
faceDisplacements = {zetaPlus == dW/2, zetaMinus == -dW/2};
outwardFaceVelocities = {vPlus == -I omega dW/2,
  vMinus == -I omega dW/2};
acousticEquations = {
  vBulk == Grad[phi[x1, x2, x3, w, t], {x1, x2, x3, w}],
  deltaP == -rhoM D[phi[x1, x2, x3, w, t], t],
  D[phi[x1, x2, x3, w, t], {t, 2}] ==
    cs0^2 Laplacian[phi[x1, x2, x3, w, t], {x1, x2, x3, w}]
};

(* Component variables used by the invariant constructor. *)
gVariables = Flatten[Table[Symbol["g" <> ToString[i] <> ToString[j]],
    {i, spatialIndices}, {j, spatialIndices}]];
gMatrix = Partition[gVariables, dimBrane];
thetaGradient = Table[Symbol["gt" <> ToString[i]], {i, spatialIndices}];
eGradient = Table[Symbol["ge" <> ToString[i]], {i, spatialIndices}];
basisVariables = Join[gVariables, {thetaField, eField}, thetaGradient,
  eGradient];
basisMonomials = Flatten[Table[basisVariables[[i]] basisVariables[[j]],
    {i, Length[basisVariables]}, {j, i, Length[basisVariables]}]];

rotationGenerators = {
  {{0, -1, 0}, {1, 0, 0}, {0, 0, 0}},
  {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},
  {{0, 0, 0}, {0, 0, -1}, {0, 1, 0}}
};

generatorVariation[generator_] := Module[{deltaG, deltaThetaGradient,
    deltaEGradient},
  deltaG = generator.gMatrix - gMatrix.generator;
  deltaThetaGradient = generator.thetaGradient;
  deltaEGradient = generator.eGradient;
  AssociationThread[basisVariables,
    Join[Flatten[deltaG], {0, 0}, deltaThetaGradient, deltaEGradient]]
];

lieDerivative[polynomial_, variation_Association] := Expand[Total[
  MapThread[D[polynomial, #1] #2 &,
    {Keys[variation], Values[variation]}]]];

constructInvariantBasis[] := Module[
  {variations, actionRows, inversionRows, constraintMatrix, nullVectors,
   invariantPolynomials, fourierRules, fourierImages, imageMonomials,
   imageMatrix, independentIndices, trialRank, currentRank, supplied,
   suppliedImages, spanRows, missing, candidate, candidateImage},

  variations = generatorVariation /@ rotationGenerators;
  actionRows = Flatten[Table[
    With[{actions = lieDerivative[#, variations[[generatorIndex]]] & /@
        basisMonomials},
      Table[Coefficient[actions[[column]], basisMonomials[[row]]],
        {row, Length[basisMonomials]}, {column, Length[basisMonomials]}]],
    {generatorIndex, Length[variations]}], 1];
  inversionRows = Table[
    coefficientVector[Expand[(basisMonomials[[column]] /.
        Thread[Join[thetaGradient, eGradient] ->
          -Join[thetaGradient, eGradient]]) - basisMonomials[[column]]],
      basisMonomials], {column, Length[basisMonomials]}];
  inversionRows = Select[inversionRows, Total[Abs[#]] != 0 &];
  constraintMatrix = Join[actionRows, inversionRows];
  nullVectors = NullSpace[constraintMatrix];
  invariantPolynomials = Factor[Dot[#, basisMonomials]] & /@ nullVectors;

  fourierRules = Join[
    Flatten[Table[gMatrix[[i, j]] -> I waveK[i] displacement[j],
      {i, spatialIndices}, {j, spatialIndices}]],
    Table[thetaGradient[[i]] -> I waveK[i] thetaField,
      {i, spatialIndices}],
    Table[eGradient[[i]] -> I waveK[i] eField,
      {i, spatialIndices}]];
  fourierImages = Expand[# /. fourierRules] & /@ invariantPolynomials;
  imageMonomials = DeleteDuplicates[Flatten[
    (First /@ CoefficientRules[#,
          Join[Table[waveK[i], {i, spatialIndices}],
            Table[displacement[i], {i, spatialIndices}],
            {thetaField, eField}]]) & /@ fourierImages, 1]];
  imageMonomials = Times @@@ (Join[Table[waveK[i], {i, spatialIndices}],
           Table[displacement[i], {i, spatialIndices}],
           {thetaField, eField}]^# & /@ imageMonomials);
  imageMatrix = Table[Coefficient[fourierImages[[row]], imageMonomials[[col]]],
    {row, Length[fourierImages]}, {col, Length[imageMonomials]}];
  independentIndices = {};
  currentRank = 0;
  Do[
    trialRank = MatrixRank[imageMatrix[[Append[independentIndices, index]]]];
    If[trialRank > currentRank, AppendTo[independentIndices, index];
      currentRank = trialRank],
    {index, Length[fourierImages]}];
  invariantPolynomials = invariantPolynomials[[independentIndices]];

  supplied = {
    Expand[Total[Flatten[Table[(gMatrix[[i, j]] - gMatrix[[j, i]])^2,
       {i, spatialIndices}, {j, i + 1, dimBrane}]]]],
    thetaField^2, thetaField eField, eField^2,
    Expand[eGradient.eGradient]
  };
  suppliedImages = Expand[# /. fourierRules] & /@ supplied;
  spanRows = suppliedImages;
  currentRank = MatrixRank[Table[coefficientVector[poly, imageMonomials],
      {poly, spanRows}]];
  missing = {};
  Do[
    candidate = invariantPolynomials[[index]];
    candidateImage = Expand[candidate /. fourierRules];
    trialRank = MatrixRank[Table[coefficientVector[poly, imageMonomials],
       {poly, Append[spanRows, candidateImage]}]];
    If[trialRank > currentRank,
      AppendTo[missing, candidate]; AppendTo[spanRows, candidateImage];
      currentRank = trialRank],
    {index, Length[invariantPolynomials]}];

  <|"VARIABLES" -> basisVariables, "MONOMIALS" -> basisMonomials,
    "ROTATION_GENERATORS" -> rotationGenerators,
    "INVARIANTS_BEFORE_DIVERGENCE_QUOTIENT" ->
      (Factor[Dot[#, basisMonomials]] & /@ nullVectors),
    "FOURIER_SUBSTITUTION" -> fourierRules,
    "INDEPENDENT_INVARIANTS" -> Join[supplied, missing],
    "SUPPLIED_INVARIANTS" -> supplied, "OMISSIONS" -> missing|>
];

energyBasisData = constructInvariantBasis[];
energyBasis = energyBasisData["INDEPENDENT_INVARIANTS"];
suppliedBasis = energyBasisData["SUPPLIED_INVARIANTS"];
omittedBasis = energyBasisData["OMISSIONS"];
extraCoefficients = Table[Symbol["chi" <> ToString[index]],
  {index, Length[omittedBasis]}];

curlSquared = suppliedBasis[[1]];
gradientESquared = suppliedBasis[[5]];
suppliedEnergyDensity = Expand[
  muR curlSquared/2 + br3 thetaField^2/2 +
  crossC w0 thetaField eField + kW w0^2 eField^2/2 +
  kappaW w0^4 gradientESquared/2];
fullEnergyDensity = Expand[suppliedEnergyDensity +
  Total[MapThread[#1 #2 &, {extraCoefficients, omittedBasis}]]];

emitShared["ENERGY_BASIS", <|
  "FIELDS_AND_FIRST_GRADIENTS" -> basisVariables,
  "CONSTRUCTION_MONOMIALS" -> basisMonomials,
  "O3_GENERATORS" -> rotationGenerators,
  "DIVERGENCE_QUOTIENT_MAP" -> energyBasisData["FOURIER_SUBSTITUTION"],
  "BASIS" -> energyBasis|>];
emitShared["ENERGY_BASIS_COUNT", Length[energyBasis]];
emitShared["ENERGY_BASIS_OMISSIONS", omittedBasis];
emitShared["ENERGY_BASIS_INDEPENDENT_TERMS",
  AssociationThread[Join[{muR, br3, crossC, kW, kappaW}, extraCoefficients],
    energyBasis]];

(* One-dimensional representative of an isotropic Fourier mode, k along x1. *)
oneDimensionalRules = Join[
  Flatten[Table[gMatrix[[i, j]] -> Which[
      i == 1 && j == 1, uLX,
      i == 1 && j == 2, uTX,
      True, 0], {i, spatialIndices}, {j, spatialIndices}]],
  {thetaField -> thetaAmp, eField -> eAmp},
  Table[thetaGradient[[i]] -> If[i == 1, thetaX, 0],
    {i, spatialIndices}],
  Table[eGradient[[i]] -> If[i == 1, eX, 0],
    {i, spatialIndices}]];
energy1D[sourceEnergy_] := Expand[sourceEnergy /. oneDimensionalRules];

jetVariables = {uL, uLX, uLXX, uT, uTX, uTXX, thetaAmp, thetaX,
  thetaXX, eAmp, eX, eXX};
jetDerivativeRules = {uL -> uLX, uLX -> uLXX, uLXX -> uLXXX,
  uT -> uTX, uTX -> uTXX, uTXX -> uTXXX,
  thetaAmp -> thetaX, thetaX -> thetaXX, thetaXX -> thetaXXX,
  eAmp -> eX, eX -> eXX, eXX -> eXXX};
totalXDerivative[expression_] := Expand[Total[MapThread[
  D[expression, #1] #2 &, {jetVariables, jetVariables /. jetDerivativeRules}]]];
planeWaveRules = {uLX -> I k uL, uLXX -> -k^2 uL,
  uLXXX -> -I k^3 uL, uTX -> I k uT, uTXX -> -k^2 uT,
  uTXXX -> -I k^3 uT, thetaX -> I k thetaAmp,
  thetaXX -> -k^2 thetaAmp, thetaXXX -> -I k^3 thetaAmp,
  eX -> I k eAmp, eXX -> -k^2 eAmp, eXXX -> -I k^3 eAmp};

functionalDerivativeFourier[density_, field_, gradient_] :=
  Expand[(D[density, field] - totalXDerivative[D[density, gradient]]) /.
    planeWaveRules];

deriveInternalOperators[sourceEnergy_] := Module[
  {density, muTheta, pW, explicitUL, explicitUT, constrainedUL,
   constrainedThickness, transverse},
  density = energy1D[sourceEnergy];
  muTheta = functionalDerivativeFourier[density, thetaAmp, thetaX];
  pW = functionalDerivativeFourier[density, eAmp, eX];
  explicitUL = functionalDerivativeFourier[density, uL, uLX];
  explicitUT = functionalDerivativeFourier[density, uT, uTX];
  constrainedUL = Expand[explicitUL + I k muTheta];
  constrainedThickness = Expand[(pW - muTheta)/w0];
  transverse = Expand[explicitUT];
  <|"ENERGY_1D" -> density, "MU_THETA" -> muTheta, "P_W" -> pW,
    "EXPLICIT_U_L" -> explicitUL, "INPLANE_INTERNAL" -> constrainedUL,
    "THICKNESS_INTERNAL" -> constrainedThickness,
    "TRANSVERSE_INTERNAL" -> transverse|>
];

internalOperators = deriveInternalOperators[fullEnergyDensity];

(* One pre-elimination face-law object feeds the solver, traction assembly,
   and kernel extraction.  Its affinity and slab velocity remain inert until
   the solve step, so their constitutive coefficients stay inspectable. *)
rawFaceLawAssociation[muThetaInput_, velocityInput_, laInput_, lvInput_,
    lxInput_, slabAffinityFactor_] := <|
  "AFFINITY_DEFINITION" ->
    faceAffinityRaw == slabAffinityFactor muThetaInput/rhoBr -
      facePressure/rhoM,
  "SLAB_VELOCITY_DEFINITION" -> faceVelocityRaw == velocityInput,
  "IMPEDANCE_LAW" -> facePressure == zFace faceBulkVelocity,
  "KINEMATIC_LAW" ->
    faceBulkVelocity == faceVelocityRaw + faceFlux/rhoM,
  "CLOSURE_LAW" ->
    faceFlux == laInput faceAffinityRaw + lvInput faceVelocityRaw,
  "TRACTION_LAW" ->
    faceTractionRaw == facePressure + lxInput faceAffinityRaw|>;

solveFace[muThetaInput_, velocityInput_, laInput_, lvInput_, lxInput_,
    slabAffinityFactor_] := Module[{rawLaws, eliminationRules, equations,
    solution, pResult, jResult, vbResult, affinityResult, tractionResult},
  rawLaws = rawFaceLawAssociation[muThetaInput, velocityInput, laInput,
    lvInput, lxInput, slabAffinityFactor];
  eliminationRules = {
    faceAffinityRaw -> Last[rawLaws["AFFINITY_DEFINITION"]],
    faceVelocityRaw -> Last[rawLaws["SLAB_VELOCITY_DEFINITION"]]};
  equations = ({rawLaws["IMPEDANCE_LAW"],
      rawLaws["KINEMATIC_LAW"], rawLaws["CLOSURE_LAW"]} /.
    eliminationRules);
  solution = First[Solve[equations,
    {facePressure, faceFlux, faceBulkVelocity}]];
  pResult = Cancel[facePressure /. solution];
  jResult = Cancel[faceFlux /. solution];
  vbResult = Cancel[faceBulkVelocity /. solution];
  affinityResult = Cancel[faceAffinityRaw /. eliminationRules /. solution];
  tractionResult = Cancel[faceTractionRaw /.
    First[Solve[rawLaws["TRACTION_LAW"], faceTractionRaw]] /.
    {faceAffinityRaw -> affinityResult, facePressure -> pResult}];
  <|"RAW_FACE_LAWS" -> rawLaws,
    "ELIMINATION_RULES" -> eliminationRules,
    "EQUATIONS" -> equations, "SOLUTION" -> solution,
    "PRESSURE" -> pResult, "FLUX" -> jResult,
    "BULK_VELOCITY" -> vbResult, "AFFINITY" -> affinityResult,
    "TRACTION" -> tractionResult|>
];

preEliminationEOMAssociation[internal_Association, tractionFactor_,
    fluxFactor_] := <|
  "INPLANE_EQUATION" -> Expand[-rhoBr omega^2 uL +
    internal["INPLANE_INTERNAL"]],
  "THICKNESS_EQUATION" -> Together[-muW omega^2 w0 eAmp +
    internal["THICKNESS_INTERNAL"] + tractionFactor
      (tractionPlusPower + tractionMinusPower)/2],
  "MASS_EQUATION_WITH_THICKNESS" -> Together[-I omega rhoBr
    (thetaAmp + eAmp + I k uL) + fluxFactor
      (jPlusPower + jMinusPower)],
  "MASS_EQUATION_FROZEN_THICKNESS" -> Together[-I omega rhoBr
    (thetaAmp + I k uL) + fluxFactor
      (jPlusPower + jMinusPower)]|>;

assembleSystem[sourceEnergy_, laInput_, lvInput_, lxInput_,
    slabAffinityFactor_, includeThickness_] := Module[
  {internal, face, vInput, preEliminationEOM, symmetricFaceRules,
   inplaneEquation, thicknessEquation, massEquation, unknowns, equations,
   matrix, determinant},
  internal = deriveInternalOperators[sourceEnergy];
  vInput = If[TrueQ[includeThickness], -I omega w0 eAmp/2, 0];
  face = solveFace[internal["MU_THETA"], vInput, laInput, lvInput,
    lxInput, slabAffinityFactor];
  preEliminationEOM = preEliminationEOMAssociation[internal, 1, 1];
  symmetricFaceRules = {tractionPlusPower -> face["TRACTION"],
    tractionMinusPower -> face["TRACTION"],
    jPlusPower -> face["FLUX"], jMinusPower -> face["FLUX"]};
  inplaneEquation = preEliminationEOM["INPLANE_EQUATION"];
  If[TrueQ[includeThickness],
    thicknessEquation = Together[
      preEliminationEOM["THICKNESS_EQUATION"] /. symmetricFaceRules];
    massEquation = Together[
      preEliminationEOM["MASS_EQUATION_WITH_THICKNESS"] /.
        symmetricFaceRules];
    unknowns = {uL, eAmp, thetaAmp};
    equations = {inplaneEquation, thicknessEquation, massEquation},
    massEquation = Together[
      preEliminationEOM["MASS_EQUATION_FROZEN_THICKNESS"] /.
        symmetricFaceRules];
    unknowns = {uL, thetaAmp};
    equations = {inplaneEquation, massEquation}
  ];
  matrix = coefficientVector[#, unknowns] & /@ equations;
  determinant = Cancel[Together[Det[matrix]]];
  <|"INTERNAL" -> internal, "FACE" -> face,
    "PRE_ELIMINATION_EOM" -> preEliminationEOM,
    "INPLANE_EQUATION" -> inplaneEquation,
    "THICKNESS_EQUATION" -> thicknessEquation,
    "MASS_EQUATION" -> massEquation, "UNKNOWNS" -> unknowns,
    "MATRIX" -> matrix, "DETERMINANT" -> determinant,
    "DISPERSION_NUMERATOR" -> Numerator[Together[determinant]]|>
];

zFace = rhoM omega/q;
fullSystem = assembleSystem[fullEnergyDensity, lambdaA, lambdaV, lambdaX,
  1, True];
placeholderSystem = assembleSystem[fullEnergyDensity, ellA, ellV, ellX,
  1, True];

premiseInventory = {
  HoldComplete[harmonicFactor], HoldComplete[qAlgebraicEquation],
  HoldComplete[acousticEquations], HoldComplete[rhoBrDefinition],
  HoldComplete[eWDefinition], HoldComplete[faceDisplacements],
  HoldComplete[outwardFaceVelocities], HoldComplete[lambdaA],
  HoldComplete[lambdaV], HoldComplete[lambdaX],
  HoldComplete[suppliedEnergyDensity], HoldComplete[$Assumptions]
};
emitLocal["PREMISE_INVENTORY", premiseInventory];

(* ---------------------------------------------------------------------- *)
(* Locus protocol.                                                        *)
(* ---------------------------------------------------------------------- *)

equationZeroForm[equal_Equal] := equal[[1]] - equal[[2]];
equationZeroForm[expression_] := expression;
constitutiveFluxExpression[closure_Equal, fluxSymbol_] := Module[
  {fluxFreeSides},
  fluxFreeSides = Select[List @@ closure, FreeQ[#, fluxSymbol] &];
  If[Length[fluxFreeSides] == 1, First[fluxFreeSides],
    Missing["ConstitutiveFluxExpression", HoldForm[closure]]]
];
swapEquationSides[equal_Equal] := Apply[Equal, Reverse[List @@ equal]];

emitLocus[baseName_String, equations_List, variables_List,
    realPremises_] := Module[
  {zeroForms, solution, identityObject, unrestrictedReduce,
   inconsistentObject, realAdmissible},
  zeroForms = Together /@ (equationZeroForm /@ equations);
  solution = Solve[equations, variables];
  identityObject = casTest[And @@ Thread[zeroForms == 0]];
  unrestrictedReduce = Reduce[And @@ equations, variables];
  inconsistentObject = casTest[unrestrictedReduce === False];
  realAdmissible = Map[
    Function[branch, Module[{testObject, status},
      testObject = Reduce[realPremises && And @@ Thread[variables ==
          (variables /. branch)], variables, Reals];
      status = Which[
        TrueQ[testObject === False], "EXCLUDED",
        TrueQ[testObject === True] || Head[testObject] =!= Reduce,
          "ADMISSIBLE",
        True, "UNDECIDED"];
      <|"BRANCH" -> branch, "STATUS_TOKEN" -> status,
        "TEST_OBJECT" -> testObject,
        "OPERANDS" -> {realPremises, branch, variables}|>
    ]], solution];
  emitShared[baseName <> "_EQUATIONS", equations];
  emitShared[baseName <> "_SOLUTION",
    <|"VARIABLES" -> variables, "SOLUTION" -> solution|>];
  emitShared[baseName <> "_IDENTICALLY_SATISFIED", <|
    "STATUS_TOKEN" -> tokenFromBoolean[identityObject],
    "TEST_OBJECT" -> identityObject, "OPERANDS" -> zeroForms|>];
  emitShared[baseName <> "_INCONSISTENT", <|
    "STATUS_TOKEN" -> tokenFromBoolean[inconsistentObject],
    "TEST_OBJECT" -> inconsistentObject,
    "OPERANDS" -> {equations, variables, unrestrictedReduce}|>];
  emitShared[baseName <> "_REAL_ADMISSIBLE", realAdmissible];
];

(* ---------------------------------------------------------------------- *)
(* B0a: projection and parity.                                            *)
(* ---------------------------------------------------------------------- *)

staticContinuityIntegrand = omegaWindow[w]
   (D[rhoField[w, x, t], t] + divX[jXField[w, x, t]] +
     D[jWField[w, x, t], w]);
finiteProjectionLeft = Inactive[Integrate][
   omegaWindow[w] (D[rhoField[w, x, t], t] + divX[jXField[w, x, t]]),
   {w, w1, w2}];
finiteProjectionSource = Expand[
  -omegaWindow[w2] jWField[w2, x, t] +
   omegaWindow[w1] jWField[w1, x, t] +
   Inactive[Integrate][D[omegaWindow[w], w] jWField[w, x, t],
     {w, w1, w2}]];
finiteProjection = finiteProjectionLeft == finiteProjectionSource;
infiniteProjection = finiteProjection /. {w1 -> -Infinity, w2 -> Infinity,
    omegaWindow[-Infinity] -> 0, omegaWindow[Infinity] -> 0};

evenCurrentRule = jWField[arg_, x, t] :> jEven[arg^2, x, t];
oddCurrentRule = jWField[arg_, x, t] :> arg jOdd[arg^2, x, t];
symmetricProjectionSource = finiteProjectionSource /. {w1 -> -ell,
    w2 -> ell};
evenParitySource = symmetricProjectionSource /. evenCurrentRule;
oddParitySource = symmetricProjectionSource /. oddCurrentRule;
evenParityReflection = Simplify[evenParitySource -
   (evenParitySource /. w -> -w)];
oddParityReflection = Simplify[oddParitySource +
   (oddParitySource /. w -> -w)];

dynamicWindowExtraTerms = Expand[
  Inactive[Integrate][D[omegaDynamic[w, x, t], t] rhoField[w, x, t] +
    gradX[omegaDynamic[w, x, t]].jXVector[w, x, t], {w, w1, w2}]];

emitShared["PROJECTION_FINITE", <|"INTEGRATED_EQUATION" ->
   Inactive[Integrate][staticContinuityIntegrand, {w, w1, w2}] == 0,
  "PROJECTED_EQUATION" -> finiteProjection|>];
emitShared["PROJECTION_INFINITE", infiniteProjection];
emitShared["PARITY_EVEN_JW", <|"SOURCE" -> evenParitySource,
  "REFLECTION_OPERANDS" -> {evenParitySource,
    evenParitySource /. w -> -w}, "REFLECTION_RESIDUAL" ->
    evenParityReflection, "TEST_OBJECT" -> casTest[evenParityReflection == 0]|>];
emitShared["PARITY_ODD_JW", <|"SOURCE" -> oddParitySource,
  "REFLECTION_OPERANDS" -> {oddParitySource,
    -(oddParitySource /. w -> -w)}, "REFLECTION_RESIDUAL" ->
    oddParityReflection, "TEST_OBJECT" -> casTest[oddParityReflection == 0]|>];
emitShared["PARITY_INTERVAL", <|"INTERVAL" -> {-ell, ell},
  "MIDPOINT" -> Simplify[Mean[{-ell, ell}]],
  "SYMMETRY_TEST_OBJECT" -> casTest[Mean[{-ell, ell}] == 0],
  "FINITE_INTERVAL_OBJECT" -> finiteProjectionSource|>];
emitShared["DYNAMIC_WINDOW_EXTRA_TERMS", dynamicWindowExtraTerms];

fluxChannelDefinitions = {
  slabAccretion == -(jPlus + jMinus),
  throughFlow == (jPlus - jMinus)/2
};
fluxChannelMatrix = coefficientVector[#,
    {jPlus, jMinus}] & /@ ({slabAccretion, throughFlow} /.
      First[Solve[fluxChannelDefinitions, {slabAccretion, throughFlow}]]);
emitShared["FLUX_CHANNELS", <|"DEFINING_RELATIONS" -> fluxChannelDefinitions,
  "TRANSFORMATION" -> fluxChannelMatrix,
  "INVERSE_TRANSFORMATION" -> Inverse[fluxChannelMatrix]|>];

(* ---------------------------------------------------------------------- *)
(* B0b and branch prescription: outgoing impermeable response.            *)
(* ---------------------------------------------------------------------- *)

upperAcousticAnsatz = aPlus Exp[I (k x + q (w - w0/2) - omega t)];
lowerAcousticAnsatz = aMinus Exp[I (k x - q (w + w0/2) - omega t)];
upperFaceObjects = {
  phiAmplitude -> aPlus,
  outwardBulkVelocity -> (D[upperAcousticAnsatz, w]/harmonicFactor /.
      w -> w0/2),
  pressureAmplitude -> (-rhoM D[upperAcousticAnsatz, t]/harmonicFactor /.
      w -> w0/2)};
lowerFaceObjects = {
  phiAmplitude -> aMinus,
  outwardBulkVelocity -> (-D[lowerAcousticAnsatz, w]/harmonicFactor /.
      w -> -w0/2),
  pressureAmplitude -> (-rhoM D[lowerAcousticAnsatz, t]/harmonicFactor /.
      w -> -w0/2)};
zUpper = Cancel[(pressureAmplitude/outwardBulkVelocity) /. upperFaceObjects];
zLower = Cancel[(pressureAmplitude/outwardBulkVelocity) /. lowerFaceObjects];
zImpermeable = Simplify[zUpper];

qPropagating = Sign[omega] Sqrt[qSquared];
qEvanescent = I Sqrt[-qSquared];
zPropagating = Refine[zImpermeable /. q -> qPropagating,
  Element[omega, Reals] && qSquared > 0 && rhoM > 0 && cs0 > 0];
zEvanescent = Refine[zImpermeable /. q -> qEvanescent,
  Element[omega, Reals] && qSquared < 0 && rhoM > 0 && cs0 > 0];
zGrazing = Limit[zImpermeable, q -> 0, Direction -> "FromAbove"];
zRegimeObjects = {
  qSquared > 0 -> <|"Q_OUT" -> qPropagating, "Z" -> zPropagating,
    "RE_IM" -> ComplexExpand[{Re[zPropagating], Im[zPropagating]}]|>,
  qSquared < 0 -> <|"Q_OUT" -> qEvanescent, "Z" -> zEvanescent,
    "RE_IM" -> ComplexExpand[{Re[zEvanescent], Im[zEvanescent]}]|>,
  qSquared == 0 -> <|"Q_OUT" -> 0, "Z_LIMIT" -> zGrazing|>
};

thicknessVelocityPair = {-I omega dW/2, -I omega dW/2};
centerVelocityPair = {-I omega zetaC, I omega zetaC};
pressureByParity = {
  thicknessChannel -> Simplify[zImpermeable thicknessVelocityPair],
  centerChannel -> Simplify[zImpermeable centerVelocityPair]
};
addedMassEquation = (zEvanescent (-I omega etaFace)) ==
  addedMass (-omega^2 etaFace);
addedMassSolution = First[Solve[addedMassEquation, addedMass]];

realAxisBranchObject = {
  omega > cs0 k -> Limit[qUpperRim[omega + I epsilon], epsilon -> 0,
    Direction -> "FromAbove"],
  -cs0 k < omega < cs0 k -> Limit[qUpperRim[omega + I epsilon],
    epsilon -> 0, Direction -> "FromAbove"],
  omega < -cs0 k -> Limit[qUpperRim[omega + I epsilon], epsilon -> 0,
    Direction -> "FromAbove"]
};
requiredRealAxisObject = {
  omega > cs0 k -> Sqrt[qSquared],
  -cs0 k < omega < cs0 k -> I Sqrt[-qSquared],
  omega < -cs0 k -> -Sqrt[qSquared]
};
branchRealAxisResidual = MapThread[
  Last[#1] - Last[#2] &, {realAxisBranchObject, requiredRealAxisObject}];

emitShared["BRANCH_REALAXIS_CHECK", <|
  "UPPER_RIM_OPERAND" -> realAxisBranchObject,
  "RADIATION_DECAY_OPERAND" -> requiredRealAxisObject,
  "RESIDUAL" -> branchRealAxisResidual,
  "TEST_OBJECT" -> casTest[And @@ Thread[branchRealAxisResidual == 0]],
  "UPPER_FACE_NORMAL_WAVENUMBER" -> q,
  "LOWER_FACE_GLOBAL_W_WAVENUMBER" -> -q|>];
emitShared["BRANCH_DEGENERATE_POINT", <|
  "POINTS" -> Solve[qSquared == 0, omega],
  "Q_OUT" -> (q /. Solve[q == 0, q]),
  "BULK_SOLUTION_LIMIT" ->
    Series[Exp[I q (w - w0/2)], {q, 0, 1}]|>];
emitShared["Z_IMPERMEABLE", <|"UPPER" -> zUpper, "LOWER" -> zLower,
  "HALF_SPACE_RESIDUAL" -> Together[zUpper - zLower]|>];
emitShared["Z_BY_REGIME", zRegimeObjects];
emitShared["Z_BY_PARITY", pressureByParity];
emitShared["ADDED_MASS", <|"DEFINING_EQUATION" -> addedMassEquation,
  "PER_FACE_OUTWARD_COORDINATE_SOLUTION" -> addedMassSolution|>];
emitShared["GRAZING_BEHAVIOUR", <|"Q_LIMIT" -> 0,
  "Z_LIMIT" -> zGrazing,
  "PRESSURE_LIMIT_AT_FIXED_FACE_VELOCITY" ->
    Limit[zImpermeable velocityFace, q -> 0],
  "PRESSURE_LIMIT_AT_FIXED_POTENTIAL" ->
    Limit[I rhoM omega potentialFace, q -> 0]|>];

(* ---------------------------------------------------------------------- *)
(* B0c: permeable response and its equation-system locus.                 *)
(* ---------------------------------------------------------------------- *)

generalFace = solveFace[muThetaFace, velocityFace, lambdaA, lambdaV,
  lambdaX, 1];
facePressureCoefficients = AssociationThread[{velocityFace, muThetaFace},
  coefficientVector[generalFace["PRESSURE"],
    {velocityFace, muThetaFace}]];
faceEquationZeroForms = equationZeroForm /@ generalFace["EQUATIONS"];
faceCoefficientMatrix = coefficientVector[#,
    {facePressure, faceFlux, faceBulkVelocity}] & /@ faceEquationZeroForms;
faceCoefficientDeterminant = Factor[Det[faceCoefficientMatrix]];
degenerateEquation = Numerator[Together[faceCoefficientDeterminant]] == 0;

emitShared["FACE_RESPONSE", <|"SUPPLIED_EQUATIONS" ->
  generalFace["EQUATIONS"], "SOLUTION" -> generalFace["SOLUTION"],
  "PRESSURE" -> generalFace["PRESSURE"],
  "FLUX" -> generalFace["FLUX"],
  "BULK_VELOCITY" -> generalFace["BULK_VELOCITY"]|>];
emitShared["FACE_RESPONSE_COEFFS", facePressureCoefficients];
emitLocus["DEGENERATE_LOCI", {degenerateEquation}, {lambdaA0},
  $Assumptions && Element[omega, Reals] && Element[q, Complexes] &&
   qAlgebraicEquation];

dissipativeCoefficientByRegime = Table[
  With[{regimeZ = Last[regime]["Z"] /. zGrazing -> DirectedInfinity[]},
    First[regime] -> AssociationThread[{thicknessChannel, centerChannel},
      Table[ComplexExpand[Re[facePressureCoefficients[velocityFace] /.
          zFace -> regimeZ]], 2]]],
  {regime, Take[zRegimeObjects, 2]}];
dissipativeGrazing = qSquared == 0 ->
  Limit[ComplexExpand[Re[facePressureCoefficients[velocityFace]]], q -> 0];
dissipationTauObjects = <|
  tauA -> <|"SMALL_ARGUMENT" -> Block[{$Assumptions = True}, Limit[
      facePressureCoefficients[velocityFace], tauA -> 0]],
    "FINITE_ARGUMENT" -> facePressureCoefficients[velocityFace],
    "LARGE_ARGUMENT" -> Block[{$Assumptions = True}, Limit[
      facePressureCoefficients[velocityFace], tauA -> Infinity]]|>,
  tauV -> <|"SMALL_ARGUMENT" -> Block[{$Assumptions = True}, Limit[
      facePressureCoefficients[velocityFace], tauV -> 0]],
    "FINITE_ARGUMENT" -> facePressureCoefficients[velocityFace],
    "LARGE_ARGUMENT" -> Block[{$Assumptions = True}, Limit[
      facePressureCoefficients[velocityFace], tauV -> Infinity]]|>,
  tauX -> <|"SMALL_ARGUMENT" -> Block[{$Assumptions = True}, Limit[
      facePressureCoefficients[velocityFace], tauX -> 0]],
    "FINITE_ARGUMENT" -> facePressureCoefficients[velocityFace],
    "LARGE_ARGUMENT" -> Block[{$Assumptions = True}, Limit[
      facePressureCoefficients[velocityFace], tauX -> Infinity]]|>
|>;
emitShared["PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY",
  Append[dissipativeCoefficientByRegime, dissipativeGrazing]];
emitShared["PERMEABLE_DISSIPATION_VS_OMEGA_TAU", dissipationTauObjects];

(* ---------------------------------------------------------------------- *)
(* B1: exact balance, linearization, rank and basis reduction.             *)
(* ---------------------------------------------------------------------- *)

sigmaExact = rhoBr (1 + eps (thetaFunction[x, t] + eFunction[x, t]));
velocityExact = eps uVelocityFunction[x, t];
fluxExact = eps (jPlusFunction[x, t] + jMinusFunction[x, t]);
exactBalanceTerms = {
  D[sigmaExact, t],
  D[sigmaExact velocityExact, x],
  fluxExact
};
linearBalanceTerms = Coefficient[Normal[Series[#, {eps, 0, 1}]], eps] & /@
  exactBalanceTerms;
linearBalance = Total[linearBalanceTerms] == 0;
fourierBalanceRules = {
  Derivative[0, 1][thetaFunction][x, t] -> -I omega thetaAmp,
  Derivative[0, 1][eFunction][x, t] -> -I omega eAmp,
  Derivative[1, 0][uVelocityFunction][x, t] -> omega k uL,
  jPlusFunction[x, t] + jMinusFunction[x, t] -> jPlus + jMinus
};
linearConstraint = Together[Total[linearBalanceTerms] /.
   fourierBalanceRules];
constraintCoefficientRow = coefficientVector[linearConstraint,
  {uL, eAmp, thetaAmp, jPlus, jMinus}];
constraintRankGeneric = MatrixRank[{constraintCoefficientRow}];
constraintRankStatic = MatrixRank[{constraintCoefficientRow /. omega -> 0}];

emitShared["CONSTRAINT", <|"EXACT_BALANCE" -> Total[exactBalanceTerms] == 0,
  "LINEAR_BALANCE" -> linearBalance,
  "FOURIER_BALANCE" -> linearConstraint == 0,
  "ASSEMBLED_FACE_BALANCE" -> fullSystem["MASS_EQUATION"] == 0|>];
emitShared["CONSTRAINT_TERM_ORIGINS", MapThread[
  <|"EXACT_TERM" -> #1, "LINEAR_TERM" -> #2,
    "DROPPED_REMAINDER" -> Expand[#1 - eps #2]|> &,
  {exactBalanceTerms, linearBalanceTerms}]];
emitShared["DOF_COUNTING_CONVENTION", <|
  "OBJECTS" -> {uL, uT, eAmp, thetaAmp},
  "COUNTED_AT_FIXED" -> {k, omega},
  "COUNT_OBJECT" -> MatrixRank[IdentityMatrix[4]]|>];
emitShared["INTERNAL_DOF_COUNT", <|
  "BEFORE_CONSTRAINT" -> MatrixRank[IdentityMatrix[4]],
  "CONSTRAINT_RANK_GENERIC" -> constraintRankGeneric,
  "AFTER_CONSTRAINT_GENERIC" -> 4 - constraintRankGeneric,
  "CONSTRAINT_RANK_AT_ZERO_FREQUENCY" -> constraintRankStatic,
  "AFTER_CONSTRAINT_AT_ZERO_FREQUENCY" -> 4 - constraintRankStatic|>];

basisOneD = energy1D /@ energyBasis;
basisFourier = Expand[# /. planeWaveRules] & /@ basisOneD;
impermeableConstraintSolution = First[Solve[
  -I omega rhoBr (thetaAmp + eAmp + I k uL) == 0, thetaAmp]];
fluxConstraintSolution = Solve[fullSystem["MASS_EQUATION"] == 0, thetaAmp];
fluxConstraintRule = If[Length[fluxConstraintSolution] > 0,
  First[fluxConstraintSolution], {}];
impermeableReducedBasis = Cancel[# /. impermeableConstraintSolution] & /@
  basisFourier;
fluxReducedBasis = Cancel[# /. fluxConstraintRule] & /@ basisFourier;
reducedMonomials = Flatten[Table[{uL, uT, eAmp}[[i]]
    {uL, uT, eAmp}[[j]], {i, 3}, {j, i, 3}]];
impermeableReductionMatrix = coefficientVector[#,
    reducedMonomials] & /@ impermeableReducedBasis;
fluxReductionMatrix = coefficientVector[#, reducedMonomials] & /@
  fluxReducedBasis;
impermeableRedundancies = NullSpace[Transpose[impermeableReductionMatrix]];
fluxRedundancies = NullSpace[Transpose[fluxReductionMatrix]];
emitShared["BASIS_REDUNDANCY_UNDER_CONSTRAINT", <|
  "BASIS_ORDER" -> energyBasis,
  "IMPERMEABLE_REDUCED_BASIS" -> impermeableReducedBasis,
  "IMPERMEABLE_REDUNDANCIES" -> impermeableRedundancies,
  "FLUX_ON_REDUCED_BASIS" -> fluxReducedBasis,
  "FLUX_ON_REDUNDANCIES" -> fluxRedundancies,
  "ROUTE_DIFFERENCE" -> comparisonRecord[
    RowReduce[impermeableReductionMatrix], RowReduce[fluxReductionMatrix]]|>];

(* ---------------------------------------------------------------------- *)
(* B2a: equations generated by the virtual constraint and balance laws.   *)
(* ---------------------------------------------------------------------- *)

lambdaXZeroSystem = assembleSystem[fullEnergyDensity, lambdaA, lambdaV, 0,
  1, True];
inplaneFixedThetaCandidate = Expand[-rhoBr omega^2 uL +
  internalOperators["EXPLICIT_U_L"]];
inplaneMuGradient = Expand[fullSystem["INPLANE_EQUATION"] -
  inplaneFixedThetaCandidate];

emitShared["INPLANE_EOM", <|
  "HELD_FIXED_IN_MU_THETA" -> {uL, eAmp},
  "MU_THETA_FUNCTIONAL_DERIVATIVE" -> internalOperators["MU_THETA"],
  "HELD_FIXED_IN_EXPLICIT_U_VARIATION" -> {thetaAmp, eAmp},
  "OPERATOR" -> fullSystem["INPLANE_EQUATION"],
  "EQUATION" -> fullSystem["INPLANE_EQUATION"] == 0|>];
emitShared["THICKNESS_EOM", <|
  "HELD_FIXED_IN_P_W" -> {uL, thetaAmp},
  "P_W_FUNCTIONAL_DERIVATIVE" -> internalOperators["P_W"],
  "VIRTUAL_CONSTRAINT_COMBINATION" ->
    internalOperators["THICKNESS_INTERNAL"],
  "OPERATOR" -> fullSystem["THICKNESS_EQUATION"],
  "EQUATION" -> fullSystem["THICKNESS_EQUATION"] == 0|>];
emitShared["BULK_FORCE_ON_THICKNESS", <|
  "PRESSURE_PART" -> fullSystem["FACE"]["PRESSURE"],
  "AFFINITY_TRACTION_PART" ->
    Expand[lambdaX fullSystem["FACE"]["AFFINITY"]],
  "TWO_FACE_VIRTUAL_WORK_COEFFICIENT" ->
    -fullSystem["FACE"]["TRACTION"]|>];
emitShared["RECIPROCAL_TRACTION_THICKNESS_EFFECT",
  comparisonRecord[fullSystem["THICKNESS_EQUATION"],
    lambdaXZeroSystem["THICKNESS_EQUATION"]]];

emitShared["CONVENTION_CHECK_INPLANE", <|
  "VIRTUAL_CONSTRAINT_OPERATOR" -> fullSystem["INPLANE_EQUATION"],
  "FIXED_THETA_CANDIDATE" -> inplaneFixedThetaCandidate,
  "DIFFERENCE" -> inplaneMuGradient,
  "MU_GRADIENT_OPERAND" -> I k internalOperators["MU_THETA"],
  "MU_GRADIENT_RESIDUAL" -> Together[inplaneMuGradient -
     I k internalOperators["MU_THETA"]],
  "TEST_OBJECT" -> casTest[Together[inplaneMuGradient -
      I k internalOperators["MU_THETA"]] == 0]|>];

(* Conservative convention check from the independently reduced energy. *)
uniformStoredEnergy = Expand[energy1D[fullEnergyDensity] /.
  {uLX -> 0, uTX -> 0, thetaX -> 0, eX -> 0,
   thetaAmp -> conservedSigma - eAmp}];
kCheck = D[uniformStoredEnergy, {eAmp, 2}];
uniformInternal = internalOperators["THICKNESS_INTERNAL"] /.
  {k -> 0, uL -> 0, thetaAmp -> conservedSigma - eAmp};
equationStiffness = Simplify[D[uniformInternal, eAmp] w0];
conservativeOmegaSquared = Simplify[equationStiffness/(muW w0^2)];
uniformEnergyHessian = D[
  energy1D[fullEnergyDensity] /. {uLX -> 0, uTX -> 0, thetaX -> 0,
    eX -> 0}, {{thetaAmp, eAmp}, 2}];
uniformPositiveDefinite = And @@ {
  uniformEnergyHessian[[1, 1]] > 0, Det[uniformEnergyHessian] > 0};
conservativeInequality = Reduce[
  conservativeOmegaSquared > 0 && w0 > 0 && muW > 0,
  {br3, crossC, kW}, Reals];

emitShared["CONVENTION_CHECK_CONSERVATIVE", <|
  "REDUCED_ENERGY" -> uniformStoredEnergy,
  "ENERGY_ROUTE_STIFFNESS" -> kCheck,
  "EQUATION_ROUTE_STIFFNESS" -> equationStiffness,
  "RESIDUAL" -> Together[equationStiffness - kCheck],
  "TEST_OBJECT" -> casTest[Together[equationStiffness - kCheck] == 0],
  "OMEGA_SQUARED" -> conservativeOmegaSquared,
  "BR3_DEPENDENCE" -> D[conservativeOmegaSquared, br3]|>];
emitShared["CONSERVATIVE_POSITIVITY_INEQUALITY", <|
  "POSITIVE_DEFINITE_OPERAND" -> uniformPositiveDefinite,
  "MODE_POSITIVITY_OPERAND" -> conservativeOmegaSquared > 0,
  "REDUCED_CONDITION" -> conservativeInequality|>];

(* ---------------------------------------------------------------------- *)
(* B2b: thickness response and bulk operator by acoustic regime.          *)
(* ---------------------------------------------------------------------- *)

thicknessDrivenSolution = First[Solve[
  fullSystem["THICKNESS_EQUATION"] == thicknessDrive, eAmp]];
thicknessResponse = Cancel[(w0 eAmp/thicknessDrive) /.
  thicknessDrivenSolution];
thicknessResponseNormalization = HoldForm[
  thicknessDisplacementAmplitude/thicknessGeneralizedPressureAmplitude];

bulkOperatorByRegime = Map[
  Function[regime, Module[{condition, zValue, pressure, force,
      accelerationRatio, velocityRatio},
    condition = First[regime];
    zValue = If[KeyExistsQ[Last[regime], "Z"], Last[regime]["Z"],
      Last[regime]["Z_LIMIT"]];
    pressure = Cancel[fullSystem["FACE"]["PRESSURE"] /. zFace -> zValue];
    force = Cancel[fullSystem["FACE"]["TRACTION"] /. zFace -> zValue];
    accelerationRatio = Cancel[force/(-omega^2 w0 eAmp)];
    velocityRatio = Cancel[force/(-I omega w0 eAmp)];
    condition -> <|"BULK_FORCE" -> force,
      "ACCELERATION_PHASE_RATIO" -> accelerationRatio,
      "VELOCITY_PHASE_RATIO" -> velocityRatio|>
  ]], zRegimeObjects];
massInterpretationTests = Map[
  Function[regime, First[regime] -> casTest[
    ComplexExpand[Re[Last[regime][If[KeyExistsQ[Last[regime], "Z"],
       "Z", "Z_LIMIT"]]]] == 0]], zRegimeObjects];

emitShared["THICKNESS_RESPONSE", thicknessResponse];
emitShared["RESPONSE_NORMALIZATION", <|
  "RATIO" -> thicknessResponseNormalization,
  "OUTPUT" -> w0 eAmp, "INPUT" -> thicknessDrive|>];
emitShared["BULK_OPERATOR_BY_REGIME", bulkOperatorByRegime];
emitShared["MASS_INTERPRETATION_VALID_WHERE", massInterpretationTests];

(* ---------------------------------------------------------------------- *)
(* B2c: eliminate thickness and density for the compression response.     *)
(* ---------------------------------------------------------------------- *)

strainRule = uL -> strainLong/(I k);
compressionEliminationEquations = {
  fullSystem["THICKNESS_EQUATION"] == 0,
  fullSystem["MASS_EQUATION"] == 0} /. strainRule;
compressionElimination = First[Solve[compressionEliminationEquations,
  {eAmp, thetaAmp}]];
longitudinalStress = Cancel[
  internalOperators["INPLANE_INTERNAL"]/(I k)];
compressionalResponse = Cancel[(longitudinalStress/strainLong) /.
  strainRule /. compressionElimination];

constraintMatrixForElimination = coefficientVector[#,
    {eAmp, thetaAmp}] & /@
  (equationZeroForm /@ compressionEliminationEquations);
constraintRankBeforeZero = MatrixRank[constraintMatrixForElimination];
constraintRankAfterZero = MatrixRank[constraintMatrixForElimination /.
  omega -> 0];
fixedKLowFrequency = Limit[compressionalResponse, omega -> 0];
fixedKHighFrequency = Limit[compressionalResponse, omega -> Infinity];
rayLowFrequency = Limit[compressionalResponse /. k -> raySlope omega,
  omega -> 0];
rayHighFrequency = Limit[compressionalResponse /. k -> raySlope omega,
  omega -> Infinity];

frozenThicknessSystem = assembleSystem[fullEnergyDensity, lambdaA, lambdaV,
  lambdaX, 1, False];
frozenThetaSolution = First[Solve[
  (frozenThicknessSystem["MASS_EQUATION"] /. strainRule) == 0, thetaAmp]];
frozenStress = Cancel[(longitudinalStress/strainLong) /. strainRule /.
  eAmp -> 0 /. frozenThetaSolution];

emitShared["COMPRESSIONAL_RESPONSE", <|
  "NORMALIZATION" -> HoldForm[longitudinalStressAmplitude/
     longitudinalStrainAmplitude],
  "STRESS_OPERAND" -> longitudinalStress,
  "DEFORMATION_OPERAND" -> strainLong,
  "ELIMINATION_SOLUTION" -> compressionElimination,
  "RESPONSE" -> compressionalResponse,
  "CONSTRAINT_RANK_GENERIC" -> constraintRankBeforeZero,
  "CONSTRAINT_RANK_AT_ZERO_FREQUENCY" -> constraintRankAfterZero,
  "DIVIDE_THEN_LIMIT" -> fixedKLowFrequency,
  "LIMIT_THEN_DIVIDE" -> Limit[
    Numerator[Together[compressionalResponse]], omega -> 0]/
     Limit[Denominator[Together[compressionalResponse]], omega -> 0]|>];
emitShared["LIMITS_AND_PATH", <|
  "FIXED_K_PATH" -> {fixedKLowFrequency, fixedKHighFrequency},
  "FIXED_OMEGA_TO_K_RATIO_PATH" -> {rayLowFrequency, rayHighFrequency},
  "PATH_RESIDUALS" -> {Together[fixedKLowFrequency - rayLowFrequency],
    Together[fixedKHighFrequency - rayHighFrequency]}|>];
emitShared["FROZEN_THICKNESS_IDENTIFICATION", <|
  "FROZEN_SOURCE_ENERGY" -> (fullEnergyDensity /. eField -> 0 /.
     Thread[eGradient -> 0]),
  "FROZEN_RESPONSE" -> frozenStress,
  "FULL_RESPONSE_RESIDUAL" -> Together[compressionalResponse - frozenStress]|>];

(* Equal-time, slab-affinity-off acceptance slice: no reference is present. *)
rawPressureSliceLaws = generalFace["RAW_FACE_LAWS"] /.
  {muThetaFace -> 0, tauA -> tauCommon, tauV -> tauCommon};
rawPressureSliceClosure = rawPressureSliceLaws["CLOSURE_LAW"] /.
  faceAffinityRaw -> Last[rawPressureSliceLaws["AFFINITY_DEFINITION"]];
rawPressureFluxExpression = constitutiveFluxExpression[
  rawPressureSliceClosure, faceFlux];
dynamicRawPressureCoefficient = Coefficient[
  rawPressureFluxExpression, facePressure];
staticRawPressureCoefficient = Block[{$Assumptions = True},
  Limit[dynamicRawPressureCoefficient, omega -> 0]];
zPermSliceMap = Solve[
  staticRawPressureCoefficient == lambdaPressureCoefficient,
  lambdaPressureCoefficient];
zPermSliceEquations = generalFace["EQUATIONS"] /.
  {muThetaFace -> 0, tauA -> tauCommon, tauV -> tauCommon};
zPermSliceSolution = First[Solve[zPermSliceEquations,
  {facePressure, faceFlux, faceBulkVelocity}]];
zPermSlice = Cancel[facePressure/velocityFace /. zPermSliceSolution];
emitShared["ZPERM_SLICE_MAP", <|
  "STATIC_RAW_PRESSURE_COEFFICIENT" -> staticRawPressureCoefficient,
  "MAPPING" -> zPermSliceMap|>];
emitShared["ZPERM_SLICE_DYNAMIC_COEFFICIENT", <|
  "FLUX_EXPRESSION" -> rawPressureFluxExpression,
  "DYNAMIC_RAW_PRESSURE_COEFFICIENT" -> dynamicRawPressureCoefficient|>];
swappedRawPressureSliceClosure = swapEquationSides[
  rawPressureSliceClosure];
oldPressureExtractorOriginal = Coefficient[
  equationZeroForm[rawPressureSliceClosure], facePressure];
oldPressureExtractorSwapped = Coefficient[
  equationZeroForm[swappedRawPressureSliceClosure], facePressure];
newPressureExtractorOriginal = Coefficient[
  constitutiveFluxExpression[rawPressureSliceClosure, faceFlux],
  facePressure];
newPressureExtractorSwapped = Coefficient[
  constitutiveFluxExpression[swappedRawPressureSliceClosure, faceFlux],
  facePressure];
emitShared["ZPERM_SLICE_MAP_REPRESENTATION_CONTROL", <|
  "ORIGINAL_CLOSURE" -> rawPressureSliceClosure,
  "SIDE_SWAPPED_CLOSURE" -> swappedRawPressureSliceClosure,
  "OLD_RESIDUAL_EXTRACTOR_OPERANDS" ->
    {oldPressureExtractorOriginal, oldPressureExtractorSwapped},
  "OLD_RESIDUAL_EXTRACTOR_DIFFERENCE" -> Together[
    oldPressureExtractorOriginal - oldPressureExtractorSwapped],
  "FLUX_EXTRACTOR_OPERANDS" ->
    {newPressureExtractorOriginal, newPressureExtractorSwapped},
  "FLUX_EXTRACTOR_RESIDUAL" -> Together[
    newPressureExtractorOriginal - newPressureExtractorSwapped],
  "STATIC_PARAMETER_DERIVATIVES" ->
    (D[staticRawPressureCoefficient, #] & /@ {omega, tauCommon})|>];
emitShared["ZPERM_SLICE", zPermSlice];

(* ---------------------------------------------------------------------- *)
(* B2d: two-port power, passivity, admissibility and reciprocity.          *)
(* ---------------------------------------------------------------------- *)

portFace = solveFace[rhoBr muSInput, velocityInput, lambdaA, lambdaV,
  lambdaX, 1];
portOutputVector = {portFace["TRACTION"], portFace["FLUX"]};
portInputVector = {velocityInput, muSInput};
portResponseMatrix = coefficientVector[#, portInputVector] & /@
  portOutputVector;
portHermitianMatrix = Simplify[(portResponseMatrix +
   ConjugateTranspose[portResponseMatrix])/2];
portPrincipalMinors = {
  portHermitianMatrix[[1, 1]], portHermitianMatrix[[2, 2]],
  Det[portHermitianMatrix]};
portDissipativityCondition = And @@ Thread[portPrincipalMinors >= 0];

twoPortLeft = 1/2 Re[(facePressurePower + lambdaX affinityPower)
    Conjugate[velocityPower] + muSPower Conjugate[fluxPower]];
bulkVelocityPower = velocityPower + fluxPower/rhoM;
affinityNormalizationPower = muSPower - facePressurePower/rhoM;
twoPortRight = 1/2 Re[facePressurePower Conjugate[bulkVelocityPower] +
    affinityNormalizationPower Conjugate[fluxPower] +
    lambdaX affinityNormalizationPower Conjugate[velocityPower]];
twoPortResidual = ComplexExpand[twoPortLeft - twoPortRight,
  {facePressurePower, affinityPower, velocityPower, muSPower, fluxPower}] /.
   affinityPower -> affinityNormalizationPower;

mixedResponseMatrix = {{lambdaA, lambdaV}, {lambdaX, 0}};
mixedHermitianMatrix = Simplify[(mixedResponseMatrix +
   ConjugateTranspose[mixedResponseMatrix])/2];
mixedPrincipalMinors = {
  mixedHermitianMatrix[[1, 1]], mixedHermitianMatrix[[2, 2]],
  Det[mixedHermitianMatrix]};

realFrequencyConjugationRules = {
  Conjugate[omega] -> omega, Conjugate[lambdaA0] -> lambdaA0,
  Conjugate[lambdaV0] -> lambdaV0, Conjugate[lambdaX0] -> lambdaX0,
  Conjugate[tauA] -> tauA, Conjugate[tauV] -> tauV,
  Conjugate[tauX] -> tauX, Conjugate[I] -> -I};
conjugateLambdaXReal = ComplexExpand[Conjugate[lambdaX]];
admissibilityCrossExpression = Cancel[Together[
  lambdaV + conjugateLambdaXReal]];
admissibilityCrossPolynomial = Numerator[admissibilityCrossExpression];
admissibilityCrossEquations = Join[
  Thread[CoefficientList[ComplexExpand[Re[admissibilityCrossPolynomial]],
      omega] == 0],
  Thread[CoefficientList[ComplexExpand[Im[admissibilityCrossPolynomial]],
      omega] == 0]];
admissibilityCoefficientRegion = Reduce[
  And @@ Join[admissibilityCrossEquations,
    {lambdaA0 >= 0, tauA >= 0, tauV >= 0, tauX >= 0}],
  {lambdaA0, lambdaV0, lambdaX0, tauA, tauV, tauX}, Reals];

generalMixedMatrix = {{onsagerA, onsagerB}, {onsagerC, onsagerEpsilon}};
allFluxMatrix = {
  {onsagerA - onsagerB onsagerC/onsagerEpsilon,
    onsagerB/onsagerEpsilon},
  {-onsagerC/onsagerEpsilon, 1/onsagerEpsilon}};
stateParities = -{-1, -1};
casimirFactor = Times @@ stateParities;
allFluxCrossEquation = Numerator[Together[
   allFluxMatrix[[1, 2]] - casimirFactor allFluxMatrix[[2, 1]]]] == 0;
allForceMatrix = {
  {1/onsagerA, -onsagerB/onsagerA},
  {onsagerC/onsagerA,
    onsagerEpsilon - onsagerC onsagerB/onsagerA}};
allForceCrossEquation = Numerator[Together[
   allForceMatrix[[1, 2]] - casimirFactor allForceMatrix[[2, 1]]]] == 0;
allFluxCleared = Factor[equationZeroForm[allFluxCrossEquation]];
allForceCleared = Factor[equationZeroForm[allForceCrossEquation]];
allFluxNormalized = normalizedPolynomial[allFluxCleared,
  {onsagerA, onsagerB, onsagerC, onsagerEpsilon}];
allForceNormalized = normalizedPolynomial[allForceCleared,
  {onsagerA, onsagerB, onsagerC, onsagerEpsilon}];
onsagerMixedRelation = allFluxNormalized /. {
  onsagerB -> lambdaV, onsagerC -> lambdaX};
onsagerPolynomial = Numerator[Together[onsagerMixedRelation]];
onsagerCoefficientEquations = Join[
  Thread[CoefficientList[ComplexExpand[Re[onsagerPolynomial]], omega] == 0],
  Thread[CoefficientList[ComplexExpand[Im[onsagerPolynomial]], omega] == 0]];
onsagerRegion = Reduce[
  And @@ Join[onsagerCoefficientEquations,
    {tauV >= 0, tauX >= 0}],
  {lambdaV0, lambdaX0, tauV, tauX}, Reals];
admissibilityWithReciprocity = Reduce[
  admissibilityCoefficientRegion && onsagerRegion,
  {lambdaA0, lambdaV0, lambdaX0, tauA, tauV, tauX}, Reals];
regionSetRelation = <|
  "ADMISSIBILITY_IMPLIES_RECIPROCITY" -> FullSimplify[
    Implies[admissibilityCoefficientRegion, onsagerRegion]],
  "RECIPROCITY_IMPLIES_ADMISSIBILITY" -> FullSimplify[
    Implies[onsagerRegion, admissibilityCoefficientRegion]],
  "INTERSECTION" -> admissibilityWithReciprocity|>;
admissibilityTimeProjection = Resolve[Exists[
  {lambdaA0, lambdaV0, lambdaX0}, admissibilityCoefficientRegion], Reals];
reciprocityTimeProjection = Resolve[Exists[
  {lambdaV0, lambdaX0}, onsagerRegion], Reals];

emitShared["TWO_PORT_POWER_IDENTITY", <|
  "SLAB_VARIABLE_OPERAND" -> twoPortLeft,
  "BULK_INTERFACE_OPERAND" -> twoPortRight,
  "AFFINITY_SUBSTITUTION" -> affinityPower == affinityNormalizationPower,
  "RESIDUAL" -> twoPortResidual,
  "TEST_OBJECT" -> casTest[twoPortResidual == 0]|>];
emitShared["PORT_DISSIPATIVITY", <|
  "INPUT_ORDER" -> portInputVector,
  "OUTPUT_ORDER" -> portOutputVector,
  "RESPONSE_MATRIX" -> portResponseMatrix,
  "HERMITIAN_POWER_MATRIX" -> portHermitianMatrix,
  "PRINCIPAL_MINORS" -> portPrincipalMinors,
  "CONDITION" -> portDissipativityCondition|>];
emitShared["PORT_CONDITION_KIND", <|
  "PARAMETER_DERIVATIVES" ->
    (D[portPrincipalMinors, #] & /@ {lambdaA0, lambdaV0, lambdaX0,
       tauA, tauV, tauX}),
  "FREQUENCY_WAVENUMBER_DERIVATIVES" ->
    (D[portPrincipalMinors, #] & /@ {omega, k, q})|>];
emitShared["ONSAGER_CONDITION", <|
  "MIXED_MATRIX" -> generalMixedMatrix,
  "ALL_FLUX_MATRIX" -> allFluxMatrix,
  "STATE_PARITIES" -> stateParities,
  "CLEARED_RELATION" -> allFluxNormalized|>];
emitShared["ONSAGER_RECIPROCITY", <|
  "ALL_FLUX_RELATION" -> allFluxNormalized,
  "ALL_FORCE_MATRIX" -> allForceMatrix,
  "ALL_FORCE_RELATION" -> allForceNormalized,
  "ROUTE_RESIDUAL" -> Together[allFluxNormalized - allForceNormalized],
  "TEST_OBJECT" -> casTest[Together[allFluxNormalized -
      allForceNormalized] == 0],
  "CONSTITUTIVE_REGION" -> onsagerRegion|>];
emitShared["ONSAGER_DETERMINABLE", <|
  "ALL_FLUX_SOLVE" -> Solve[allFluxCrossEquation, onsagerC],
  "ALL_FORCE_SOLVE" -> Solve[allForceCrossEquation, onsagerC],
  "RESIDUAL" -> Together[allFluxNormalized - allForceNormalized]|>];
emitShared["RELAXATION_TIME_RELATIONS", <|
  "ADMISSIBILITY" -> admissibilityTimeProjection,
  "RECIPROCITY" -> reciprocityTimeProjection|>];
emitShared["COEFFICIENT_ADMISSIBILITY", <|
  "HERMITIAN_POWER_MATRIX" -> mixedHermitianMatrix,
  "PRINCIPAL_MINORS" -> mixedPrincipalMinors,
  "UNCONDITIONAL_REGION" -> admissibilityCoefficientRegion,
  "RECIPROCITY_REGION" -> onsagerRegion,
  "WITH_RECIPROCITY" -> admissibilityWithReciprocity,
  "SET_RELATION" -> regionSetRelation|>];

(* ---------------------------------------------------------------------- *)
(* Causality diagnostic: extraction, inert propagation and pole inventory. *)
(* ---------------------------------------------------------------------- *)

extractFaceKernels[rawLaws_Association] := Module[
  {fluxExpression, tractionExpression},
  fluxExpression = constitutiveFluxExpression[
    rawLaws["CLOSURE_LAW"], faceFlux];
  tractionExpression = constitutiveFluxExpression[
    rawLaws["TRACTION_LAW"], faceTractionRaw];
  <|"A" -> Coefficient[fluxExpression /. faceVelocityRaw -> 0,
      faceAffinityRaw],
    "V" -> Coefficient[fluxExpression /. faceAffinityRaw -> 0,
      faceVelocityRaw],
    "X" -> Coefficient[tractionExpression - facePressure,
      faceAffinityRaw]|>
];
assembledRawFaceLaws = fullSystem["FACE"]["RAW_FACE_LAWS"];
assembledKernelExtraction = extractFaceKernels[assembledRawFaceLaws];
kernelExtracted = <|
  "A_PLUS" -> assembledKernelExtraction["A"],
  "V_PLUS" -> assembledKernelExtraction["V"],
  "X_PLUS" -> assembledKernelExtraction["X"],
  "A_MINUS" -> assembledKernelExtraction["A"],
  "V_MINUS" -> assembledKernelExtraction["V"],
  "X_MINUS" -> assembledKernelExtraction["X"]
|>;
kernelOrientationIdentities = <|
  "A_PLUS" -> rationalIdentityRecord[kernelExtracted["A_PLUS"],
    lambdaA, omega],
  "V_PLUS" -> rationalIdentityRecord[kernelExtracted["V_PLUS"],
    lambdaV, omega],
  "X_PLUS" -> rationalIdentityRecord[kernelExtracted["X_PLUS"],
    lambdaX, omega],
  "A_MINUS" -> rationalIdentityRecord[kernelExtracted["A_MINUS"],
    lambdaA, omega],
  "V_MINUS" -> rationalIdentityRecord[kernelExtracted["V_MINUS"],
    lambdaV, omega],
  "X_MINUS" -> rationalIdentityRecord[kernelExtracted["X_MINUS"],
    lambdaX, omega]
|>;
kernelReplacement = {ellA -> kernelExtracted["A_PLUS"],
  ellV -> kernelExtracted["V_PLUS"], ellX -> kernelExtracted["X_PLUS"]};
propagatedFace = AssociationMap[
  Cancel[placeholderSystem["FACE"][#] /. kernelReplacement] &,
  {"PRESSURE", "FLUX", "BULK_VELOCITY", "AFFINITY", "TRACTION"}];
kernelPropagationResiduals = <|
  "FACE_PRESSURE" -> comparisonRecord[propagatedFace["PRESSURE"],
    fullSystem["FACE"]["PRESSURE"]],
  "FACE_FLUX" -> comparisonRecord[propagatedFace["FLUX"],
    fullSystem["FACE"]["FLUX"]],
  "INPLANE_EQUATION" -> comparisonRecord[
    placeholderSystem["INPLANE_EQUATION"] /. kernelReplacement,
    fullSystem["INPLANE_EQUATION"]],
  "THICKNESS_EQUATION" -> comparisonRecord[
    placeholderSystem["THICKNESS_EQUATION"] /. kernelReplacement,
    fullSystem["THICKNESS_EQUATION"]],
  "MASS_EQUATION" -> comparisonRecord[
    placeholderSystem["MASS_EQUATION"] /. kernelReplacement,
    fullSystem["MASS_EQUATION"]],
  "DISPERSION_DETERMINANT" -> comparisonRecord[
    placeholderSystem["DETERMINANT"] /. kernelReplacement,
    fullSystem["DETERMINANT"]]|>;

bareKernelDenominators = <|"A" -> Denominator[Together[lambdaA]],
  "V" -> Denominator[Together[lambdaV]],
  "X" -> Denominator[Together[lambdaX]]|>;
downstreamKernelObjects = <|
  "FACE_PRESSURE" -> fullSystem["FACE"]["PRESSURE"],
  "FACE_FLUX" -> fullSystem["FACE"]["FLUX"],
  "THICKNESS_EQUATION" -> fullSystem["THICKNESS_EQUATION"],
  "MASS_EQUATION" -> fullSystem["MASS_EQUATION"],
  "DISPERSION_DETERMINANT" -> fullSystem["DETERMINANT"]|>;
kernelPoleInventory = AssociationMap[
  Function[channel, Module[{bareDenominator, bareLocation},
    bareDenominator = bareKernelDenominators[channel];
    bareLocation = Solve[bareDenominator == 0, omega];
    <|"BARE_DENOMINATOR" -> bareDenominator,
      "BARE_LOCATION" -> bareLocation,
      "DOWNSTREAM" -> AssociationMap[
        Function[objectName, Module[{cancelled, denominator, residues},
          cancelled = Cancel[downstreamKernelObjects[objectName]];
          denominator = Factor[Denominator[Together[cancelled]]];
          residues = Limit[cancelled (omega - (omega /. First[bareLocation])),
              omega -> (omega /. First[bareLocation])];
          <|"CANCELLED_OBJECT" -> cancelled,
            "REDUCED_DENOMINATOR" -> denominator,
            "BARE_LOCATION_RESIDUE" -> residues,
            "REDUCED_POLE_EQUATION" -> denominator == 0|>
        ]], Keys[downstreamKernelObjects]]|>
  ]], Keys[bareKernelDenominators]];

requiredOrientationRecords = {
  "A_PLUS", "V_PLUS", "X_PLUS", "A_MINUS", "V_MINUS", "X_MINUS"};
requiredPropagationRecords = {"FACE_PRESSURE", "FACE_FLUX",
  "INPLANE_EQUATION", "THICKNESS_EQUATION", "MASS_EQUATION",
  "DISPERSION_DETERMINANT"};
aggregateNamedTestRecords[records_Association, requiredKeys_List] := Module[
  {lookups, presentRecords, testObjects, presenceTest, aggregateTest},
  lookups = Lookup[records, requiredKeys,
    Missing["RequiredRecord"]];
  presentRecords = Select[lookups,
    AssociationQ[#] && Length[#] > 0 && KeyExistsQ[#, "TEST_OBJECT"] &];
  testObjects = Lookup[presentRecords, "TEST_OBJECT"];
  presenceTest = Length[presentRecords] == Length[requiredKeys];
  aggregateTest = If[TrueQ[presenceTest], And @@ testObjects, False];
  <|"REQUIRED_RECORDS" -> requiredKeys,
    "AVAILABLE_RECORDS" -> Keys[records],
    "PRESENT_RECORD_COUNT" -> Length[presentRecords],
    "REQUIRED_RECORD_COUNT" -> Length[requiredKeys],
    "PRESENCE_TEST_OBJECT" -> presenceTest,
    "TEST_OBJECTS" -> testObjects,
    "AGGREGATE_TEST_OBJECT" -> casTest[aggregateTest]|>
];
orientationAggregation = aggregateNamedTestRecords[
  kernelOrientationIdentities, requiredOrientationRecords];
propagationAggregation = aggregateNamedTestRecords[
  kernelPropagationResiduals, requiredPropagationRecords];
causalityAggregate = casTest[
  orientationAggregation["AGGREGATE_TEST_OBJECT"] &&
  propagationAggregation["AGGREGATE_TEST_OBJECT"]];

kernelAblatedSystem = assembleSystem[fullEnergyDensity,
  lambdaA /. lambdaA0 -> lambdaA0 + kernelAblation, lambdaV, lambdaX,
  1, True];
kernelAblatedExtraction = extractFaceKernels[
  kernelAblatedSystem["FACE"]["RAW_FACE_LAWS"]];
kernelUnrelatedSystem = assembleSystem[
  fullEnergyDensity /. crossC -> crossC + kernelUnrelatedAblation,
  lambdaA, lambdaV, lambdaX, 1, True];
kernelUnrelatedExtraction = extractFaceKernels[
  kernelUnrelatedSystem["FACE"]["RAW_FACE_LAWS"]];
removedOrientationAggregation = aggregateNamedTestRecords[
  KeyDrop[kernelOrientationIdentities, First[requiredOrientationRecords]],
  requiredOrientationRecords];
unrelatedCausalityAggregation = aggregateNamedTestRecords[
  kernelOrientationIdentities, requiredOrientationRecords];
emitShared["KERNEL_ORIENTATION_IDENTITIES", kernelOrientationIdentities];
emitShared["KERNEL_ORIENTATION_CONTROLS", <|
  "ASSEMBLED_KERNEL_ABLATION" -> rationalIdentityRecord[
    kernelAblatedExtraction["A"], lambdaA, omega],
  "UNRELATED_SOURCE_OPERANDS" ->
    {assembledKernelExtraction, kernelUnrelatedExtraction},
  "UNRELATED_SOURCE_RESIDUALS" -> AssociationMap[
    Together[assembledKernelExtraction[#] -
      kernelUnrelatedExtraction[#]] &, Keys[assembledKernelExtraction]]|>];
emitShared["KERNEL_PROPAGATION_RESIDUALS", kernelPropagationResiduals];
emitShared["KERNEL_POLE_LOCATIONS", kernelPoleInventory];
emitShared["CAUSALITY_CHECK", <|
  "ORIENTATION_RECORDS" -> orientationAggregation,
  "PROPAGATION_RECORDS" -> propagationAggregation,
  "TEST_OBJECT" -> causalityAggregate,
  "REMOVED_RECORD_CONTROL" -> removedOrientationAggregation,
  "UNRELATED_OBJECT_CONTROL" -> unrelatedCausalityAggregation,
  "REMOVED_RECORD_INDICATOR_RESIDUAL" ->
    Boole[TrueQ[causalityAggregate]] - Boole[TrueQ[
      removedOrientationAggregation["AGGREGATE_TEST_OBJECT"]]],
  "UNRELATED_OBJECT_INDICATOR_RESIDUAL" ->
    Boole[TrueQ[causalityAggregate]] - Boole[TrueQ[
      unrelatedCausalityAggregation["AGGREGATE_TEST_OBJECT"] &&
      propagationAggregation["AGGREGATE_TEST_OBJECT"]]],
  "COVERAGE_BOUNDARIES" -> {
    Thread[{lambdaA0, lambdaV0, lambdaX0} == 0],
    Thread[{tauA, tauV, tauX} == 0]}|>];

(* ---------------------------------------------------------------------- *)
(* Energy accounting and the two independent power discriminators.        *)
(* ---------------------------------------------------------------------- *)

preEliminationEnergyEOM = fullSystem["PRE_ELIMINATION_EOM"];
slabSpecificChemicalPotential =
  fullSystem["INTERNAL"]["MU_THETA"]/rhoBr;
slabPowerPairingTerms = <|
  "INPLANE_EOM_TIMES_CONJUGATE_U_DOT" -> 1/2 Re[
    preEliminationEnergyEOM["INPLANE_EQUATION"]
      Conjugate[-I omega uL]],
  "THICKNESS_EOM_TIMES_CONJUGATE_DELTA_W_DOT" -> 1/2 Re[
    preEliminationEnergyEOM["THICKNESS_EQUATION"]
      Conjugate[-I omega w0 eAmp]],
  "MASS_BALANCE_TIMES_MU_S" -> 1/2 Re[
    slabSpecificChemicalPotential Conjugate[
      preEliminationEnergyEOM["MASS_EQUATION_WITH_THICKNESS"]]]|>;
slabVirtualPairing = Total[Values[slabPowerPairingTerms]];
slabBoundaryPairing = Expand[slabVirtualPairing -
  (slabVirtualPairing /. {tractionPlusPower -> 0,
    tractionMinusPower -> 0, jPlusPower -> 0, jMinusPower -> 0})];

thicknessVelocityAmplitude = -I omega dW;
faceVelocityScale = Cancel[thicknessVelocityAmplitude/
  Last[First[outwardFaceVelocities]]];
powerFaceRawLaws = rawFaceLawAssociation[rhoBr muSFree,
    velocityFree, lambdaA, lambdaV, lambdaX, 1] /.
  {facePressure -> facePressureFree, faceFlux -> fluxFree,
    faceBulkVelocity -> bulkVelocityFree,
    faceAffinityRaw -> affinityFree,
    faceVelocityRaw -> velocityFree,
    faceTractionRaw -> tractionFree};

slabFacePowerFromEOM[eom_Association, rawLaws_Association] := Module[
  {tractionCoefficient, fluxCoefficient, tractionExpression},
  tractionCoefficient = faceVelocityScale Coefficient[
    eom["THICKNESS_EQUATION"], tractionPlusPower];
  fluxCoefficient = Coefficient[
    eom["MASS_EQUATION_WITH_THICKNESS"], jPlusPower];
  tractionExpression = constitutiveFluxExpression[
    rawLaws["TRACTION_LAW"], tractionFree];
  -1/2 Re[tractionCoefficient tractionExpression
      Conjugate[velocityFree] +
    fluxCoefficient muSFree Conjugate[fluxFree]]
];

routeASlabFacePower = slabFacePowerFromEOM[
  preEliminationEnergyEOM, powerFaceRawLaws];
powerFluxExpression = constitutiveFluxExpression[
  powerFaceRawLaws["CLOSURE_LAW"], fluxFree];
powerBulkVelocityExpression = constitutiveFluxExpression[
  powerFaceRawLaws["KINEMATIC_LAW"], bulkVelocityFree];
powerAffinityExpression = affinityFree;
powerMuSExpression = Cancel[muSFree /.
  First[Solve[powerFaceRawLaws["AFFINITY_DEFINITION"], muSFree]]];
powerTractionExpression = constitutiveFluxExpression[
  powerFaceRawLaws["TRACTION_LAW"], tractionFree];
powerStateRules = {muSFree -> powerMuSExpression};
routeASlabPowerOnFaceLaws = routeASlabFacePower /. powerStateRules;
routeBOutgoingBulkPower = 1/2 Re[
  facePressureFree Conjugate[powerBulkVelocityExpression]];
routeBConstitutivePower = 1/2 Re[
  powerAffinityExpression Conjugate[fluxFree] +
  (powerTractionExpression - facePressureFree)
    Conjugate[velocityFree]];
routeBTotalOutgoingPower = routeBOutgoingBulkPower +
  routeBConstitutivePower;
energyExchangeResidual = FullSimplify[ComplexExpand[
  routeASlabPowerOnFaceLaws + routeBTotalOutgoingPower,
  {facePressureFree, affinityFree, velocityFree, fluxFree}]];
signedExchangeTerms = {-routeBOutgoingBulkPower,
  -(routeBConstitutivePower /. fluxFree -> powerFluxExpression)};
routeBB0bImpedanceOperand = facePressureFree ==
  zImpermeable powerBulkVelocityExpression;

emitShared["ENERGY_SINKS", <|
  "ROUTE_A_SLAB_PAIRING" -> slabPowerPairingTerms,
  "ROUTE_A_BOUNDARY_OPERAND" -> routeASlabPowerOnFaceLaws,
  "ROUTE_B_OUTGOING_BULK_OPERAND" -> routeBOutgoingBulkPower,
  "ROUTE_B_B0B_IMPEDANCE_OPERAND" -> routeBB0bImpedanceOperand,
  "ROUTE_B_CONSTITUTIVE_OPERAND" -> routeBConstitutivePower,
  "SIGNED_EXCHANGE_TERMS" -> signedExchangeTerms,
  "REAL_POWER_SIGNS" -> (Sign[Re[#]] & /@ signedExchangeTerms)|>];
emitShared["ENERGY_SOURCES", <|
  "ROUTE_A_SLAB_PAIRING" -> slabPowerPairingTerms,
  "SIGNED_EXCHANGE_TERMS" -> signedExchangeTerms,
  "NEGATED_REAL_POWER_SIGNS" -> (Sign[-Re[#]] & /@
    signedExchangeTerms)|>];
emitShared["UNATTRIBUTED_SINK_TERMS", <|
  "ROUTE_A_SLAB_OPERAND" -> routeASlabPowerOnFaceLaws,
  "ROUTE_B_BULK_PLUS_CONSTITUTIVE_OPERAND" -> routeBTotalOutgoingPower,
  "RESIDUAL" -> energyExchangeResidual|>];
emitShared["UNATTRIBUTED_EXCHANGE_TERMS", comparisonRecord[
  routeASlabPowerOnFaceLaws, -routeBTotalOutgoingPower]];

upperPressureAmplitude = -rhoM D[upperAcousticAnsatz, t]/harmonicFactor /.
  w -> w0/2;
upperVelocityAmplitude = D[upperAcousticAnsatz, w]/harmonicFactor /.
  w -> w0/2;
lowerPressureAmplitude = -rhoM D[lowerAcousticAnsatz, t]/harmonicFactor /.
  w -> -w0/2;
lowerVelocityAmplitude = -D[lowerAcousticAnsatz, w]/harmonicFactor /.
  w -> -w0/2;
outgoingBulkPower = 1/2 Re[
  upperPressureAmplitude Conjugate[upperVelocityAmplitude] +
  lowerPressureAmplitude Conjugate[lowerVelocityAmplitude]];
routeAPressureFacePower = routeASlabFacePower /.
  {lambdaA0 -> 0, lambdaV0 -> 0, lambdaX0 -> 0, fluxFree -> 0};
slabPressurePower = Total[{
  routeAPressureFacePower /. {facePressureFree -> upperPressureAmplitude,
    velocityFree -> upperVelocityAmplitude},
  routeAPressureFacePower /. {facePressureFree -> lowerPressureAmplitude,
    velocityFree -> lowerVelocityAmplitude}}];
pressurePowerResidual = ComplexExpand[
  slabPressurePower + outgoingBulkPower,
  {aPlus, aMinus}];

routeAEOMAblated = preEliminationEOMAssociation[
  fullSystem["INTERNAL"], 1 + routeAEOMAblation, 1];
routeAEOMAblatedFacePower = slabFacePowerFromEOM[
  routeAEOMAblated, powerFaceRawLaws];
routeAEOMAblatedPressureFacePower = routeAEOMAblatedFacePower /.
  {lambdaA0 -> 0, lambdaV0 -> 0, lambdaX0 -> 0, fluxFree -> 0};
routeAEOMAblatedPressurePower = Total[{
  routeAEOMAblatedPressureFacePower /.
    {facePressureFree -> upperPressureAmplitude,
      velocityFree -> upperVelocityAmplitude},
  routeAEOMAblatedPressureFacePower /.
    {facePressureFree -> lowerPressureAmplitude,
      velocityFree -> lowerVelocityAmplitude}}];
pressurePowerEOMAblationResidual = ComplexExpand[
  routeAEOMAblatedPressurePower + outgoingBulkPower,
  {aPlus, aMinus}];
routeAUnrelatedEOM = preEliminationEOMAssociation[
  deriveInternalOperators[
    fullEnergyDensity /. kappaW -> kappaW + energyUnrelatedAblation],
  1, 1];
routeAUnrelatedPressureFacePower = slabFacePowerFromEOM[
    routeAUnrelatedEOM, powerFaceRawLaws] /.
  {lambdaA0 -> 0, lambdaV0 -> 0, lambdaX0 -> 0, fluxFree -> 0};
routeAUnrelatedPressurePower = Total[{
  routeAUnrelatedPressureFacePower /.
    {facePressureFree -> upperPressureAmplitude,
      velocityFree -> upperVelocityAmplitude},
  routeAUnrelatedPressureFacePower /.
    {facePressureFree -> lowerPressureAmplitude,
      velocityFree -> lowerVelocityAmplitude}}];
pressurePowerUnrelatedResidual = ComplexExpand[
  routeAUnrelatedPressurePower + outgoingBulkPower,
  {aPlus, aMinus}];
emitShared["PRESSURE_WORK_SIGN_CHECK", <|
  "ROUTE_A_PAIRING_OPERANDS" -> slabPowerPairingTerms,
  "ROUTE_A_BOUNDARY_PAIRING" -> slabBoundaryPairing,
  "SLAB_OFF_SHELL_OPERAND" -> slabPressurePower,
  "OUTGOING_BULK_OPERAND" -> outgoingBulkPower,
  "RESIDUAL" -> pressurePowerResidual,
  "TEST_OBJECT" -> casTest[pressurePowerResidual == 0],
  "ONE_SIDED_EOM_ABLATION" -> <|
    "ROUTE_A_OPERAND" -> routeAEOMAblatedPressurePower,
    "ROUTE_B_OPERAND" -> outgoingBulkPower,
    "RESIDUAL" -> pressurePowerEOMAblationResidual|>,
  "UNRELATED_ENERGY_ABLATION" -> <|
    "ROUTE_A_OPERAND" -> routeAUnrelatedPressurePower,
    "ROUTE_B_OPERAND" -> outgoingBulkPower,
    "RESIDUAL" -> pressurePowerUnrelatedResidual|>,
  "RESTRICTIONS" -> {Element[omega, Reals], qSquared > 0,
    lambdaA0 == 0, lambdaV0 == 0, lambdaX0 == 0}|>];

routeAEOMAblatedPowerOnFaceLaws = routeAEOMAblatedFacePower /.
  powerStateRules;
routeAClosureAblatedLaws = rawFaceLawAssociation[rhoBr muSFree,
    velocityFree, lambdaA + routeAClosureAblation, lambdaV, lambdaX, 1] /.
  {facePressure -> facePressureFree, faceFlux -> fluxFree,
    faceBulkVelocity -> bulkVelocityFree,
    faceAffinityRaw -> affinityFree,
    faceVelocityRaw -> velocityFree,
    faceTractionRaw -> tractionFree};
routeAClosureAblatedFlux = constitutiveFluxExpression[
  routeAClosureAblatedLaws["CLOSURE_LAW"], fluxFree];
routeAClosureAblatedPower = slabFacePowerFromEOM[
    preEliminationEnergyEOM, routeAClosureAblatedLaws] /.
  {fluxFree -> routeAClosureAblatedFlux,
    muSFree -> powerMuSExpression};
routeASlabPowerOnClosure = routeASlabPowerOnFaceLaws /.
  fluxFree -> powerFluxExpression;
routeBTotalOutgoingPowerOnClosure = routeBTotalOutgoingPower /.
  fluxFree -> powerFluxExpression;
routeAUnrelatedFacePower = slabFacePowerFromEOM[
    routeAUnrelatedEOM, powerFaceRawLaws] /. powerStateRules;
twoPortExpandedA = ComplexExpand[routeASlabPowerOnFaceLaws,
  {facePressureFree, affinityFree, velocityFree, fluxFree}];
twoPortExpandedB = ComplexExpand[-routeBTotalOutgoingPower,
  {facePressureFree, affinityFree, velocityFree, fluxFree}];
twoPortOrders = Range[0, Exponent[twoPortExpandedA, lambdaX0]];
twoPortOrderRecords = Table[
  With[{operandA = Coefficient[twoPortExpandedA, lambdaX0, order],
      operandB = Coefficient[twoPortExpandedB, lambdaX0, order]},
    order -> comparisonRecord[operandA, operandB]],
  {order, twoPortOrders}];
emitShared["FULL_TWO_PORT_BALANCE_CHECK", <|
  "ROUTE_A_SLAB_OPERAND" -> routeASlabPowerOnFaceLaws,
  "ROUTE_B_OUTGOING_BULK_OPERAND" -> routeBOutgoingBulkPower,
  "ROUTE_B_CONSTITUTIVE_OPERAND" -> routeBConstitutivePower,
  "RESIDUAL" -> energyExchangeResidual,
  "ORDERS_IN_RECIPROCAL_COEFFICIENT" -> twoPortOrderRecords,
  "ONE_SIDED_EOM_ABLATION" -> <|
    "ROUTE_A_OPERAND" -> routeAEOMAblatedPowerOnFaceLaws,
    "ROUTE_B_OPERAND" -> routeBTotalOutgoingPower,
    "RESIDUAL" -> FullSimplify[ComplexExpand[
      routeAEOMAblatedPowerOnFaceLaws + routeBTotalOutgoingPower,
      {facePressureFree, affinityFree, velocityFree, fluxFree}]]|>,
  "ONE_SIDED_CLOSURE_ABLATION" -> <|
    "ROUTE_A_OPERAND" -> routeAClosureAblatedPower,
    "ROUTE_B_OPERAND" -> routeBTotalOutgoingPowerOnClosure,
    "RESIDUAL" -> FullSimplify[ComplexExpand[
      routeAClosureAblatedPower + routeBTotalOutgoingPowerOnClosure,
      {facePressureFree, affinityFree, velocityFree}]]|>,
  "UNRELATED_ENERGY_ABLATION" -> <|
    "ROUTE_A_OPERAND" -> routeAUnrelatedFacePower,
    "ROUTE_B_OPERAND" -> routeBTotalOutgoingPower,
    "RESIDUAL" -> FullSimplify[ComplexExpand[
      routeAUnrelatedFacePower + routeBTotalOutgoingPower,
      {facePressureFree, affinityFree, velocityFree, fluxFree}]]|>|>];

(* ---------------------------------------------------------------------- *)
(* B4: the k=0 impermeable breathing slice with the bulk retained.         *)
(* ---------------------------------------------------------------------- *)

breathingInternalEnergy = energy1D[fullEnergyDensity] /.
  {uL -> 0, uLX -> 0, uT -> 0, uTX -> 0, thetaX -> 0, eX -> 0,
   thetaAmp -> conservedSigma - eAmp};
breathingStiffness = D[breathingInternalEnergy, {eAmp, 2}];
qBreathing = omega/cs0;
zBreathing = Cancel[zImpermeable /. {k -> 0, q -> qBreathing}];
breathingPressure = Cancel[zBreathing (-I omega w0 eAmp/2)];
breathingEquation = Expand[-muW omega^2 w0 eAmp +
  breathingStiffness eAmp/w0 + breathingPressure];
breathingDispersionPolynomial = Cancel[breathingEquation/eAmp];
breathingRoots = Solve[breathingDispersionPolynomial == 0, omega];
breathingRootImaginaryParts = ComplexExpand[
  Im[omega /. breathingRoots],
  {Sqrt[-4 breathingStiffness muW + rhoM^2 cs0^2 w0^2]}];
breathingGrowthCondition = LogicalExpand[Exists[omega,
  breathingDispersionPolynomial == 0 && Im[omega] > 0]];

emitShared["BREATHING_SLICE_DISPERSION", <|
  "BULK_IMPEDANCE" -> zBreathing,
  "EQUATION" -> breathingDispersionPolynomial == 0,
  "ROOTS" -> breathingRoots,
  "IMAGINARY_PARTS" -> breathingRootImaginaryParts|>];
emitShared["BREATHING_STABILITY_CONDITION", <|
  "STIFFNESS_FROM_ENERGY" -> breathingStiffness,
  "ROOT_GROWTH_CONDITION" -> breathingGrowthCondition,
  "MODULUS_DEPENDENCE" ->
    (D[breathingStiffness, #] & /@ {br3, crossC, kW})|>];

(* ---------------------------------------------------------------------- *)
(* B5: longitudinal dispersion, exact implicit roots and sheet audit.      *)
(* ---------------------------------------------------------------------- *)

longitudinalDispersion = fullSystem["DETERMINANT"];
longitudinalPolynomial = fullSystem["DISPERSION_NUMERATOR"];
longitudinalRootSolve = Solve[
  longitudinalPolynomial == 0, omega, Cubics -> False, Quartics -> False];
longitudinalRootValues = omega /. longitudinalRootSolve;

complexRootSubstitution = {omega -> omegaR + I omegaI,
  q -> qR + I qI};
rootEquationComplex = Together[longitudinalPolynomial /.
  complexRootSubstitution];
qEquationComplex = Together[(q^2 - qSquared) /.
  complexRootSubstitution];
rootRealEquations = Thread[ComplexExpand[{
    Re[Numerator[rootEquationComplex]], Im[Numerator[rootEquationComplex]],
    Re[Numerator[qEquationComplex]], Im[Numerator[qEquationComplex]]}] == 0];
rootVariables = {omegaR, omegaI, qR, qI};
rootImplicitRegion = ImplicitRegion[And @@ rootRealEquations, rootVariables];
growthRootRegion = ImplicitRegion[And @@ Join[rootRealEquations,
   {omegaI > 0}], rootVariables];
decayRootRegion = ImplicitRegion[And @@ Join[rootRealEquations,
   {omegaI < 0}], rootVariables];
realRootRegion = ImplicitRegion[And @@ Join[rootRealEquations,
   {omegaI == 0}], rootVariables];
stabilityConditionObject = Exists[rootVariables,
  And @@ Join[rootRealEquations, {omegaI > 0}]];

physicalSheetEquation = q == qContinued[omega];
oppositeSheetEquation = q == -qContinued[omega];
physicalSheetRootSet = ImplicitRegion[
  longitudinalPolynomial == 0 && physicalSheetEquation, {omega, q}];
oppositeSheetRootSet = ImplicitRegion[
  (longitudinalPolynomial /. q -> -q) == 0 && oppositeSheetEquation,
  {omega, q}];
branchSensitivityRatio = ConditionalExpression[
  Im[omegaOpposite]/Im[omegaPhysical],
  (longitudinalPolynomial /. {omega -> omegaPhysical,
      q -> qContinued[omegaPhysical]}) == 0 &&
   (longitudinalPolynomial /. {omega -> omegaOpposite,
      q -> -qContinued[omegaOpposite]}) == 0 &&
   Im[omegaPhysical] != 0];

interfaceZeroSystem = assembleSystem[fullEnergyDensity, 0, 0, 0, 1, True];
impermeableSystem = assembleSystem[fullEnergyDensity, 0, 0, lambdaX, 1, True];
bulkEvanescentSystem = fullSystem["DETERMINANT"] /. q -> qEvanescent;
bulkPropagatingSystem = fullSystem["DETERMINANT"] /. q -> qPropagating;

thresholdMatrixQ0 = Map[
  SeriesCoefficient[Together[#], {q, 0, 0}] &,
  fullSystem["MATRIX"], {2}];
thresholdSoundConeBlocks = <|
  "POSITIVE_SOUND_CONE" -> Simplify[
    thresholdMatrixQ0 /. omega -> cs0 k],
  "NEGATIVE_SOUND_CONE" -> Simplify[
    thresholdMatrixQ0 /. omega -> -cs0 k]|>;
thresholdLeadingOrderAvailable =
  AllTrue[Values[thresholdSoundConeBlocks], MatrixQ] &&
  FreeQ[thresholdSoundConeBlocks, _Missing | _SeriesCoefficient |
    Indeterminate | ComplexInfinity];
modeCategoryByUnknown = AssociationThread[
  fullSystem["UNKNOWNS"],
  {longitudinalModeCategory, thicknessModeCategory,
    densityModeCategory}];
classifyLeadingBlock[block_?MatrixQ] := Module[
  {rank, nullspace, strata},
  rank = MatrixRank[block];
  nullspace = NullSpace[block];
  strata = Map[Function[nullVector, <|
      "LEADING_VECTOR" -> nullVector,
      "COMPONENT_SUPPORT" -> AssociationThread[
        fullSystem["UNKNOWNS"], Not[TrueQ[PossibleZeroQ[#]]] & /@
          nullVector],
      "MODE_CATEGORIES" -> Pick[Values[modeCategoryByUnknown],
        Not[TrueQ[PossibleZeroQ[#]]] & /@ nullVector]|>], nullspace];
  <|"RANK" -> rank, "NULLSPACE" -> nullspace,
    "NULLITY" -> Length[nullspace], "MODE_STRATA" -> strata|>
];
thresholdLeadingClassification = If[TrueQ[thresholdLeadingOrderAvailable],
  AssociationMap[Function[branch, Module[
    {block, dispersion, stratumRule, stratumBlock},
    block = thresholdSoundConeBlocks[branch];
    dispersion = Factor[Together[Det[block]]];
    stratumRule = First[Solve[dispersion == 0, chi5]];
    stratumBlock = Map[Cancel[# /. stratumRule] &, block, {2}];
    <|"LAURENT_Q0_BLOCK" -> block,
      "LEADING_DISPERSION_EQUATION" -> dispersion == 0,
      "STRATUM_PARAMETER_RULE" -> stratumRule,
      "STRATUM_BLOCK" -> stratumBlock,
      "CLASSIFICATION" -> classifyLeadingBlock[stratumBlock]|>
  ]], Keys[thresholdSoundConeBlocks]], Missing["LeadingOrderBlock"]];
thresholdBulkNormFinite = Integrate[Abs[thresholdBulkAmplitude]^2,
  {w, w0/2, radialCutoff}, Assumptions -> radialCutoff > w0/2];
thresholdBulkNorm = Block[{$Assumptions = True},
  Limit[thresholdBulkNormFinite, radialCutoff -> Infinity]];
grazingPressureLaurentCoefficient = SeriesCoefficient[
  zImpermeable velocityFace, {q, 0, -1}];
ablateLeadingStratum[payload_Association] := Module[
  {block, cofactors, position, ablatedBlock, ablatedClassification},
  block = payload["STRATUM_BLOCK"];
  cofactors = Table[(-1)^(row + column) Det[
    Delete[Map[Delete[#, column] &, block], row]],
    {row, Length[block]}, {column, Length[First[block]]}];
  position = First[Position[cofactors,
    entry_ /; !TrueQ[PossibleZeroQ[entry]], {2}, Heads -> False],
    Missing["NonzeroCofactor"]];
  ablatedBlock = ReplacePart[block, position ->
    (Extract[block, position] + grazingBlockAblation)];
  ablatedClassification = classifyLeadingBlock[ablatedBlock];
  <|"POSITION" -> position, "BLOCK" -> ablatedBlock,
    "CLASSIFICATION" -> ablatedClassification,
    "RANK_RESIDUAL" ->
      ablatedClassification["RANK"] -
        payload["CLASSIFICATION"]["RANK"],
    "NULLITY_RESIDUAL" ->
      ablatedClassification["NULLITY"] -
        payload["CLASSIFICATION"]["NULLITY"]|>
];
thresholdAblationPayload = If[TrueQ[thresholdLeadingOrderAvailable],
  AssociationMap[ablateLeadingStratum[
    thresholdLeadingClassification[#]] &,
    Keys[thresholdLeadingClassification]],
  Missing["LeadingOrderBlock"]];
thresholdUnrelatedClassification = thresholdLeadingClassification;
grazingCriteria = If[TrueQ[thresholdLeadingOrderAvailable], <|
    "LAURENT_Q0_BLOCKS" -> thresholdSoundConeBlocks,
    "LEADING_ORDER_STRATA" -> thresholdLeadingClassification,
    "NONDEGENERATE_LEADING_PAYLOAD" -> casTest[
      AllTrue[Values[thresholdLeadingClassification],
        !AllTrue[Flatten[#1["LAURENT_Q0_BLOCK"]],
          TrueQ[PossibleZeroQ[#]] &] &&
        Length[#1["CLASSIFICATION"]["NULLSPACE"]] > 0 &]],
    "MODE_CATEGORY_MAP" -> modeCategoryByUnknown,
    "COALESCENCE_SERIES" ->
      Series[Exp[I q (w - w0/2)], {q, 0, 1}],
    "BULK_NORM" -> thresholdBulkNorm,
    "ZERO_BULK_AMPLITUDE_CONDITION" ->
      Solve[thresholdBulkAmplitude == 0, thresholdBulkAmplitude],
    "FINITE_PRESSURE_CONDITION" -> Solve[
      grazingPressureLaurentCoefficient == 0, velocityFace],
    "LEADING_BLOCK_ABLATION" -> thresholdAblationPayload,
    "UNRELATED_BULK_AMPLITUDE_ABLATION" -> <|
      "OPERAND" -> thresholdBulkNorm /.
        thresholdBulkAmplitude ->
          thresholdBulkAmplitude + grazingUnrelatedAblation,
      "CLASSIFICATION" -> thresholdUnrelatedClassification,
      "RANK_RESIDUALS" -> AssociationMap[
        thresholdUnrelatedClassification[#]["CLASSIFICATION"]["RANK"] -
          thresholdLeadingClassification[#]["CLASSIFICATION"]["RANK"] &,
        Keys[thresholdLeadingClassification]],
      "NULLITY_RESIDUALS" -> AssociationMap[
        thresholdUnrelatedClassification[#]["CLASSIFICATION"]["NULLITY"] -
          thresholdLeadingClassification[#]["CLASSIFICATION"]["NULLITY"] &,
        Keys[thresholdLeadingClassification]]|>|>,
  <|"AVAILABLE_OPERANDS" -> <|
      "MATRIX_PENCIL" -> fullSystem["MATRIX"],
      "LAURENT_Q0_ATTEMPT" -> thresholdMatrixQ0|>,
    "STATUS" -> "MISSING_OPERAND"|>];

closedFormRootTest = casTest[FreeQ[longitudinalRootSolve, _Solve]];
rootStabilityClasses = AssociationThread[
  longitudinalRootValues, Sign[Im[#]] & /@ longitudinalRootValues];

emitShared["LONGITUDINAL_DISPERSION", <|
  "UNKNOWN_ORDER" -> fullSystem["UNKNOWNS"],
  "MATRIX" -> fullSystem["MATRIX"],
  "RATIONAL_DETERMINANT" -> longitudinalDispersion,
  "NUMERATOR_EQUATION" -> longitudinalPolynomial == 0,
  "BULK_EQUATION" -> qAlgebraicEquation|>];
emitShared["ROOTS", <|
  "CAS_SOLUTION" -> longitudinalRootSolve,
  "CLOSED_FORM_TEST_OBJECT" -> closedFormRootTest,
  "EXACT_IMPLICIT_REGION" -> rootImplicitRegion|>];
emitShared["IMAGINARY_PART", <|
  "EXPLICIT_ROOT_OBJECTS" -> (Im /@ longitudinalRootValues),
  "IMPLICIT_COORDINATE" -> omegaI,
  "GROWTH_REGION" -> growthRootRegion,
  "DECAY_REGION" -> decayRootRegion,
  "REAL_REGION" -> realRootRegion|>];
emitShared["DISSIPATION_ORIGIN", <|
  "FULL_DETERMINANT" -> longitudinalDispersion,
  "INTERFACE_COEFFICIENT_CUT" -> interfaceZeroSystem["DETERMINANT"],
  "INTERFACE_DIFFERENCE" -> Together[longitudinalDispersion -
    interfaceZeroSystem["DETERMINANT"]],
  "IMPERMEABLE_RADIATION_OPERAND" -> impermeableSystem["DETERMINANT"],
  "PROPAGATING_BULK_OPERAND" -> bulkPropagatingSystem,
  "EVANESCENT_BULK_OPERAND" -> bulkEvanescentSystem|>];
emitShared["ROOT_STABILITY_CLASS", rootStabilityClasses];
emitShared["STABILITY_CONDITION", <|
  "QUANTIFIED_ROOT_CONDITION" -> stabilityConditionObject,
  "MODULUS_AND_C_DEPENDENCE" ->
    (D[longitudinalPolynomial, #] & /@
      Join[{br3, crossC, kW, kappaW, muR}, extraCoefficients])|>];
emitShared["GRAZING_MODE_CLASSIFICATION", grazingCriteria];
emitShared["BRANCH_SENSITIVITY", <|
  "PHYSICAL_SHEET_ROOT_SET" -> physicalSheetRootSet,
  "OPPOSITE_SHEET_ROOT_SET" -> oppositeSheetRootSet,
  "IMAGINARY_PART_RATIO" -> branchSensitivityRatio|>];
emitShared["SHEET_OF_EACH_ROOT", <|
  "ROOT_OBJECTS" -> longitudinalRootValues,
  "PHYSICAL_SHEET_DEFINITION" -> physicalSheetEquation,
  "OPPOSITE_SHEET_DEFINITION" -> oppositeSheetEquation,
  "UPPER_RIM_CONTINUATION" -> qContinued[omegaR + I omegaI],
  "LOWER_CUT_JUMPS" -> {
    Limit[qContinued[cs0 k + horizontalOffset + I omegaI],
      horizontalOffset -> 0, Direction -> "FromAbove"] -
     Limit[qContinued[cs0 k + horizontalOffset + I omegaI],
      horizontalOffset -> 0, Direction -> "FromBelow"],
    Limit[qContinued[-cs0 k + horizontalOffset + I omegaI],
      horizontalOffset -> 0, Direction -> "FromAbove"] -
     Limit[qContinued[-cs0 k + horizontalOffset + I omegaI],
      horizontalOffset -> 0, Direction -> "FromBelow"]},
  "COMPLEX_ROOT_SPATIAL_RESELECTION_OPERAND" ->
    {physicalSheetEquation, qAlgebraicEquation},
  "CUT_CROSSING_LOCUS" ->
    {reOmega == cs0 k, reOmega == -cs0 k}|>];

rootAdmissibilityObjects = <|
  "ROOT_REGION" -> rootImplicitRegion,
  "PORT_REGION" -> portDissipativityCondition,
  "THERMODYNAMIC_REGION" -> admissibilityCoefficientRegion,
  "RECIPROCITY_REGION" -> onsagerRegion,
  "REGION_SET_RELATION" -> regionSetRelation|>;
emitShared["GROWTH_INSIDE_ADMISSIBLE", Join[rootAdmissibilityObjects,
  <|"ROOT_SIGN_CONDITION" -> omegaI > 0|>]];
emitShared["DECAY_INSIDE_ADMISSIBLE", Join[rootAdmissibilityObjects,
  <|"ROOT_SIGN_CONDITION" -> omegaI < 0|>]];
emitShared["ROOT_POWER_SOURCE", <|
  "ROOT_OUTSIDE_PORT_REGION" -> Not[portDissipativityCondition],
  "SIGNED_EXCHANGE_TERMS" -> signedExchangeTerms,
  "TRANSPORT_RESIDUAL" -> energyExchangeResidual|>];
emitShared["GROWTH_ARTIFACT_DIAGNOSTICS", <|
  "ROOT_REGION" -> growthRootRegion,
  "KERNEL_ORIENTATION" -> kernelOrientationIdentities,
  "KERNEL_PROPAGATION" -> kernelPropagationResiduals,
  "KERNEL_POLES" -> kernelPoleInventory,
  "SHEET" -> physicalSheetEquation,
  "ADMISSIBILITY" -> rootAdmissibilityObjects|>];
emitShared["DECAY_ARTIFACT_DIAGNOSTICS", <|
  "ROOT_REGION" -> decayRootRegion,
  "KERNEL_ORIENTATION" -> kernelOrientationIdentities,
  "KERNEL_PROPAGATION" -> kernelPropagationResiduals,
  "KERNEL_POLES" -> kernelPoleInventory,
  "SHEET" -> physicalSheetEquation,
  "ADMISSIBILITY" -> rootAdmissibilityObjects|>];
emitShared["RECIPROCAL_TRACTION_ROOT_EFFECT", <|
  "FULL_ROOT_SET" -> rootImplicitRegion,
  "FORM_CUT_DETERMINANT" -> lambdaXZeroSystem["DETERMINANT"],
  "DETERMINANT_DIFFERENCE" -> Together[longitudinalDispersion -
    lambdaXZeroSystem["DETERMINANT"]],
  "FORM_CUT_ROOT_SET" -> ImplicitRegion[
    lambdaXZeroSystem["DISPERSION_NUMERATOR"] == 0 &&
      qAlgebraicEquation, {omega, q}]|>];

(* ---------------------------------------------------------------------- *)
(* B6: transverse coupling generated from the same basis and constraint.   *)
(* ---------------------------------------------------------------------- *)

transverseInternal = internalOperators["TRANSVERSE_INTERNAL"];
transverseCouplingCoefficients = coefficientVector[transverseInternal,
  {eAmp, thetaAmp}];
transverseOperator = Coefficient[
  -rhoBr omega^2 uT + transverseInternal, uT];
transverseDispersion = transverseOperator == 0;
transverseRoots = Solve[transverseDispersion, omega];
transverseNormalization = HoldForm[
  transverseInPlaneEquationCoefficientOfThicknessAmplitude];
transverseParameterDependence = D[transverseOperator, #] & /@
  {lambdaA0, lambdaV0, lambdaX0, tauA, tauV, tauX};

emitShared["TRANSVERSE_COUPLING", <|
  "NORMALIZATION" -> transverseNormalization,
  "OUTPUT_EQUATION" -> (-rhoBr omega^2 uT + transverseInternal),
  "INPUT_AMPLITUDES" -> {eAmp, thetaAmp},
  "COEFFICIENTS" -> transverseCouplingCoefficients,
  "NORMALIZATION_DETERMINACY" ->
    casTest[transverseCouplingCoefficients != ConstantArray[0, 2]]|>];
emitShared["TRANSVERSE_DISPERSION", <|
  "EQUATION" -> transverseDispersion,
  "ROOTS" -> transverseRoots|>];
emitShared["TRANSVERSE_DISSIPATION", <|
  "ROOT_IMAGINARY_PARTS" -> (Im[omega /. #] & /@ transverseRoots),
  "INTERFACE_PARAMETER_DEPENDENCE" -> transverseParameterDependence,
  "SLAB_AFFINITY_DEPENDENCE" -> D[transverseOperator, slabAffinityFactor]|>];

(* ---------------------------------------------------------------------- *)
(* B8: source-level form controls.                                        *)
(* ---------------------------------------------------------------------- *)

omegaShifted[w_] := Sech[(w - windowShift)/windowWidth]^2;
omegaCentered[w_] := Sech[w/windowWidth]^2;
controlWSource[current_] := Expand[
  -omegaShifted[ell] current[ell] + omegaShifted[-ell] current[-ell] +
   Inactive[Integrate][D[omegaShifted[w], w] current[w],
     {w, -ell, ell}]];
controlISource[current_] := Expand[
  -omegaCentered[ell + intervalShift] current[ell + intervalShift] +
   omegaCentered[-ell] current[-ell] +
   Inactive[Integrate][D[omegaCentered[w], w] current[w],
     {w, -ell, ell + intervalShift}]];
evenControlCurrent[arg_] := evenProfile[arg^2];
oddControlCurrent[arg_] := arg oddProfile[arg^2];
windowControlObjects = {
  evenCurrent -> controlWSource[evenControlCurrent],
  oddCurrent -> controlWSource[oddControlCurrent]};
intervalControlObjects = {
  evenCurrent -> controlISource[evenControlCurrent],
  oddCurrent -> controlISource[oddControlCurrent]};

emitShared["CONTROL_WINDOW_PARITY", <|
  "FAMILY" -> omegaShifted[w], "INTERVAL" -> {-ell, ell},
  "SOURCES" -> windowControlObjects,
  "PARAMETER_DERIVATIVES" ->
    (D[Last[#], windowShift] & /@ windowControlObjects),
  "IDENTICAL_INDEPENDENCE_TESTS" ->
    (casTest[D[Last[#], windowShift] == 0] & /@ windowControlObjects)|>];
emitShared["CONTROL_INTERVAL_SYMMETRY", <|
  "WINDOW" -> omegaCentered[w],
  "INTERVAL_FAMILY" -> {-ell, ell + intervalShift},
  "SOURCES" -> intervalControlObjects,
  "PARAMETER_DERIVATIVES" ->
    (D[Last[#], intervalShift] & /@ intervalControlObjects),
  "IDENTICAL_INDEPENDENCE_TESTS" ->
    (casTest[D[Last[#], intervalShift] == 0] & /@ intervalControlObjects)|>];

noThicknessEnergy = Expand[fullEnergyDensity /.
  Join[{eField -> 0}, Thread[eGradient -> 0]]];
noGradientEnergy = fullEnergyDensity /. kappaW -> 0;
noCrossEnergy = fullEnergyDensity /. crossC -> 0;
controlNoThicknessSystem = assembleSystem[noThicknessEnergy, lambdaA,
  lambdaV, lambdaX, 1, False];
controlNoGradientSystem = assembleSystem[noGradientEnergy, lambdaA,
  lambdaV, lambdaX, 1, True];
controlImpermeableSystem = assembleSystem[fullEnergyDensity, 0, 0,
  lambdaX, 1, True];
controlNoCrossSystem = assembleSystem[noCrossEnergy, lambdaA, lambdaV,
  lambdaX, 1, True];
controlNoMuSystem = assembleSystem[fullEnergyDensity, lambdaA, lambdaV,
  lambdaX, 0, True];
controlNoReciprocalSystem = assembleSystem[fullEnergyDensity, lambdaA,
  lambdaV, 0, 1, True];

controlNoThicknessCompression = frozenStress;
controlNoGradientCompression = Module[{eqs, sol},
  eqs = {controlNoGradientSystem["THICKNESS_EQUATION"] == 0,
      controlNoGradientSystem["MASS_EQUATION"] == 0} /. strainRule;
  sol = First[Solve[eqs, {eAmp, thetaAmp}]];
  Cancel[(controlNoGradientSystem["INTERNAL"]["INPLANE_INTERNAL"]/(I k)
      /strainLong) /. strainRule /. sol]
];
controlNoCrossCompression = Module[{eqs, sol},
  eqs = {controlNoCrossSystem["THICKNESS_EQUATION"] == 0,
      controlNoCrossSystem["MASS_EQUATION"] == 0} /. strainRule;
  sol = First[Solve[eqs, {eAmp, thetaAmp}]];
  Cancel[(controlNoCrossSystem["INTERNAL"]["INPLANE_INTERNAL"]/(I k)
      /strainLong) /. strainRule /. sol]
];
thicknessResponseForSystem[system_] := Module[{solution},
  solution = First[Solve[system["THICKNESS_EQUATION"] == thicknessDrive,
    eAmp]];
  Cancel[(w0 eAmp/thicknessDrive) /. solution]
];
compressionResponseForSystem[system_] := Module[{eqs, solution, stress},
  eqs = {system["THICKNESS_EQUATION"] == 0,
      system["MASS_EQUATION"] == 0} /. strainRule;
  solution = First[Solve[eqs, {eAmp, thetaAmp}]];
  stress = system["INTERNAL"]["INPLANE_INTERNAL"]/(I k);
  Cancel[(stress/strainLong) /. strainRule /. solution]
];

emitShared["CONTROL_NO_THICKNESS", <|
  "SOURCE_ENERGY" -> noThicknessEnergy,
  "COMPRESSIONAL_RESPONSE" -> controlNoThicknessCompression,
  "LONGITUDINAL_MATRIX" -> controlNoThicknessSystem["MATRIX"],
  "LONGITUDINAL_DISPERSION" -> controlNoThicknessSystem["DETERMINANT"],
  "FULL_DISPERSION_RESIDUAL" -> Together[longitudinalDispersion -
    controlNoThicknessSystem["DETERMINANT"]]|>];
thicknessSimultaneousCuts = {
  eAmp -> 0, velocityInput -> 0, lambdaV velocityInput -> 0,
  facePressurePower velocityPower -> 0,
  lambdaX affinityPower velocityPower -> 0};
emitShared["CONTROL_A_ATTRIBUTION", <|
  "SIMULTANEOUS_SOURCE_CUTS" -> thicknessSimultaneousCuts,
  "CUT_COUNT" -> Length[thicknessSimultaneousCuts],
  "SINGLE_CHANNEL_ATTRIBUTION_RANK" ->
    MatrixRank[IdentityMatrix[Length[thicknessSimultaneousCuts]]],
  "REMAINING_TRANSFER" ->
    Cancel[controlNoThicknessSystem["FACE"]["FLUX"]]|>];
emitShared["CONTROL_NO_GRADIENT_STIFFNESS", <|
  "SOURCE_ENERGY" -> noGradientEnergy,
  "THICKNESS_RESPONSE" ->
    thicknessResponseForSystem[controlNoGradientSystem],
  "COMPRESSIONAL_RESPONSE" -> controlNoGradientCompression,
  "DISPERSION" -> controlNoGradientSystem["DETERMINANT"],
  "DISPERSION_RESIDUAL" -> Together[longitudinalDispersion -
    controlNoGradientSystem["DETERMINANT"]]|>];
emitShared["CONTROL_IMPERMEABLE", <|
  "SOURCE_LAW_COEFFICIENTS" -> {0, 0, lambdaX},
  "DISPERSION" -> controlImpermeableSystem["DETERMINANT"],
  "FULL_DISPERSION_RESIDUAL" -> Together[longitudinalDispersion -
    controlImpermeableSystem["DETERMINANT"]]|>];
emitShared["CONTROL_NO_CROSS_TERM", <|
  "SOURCE_ENERGY" -> noCrossEnergy,
  "COMPRESSIONAL_RESPONSE" -> controlNoCrossCompression,
  "DISPERSION" -> controlNoCrossSystem["DETERMINANT"],
  "DISPERSION_RESIDUAL" -> Together[longitudinalDispersion -
    controlNoCrossSystem["DETERMINANT"]]|>];

noMuPortFace = solveFace[rhoBr muSInput, velocityInput, lambdaA, lambdaV,
  lambdaX, 0];
noMuPortMatrix = coefficientVector[#,
    portInputVector] & /@ {noMuPortFace["TRACTION"],
      noMuPortFace["FLUX"]};
noMuHermitian = Simplify[(noMuPortMatrix +
  ConjugateTranspose[noMuPortMatrix])/2];
emitShared["CONTROL_NO_MU_COUPLING", <|
  "SOURCE_AFFINITY_FACTOR" -> 0,
  "DISPERSION" -> controlNoMuSystem["DETERMINANT"],
  "PORT_HERMITIAN_MATRIX" -> noMuHermitian,
  "THERMODYNAMIC_MATRIX" -> mixedHermitianMatrix,
  "DISPERSION_RESIDUAL" -> Together[longitudinalDispersion -
    controlNoMuSystem["DETERMINANT"]]|>];

noReciprocalPortFace = solveFace[rhoBr muSInput, velocityInput,
  lambdaA, lambdaV, 0, 1];
noReciprocalPortMatrix = coefficientVector[#, portInputVector] & /@
  {noReciprocalPortFace["PRESSURE"], noReciprocalPortFace["FLUX"]};
noReciprocalHermitian = Simplify[(noReciprocalPortMatrix +
  ConjugateTranspose[noReciprocalPortMatrix])/2];
noReciprocalMixedSource = {{lambdaA, lambdaV}, {0, 0}};
noReciprocalMixedHermitian = Simplify[(noReciprocalMixedSource +
  ConjugateTranspose[noReciprocalMixedSource])/2];
noReciprocalPowerIdentity = twoPortLeft /. lambdaX0 -> 0;
emitShared["CONTROL_NO_RECIPROCAL_TRACTION", <|
  "RECOMPUTED_INPLANE" -> controlNoReciprocalSystem["INPLANE_EQUATION"],
  "RECOMPUTED_THICKNESS" -> controlNoReciprocalSystem["THICKNESS_EQUATION"],
  "RECOMPUTED_MASS" -> controlNoReciprocalSystem["MASS_EQUATION"],
  "RECOMPUTED_THICKNESS_RESPONSE" ->
    thicknessResponseForSystem[controlNoReciprocalSystem],
  "RECOMPUTED_COMPRESSIONAL_RESPONSE" ->
    compressionResponseForSystem[controlNoReciprocalSystem],
  "RECOMPUTED_DISPERSION" -> controlNoReciprocalSystem["DETERMINANT"],
  "RECOMPUTED_ROOT_SET" -> ImplicitRegion[
    controlNoReciprocalSystem["DISPERSION_NUMERATOR"] == 0 &&
      qAlgebraicEquation, {omega, q}],
  "MECHANICAL_OPERATOR_COMPARISON" -> comparisonRecord[
    controlNoReciprocalSystem["THICKNESS_EQUATION"],
    fullSystem["THICKNESS_EQUATION"] /. lambdaX0 -> 0],
  "POWER_MATRIX" -> noReciprocalHermitian,
  "THERMODYNAMIC_POWER_MATRIX" -> noReciprocalMixedHermitian,
  "POWER_IDENTITY_COMPARISON" -> comparisonRecord[
    noReciprocalPowerIdentity, twoPortLeft /. lambdaX0 -> 0],
  "POWER_FORM_COMPARISON" -> comparisonRecord[
    noReciprocalMixedHermitian, mixedHermitianMatrix /. lambdaX0 -> 0]|>];

transverseCouplingFromEnergy[sourceEnergy_] := Module[{operator},
  operator = deriveInternalOperators[sourceEnergy]["TRANSVERSE_INTERNAL"];
  coefficientVector[operator, {eAmp, thetaAmp}]
];
transverseControlObjects = <|
  "A" -> transverseCouplingFromEnergy[noThicknessEnergy],
  "B" -> transverseCouplingFromEnergy[noGradientEnergy],
  "C" -> transverseCouplingFromEnergy[fullEnergyDensity],
  "D" -> transverseCouplingFromEnergy[noCrossEnergy]|>;
emitShared["CONTROLS_ON_TRANSVERSE", <|
  "FULL" -> transverseCouplingCoefficients,
  "CONTROLS" -> transverseControlObjects,
  "RESIDUALS" -> AssociationThread[Keys[transverseControlObjects],
    (Together[transverseCouplingCoefficients - #] & /@
      Values[transverseControlObjects])],
  "CHANNEL_PARAMETER_DERIVATIVES" -> transverseParameterDependence|>];

(* ---------------------------------------------------------------------- *)
(* B9: validity and the discarded convective correction.                  *)
(* ---------------------------------------------------------------------- *)

convectiveDispersion = Expand[(omega - vDr q)^2 -
  cs0^2 (k^2 + q^2)];
restDispersion = Expand[omega^2 - cs0^2 (k^2 + q^2)];
discardedConvectiveCorrection = Expand[convectiveDispersion - restDispersion];
leadingConvectiveCorrection = Coefficient[
  Normal[Series[discardedConvectiveCorrection, {vDr, 0, 1}]], vDr] vDr;
relativeConvectiveCorrection = Cancel[
  leadingConvectiveCorrection/omega^2];
timescaleMeasure = Abs[vDr q/omega];
speedMeasure = Abs[vDr/cs0];

emitShared["VALIDITY_CONDITIONS", <|
  "LINEAR_AMPLITUDES" -> {Abs[thetaAmp] < 1, Abs[eAmp] < 1,
    Abs[k uL] < 1},
  "ACOUSTIC_BRANCH_DOMAIN" -> qAlgebraicEquation,
  "INDEPENDENT_KERNEL_ARGUMENTS" ->
    {Abs[omega tauA], Abs[omega tauV], Abs[omega tauX]},
  "KERNEL_REGIMES" -> Table[
    {Abs[omega tau] < 1, Abs[omega tau] == 1,
      Abs[omega tau] > 1}, {tau, {tauA, tauV, tauX}}]|>];
emitShared["VALIDITY_TIMESCALE", <|
  "MEASURE" -> timescaleMeasure,
  "SMALLNESS_CONDITION" -> timescaleMeasure < 1,
  "COMPLEX_MEASURE" -> ComplexExpand[Abs[vDr (qR + I qI)/
      (omegaR + I omegaI)]]|>];
emitShared["VALIDITY_FLOW_SPEED", <|
  "MEASURE" -> speedMeasure,
  "SMALLNESS_CONDITION" -> speedMeasure < 1|>];
emitShared["DISCARDED_CONVECTIVE_CORRECTION", <|
  "CONVECTIVE_OPERATOR" -> convectiveDispersion,
  "REST_OPERATOR" -> restDispersion,
  "DIFFERENCE" -> discardedConvectiveCorrection,
  "LEADING_TERM" -> leadingConvectiveCorrection,
  "RELATIVE_TERM" -> relativeConvectiveCorrection|>];
emitShared["VALIDITY_FAILURE_REGION", <|
  "TIMESCALE_REGION" -> timescaleMeasure >= 1,
  "SPEED_REGION" -> speedMeasure >= 1,
  "COMPLEX_NORM_DEFINED_WHERE" ->
    Reduce[omegaR^2 + omegaI^2 > 0, {omegaR, omegaI}, Reals]|>];

(* ---------------------------------------------------------------------- *)
(* B7: dimensions derived from the supplied equations and energy terms.   *)
(* ---------------------------------------------------------------------- *)

dimL = {1, 0, 0};
dimT = {0, 1, 0};
dimM = {0, 0, 1};
dimZero = {0, 0, 0};
dimEnergyPerBraneVolume = dimM - dimL - 2 dimT;
dimBulkMassDensity = dimM - 4 dimL;
dimBraneMassDensity = dimM - 3 dimL;
dimVelocity = dimL - dimT;
dimAcceleration = dimL - 2 dimT;
dimBulkPressure = dimM - 2 dimL - 2 dimT;
dimMassFlux = dimBulkMassDensity + dimVelocity;
dimAffinity = dimBulkPressure - dimBulkMassDensity;

basisAtomDimensions = AssociationThread[basisVariables,
  Join[ConstantArray[dimZero, Length[gVariables] + 2],
    ConstantArray[-dimL, Length[thetaGradient] + Length[eGradient]]]];
basisTermDimension[polynomial_] := Module[{rules, exponents},
  rules = CoefficientRules[polynomial, basisVariables];
  exponents = First[First[rules]];
  Total[MapThread[#1 #2 &, {exponents, Values[basisAtomDimensions]}]]
];
extraCoefficientDimensions = dimEnergyPerBraneVolume -
    basisTermDimension[#] & /@ omittedBasis;

transverseZeroTest = casTest[And @@ Thread[
   transverseCouplingCoefficients == 0]];
transverseCouplingDimension = If[TrueQ[transverseZeroTest],
  Missing["Undetermined", transverseNormalization],
  dimBulkPressure];

dimensionObjects = <|
  "Z" -> (dimBulkPressure - dimVelocity),
  "M_ADD" -> (dimBulkPressure - dimAcceleration),
  "RHO_M" -> dimBulkMassDensity,
  "C_S0" -> dimVelocity,
  "V_DR" -> dimVelocity,
  "RHO_BR0" -> dimBraneMassDensity,
  "B_RHO" -> (dimEnergyPerBraneVolume - dimL),
  "B_RHO_3" -> dimEnergyPerBraneVolume,
  "MU_W" -> (dimEnergyPerBraneVolume - 2 (dimL - dimT)),
  "K_W" -> (dimEnergyPerBraneVolume - 2 dimL),
  "KAPPA_W" -> (dimEnergyPerBraneVolume - 2 dimL),
  "C" -> (dimEnergyPerBraneVolume - dimL),
  "MU_R" -> dimEnergyPerBraneVolume,
  "THICKNESS_RESPONSE" -> (dimL - dimBulkPressure),
  "COMPRESSIONAL_RESPONSE" -> dimEnergyPerBraneVolume,
  "TRANSVERSE_COUPLING" -> transverseCouplingDimension,
  "LAMBDA_A0" -> (dimMassFlux - dimAffinity),
  "LAMBDA_V0" -> (dimMassFlux - dimVelocity),
  "LAMBDA_X0" -> (dimBulkPressure - dimAffinity),
  "TAU_A" -> dimT, "TAU_V" -> dimT, "TAU_X" -> dimT,
  "AFFINITY" -> dimAffinity,
  "MU_THETA" -> dimEnergyPerBraneVolume,
  "MU_S" -> (dimEnergyPerBraneVolume - dimBraneMassDensity),
  "PROJECTION_SOURCE" -> (dimBraneMassDensity - dimT),
  "FACE_RESPONSE_VELOCITY_COEFFICIENT" ->
    (dimBulkPressure - dimVelocity),
  "FACE_RESPONSE_MU_COEFFICIENT" ->
    (dimBulkPressure - dimEnergyPerBraneVolume)
|>;
Do[AssociateTo[dimensionObjects,
  "CHI" <> ToString[index] -> extraCoefficientDimensions[[index]]],
  {index, Length[extraCoefficientDimensions]}];

dimensionRoutes = <|
  "Z" -> HoldForm[pressureDimension - velocityDimension],
  "M_ADD" -> HoldForm[pressureDimension - accelerationDimension],
  "RHO_M" -> HoldForm[massDimension - 4 lengthDimension],
  "C_S0" -> HoldForm[frequencyDimension - wavenumberDimension],
  "V_DR" -> HoldForm[normalDistanceDimension/timeDimension],
  "RHO_BR0" -> HoldForm[rho4Dimension + w0Dimension],
  "B_RHO" -> HoldForm[braneEnergyDensityDimension - w0Dimension],
  "B_RHO_3" -> HoldForm[thetaQuadraticEnergyDimension],
  "MU_W" -> HoldForm[braneEnergyDensityDimension -
     2 thicknessVelocityDimension],
  "K_W" -> HoldForm[braneEnergyDensityDimension - 2 w0Dimension],
  "KAPPA_W" -> HoldForm[braneEnergyDensityDimension - 2 w0Dimension -
     2 gradientThicknessDimension],
  "C" -> HoldForm[braneEnergyDensityDimension - w0Dimension],
  "MU_R" -> HoldForm[braneEnergyDensityDimension -
     2 curlDisplacementDimension],
  "THICKNESS_RESPONSE" -> HoldForm[
     thicknessDisplacementDimension - generalizedPressureDimension],
  "COMPRESSIONAL_RESPONSE" -> HoldForm[
     longitudinalStressDimension - longitudinalStrainDimension],
  "TRANSVERSE_COUPLING" -> HoldForm[
     transverseEquationDimension - thicknessAmplitudeDimension],
  "LAMBDA_A0" -> HoldForm[massFluxDimension - affinityDimension],
  "LAMBDA_V0" -> HoldForm[massFluxDimension - faceVelocityDimension],
  "LAMBDA_X0" -> HoldForm[tractionDimension - affinityDimension],
  "TAU_A" -> HoldForm[omegaDimension tauADimension == 0],
  "TAU_V" -> HoldForm[omegaDimension tauVDimension == 0],
  "TAU_X" -> HoldForm[omegaDimension tauXDimension == 0],
  "AFFINITY" -> HoldForm[muSDimension == pressureDimension -
     rhoMDimension],
  "MU_THETA" -> HoldForm[functionalDerivativeOfEnergyWithRespectToTheta],
  "MU_S" -> HoldForm[muThetaDimension - rhoBrDimension],
  "PROJECTION_SOURCE" -> HoldForm[integratedMassDensityDimension -
     timeDimension],
  "FACE_RESPONSE_VELOCITY_COEFFICIENT" -> HoldForm[
     facePressureDimension - faceVelocityDimension],
  "FACE_RESPONSE_MU_COEFFICIENT" -> HoldForm[
     facePressureDimension - muThetaDimension]
|>;
Do[AssociateTo[dimensionRoutes, "CHI" <> ToString[index] ->
   HoldForm[braneEnergyDensityDimension -
     omittedInvariantDimension[index]]],
  {index, Length[extraCoefficientDimensions]}];

dimensionRouteKinds = AssociationMap[
  If[MemberQ[{"Z", "M_ADD", "RHO_M", "RHO_BR0",
      "THICKNESS_RESPONSE", "COMPRESSIONAL_RESPONSE",
      "TRANSVERSE_COUPLING"}, #], "DEFINITIONAL", "INDEPENDENT"] &,
  Keys[dimensionObjects]];

Do[
  emitShared["DIM_" <> name, <|"VECTOR_L_T_M" -> dimensionObjects[name],
    "ROUTE" -> dimensionRoutes[name]|>];
  emitShared["DIM_ROUTE_KIND_" <> name, dimensionRouteKinds[name]],
  {name, Keys[dimensionObjects]}];

atomDimensions = <|
  omega -> -dimT, k -> -dimL, q -> -dimL, w0 -> dimL,
  rhoM -> dimBulkMassDensity, rhoBr -> dimBraneMassDensity,
  rho4 -> dimM - 4 dimL, cs0 -> dimVelocity, vDr -> dimVelocity,
  uL -> dimL, uT -> dimL, eAmp -> dimZero, thetaAmp -> dimZero,
  muW -> dimensionObjects["MU_W"], muR -> dimensionObjects["MU_R"],
  br3 -> dimensionObjects["B_RHO_3"], bRho -> dimensionObjects["B_RHO"],
  crossC -> dimensionObjects["C"], kW -> dimensionObjects["K_W"],
  kappaW -> dimensionObjects["KAPPA_W"],
  lambdaA0 -> dimensionObjects["LAMBDA_A0"],
  lambdaV0 -> dimensionObjects["LAMBDA_V0"],
  lambdaX0 -> dimensionObjects["LAMBDA_X0"],
  tauA -> dimT, tauV -> dimT, tauX -> dimT,
  ellA -> dimensionObjects["LAMBDA_A0"],
  ellV -> dimensionObjects["LAMBDA_V0"],
  ellX -> dimensionObjects["LAMBDA_X0"],
  facePressure -> dimBulkPressure, facePressurePower -> dimBulkPressure,
  facePressureFree -> dimBulkPressure, deltaPExtract -> dimBulkPressure,
  faceFlux -> dimMassFlux, fluxPower -> dimMassFlux,
  fluxFree -> dimMassFlux, velocityInput -> dimVelocity,
  velocityFace -> dimVelocity, velocityPower -> dimVelocity,
  velocityFree -> dimVelocity, faceBulkVelocity -> dimVelocity,
  muThetaFace -> dimEnergyPerBraneVolume,
  muThetaPower -> dimEnergyPerBraneVolume,
  muSInput -> dimAffinity, muSPower -> dimAffinity, muSFree -> dimAffinity,
  affinityExtract -> dimAffinity, affinityPower -> dimAffinity,
  affinityFree -> dimAffinity, tractionExtract -> dimBulkPressure,
  pressureAmplitude -> dimBulkPressure,
  outwardBulkVelocity -> dimVelocity,
  thicknessDrive -> dimBulkPressure, strainLong -> dimZero,
  raySlope -> -dimVelocity
|>;
Do[AssociateTo[atomDimensions, extraCoefficients[[index]] ->
   extraCoefficientDimensions[[index]]],
  {index, Length[extraCoefficients]}];

Clear[dimensionOfExpression];
dimensionOfExpression[expression_?NumericQ] := dimZero;
dimensionOfExpression[s_Symbol] := Lookup[atomDimensions, s,
  Missing["DimensionNotEstablished", s]];
dimensionOfExpression[expression_Plus] := Module[{dimensions},
  dimensions = dimensionOfExpression /@ (List @@ expression);
  If[SameQ @@ dimensions, First[dimensions],
    Missing["Inhomogeneous", dimensions]]
];
dimensionOfExpression[expression_Times] := Module[{dimensions},
  dimensions = dimensionOfExpression /@ (List @@ expression);
  If[AnyTrue[dimensions, MissingQ],
    Missing["FactorDimensionNotEstablished", dimensions], Total[dimensions]]
];
dimensionOfExpression[expression_Power] := Module[{base, exponent},
  base = First[expression]; exponent = Last[expression];
  If[NumericQ[exponent], exponent dimensionOfExpression[base],
    Missing["ExponentDimensionNotEstablished", HoldForm[expression]]]
];
dimensionOfExpression[Conjugate[argument_]] := dimensionOfExpression[argument];
dimensionOfExpression[Re[argument_]] := dimensionOfExpression[argument];
dimensionOfExpression[Im[argument_]] := dimensionOfExpression[argument];
dimensionOfExpression[Abs[argument_]] := dimensionOfExpression[argument];
dimensionOfExpression[Sign[argument_]] := dimZero;
dimensionOfExpression[expression_] := Missing[
  "DimensionNotEstablished", HoldForm[expression]];

homogeneityRecord[equationName_, expression_] := Module[
  {terms, dimensions, testObject},
  terms = If[Head[Expand[expression]] === Plus,
    List @@ Expand[expression], {Expand[expression]}];
  dimensions = dimensionOfExpression /@ terms;
  testObject = casTest[SameQ @@ dimensions && !AnyTrue[dimensions, MissingQ]];
  <|"EQUATION" -> equationName, "TERMS" -> terms,
    "TERM_DIMENSIONS" -> dimensions, "TEST_OBJECT" -> testObject|>
];

homogeneityInputs = <|
  "INPLANE_EOM" -> fullSystem["INPLANE_EQUATION"],
  "THICKNESS_EOM" -> fullSystem["THICKNESS_EQUATION"],
  "MASS_BALANCE" -> fullSystem["MASS_EQUATION"],
  "AFFINITY" -> (muThetaFace/rhoBr - facePressure/rhoM),
  "CLOSURE" -> equationZeroForm[generalFace["EQUATIONS"][[3]]],
  "FACE_RESPONSE" ->
    (facePressure - generalFace["PRESSURE"]),
  "TWO_PORT_POWER_IDENTITY" ->
    (twoPortLeft - twoPortRight /. affinityPower ->
       affinityNormalizationPower),
  "DISPERSION_DETERMINANT" -> longitudinalDispersion
|>;
homogeneityRecords = AssociationMap[
  homogeneityRecord[#, homogeneityInputs[#]] &, Keys[homogeneityInputs]];
Do[emitShared["HOMOGENEITY_" <> name, homogeneityRecords[name]],
  {name, Keys[homogeneityRecords]}];

ablationBaseDimensions = homogeneityRecords["MASS_BALANCE"][
  "TERM_DIMENSIONS"];
ablationCorruptedDimensions = ReplacePart[ablationBaseDimensions,
  1 -> (First[ablationBaseDimensions] + dimL)];
ablationCorruptTest = casTest[SameQ @@ ablationCorruptedDimensions];
ablationRestoreTest = casTest[SameQ @@ ablationBaseDimensions];
emitShared["HOMOGENEITY_ABLATION_DEMO", <|
  "SOURCE_EQUATION" -> homogeneityInputs["MASS_BALANCE"],
  "UNMODIFIED_DIMENSIONS" -> ablationBaseDimensions,
  "CORRUPTED_DIMENSIONS" -> ablationCorruptedDimensions,
  "CORRUPTED_TEST_OBJECT" -> ablationCorruptTest,
  "RESTORED_DIMENSIONS" -> ablationBaseDimensions,
  "RESTORED_TEST_OBJECT" -> ablationRestoreTest|>];

localInventoryName = standardEmissionName[localObject["TAG_NAMES"]];
emitLocal["TAG_NAMES", Append[localNames, localInventoryName]];
