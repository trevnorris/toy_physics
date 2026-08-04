$HistoryLength = 0;
ClearAll["Global`*"];

(* Supplied premise P1: the medium has GNLS-superfluid substructure. *)
(* Supplied premise P2: a scalar one-component superfluid has no spin-1 transverse excitation. *)
(* Supplied premise P3: the brane is a symbolic D-dimensional sheet in a higher-dimensional bulk. *)
(* Supplied premise P4: the bulk carries no shear and is absent from the action. *)
(* Supplied premise P5: the brane stiffness has the MacCullagh curl-only form. *)
(* Supplied premise P6: response is linear about an unstrained brane at rest, without dissipation. *)
(* The wave-vector domain supplied to the mode tests is k.k != 0. *)

emittedTags = {};
emit[tag_, value_] := Module[{},
  If[MemberQ[emittedTags, tag], Exit[91]];
  AppendTo[emittedTags, tag];
  Print[tag <> ": " <> ToString[value, InputForm]]
];

fieldHeads = {u1, u2, u3};
coordinates = {t, x, y, z};
spaceCoordinates = Rest[coordinates];
fieldVector = Through[fieldHeads[t, x, y, z]];
velocityVector = D[fieldVector, t];
curlVector = Curl[fieldVector, spaceCoordinates];
divergenceScalar = Div[fieldVector, spaceCoordinates];
gradientMatrix = Table[D[fieldVector[[j]], spaceCoordinates[[i]]], {i, 3}, {j, 3}];

(* Every physical parameter combination below constructs an action. *)
mainLagrangian = rhoBr velocityVector.velocityVector/2 - muR curlVector.curlVector/2;
x1Lagrangian = lambdaRho rhoBr velocityVector.velocityVector/2 - muR curlVector.curlVector/2;
x2Lagrangian = rhoBr velocityVector.velocityVector/2 - lambdaMu muR curlVector.curlVector/2;
x3Lagrangian = rhoBr velocityVector.velocityVector/2 - muR divergenceScalar^2/2;
x4Lagrangian = rhoBr velocityVector.velocityVector/2 - muR Total[Flatten[gradientMatrix]^2]/2;
x5Lagrangian = rhoBr velocityVector.velocityVector/2 + muR curlVector.curlVector/2;
x6Lagrangian = velocityVector.DiagonalMatrix[{rhoBr, rhoBr, rhoZ}].velocityVector/2 -
  muR curlVector.curlVector/2;

(* The only other hand construction involving wave symbols is the supplied ansatz. *)
amplitudeVector = {a1, a2, a3};
waveVector = {kx, ky, kz};
wavePhase = Exp[I (waveVector.spaceCoordinates - omega t)];
planeWaveRules = MapThread[
  #1 -> Function[{t, x, y, z}, Evaluate[#2 wavePhase]] &,
  {fieldHeads, amplitudeVector}
];
planeWaveAnsatzLeft = fieldVector;
planeWaveAnsatzRight = fieldVector /. planeWaveRules;
planeWaveAnsatzResidual = planeWaveAnsatzLeft - planeWaveAnsatzRight;
planeWaveAnsatz = Thread[planeWaveAnsatzLeft == planeWaveAnsatzRight];
waveScalar = waveVector.waveVector;
transverseGenerator = Expand[waveScalar IdentityMatrix[Length[waveVector]] -
  Outer[Times, waveVector, waveVector]];
longitudinalGenerator = Outer[Times, waveVector, waveVector];

evenPolynomial[polynomial_, variable_, squaredVariable_] := Module[
  {expanded = Expand[polynomial], degree},
  degree = Exponent[expanded, variable];
  Sum[Coefficient[expanded, variable, 2 index] squaredVariable^index,
    {index, 0, Floor[degree/2]}]
];

rootMultisetFromPolynomial[polynomial_] := Module[
  {factorPairs, repeatedRoots},
  factorPairs = Select[Rest[FactorList[Factor[polynomial]]],
    Not[FreeQ[First[#], omegaSquared]] &];
  repeatedRoots = Map[
    Function[pair,
      Flatten[ConstantArray[
        omegaSquared /. Solve[pair[[1]] == 0, omegaSquared], pair[[2]]], 1]
    ],
    factorPairs
  ];
  Flatten[repeatedRoots, 1]
];

rootDiagnostic[matrixInOmegaSquared_, root_] := Module[
  {matrixAtRoot, stackedMatrix, stackedRank, dimensionOperand,
   e1Residual, e1Test, e2Count, e3Operand, e3Reference, e3Residual,
   e3Test, e4Operand, e4Reference, e4Residual, e4Test},
  matrixAtRoot = FullSimplify[matrixInOmegaSquared /. omegaSquared -> root];
  stackedMatrix = ArrayFlatten[{{matrixAtRoot},
    {Transpose[Transpose[{waveVector}]]}}];
  stackedRank = MatrixRank[stackedMatrix];
  dimensionOperand = Length[waveVector];
  e1Residual = FullSimplify[stackedRank - dimensionOperand];
  e1Test = FullSimplify[stackedRank < dimensionOperand];
  e2Count = FullSimplify[dimensionOperand - stackedRank];
  e3Operand = FullSimplify[matrixAtRoot.waveVector];
  e3Reference = ConstantArray[0, Length[waveVector]];
  e3Residual = FullSimplify[e3Operand - e3Reference];
  e3Test = SameQ[e3Residual, e3Reference];
  e4Operand = FullSimplify[matrixAtRoot.transverseGenerator];
  e4Reference = ConstantArray[0, Dimensions[e4Operand]];
  e4Residual = FullSimplify[e4Operand - e4Reference];
  e4Test = SameQ[e4Residual, e4Reference];
  <|
    "Root" -> root, "Matrix" -> matrixAtRoot, "Stack" -> stackedMatrix,
    "E1Left" -> stackedRank, "E1Right" -> dimensionOperand,
    "E1Residual" -> e1Residual, "E1Test" -> e1Test, "E2" -> e2Count,
    "E3Left" -> e3Operand, "E3Right" -> e3Reference,
    "E3Residual" -> e3Residual, "E3Test" -> e3Test,
    "E4Left" -> e4Operand, "E4Right" -> e4Reference,
    "E4Residual" -> e4Residual, "E4Test" -> e4Test
  |>
];

selectRootsByTest[roots_, diagnostics_, key_] := Select[
  roots,
  Function[currentRoot,
    TrueQ[diagnostics[[First[FirstPosition[roots, currentRoot]], key]]]
  ]
];

deriveFromAction[lagrangian_] := Module[
  {eulerResidual, equationReference, equationDifference, equations,
   substitutedResidual, planeResidual,
   dynamicalMatrix, determinant, determinantOmegaSquared, rootRules,
   roots, matrixOmegaSquared, diagnostics, transverseRoots,
   longitudinalRoots, transverseSubspaceRoots, rootMultiset, rootNullities},
  eulerResidual = Table[
    D[lagrangian, fieldVector[[component]]] -
      Total[Table[
        D[D[lagrangian, D[fieldVector[[component]], coordinate]], coordinate],
        {coordinate, coordinates}
      ]],
    {component, Length[fieldVector]}
  ];
  equationReference = ConstantArray[0, Length[fieldVector]];
  equationDifference = eulerResidual - equationReference;
  equations = Thread[eulerResidual == equationReference];
  substitutedResidual = eulerResidual /. planeWaveRules;
  planeResidual = FullSimplify[substitutedResidual/wavePhase];
  dynamicalMatrix = Table[
    Coefficient[Expand[planeResidual[[row]]], amplitudeVector[[column]]],
    {row, Length[fieldVector]}, {column, Length[amplitudeVector]}
  ];
  determinant = Factor[Det[dynamicalMatrix]];
  determinantOmegaSquared = Factor[
    evenPolynomial[determinant, omega, omegaSquared]
  ];
  rootRules = DeleteDuplicates[Solve[determinantOmegaSquared == 0, omegaSquared]];
  roots = omegaSquared /. rootRules;
  matrixOmegaSquared = Map[
    evenPolynomial[#, omega, omegaSquared] &,
    dynamicalMatrix,
    {2}
  ];
  diagnostics = rootDiagnostic[matrixOmegaSquared, #] & /@ roots;
  transverseRoots = selectRootsByTest[roots, diagnostics, "E1Test"];
  longitudinalRoots = selectRootsByTest[roots, diagnostics, "E3Test"];
  transverseSubspaceRoots = selectRootsByTest[roots, diagnostics, "E4Test"];
  rootMultiset = rootMultisetFromPolynomial[determinantOmegaSquared];
  rootNullities = ({#, Length[waveVector] - MatrixRank[
      FullSimplify[matrixOmegaSquared /. omegaSquared -> #]]} &) /@ rootMultiset;
  <|
    "Lagrangian" -> lagrangian, "EulerResidual" -> eulerResidual,
    "EquationReference" -> equationReference, "EquationDifference" -> equationDifference,
    "Equations" -> equations, "SubstitutedResidual" -> substitutedResidual,
    "PlaneResidual" -> planeResidual, "Matrix" -> dynamicalMatrix,
    "Determinant" -> determinant, "DeterminantOmegaSquared" -> determinantOmegaSquared,
    "RootRules" -> rootRules, "Roots" -> roots,
    "MatrixOmegaSquared" -> matrixOmegaSquared, "Diagnostics" -> diagnostics,
    "TransverseRoots" -> transverseRoots, "LongitudinalRoots" -> longitudinalRoots,
    "TransverseSubspaceRoots" -> transverseSubspaceRoots,
    "RootMultiset" -> rootMultiset, "RootNullities" -> rootNullities
  |>
];

emitRootDiagnostics[prefix_, diagnostics_] := Module[{index, diagnostic, base},
  Do[
    diagnostic = diagnostics[[index]];
    base = prefix <> "ROOT" <> ToString[index];
    emit[base, diagnostic["Root"]];
    emit[base <> "_E1_STACK", diagnostic["Stack"]];
    emit[base <> "_E1_LEFT", diagnostic["E1Left"]];
    emit[base <> "_E1_RIGHT", diagnostic["E1Right"]];
    emit[base <> "_E1_RESIDUAL", diagnostic["E1Residual"]];
    emit[base <> "_E1_TEST", diagnostic["E1Test"]];
    emit[base <> "_E2", diagnostic["E2"]];
    emit[base <> "_E3_LEFT", diagnostic["E3Left"]];
    emit[base <> "_E3_RIGHT", diagnostic["E3Right"]];
    emit[base <> "_E3_RESIDUAL", diagnostic["E3Residual"]];
    emit[base <> "_E3_TEST", diagnostic["E3Test"]];
    emit[base <> "_E4_LEFT", diagnostic["E4Left"]];
    emit[base <> "_E4_RIGHT", diagnostic["E4Right"]];
    emit[base <> "_E4_RESIDUAL", diagnostic["E4Residual"]];
    emit[base <> "_E4_TEST", diagnostic["E4Test"]],
    {index, Length[diagnostics]}
  ]
];

mainData = deriveFromAction[mainLagrangian];

emit["WL_S9_LAGRANGIAN", mainData["Lagrangian"]];
emit["WL_S9_EULER_LAGRANGE_RESIDUAL", mainData["EulerResidual"]];
emit["WL_S9_EQUATION_REFERENCE", mainData["EquationReference"]];
emit["WL_S9_EQUATION_DIFFERENCE", mainData["EquationDifference"]];
emit["WL_S9_EQUATIONS", mainData["Equations"]];
emit["WL_S9_PLANE_WAVE_ANSATZ_LEFT", planeWaveAnsatzLeft];
emit["WL_S9_PLANE_WAVE_ANSATZ_RIGHT", planeWaveAnsatzRight];
emit["WL_S9_PLANE_WAVE_ANSATZ_RESIDUAL", planeWaveAnsatzResidual];
emit["WL_S9_PLANE_WAVE_ANSATZ", planeWaveAnsatz];
emit["WL_S9_PLANE_WAVE_SUBSTITUTED_RESIDUAL", mainData["SubstitutedResidual"]];
emit["WL_S9_PLANE_WAVE_RESIDUAL", mainData["PlaneResidual"]];
emit["WL_S9_DYNAMICAL_MATRIX", mainData["Matrix"]];
emit["WL_S9_DETERMINANT", mainData["Determinant"]];
emit["WL_S9_OMEGA_SQUARED_SOLUTIONS", mainData["RootRules"]];
emit["WL_S9_TRANSVERSE_GENERATOR", transverseGenerator];
emit["WL_S9_LONGITUDINAL_GENERATOR", longitudinalGenerator];
emitRootDiagnostics["WL_S9_", mainData["Diagnostics"]];
emit["WL_S9_E1_ROOT_SUBSET", mainData["TransverseRoots"]];
emit["WL_S9_E3_ROOT_SUBSET", mainData["LongitudinalRoots"]];
emit["WL_S9_E4_ROOT_SUBSET", mainData["TransverseSubspaceRoots"]];
emit["WL_S9_ROOT_MULTISET", mainData["RootMultiset"]];
emit["WL_S9_ROOT_NULLITIES", mainData["RootNullities"]];

toWaveScalar[expression_] := FullSimplify[Factor[expression] /. waveScalar -> q];
transverseRootsInQ = toWaveScalar /@ mainData["TransverseRoots"];
candidateSpeedsSquared = FullSimplify[#/q] & /@ transverseRootsInQ;
homogeneityDefects = FullSimplify[q D[#, q] - #] & /@ transverseRootsInQ;
speedVariations = FullSimplify[D[#, q]] & /@ candidateSpeedsSquared;
numeratorDegrees = Exponent[Numerator[Together[#]], q] & /@ transverseRootsInQ;
denominatorDegrees = Exponent[Denominator[Together[#]], q] & /@ transverseRootsInQ;

Do[
  emit["WL_S9_E1_ROOT_Q" <> ToString[index], transverseRootsInQ[[index]]];
  emit["WL_S9_CANDIDATE_SPEED_SQUARED" <> ToString[index], candidateSpeedsSquared[[index]]];
  emit["WL_S9_HOMOGENEITY_DEFECT" <> ToString[index], homogeneityDefects[[index]]];
  emit["WL_S9_SPEED_VARIATION" <> ToString[index], speedVariations[[index]]];
  emit["WL_S9_NUMERATOR_Q_DEGREE" <> ToString[index], numeratorDegrees[[index]]];
  emit["WL_S9_DENOMINATOR_Q_DEGREE" <> ToString[index], denominatorDegrees[[index]]],
  {index, Length[transverseRootsInQ]}
];

lengthDimension = {1, 0, 0};
timeDimension = {0, 1, 0};
massDimension = {0, 0, 1};
accelerationDimension = lengthDimension - 2 timeDimension;
forceDimension = massDimension + accelerationDimension;
energyDimension = forceDimension + lengthDimension;
energyDensityDimension = energyDimension - D lengthDimension;
displacementDimension = lengthDimension;
timeDerivativeDimension = displacementDimension - timeDimension;
curlDimension = displacementDimension - lengthDimension;
rhoUnknownDimension = Array[rhoExponent, 3];
muUnknownDimension = Array[muExponent, 3];
dimensionEquationLeft = Join[
  rhoUnknownDimension + 2 timeDerivativeDimension,
  muUnknownDimension + 2 curlDimension
];
dimensionEquationRight = Flatten[ConstantArray[energyDensityDimension, 2], 1];
dimensionEquationResidual = dimensionEquationLeft - dimensionEquationRight;
dimensionEquations = Thread[dimensionEquationLeft == dimensionEquationRight];
dimensionSolution = First[Solve[dimensionEquations,
  Join[rhoUnknownDimension, muUnknownDimension]]];
rhoSolvedDimension = FullSimplify[rhoUnknownDimension /. dimensionSolution];
muSolvedDimension = FullSimplify[muUnknownDimension /. dimensionSolution];
modulusDensityDimensionDifference = FullSimplify[muSolvedDimension - rhoSolvedDimension];
parameterPower[expression_, parameter_] :=
  Exponent[Numerator[Together[expression]], parameter] -
  Exponent[Denominator[Together[expression]], parameter];
speedSquaredParameterPowers = ({parameterPower[#, rhoBr],
    parameterPower[#, muR]} &) /@ candidateSpeedsSquared;
speedSquaredImpliedDimensions = FullSimplify[
  #[[1]] rhoSolvedDimension + #[[2]] muSolvedDimension] & /@
  speedSquaredParameterPowers;
velocitySquaredDimension = 2 (lengthDimension - timeDimension);
speedSquaredDimensionDifferences = FullSimplify[# - velocitySquaredDimension] & /@
  speedSquaredImpliedDimensions;

emit["WL_S9_ENERGY_DENSITY_DIMENSION", energyDensityDimension];
emit["WL_S9_DIMENSION_EQUATION_LEFT", dimensionEquationLeft];
emit["WL_S9_DIMENSION_EQUATION_RIGHT", dimensionEquationRight];
emit["WL_S9_DIMENSION_EQUATION_RESIDUAL", dimensionEquationResidual];
emit["WL_S9_DIMENSION_EQUATIONS", dimensionEquations];
emit["WL_S9_RHO_DIMENSION", rhoSolvedDimension];
emit["WL_S9_MU_DIMENSION", muSolvedDimension];
emit["WL_S9_MU_RHO_DIMENSION_DIFFERENCE", modulusDensityDimensionDifference];
emit["WL_S9_SPEED_SQUARED_PARAMETER_POWERS", speedSquaredParameterPowers];
emit["WL_S9_SPEED_SQUARED_IMPLIED_DIMENSION", speedSquaredImpliedDimensions];
emit["WL_S9_VELOCITY_SQUARED_DIMENSION", velocitySquaredDimension];
emit["WL_S9_SPEED_SQUARED_DIMENSION_DIFFERENCE", speedSquaredDimensionDifferences];
emit["WL_S9_RHO_DIMENSION_D3", rhoSolvedDimension /. D -> 3];
emit["WL_S9_MU_DIMENSION_D3", muSolvedDimension /. D -> 3];
emit["WL_S9_MU_RHO_DIMENSION_DIFFERENCE_D3",
  modulusDensityDimensionDifference /. D -> 3];

emitControl[identifier_, lagrangian_] := Module[
  {data, prefix, multisetDifference},
  data = deriveFromAction[lagrangian];
  prefix = "WL_S9_" <> identifier <> "_";
  multisetDifference = FullSimplify[data["RootMultiset"] - mainData["RootMultiset"]];
  emit[prefix <> "LAGRANGIAN", data["Lagrangian"]];
  emit[prefix <> "DETERMINANT", data["Determinant"]];
  emit[prefix <> "ROOT_MULTISET", data["RootMultiset"]];
  emit[prefix <> "ROOT_NULLITIES", data["RootNullities"]];
  emitRootDiagnostics[prefix, data["Diagnostics"]];
  emit[prefix <> "E1_ROOT_SUBSET", data["TransverseRoots"]];
  emit[prefix <> "E3_ROOT_SUBSET", data["LongitudinalRoots"]];
  emit[prefix <> "E4_ROOT_SUBSET", data["TransverseSubspaceRoots"]];
  emit[prefix <> "ROOT_MULTISET_LEFT", data["RootMultiset"]];
  emit[prefix <> "ROOT_MULTISET_RIGHT", mainData["RootMultiset"]];
  emit[prefix <> "ROOT_MULTISET_RESIDUAL", multisetDifference];
  data
];

x1Data = emitControl["X1", x1Lagrangian];
x2Data = emitControl["X2", x2Lagrangian];
x3Data = emitControl["X3", x3Lagrangian];
x4Data = emitControl["X4", x4Lagrangian];
x5Data = emitControl["X5", x5Lagrangian];
x6Data = emitControl["X6", x6Lagrangian];

x4TransverseRootOperand = First[x4Data["TransverseRoots"]];
x4NonTransverseRootOperand = First[x4Data["LongitudinalRoots"]];
x4RootDifference = FullSimplify[x4TransverseRootOperand - x4NonTransverseRootOperand];
emit["WL_S9_X4_TRANSVERSE_ROOT_OPERAND", x4TransverseRootOperand];
emit["WL_S9_X4_E3_ROOT_OPERAND", x4NonTransverseRootOperand];
emit["WL_S9_X4_ROOT_DIFFERENCE", x4RootDifference];
emit["WL_S9_X6_FULL_ROOT_MULTISET", x6Data["RootMultiset"]];

Exit[0];
