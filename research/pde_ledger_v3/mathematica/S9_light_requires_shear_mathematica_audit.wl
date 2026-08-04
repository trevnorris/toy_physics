$HistoryLength = 0;
ClearAll["Global`*"];

(* Supplied premise P1: the medium has GNLS-superfluid substructure. *)
(* Supplied premise P2: a scalar one-component superfluid has no spin-1 transverse excitation. *)
(* Supplied premise P3: the brane is a symbolic D-dimensional sheet in a higher-dimensional bulk. *)
(* Supplied premise P4: the bulk carries no shear and is absent from the action. *)
(* Supplied premise P5: the brane stiffness has the MacCullagh curl-only form. *)
(* Supplied premise P6: response is linear about an unstrained brane at rest, without dissipation. *)
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
laplacianVector = Laplacian[fieldVector, spaceCoordinates];

(* Every physical parameter combination below constructs an action. *)
mainLagrangian = rhoBr velocityVector.velocityVector/2 - muR curlVector.curlVector/2;
x1Lagrangian = lambdaRho rhoBr velocityVector.velocityVector/2 - muR curlVector.curlVector/2;
x2Lagrangian = rhoBr velocityVector.velocityVector/2 - lambdaMu muR curlVector.curlVector/2;
x3Lagrangian = rhoBr velocityVector.velocityVector/2 - muR divergenceScalar^2/2;
x4Lagrangian = rhoBr velocityVector.velocityVector/2 - muR Total[Flatten[gradientMatrix]^2]/2;
x5Lagrangian = rhoBr velocityVector.velocityVector/2 + muR curlVector.curlVector/2;
x6Lagrangian = velocityVector.DiagonalMatrix[{rhoBr, rhoBr, rhoZ}].velocityVector/2 -
  muR curlVector.curlVector/2;
x7Lagrangian = rhoBr velocityVector.velocityVector/2 -
  muF laplacianVector.laplacianVector/2;
x8Lagrangian = rhoBr velocityVector.velocityVector/2 -
  muR curlVector.curlVector/2 - muG fieldVector.fieldVector/2;

(* The only other hand construction involving wave symbols is the supplied ansatz. *)
amplitudeVector = {a1, a2, a3};
pairedAmplitudeVector = {b1, b2, b3};
waveVector = {kx, ky, kz};
wavePhase = Exp[I (waveVector.spaceCoordinates - omega t)];
pairedWavePhase = Exp[-I (waveVector.spaceCoordinates - omega t)];
planeWaveRules = MapThread[
  #1 -> Function[{t, x, y, z}, Evaluate[#2 wavePhase]] &,
  {fieldHeads, amplitudeVector}
];
planeWaveAnsatz = Thread[fieldVector == (fieldVector /. planeWaveRules)];
quadraticPlaneWaveRules = MapThread[
  #1 -> Function[{t, x, y, z}, Evaluate[#2 wavePhase + #3 pairedWavePhase]] &,
  {fieldHeads, amplitudeVector, pairedAmplitudeVector}
];
quadraticPlaneWaveAnsatz =
  Thread[fieldVector == (fieldVector /. quadraticPlaneWaveRules)];
waveCovector = FullSimplify[
  (-I D[wavePhase, #]/wavePhase) & /@ spaceCoordinates];
waveScalar = Expand[waveCovector.waveCovector];
longitudinalGenerator = Outer[Times, waveCovector, waveCovector];
transverseGenerator = Expand[
  waveScalar IdentityMatrix[Length[waveCovector]] - longitudinalGenerator];

positiveAssumptionQuantities =
  {rhoBr, muR, muF, muG, rhoZ, lambdaRho, lambdaMu, waveScalar};
realWavevectorComponents = waveCovector;
positiveAssumptions = Thread[positiveAssumptionQuantities > 0];
realWavevectorAssumptions = Element[#, Reals] & /@ realWavevectorComponents;
assumptionList = Join[positiveAssumptions, realWavevectorAssumptions];
assumptionSet = And @@ assumptionList;

casSimplify[expression_, assumptions_: assumptionSet] :=
  FullSimplify[expression, Assumptions -> assumptions];

evenPolynomial[polynomial_, variable_, squaredVariable_] := Module[
  {expanded = Expand[polynomial], degree},
  degree = Exponent[expanded, variable];
  Sum[Coefficient[expanded, variable, 2 index] squaredVariable^index,
    {index, 0, Floor[degree/2]}]
];

discardedOddPart[polynomial_, variable_, squaredVariable_] :=
  casSimplify[Expand[polynomial -
    (evenPolynomial[polynomial, variable, squaredVariable] /.
      squaredVariable -> variable^2)]];

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

rootDiagnostic[matrixInOmegaSquared_, root_, diagnosticWaveVector_,
    diagnosticTransverseGenerator_, diagnosticAssumptions_] := Module[
  {omegaSquaredSign, matrixAtRoot, stackedMatrix, stackedRank, dimensionOperand,
   e1Residual, e1Test, e2Count, e3Operand, e3Reference, e3Residual,
   e3Test, e4Operand, e4Reference, e4Residual, e4Test},
  omegaSquaredSign = casSimplify[Sign[root], diagnosticAssumptions];
  matrixAtRoot = casSimplify[
    matrixInOmegaSquared /. omegaSquared -> root, diagnosticAssumptions];
  stackedMatrix = ArrayFlatten[{{matrixAtRoot},
    {Transpose[Transpose[{diagnosticWaveVector}]]}}];
  stackedRank = Assuming[diagnosticAssumptions,
    MatrixRank[casSimplify[stackedMatrix, diagnosticAssumptions]]];
  dimensionOperand = Length[diagnosticWaveVector];
  e1Residual = casSimplify[stackedRank - dimensionOperand,
    diagnosticAssumptions];
  e1Test = casSimplify[stackedRank < dimensionOperand,
    diagnosticAssumptions];
  e2Count = casSimplify[dimensionOperand - stackedRank,
    diagnosticAssumptions];
  e3Operand = casSimplify[matrixAtRoot.diagnosticWaveVector,
    diagnosticAssumptions];
  e3Reference = ConstantArray[0, Length[diagnosticWaveVector]];
  e3Residual = casSimplify[e3Operand - e3Reference, diagnosticAssumptions];
  e3Test = SameQ[e3Residual, e3Reference];
  e4Operand = casSimplify[matrixAtRoot.diagnosticTransverseGenerator,
    diagnosticAssumptions];
  e4Reference = ConstantArray[0, Dimensions[e4Operand]];
  e4Residual = casSimplify[e4Operand - e4Reference, diagnosticAssumptions];
  e4Test = SameQ[e4Residual, e4Reference];
  <|
    "Root" -> root, "OmegaSquaredSign" -> omegaSquaredSign,
    "Matrix" -> matrixAtRoot, "Stack" -> stackedMatrix,
    "E1Left" -> stackedRank, "E1Right" -> dimensionOperand,
    "E1Residual" -> e1Residual, "E1Test" -> e1Test, "E2" -> e2Count,
    "E3Left" -> e3Operand, "E3Right" -> e3Reference,
    "E3Residual" -> e3Residual, "E3Test" -> e3Test,
    "E4Left" -> e4Operand, "E4Right" -> e4Reference,
    "E4Residual" -> e4Residual, "E4Test" -> e4Test
  |>
];

selectRootsByTest[roots_, diagnostics_, key_] := MapThread[
  If[TrueQ[#2[key]], #1, Nothing] &,
  {roots, diagnostics}
];

derivativeAtomPattern = HoldPattern[
  Derivative[dt_, dx_, dy_, dz_][field_][t, x, y, z]
];

derivativeOrder[atom_] := atom /. derivativeAtomPattern :> {dt, dx, dy, dz};

fieldDerivativeAtoms[expression_, fieldHead_] := DeleteDuplicates[Cases[
  expression,
  atom : HoldPattern[Derivative[_, _, _, _][fieldHead][t, x, y, z]] :> atom,
  Infinity
]];

applyDerivativeMultiOrder[expression_, orders_] := Fold[
  D[#1, #2] &,
  expression,
  Flatten[MapThread[ConstantArray, {coordinates, orders}]]
];

eulerResidualFromAction[lagrangian_] := Map[
  Function[fieldHead,
    Module[{bareField, derivativeAtoms, variationalAtoms},
      bareField = fieldHead[t, x, y, z];
      derivativeAtoms = fieldDerivativeAtoms[lagrangian, fieldHead];
      variationalAtoms = Join[
        If[FreeQ[lagrangian, bareField], {}, {{{0, 0, 0, 0}, bareField}}],
        ({derivativeOrder[#], #} & /@ derivativeAtoms)
      ];
      casSimplify[Total[
        (-1)^Total[#[[1]]] applyDerivativeMultiOrder[
          D[lagrangian, #[[2]]], #[[1]]
        ] & /@ variationalAtoms
      ]]
    ]
  ],
  fieldHeads
];

deriveFromAction[lagrangian_] := Module[
  {eulerResidual, equations, substitutedResidual, planeResidual,
   dynamicalMatrixA, pairedPlaneWaveLagrangian, pairedAmplitudeKernel,
   planeWaveQuadratic, dynamicalMatrixB,
   dynamicalMatrixResidual, determinant, determinantOddPart,
   determinantOmegaSquared, rootRules, roots, matrixOddPart,
   matrixOmegaSquared, diagnostics, transverseRoots, longitudinalRoots,
   transverseSubspaceRoots, rootMultiset, rootNullities},
  eulerResidual = eulerResidualFromAction[lagrangian];
  equations = Thread[eulerResidual == ConstantArray[0, Length[fieldVector]]];
  substitutedResidual = eulerResidual /. planeWaveRules;
  planeResidual = casSimplify[substitutedResidual/wavePhase];
  dynamicalMatrixA = casSimplify[Table[
    Coefficient[Expand[planeResidual[[row]]], amplitudeVector[[column]]],
    {row, Length[fieldVector]}, {column, Length[amplitudeVector]}
  ]];
  pairedPlaneWaveLagrangian = Expand[lagrangian /. quadraticPlaneWaveRules];
  pairedAmplitudeKernel = casSimplify[Table[
    D[pairedPlaneWaveLagrangian, amplitudeVector[[row]],
      pairedAmplitudeVector[[column]]],
    {row, Length[amplitudeVector]}, {column, Length[pairedAmplitudeVector]}
  ]];
  planeWaveQuadratic = casSimplify[
    amplitudeVector.pairedAmplitudeKernel.amplitudeVector/2];
  dynamicalMatrixB = casSimplify[Table[
    D[planeWaveQuadratic, amplitudeVector[[row]], amplitudeVector[[column]]],
    {row, Length[amplitudeVector]}, {column, Length[amplitudeVector]}
  ]];
  dynamicalMatrixResidual = casSimplify[dynamicalMatrixA - dynamicalMatrixB];
  determinant = Factor[casSimplify[Det[dynamicalMatrixA]]];
  determinantOddPart = discardedOddPart[determinant, omega, omegaSquared];
  determinantOmegaSquared = Factor[
    evenPolynomial[determinant, omega, omegaSquared]
  ];
  rootRules = DeleteDuplicates[
    Solve[determinantOmegaSquared == 0, omegaSquared]];
  roots = casSimplify[omegaSquared /. rootRules];
  matrixOddPart = Map[
    discardedOddPart[#, omega, omegaSquared] &,
    dynamicalMatrixA,
    {2}
  ];
  matrixOmegaSquared = Map[
    evenPolynomial[#, omega, omegaSquared] &,
    dynamicalMatrixA,
    {2}
  ];
  diagnostics = rootDiagnostic[matrixOmegaSquared, #, waveCovector,
    transverseGenerator, assumptionSet] & /@ roots;
  transverseRoots = selectRootsByTest[roots, diagnostics, "E1Test"];
  longitudinalRoots = selectRootsByTest[roots, diagnostics, "E3Test"];
  transverseSubspaceRoots = selectRootsByTest[roots, diagnostics, "E4Test"];
  rootMultiset = rootMultisetFromPolynomial[determinantOmegaSquared];
  rootNullities = ({#, Length[waveVector] - Assuming[assumptionSet,
      MatrixRank[casSimplify[matrixOmegaSquared /. omegaSquared -> #]]]} &) /@
    rootMultiset;
  <|
    "Lagrangian" -> lagrangian, "EulerResidual" -> eulerResidual,
    "Equations" -> equations, "SubstitutedResidual" -> substitutedResidual,
    "PlaneResidual" -> planeResidual, "PlaneWaveQuadratic" -> planeWaveQuadratic,
    "MatrixA" -> dynamicalMatrixA, "MatrixB" -> dynamicalMatrixB,
    "MatrixResidual" -> dynamicalMatrixResidual,
    "Determinant" -> determinant, "DeterminantOddPart" -> determinantOddPart,
    "DeterminantOmegaSquared" -> determinantOmegaSquared,
    "RootRules" -> rootRules, "Roots" -> roots,
    "MatrixOddPart" -> matrixOddPart,
    "MatrixOmegaSquared" -> matrixOmegaSquared, "Diagnostics" -> diagnostics,
    "TransverseRoots" -> transverseRoots, "LongitudinalRoots" -> longitudinalRoots,
    "TransverseSubspaceRoots" -> transverseSubspaceRoots,
    "RootMultiset" -> rootMultiset, "RootNullities" -> rootNullities
  |>
];

emptyComputedObject[objects_] := Select[objects, False &];

emitRootDiagnostics[prefix_, diagnostics_, slotCount_] := Module[
  {index, diagnosticSlot, diagnostic, base, keys},
  keys = {
    {"", "Root"}, {"_OMEGA_SQUARED_SIGN", "OmegaSquaredSign"},
    {"_E1_STACK", "Stack"}, {"_E1_LEFT", "E1Left"},
    {"_E1_RIGHT", "E1Right"}, {"_E1_RESIDUAL", "E1Residual"},
    {"_E1_TEST", "E1Test"}, {"_E2", "E2"},
    {"_E3_LEFT", "E3Left"}, {"_E3_RIGHT", "E3Right"},
    {"_E3_RESIDUAL", "E3Residual"}, {"_E3_TEST", "E3Test"},
    {"_E4_LEFT", "E4Left"}, {"_E4_RIGHT", "E4Right"},
    {"_E4_RESIDUAL", "E4Residual"}, {"_E4_TEST", "E4Test"}
  };
  Do[
    diagnosticSlot = Pick[diagnostics, Range[Length[diagnostics]], index];
    base = prefix <> "ROOT" <> ToString[index];
    If[diagnosticSlot === {},
      Do[emit[base <> key[[1]], emptyComputedObject[diagnosticSlot]], {key, keys}],
      diagnostic = First[diagnosticSlot];
      Do[emit[base <> key[[1]], diagnostic[key[[2]]]], {key, keys}]
    ],
    {index, slotCount}
  ]
];

lengthDimension = {1, 0, 0};
timeDimension = {0, 1, 0};
massDimension = {0, 0, 1};
accelerationDimension = lengthDimension - 2 timeDimension;
forceDimension = massDimension + accelerationDimension;
energyDimension = forceDimension + lengthDimension;
energyDensityDimension = energyDimension - D lengthDimension;
displacementDimension = lengthDimension;

physicalCoefficientSymbols = {rhoBr, rhoZ, muR, muF, muG};
coefficientUnknownDimensions = <|
  rhoBr -> Array[rhoExponent, 3],
  rhoZ -> Array[rhoZExponent, 3],
  muR -> Array[muExponent, 3],
  muF -> Array[muFExponent, 3],
  muG -> Array[muGExponent, 3]
|>;

parameterPower[expression_, parameter_] :=
  Exponent[Numerator[Together[expression]], parameter] -
  Exponent[Denominator[Together[expression]], parameter];

lagrangianTerms[lagrangian_] := Module[{expanded = Expand[lagrangian]},
  If[Head[expanded] === Plus, List @@ expanded, {expanded}]
];

fieldAtomsForHead[expression_, fieldHead_] := Module[{bareField},
  bareField = fieldHead[t, x, y, z];
  Join[
    If[FreeQ[expression, bareField], {}, {bareField}],
    fieldDerivativeAtoms[expression, fieldHead]
  ]
];

fieldAtoms[expression_] := DeleteDuplicates[Flatten[
  fieldAtomsForHead[expression, #] & /@ fieldHeads
]];

fieldMultiOrder[atom_] := If[
  MatchQ[atom, derivativeAtomPattern], derivativeOrder[atom],
  ConstantArray[0, Length[coordinates]]
];

termFieldOrders[term_] := Module[{atoms},
  atoms = fieldAtoms[term];
  Flatten[
    ConstantArray[fieldMultiOrder[#], Exponent[term, #]] & /@ atoms,
    1
  ]
];

termCoefficient[term_] := Module[{atoms},
  atoms = fieldAtoms[term];
  casSimplify[term/Times @@ (#^Exponent[term, #] & /@ atoms)]
];

dimensionWalk[expression_, symbolDimensionMap_] := Module[
  {walk, result},
  walk[node_] := Module[
    {parts, dimensions, unknowns, mismatches, unsupported, reference,
     differing},
    Which[
      NumberQ[node],
        <|"Dimension" -> ConstantArray[0, 3], "UnknownSymbols" -> {},
          "SumMismatches" -> {}, "Unsupported" -> {}|>,
      Head[node] === Symbol,
        If[KeyExistsQ[symbolDimensionMap, node],
          <|"Dimension" -> symbolDimensionMap[node], "UnknownSymbols" -> {},
            "SumMismatches" -> {}, "Unsupported" -> {}|>,
          <|"Dimension" -> Missing["UnknownSymbol", node],
            "UnknownSymbols" -> {node}, "SumMismatches" -> {},
            "Unsupported" -> {}|>
        ],
      Head[node] === Plus,
        parts = walk /@ (List @@ node);
        dimensions = Lookup[parts, "Dimension"];
        unknowns = DeleteDuplicates[Flatten[Lookup[parts, "UnknownSymbols"]]];
        mismatches = Flatten[Lookup[parts, "SumMismatches"], 1];
        unsupported = Flatten[Lookup[parts, "Unsupported"]];
        If[unknowns === {} && unsupported === {} &&
            FreeQ[dimensions, _Missing],
          reference = First[dimensions];
          differing = Select[Range[Length[dimensions]],
            !SameQ[casSimplify[dimensions[[#]] - reference],
              ConstantArray[0, 3]] &];
          If[differing =!= {},
            AppendTo[mismatches, <|"Summands" -> List @@ node,
              "Dimensions" -> dimensions|>]
          ]
        ];
        <|"Dimension" -> If[mismatches === {} && unknowns === {} &&
              unsupported === {}, First[dimensions], Missing["DimensionFailure"]],
          "UnknownSymbols" -> unknowns, "SumMismatches" -> mismatches,
          "Unsupported" -> unsupported|>,
      Head[node] === Times,
        parts = walk /@ (List @@ node);
        dimensions = Lookup[parts, "Dimension"];
        unknowns = DeleteDuplicates[Flatten[Lookup[parts, "UnknownSymbols"]]];
        mismatches = Flatten[Lookup[parts, "SumMismatches"], 1];
        unsupported = Flatten[Lookup[parts, "Unsupported"]];
        <|"Dimension" -> If[unknowns === {} && mismatches === {} &&
              unsupported === {} && FreeQ[dimensions, _Missing],
            casSimplify[Total[dimensions]], Missing["DimensionFailure"]],
          "UnknownSymbols" -> unknowns, "SumMismatches" -> mismatches,
          "Unsupported" -> unsupported|>,
      Head[node] === Power && MatchQ[node[[2]], _Integer | _Rational],
        result = walk[node[[1]]];
        <|"Dimension" -> If[FreeQ[result["Dimension"], _Missing],
              casSimplify[node[[2]] result["Dimension"]],
              Missing["DimensionFailure"]],
          "UnknownSymbols" -> result["UnknownSymbols"],
          "SumMismatches" -> result["SumMismatches"],
          "Unsupported" -> result["Unsupported"]|>,
      True,
        <|"Dimension" -> Missing["UnsupportedExpression", node],
          "UnknownSymbols" -> DeleteDuplicates[Cases[node, _Symbol, Infinity]],
          "SumMismatches" -> {}, "Unsupported" -> {node}|>
    ]
  ];
  walk[expression]
];

deriveDimensionBlock[lagrangian_, candidateSpeedsSquared_] := Module[
  {terms, derivativeOrders, coefficients, presentCoefficients,
   coefficientDimensions, termDimensions, equationLeft, equationRight,
   equationResidual, equations, unknowns, solutionSet, solution,
   solvedDimension, rhoDimension, rhoZDimension, muDimension,
   muFDimension, muGDimension, modulusDensityDifference, parameterPowers,
   symbolDimensionMap, dimensionWalks, impliedDimensions,
   dimensionUnknownSymbols, dimensionSumMismatches,
   dimensionUnsupportedExpressions, velocitySquaredDimension,
   dimensionDifferences},
  terms = lagrangianTerms[lagrangian];
  derivativeOrders = termFieldOrders /@ terms;
  coefficients = termCoefficient /@ terms;
  presentCoefficients = Select[physicalCoefficientSymbols,
    Function[coefficientSymbol,
      AnyTrue[coefficients,
        Function[coefficient, !FreeQ[coefficient, coefficientSymbol]]]
    ]
  ];
  coefficientDimensions = AssociationMap[
    Function[coefficientSymbol,
      coefficientUnknownDimensions[coefficientSymbol]],
    presentCoefficients
  ];
  termDimensions = MapThread[
    Function[{coefficient, orders},
      casSimplify[
        Total[(parameterPower[coefficient, #] coefficientDimensions[#]) & /@
          presentCoefficients] +
        Total[(displacementDimension - #[[1]] timeDimension -
          Total[Rest[#]] lengthDimension) & /@ orders]
      ]
    ],
    {coefficients, derivativeOrders}
  ];
  equationLeft = Flatten[termDimensions];
  equationRight = Flatten[ConstantArray[energyDensityDimension, Length[terms]], 1];
  equationResidual = casSimplify[equationLeft - equationRight];
  equations = Thread[equationLeft == equationRight];
  unknowns = Flatten[Values[coefficientDimensions]];
  solutionSet = Solve[equations, unknowns];
  If[solutionSet === {},
    Return[<|
      "DerivativeOrders" -> derivativeOrders,
      "TermDimensions" -> termDimensions,
      "EnergyDensityDimension" -> energyDensityDimension,
      "EquationLeft" -> equationLeft, "EquationRight" -> equationRight,
      "EquationResidual" -> equationResidual, "Equations" -> equations,
      "SolutionSet" -> solutionSet
    |>]
  ];
  solution = First[solutionSet];
  solvedDimension[coefficientSymbol_] := If[
    MemberQ[presentCoefficients, coefficientSymbol],
    casSimplify[coefficientDimensions[coefficientSymbol] /. solution],
    emptyComputedObject[presentCoefficients]
  ];
  rhoDimension = solvedDimension[rhoBr];
  rhoZDimension = solvedDimension[rhoZ];
  muDimension = solvedDimension[muR];
  muFDimension = solvedDimension[muF];
  muGDimension = solvedDimension[muG];
  modulusDensityDifference = If[
    Length[rhoDimension] == 3 && Length[muDimension] == 3,
    casSimplify[muDimension - rhoDimension],
    emptyComputedObject[presentCoefficients]
  ];
  parameterPowers = (If[TrueQ[# == 0],
      emptyComputedObject[physicalCoefficientSymbols],
      {parameterPower[#, rhoBr], parameterPower[#, rhoZ],
        parameterPower[#, muR], parameterPower[#, muF],
        parameterPower[#, muG]}
    ] &) /@ candidateSpeedsSquared;
  symbolDimensionMap = Join[
    AssociationMap[
      casSimplify[coefficientDimensions[#] /. solution] &,
      presentCoefficients
    ],
    <|kx -> -lengthDimension, ky -> -lengthDimension,
      kz -> -lengthDimension, q -> -2 lengthDimension,
      omega -> -timeDimension, lambdaRho -> ConstantArray[0, 3],
      lambdaMu -> ConstantArray[0, 3], D -> ConstantArray[0, 3]|>
  ];
  dimensionWalks = dimensionWalk[#, symbolDimensionMap] & /@
    candidateSpeedsSquared;
  impliedDimensions = Lookup[dimensionWalks, "Dimension"];
  dimensionUnknownSymbols = Lookup[dimensionWalks, "UnknownSymbols"];
  dimensionSumMismatches = Lookup[dimensionWalks, "SumMismatches"];
  dimensionUnsupportedExpressions = Lookup[dimensionWalks, "Unsupported"];
  velocitySquaredDimension = 2 (lengthDimension - timeDimension);
  dimensionDifferences = Map[
    If[Length[#] == 3, casSimplify[# - velocitySquaredDimension],
      emptyComputedObject[#]] &,
    impliedDimensions
  ];
  <|
    "DerivativeOrders" -> derivativeOrders,
    "TermDimensions" -> termDimensions,
    "EnergyDensityDimension" -> energyDensityDimension,
    "EquationLeft" -> equationLeft, "EquationRight" -> equationRight,
    "EquationResidual" -> equationResidual, "Equations" -> equations,
    "SolutionSet" -> solutionSet, "RhoDimension" -> rhoDimension,
    "RhoZDimension" -> rhoZDimension, "MuDimension" -> muDimension,
    "MuFDimension" -> muFDimension, "MuGDimension" -> muGDimension,
    "MuRhoDifference" -> modulusDensityDifference,
    "ParameterPowers" -> parameterPowers,
    "SymbolDimensionMap" -> symbolDimensionMap,
    "ImpliedDimensions" -> impliedDimensions,
    "DimensionUnknownSymbols" -> dimensionUnknownSymbols,
    "DimensionSumMismatches" -> dimensionSumMismatches,
    "DimensionUnsupportedExpressions" -> dimensionUnsupportedExpressions,
    "VelocitySquaredDimension" -> velocitySquaredDimension,
    "DimensionDifferences" -> dimensionDifferences
  |>
];

toWaveScalar[expression_] := casSimplify[Factor[expression] /. waveScalar -> q,
  assumptionSet /. waveScalar -> q];

deriveSpeedBlock[data_] := Module[
  {rootsInQ, wavevectorOccurrences, candidateSpeedsSquared,
   homogeneityDefects, speedVariations, numeratorDegrees, denominatorDegrees},
  rootsInQ = toWaveScalar /@ data["TransverseRoots"];
  wavevectorOccurrences = Function[rootInQ,
    Thread[waveCovector ->
      (Function[component, !FreeQ[rootInQ, component]] /@ waveCovector)]
  ] /@ rootsInQ;
  candidateSpeedsSquared = casSimplify[#/q,
      assumptionSet /. waveScalar -> q] & /@ rootsInQ;
  homogeneityDefects = casSimplify[q D[#, q] - #,
      assumptionSet /. waveScalar -> q] & /@ rootsInQ;
  speedVariations = casSimplify[D[#, q],
      assumptionSet /. waveScalar -> q] & /@ candidateSpeedsSquared;
  numeratorDegrees = If[TrueQ[# == 0], emptyComputedObject[rootsInQ],
      Exponent[Numerator[Together[#]], q]] & /@ rootsInQ;
  denominatorDegrees = Exponent[Denominator[Together[#]], q] & /@ rootsInQ;
  <|
    "RootsInQ" -> rootsInQ, "WavevectorOccurrences" -> wavevectorOccurrences,
    "CandidateSpeedsSquared" -> candidateSpeedsSquared,
    "HomogeneityDefects" -> homogeneityDefects,
    "SpeedVariations" -> speedVariations,
    "NumeratorDegrees" -> numeratorDegrees,
    "DenominatorDegrees" -> denominatorDegrees
  |>
];

deriveScalingBlock[data_] := Module[
  {scalingRules, scalingAssumptions, scaledOmegaSquared,
   lambdaSquaredOmegaSquared, scalingResiduals},
  scalingRules = Thread[waveCovector -> lambdaScale waveCovector];
  scalingAssumptions = assumptionSet && lambdaScale > 0;
  scaledOmegaSquared = casSimplify[# /. scalingRules,
      scalingAssumptions] & /@ data["Roots"];
  lambdaSquaredOmegaSquared = casSimplify[lambdaScale^2 #,
      scalingAssumptions] & /@ data["Roots"];
  scalingResiduals = casSimplify[(# /. scalingRules) - lambdaScale^2 #,
      scalingAssumptions] & /@ data["Roots"];
  <|
    "ScaledOmegaSquared" -> scaledOmegaSquared,
    "LambdaSquaredOmegaSquared" -> lambdaSquaredOmegaSquared,
    "ScalingResiduals" -> scalingResiduals
  |>
];

directionSpecifications = {
  {"GENERIC", {}},
  {"KX_KY_ZERO", {kz -> 0}},
  {"ZERO_ZERO_KZ", {kx -> 0, ky -> 0}},
  {"KX_ZERO_ZERO", {ky -> 0, kz -> 0}},
  {"KX_KX_KX", {ky -> kx, kz -> kx}}
};

specialDirectionDiagnostics[data_, rules_] := Module[
  {specialMatrix, specialWaveVector, specialGenerator, specialAssumptions},
  specialMatrix = data["MatrixOmegaSquared"] /. rules;
  specialWaveVector = waveCovector /. rules;
  specialGenerator = transverseGenerator /. rules;
  specialAssumptions = assumptionSet /. rules;
  rootDiagnostic[specialMatrix, # /. rules, specialWaveVector,
      specialGenerator, specialAssumptions] & /@ data["Roots"]
];

emitSpecialDirectionDiagnostics[prefix_, data_, slotCount_] := Module[
  {label, rules, diagnostics, index, slot, diagnostic, base},
  Do[
    label = specification[[1]];
    rules = specification[[2]];
    diagnostics = specialDirectionDiagnostics[data, rules];
    Do[
      slot = Pick[diagnostics, Range[Length[diagnostics]], index];
      base = prefix <> "ROOT" <> ToString[index] <> "_" <> label;
      If[slot === {},
        emit[base <> "_E1", emptyComputedObject[slot]];
        emit[base <> "_E2", emptyComputedObject[slot]];
        emit[base <> "_E3", emptyComputedObject[slot]];
        emit[base <> "_E4", emptyComputedObject[slot]],
        diagnostic = First[slot];
        emit[base <> "_E1", diagnostic["E1Test"]];
        emit[base <> "_E2", diagnostic["E2"]];
        emit[base <> "_E3", diagnostic["E3Test"]];
        emit[base <> "_E4", diagnostic["E4Test"]]
      ],
      {index, slotCount}
    ],
    {specification, directionSpecifications}
  ]
];

locusDiagnostics[data_, rules_, locusAssumptions_] := Module[
  {locusMatrix, locusWaveVector, locusGenerator},
  locusMatrix = data["MatrixOmegaSquared"] /. rules;
  locusWaveVector = waveCovector /. rules;
  locusGenerator = transverseGenerator /. rules;
  rootDiagnostic[locusMatrix, # /. rules, locusWaveVector,
      locusGenerator, locusAssumptions] & /@ data["Roots"]
];

emitLocusDiagnostics[prefix_, data_, slotCount_, label_, rules_,
    locusAssumptions_] := Module[
  {diagnostics, index, slot, diagnostic, base},
  diagnostics = locusDiagnostics[data, rules, locusAssumptions];
  Do[
    slot = Pick[diagnostics, Range[Length[diagnostics]], index];
    base = prefix <> "ROOT" <> ToString[index] <> "_" <> label;
    If[slot === {},
      emit[base <> "_E1", emptyComputedObject[slot]];
      emit[base <> "_E2", emptyComputedObject[slot]];
      emit[base <> "_E3", emptyComputedObject[slot]];
      emit[base <> "_E4", emptyComputedObject[slot]],
      diagnostic = First[slot];
      emit[base <> "_E1", diagnostic["E1Test"]];
      emit[base <> "_E2", diagnostic["E2"]];
      emit[base <> "_E3", diagnostic["E3Test"]];
      emit[base <> "_E4", diagnostic["E4Test"]]
    ],
    {index, slotCount}
  ]
];

valueAtSlot[values_, index_] := Module[{slot},
  slot = Pick[values, Range[Length[values]], index];
  If[slot === {}, emptyComputedObject[slot], First[slot]]
];

emitPackage[identifier_, data_, rootSlotCount_, speedSlotCount_] := Module[
  {prefix, speedData, scalingData, dimensionData, zeroMatrix, index,
   zeroWavevectorRules, zeroWavevectorAssumptions, equalDensityRules,
   equalDensityAssumptions},
  prefix = If[identifier === "", "WL_S9_", "WL_S9_" <> identifier <> "_"];
  speedData = deriveSpeedBlock[data];
  scalingData = deriveScalingBlock[data];
  dimensionData = deriveDimensionBlock[
    data["Lagrangian"], speedData["CandidateSpeedsSquared"]];
  zeroMatrix = ConstantArray[0, Dimensions[data["MatrixResidual"]]];

  emit[prefix <> "ASSUMPTIONS", assumptionSet];
  emit[prefix <> "LAGRANGIAN", data["Lagrangian"]];
  emit[prefix <> "EULER_LAGRANGE_RESIDUAL", data["EulerResidual"]];
  emit[prefix <> "EQUATIONS", data["Equations"]];
  emit[prefix <> "PLANE_WAVE_ANSATZ", planeWaveAnsatz];
  emit[prefix <> "QUADRATIC_PLANE_WAVE_ANSATZ", quadraticPlaneWaveAnsatz];
  emit[prefix <> "PLANE_WAVE_SUBSTITUTED_RESIDUAL", data["SubstitutedResidual"]];
  emit[prefix <> "PLANE_WAVE_RESIDUAL", data["PlaneResidual"]];
  emit[prefix <> "PLANE_WAVE_LAGRANGIAN_QUADRATIC", data["PlaneWaveQuadratic"]];
  emit[prefix <> "DYNAMICAL_MATRIX_A", data["MatrixA"]];
  emit[prefix <> "DYNAMICAL_MATRIX_B", data["MatrixB"]];
  emit[prefix <> "DYNAMICAL_MATRIX_RESIDUAL", data["MatrixResidual"]];
  If[!SameQ[data["MatrixResidual"], zeroMatrix], Exit[92]];

  emit[prefix <> "DETERMINANT", data["Determinant"]];
  emit[prefix <> "DETERMINANT_DISCARDED_ODD_PART", data["DeterminantOddPart"]];
  emit[prefix <> "DYNAMICAL_MATRIX_DISCARDED_ODD_PART", data["MatrixOddPart"]];
  If[!SameQ[data["DeterminantOddPart"], 0] ||
      !SameQ[data["MatrixOddPart"], zeroMatrix], Exit[93]];

  emit[prefix <> "OMEGA_SQUARED_SOLUTIONS", data["RootRules"]];
  emit[prefix <> "TRANSVERSE_GENERATOR", transverseGenerator];
  emit[prefix <> "LONGITUDINAL_GENERATOR", longitudinalGenerator];
  emitRootDiagnostics[prefix, data["Diagnostics"], rootSlotCount];
  emit[prefix <> "E1_ROOT_SUBSET", data["TransverseRoots"]];
  emit[prefix <> "E3_ROOT_SUBSET", data["LongitudinalRoots"]];
  emit[prefix <> "E4_ROOT_SUBSET", data["TransverseSubspaceRoots"]];
  emit[prefix <> "ROOT_MULTISET", data["RootMultiset"]];
  emit[prefix <> "ROOT_NULLITIES", data["RootNullities"]];
  emitSpecialDirectionDiagnostics[prefix, data, rootSlotCount];
  If[identifier === "",
    zeroWavevectorRules = Thread[waveCovector -> ConstantArray[0,
      Length[waveCovector]]];
    zeroWavevectorAssumptions =
      (assumptionSet /. (Last[positiveAssumptions] -> True)) /.
        zeroWavevectorRules;
    emitLocusDiagnostics[prefix, data, Length[data["Diagnostics"]],
      "ZERO_ZERO_ZERO", zeroWavevectorRules, zeroWavevectorAssumptions]
  ];
  If[identifier === "X6",
    equalDensityRules = {rhoZ -> rhoBr};
    equalDensityAssumptions = assumptionSet /. equalDensityRules;
    emitLocusDiagnostics[prefix, data, Length[data["Diagnostics"]],
      "RHO_Z_RHO_BR", equalDensityRules, equalDensityAssumptions]
  ];

  Do[
    emit[prefix <> "ROOT" <> ToString[index] <>
      "_WAVEVECTOR_SCALED_OMEGA_SQUARED",
      valueAtSlot[scalingData["ScaledOmegaSquared"], index]];
    emit[prefix <> "ROOT" <> ToString[index] <>
      "_LAMBDA_SCALE_SQUARED_OMEGA_SQUARED",
      valueAtSlot[scalingData["LambdaSquaredOmegaSquared"], index]];
    emit[prefix <> "ROOT" <> ToString[index] <> "_SCALING_RESIDUAL",
      valueAtSlot[scalingData["ScalingResiduals"], index]],
    {index, rootSlotCount}
  ];

  Do[
    emit[prefix <> "E1_ROOT_Q" <> ToString[index],
      valueAtSlot[speedData["RootsInQ"], index]];
    emit[prefix <> "E1_ROOT_Q_WAVEVECTOR_OCCURRENCE" <> ToString[index],
      valueAtSlot[speedData["WavevectorOccurrences"], index]];
    emit[prefix <> "CANDIDATE_SPEED_SQUARED" <> ToString[index],
      valueAtSlot[speedData["CandidateSpeedsSquared"], index]];
    emit[prefix <> "HOMOGENEITY_DEFECT" <> ToString[index],
      valueAtSlot[speedData["HomogeneityDefects"], index]];
    emit[prefix <> "SPEED_VARIATION" <> ToString[index],
      valueAtSlot[speedData["SpeedVariations"], index]];
    emit[prefix <> "NUMERATOR_Q_DEGREE" <> ToString[index],
      valueAtSlot[speedData["NumeratorDegrees"], index]];
    emit[prefix <> "DENOMINATOR_Q_DEGREE" <> ToString[index],
      valueAtSlot[speedData["DenominatorDegrees"], index]],
    {index, speedSlotCount}
  ];

  emit[prefix <> "LAGRANGIAN_TERM_DERIVATIVE_ORDERS",
    dimensionData["DerivativeOrders"]];
  emit[prefix <> "LAGRANGIAN_TERM_DIMENSION_EXPRESSIONS",
    dimensionData["TermDimensions"]];
  emit[prefix <> "ENERGY_DENSITY_DIMENSION",
    dimensionData["EnergyDensityDimension"]];
  emit[prefix <> "DIMENSION_EQUATION_LEFT", dimensionData["EquationLeft"]];
  emit[prefix <> "DIMENSION_EQUATION_RIGHT", dimensionData["EquationRight"]];
  emit[prefix <> "DIMENSION_EQUATION_RESIDUAL", dimensionData["EquationResidual"]];
  emit[prefix <> "DIMENSION_EQUATIONS", dimensionData["Equations"]];
  emit[prefix <> "DIMENSION_SOLUTION", dimensionData["SolutionSet"]];
  If[dimensionData["SolutionSet"] === {}, Exit[94]];
  emit[prefix <> "RHO_DIMENSION", dimensionData["RhoDimension"]];
  emit[prefix <> "RHO_Z_DIMENSION", dimensionData["RhoZDimension"]];
  emit[prefix <> "MU_DIMENSION", dimensionData["MuDimension"]];
  emit[prefix <> "MU_F_DIMENSION", dimensionData["MuFDimension"]];
  emit[prefix <> "MU_G_DIMENSION", dimensionData["MuGDimension"]];
  emit[prefix <> "MU_RHO_DIMENSION_DIFFERENCE",
    dimensionData["MuRhoDifference"]];
  emit[prefix <> "SPEED_SQUARED_PARAMETER_SYMBOLS", physicalCoefficientSymbols];
  emit[prefix <> "SPEED_SQUARED_PARAMETER_POWERS",
    dimensionData["ParameterPowers"]];
  emit[prefix <> "SPEED_SQUARED_SYMBOL_DIMENSION_MAP",
    dimensionData["SymbolDimensionMap"]];
  emit[prefix <> "SPEED_SQUARED_IMPLIED_DIMENSION",
    dimensionData["ImpliedDimensions"]];
  emit[prefix <> "VELOCITY_SQUARED_DIMENSION",
    dimensionData["VelocitySquaredDimension"]];
  emit[prefix <> "SPEED_SQUARED_DIMENSION_DIFFERENCE",
    dimensionData["DimensionDifferences"]];
  emit[prefix <> "SPEED_SQUARED_DIMENSION_UNKNOWN_SYMBOLS",
    dimensionData["DimensionUnknownSymbols"]];
  emit[prefix <> "SPEED_SQUARED_DIMENSION_SUM_MISMATCHES",
    dimensionData["DimensionSumMismatches"]];
  emit[prefix <> "SPEED_SQUARED_DIMENSION_UNSUPPORTED_EXPRESSIONS",
    dimensionData["DimensionUnsupportedExpressions"]];
  If[!AllTrue[dimensionData["DimensionUnknownSymbols"], # === {} &],
    Exit[95]
  ];
  If[!AllTrue[dimensionData["DimensionSumMismatches"], # === {} &],
    Exit[96]
  ];
  If[!AllTrue[dimensionData["DimensionUnsupportedExpressions"], # === {} &],
    Exit[97]
  ];
  emit[prefix <> "RHO_DIMENSION_D3", dimensionData["RhoDimension"] /. D -> 3];
  emit[prefix <> "RHO_Z_DIMENSION_D3", dimensionData["RhoZDimension"] /. D -> 3];
  emit[prefix <> "MU_DIMENSION_D3", dimensionData["MuDimension"] /. D -> 3];
  emit[prefix <> "MU_F_DIMENSION_D3", dimensionData["MuFDimension"] /. D -> 3];
  emit[prefix <> "MU_G_DIMENSION_D3", dimensionData["MuGDimension"] /. D -> 3];
  emit[prefix <> "MU_RHO_DIMENSION_DIFFERENCE_D3",
    dimensionData["MuRhoDifference"] /. D -> 3];
];

actionRecords = {
  {"", mainLagrangian}, {"X1", x1Lagrangian}, {"X2", x2Lagrangian},
  {"X3", x3Lagrangian}, {"X4", x4Lagrangian}, {"X5", x5Lagrangian},
  {"X6", x6Lagrangian}, {"X7", x7Lagrangian}, {"X8", x8Lagrangian}
};
derivedRecords = {#[[1]], deriveFromAction[#[[2]]]} & /@ actionRecords;
rootSlotCount = Max[Length[#[[2]]["Diagnostics"]] & /@ derivedRecords];
speedSlotCount = Max[Length[#[[2]]["TransverseRoots"]] & /@ derivedRecords];

emitPackage[#[[1]], #[[2]], rootSlotCount, speedSlotCount] & /@ derivedRecords;

Exit[0];
