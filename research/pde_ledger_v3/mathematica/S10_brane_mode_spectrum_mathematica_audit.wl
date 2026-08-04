$HistoryLength = 0;

ClearAll["Global`*"];
$HistoryLength = 0;

emittedTagCount = 0;
localTagNames = {};

emit[tag_String, expression_] := Module[{payload},
  payload = ToString[expression, InputForm, PageWidth -> Infinity];
  Print[tag <> ": " <> payload];
  emittedTagCount = emittedTagCount + 1;
];

emitLocal[tag_String, expression_] := Module[{},
  localTagNames = Append[localTagNames, tag];
  emit[tag, expression];
];

localizeTag[tag_String] :=
  StringReplace[tag, StartOfString ~~ "WL_S10_" -> "WL_S10_LOCAL_"];

emitSolverObject[tag_String, expression_] := Module[
  {held, conditionalParts, localBase},
  held = expression;
  localBase = localizeTag[tag];
  If[Head[held] === ConditionalExpression,
    emit[tag, held[[1]]];
    emitLocal[localBase <> "_SOLVER_CONDITION", held[[2]]],
    emit[tag, held];
    conditionalParts = Cases[
      held,
      ConditionalExpression[value_, condition_] :> {value, condition},
      Infinity
    ];
    Do[
      emitLocal[
        localBase <> "_SOLVER_CONDITION" <> ToString[index],
        conditionalParts[[index, 2]]
      ],
      {index, Length[conditionalParts]}
    ];
  ];
];

fullSimplifyWith[expression_, assumptions_] :=
  FullSimplify[expression, Assumptions -> assumptions];

refineWith[expression_, assumptions_] :=
  Refine[expression, Assumptions -> assumptions];

reduceAllowed[formula_, variables_List, assumptions_] :=
  Reduce[And[formula, assumptions], variables, Reals];

realLocusSolve[formula_, variables_List] :=
  Quiet[Solve[formula, variables, Reals], Solve::svars];

braneDimension = Symbol["braneDimension"];
rhoBr = Symbol["rhoBr"];
muR = Symbol["muR"];
sRho = Symbol["sRho"];
coefficientScale = Symbol["coefficientScale"];
omegaSquared = Symbol["omegaSquared"];
lambdaScale = Symbol["lambdaScale"];
timeCoordinate = Symbol["t"];

packagePrefix[package_String, dimension_Integer] :=
  "WL_S10_" <> package <> "_D" <> ToString[dimension];

rootPrefix[prefix_String, rootIndex_Integer] :=
  prefix <> "_ROOT" <> ToString[rootIndex];

packageControlAssumptions[package_String] := Switch[
  package,
  "XFORM_ANISO",
    And[sRho > 0, Element[sRho, Reals], Unequal[sRho, 1]],
  "XCOEF_SCALE",
    And[coefficientScale > 0, Element[coefficientScale, Reals]],
  _,
    True
];

jointAssumptionExpression[
  package_String,
  wavevector_List,
  amplitudes_List
] := Apply[
  And,
  Join[
    {
      rhoBr > 0,
      muR > 0,
      Total[wavevector^2] > 0
    },
    Element[#, Reals] & /@ wavevector,
    Element[#, Reals] & /@ amplitudes,
    {
      Element[braneDimension, Integers],
      braneDimension > 0,
      packageControlAssumptions[package]
    }
  ]
];

buildPackage[package_String, dimension_Integer] := Module[
  {
    coordinates, fieldHeads, fields, amplitudes, wavevector, arguments,
    velocities, gradient, curlStiffness, fullGradientStiffness,
    divergenceStiffness, stiffnessDensity, inertialCoefficients,
    stiffnessCoefficient, stiffnessMultiplier, kineticDensity,
    lagrangian, controlData, dimensionPremises
  },
  coordinates = Table[Symbol["x" <> ToString[index]], {index, dimension}];
  fieldHeads = Table[Symbol["u" <> ToString[index]], {index, dimension}];
  amplitudes = Table[Symbol["a" <> ToString[index]], {index, dimension}];
  wavevector = Table[Symbol["k" <> ToString[index]], {index, dimension}];
  arguments = Join[coordinates, {timeCoordinate}];
  fields = (#[Sequence @@ arguments] &) /@ fieldHeads;
  velocities = D[#, timeCoordinate] & /@ fields;
  gradient = Table[
    D[fields[[fieldIndex]], coordinates[[coordinateIndex]]],
    {coordinateIndex, dimension},
    {fieldIndex, dimension}
  ];

  curlStiffness = (1/2) Sum[
    (gradient[[coordinateIndex, fieldIndex]] -
       gradient[[fieldIndex, coordinateIndex]])^2,
    {coordinateIndex, dimension},
    {fieldIndex, dimension}
  ];
  fullGradientStiffness = Sum[
    gradient[[coordinateIndex, fieldIndex]]^2,
    {coordinateIndex, dimension},
    {fieldIndex, dimension}
  ];
  divergenceStiffness = Total[Diagonal[gradient]]^2;

  inertialCoefficients = Switch[
    package,
    "XFORM_ANISO",
      Join[{sRho rhoBr}, ConstantArray[rhoBr, dimension - 1]],
    _,
      ConstantArray[rhoBr, dimension]
  ];
  stiffnessDensity = Switch[
    package,
    "XFORM_FULLGRAD", fullGradientStiffness,
    "XFORM_DIVONLY", divergenceStiffness,
    _, curlStiffness
  ];
  stiffnessCoefficient = Switch[
    package,
    "XCOEF_SCALE", coefficientScale muR,
    _, muR
  ];
  stiffnessMultiplier = Switch[
    package,
    "XFORM_SIGNFLIP", 1,
    _, -1
  ];
  kineticDensity = (1/2) Total[
    inertialCoefficients velocities^2
  ];
  lagrangian = kineticDensity +
    stiffnessMultiplier (stiffnessCoefficient/2) stiffnessDensity;

  controlData = {
    actionInertialCoefficients -> inertialCoefficients,
    actionStiffnessDensity -> stiffnessDensity,
    actionStiffnessCoefficient -> stiffnessCoefficient,
    actionStiffnessMultiplier -> stiffnessMultiplier
  };
  dimensionPremises = Switch[
    package,
    "XFORM_ANISO", {dimensionOf[sRho] == {0, 0, 0}},
    "XCOEF_SCALE", {dimensionOf[coefficientScale] == {0, 0, 0}},
    _, {}
  ];

  <|
    "Coordinates" -> coordinates,
    "FieldHeads" -> fieldHeads,
    "Fields" -> fields,
    "Amplitudes" -> amplitudes,
    "Wavevector" -> wavevector,
    "Velocities" -> velocities,
    "Gradient" -> gradient,
    "CurlStiffness" -> curlStiffness,
    "StiffnessDensity" -> stiffnessDensity,
    "InertialCoefficients" -> inertialCoefficients,
    "StiffnessCoefficients" -> {stiffnessCoefficient},
    "Lagrangian" -> lagrangian,
    "ControlData" -> controlData,
    "DimensionPremises" -> dimensionPremises
  |>
];

eulerLagrangeSystem[data_Association] := Module[
  {lagrangian, fields, velocities, gradient, coordinates, dimension},
  lagrangian = data["Lagrangian"];
  fields = data["Fields"];
  velocities = data["Velocities"];
  gradient = data["Gradient"];
  coordinates = data["Coordinates"];
  dimension = Length[fields];
  Expand@Table[
    D[D[lagrangian, velocities[[fieldIndex]]], timeCoordinate] +
      Sum[
        D[
          D[lagrangian, gradient[[coordinateIndex, fieldIndex]]],
          coordinates[[coordinateIndex]]
        ],
        {coordinateIndex, dimension}
      ] - D[lagrangian, fields[[fieldIndex]]],
    {fieldIndex, dimension}
  ]
];

ansatzRules[data_Association] := Module[
  {
    dimension, dummyArguments, dummyCoordinates, dummyTime,
    amplitudes, wavevector
  },
  dimension = Length[data["Fields"]];
  dummyArguments = Table[Unique["ansatzArgument"], {dimension + 1}];
  dummyCoordinates = Take[dummyArguments, dimension];
  dummyTime = Last[dummyArguments];
  amplitudes = data["Amplitudes"];
  wavevector = data["Wavevector"];
  MapThread[
    Rule[
      #1,
      Function[
        Evaluate[dummyArguments],
        Evaluate[#2 Cos[
          Total[wavevector dummyCoordinates] -
            Sqrt[omegaSquared] dummyTime
        ]]
      ]
    ] &,
    {data["FieldHeads"], amplitudes}
  ]
];

periodAverage[expression_, assumptions_] := Module[{integral},
  integral = Integrate[
    expression,
    {timeCoordinate, 0, 2 Pi/Sqrt[omegaSquared]},
    Assumptions -> And[assumptions, omegaSquared > 0],
    GenerateConditions -> False
  ];
  fullSimplifyWith[
    (Sqrt[omegaSquared]/(2 Pi)) integral,
    assumptions
  ]
];

deriveMatrices[data_Association, eom_List, assumptions_] := Module[
  {
    rules, coordinates, amplitudes, wavevector, phase, commonFactor,
    eomOnAnsatz, strippedEquations, matrixA, lagrangianOnAnsatz,
    averagedLagrangian, matrixB, dimension
  },
  rules = ansatzRules[data];
  coordinates = data["Coordinates"];
  amplitudes = data["Amplitudes"];
  wavevector = data["Wavevector"];
  dimension = Length[amplitudes];
  phase = Total[wavevector coordinates] -
    Sqrt[omegaSquared] timeCoordinate;
  commonFactor = Cos[phase];
  eomOnAnsatz = Expand[eom /. rules];
  strippedEquations =
    fullSimplifyWith[Together[#/commonFactor], assumptions] & /@
      eomOnAnsatz;
  matrixA = Table[
    fullSimplifyWith[
      Coefficient[strippedEquations[[row]], amplitudes[[column]]],
      assumptions
    ],
    {row, dimension},
    {column, dimension}
  ];

  lagrangianOnAnsatz = Expand[data["Lagrangian"] /. rules];
  averagedLagrangian = periodAverage[lagrangianOnAnsatz, assumptions];
  matrixB = Table[
    fullSimplifyWith[
      D[averagedLagrangian, amplitudes[[row]], amplitudes[[column]]],
      assumptions
    ],
    {row, dimension},
    {column, dimension}
  ];

  <|
    "AnsatzRules" -> rules,
    "EOMOnAnsatz" -> eomOnAnsatz,
    "StrippedEquations" -> strippedEquations,
    "MatrixA" -> matrixA,
    "LagrangianOnAnsatz" -> lagrangianOnAnsatz,
    "AveragedLagrangian" -> averagedLagrangian,
    "MatrixB" -> matrixB
  |>
];

distinctUnderAssumptions[expressions_List, assumptions_] :=
  DeleteDuplicates[
    expressions,
    TrueQ[fullSimplifyWith[#1 == #2, assumptions]] &
  ];

removeConditionalWrappers[expression_] :=
  expression /. ConditionalExpression[value_, condition_] :> value;

runSpectrum[
  prefix_String,
  matrix_List,
  wavevector_List,
  assumptions_
] := Module[
  {
    determinant, solutions, rootsFromSolutions, roots, pairs,
    coincidenceRecords, pair, pairPrefix, equation, locus, allowed,
    intersects
  },
  determinant = Factor[Det[matrix]];
  emit[prefix <> "_Q3_DETERMINANT", determinant];
  solutions = Solve[determinant == 0, omegaSquared];
  solutions = fullSimplifyWith[solutions, assumptions];
  emitSolverObject[prefix <> "_Q3_SOLUTIONS", solutions];
  rootsFromSolutions = Flatten[{
    omegaSquared /. removeConditionalWrappers[solutions]
  }];
  rootsFromSolutions = Select[
    rootsFromSolutions,
    FreeQ[#, omegaSquared] &
  ];
  roots = distinctUnderAssumptions[rootsFromSolutions, assumptions];
  emit[prefix <> "_Q3_DISTINCT_ROOTS", roots];
  emit[prefix <> "_Q3_ROOT_COUNT", Length[roots]];
  emit[prefix <> "_ROOT_ORDERING", roots];
  Do[
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q3_SIGN",
      refineWith[Sign[roots[[rootIndex]]], assumptions]
    ],
    {rootIndex, Length[roots]}
  ];

  pairs = Subsets[Range[Length[roots]], {2}];
  coincidenceRecords = {};
  Do[
    pair = pairs[[pairIndex]];
    pairPrefix = prefix <> "_ROOT" <> ToString[pair[[1]]] <>
      "_ROOT" <> ToString[pair[[2]]] <> "_Q3_COINCIDENCE";
    equation = Together[
      roots[[pair[[1]]]] - roots[[pair[[2]]]]
    ] == 0;
    locus = realLocusSolve[equation, wavevector];
    allowed = reduceAllowed[equation, wavevector, assumptions];
    intersects = Not[SameQ[allowed, False]];
    emit[pairPrefix <> "_OPERANDS", {equation, assumptions}];
    emitSolverObject[pairPrefix <> "_LOCUS", locus];
    emit[pairPrefix <> "_ALLOWED_INTERSECTION", allowed];
    emit[pairPrefix <> "_INTERSECTION_TEST", intersects];
    coincidenceRecords = Append[
      coincidenceRecords,
      <|
        "Pair" -> pair,
        "Equation" -> equation,
        "Locus" -> locus,
        "Allowed" -> allowed,
        "Intersects" -> intersects
      |>
    ],
    {pairIndex, Length[pairs]}
  ];
  emit[prefix <> "_Q3_ROOT_COINCIDENCE_LOCI", coincidenceRecords];

  <|
    "Determinant" -> determinant,
    "Solutions" -> solutions,
    "Roots" -> roots,
    "CoincidenceRecords" -> coincidenceRecords
  |>
];

computeModeData[
  prefix_String,
  matrix_List,
  wavevector_List,
  root_,
  rootIndex_Integer,
  assumptions_
] := Module[
  {
    dimension, zeroTest, rootMatrix, rank, nullity, stackedMatrix,
    stackedRank, transverseNullity, nullityDifference,
    wavevectorProduct, basis, basisDots, basisResiduals, basisCount,
    basisCountResidual, squaredWavevector, currentPrefix
  },
  dimension = Length[wavevector];
  zeroTest = Function[
    entry,
    TrueQ[fullSimplifyWith[entry == 0, assumptions]]
  ];
  rootMatrix = Map[
    fullSimplifyWith[#, assumptions] &,
    matrix /. omegaSquared -> root,
    {2}
  ];
  rank = MatrixRank[rootMatrix, ZeroTest -> zeroTest];
  nullity = dimension - rank;
  stackedMatrix = Join[rootMatrix, {wavevector}];
  stackedRank = MatrixRank[stackedMatrix, ZeroTest -> zeroTest];
  transverseNullity = dimension - stackedRank;
  nullityDifference = nullity - transverseNullity;
  wavevectorProduct = fullSimplifyWith[rootMatrix . wavevector, assumptions];
  basis = NullSpace[rootMatrix, ZeroTest -> zeroTest];
  basisDots = fullSimplifyWith[# . wavevector, assumptions] & /@ basis;
  squaredWavevector = Total[wavevector^2];
  basisResiduals = MapThread[
    fullSimplifyWith[
      squaredWavevector #1 - #2 wavevector,
      assumptions
    ] &,
    {basis, basisDots}
  ];
  basisCount = Length[basis];
  basisCountResidual = basisCount - nullity;
  currentPrefix = rootPrefix[prefix, rootIndex];

  emit[currentPrefix <> "_N1_MATRIX", rootMatrix];
  emit[currentPrefix <> "_N2_RANK", rank];
  emit[currentPrefix <> "_N2_NULLITY", nullity];
  emit[currentPrefix <> "_N3_STACKED_MATRIX", stackedMatrix];
  emit[currentPrefix <> "_N3_STACKED_RANK", stackedRank];
  emit[currentPrefix <> "_N3_TRANSVERSE_NULLITY", transverseNullity];
  emit[currentPrefix <> "_N4_NULLITY_DIFFERENCE", nullityDifference];
  emit[currentPrefix <> "_N5_WAVEVECTOR_PRODUCT", wavevectorProduct];
  emit[currentPrefix <> "_N6_NULLSPACE_BASIS", basis];
  emit[currentPrefix <> "_N6_BASIS_DOTS", basisDots];
  emit[currentPrefix <> "_N6_BASIS_RESIDUALS", basisResiduals];
  emit[currentPrefix <> "_N7_BASIS_COUNT", basisCount];
  emit[currentPrefix <> "_N7_COUNT_RESIDUAL", basisCountResidual];

  <|
    "Root" -> root,
    "RootMatrix" -> rootMatrix,
    "Rank" -> rank,
    "Nullity" -> nullity,
    "StackedMatrix" -> stackedMatrix,
    "StackedRank" -> stackedRank,
    "TransverseNullity" -> transverseNullity,
    "Basis" -> basis
  |>
];

runModeSet[
  prefix_String,
  matrix_List,
  wavevector_List,
  roots_List,
  assumptions_
] := Table[
  computeModeData[
    prefix,
    matrix,
    wavevector,
    roots[[rootIndex]],
    rootIndex,
    assumptions
  ],
  {rootIndex, Length[roots]}
];

runScaling[
  prefix_String,
  roots_List,
  wavevector_List,
  assumptions_
] := Module[
  {
    scaledRoots, ratios, ratioRecords, root, scaledRoot, zeroTest,
    ratio, numeratorExponent, denominatorExponent, exponentCandidate,
    purePowerTest, currentPrefix
  },
  scaledRoots = {};
  ratios = {};
  ratioRecords = {};
  Do[
    root = roots[[rootIndex]];
    currentPrefix = rootPrefix[prefix, rootIndex];
    scaledRoot = fullSimplifyWith[
      root /. Thread[wavevector -> lambdaScale wavevector],
      assumptions
    ];
    zeroTest = TrueQ[fullSimplifyWith[root == 0, assumptions]];
    emit[currentPrefix <> "_Q5_SCALED_ROOT", scaledRoot];
    emit[currentPrefix <> "_Q5_ORIGINAL_ROOT", root];
    emit[currentPrefix <> "_Q5_RATIO_DOMAIN", Not[zeroTest]];
    scaledRoots = Append[scaledRoots, scaledRoot];
    If[Not[zeroTest],
      ratio = fullSimplifyWith[
        Cancel[Together[scaledRoot/root]],
        assumptions
      ];
      numeratorExponent = Exponent[
        Numerator[Together[ratio]],
        lambdaScale
      ];
      denominatorExponent = Exponent[
        Denominator[Together[ratio]],
        lambdaScale
      ];
      exponentCandidate = numeratorExponent - denominatorExponent;
      purePowerTest = FreeQ[
        fullSimplifyWith[
          Cancel[ratio/lambdaScale^exponentCandidate],
          assumptions
        ],
        lambdaScale
      ];
      emit[currentPrefix <> "_Q5_SCALING_RATIO", ratio];
      emit[currentPrefix <> "_Q5_PURE_POWER_TEST", purePowerTest];
      emit[currentPrefix <> "_Q5_SCALING_EXPONENT_CANDIDATE", exponentCandidate];
      If[purePowerTest,
        emit[currentPrefix <> "_Q5_SCALING_EXPONENT", exponentCandidate]
      ];
      ratios = Append[ratios, ratio];
      ratioRecords = Append[
        ratioRecords,
        <|
          "RootIndex" -> rootIndex,
          "Ratio" -> ratio,
          "ExponentCandidate" -> exponentCandidate,
          "PurePowerTest" -> purePowerTest
        |>
      ],
      ratioRecords = Append[
        ratioRecords,
        <|
          "RootIndex" -> rootIndex,
          "RatioDomain" -> Not[zeroTest]
        |>
      ]
    ],
    {rootIndex, Length[roots]}
  ];
  <|
    "ScaledRoots" -> scaledRoots,
    "Ratios" -> ratios,
    "RatioRecords" -> ratioRecords
  |>
];

makeDimensionAnalyzer[
  fieldHeads_List,
  symbolDimensions_Association
] := Module[{fieldDimension, zeroDimension},
  fieldDimension = {1, 0, 0};
  zeroDimension = {0, 0, 0};
  Function[expression,
    Module[{walk},
      walk[current_] := Module[
        {
          children, dimensions, groups, constraints, orders,
          differentiatedField, firstDimension
        },
        Which[
          NumericQ[current],
            <|
              "Dimension" -> zeroDimension,
              "AddendDimensions" -> {},
              "Constraints" -> {}
            |>,
          KeyExistsQ[symbolDimensions, current],
            <|
              "Dimension" -> symbolDimensions[current],
              "AddendDimensions" -> {},
              "Constraints" -> {}
            |>,
          MemberQ[fieldHeads, Head[current]],
            <|
              "Dimension" -> fieldDimension,
              "AddendDimensions" -> {},
              "Constraints" -> {}
            |>,
          Head[Head[Head[current]]] === Derivative &&
              MemberQ[fieldHeads, First[List @@ Head[current]]],
            orders = List @@ Head[Head[current]];
            differentiatedField = First[List @@ Head[current]];
            <|
              "Dimension" ->
                fieldDimension +
                  {-Total[Most[orders]], -Last[orders], 0},
              "AddendDimensions" -> {},
              "Constraints" -> {}
            |>,
          Head[current] === Plus,
            children = walk /@ (List @@ current);
            dimensions = Lookup[children, "Dimension"];
            firstDimension = First[dimensions];
            groups = Join[
              Flatten[Lookup[children, "AddendDimensions"], 1],
              {dimensions}
            ];
            constraints = Join[
              Flatten[Lookup[children, "Constraints"], 1],
              Flatten[Thread[# == firstDimension] & /@ Rest[dimensions]]
            ];
            <|
              "Dimension" -> firstDimension,
              "AddendDimensions" -> groups,
              "Constraints" -> constraints
            |>,
          Head[current] === Times,
            children = walk /@ (List @@ current);
            <|
              "Dimension" -> Total[Lookup[children, "Dimension"]],
              "AddendDimensions" ->
                Flatten[Lookup[children, "AddendDimensions"], 1],
              "Constraints" ->
                Flatten[Lookup[children, "Constraints"], 1]
            |>,
          Head[current] === Power,
            children = walk[First[current]];
            <|
              "Dimension" -> Last[current] children["Dimension"],
              "AddendDimensions" -> children["AddendDimensions"],
              "Constraints" -> children["Constraints"]
            |>,
          Head[current] === Abs,
            walk[First[current]],
          Head[current] === ConditionalExpression,
            walk[First[current]],
          True,
            <|
              "Dimension" -> zeroDimension,
              "AddendDimensions" -> {},
              "Constraints" -> {}
            |>
        ]
      ];
      walk[expression]
    ]
  ]
];

resolveDimensionAnalysis[analysis_Association, rules_List, assumptions_] := Module[
  {dimension, termDimensions, constraints, homogeneity},
  dimension = fullSimplifyWith[analysis["Dimension"] /. rules, assumptions];
  termDimensions = fullSimplifyWith[
    analysis["AddendDimensions"] /. rules,
    assumptions
  ];
  constraints = analysis["Constraints"];
  homogeneity = And @@ (
    TrueQ[fullSimplifyWith[# /. rules, assumptions]] & /@ constraints
  );
  <|
    "Dimension" -> dimension,
    "AddendDimensions" -> termDimensions,
    "Homogeneity" -> homogeneity
  |>
];

resolvedExpressionDimension[
  expression_,
  analyzer_,
  rules_List,
  assumptions_
] := Module[{zeroTest, analysis},
  zeroTest = TrueQ[fullSimplifyWith[expression == 0, assumptions]];
  If[zeroTest,
    <|
      "Dimension" -> Indeterminate,
      "AddendDimensions" -> {},
      "Homogeneity" -> True
    |>,
    analysis = analyzer[expression];
    resolveDimensionAnalysis[analysis, rules, assumptions]
  ]
];

runDimensions[
  prefix_String,
  data_Association,
  eom_List,
  matrixData_Association,
  spectrumData_Association,
  scalingData_Association,
  assumptions_
] := Module[
  {
    dimensionNames, rhoDimension, muDimension, sRhoDimension,
    scaleDimension, coefficientDimensionMap, actionCoefficientSymbols,
    controlSymbols, unknowns, controlEquations, energyDensityDimension,
    expandedLagrangian, actionTerms, analyzer, actionTermAnalyses,
    actionTermDimensions, dimensionEquations, dimensionSolutions,
    dimensionRules, inertialAnalyses, stiffnessAnalyses,
    inertialResolved, stiffnessResolved, squaredWavevector,
    rateCoefficientExpressions, rateResolved, auditGroups, groupNames,
    groupExpressions, groupResolved
  },
  dimensionNames = {"Length", "Time", "Mass"};
  rhoDimension = Table[Symbol["dimRho" <> name], {name, dimensionNames}];
  muDimension = Table[Symbol["dimMu" <> name], {name, dimensionNames}];
  sRhoDimension = Table[Symbol["dimSRho" <> name], {name, dimensionNames}];
  scaleDimension = Table[Symbol["dimScale" <> name], {name, dimensionNames}];
  coefficientDimensionMap = <|
    rhoBr -> rhoDimension,
    muR -> muDimension,
    sRho -> sRhoDimension,
    coefficientScale -> scaleDimension,
    lambdaScale -> {0, 0, 0},
    omegaSquared -> {0, -2, 0}
  |>;
  Do[
    AssociateTo[
      coefficientDimensionMap,
      {
        data["Amplitudes"][[index]] -> {1, 0, 0},
        data["Wavevector"][[index]] -> {-1, 0, 0}
      }
    ],
    {index, Length[data["Amplitudes"]]}
  ];
  analyzer = makeDimensionAnalyzer[
    data["FieldHeads"],
    coefficientDimensionMap
  ];

  energyDensityDimension = {2 - braneDimension, -2, 1};
  expandedLagrangian = Expand[data["Lagrangian"]];
  actionTerms = If[
    Head[expandedLagrangian] === Plus,
    List @@ expandedLagrangian,
    {expandedLagrangian}
  ];
  actionTermAnalyses = analyzer /@ actionTerms;
  actionTermDimensions = Lookup[actionTermAnalyses, "Dimension"];
  actionCoefficientSymbols = Select[
    {rhoBr, muR, sRho, coefficientScale},
    Not[FreeQ[expandedLagrangian, #]] &
  ];
  controlSymbols = Intersection[
    actionCoefficientSymbols,
    {sRho, coefficientScale}
  ];
  unknowns = Flatten[Lookup[coefficientDimensionMap, actionCoefficientSymbols]];
  controlEquations = Flatten[
    Table[
      Thread[coefficientDimensionMap[symbol] == {0, 0, 0}],
      {symbol, controlSymbols}
    ]
  ];
  dimensionEquations = Join[
    Flatten[
      Thread[# == energyDensityDimension] & /@ actionTermDimensions
    ],
    controlEquations
  ];
  dimensionSolutions = Solve[dimensionEquations, unknowns];
  emit[prefix <> "_Q6_ENERGY_DENSITY_DIMENSION", energyDensityDimension];
  emit[prefix <> "_Q6_DIMENSION_PREMISES", data["DimensionPremises"]];
  emit[prefix <> "_Q6_DIMENSION_EQUATIONS", dimensionEquations];
  emitSolverObject[prefix <> "_Q6_DIMENSION_SOLUTIONS", dimensionSolutions];
  dimensionRules = If[Length[dimensionSolutions] > 0, First[dimensionSolutions], {}];

  inertialAnalyses = analyzer /@ data["InertialCoefficients"];
  stiffnessAnalyses = analyzer /@ data["StiffnessCoefficients"];
  inertialResolved =
    resolveDimensionAnalysis[#, dimensionRules, assumptions] & /@
      inertialAnalyses;
  stiffnessResolved =
    resolveDimensionAnalysis[#, dimensionRules, assumptions] & /@
      stiffnessAnalyses;
  emit[prefix <> "_Q6_INERTIAL_COEFFICIENTS", data["InertialCoefficients"]];
  emit[
    prefix <> "_Q6_INERTIAL_COEFFICIENT_DIMENSIONS",
    Lookup[inertialResolved, "Dimension"]
  ];
  emit[
    prefix <> "_Q6_INERTIAL_COEFFICIENT_TERM_DIMENSIONS",
    Lookup[inertialResolved, "AddendDimensions"]
  ];
  emit[
    prefix <> "_Q6_INERTIAL_COEFFICIENT_HOMOGENEITY",
    Lookup[inertialResolved, "Homogeneity"]
  ];
  emit[prefix <> "_Q6_STIFFNESS_COEFFICIENTS", data["StiffnessCoefficients"]];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_DIMENSIONS",
    Lookup[stiffnessResolved, "Dimension"]
  ];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_TERM_DIMENSIONS",
    Lookup[stiffnessResolved, "AddendDimensions"]
  ];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_HOMOGENEITY",
    Lookup[stiffnessResolved, "Homogeneity"]
  ];

  squaredWavevector = Total[data["Wavevector"]^2];
  rateCoefficientExpressions = Cancel[
      Together[#/squaredWavevector]
    ] & /@ spectrumData["Roots"];
  rateResolved = resolvedExpressionDimension[
      #,
      analyzer,
      dimensionRules,
      assumptions
    ] & /@ rateCoefficientExpressions;
  Do[
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q6_RATE_COEFFICIENT",
      rateCoefficientExpressions[[rootIndex]]
    ];
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q6_RATE_COEFFICIENT_DIMENSION",
      rateResolved[[rootIndex, "Dimension"]]
    ];
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q6_RATE_COEFFICIENT_TERM_DIMENSIONS",
      rateResolved[[rootIndex, "AddendDimensions"]]
    ];
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q6_RATE_COEFFICIENT_HOMOGENEITY",
      rateResolved[[rootIndex, "Homogeneity"]]
    ],
    {rootIndex, Length[rateCoefficientExpressions]}
  ];

  auditGroups = <|
    "LAGRANGIAN" -> {data["Lagrangian"]},
    "EOM" -> eom,
    "AVERAGED_LAGRANGIAN" -> {matrixData["AveragedLagrangian"]},
    "MATRIX_A_ENTRIES" -> Flatten[matrixData["MatrixA"]],
    "MATRIX_B_ENTRIES" -> Flatten[matrixData["MatrixB"]],
    "DETERMINANT" -> {spectrumData["Determinant"]},
    "ROOTS" -> spectrumData["Roots"],
    "SCALED_ROOTS" -> scalingData["ScaledRoots"],
    "SCALING_RATIOS" -> scalingData["Ratios"]
  |>;
  groupNames = Keys[auditGroups];
  Do[
    groupExpressions = auditGroups[groupNames[[groupIndex]]];
    groupResolved = resolvedExpressionDimension[
        #,
        analyzer,
        dimensionRules,
        assumptions
      ] & /@ groupExpressions;
    emit[
      prefix <> "_Q6_" <> groupNames[[groupIndex]] <> "_DIMENSIONS",
      Lookup[groupResolved, "Dimension"]
    ];
    emit[
      prefix <> "_Q6_" <> groupNames[[groupIndex]] <> "_TERM_DIMENSIONS",
      Lookup[groupResolved, "AddendDimensions"]
    ];
    emit[
      prefix <> "_Q6_" <> groupNames[[groupIndex]] <> "_HOMOGENEITY",
      Lookup[groupResolved, "Homogeneity"]
    ],
    {groupIndex, Length[groupNames]}
  ];

  <|
    "DimensionSolutions" -> dimensionSolutions,
    "DimensionRules" -> dimensionRules
  |>
];

runCurlComparison[prefix_String, dimension_Integer] := Module[
  {gradientSymbols, curlStiffness, curlVector, curlNorm},
  If[dimension == 3,
    gradientSymbols = Table[
      Symbol["g" <> ToString[row] <> "x" <> ToString[column]],
      {row, 3},
      {column, 3}
    ];
    curlStiffness = (1/2) Sum[
      (gradientSymbols[[row, column]] -
         gradientSymbols[[column, row]])^2,
      {row, 3},
      {column, 3}
    ];
    curlVector = {
      gradientSymbols[[2, 3]] - gradientSymbols[[3, 2]],
      gradientSymbols[[3, 1]] - gradientSymbols[[1, 3]],
      gradientSymbols[[1, 2]] - gradientSymbols[[2, 1]]
    };
    curlNorm = Expand[curlVector . curlVector];
    emit[prefix <> "_Q7_STIFFNESS", Expand[curlStiffness]];
    emit[prefix <> "_Q7_CURL_NORM", curlNorm];
    emit[prefix <> "_Q7_RESIDUAL", Expand[curlStiffness - curlNorm]];
  ];
];

pointSpectrum[
  prefix_String,
  matrix_List,
  wavevector_List,
  assumptions_
] := Module[
  {determinant, solutions, rootValues, roots, pairs, coincidenceTests},
  determinant = Factor[Det[matrix]];
  emit[prefix <> "_Q3_DETERMINANT", determinant];
  solutions = fullSimplifyWith[
    Solve[determinant == 0, omegaSquared],
    assumptions
  ];
  emitSolverObject[prefix <> "_Q3_SOLUTIONS", solutions];
  rootValues = Flatten[{
    omegaSquared /. removeConditionalWrappers[solutions]
  }];
  rootValues = Select[rootValues, FreeQ[#, omegaSquared] &];
  roots = distinctUnderAssumptions[rootValues, assumptions];
  emit[prefix <> "_Q3_DISTINCT_ROOTS", roots];
  emit[prefix <> "_Q3_ROOT_COUNT", Length[roots]];
  emit[prefix <> "_ROOT_ORDERING", roots];
  Do[
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q3_SIGN",
      refineWith[Sign[roots[[rootIndex]]], assumptions]
    ],
    {rootIndex, Length[roots]}
  ];
  pairs = Subsets[Range[Length[roots]], {2}];
  coincidenceTests = Table[
    {
      pair,
      fullSimplifyWith[
        roots[[pair[[1]]]] == roots[[pair[[2]]]],
        assumptions
      ]
    },
    {pair, pairs}
  ];
  emit[prefix <> "_Q3_ROOT_COINCIDENCE_LOCI", coincidenceTests];
  runModeSet[prefix, matrix, wavevector, roots, assumptions];
  <|"Roots" -> roots, "CoincidenceTests" -> coincidenceTests|>
];

runRankStrata[
  prefix_String,
  data_Association,
  matrix_List,
  spectrumData_Association,
  modeData_List,
  assumptions_
] := Module[
  {
    dimension, wavevector, amplitudes, allIndices, stratumRecords,
    rootMatrix, genericRank, rowSelections, columnSelections, minors,
    equations, locus, allowed, intersects, currentPrefix,
    coincidenceRecords, gatheredRecords, uniqueRecords, record,
    pointVariables, pointFormula, pointInstances, point,
    pointEqualities, pointAssumptions, pointMatrix, pointWavevector,
    stratumPrefix, controlVariables
  },
  dimension = Length[data["Wavevector"]];
  wavevector = data["Wavevector"];
  amplitudes = data["Amplitudes"];
  allIndices = Range[dimension];
  coincidenceRecords = spectrumData["CoincidenceRecords"];
  stratumRecords = {};

  Do[
    rootMatrix = modeData[[rootIndex, "RootMatrix"]];
    genericRank = modeData[[rootIndex, "Rank"]];
    rowSelections = Subsets[allIndices, {genericRank}];
    columnSelections = Subsets[allIndices, {genericRank}];
    minors = If[
      genericRank == 0,
      {Apply[Times, {}]},
      Flatten[
        Table[
          Factor[Together[Det[rootMatrix[[rows, columns]]]]],
          {rows, rowSelections},
          {columns, columnSelections}
        ]
      ]
    ];
    equations = (# == 0) & /@ minors;
    locus = realLocusSolve[equations, wavevector];
    allowed = reduceAllowed[And @@ equations, wavevector, assumptions];
    intersects = Not[SameQ[allowed, False]];
    currentPrefix = rootPrefix[prefix, rootIndex] <> "_Q8_RANK_DROP";
    emit[currentPrefix <> "_MINORS", minors];
    emit[currentPrefix <> "_OPERANDS", {equations, assumptions}];
    emitSolverObject[currentPrefix <> "_LOCUS", locus];
    emit[currentPrefix <> "_ALLOWED_INTERSECTION", allowed];
    emit[currentPrefix <> "_INTERSECTION_TEST", intersects];
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q8_ROOT_COINCIDENCE_LOCI",
      coincidenceRecords
    ];
    If[intersects,
      stratumRecords = Append[
        stratumRecords,
        <|
          "Source" -> {rankDrop, rootIndex},
          "Equations" -> equations,
          "Allowed" -> allowed
        |>
      ]
    ],
    {rootIndex, Length[modeData]}
  ];

  Do[
    If[TrueQ[coincidenceRecords[[recordIndex, "Intersects"]]],
      stratumRecords = Append[
        stratumRecords,
        <|
          "Source" -> {
            rootCoincidence,
            coincidenceRecords[[recordIndex, "Pair"]]
          },
          "Equations" -> {coincidenceRecords[[recordIndex, "Equation"]]},
          "Allowed" -> coincidenceRecords[[recordIndex, "Allowed"]]
        |>
      ]
    ],
    {recordIndex, Length[coincidenceRecords]}
  ];

  gatheredRecords = GatherBy[
    stratumRecords,
    ToString[#1["Allowed"], InputForm, PageWidth -> Infinity] &
  ];
  uniqueRecords = Map[
    Function[group,
      <|
        "Sources" -> Lookup[group, "Source"],
        "Equations" -> First[group]["Equations"],
        "Allowed" -> First[group]["Allowed"]
      |>
    ],
    gatheredRecords
  ];
  emit[prefix <> "_Q8_ALLOWED_STRATA", uniqueRecords];

  controlVariables = Select[
    {sRho, coefficientScale},
    Not[FreeQ[assumptions, #]] &
  ];
  pointVariables = Join[
    wavevector,
    amplitudes,
    {rhoBr, muR, braneDimension},
    controlVariables
  ];
  Do[
    record = uniqueRecords[[stratumIndex]];
    pointFormula = And[
      And @@ record["Equations"],
      assumptions,
      braneDimension == dimension
    ];
    pointInstances = FindInstance[
      pointFormula,
      pointVariables,
      Reals,
      1
    ];
    If[Length[pointInstances] == 0, Quit[2]];
    point = First[pointInstances];
    stratumPrefix = prefix <> "_STRATUM" <> ToString[stratumIndex];
    emit[stratumPrefix <> "_Q8_SOURCE", record["Sources"]];
    emit[stratumPrefix <> "_Q8_POINT", point];
    pointEqualities = And @@ (point /. Rule[left_, right_] :> left == right);
    pointAssumptions = And[assumptions, pointEqualities];
    pointMatrix = Map[
      fullSimplifyWith[#, pointAssumptions] &,
      matrix /. point,
      {2}
    ];
    pointWavevector = wavevector /. point;
    pointSpectrum[
      stratumPrefix,
      pointMatrix,
      pointWavevector,
      pointAssumptions
    ],
    {stratumIndex, Length[uniqueRecords]}
  ];
];

runPackage[package_String, dimension_Integer] := Module[
  {
    prefix, data, assumptions, expandedLagrangian, eom, matrixData,
    matrixA, matrixB, matrixResidual, matrixEntryRatio, spectrumData,
    modeData, scalingData
  },
  prefix = packagePrefix[package, dimension];
  data = buildPackage[package, dimension];
  assumptions = jointAssumptionExpression[
    package,
    data["Wavevector"],
    data["Amplitudes"]
  ];
  expandedLagrangian = Expand[data["Lagrangian"]];

  emit[prefix <> "_PREMISE_FIELD_CONTENT", suppliedInPlaneField[
    data["Fields"], braneDimension, separateSector[h]
  ]];
  emit[prefix <> "_PREMISE_BACKGROUND_STATE", suppliedBackground[v0 == 0]];
  emit[prefix <> "_PREMISE_DISSIPATION", suppliedNoDissipation[
    data["Lagrangian"]
  ]];
  emit[prefix <> "_PREMISE_RESPONSE", suppliedQuadraticResponse[
    data["Lagrangian"]
  ]];
  emit[prefix <> "_PREMISE_BASELINE_STIFFNESS", data["CurlStiffness"]];
  emit[prefix <> "_JOINT_ASSUMPTIONS", assumptions];
  emit[prefix <> "_ACTION_CONTROL", data["ControlData"]];
  emit[prefix <> "_SYMBOLIC_SIMPLIFIERS", {
    Together, Cancel, Factor, FullSimplify
  }];

  emit[prefix <> "_Q1_LAGRANGIAN", expandedLagrangian];
  eom = eulerLagrangeSystem[data];
  emit[prefix <> "_Q1_EULER_LAGRANGE_SYSTEM", Thread[eom == 0]];

  matrixData = deriveMatrices[data, eom, assumptions];
  matrixA = matrixData["MatrixA"];
  matrixB = matrixData["MatrixB"];
  matrixResidual = Map[
    fullSimplifyWith[#, assumptions] &,
    matrixA - matrixB,
    {2}
  ];
  matrixEntryRatio = fullSimplifyWith[
    Cancel[Together[matrixA[[1, 1]]/matrixB[[1, 1]]]],
    assumptions
  ];
  emit[prefix <> "_Q2_MATRIX_A", matrixA];
  emit[prefix <> "_Q2_MATRIX_B", matrixB];
  emit[prefix <> "_Q2_MATRIX_RESIDUAL", matrixResidual];
  emit[prefix <> "_Q2_MATRIX_ENTRY_RATIO", matrixEntryRatio];
  emit[prefix <> "_Q2_ROUTE_RESIDUAL_SCOPE", codingConsistencyOnly];
  emit[prefix <> "_Q2_DOWNSTREAM_ROUTE", quadraticFormRoute];

  spectrumData = runSpectrum[
    prefix,
    matrixB,
    data["Wavevector"],
    assumptions
  ];
  modeData = runModeSet[
    prefix,
    matrixB,
    data["Wavevector"],
    spectrumData["Roots"],
    assumptions
  ];
  scalingData = runScaling[
    prefix,
    spectrumData["Roots"],
    data["Wavevector"],
    assumptions
  ];
  runDimensions[
    prefix,
    data,
    eom,
    matrixData,
    spectrumData,
    scalingData,
    assumptions
  ];
  runCurlComparison[prefix, dimension];
  runRankStrata[
    prefix,
    data,
    matrixB,
    spectrumData,
    modeData,
    assumptions
  ];
];

packageRuns = {
  {"MAIN", {2, 3, 4, 5}},
  {"XFORM_FULLGRAD", {3, 4}},
  {"XFORM_DIVONLY", {3, 4}},
  {"XFORM_SIGNFLIP", {3, 4}},
  {"XFORM_ANISO", {3, 4}},
  {"XCOEF_SCALE", {3}}
};

declaredRunPairs = Flatten[
  Table[
    Table[
      {packageRuns[[packageIndex, 1]], dimension},
      {dimension, packageRuns[[packageIndex, 2]]}
    ],
    {packageIndex, Length[packageRuns]}
  ],
  1
];
actualRunPairs = declaredRunPairs;
skippedRunPairs = Complement[declaredRunPairs, actualRunPairs];

emit["WL_S10_RUN_PAIRS", actualRunPairs];
emit["WL_S10_SKIPPED_PAIRS", skippedRunPairs];

totalRuntime = First@AbsoluteTiming[
  Do[
    runPackage[
      packageRuns[[packageIndex, 1]],
      dimension
    ],
    {packageIndex, Length[packageRuns]},
    {dimension, packageRuns[[packageIndex, 2]]}
  ]
];

emit["WL_S10_RUNTIME_SECONDS", totalRuntime];
localListingTag = "WL_S10_LOCAL_TAG_NAMES";
localTagNames = Append[DeleteDuplicates[localTagNames], localListingTag];
emit[localListingTag, localTagNames];
