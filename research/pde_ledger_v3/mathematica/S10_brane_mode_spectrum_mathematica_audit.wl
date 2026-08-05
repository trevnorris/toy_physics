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

realLocusSolve[formula_, variables_List] := Module[
  {result, messages},
  Block[{$MessageList = {}},
    result = Solve[formula, variables, Reals];
    messages = $MessageList;
  ];
  <|
    "Result" -> result,
    "SolveSvarsMessageFired" -> Not[FreeQ[
      messages,
      HoldPattern[MessageName[Solve, "svars"]]
    ]]
  |>
];

classifyReduceOutcome[result_] := Which[
  SameQ[result, False], decidedEmpty,
  Not[FreeQ[HoldComplete[result], _Reduce | $Failed | $Aborted]], undecided,
  True, decidedNonempty
];

intersectionValue[outcome_] := Switch[
  outcome,
  decidedEmpty, False,
  decidedNonempty, True,
  _, Indeterminate
];

braneDimension = Symbol["braneDimension"];
rhoBr = Symbol["rhoBr"];
muR = Symbol["muR"];
sRho = Symbol["sRho"];
coefficientScale = Symbol["coefficientScale"];
omegaSquared = Symbol["omegaSquared"];
lambdaScale = Symbol["lambdaScale"];
timeCoordinate = Symbol["t"];
phaseVariable = Symbol["phi"];

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
      lambdaScale > 0,
      Total[wavevector^2] > 0
    },
    Element[#, Reals] & /@ wavevector,
    Element[#, Reals] & /@ amplitudes,
    {
      Element[braneDimension, Integers],
      Element[lambdaScale, Reals],
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
    amplitudes, wavevector, frequencyBranch, phaseExpression,
    phaseFunction, rules
  },
  dimension = Length[data["Fields"]];
  dummyArguments = Table[Unique["ansatzArgument"], {dimension + 1}];
  dummyCoordinates = Take[dummyArguments, dimension];
  dummyTime = Last[dummyArguments];
  amplitudes = data["Amplitudes"];
  wavevector = data["Wavevector"];
  frequencyBranch = Sqrt[omegaSquared];
  phaseExpression = Total[wavevector dummyCoordinates] -
    frequencyBranch dummyTime;
  phaseFunction = Function[
    Evaluate[dummyArguments],
    Evaluate[phaseExpression]
  ];
  rules = MapThread[
    Rule[
      #1,
      Function[
        Evaluate[dummyArguments],
        Evaluate[#2 Cos[phaseExpression]]
      ]
    ] &,
    {data["FieldHeads"], amplitudes}
  ];
  <|
    "Rules" -> rules,
    "FrequencyBranch" -> frequencyBranch,
    "PhaseFunction" -> phaseFunction
  |>
];

periodAverage[expression_, phase_, assumptions_] := Module[
  {
    averageWindow, normalization, phaseAverageSpecification,
    phaseExpression, integral, averagedExpression, conditions, value
  },
  averageWindow = {0, 2 Pi};
  normalization = 1/(averageWindow[[2]] - averageWindow[[1]]);
  phaseAverageSpecification = suppliedPhaseAverage[
    phaseVariable,
    averageWindow,
    normalization
  ];
  phaseExpression = expression /. {
    -phase -> -phaseVariable,
    phase -> phaseVariable
  };
  integral = Integrate[
    phaseExpression,
    Prepend[averageWindow, phaseVariable],
    Assumptions -> assumptions
  ];
  averagedExpression = fullSimplifyWith[
    normalization integral,
    assumptions
  ];
  conditions = DeleteDuplicates@Cases[
    averagedExpression,
    ConditionalExpression[_, condition_] :> condition,
    Infinity
  ];
  value = removeConditionalWrappers[averagedExpression];
  <|
    "Integral" -> integral,
    "Expression" -> averagedExpression,
    "Conditions" -> conditions,
    "Value" -> value,
    "Specification" -> phaseAverageSpecification
  |>
];

deriveMatrices[data_Association, eom_List, assumptions_] := Module[
  {
    ansatzData, rules, coordinates, amplitudes, wavevector, phase, commonFactor,
    eomOnAnsatz, strippedEquations, matrixA, lagrangianOnAnsatz,
    averageData, averagedLagrangian, matrixB, dimension
  },
  ansatzData = ansatzRules[data];
  rules = ansatzData["Rules"];
  coordinates = data["Coordinates"];
  amplitudes = data["Amplitudes"];
  wavevector = data["Wavevector"];
  dimension = Length[amplitudes];
  phase = Apply[
    ansatzData["PhaseFunction"],
    Join[coordinates, {timeCoordinate}]
  ];
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
  averageData = periodAverage[lagrangianOnAnsatz, phase, assumptions];
  averagedLagrangian = averageData["Value"];
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
    "AnsatzFrequencyBranch" -> ansatzData["FrequencyBranch"],
    "EOMOnAnsatz" -> eomOnAnsatz,
    "StrippedEquations" -> strippedEquations,
    "MatrixA" -> matrixA,
    "MatrixARoute" -> matrixRoute[equationsOfMotionRoute, matrixA],
    "LagrangianOnAnsatz" -> lagrangianOnAnsatz,
    "PeriodAverageIntegral" -> averageData["Integral"],
    "PeriodAverageExpression" -> averageData["Expression"],
    "PeriodAverageConditions" -> averageData["Conditions"],
    "PeriodAverageSpecification" -> averageData["Specification"],
    "AveragedLagrangian" -> averagedLagrangian,
    "MatrixB" -> matrixB,
    "MatrixBRoute" -> matrixRoute[quadraticFormRoute, matrixB]
  |>
];

distinctUnderAssumptions[expressions_List, assumptions_] :=
  DeleteDuplicates[
    expressions,
    TrueQ[fullSimplifyWith[#1 == #2, assumptions]] &
  ];

removeConditionalWrappers[expression_] :=
  expression /. ConditionalExpression[value_, condition_] :> value;

solutionRootRecords[expression_] := Module[{walk},
  walk[current_, inheritedCondition_] := Which[
    Head[current] === ConditionalExpression,
      walk[
        current[[1]],
        And[inheritedCondition, current[[2]]]
      ],
    MatchQ[current, HoldPattern[Rule[omegaSquared, _]]],
      Module[{rightHandSide, localConditions},
        rightHandSide = current[[2]];
        localConditions = Cases[
          rightHandSide,
          ConditionalExpression[_, condition_] :> condition,
          Infinity
        ];
        {
          <|
            "Root" -> removeConditionalWrappers[rightHandSide],
            "Condition" -> And[
              inheritedCondition,
              And @@ localConditions
            ]
          |>
        }
      ],
    ListQ[current],
      Flatten[walk[#, inheritedCondition] & /@ current],
    True,
      {}
  ];
  walk[expression, True]
];

conditionsForRoot[root_, records_List, assumptions_] := DeleteCases[
  DeleteDuplicates@Lookup[
    Select[
      records,
      TrueQ[
        fullSimplifyWith[#1["Root"] == root, assumptions]
      ] &
    ],
    "Condition"
  ],
  True
];

runSpectrum[
  prefix_String,
  matrixRouteObject_matrixRoute,
  wavevector_List,
  assumptions_
] := Module[
  {
    matrix,
    determinant, solutions, solutionRecords, rootsFromSolutions, filteredRoots,
    discardedRoots, roots, pairs,
    coincidenceRecords, coincidencePayload, pair, pairPrefix, equation,
    locusData, locus, allowed,
    outcome, intersects
  },
  matrix = matrixRouteObject[[2]];
  emit[prefix <> "_Q2_DOWNSTREAM_ROUTE", matrixRouteObject[[1]]];
  determinant = Factor[Det[matrix]];
  emit[prefix <> "_Q3_DETERMINANT", determinant];
  solutions = Solve[determinant == 0, omegaSquared];
  solutions = fullSimplifyWith[solutions, assumptions];
  emitSolverObject[prefix <> "_Q3_SOLUTIONS", solutions];
  solutionRecords = solutionRootRecords[solutions];
  rootsFromSolutions = Lookup[solutionRecords, "Root"];
  filteredRoots = Select[
    rootsFromSolutions,
    FreeQ[#, omegaSquared] &
  ];
  discardedRoots = Select[
    rootsFromSolutions,
    Not[FreeQ[#, omegaSquared]] &
  ];
  emit[
    prefix <> "_Q3_ROOT_CANDIDATE_COUNT_BEFORE_FILTER",
    Length[rootsFromSolutions]
  ];
  emit[prefix <> "_Q3_ROOTS_DISCARDED_BY_FILTER", discardedRoots];
  emit[
    prefix <> "_Q3_ROOT_CANDIDATE_COUNT_AFTER_FILTER",
    Length[filteredRoots]
  ];
  roots = distinctUnderAssumptions[filteredRoots, assumptions];
  emit[prefix <> "_Q3_DISTINCT_ROOTS", roots];
  emit[
    prefix <> "_Q3_ROOT_LIST_COUNTS",
    <|
      "FilteredCandidateCount" -> Length[filteredRoots],
      "DistinctRootCount" -> Length[roots]
    |>
  ];
  emit[prefix <> "_Q3_ROOT_COUNT", Length[roots]];
  emit[prefix <> "_ROOT_ORDERING", roots];
  Do[
    Module[{rootConditions},
      rootConditions = conditionsForRoot[
        roots[[rootIndex]],
        solutionRecords,
        assumptions
      ];
      emit[
        rootPrefix[prefix, rootIndex] <> "_Q3_SOLVER_CONDITIONS",
        rootConditions
      ];
    ];
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
    locusData = realLocusSolve[equation, wavevector];
    locus = locusData["Result"];
    allowed = reduceAllowed[equation, wavevector, assumptions];
    outcome = classifyReduceOutcome[allowed];
    intersects = intersectionValue[outcome];
    emit[pairPrefix <> "_OPERANDS", {equation, assumptions}];
    emitSolverObject[pairPrefix <> "_LOCUS", locus];
    emitLocal[
      localizeTag[pairPrefix <> "_LOCUS_SOLVE_SVARS_MESSAGE"],
      locusData["SolveSvarsMessageFired"]
    ];
    emit[pairPrefix <> "_ALLOWED_INTERSECTION", allowed];
    emit[pairPrefix <> "_INTERSECTION_OUTCOME", outcome];
    emit[pairPrefix <> "_INTERSECTION_TEST", intersects];
    coincidenceRecords = Append[
      coincidenceRecords,
      <|
        "Pair" -> pair,
        "Equation" -> equation,
        "Locus" -> locus,
        "Allowed" -> allowed,
        "IntersectionOutcome" -> outcome,
        "Intersects" -> intersects
      |>
    ],
    {pairIndex, Length[pairs]}
  ];
  coincidencePayload = If[
    Length[pairs] == 0,
    noRootPairs[Length[roots]],
    coincidenceRecords
  ];
  emit[prefix <> "_Q3_ROOT_COINCIDENCE_LOCI", coincidencePayload];

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
    "WavevectorProduct" -> wavevectorProduct,
    "Basis" -> basis,
    "BasisDots" -> basisDots,
    "BasisResiduals" -> basisResiduals
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
      emit[
        currentPrefix <> "_Q5_SCALING_EXPONENT",
        If[purePowerTest, exponentCandidate, notAPurePower]
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
          "RatioDomain" -> Not[zeroTest],
          "Ratio" -> undefinedRatio[root],
          "ExponentCandidate" -> notApplicable[undefinedRatio[root]],
          "PurePowerTest" -> notApplicable[undefinedRatio[root]]
        |>
      ];
      emit[
        currentPrefix <> "_Q5_SCALING_RATIO",
        undefinedRatio[root]
      ];
      emit[
        currentPrefix <> "_Q5_PURE_POWER_TEST",
        notApplicable[undefinedRatio[root]]
      ];
      emit[
        currentPrefix <> "_Q5_SCALING_EXPONENT_CANDIDATE",
        notApplicable[undefinedRatio[root]]
      ];
      emit[
        currentPrefix <> "_Q5_SCALING_EXPONENT",
        notApplicable[undefinedRatio[root]]
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
              "Constraints" -> {},
              "UnhandledHeads" -> {}
            |>,
          KeyExistsQ[symbolDimensions, current],
            <|
              "Dimension" -> symbolDimensions[current],
              "AddendDimensions" -> {},
              "Constraints" -> {},
              "UnhandledHeads" -> {}
            |>,
          MemberQ[fieldHeads, Head[current]],
            <|
              "Dimension" -> fieldDimension,
              "AddendDimensions" -> {},
              "Constraints" -> {},
              "UnhandledHeads" -> {}
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
              "Constraints" -> {},
              "UnhandledHeads" -> {}
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
              "Constraints" -> constraints,
              "UnhandledHeads" -> DeleteDuplicates@Flatten[
                Lookup[children, "UnhandledHeads"],
                1
              ]
            |>,
          Head[current] === Times,
            children = walk /@ (List @@ current);
            <|
              "Dimension" -> Total[Lookup[children, "Dimension"]],
              "AddendDimensions" ->
                Flatten[Lookup[children, "AddendDimensions"], 1],
              "Constraints" ->
                Flatten[Lookup[children, "Constraints"], 1],
              "UnhandledHeads" -> DeleteDuplicates@Flatten[
                Lookup[children, "UnhandledHeads"],
                1
              ]
            |>,
          Head[current] === Power,
            children = walk[First[current]];
            <|
              "Dimension" -> Last[current] children["Dimension"],
              "AddendDimensions" -> children["AddendDimensions"],
              "Constraints" -> children["Constraints"],
              "UnhandledHeads" -> children["UnhandledHeads"]
            |>,
          Head[current] === Abs,
            walk[First[current]],
          Head[current] === ConditionalExpression,
            walk[First[current]],
          True,
            <|
              "Dimension" -> dimensionIndeterminate[
                HoldForm[Head[current]]
              ],
              "AddendDimensions" -> {},
              "Constraints" -> {},
              "UnhandledHeads" -> {HoldForm[Head[current]]}
            |>
        ]
      ];
      walk[expression]
    ]
  ]
];

resolveDimensionAnalysis[analysis_Association, rules_List, assumptions_] := Module[
  {
    dimension, termDimensions, pairwiseDifferences, homogeneity,
    unhandledHeads
  },
  dimension = fullSimplifyWith[analysis["Dimension"] /. rules, assumptions];
  termDimensions = fullSimplifyWith[
    analysis["AddendDimensions"] /. rules,
    assumptions
  ];
  unhandledHeads = analysis["UnhandledHeads"];
  pairwiseDifferences = Function[group,
      fullSimplifyWith[#[[1]] - #[[2]], assumptions] & /@
        Subsets[group, {2}]
    ] /@ termDimensions;
  homogeneity = If[
    Length[unhandledHeads] > 0,
    False,
    And @@ Flatten[
      Function[group,
        TrueQ[
          fullSimplifyWith[
            # == ConstantArray[0, Length[#]],
            assumptions
          ]
        ] & /@ group
      ] /@ pairwiseDifferences
    ]
  ];
  <|
    "Dimension" -> dimension,
    "AddendDimensions" -> termDimensions,
    "PairwiseDifferences" -> pairwiseDifferences,
    "Homogeneity" -> homogeneity,
    "UnhandledHeads" -> unhandledHeads
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
      "PairwiseDifferences" -> {},
      "Homogeneity" -> True,
      "UnhandledHeads" -> {}
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
  modeData_List,
  curlData_Association,
  rankData_Association,
  assumptions_
] := Module[
  {
    dimensionNames, rhoDimension, muDimension, sRhoDimension,
    scaleDimension, coefficientDimensionMap, actionCoefficientSymbols,
    controlSymbols, unknowns, controlEquations, energyDensityDimension,
    expandedLagrangian, actionTerms, analyzer, actionTermAnalyses,
    actionTermDimensions, dimensionEquations, dimensionSolutions,
    dimensionEquationResiduals, dimensionCoefficientArrays,
    dimensionAugmentedMatrix, independentDimensionEquationCount,
    unknownCoefficientDimensionCount, dimensionEquationCountDifference,
    dimensionSystemDetermination, actionHomogeneityVacuity,
    concreteDimension, braneDimensionRelation,
    dimensionSolutionsSpecialized, dimensionRules, inertialAnalyses,
    stiffnessAnalyses, solvedActionResolved,
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
  Do[
    AssociateTo[
      coefficientDimensionMap,
      curlData["GradientSymbols"][[row, column]] -> {0, 0, 0}
    ],
    {row, Length[curlData["GradientSymbols"]]},
    {column, Length[curlData["GradientSymbols"]]}
  ];
  analyzer = makeDimensionAnalyzer[
    data["FieldHeads"],
    coefficientDimensionMap
  ];
  emit[
    prefix <> "_PREMISE_FIELD_DIMENSION",
    suppliedFieldDimension[
      u,
      analyzer[First[data["Fields"]]]["Dimension"]
    ]
  ];

  concreteDimension = Length[data["Wavevector"]];
  braneDimensionRelation = braneDimension == concreteDimension;
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
  dimensionEquationResiduals = dimensionEquations /. Equal -> Subtract;
  dimensionCoefficientArrays = CoefficientArrays[
    dimensionEquationResiduals,
    unknowns
  ];
  dimensionAugmentedMatrix = MapThread[
    Append,
    {
      Normal[dimensionCoefficientArrays[[2]]],
      Normal[dimensionCoefficientArrays[[1]]]
    }
  ];
  independentDimensionEquationCount = MatrixRank[
    dimensionAugmentedMatrix
  ];
  unknownCoefficientDimensionCount = Length[unknowns];
  dimensionEquationCountDifference =
    independentDimensionEquationCount - unknownCoefficientDimensionCount;
  dimensionSystemDetermination = Which[
    dimensionEquationCountDifference > 0, overdetermined,
    dimensionEquationCountDifference == 0, exactlyDetermined,
    True, underdetermined
  ];
  actionHomogeneityVacuity = TrueQ[
    dimensionEquationCountDifference <= 0
  ];
  dimensionSolutionsSpecialized = fullSimplifyWith[
    dimensionSolutions /. braneDimension -> concreteDimension,
    assumptions
  ];
  emit[
    prefix <> "_Q6_BRANE_DIMENSION_RELATION",
    braneDimensionRelation
  ];
  emit[prefix <> "_Q6_ENERGY_DENSITY_DIMENSION", energyDensityDimension];
  emit[prefix <> "_Q6_DIMENSION_PREMISES", data["DimensionPremises"]];
  emit[
    prefix <> "_Q6_ACTION_TERM_DIMENSIONS_SYMBOLIC",
    actionTermDimensions
  ];
  emit[
    prefix <> "_Q6_ACTION_TERM_PAIRWISE_DIFFERENCES_SYMBOLIC",
    (#[[1]] - #[[2]]) & /@ Subsets[actionTermDimensions, {2}]
  ];
  emit[
    prefix <> "_Q6_ACTION_TERM_UNHANDLED_HEADS",
    Lookup[actionTermAnalyses, "UnhandledHeads"]
  ];
  emit[prefix <> "_Q6_DIMENSION_EQUATIONS", dimensionEquations];
  emit[
    prefix <> "_Q6_DIMENSION_EQUATION_INDEPENDENT_COUNT",
    independentDimensionEquationCount
  ];
  emit[
    prefix <> "_Q6_UNKNOWN_COEFFICIENT_DIMENSION_COUNT",
    unknownCoefficientDimensionCount
  ];
  emit[
    prefix <> "_Q6_DIMENSION_EQUATION_COUNT_DIFFERENCE",
    dimensionEquationCountDifference
  ];
  emit[
    prefix <> "_Q6_DIMENSION_SYSTEM_DETERMINATION",
    dimensionSystemDetermination
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_HOMOGENEITY_VACUITY",
    <|
      "EquationCountDifference" -> dimensionEquationCountDifference,
      "Vacuous" -> actionHomogeneityVacuity
    |>
  ];
  emitSolverObject[prefix <> "_Q6_DIMENSION_SOLUTIONS", dimensionSolutions];
  emitSolverObject[
    prefix <> "_Q6_DIMENSION_SOLUTIONS_SPECIALIZED",
    dimensionSolutionsSpecialized
  ];
  dimensionRules = If[Length[dimensionSolutions] > 0, First[dimensionSolutions], {}];
  solvedActionResolved = resolvedExpressionDimension[
    data["Lagrangian"],
    analyzer,
    dimensionRules,
    assumptions
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_DIMENSION",
    solvedActionResolved["Dimension"]
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_TERM_DIMENSIONS",
    solvedActionResolved["AddendDimensions"]
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_TERM_PAIRWISE_DIFFERENCES",
    solvedActionResolved["PairwiseDifferences"]
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_TERM_HOMOGENEITY",
    solvedActionResolved["Homogeneity"]
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_TERM_UNHANDLED_HEADS",
    solvedActionResolved["UnhandledHeads"]
  ];

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
    prefix <> "_Q6_INERTIAL_COEFFICIENT_DIMENSIONS_SPECIALIZED",
    Lookup[inertialResolved, "Dimension"] /.
      braneDimension -> concreteDimension
  ];
  emit[
    prefix <> "_Q6_INERTIAL_COEFFICIENT_TERM_DIMENSIONS",
    Lookup[inertialResolved, "AddendDimensions"]
  ];
  emit[
    prefix <> "_Q6_INERTIAL_COEFFICIENT_PAIRWISE_DIFFERENCES",
    Lookup[inertialResolved, "PairwiseDifferences"]
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_INERTIAL_COEFFICIENT_HOMOGENEITY",
    Lookup[inertialResolved, "Homogeneity"]
  ];
  emit[
    prefix <> "_Q6_INERTIAL_COEFFICIENT_UNHANDLED_HEADS",
    Lookup[inertialResolved, "UnhandledHeads"]
  ];
  emit[prefix <> "_Q6_STIFFNESS_COEFFICIENTS", data["StiffnessCoefficients"]];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_DIMENSIONS",
    Lookup[stiffnessResolved, "Dimension"]
  ];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_DIMENSIONS_SPECIALIZED",
    Lookup[stiffnessResolved, "Dimension"] /.
      braneDimension -> concreteDimension
  ];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_TERM_DIMENSIONS",
    Lookup[stiffnessResolved, "AddendDimensions"]
  ];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_PAIRWISE_DIFFERENCES",
    Lookup[stiffnessResolved, "PairwiseDifferences"]
  ];
  emit[
    prefix <> "_Q6_SOLVED_ACTION_STIFFNESS_COEFFICIENT_HOMOGENEITY",
    Lookup[stiffnessResolved, "Homogeneity"]
  ];
  emit[
    prefix <> "_Q6_STIFFNESS_COEFFICIENT_UNHANDLED_HEADS",
    Lookup[stiffnessResolved, "UnhandledHeads"]
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
      rootPrefix[prefix, rootIndex] <>
        "_Q6_RATE_COEFFICIENT_PAIRWISE_DIFFERENCES",
      rateResolved[[rootIndex, "PairwiseDifferences"]]
    ];
    emit[
      rootPrefix[prefix, rootIndex] <>
        "_Q6_DERIVED_RATE_COEFFICIENT_HOMOGENEITY",
      rateResolved[[rootIndex, "Homogeneity"]]
    ];
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q6_RATE_COEFFICIENT_UNHANDLED_HEADS",
      rateResolved[[rootIndex, "UnhandledHeads"]]
    ],
    {rootIndex, Length[rateCoefficientExpressions]}
  ];

  auditGroups = <|
    "EOM" -> eom,
    "AVERAGED_LAGRANGIAN" -> {matrixData["AveragedLagrangian"]},
    "MATRIX_A_ENTRIES" -> Flatten[matrixData["MatrixA"]],
    "MATRIX_B_ENTRIES" -> Flatten[matrixData["MatrixB"]],
    "DETERMINANT" -> {spectrumData["Determinant"]},
    "ROOTS" -> spectrumData["Roots"],
    "SCALED_ROOTS" -> scalingData["ScaledRoots"],
    "SCALING_RATIOS" -> scalingData["Ratios"],
    "N1_MATRIX_ENTRIES" -> Flatten[
      Lookup[modeData, "RootMatrix", {}]
    ],
    "N5_WAVEVECTOR_PRODUCT_ENTRIES" -> Flatten[
      Lookup[modeData, "WavevectorProduct", {}]
    ],
    "N6_BASIS_ENTRIES" -> Flatten[Lookup[modeData, "Basis", {}]],
    "N6_BASIS_DOTS" -> Flatten[Lookup[modeData, "BasisDots", {}]],
    "N6_BASIS_RESIDUAL_ENTRIES" -> Flatten[
      Lookup[modeData, "BasisResiduals", {}]
    ],
    "Q8_MINOR_EXPRESSIONS" -> Flatten[rankData["Minors"]]
  |>;
  If[Length[curlData["Expressions"]] > 0,
    AssociateTo[auditGroups, "Q7_EXPRESSIONS" -> curlData["Expressions"]]
  ];
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
      prefix <> "_Q6_" <> groupNames[[groupIndex]] <>
        "_PAIRWISE_DIFFERENCES",
      Lookup[groupResolved, "PairwiseDifferences"]
    ];
    emit[
      prefix <> "_Q6_DERIVED_" <> groupNames[[groupIndex]] <>
        "_HOMOGENEITY",
      Lookup[groupResolved, "Homogeneity"]
    ];
    emit[
      prefix <> "_Q6_" <> groupNames[[groupIndex]] <> "_UNHANDLED_HEADS",
      Lookup[groupResolved, "UnhandledHeads"]
    ],
    {groupIndex, Length[groupNames]}
  ];

  <|
    "DimensionSolutions" -> dimensionSolutions,
    "DimensionRules" -> dimensionRules
  |>
];

runCurlComparison[prefix_String, data_Association] := Module[
  {
    dimension, gradientSymbols, gradientRules, actionStiffness,
    curlVector, curlNorm, comparisonExpressions
  },
  dimension = Length[data["Gradient"]];
  gradientSymbols = {};
  comparisonExpressions = {};
  If[dimension == 3,
    gradientSymbols = Table[
      Symbol["g" <> ToString[row] <> "x" <> ToString[column]],
      {row, 3},
      {column, 3}
    ];
    gradientRules = Thread[
      Flatten[data["Gradient"]] -> Flatten[gradientSymbols]
    ];
    actionStiffness = Expand[
      data["StiffnessDensity"] /. gradientRules
    ];
    curlVector = {
      gradientSymbols[[2, 3]] - gradientSymbols[[3, 2]],
      gradientSymbols[[3, 1]] - gradientSymbols[[1, 3]],
      gradientSymbols[[1, 2]] - gradientSymbols[[2, 1]]
    };
    curlNorm = Expand[curlVector . curlVector];
    emit[prefix <> "_Q7_PACKAGE_STIFFNESS_DENSITY", actionStiffness];
    emit[prefix <> "_Q7_ORDINARY_CURL_NORM", curlNorm];
    emit[
      prefix <> "_Q7_PACKAGE_STIFFNESS_VS_ORDINARY_CURL_RESIDUAL",
      Expand[actionStiffness - curlNorm]
    ];
    comparisonExpressions = {
      actionStiffness,
      curlNorm,
      Expand[actionStiffness - curlNorm]
    };
  ];
  <|
    "GradientSymbols" -> gradientSymbols,
    "Expressions" -> comparisonExpressions
  |>
];

pointSpectrum[
  prefix_String,
  matrix_List,
  wavevector_List,
  parameterVariables_List,
  assumptions_
] := Module[
  {
    determinant, solutions, solutionRecords, rootValues, filteredRoots,
    discardedRoots,
    roots, pairs, coincidenceTests, pair, pairPrefix,
    equation, locus,
    allowed, outcome, intersects, locusData
  },
  determinant = Factor[Det[matrix]];
  emit[prefix <> "_Q3_DETERMINANT", determinant];
  solutions = fullSimplifyWith[
    Solve[determinant == 0, omegaSquared],
    assumptions
  ];
  emitSolverObject[prefix <> "_Q3_SOLUTIONS", solutions];
  solutionRecords = solutionRootRecords[solutions];
  rootValues = Lookup[solutionRecords, "Root"];
  filteredRoots = Select[rootValues, FreeQ[#, omegaSquared] &];
  discardedRoots = Select[rootValues, Not[FreeQ[#, omegaSquared]] &];
  emit[
    prefix <> "_Q3_ROOT_CANDIDATE_COUNT_BEFORE_FILTER",
    Length[rootValues]
  ];
  emit[prefix <> "_Q3_ROOTS_DISCARDED_BY_FILTER", discardedRoots];
  emit[
    prefix <> "_Q3_ROOT_CANDIDATE_COUNT_AFTER_FILTER",
    Length[filteredRoots]
  ];
  roots = distinctUnderAssumptions[filteredRoots, assumptions];
  emit[prefix <> "_Q3_DISTINCT_ROOTS", roots];
  emit[
    prefix <> "_Q3_ROOT_LIST_COUNTS",
    <|
      "FilteredCandidateCount" -> Length[filteredRoots],
      "DistinctRootCount" -> Length[roots]
    |>
  ];
  emit[prefix <> "_Q3_ROOT_COUNT", Length[roots]];
  emit[prefix <> "_ROOT_ORDERING", roots];
  Do[
    Module[{rootConditions},
      rootConditions = conditionsForRoot[
        roots[[rootIndex]],
        solutionRecords,
        assumptions
      ];
      emit[
        rootPrefix[prefix, rootIndex] <> "_Q3_SOLVER_CONDITIONS",
        rootConditions
      ];
    ];
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q3_SIGN",
      refineWith[Sign[roots[[rootIndex]]], assumptions]
    ],
    {rootIndex, Length[roots]}
  ];
  pairs = Subsets[Range[Length[roots]], {2}];
  coincidenceTests = {};
  Do[
    pair = pairs[[pairIndex]];
    pairPrefix = prefix <> "_ROOT" <> ToString[pair[[1]]] <>
      "_ROOT" <> ToString[pair[[2]]] <> "_Q3_COINCIDENCE";
    equation = Together[
      roots[[pair[[1]]]] - roots[[pair[[2]]]]
    ] == 0;
    locusData = realLocusSolve[equation, parameterVariables];
    locus = locusData["Result"];
    allowed = reduceAllowed[equation, parameterVariables, assumptions];
    outcome = classifyReduceOutcome[allowed];
    intersects = intersectionValue[outcome];
    emit[pairPrefix <> "_OPERANDS", {equation, assumptions}];
    emitSolverObject[pairPrefix <> "_PARAMETER_LOCUS", locus];
    emitLocal[
      localizeTag[
        pairPrefix <> "_PARAMETER_LOCUS_SOLVE_SVARS_MESSAGE"
      ],
      locusData["SolveSvarsMessageFired"]
    ];
    emit[pairPrefix <> "_ALLOWED_PARAMETER_REGION", allowed];
    emit[pairPrefix <> "_PARAMETER_INTERSECTION_OUTCOME", outcome];
    emit[pairPrefix <> "_PARAMETER_INTERSECTION_TEST", intersects];
    coincidenceTests = Append[
      coincidenceTests,
      <|
        "Pair" -> pair,
        "Equation" -> equation,
        "ParameterLocus" -> locus,
        "AllowedParameterRegion" -> allowed,
        "IntersectionOutcome" -> outcome,
        "Intersects" -> intersects
      |>
    ],
    {pairIndex, Length[pairs]}
  ];
  emit[
    prefix <> "_Q3_STRATUM_ROOT_COINCIDENCE_RECORDS",
    If[
      Length[pairs] == 0,
      noRootPairs[Length[roots]],
      coincidenceTests
    ]
  ];
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
    dimension, wavevector, allIndices, stratumRecords,
    rootMatrix, genericRank, rowSelections, columnSelections, minors,
    equations, locusData, locus, allowed, outcome, intersects, currentPrefix,
    coincidenceRecords, gatheredRecords, uniqueRecords, record,
    sameStratumRegionQ,
    pointVariables, pointFormula, pointInstances, point, pointLocated,
    pointSearchRecord, spectrumParameterVariables,
    parameterSpecializationRules, parameterSpecializations,
    pointAssumptions, pointMatrix, pointWavevector,
    pointProjectionVariables, pointProjectionFormula, stratumPrefix,
    controlVariables, allMinors
  },
  dimension = Length[data["Wavevector"]];
  wavevector = data["Wavevector"];
  allIndices = Range[dimension];
  coincidenceRecords = spectrumData["CoincidenceRecords"];
  stratumRecords = {};
  allMinors = {};

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
    allMinors = Join[allMinors, minors];
    locusData = realLocusSolve[equations, wavevector];
    locus = locusData["Result"];
    allowed = reduceAllowed[And @@ equations, wavevector, assumptions];
    outcome = classifyReduceOutcome[allowed];
    intersects = intersectionValue[outcome];
    currentPrefix = rootPrefix[prefix, rootIndex] <> "_Q8_RANK_DROP";
    emit[currentPrefix <> "_MINORS", minors];
    emit[currentPrefix <> "_OPERANDS", {equations, assumptions}];
    emitSolverObject[currentPrefix <> "_LOCUS", locus];
    emitLocal[
      localizeTag[currentPrefix <> "_LOCUS_SOLVE_SVARS_MESSAGE"],
      locusData["SolveSvarsMessageFired"]
    ];
    emit[currentPrefix <> "_ALLOWED_INTERSECTION", allowed];
    emit[currentPrefix <> "_INTERSECTION_OUTCOME", outcome];
    emit[currentPrefix <> "_INTERSECTION_TEST", intersects];
    emit[
      rootPrefix[prefix, rootIndex] <> "_Q8_ROOT_COINCIDENCE_LOCI",
      coincidenceRecords
    ];
    If[TrueQ[intersects],
      stratumRecords = Append[
        stratumRecords,
        <|
          "Source" -> {rankDrop, rootIndex},
          "Equations" -> equations,
          "Allowed" -> allowed,
          "IntersectionOutcome" -> outcome
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

  sameStratumRegionQ = Function[{left, right},
    SameQ[
      Reduce[
        And[
          assumptions,
          Xor[left["Allowed"], right["Allowed"]]
        ],
        wavevector,
        Reals
      ],
      False
    ]
  ];
  gatheredRecords = Gather[
    stratumRecords,
    sameStratumRegionQ
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
  spectrumParameterVariables = Join[{rhoBr, muR}, controlVariables];
  pointVariables = wavevector;
  Do[
    record = uniqueRecords[[stratumIndex]];
    pointFormula = And[
      And @@ record["Equations"],
      assumptions,
      braneDimension == dimension
    ];
    pointProjectionVariables = Join[
      data["Amplitudes"],
      {rhoBr, muR, braneDimension, lambdaScale},
      controlVariables
    ];
    pointProjectionFormula = Reduce[
      Apply[Exists, {pointProjectionVariables, pointFormula}],
      wavevector,
      Reals
    ];
    stratumPrefix = prefix <> "_STRATUM" <> ToString[stratumIndex];
    pointInstances = Quiet[
      FindInstance[pointProjectionFormula, pointVariables, Reals, 1]
    ];
    pointLocated = ListQ[pointInstances] && Length[pointInstances] > 0;
    point = If[
      pointLocated,
      First[pointInstances],
      pointUnavailable[pointInstances]
    ];
    parameterSpecializationRules = Cases[
      Flatten[{point}],
      Rule[symbol_, value_] /; MemberQ[spectrumParameterVariables, symbol]
    ];
    parameterSpecializations = If[
      Length[parameterSpecializationRules] == 0,
      noParametersSpecialized[spectrumParameterVariables],
      parameterSpecializationRules
    ];
    pointSearchRecord = <|
      "AllowedRegion" -> record["Allowed"],
      "ProjectionFormula" -> pointProjectionFormula,
      "Variables" -> pointVariables,
      "Output" -> pointInstances,
      "Located" -> pointLocated
    |>;
    emit[stratumPrefix <> "_Q8_SOURCE", record["Sources"]];
    emit[stratumPrefix <> "_Q8_POINT_SEARCH", pointSearchRecord];
    emit[stratumPrefix <> "_Q8_POINT", point];
    emit[
      stratumPrefix <> "_Q8_PARAMETER_SPECIALIZATIONS",
      parameterSpecializations
    ];
    If[pointLocated,
      pointAssumptions = And[
          assumptions,
          And @@ record["Equations"],
          braneDimension == dimension
        ] /. point;
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
        spectrumParameterVariables,
        pointAssumptions
      ]
    ],
    {stratumIndex, Length[uniqueRecords]}
  ];
  <|"Minors" -> allMinors|>
];

runPackage[package_String, dimension_Integer] := Module[
  {
    prefix, data, assumptions, expandedLagrangian, eom, matrixData,
    matrixA, matrixB, matrixResidual, matrixEntryRatio,
    downstreamMatrixRoute, downstreamMatrix, spectrumData, modeData,
    scalingData, curlData, rankData
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

  emit[prefix <> "_Q1_LAGRANGIAN", expandedLagrangian];
  eom = eulerLagrangeSystem[data];
  emit[prefix <> "_Q1_EULER_LAGRANGE_SYSTEM", Thread[eom == 0]];

  matrixData = deriveMatrices[data, eom, assumptions];
  emit[
    prefix <> "_PREMISE_PERIOD_AVERAGE",
    matrixData["PeriodAverageSpecification"]
  ];
  emitLocal[
    localizeTag[prefix <> "_ANSATZ_FREQUENCY_BRANCH"],
    matrixData["AnsatzFrequencyBranch"]
  ];
  matrixA = matrixData["MatrixA"];
  matrixB = matrixData["MatrixB"];
  downstreamMatrixRoute = matrixData["MatrixBRoute"];
  downstreamMatrix = downstreamMatrixRoute[[2]];
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
  emit[
    prefix <> "_Q2_PERIOD_AVERAGE_CONDITIONS",
    matrixData["PeriodAverageConditions"]
  ];
  emit[prefix <> "_Q2_MATRIX_B", matrixB];
  emit[prefix <> "_Q2_MATRIX_RESIDUAL", matrixResidual];
  emit[prefix <> "_Q2_MATRIX_ENTRY_RATIO", matrixEntryRatio];
  emit[prefix <> "_Q2_ROUTE_RESIDUAL_SCOPE", codingConsistencyOnly];

  spectrumData = runSpectrum[
    prefix,
    downstreamMatrixRoute,
    data["Wavevector"],
    assumptions
  ];
  modeData = runModeSet[
    prefix,
    downstreamMatrix,
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
  curlData = runCurlComparison[prefix, data];
  rankData = runRankStrata[
    prefix,
    data,
    downstreamMatrix,
    spectrumData,
    modeData,
    assumptions
  ];
  runDimensions[
    prefix,
    data,
    eom,
    matrixData,
    spectrumData,
    scalingData,
    modeData,
    curlData,
    rankData,
    assumptions
  ];
  actualRunPairs = Append[actualRunPairs, {package, dimension}];
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
actualRunPairs = {};

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

skippedRunPairs = Complement[declaredRunPairs, actualRunPairs];
emit["WL_S10_RUN_PAIRS", actualRunPairs];
emit["WL_S10_SKIPPED_PAIRS", skippedRunPairs];
emit["WL_S10_RUNTIME_SECONDS", totalRuntime];
localListingTag = "WL_S10_LOCAL_TAG_NAMES";
localTagNames = Append[DeleteDuplicates[localTagNames], localListingTag];
emit[localListingTag, localTagNames];
emit["WL_S10_EMITTED_TAG_COUNT", emittedTagCount + 1];
