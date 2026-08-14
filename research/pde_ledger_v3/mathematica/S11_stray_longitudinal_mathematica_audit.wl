$HistoryLength = 0;
ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};

tokenProvedTrue = "PROVED_TRUE";
tokenProvedFalse = "PROVED_FALSE";
tokenUndecided = "UNDECIDED";
tokenAdmissible = "ADMISSIBLE";
tokenExcluded = "EXCLUDED";
tokenPositive = "POSITIVE";
tokenNegative = "NEGATIVE";
tokenZero = "ZERO";
tokenNotApplicable = "NOT_APPLICABLE";
tokenUndefinedRatio = "UNDEFINED_RATIO";
tokenNotPurePower = "NOT_A_PURE_POWER";
tokenOverDetermined = "OVER_DETERMINED";
tokenExactlyDetermined = "EXACTLY_DETERMINED";
tokenUnderDetermined = "UNDER_DETERMINED";
tokenConstant = "CONSTANT";
tokenVaries = "VARIES";
tokenCompleteCoverage = "COMPLETE_COVERAGE";
tokenIncompleteCoverage = "INCOMPLETE_COVERAGE";
tokenRankDrop = "RANK_DROP";
tokenStackedDrop = "STACKED_DROP";
tokenRootCoincidence = "ROOT_COINCIDENCE";

(* Internal emission keys are mapped to the shared §8 names only here. *)
cellEmission[package_String, dimension_Integer, quantity_String] :=
  HoldComplete[CellObject, package, dimension, quantity];
localCellEmission[package_String, dimension_Integer, quantity_String] :=
  HoldComplete[LocalCellObject, package, dimension, quantity];
runEmission[quantity_String] := HoldComplete[RunObject, quantity];
localInventoryEmission[] := HoldComplete[LocalInventoryObject];

standardEmissionName[HoldComplete[CellObject, package_String, dimension_Integer,
    quantity_String]] := "WL_S11_" <> package <> "_D" <> ToString[dimension] <>
  "_" <> quantity;
standardEmissionName[HoldComplete[LocalCellObject, package_String,
    dimension_Integer, quantity_String]] := "WL_S11_LOCAL_" <> package <> "_D" <>
  ToString[dimension] <> "_" <> quantity;
standardEmissionName[HoldComplete[RunObject, quantity_String]] :=
  "WL_S11_" <> quantity;
standardEmissionName[HoldComplete[LocalInventoryObject]] :=
  "WL_S11_LOCAL_TAG_NAMES";

emittedNames = <||>;
localNames = {};
emit[internalTag_, payload_] := Module[{emittedName, rendered, stream},
  emittedName = standardEmissionName[internalTag];
  If[! StringQ[emittedName], Quit[90]];
  If[KeyExistsQ[emittedNames, emittedName], Quit[91]];
  AssociateTo[emittedNames, emittedName -> True];
  If[StringStartsQ[emittedName, "WL_S11_LOCAL_"],
    AppendTo[localNames, emittedName]];
  rendered = ToString[payload, InputForm, PageWidth -> Infinity];
  stream = First[$Output];
  WriteString[stream, emittedName <> ": " <> rendered <> "\n"];
  Flush[stream];
];

emitCell[package_, dimension_, quantity_, payload_] :=
  emit[cellEmission[package, dimension, quantity], payload];
emitLocalCell[package_, dimension_, quantity_, payload_] :=
  emit[localCellEmission[package, dimension, quantity], payload];

engineSimplify[expression_, assumptions_] :=
  FullSimplify[expression, Assumptions -> assumptions];
unrestrictedSimplify[expression_] := FullSimplify[expression];
zeroTest[assumptions_] := Function[value,
  TrueQ[engineSimplify[value == 0, assumptions]]];
assumedRank[matrix_, assumptions_] :=
  MatrixRank[matrix, ZeroTest -> zeroTest[assumptions]];
assumedNullSpace[matrix_, assumptions_] :=
  NullSpace[matrix, ZeroTest -> zeroTest[assumptions]];

relationResidual[relation_] := Which[
  TrueQ[relation], 0,
  SameQ[relation, False], 1,
  Head[relation] === Equal, Subtract @@ List @@ relation,
  True, relation
];

globalSymbolsIn[expression_] := DeleteDuplicates[Cases[
  expression, symbol_Symbol /; Context[symbol] === "Global`",
  {0, Infinity}]];

orderedGlobalSymbols[expression_] :=
  SortBy[globalSymbolsIn[expression], SymbolName];

truthStatus[testObject_] := Which[
  TrueQ[testObject], tokenProvedTrue,
  SameQ[testObject, False], tokenProvedFalse,
  True, tokenUndecided
];

admissibilityStatus[testObject_] := Which[
  TrueQ[testObject], tokenAdmissible,
  SameQ[testObject, False], tokenExcluded,
  True, tokenUndecided
];

fourWaySign[operand_, assumptions_] := Module[{positiveTest, negativeTest, zeroTestObject},
  positiveTest = engineSimplify[operand > 0, assumptions];
  negativeTest = engineSimplify[operand < 0, assumptions];
  zeroTestObject = engineSimplify[operand == 0, assumptions];
  Which[
    TrueQ[positiveTest], tokenPositive,
    TrueQ[negativeTest], tokenNegative,
    TrueQ[zeroTestObject], tokenZero,
    True, tokenUndecided
  ]
];

signPayload[operand_, assumptions_] := {
  "SIGN_TOKEN" -> fourWaySign[operand, assumptions],
  "OPERAND" -> operand
};

realExistenceTest[condition_] := Module[{variables},
  variables = orderedGlobalSymbols[condition];
  Quiet[If[variables === {}, unrestrictedSimplify[condition],
    With[{quantifiedVariables = variables, predicate = condition},
      Resolve[Exists[quantifiedVariables, predicate], Reals]]]]
];

exactProjectedInstance[condition_, projectedVariables_List] := Module[
  {allVariables, instances, instance, values},
  allVariables = DeleteDuplicates[Join[orderedGlobalSymbols[condition], projectedVariables]];
  instances = Quiet[If[allVariables === {},
    If[TrueQ[unrestrictedSimplify[condition]], {{}}, {}],
    With[{instanceVariables = allVariables, predicate = condition},
      FindInstance[predicate, instanceVariables, Reals, 1]]]];
  If[! ListQ[instances] || instances === {}, Return[tokenNotApplicable]];
  instance = First[instances];
  values = Table[variable /. instance, {variable, projectedVariables}];
  If[AnyTrue[MapThread[SameQ[#1, #2] &, {values, projectedVariables}], TrueQ],
    tokenNotApplicable,
    MapThread[(#1 -> #2) &, {projectedVariables, values}]
  ]
];

plainBranchRules[branch_] := Cases[branch, _Rule, {1}] /.
  Rule[left_, ConditionalExpression[right_, _]] :> Rule[left, right];

branchConditions[branch_] := Cases[
  branch, Rule[_, ConditionalExpression[_, condition_]] :> condition, {1}];

branchDefiningEquations[branch_] :=
  (#[[1]] == #[[2]] &) /@ plainBranchRules[branch];

conditionTerms[condition_] :=
  If[Head[condition] === And, List @@ condition, {condition}];

branchRelations[branch_] := Join[
  branchDefiningEquations[branch],
  Flatten[conditionTerms /@ branchConditions[branch], 1]
];

polynomialLocusQ[residuals_List] := And @@ Map[
  Function[residual,
    PolynomialQ[residual, orderedGlobalSymbols[residual]]],
  residuals
];

canonicalLocus[residuals_List, variables_List] := Module[{},
  If[! polynomialLocusQ[residuals], Return[tokenNotApplicable]];
  Quiet[GroebnerBasis[residuals, variables,
    MonomialOrder -> Lexicographic,
    CoefficientDomain -> RationalFunctions]]
];

emitLocus[package_String, dimension_Integer, baseQuantity_String,
    equations_List, variables_List, assumptions_] := Module[
  {solution, identityResiduals, identityOperands, identityTest,
   unrestrictedReduction, inconsistentOperands, inconsistentTest,
   realAdmissible, admittedBranches, branch, branchOperands, branchTest,
   branchStatus, point, canonical, fullCondition, realTest, realStatus,
   realWitness, realStatusOperands},
  emitCell[package, dimension, baseQuantity <> "_EQUATIONS", equations];
  solution = Quiet[Solve[And @@ equations, variables]];
  emitCell[package, dimension, baseQuantity <> "_SOLUTION", {
    "SOLVE_VARIABLES" -> variables,
    "SOLUTION_SET" -> solution
  }];

  identityResiduals = relationResidual /@ equations;
  identityOperands = {"RESIDUALS" -> identityResiduals,
    "SOLVE_VARIABLES" -> variables};
  identityTest = unrestrictedSimplify[
    And @@ (# == 0 & /@ identityResiduals)];
  emitCell[package, dimension, baseQuantity <> "_IDENTICALLY_SATISFIED", {
    "STATUS_TOKEN" -> truthStatus[identityTest],
    "TEST_OBJECT" -> identityTest,
    "OPERANDS" -> identityOperands
  }];

  unrestrictedReduction = Quiet[Reduce[And @@ equations, variables, Complexes]];
  inconsistentOperands = {"EQUATIONS" -> equations,
    "SOLVE_VARIABLES" -> variables,
    "UNRESTRICTED_REDUCTION" -> unrestrictedReduction};
  inconsistentTest = SameQ[unrestrictedReduction, False];
  emitCell[package, dimension, baseQuantity <> "_INCONSISTENT", {
    "STATUS_TOKEN" -> truthStatus[inconsistentTest],
    "TEST_OBJECT" -> inconsistentTest,
    "OPERANDS" -> inconsistentOperands
  }];

  admittedBranches = {};
  realAdmissible = If[ListQ[solution],
    Table[
      branch = solution[[branchIndex]];
      branchOperands = {"BRANCH" -> branch, "PREMISES" -> assumptions};
      branchTest = realExistenceTest[And[assumptions, And @@ branchRelations[branch]]];
      branchStatus = admissibilityStatus[branchTest];
      point = If[SameQ[branchStatus, tokenAdmissible],
        exactProjectedInstance[
          And[assumptions, And @@ branchRelations[branch]], variables],
        tokenNotApplicable];
      If[SameQ[branchStatus, tokenAdmissible] && ListQ[point],
        AppendTo[admittedBranches, <|
          "Branch" -> branch,
          "Rules" -> plainBranchRules[branch],
          "Equations" -> branchDefiningEquations[branch],
          "Variables" -> variables,
          "Point" -> point|>]];
      {
        "BRANCH" -> branch,
        "STATUS_TOKEN" -> branchStatus,
        "TEST_OBJECT" -> branchTest,
        "OPERANDS" -> branchOperands
      },
      {branchIndex, Length[solution]}],
    {{
      "BRANCH" -> solution,
      "STATUS_TOKEN" -> tokenUndecided,
      "TEST_OBJECT" -> Failure["OpaqueSolverResult", <|
        "SolverResult" -> solution|>],
      "OPERANDS" -> {"BRANCH" -> solution, "PREMISES" -> assumptions}
    }}
  ];
  emitCell[package, dimension, baseQuantity <> "_REAL_ADMISSIBLE", realAdmissible];

  canonical = canonicalLocus[identityResiduals, variables];
  emitCell[package, dimension, baseQuantity <> "_CANONICAL_LOCUS", canonical];

  fullCondition = And[assumptions, And @@ equations];
  realTest = realExistenceTest[fullCondition];
  realStatus = Which[
    TrueQ[realTest], "PROVED_NONEMPTY",
    SameQ[realTest, False], "PROVED_EMPTY",
    True, tokenUndecided
  ];
  realWitness = If[SameQ[realStatus, "PROVED_NONEMPTY"],
    exactProjectedInstance[fullCondition, variables], tokenNotApplicable];
  If[SameQ[realStatus, "PROVED_NONEMPTY"] && ! ListQ[realWitness],
    realStatus = tokenUndecided;
    realWitness = tokenNotApplicable];
  emitCell[package, dimension, baseQuantity <> "_REAL_STATUS", realStatus];
  emitCell[package, dimension, baseQuantity <> "_REAL_WITNESS", realWitness];
  realStatusOperands = {
    "EQUATIONS" -> equations,
    "SOLVE_VARIABLES" -> variables,
    "PREMISES" -> assumptions,
    "TEST_OBJECT" -> realTest
  };
  emitCell[package, dimension, baseQuantity <> "_REAL_STATUS_OPERANDS",
    realStatusOperands];
  <|
    "BaseQuantity" -> baseQuantity,
    "Equations" -> equations,
    "Variables" -> variables,
    "Solution" -> solution,
    "RealAdmissible" -> realAdmissible,
    "AdmittedBranches" -> admittedBranches,
    "RealStatus" -> realStatus,
    "RealWitness" -> realWitness|>
];

coordinateFrame[dimension_Integer] := Module[
  {spaceCoordinates, timeCoordinate, fieldHeads, arguments, fields, gradient},
  spaceCoordinates = Table[Symbol["x" <> ToString[index]], {index, dimension}];
  timeCoordinate = Symbol["t"];
  fieldHeads = Table[Symbol["u" <> ToString[index]], {index, dimension}];
  arguments = Join[spaceCoordinates, {timeCoordinate}];
  fields = (#[Sequence @@ arguments] &) /@ fieldHeads;
  gradient = Table[D[fields[[column]], spaceCoordinates[[row]]],
    {row, dimension}, {column, dimension}];
  <|
    "X" -> spaceCoordinates,
    "T" -> timeCoordinate,
    "Heads" -> fieldHeads,
    "Arguments" -> arguments,
    "Fields" -> fields,
    "Gradient" -> gradient|>
];

curlDensity[gradient_] := Module[{dimension = Length[gradient]},
  1/2 Total[Flatten[Table[
    (gradient[[row, column]] - gradient[[column, row]])^2,
    {row, dimension}, {column, dimension}]]]
];

divDensity[gradient_] := Tr[gradient]^2;

symtlDensity[gradient_] := Module[{dimension = Length[gradient], tracelessPart},
  tracelessPart = (gradient + Transpose[gradient])/2 -
    Tr[gradient] IdentityMatrix[dimension]/dimension;
  Total[Flatten[tracelessPart^2]]
];

quadraticCoordinates[polynomial_, variables_, pairs_] := Map[
  Function[pair,
    If[pair[[1]] === pair[[2]],
      Coefficient[Expand[polynomial], variables[[pair[[1]]]], 2],
      Coefficient[
        Coefficient[Expand[polynomial], variables[[pair[[1]]]], 1],
        variables[[pair[[2]]]], 1]
    ]],
  pairs
];

polynomialFromCoordinates[row_, monomials_] :=
  Expand[Total[MapThread[Times, {row, monomials}]]];

eulerOperatorForGradientDensity[polynomial_, gradientVariables_, frame_] := Module[
  {dimension, gradient, spaceCoordinates},
  dimension = Length[frame["X"]];
  gradient = frame["Gradient"];
  spaceCoordinates = frame["X"];
  Table[
    Expand[Total[Table[
      D[D[polynomial, gradientVariables[[(row - 1) dimension + column]]] /.
        Thread[gradientVariables -> Flatten[gradient]], spaceCoordinates[[row]]],
      {row, dimension}]]],
    {column, dimension}]
];

buildInvariantCensus[dimension_Integer, package_String] := Module[
  {variableCount, gradientVariables, gradientMatrix, pairs, monomials,
   generators, variationVariables, actionRows, connectedConstraints,
   v1Raw, v1Basis, v1Dimension, reflection, reflectedGradient,
   reflectionRows, fullConstraints, v2Raw, v2Basis, v2Dimension,
   v1Polynomials, reflectedPolynomials, v4Residual, v4Sum, frame,
   v5Euler, pivots, coordinateActionRows, v6Operator, v6Coordinates,
   v6Ambient, v6Basis, v6Dimension, pdPolynomial},
  variableCount = dimension^2;
  gradientVariables = Table[Symbol["g" <> ToString[index]],
    {index, variableCount}];
  gradientMatrix = Partition[gradientVariables, dimension];
  pairs = Flatten[Table[{left, right}, {left, variableCount},
    {right, left, variableCount}], 1];
  monomials = (gradientVariables[[#[[1]]]] gradientVariables[[#[[2]]]] &) /@ pairs;
  emitCell[package, dimension, "MONOMIAL_ORDERING", monomials];

  generators = Flatten[Table[
    SparseArray[{{row, column} -> 1, {column, row} -> -1},
      {dimension, dimension}],
    {row, 1, dimension - 1}, {column, row + 1, dimension}], 1];
  connectedConstraints = If[generators === {},
    SparseArray[{}, {0, Length[monomials]}],
    SparseArray[Join @@ Table[
      variationVariables = Flatten[
        generators[[generatorIndex]].gradientMatrix -
          gradientMatrix.generators[[generatorIndex]]];
      actionRows = Table[
        quadraticCoordinates[
          Expand[Total[MapThread[
            (D[monomials[[monomialIndex]], #1] #2) &,
            {gradientVariables, variationVariables}]]],
          gradientVariables, pairs],
        {monomialIndex, Length[monomials]}];
      Transpose[actionRows],
      {generatorIndex, Length[generators]}]]];
  v1Raw = NullSpace[connectedConstraints];
  v1Basis = If[v1Raw === {}, {}, RowReduce[v1Raw]];
  v1Dimension = Length[v1Basis];
  emitCell[package, dimension, "V1_BASIS", v1Basis];
  emitCell[package, dimension, "V1_DIM", v1Dimension];

  reflection = DiagonalMatrix[Join[{-1}, ConstantArray[1, dimension - 1]]];
  reflectedGradient = Expand[reflection.gradientMatrix.Transpose[reflection]];
  reflectionRows = Table[
    quadraticCoordinates[
      monomials[[monomialIndex]] /.
        Thread[gradientVariables -> Flatten[reflectedGradient]],
      gradientVariables, pairs],
    {monomialIndex, Length[monomials]}];
  fullConstraints = Join[connectedConstraints,
    SparseArray[Transpose[reflectionRows - IdentityMatrix[Length[monomials]]]]];
  v2Raw = NullSpace[fullConstraints];
  v2Basis = If[v2Raw === {}, {}, RowReduce[v2Raw]];
  v2Dimension = Length[v2Basis];
  emitCell[package, dimension, "V2_BASIS", v2Basis];
  emitCell[package, dimension, "V2_DIM", v2Dimension];
  emitCell[package, dimension, "V3_DIFFERENCE", v1Dimension - v2Dimension];

  emitCell[package, dimension, "R0_MATRIX", reflection];
  emitCell[package, dimension, "R0_ORTHOGONALITY_RESIDUAL",
    Transpose[reflection].reflection - IdentityMatrix[dimension]];
  emitCell[package, dimension, "R0_DETERMINANT", Det[reflection]];

  v1Polynomials = polynomialFromCoordinates[#, monomials] & /@ v1Basis;
  reflectedPolynomials = Expand[# /.
      Thread[gradientVariables -> Flatten[reflectedGradient]]] & /@ v1Polynomials;
  v4Residual = MapThread[Expand[#1 - #2] &,
    {reflectedPolynomials, v1Polynomials}];
  v4Sum = MapThread[Expand[#1 + #2] &,
    {reflectedPolynomials, v1Polynomials}];
  emitCell[package, dimension, "V4_REFLECTED", reflectedPolynomials];
  emitCell[package, dimension, "V4_RESIDUAL", v4Residual];
  emitCell[package, dimension, "V4_SUM", v4Sum];

  frame = coordinateFrame[dimension];
  v5Euler = eulerOperatorForGradientDensity[#, gradientVariables, frame] & /@
    v1Polynomials;
  emitCell[package, dimension, "V5_EULER_LAGRANGE", v5Euler];

  If[v1Dimension == 0,
    v6Operator = {};
    v6Basis = {},
    pivots = FirstPosition[#, 1][[1]] & /@ v1Basis;
    coordinateActionRows = Table[
      quadraticCoordinates[reflectedPolynomials[[basisIndex]],
        gradientVariables, pairs][[pivots]],
      {basisIndex, v1Dimension}];
    v6Operator = Transpose[coordinateActionRows];
    v6Coordinates = NullSpace[v6Operator + IdentityMatrix[v1Dimension]];
    v6Ambient = (#.v1Basis) & /@ v6Coordinates;
    v6Basis = If[v6Ambient === {}, {}, RowReduce[v6Ambient]]
  ];
  v6Dimension = Length[v6Basis];
  emitCell[package, dimension, "V6_OPERATOR", v6Operator];
  emitCell[package, dimension, "V6_BASIS", v6Basis];
  emitCell[package, dimension, "V6_DIM", v6Dimension];
  emitCell[package, dimension, "V7_RESIDUAL",
    v6Dimension - (v1Dimension - v2Dimension)];

  pdPolynomial = Expand[Total[
    polynomialFromCoordinates[#, monomials] & /@ v6Basis]];
  <|
    "D" -> dimension,
    "GVariables" -> gradientVariables,
    "GMatrix" -> gradientMatrix,
    "Monomials" -> monomials,
    "V1Basis" -> v1Basis,
    "V1Dim" -> v1Dimension,
    "V2Basis" -> v2Basis,
    "V2Dim" -> v2Dimension,
    "R0" -> reflection,
    "V1Polynomials" -> v1Polynomials,
    "ReflectedPolynomials" -> reflectedPolynomials,
    "V4Residual" -> v4Residual,
    "V4Sum" -> v4Sum,
    "V5Euler" -> v5Euler,
    "V6Operator" -> v6Operator,
    "V6Basis" -> v6Basis,
    "V6Dim" -> v6Dimension,
    "PDPolynomial" -> pdPolynomial|>
];

emitInvariantCensus[data_, package_String, dimension_Integer] := Module[{},
  emitCell[package, dimension, "MONOMIAL_ORDERING", data["Monomials"]];
  emitCell[package, dimension, "V1_BASIS", data["V1Basis"]];
  emitCell[package, dimension, "V1_DIM", Length[data["V1Basis"]]];
  emitCell[package, dimension, "V2_BASIS", data["V2Basis"]];
  emitCell[package, dimension, "V2_DIM", Length[data["V2Basis"]]];
  emitCell[package, dimension, "V3_DIFFERENCE",
    Length[data["V1Basis"]] - Length[data["V2Basis"]]];
  emitCell[package, dimension, "R0_MATRIX", data["R0"]];
  emitCell[package, dimension, "R0_ORTHOGONALITY_RESIDUAL",
    Transpose[data["R0"]].data["R0"] - IdentityMatrix[dimension]];
  emitCell[package, dimension, "R0_DETERMINANT", Det[data["R0"]]];
  emitCell[package, dimension, "V4_REFLECTED", data["ReflectedPolynomials"]];
  emitCell[package, dimension, "V4_RESIDUAL", data["V4Residual"]];
  emitCell[package, dimension, "V4_SUM", data["V4Sum"]];
  emitCell[package, dimension, "V5_EULER_LAGRANGE", data["V5Euler"]];
  emitCell[package, dimension, "V6_OPERATOR", data["V6Operator"]];
  emitCell[package, dimension, "V6_BASIS", data["V6Basis"]];
  emitCell[package, dimension, "V6_DIM", Length[data["V6Basis"]]];
  emitCell[package, dimension, "V7_RESIDUAL",
    Length[data["V6Basis"]] -
      (Length[data["V1Basis"]] - Length[data["V2Basis"]])];
];

stiffnessBlueprint[package_String] := Switch[package,
  "MAIN" | "XKIN_ANISO", {{muR/2, "curl"}, {bComp/2, "div"}},
  "XFORM_CURLONLY", {{muR/2, "curl"}},
  "XFORM_DIVONLY", {{bComp/2, "div"}},
  "XFORM_TRACELESS", {{muR/2, "curl"}, {muBr, "symtl"}},
  "XFORM_EXTRA", {{muR/2, "curl"}, {bComp/2, "div"}, {beta/2, "pd"}},
  "XCOEF_BSCALE", {{muR/2, "curl"}, {s bComp/2, "div"}},
  "XCOEF_BSIGN", {{muR/2, "curl"}, {-bComp/2, "div"}}
];

kineticRecords[package_String, velocityJet_List] := If[package === "XKIN_ANISO",
  {
    <|"Factor" -> rhoBr/2, "DensityJet" -> Total[Rest[velocityJet]^2]|>,
    <|"Factor" -> sRho rhoBr/2, "DensityJet" -> First[velocityJet]^2|>
  },
  {<|"Factor" -> rhoBr/2, "DensityJet" -> Total[velocityJet^2]|>}
];

stiffnessRecords[package_String, densityJet_, densityCoordinate_] := Map[
  Function[item, <|
    "Factor" -> item[[1]],
    "Kind" -> item[[2]],
    "DensityJet" -> densityJet[item[[2]]],
    "DensityCoordinate" -> densityCoordinate[item[[2]]]|>],
  stiffnessBlueprint[package]
];

makePremises[package_String, dimension_Integer, wavevector_List,
    amplitudes_List] := Module[{symbolicEntries, packageEntries, symbolicJoint},
  symbolicEntries = {
    rhoBr > 0,
    muR > 0,
    bComp > 0,
    Total[wavevector^2] > 0,
    And @@ (Element[#, Reals] & /@ wavevector),
    And @@ (Element[#, Reals] & /@ amplitudes),
    Element[dimensionSymbol, Integers] && dimensionSymbol > 0
  };
  packageEntries = Switch[package,
    "XFORM_TRACELESS", {muBr > 0},
    "XFORM_EXTRA", {Element[beta, Reals]},
    "XCOEF_BSCALE", {s > 0, s != 1},
    "XKIN_ANISO", {sRho > 0, sRho != 1},
    _, {}
  ];
  symbolicJoint = And @@ Join[symbolicEntries, packageEntries];
  <|
    "SymbolicEntries" -> Join[symbolicEntries, packageEntries],
    "SymbolicJoint" -> symbolicJoint,
    "AtDimension" -> (symbolicJoint /. dimensionSymbol -> dimension)|>
];

premiseInventory[premises_, package_String, dimension_Integer,
    amplitudes_List] := {
  "JOINT_SYMBOLIC_ASSUMPTIONS" -> premises["SymbolicJoint"],
  "PACKAGE" -> package,
  "BRANE_DIMENSION" -> dimension,
  "DISPLACEMENT_DIMENSION" -> {1, 0, 0},
  "IN_PLANE_COMPONENTS" -> amplitudes,
  "OUT_OF_PLANE_FIELD_EXCLUDED" -> hBranon,
  "BACKGROUND_VELOCITY" -> 0,
  "DISSIPATIVE_TERMS" -> {},
  "MAXIMUM_FIELD_DEGREE" -> 2,
  "WALL_WIDTH_FIELDS" -> {},
  "REAL_PHASE_ANSATZ" -> Cos,
  "PHASE_AVERAGE_DOMAIN" -> {phaseVariable, 0, 2 Pi},
  "DIMENSIONLESS_COEFFICIENTS" ->
    If[package === "XCOEF_BSCALE", {s}, {}],
  "BULK_MODE_CONTENT" -> {ScalarSoundMode},
  "BULK_AMPLITUDE_PREMISE" -> Element[A, Reals],
  "BULK_SPEED_PREMISE" -> cs0 > 0,
  "INTERFACE_EQUATIONS_SUPPLIED" -> {}
};

multiplicativeFactors[expression_] := Module[{factored = Factor[expression]},
  If[Head[factored] === Times, List @@ factored, {factored}]
];

routeAStrippedFactor[planeExpressions_, amplitudes_, phaseCoordinates_List] := Module[
  {coefficientEntries, nonzeroEntries, commonFactors, phaseFactors},
  coefficientEntries = Flatten[Table[
    Coefficient[planeExpressions[[row]], amplitudes[[column]]],
    {row, Length[planeExpressions]}, {column, Length[amplitudes]}]];
  nonzeroEntries = Select[coefficientEntries, ! TrueQ[# === 0] &];
  commonFactors = If[nonzeroEntries === {}, {},
    Fold[Intersection[#1, #2, SameTest -> SameQ] &,
      multiplicativeFactors[First[nonzeroEntries]],
      multiplicativeFactors /@ Rest[nonzeroEntries]]];
  phaseFactors = Select[commonFactors,
    ! FreeQ[#, Alternatives @@ phaseCoordinates] &];
  If[phaseFactors === {}, Missing["NoCommonPhaseFactor"], Times @@ phaseFactors]
];

coefficientOrderingFromRecords[kineticRecordsList_, stiffnessRecordsList_] := Module[
  {factors, inertialCoefficients},
  factors = Join[Lookup[kineticRecordsList, "Factor"],
    Lookup[stiffnessRecordsList, "Factor"]];
  inertialCoefficients = Flatten[globalSymbolsIn /@ Lookup[kineticRecordsList, "Factor"]];
  SortBy[DeleteDuplicates[Join[
    Flatten[globalSymbolsIn /@ factors], inertialCoefficients]], SymbolName]
];

mergeMultiplicityPairs[pairs_List, assumptions_] := Fold[
  Function[{accumulator, pair}, Module[{position, index},
    position = FirstPosition[accumulator,
      existing_List /; Length[existing] == 2 &&
        TrueQ[engineSimplify[existing[[1]] == pair[[1]], assumptions]],
      Missing["NotFound"], {1}];
    If[MissingQ[position], Append[accumulator, pair],
      index = First[position];
      ReplacePart[accumulator, position -> {
        accumulator[[index, 1]],
        accumulator[[index, 2]] + pair[[2]]}]]
  ]],
  {}, pairs
];

rootMultiplicityPairs[determinant_, assumptions_] := Module[
  {polynomial, factorPairs, rawPairs, solutions, roots},
  polynomial = Numerator[Together[determinant]];
  factorPairs = Select[Rest[FactorList[polynomial]],
    ! FreeQ[First[#], omegaSquared] &];
  rawPairs = Flatten[Map[
    Function[factorPair,
      solutions = Quiet[Solve[factorPair[[1]] == 0, omegaSquared]];
      roots = omegaSquared /. solutions;
      ({#, factorPair[[2]]} &) /@ roots],
    factorPairs], 1];
  mergeMultiplicityPairs[rawPairs, assumptions]
];

distinctUnderAssumptions[values_List, assumptions_] := DeleteDuplicates[
  values, TrueQ[engineSimplify[#1 == #2, assumptions]] &];

allMaximalMinors[matrix_, rank_Integer] := Module[{rowSets, columnSets},
  If[rank == 0, Return[{}]];
  rowSets = Subsets[Range[Length[matrix]], {rank}];
  columnSets = Subsets[Range[Length[First[matrix]]], {rank}];
  Flatten[Table[
    Factor[Det[matrix[[rowSet, columnSet]]]],
    {rowSet, rowSets}, {columnSet, columnSets}]]
];

unscopedQuantity[scope_String, pointEvidence_, quantity_String] :=
  scope <> If[TrueQ[pointEvidence], "POINT_EVIDENCE_", ""] <> quantity;

rootQuantity[scope_String, pointEvidence_, rootIndex_Integer, quantity_String] :=
  scope <> "ROOT" <> ToString[rootIndex] <> "_" <>
    If[TrueQ[pointEvidence], "POINT_EVIDENCE_", ""] <> quantity;

changeLocusSuffixes = {
  "EQUATIONS", "SOLUTION", "IDENTICALLY_SATISFIED", "INCONSISTENT",
  "REAL_ADMISSIBLE", "CANONICAL_LOCUS", "REAL_STATUS", "REAL_WITNESS",
  "REAL_STATUS_OPERANDS"};

emitComponentCount[package_, dimension_, baseQuantity_String, value_,
    freeParameters_List] := Module[{status, certificate},
  status = If[freeParameters === {}, tokenConstant, tokenUndecided];
  certificate = If[SameQ[status, tokenConstant], {
      "FREE_PARAMETERS" -> freeParameters,
      "TEST_OBJECT" -> True,
      "OPERANDS" -> {"VALUE" -> value, "FREE_PARAMETERS" -> freeParameters}
    }, tokenNotApplicable];
  emitCell[package, dimension, baseQuantity, {
    "STATUS_TOKEN" -> status,
    "VALUE" -> value
  }];
  emitCell[package, dimension, baseQuantity <> "_STATUS", status];
  emitCell[package, dimension, baseQuantity <> "_CONSTANCY_CERTIFICATE", certificate];
  Scan[Function[suffix,
    emitCell[package, dimension,
      baseQuantity <> "_CHANGE_LOCUS_" <> suffix, tokenNotApplicable]],
    changeLocusSuffixes];
  status
];

emitCountObject[package_, dimension_, baseQuantity_String, value_,
    componentCounts_, freeParameters_List] := If[TrueQ[componentCounts],
  emitComponentCount[package, dimension, baseQuantity, value, freeParameters],
  emitCell[package, dimension, baseQuantity, value];
  tokenNotApplicable
];

candidateRecords[locus_, sourceToken_, sourceLocusTag_] := Map[
  Function[branchRecord, Join[branchRecord, <|
    "Source" -> sourceToken,
    "SourceLocusTag" -> sourceLocusTag|>]],
  locus["AdmittedBranches"]
];

computeSpectrumAndModes[matrix_, coefficients_List, wavevector_List, assumptions_,
    package_String, dimension_Integer, scope_String, pointEvidence_,
    includeRankLoci_, componentCounts_, freeParameters_List,
    activeKVariables_List, activeCoefficientVariables_List] := Module[
  {determinant, multiplicityPairs, solverSolution, degree, rootCountAll,
   roots, rootCountDistinct, rootCountStatuses = {}, pairIndices,
   coincidenceEquations, coincidenceK, coincidenceCoefficient,
   coincidenceCandidates = {}, rankCandidates = {}, stackedCandidates = {},
   rootObjects, root, rootBase, matrixAtRoot, rank, nullity, stacked,
   stackedRank, transverseNullity, basis, basisDots, basisResiduals,
   basisCount, mDotK, rankMinors, stackedMinors, rankEquations, stackedEquations,
   rankKLocus, rankCoefficientLocus, rankJointLocus, stackedKLocus,
   stackedCoefficientLocus, stackedJointLocus, rootCountBase,
   countStatus, q4Derived = {}, rootPrefix},
  determinant = Factor[Det[matrix]];
  emitCell[package, dimension,
    unscopedQuantity[scope, pointEvidence, "DET_M"], determinant];
  multiplicityPairs = rootMultiplicityPairs[determinant, assumptions];
  emitCell[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_MULTIPLICITY_PAIRS"],
    multiplicityPairs];
  solverSolution = Quiet[Solve[determinant == 0, omegaSquared]];
  emitCell[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_SOLUTION_SET"],
    solverSolution];
  rootCountAll = Total[Last /@ multiplicityPairs];
  rootCountBase = unscopedQuantity[scope, pointEvidence, "ROOT_COUNT_ALL"];
  countStatus = emitCountObject[package, dimension, rootCountBase, rootCountAll,
    componentCounts, freeParameters];
  If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
  degree = Exponent[Numerator[Together[determinant]], omegaSquared];
  countStatus = emitCountObject[package, dimension,
    unscopedQuantity[scope, pointEvidence, "DET_M_DEGREE"], degree,
    componentCounts, freeParameters];
  If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
  emitCell[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_DEGREE_RESIDUAL"],
    degree - rootCountAll];
  roots = distinctUnderAssumptions[First /@ multiplicityPairs, assumptions];
  emitCell[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_DISTINCT"], roots];
  rootCountDistinct = Length[roots];
  countStatus = emitCountObject[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_COUNT_DISTINCT"],
    rootCountDistinct, componentCounts, freeParameters];
  If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
  emitCell[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_ORDERING"], roots];

  Do[
    rootBase = rootQuantity[scope, pointEvidence, rootIndex, "VALUE"];
    emitCell[package, dimension, rootBase, roots[[rootIndex]]];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "SIGN"],
      signPayload[roots[[rootIndex]], assumptions]],
    {rootIndex, Length[roots]}];

  pairIndices = Subsets[Range[Length[roots]], {2}];
  coincidenceEquations = (engineSimplify[
      roots[[#[[1]]]] - roots[[#[[2]]]], assumptions] == 0) & /@ pairIndices;
  coincidenceK = emitLocus[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_COINCIDENCE_K"],
    coincidenceEquations, activeKVariables, assumptions];
  coincidenceCoefficient = emitLocus[package, dimension,
    unscopedQuantity[scope, pointEvidence, "ROOT_COINCIDENCE_COEFF"],
    coincidenceEquations, activeCoefficientVariables, assumptions];
  If[TrueQ[includeRankLoci],
    coincidenceCandidates = Join[
      candidateRecords[coincidenceK, tokenRootCoincidence,
        standardEmissionName[cellEmission[package, dimension,
          unscopedQuantity[scope, pointEvidence, "ROOT_COINCIDENCE_K"]]]],
      candidateRecords[coincidenceCoefficient, tokenRootCoincidence,
        standardEmissionName[cellEmission[package, dimension,
          unscopedQuantity[scope, pointEvidence, "ROOT_COINCIDENCE_COEFF"]]]]
    ]];

  rootObjects = Table[
    root = roots[[rootIndex]];
    matrixAtRoot = Map[engineSimplify[#, assumptions] &,
      matrix /. omegaSquared -> root, {2}];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N1"], matrixAtRoot];
    rank = assumedRank[matrixAtRoot, assumptions];
    countStatus = emitCountObject[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N2_RANK"], rank,
      componentCounts, freeParameters];
    If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
    nullity = dimension - rank;
    countStatus = emitCountObject[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N2_NULLITY"], nullity,
      componentCounts, freeParameters];
    If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
    stacked = Join[matrixAtRoot, {wavevector}];
    stackedRank = assumedRank[stacked, assumptions];
    countStatus = emitCountObject[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N3_STACKED_RANK"],
      stackedRank, componentCounts, freeParameters];
    If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
    transverseNullity = dimension - stackedRank;
    countStatus = emitCountObject[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N3_TRANSVERSE_NULLITY"],
      transverseNullity, componentCounts, freeParameters];
    If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N4_NULLITY_DIFFERENCE"],
      nullity - transverseNullity];
    mDotK = Map[engineSimplify[#, assumptions] &, matrixAtRoot.wavevector];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N5_M_DOT_K"], mDotK];
    basis = assumedNullSpace[matrixAtRoot, assumptions];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N6_BASIS"], basis];
    basisDots = Map[engineSimplify[#.wavevector, assumptions] &, basis];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N6_DOT_K"], basisDots];
    basisResiduals = MapThread[
      Function[{basisVector, dotProduct},
        Map[engineSimplify[#, assumptions] &,
          Total[wavevector^2] basisVector - dotProduct wavevector]],
      {basis, basisDots}];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N6_RESIDUAL"],
      basisResiduals];
    basisCount = Length[basis];
    countStatus = emitCountObject[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N7_BASIS_COUNT"],
      basisCount, componentCounts, freeParameters];
    If[TrueQ[componentCounts], AppendTo[rootCountStatuses, countStatus]];
    emitCell[package, dimension,
      rootQuantity[scope, pointEvidence, rootIndex, "N7_RESIDUAL"],
      basisCount - nullity];

    If[TrueQ[includeRankLoci],
      rankMinors = allMaximalMinors[matrixAtRoot, rank];
      emitCell[package, dimension,
        rootQuantity[scope, pointEvidence, rootIndex, "RANK_DROP_MINORS"],
        rankMinors];
      rankEquations = (# == 0) & /@ rankMinors;
      rootPrefix = rootQuantity[scope, pointEvidence, rootIndex, "RANK_DROP_K"];
      rankKLocus = emitLocus[package, dimension, rootPrefix, rankEquations,
        activeKVariables, assumptions];
      rankCandidates = Join[rankCandidates,
        candidateRecords[rankKLocus, tokenRankDrop,
          standardEmissionName[cellEmission[package, dimension, rootPrefix]]]];
      rootPrefix = rootQuantity[scope, pointEvidence, rootIndex, "RANK_DROP_COEFF"];
      rankCoefficientLocus = emitLocus[package, dimension, rootPrefix,
        rankEquations, activeCoefficientVariables, assumptions];
      rankCandidates = Join[rankCandidates,
        candidateRecords[rankCoefficientLocus, tokenRankDrop,
          standardEmissionName[cellEmission[package, dimension, rootPrefix]]]];
      rootPrefix = rootQuantity[scope, pointEvidence, rootIndex, "RANK_DROP_JOINT"];
      rankJointLocus = emitLocus[package, dimension, rootPrefix, rankEquations,
        Join[activeKVariables, activeCoefficientVariables], assumptions];
      rankCandidates = Join[rankCandidates,
        candidateRecords[rankJointLocus, tokenRankDrop,
          standardEmissionName[cellEmission[package, dimension, rootPrefix]]]];

      stackedMinors = allMaximalMinors[stacked, stackedRank];
      emitCell[package, dimension,
        rootQuantity[scope, pointEvidence, rootIndex, "STACKED_DROP_MINORS"],
        stackedMinors];
      stackedEquations = (# == 0) & /@ stackedMinors;
      rootPrefix = rootQuantity[scope, pointEvidence, rootIndex, "STACKED_DROP_K"];
      stackedKLocus = emitLocus[package, dimension, rootPrefix,
        stackedEquations, activeKVariables, assumptions];
      stackedCandidates = Join[stackedCandidates,
        candidateRecords[stackedKLocus, tokenStackedDrop,
          standardEmissionName[cellEmission[package, dimension, rootPrefix]]]];
      rootPrefix = rootQuantity[scope, pointEvidence, rootIndex, "STACKED_DROP_COEFF"];
      stackedCoefficientLocus = emitLocus[package, dimension, rootPrefix,
        stackedEquations, activeCoefficientVariables, assumptions];
      stackedCandidates = Join[stackedCandidates,
        candidateRecords[stackedCoefficientLocus, tokenStackedDrop,
          standardEmissionName[cellEmission[package, dimension, rootPrefix]]]];
      rootPrefix = rootQuantity[scope, pointEvidence, rootIndex, "STACKED_DROP_JOINT"];
      stackedJointLocus = emitLocus[package, dimension, rootPrefix,
        stackedEquations, Join[activeKVariables, activeCoefficientVariables],
        assumptions];
      stackedCandidates = Join[stackedCandidates,
        candidateRecords[stackedJointLocus, tokenStackedDrop,
          standardEmissionName[cellEmission[package, dimension, rootPrefix]]]]
    ];
    AppendTo[q4Derived, <|
      "MDotK" -> mDotK,
      "BasisResiduals" -> basisResiduals|>];
    <|
      "Root" -> root,
      "Matrix" -> matrixAtRoot,
      "Rank" -> rank,
      "Nullity" -> nullity,
      "Stacked" -> stacked,
      "StackedRank" -> stackedRank,
      "Basis" -> basis,
      "BasisResiduals" -> basisResiduals|>,
    {rootIndex, Length[roots]}];
  <|
    "Determinant" -> determinant,
    "MultiplicityPairs" -> multiplicityPairs,
    "Roots" -> roots,
    "RootObjects" -> rootObjects,
    "Q4Derived" -> q4Derived,
    "Candidates" -> Join[rankCandidates, stackedCandidates,
      coincidenceCandidates],
    "CountStatuses" -> rootCountStatuses|>
];

dimensionVectorsEqualQ[left_, right_, assumptions_] :=
  ListQ[left] && ListQ[right] && Length[left] == Length[right] &&
    TrueQ[engineSimplify[And @@ Thread[left == right], assumptions]];

dimensionOfScalar[expression_, atomDimensions_, assumptions_] := Module[
  {pieces, dimensions, uniqueDimensions, baseDimension},
  Which[
    NumberQ[expression], {0, 0, 0},
    KeyExistsQ[atomDimensions, Unevaluated[expression]],
      atomDimensions[Unevaluated[expression]],
    Head[expression] === ConditionalExpression,
      dimensionOfScalar[expression[[1]], atomDimensions, assumptions],
    Head[expression] === Plus,
      pieces = List @@ expression;
      dimensions = dimensionOfScalar[#, atomDimensions, assumptions] & /@ pieces;
      uniqueDimensions = DeleteDuplicates[dimensions,
        dimensionVectorsEqualQ[#1, #2, assumptions] &];
      If[Length[uniqueDimensions] == 1, First[uniqueDimensions],
        DimensionalAlternatives @@ uniqueDimensions],
    Head[expression] === Times,
      dimensions = dimensionOfScalar[#, atomDimensions, assumptions] & /@
        (List @@ expression);
      If[And @@ (ListQ /@ dimensions), Total[dimensions],
        DimensionalProduct[dimensions]],
    Head[expression] === Power &&
        (IntegerQ[expression[[2]]] || RationalQ[expression[[2]]]),
      baseDimension = dimensionOfScalar[expression[[1]], atomDimensions, assumptions];
      If[ListQ[baseDimension], expression[[2]] baseDimension,
        DimensionalPower[baseDimension, expression[[2]]]],
    True, UnknownDimension[expression]
  ]
];

dimensionOfDeclaredTerm[record_, atomDimensions_, assumptions_] := Module[
  {density, factorDimension, densityDimension},
  density = record["DensityJet"];
  If[TrueQ[density === 0], Return[ZeroFormDimensionUndetermined[density]]];
  factorDimension = dimensionOfScalar[record["Factor"], atomDimensions, assumptions];
  densityDimension = dimensionOfScalar[density, atomDimensions, assumptions];
  If[ListQ[factorDimension] && ListQ[densityDimension],
    factorDimension + densityDimension,
    DimensionalProduct[{factorDimension, densityDimension}]]
];

expandedTerms[expression_] := If[TrueQ[expression === 0], {},
  If[Head[Expand[expression]] === Plus, List @@ Expand[expression],
    {Expand[expression]}]];

homogeneityRecord[objectName_String, expression_, atomDimensions_, assumptions_] := Module[
  {scalars, termDimensions, referenceDimensions, testObject},
  scalars = If[ListQ[expression], Flatten[expression], {expression}];
  termDimensions = Map[
    Function[scalar,
      dimensionOfScalar[#, atomDimensions, assumptions] & /@ expandedTerms[scalar]],
    scalars];
  referenceDimensions = DeleteDuplicates[Flatten[termDimensions, 1],
    dimensionVectorsEqualQ[#1, #2, assumptions] &];
  testObject = If[referenceDimensions === {}, True,
    And[And @@ (ListQ /@ referenceDimensions), Length[referenceDimensions] == 1]];
  {
    "OBJECT_NAME" -> objectName,
    "EXPRESSION" -> expression,
    "TERM_DIMENSIONS" -> termDimensions,
    "TEST_OBJECT" -> testObject
  }
];

buildDimensions[coefficientOrdering_List, actionRecords_List, velocityJet_List,
    gradientVariables_List, amplitudes_List, wavevector_List, assumptions_,
    package_String, dimension_Integer] := Module[
  {dimensionlessCoefficients, unknownCoefficients, dimensionVariables,
   coefficientDimensions, baseAtomDimensions, targetDimension, termDimensions,
   equationVectors, dimensionEquations, flatUnknowns, solution, solutionRules,
   finalCoefficientDimensions, firstSlotVariables, firstSlotExpressions,
   coefficientMatrix, equationCount, unknownCount, determinacy,
   finalAtomDimensions, actionHomogeneity, soundSpeedDimensionVariables,
   soundSpeedDimensionEquations, soundSpeedDimensionSolution,
   soundSpeedDimension},
  dimensionlessCoefficients = If[MemberQ[coefficientOrdering, s], {s}, {}];
  unknownCoefficients = Select[coefficientOrdering,
    ! MemberQ[dimensionlessCoefficients, #] &];
  dimensionVariables = Association@Table[
    coefficient -> Table[
      Symbol["dim" <> SymbolName[coefficient] <> "Slot" <> ToString[slot]],
      {slot, 3}],
    {coefficient, unknownCoefficients}];
  coefficientDimensions = Association@Table[
    coefficient -> If[MemberQ[dimensionlessCoefficients, coefficient],
      {0, 0, 0}, dimensionVariables[coefficient]],
    {coefficient, coefficientOrdering}];
  baseAtomDimensions = Join[coefficientDimensions,
    AssociationThread[velocityJet,
      ConstantArray[{1, -1, 0}, Length[velocityJet]]],
    AssociationThread[gradientVariables,
      ConstantArray[{0, 0, 0}, Length[gradientVariables]]],
    AssociationThread[amplitudes,
      ConstantArray[{1, 0, 0}, Length[amplitudes]]],
    AssociationThread[wavevector,
      ConstantArray[{-1, 0, 0}, Length[wavevector]]],
    <|omegaSquared -> {0, -2, 0}, lambdaScale -> {0, 0, 0},
      kwSquared -> {-2, 0, 0}|>];
  targetDimension = {2 - dimensionSymbol, -2, 1};
  termDimensions = dimensionOfDeclaredTerm[#, baseAtomDimensions, assumptions] & /@
    actionRecords;
  equationVectors = Select[termDimensions, ListQ];
  dimensionEquations = Flatten[Thread[# == targetDimension] & /@ equationVectors];
  flatUnknowns = Flatten[Values[dimensionVariables]];
  solution = Quiet[Solve[dimensionEquations, flatUnknowns]];
  emitCell[package, dimension, "DIM_EQUATIONS", dimensionEquations];
  emitCell[package, dimension, "DIM_SOLUTION", solution];
  solutionRules = If[ListQ[solution] && solution =!= {}, First[solution], {}];
  finalCoefficientDimensions = Table[
    {coefficient,
      engineSimplify[coefficientDimensions[coefficient] /. solutionRules,
        assumptions]},
    {coefficient, coefficientOrdering}];
  emitCell[package, dimension, "DIM_COEFFICIENTS", finalCoefficientDimensions];

  firstSlotVariables = If[unknownCoefficients === {}, {},
    dimensionVariables[#][[1]] & /@ unknownCoefficients];
  firstSlotExpressions = If[equationVectors === {}, {},
    (#[[1]] - targetDimension[[1]]) & /@ equationVectors];
  coefficientMatrix = If[firstSlotVariables === {} || firstSlotExpressions === {}, {},
    Normal[Last[CoefficientArrays[firstSlotExpressions, firstSlotVariables]]]];
  equationCount = If[firstSlotVariables === {} || firstSlotExpressions === {}, 0,
    assumedRank[coefficientMatrix, assumptions]];
  unknownCount = Length[unknownCoefficients];
  determinacy = Which[
    ! ListQ[solution] || solution === {}, tokenOverDetermined,
    equationCount > unknownCount, tokenOverDetermined,
    equationCount < unknownCount, tokenUnderDetermined,
    True, tokenExactlyDetermined
  ];
  emitCell[package, dimension, "DIM_EQUATION_COUNT", equationCount];
  emitCell[package, dimension, "DIM_UNKNOWN_COUNT", unknownCount];
  emitCell[package, dimension, "DIM_COUNT_DIFFERENCE", equationCount - unknownCount];
  emitCell[package, dimension, "DIM_DETERMINACY", determinacy];

  soundSpeedDimensionVariables = Table[
    Symbol["dimCs0Slot" <> ToString[slot]], {slot, 3}];
  soundSpeedDimensionEquations = Thread[
    2 soundSpeedDimensionVariables + {-2, 0, 0} == {0, -2, 0}];
  soundSpeedDimensionSolution = Quiet[Solve[
    soundSpeedDimensionEquations, soundSpeedDimensionVariables]];
  soundSpeedDimension = soundSpeedDimensionVariables /.
    First[soundSpeedDimensionSolution];
  finalAtomDimensions = Join[baseAtomDimensions /. solutionRules,
    <|cs0 -> soundSpeedDimension|>];
  actionHomogeneity = Table[
    homogeneityRecord["ACTION_TERM_" <> ToString[index],
      actionRecords[[index]]["Factor"] actionRecords[[index]]["DensityJet"],
      finalAtomDimensions, assumptions],
    {index, Length[actionRecords]}];
  emitCell[package, dimension, "DIM_HOMOGENEITY_ACTION", actionHomogeneity];
  <|
    "CoefficientDimensions" -> finalCoefficientDimensions,
    "AtomDimensions" -> finalAtomDimensions,
    "Solution" -> solution,
    "SolutionRules" -> solutionRules|>
];

emitScaledRoots[spectrum_, wavevector_List, assumptions_, package_String,
    dimension_Integer] := Module[
  {root, scaled, ratio, numeratorExponent, denominatorExponent, candidate,
   exponent, scaleAssumptions},
  scaleAssumptions = And[assumptions, lambdaScale > 0];
  Do[
    root = spectrum["Roots"][[rootIndex]];
    scaled = engineSimplify[root /. Thread[wavevector -> lambdaScale wavevector],
      scaleAssumptions];
    emitCell[package, dimension, "ROOT" <> ToString[rootIndex] <> "_SCALED", scaled];
    emitCell[package, dimension, "ROOT" <> ToString[rootIndex] <> "_UNSCALED", root];
    If[TrueQ[engineSimplify[root == 0, assumptions]],
      ratio = tokenUndefinedRatio;
      exponent = tokenUndefinedRatio,
      ratio = engineSimplify[scaled/root, scaleAssumptions];
      numeratorExponent = Exponent[Numerator[Together[ratio]], lambdaScale];
      denominatorExponent = Exponent[Denominator[Together[ratio]], lambdaScale];
      candidate = numeratorExponent - denominatorExponent;
      exponent = If[TrueQ[engineSimplify[ratio == lambdaScale^candidate,
          scaleAssumptions]], candidate, tokenNotPurePower]
    ];
    emitCell[package, dimension,
      "ROOT" <> ToString[rootIndex] <> "_SCALE_RATIO", ratio];
    emitCell[package, dimension,
      "ROOT" <> ToString[rootIndex] <> "_SCALE_EXPONENT", exponent],
    {rootIndex, Length[spectrum["Roots"]]}]
];

emitQ7[stiffnessRecordList_, coordinateGradient_, independentGradient_,
    package_String, dimension_Integer] := Module[
  {gradientRules, fullStiffness, curlTerm, unweightedCurl, leviCivita,
   curlVector, curlReference, residual},
  gradientRules = Thread[Flatten[coordinateGradient] -> Flatten[independentGradient]];
  fullStiffness = Total[(#["Factor"] #["DensityCoordinate"]) & /@
      stiffnessRecordList] /. gradientRules;
  curlTerm = Total[Map[
      Function[record, If[record["Kind"] === "curl",
        record["Factor"] record["DensityCoordinate"], 0]],
      stiffnessRecordList]] /. gradientRules;
  unweightedCurl = curlDensity[coordinateGradient] /. gradientRules;
  leviCivita = LeviCivitaTensor[dimension];
  curlVector = Table[Total[Flatten[Table[
      leviCivita[[component, row, column]] independentGradient[[row, column]],
      {row, dimension}, {column, dimension}]]],
    {component, dimension}];
  curlReference = Expand[curlVector.curlVector];
  residual = Expand[unweightedCurl - curlReference];
  emitCell[package, dimension, "Q7_W_FULL", Expand[fullStiffness]];
  emitCell[package, dimension, "Q7_CURL_TERM", Expand[curlTerm]];
  emitCell[package, dimension, "Q7_CURL_DENSITY", Expand[unweightedCurl]];
  emitCell[package, dimension, "Q7_CURL_REFERENCE", curlReference];
  emitCell[package, dimension, "Q7_RESIDUAL", residual];
  residual
];

sourceOrder = {tokenRankDrop, tokenStackedDrop, tokenRootCoincidence};
sourcePosition[source_] := FirstPosition[sourceOrder, source][[1]];

mergeCandidates[candidates_List] := Module[{merged = {}, key, position, record},
  Do[
    record = candidates[[candidateIndex]];
    key = ToString[record["Branch"],
      InputForm, PageWidth -> Infinity];
    position = FirstPosition[merged,
      existing_Association /; SameQ[existing["Key"], key],
      Missing["NotFound"], {1}];
    If[MissingQ[position],
      AppendTo[merged, <|
        "Key" -> key,
        "Branch" -> record["Branch"],
        "Rules" -> record["Rules"],
        "Equations" -> record["Equations"],
        "Variables" -> record["Variables"],
        "Point" -> record["Point"],
        "Contributors" -> {{record["Source"], record["SourceLocusTag"]}}|>],
      merged[[First[position], "Contributors"]] = DeleteDuplicates[
        Append[merged[[First[position], "Contributors"]],
          {record["Source"], record["SourceLocusTag"]}]];
      merged[[First[position], "Variables"]] = DeleteDuplicates[Join[
        merged[[First[position], "Variables"]], record["Variables"]]]
    ],
    {candidateIndex, Length[candidates]}];
  merged
];

componentFreeParameters[restrictedMatrix_, wavevector_List,
    coefficients_List] := Select[DeleteDuplicates[Join[wavevector, coefficients]],
  ! FreeQ[restrictedMatrix, #] &];

failureJacobianPayload[componentRoots_, coefficientOrdering_, equations_, point_] := {
  "FAILURE_TOKEN" -> "MISSING_TANGENT_COORDINATES_AND_OFF_STRATUM_EXTENSION",
  "RECOMPUTED_ROOTS" -> componentRoots,
  "COEFFICIENT_ORDERING" -> coefficientOrdering,
  "DEFINING_EQUATIONS" -> equations,
  "EVALUATION_POINT" -> point
};

emitStrata[candidates_List, matrix_, genericRootJacobian_, coefficients_List,
    wavevector_List, assumptions_, package_String, dimension_Integer] := Module[
  {strata, ordering, stratum, contributors, sources, sourceTags, scope,
   restrictionRules, restrictedMatrix, restrictedAssumptions, freeParameters,
   activeK, activeCoefficients, componentSpectrum, coverage, point,
   pointMatrix, pointAssumptions, pointActiveK, pointActiveCoefficients,
   pointSpectrum, restrictedJacobian},
  strata = Map[
    Function[component,
      Append[component, "Point" -> exactProjectedInstance[
        And[assumptions, And @@ component["Equations"]],
        component["Variables"]]]],
    mergeCandidates[candidates]];
  strata = Select[strata, ListQ[#["Point"]] &];
  ordering = Lookup[strata, "Branch", {}];
  emitCell[package, dimension, "STRATUM_ORDERING", ordering];
  Do[
    stratum = strata[[stratumIndex]];
    scope = "STRATUM" <> ToString[stratumIndex] <> "_";
    contributors = SortBy[stratum["Contributors"],
      sourcePosition[First[#]] &];
    sources = DeleteDuplicates[First /@ contributors];
    sourceTags = Last /@ contributors;
    restrictionRules = stratum["Rules"];
    restrictedMatrix = Map[engineSimplify[#,
          assumptions /. restrictionRules] &,
      matrix /. restrictionRules, {2}];
    restrictedAssumptions = engineSimplify[assumptions /. restrictionRules, True];
    freeParameters = componentFreeParameters[restrictedMatrix, wavevector, coefficients];
    point = stratum["Point"];

    emitCell[package, dimension, scope <> "SOURCES", sources];
    emitCell[package, dimension, scope <> "SOURCE_LOCUS_TAGS", sourceTags];
    emitCell[package, dimension, scope <> "DEFINING_EQUATIONS",
      stratum["Equations"]];
    emitCell[package, dimension, scope <> "FREE_PARAMETERS", freeParameters];
    emitCell[package, dimension, scope <> "POINT", point];

    activeK = Select[wavevector, MemberQ[freeParameters, #] &];
    activeCoefficients = Select[coefficients, MemberQ[freeParameters, #] &];
    componentSpectrum = computeSpectrumAndModes[restrictedMatrix, coefficients,
      wavevector, restrictedAssumptions, package, dimension, scope, False,
      False, True, freeParameters, activeK, activeCoefficients];
    coverage = If[And @@ (MemberQ[{tokenConstant, tokenVaries}, #] & /@
          componentSpectrum["CountStatuses"]),
      tokenCompleteCoverage, tokenIncompleteCoverage];
    emitCell[package, dimension, scope <> "COMPONENT_Q3_Q4_COVERAGE", coverage];
    ClearSystemCache[];

    pointMatrix = Map[engineSimplify[#, assumptions /. point] &,
      matrix /. point, {2}];
    pointAssumptions = engineSimplify[assumptions /. point, True];
    pointActiveK = Select[wavevector, ! FreeQ[pointMatrix, #] &];
    pointActiveCoefficients = Select[coefficients, ! FreeQ[pointMatrix, #] &];
    pointSpectrum = computeSpectrumAndModes[pointMatrix, coefficients, wavevector,
      pointAssumptions, package, dimension, scope, True, False, False, {},
      pointActiveK, pointActiveCoefficients];

    restrictedJacobian = Quiet[genericRootJacobian /. point];
    emitCell[package, dimension,
      scope <> "POINT_EVIDENCE_ROOT_COEFFICIENT_JACOBIAN_RESTRICTED",
      restrictedJacobian];
    emitCell[package, dimension,
      scope <> "ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED",
      failureJacobianPayload[componentSpectrum["Roots"], coefficients,
        stratum["Equations"], point]];
    ClearSystemCache[],
    {stratumIndex, Length[strata]}];
  strata
];

closureRankFromMatrix[matrix_List, assumptions_] := Module[
  {rowCount, columnCount, possibleRanks, realisedRanks},
  rowCount = Length[matrix];
  columnCount = If[rowCount == 0, 0, Length[First[matrix]]];
  possibleRanks = Range[Min[rowCount, columnCount]];
  realisedRanks = Select[possibleRanks,
    Function[rank,
      AnyTrue[allMaximalMinors[matrix, rank],
        ! TrueQ[engineSimplify[# == 0, assumptions]] &]]];
  Max[Prepend[realisedRanks, 0]]
];

emitQ11[spectrum_, coefficients_List, wavevector_List, amplitudes_List,
    assumptions_, package_String, dimension_Integer] := Module[
  {bulkAssumptions, bulkFieldRecord, bulkDispersion, kwObjects, root,
   kwEquation, kwSolutions, kwExpression, closureSourceEquations,
   closureEquations, closureUnknowns, closureMatrix, closureRank},
  bulkAssumptions = And[assumptions, Element[A, Reals], cs0 > 0];
  bulkFieldRecord = <|
    "Amplitude" -> A,
    "FieldEquation" -> (bulkPhi[Sequence @@ wavevector, w, t] ==
      A Cos[Total[MapThread[Times, {wavevector,
          Table[Symbol["x" <> ToString[index]], {index, dimension}]}]] +
        kw w - omegaLinear t]),
    "Dispersion" -> (omegaSquared == cs0^2 (Total[wavevector^2] + kwSquared))|>;
  bulkDispersion = bulkFieldRecord["Dispersion"];
  kwObjects = Table[
    root = spectrum["Roots"][[rootIndex]];
    kwEquation = bulkDispersion /. omegaSquared -> root;
    emitCell[package, dimension,
      "ROOT" <> ToString[rootIndex] <> "_KW_EQUATION", kwEquation];
    kwSolutions = Quiet[Solve[kwEquation, kwSquared]];
    kwExpression = If[ListQ[kwSolutions] && kwSolutions =!= {},
      kwSquared /. First[kwSolutions],
      Failure["UnsolvedBulkNormalWavevector", <|
        "Equation" -> kwEquation, "SolverResult" -> kwSolutions|>]];
    emitCell[package, dimension,
      "ROOT" <> ToString[rootIndex] <> "_KW_SQUARED", kwExpression];
    emitCell[package, dimension,
      "ROOT" <> ToString[rootIndex] <> "_KW_SIGN",
      signPayload[kwExpression, bulkAssumptions]];
    emitLocus[package, dimension,
      "ROOT" <> ToString[rootIndex] <> "_KW_ZERO_LOCUS",
      {kwExpression == 0}, Join[coefficients, {cs0}], bulkAssumptions];
    kwExpression,
    {rootIndex, Length[spectrum["Roots"]]}];

  closureSourceEquations = Join[
    Table[braneField[component] == amplitudes[[component]] Cos[phaseVariable],
      {component, dimension}],
    {bulkFieldRecord["FieldEquation"], bulkFieldRecord["Dispersion"]}];
  closureEquations = Select[closureSourceEquations,
    Function[equation,
      ! FreeQ[equation, bulkFieldRecord["Amplitude"]] &&
        AnyTrue[amplitudes, ! FreeQ[equation, #] &]]];
  closureUnknowns = Join[amplitudes, {bulkFieldRecord["Amplitude"]}];
  closureMatrix = Table[
    Coefficient[relationResidual[equation], unknown],
    {equation, closureEquations}, {unknown, closureUnknowns}];
  closureRank = closureRankFromMatrix[closureMatrix, bulkAssumptions];
  emitCell[package, dimension, "C1_EQUATIONS", closureEquations];
  emitCell[package, dimension, "C2_UNKNOWNS", closureUnknowns];
  emitCell[package, dimension, "C2_COUNT", Length[closureUnknowns]];
  emitCell[package, dimension, "C3_RANK", closureRank];
  emitCell[package, dimension, "C4_DIFFERENCE",
    Length[closureUnknowns] - closureRank];
  kwObjects
];

runCell[package_String, dimension_Integer, invariantData_] := Module[
  {frame, spaceCoordinates, timeCoordinate, fields, coordinateGradient,
   gradientVariables, gradientJet, velocityJet, bareFieldJet, amplitudes,
   wavevector, phase, premises, assumptions, bulkAssumptions,
   pdJet, pdCoordinate, densityJet, densityCoordinate, kineticRecordList,
   stiffnessRecordList, jetRules, kineticTermsJet, kineticTermsCoordinate,
   stiffnessTermsJet, stiffnessTermsCoordinate, stiffnessFunctional,
   lagrangianJet, lagrangianCoordinate, eulerExpressions, eulerSystem,
   ansatzRules, planeEuler, strippedFactor, amplitudeEquations, matrixA,
   planeGradientRules, planeVelocityRules, planeLagrangian,
   averagedLagrangian, matrixB, matrixResidual, routeSelection,
   downstreamMatrix, coefficientOrdering, matrixCoefficientJacobian,
   actionRecords, dimensionData, spectrum, rootDimensionObjects,
   scaledObjects, genericRootJacobian, strata, q7Residual, kwObjects,
   derivedHomogeneity, independentGradient},
  frame = coordinateFrame[dimension];
  spaceCoordinates = frame["X"];
  timeCoordinate = frame["T"];
  fields = frame["Fields"];
  coordinateGradient = frame["Gradient"];
  gradientVariables = invariantData["GVariables"];
  gradientJet = Partition[gradientVariables, dimension];
  velocityJet = Table[Symbol["v" <> ToString[index]], {index, dimension}];
  bareFieldJet = Table[Symbol["uBare" <> ToString[index]], {index, dimension}];
  amplitudes = Table[Symbol["a" <> ToString[index]], {index, dimension}];
  wavevector = Table[Symbol["k" <> ToString[index]], {index, dimension}];
  phase = Total[MapThread[Times, {wavevector, spaceCoordinates}]] -
    omegaLinear timeCoordinate;

  premises = makePremises[package, dimension, wavevector, amplitudes];
  assumptions = premises["AtDimension"];
  bulkAssumptions = And[assumptions, Element[A, Reals], cs0 > 0];
  emitLocalCell[package, dimension, "PREMISE_INVENTORY",
    premiseInventory[premises, package, dimension, amplitudes]];

  pdJet = invariantData["PDPolynomial"];
  pdCoordinate = pdJet /. Thread[gradientVariables -> Flatten[coordinateGradient]];
  emitCell[package, dimension, "PD_TERM", Expand[pdCoordinate]];
  densityJet = <|
    "curl" -> curlDensity[gradientJet],
    "div" -> divDensity[gradientJet],
    "symtl" -> symtlDensity[gradientJet],
    "pd" -> pdJet|>;
  densityCoordinate = <|
    "curl" -> curlDensity[coordinateGradient],
    "div" -> divDensity[coordinateGradient],
    "symtl" -> symtlDensity[coordinateGradient],
    "pd" -> pdCoordinate|>;
  kineticRecordList = kineticRecords[package, velocityJet];
  stiffnessRecordList = stiffnessRecords[package, densityJet, densityCoordinate];
  jetRules = Join[
    Thread[velocityJet -> (D[#, timeCoordinate] & /@ fields)],
    Thread[gradientVariables -> Flatten[coordinateGradient]],
    Thread[bareFieldJet -> fields]];
  kineticTermsJet = (#["Factor"] #["DensityJet"]) & /@ kineticRecordList;
  kineticTermsCoordinate = kineticTermsJet /. jetRules;
  stiffnessTermsJet = (#["Factor"] #["DensityJet"]) & /@ stiffnessRecordList;
  stiffnessTermsCoordinate = (#["Factor"] #["DensityCoordinate"]) & /@
    stiffnessRecordList;
  stiffnessFunctional = Total[stiffnessTermsCoordinate];
  lagrangianJet = Total[kineticTermsJet] - Total[stiffnessTermsJet];
  lagrangianCoordinate = Expand[lagrangianJet /. jetRules];
  emitCell[package, dimension, "LAGRANGIAN", lagrangianCoordinate];
  emitCell[package, dimension, "KINETIC_TERMS", kineticTermsCoordinate];
  emitCell[package, dimension, "STIFFNESS_TERMS", stiffnessTermsCoordinate];

  eulerExpressions = Table[Expand[
      D[D[lagrangianJet, velocityJet[[component]]] /. jetRules,
        timeCoordinate] +
      Total[Table[
        D[D[lagrangianJet, gradientJet[[row, component]]] /. jetRules,
          spaceCoordinates[[row]]],
        {row, dimension}]] -
      (D[lagrangianJet, bareFieldJet[[component]]] /. jetRules)],
    {component, dimension}];
  eulerSystem = (# == 0) & /@ eulerExpressions;
  emitCell[package, dimension, "EULER_LAGRANGE_SYSTEM", eulerSystem];

  ansatzRules = Table[With[
      {head = frame["Heads"][[component]],
       amplitude = amplitudes[[component]],
       arguments = frame["Arguments"], phaseExpression = phase},
      head -> Function[Evaluate[arguments],
        Evaluate[amplitude Cos[phaseExpression]]]],
    {component, dimension}];
  planeEuler = TrigFactor[eulerExpressions /. ansatzRules] /.
    omegaLinear^2 -> omegaSquared;
  strippedFactor = routeAStrippedFactor[planeEuler, amplitudes,
    Join[spaceCoordinates, {timeCoordinate}]];
  If[MissingQ[strippedFactor],
    Throw[Failure["RouteAStrippedFactorUnavailable", <|
      "PlaneEulerExpressions" -> planeEuler|>], operationalFailure]];
  amplitudeEquations = engineSimplify[#/strippedFactor, assumptions] & /@
    planeEuler;
  matrixA = Table[engineSimplify[
      Coefficient[amplitudeEquations[[row]], amplitudes[[column]]], assumptions],
    {row, dimension}, {column, dimension}];
  emitCell[package, dimension, "M_ROUTE_A_STRIPPED_FACTOR", strippedFactor];
  emitCell[package, dimension, "M_A", matrixA];

  planeGradientRules = Thread[gradientVariables -> Flatten[Table[
      D[amplitudes[[column]] Cos[phaseVariable], phaseVariable] wavevector[[row]],
      {row, dimension}, {column, dimension}]]];
  planeVelocityRules = Thread[velocityJet -> Table[
      D[amplitudes[[component]] Cos[phaseVariable], phaseVariable] (-omegaLinear),
      {component, dimension}]];
  planeLagrangian = lagrangianJet /.
    Join[planeGradientRules, planeVelocityRules];
  averagedLagrangian = engineSimplify[
    Integrate[Expand[planeLagrangian], {phaseVariable, 0, 2 Pi}]/(2 Pi) /.
      omegaLinear^2 -> omegaSquared, assumptions];
  matrixB = Table[engineSimplify[
      D[averagedLagrangian, amplitudes[[row]], amplitudes[[column]]],
      assumptions],
    {row, dimension}, {column, dimension}];
  emitCell[package, dimension, "M_B", matrixB];
  matrixResidual = Map[engineSimplify[#, assumptions] &, matrixA - matrixB, {2}];
  emitCell[package, dimension, "M_RESIDUAL", matrixResidual];
  emitCell[package, dimension, "M_RATIO",
    engineSimplify[matrixA[[1, 1]]/matrixB[[1, 1]], assumptions]];
  emitCell[package, dimension, "M_ROUTE_RESIDUAL_SCOPE",
    "CODING_CONSISTENCY_ONLY"];
  routeSelection = <|"TOKEN" -> "M_B", "MATRIX" -> matrixB|>;
  downstreamMatrix = routeSelection["MATRIX"];
  emitCell[package, dimension, "M_ROUTE_USED", routeSelection["TOKEN"]];

  coefficientOrdering = coefficientOrderingFromRecords[
    kineticRecordList, stiffnessRecordList];
  emitCell[package, dimension, "COEFFICIENT_ORDERING", coefficientOrdering];
  matrixCoefficientJacobian = Table[
    Map[engineSimplify[#, assumptions] &,
      D[downstreamMatrix, coefficient], {2}],
    {coefficient, coefficientOrdering}];
  emitCell[package, dimension, "M_COEFFICIENT_JACOBIAN",
    matrixCoefficientJacobian];

  actionRecords = Join[kineticRecordList,
    Map[Append[#, "Factor" -> -#["Factor"]] &, stiffnessRecordList]];
  dimensionData = buildDimensions[coefficientOrdering, actionRecords,
    velocityJet, gradientVariables, amplitudes, wavevector, assumptions,
    package, dimension];

  spectrum = computeSpectrumAndModes[downstreamMatrix, coefficientOrdering,
    wavevector, assumptions, package, dimension, "", False, True, False, {},
    wavevector, coefficientOrdering];
  rootDimensionObjects = Table[
    dimensionOfScalar[spectrum["Roots"][[rootIndex]]/Total[wavevector^2],
      dimensionData["AtomDimensions"], assumptions],
    {rootIndex, Length[spectrum["Roots"]]}];
  Do[emitCell[package, dimension,
      "ROOT" <> ToString[rootIndex] <> "_DIM_OVER_KSQ",
      rootDimensionObjects[[rootIndex]]],
    {rootIndex, Length[rootDimensionObjects]}];
  emitScaledRoots[spectrum, wavevector, assumptions, package, dimension];
  genericRootJacobian = Table[
    engineSimplify[D[spectrum["Roots"][[rootIndex]], coefficient], assumptions],
    {rootIndex, Length[spectrum["Roots"]]},
    {coefficient, coefficientOrdering}];
  emitCell[package, dimension, "ROOT_COEFFICIENT_JACOBIAN",
    genericRootJacobian];

  strata = emitStrata[spectrum["Candidates"], downstreamMatrix,
    genericRootJacobian, coefficientOrdering, wavevector, assumptions,
    package, dimension];

  q7Residual = If[dimension == 3,
    independentGradient = Partition[invariantData["GVariables"], dimension];
    emitQ7[stiffnessRecordList, coordinateGradient, independentGradient,
      package, dimension],
    tokenNotApplicable];

  kwObjects = emitQ11[spectrum, coefficientOrdering, wavevector, amplitudes,
    assumptions, package, dimension];

  derivedHomogeneity = Join[
    {homogeneityRecord["DET_M", spectrum["Determinant"],
      dimensionData["AtomDimensions"], assumptions]},
    Table[homogeneityRecord["ROOT" <> ToString[rootIndex] <> "_VALUE",
      spectrum["Roots"][[rootIndex]], dimensionData["AtomDimensions"],
      assumptions], {rootIndex, Length[spectrum["Roots"]]}],
    Table[homogeneityRecord["ROOT" <> ToString[rootIndex] <> "_N5_M_DOT_K",
      spectrum["Q4Derived"][[rootIndex]]["MDotK"],
      dimensionData["AtomDimensions"], assumptions],
      {rootIndex, Length[spectrum["Q4Derived"]]}],
    Table[homogeneityRecord["ROOT" <> ToString[rootIndex] <> "_N6_RESIDUAL",
      spectrum["Q4Derived"][[rootIndex]]["BasisResiduals"],
      dimensionData["AtomDimensions"], assumptions],
      {rootIndex, Length[spectrum["Q4Derived"]]}],
    If[dimension == 3,
      {homogeneityRecord["Q7_RESIDUAL", q7Residual,
        dimensionData["AtomDimensions"], assumptions]}, {}],
    Table[homogeneityRecord["ROOT" <> ToString[rootIndex] <> "_KW_SQUARED",
      kwObjects[[rootIndex]], dimensionData["AtomDimensions"], bulkAssumptions],
      {rootIndex, Length[kwObjects]}]
  ];
  emitCell[package, dimension, "DIM_HOMOGENEITY_DERIVED", derivedHomogeneity];
];

packageOrder = {
  "MAIN", "XFORM_CURLONLY", "XFORM_EXTRA", "XFORM_DIVONLY",
  "XFORM_TRACELESS", "XCOEF_BSCALE", "XCOEF_BSIGN", "XKIN_ANISO"};
declaredSweep = <|
  "MAIN" -> {2, 3, 4, 5},
  "XFORM_CURLONLY" -> {2, 3, 4, 5},
  "XFORM_EXTRA" -> {2, 3, 4, 5},
  "XFORM_DIVONLY" -> {3, 4},
  "XFORM_TRACELESS" -> {3, 4},
  "XCOEF_BSCALE" -> {3},
  "XCOEF_BSIGN" -> {3},
  "XKIN_ANISO" -> {2, 3, 4, 5}|>;
sweepDimensions = Sort[DeleteDuplicates[Flatten[Values[declaredSweep]]]];
declaredPairs = Reap[
    Do[If[MemberQ[declaredSweep[package], dimension],
      Sow[{package, dimension}]],
      {dimension, sweepDimensions}, {package, packageOrder}]][[2, 1]];

scriptArguments = Rest[$ScriptCommandLine];
selectedPairs = Which[
  scriptArguments === {}, declaredPairs,
  Length[scriptArguments] == 2 && MemberQ[packageOrder, First[scriptArguments]] &&
      StringMatchQ[Last[scriptArguments], DigitCharacter ..],
    selectedDimension = FromDigits[Last[scriptArguments]];
    If[MemberQ[declaredSweep[First[scriptArguments]], selectedDimension],
      {{First[scriptArguments], selectedDimension}}, Quit[64]],
  True, Quit[64]
];

runPairs = {};
invariantCache = <||>;
operationalResult = Catch[
  Do[
    package = pair[[1]];
    dimension = pair[[2]];
    If[! KeyExistsQ[invariantCache, dimension],
      invariant = buildInvariantCensus[dimension, package];
      AssociateTo[invariantCache, dimension -> invariant],
      invariant = invariantCache[dimension];
      emitInvariantCensus[invariant, package, dimension]
    ];
    runCell[package, dimension, invariant];
    AppendTo[runPairs, {package, dimension}],
    {pair, selectedPairs}];
  completedNormally,
  operationalFailure
];

If[Head[operationalResult] === Failure, Quit[70]];
emit[runEmission["RUN_PAIRS"], runPairs];
emit[runEmission["SKIPPED_PAIRS"], Complement[declaredPairs, runPairs]];
localInventoryName = standardEmissionName[localInventoryEmission[]];
emit[localInventoryEmission[], Append[localNames, localInventoryName]];
