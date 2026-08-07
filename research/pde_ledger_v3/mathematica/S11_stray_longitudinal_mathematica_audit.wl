$HistoryLength = 0;

ClearAll["Global`*"];
$HistoryLength = 0;
$Messages = {OutputStream["stderr", 2]};

sScale = s;
bulkAmplitude = A;
ansatzProfile = Cos;
bulkSpeedPremise = cs0 > 0;

emittedNames = <||>;
localNames = {};
localTagNamesTag = "WL_S11_LOCAL_TAG_NAMES";
engineMethods = <|
  "Simplifier" -> FullSimplify,
  "CoefficientSortKey" -> SymbolName
|>;

emit[tag_String, payload_] := Module[{rendered},
  If[KeyExistsQ[emittedNames, tag], Quit[2]];
  AssociateTo[emittedNames, tag -> True];
  If[StringStartsQ[tag, "WL_S11_LOCAL_"], AppendTo[localNames, tag]];
  rendered = ToString[payload, InputForm];
  WriteString[$Output[[1]], tag <> ": " <> rendered <> "\n"];
  Flush[$Output[[1]]];
];

simp[expr_, assumptions_] :=
  engineMethods["Simplifier"][expr, Assumptions -> assumptions];
simpUnrestricted[expr_] := engineMethods["Simplifier"][expr];
heldMethod[method_] := HoldForm[method];
zeroTest[assumptions_] := Function[z, TrueQ[simp[z == 0, assumptions]]];
assumedRank[matrix_, assumptions_] :=
  MatrixRank[matrix, ZeroTest -> zeroTest[assumptions]];
assumedNullSpace[matrix_, assumptions_] :=
  NullSpace[matrix, ZeroTest -> zeroTest[assumptions]];

relationResidual[relation_] := Which[
  TrueQ[relation], 0,
  TrueQ[relation === False], 1,
  Head[relation] === Equal, Subtract @@ List @@ relation,
  True, relation
];

realConditionFromRelations[relations_List] := Module[{residuals},
  residuals = relationResidual /@ relations;
  And @@ Flatten[({ComplexExpand[Re[#]] == 0, ComplexExpand[Im[#]] == 0} &) /@ residuals]
];

multiplicativeFactors[expr_] := Module[{factored = Factor[expr]},
  If[Head[factored] === Times, List @@ factored, {factored}]
];

routeAStrippedFactor[planeEquations_, amplitudes_, phaseCoordinates_List] := Module[
  {coefficientEntries, nonzeroEntries, commonFactors, phaseFactors},
  coefficientEntries = Flatten[Table[
    Coefficient[planeEquations[[i]], amplitudes[[j]]],
    {i, Length[planeEquations]}, {j, Length[amplitudes]}]];
  nonzeroEntries = Select[coefficientEntries, ! TrueQ[# === 0] &];
  commonFactors = If[nonzeroEntries === {}, {},
    Fold[Intersection[#1, #2, SameTest -> SameQ] &,
      multiplicativeFactors[First[nonzeroEntries]],
      multiplicativeFactors /@ Rest[nonzeroEntries]]];
  phaseFactors = Select[commonFactors,
    ! FreeQ[#, Alternatives @@ phaseCoordinates] &];
  If[phaseFactors === {}, Missing["NoCommonPhaseFactor"], Times @@ phaseFactors]
];

stratumDerivativeFailure[roots_, coefficients_, equations_, point_] :=
  Failure["UndefinedStratumDerivative", <|
    "RootsRecomputedOnStratum" -> roots,
    "CoefficientOrdering" -> coefficients,
    "DefiningEquations" -> equations,
    "EvaluationPoint" -> point,
    "MissingConstruction" -> {TangentCoordinates, OffStratumExtension}|>];

conditionalParts[expr_] := If[
  Head[expr] === ConditionalExpression,
  {expr[[1]], expr[[2]]},
  {expr, True}
];

fourWaySign[operand_, assumptions_] := Which[
  TrueQ[simp[operand > 0, assumptions]], POSITIVE,
  TrueQ[simp[operand < 0, assumptions]], NEGATIVE,
  TrueQ[simp[operand == 0, assumptions]], ZERO,
  True, UNDECIDED
];

globalSymbolsIn[expr_] := DeleteDuplicates[
  Cases[expr, z_Symbol /; Context[z] === "Global`", {0, Infinity}]
];

coordinateFrame[d_Integer] := Module[{x, t, heads, args, fields, gradient},
  x = Table[Symbol["x" <> ToString[i]], {i, d}];
  t = Symbol["t"];
  heads = Table[Symbol["u" <> ToString[j]], {j, d}];
  args = Join[x, {t}];
  fields = (#[Sequence @@ args] &) /@ heads;
  gradient = Table[D[fields[[j]], x[[i]]], {i, d}, {j, d}];
  <|"X" -> x, "T" -> t, "Heads" -> heads, "Arguments" -> args,
    "Fields" -> fields, "Gradient" -> gradient|>
];

curlDensity[gradient_] := Module[{d = Length[gradient]},
  1/2 Total[Flatten[Table[(gradient[[i, j]] - gradient[[j, i]])^2,
    {i, d}, {j, d}]]]
];

divDensity[gradient_] := Tr[gradient]^2;

symtlDensity[gradient_] := Module[{d = Length[gradient], stl},
  stl = (gradient + Transpose[gradient])/2 -
    Tr[gradient] IdentityMatrix[d]/d;
  Total[Flatten[stl^2]]
];

quadCoordinates[poly_, variables_, pairs_] := Map[
  Function[pair,
    If[pair[[1]] === pair[[2]],
      Coefficient[Expand[poly], variables[[pair[[1]]]], 2],
      Coefficient[
        Coefficient[Expand[poly], variables[[pair[[1]]]], 1],
        variables[[pair[[2]]]], 1]
    ]
  ],
  pairs
];

polyFromCoordinates[row_, monomials_] := Expand[Total[MapThread[Times, {row, monomials}]]];

eulerOperatorForDensity[poly_, gVariables_, frame_] := Module[
  {d = Length[frame["X"]], gradient = frame["Gradient"], x = frame["X"]},
  Table[
    Expand[Total[Table[
      D[(D[poly, gVariables[[(i - 1) d + j]]] /. Thread[gVariables -> Flatten[gradient]]), x[[i]]],
      {i, d}]]],
    {j, d}
  ]
];

buildInvariantCensus[d_Integer] := Module[
  {n, gVariables, gMatrix, pairs, monomials, generators, deltaVariables,
   actionRows, connectedConstraints, v1Raw, v1Basis, r0, reflectedG,
   reflectionRows, fullConstraints, v2Raw, v2Basis, v1Dim, v2Dim,
   basisPolynomials, reflectedPolynomials, pivots, coordinateActionRows,
   operator, v6Coordinates, v6Ambient, v6Basis, v6Dim, frame,
   eulerOperators, pdPolynomial},
  n = d^2;
  gVariables = Table[Symbol["g" <> ToString[i]], {i, n}];
  gMatrix = Partition[gVariables, d];
  pairs = Flatten[Table[{p, q}, {p, n}, {q, p, n}], 1];
  monomials = (gVariables[[#[[1]]]] gVariables[[#[[2]]]] &) /@ pairs;
  generators = Flatten[Table[
    SparseArray[{{i, j} -> 1, {j, i} -> -1}, {d, d}],
    {i, 1, d - 1}, {j, i + 1, d}], 1];
  connectedConstraints = If[generators === {},
    SparseArray[{}, {0, Length[monomials]}],
    SparseArray[Join @@ Table[
      deltaVariables = Flatten[generators[[q]].gMatrix - gMatrix.generators[[q]]];
      actionRows = Table[
        quadCoordinates[
          Expand[Total[MapThread[(D[monomials[[p]], #1] #2) &,
            {gVariables, deltaVariables}]]],
          gVariables, pairs],
        {p, Length[monomials]}];
      Transpose[actionRows],
      {q, Length[generators]}]]
  ];
  v1Raw = NullSpace[connectedConstraints];
  v1Basis = RowReduce[v1Raw];
  r0 = DiagonalMatrix[Join[{-1}, ConstantArray[1, d - 1]]];
  reflectedG = Expand[r0.gMatrix.Transpose[r0]];
  reflectionRows = Table[
    quadCoordinates[monomials[[p]] /. Thread[gVariables -> Flatten[reflectedG]],
      gVariables, pairs],
    {p, Length[monomials]}];
  fullConstraints = Join[connectedConstraints,
    SparseArray[Transpose[reflectionRows - IdentityMatrix[Length[monomials]]]]];
  v2Raw = NullSpace[fullConstraints];
  v2Basis = RowReduce[v2Raw];
  v1Dim = Length[v1Basis];
  v2Dim = Length[v2Basis];
  basisPolynomials = polyFromCoordinates[#, monomials] & /@ v1Basis;
  reflectedPolynomials = Expand[# /. Thread[gVariables -> Flatten[reflectedG]]] & /@
    basisPolynomials;
  If[v1Dim == 0,
    operator = {};
    v6Basis = {};
    v6Dim = 0,
    pivots = FirstPosition[#, 1][[1]] & /@ v1Basis;
    coordinateActionRows = Table[
      quadCoordinates[reflectedPolynomials[[i]], gVariables, pairs][[pivots]],
      {i, v1Dim}];
    operator = Transpose[coordinateActionRows];
    v6Coordinates = NullSpace[operator + IdentityMatrix[v1Dim]];
    v6Ambient = (#.v1Basis) & /@ v6Coordinates;
    v6Basis = If[v6Ambient === {}, {}, RowReduce[v6Ambient]];
    v6Dim = Length[v6Basis]
  ];
  frame = coordinateFrame[d];
  eulerOperators = eulerOperatorForDensity[#, gVariables, frame] & /@ basisPolynomials;
  pdPolynomial = polyFromCoordinates[#, monomials] & /@ v6Basis // Total // Expand;
  <|
    "D" -> d, "GVariables" -> gVariables, "GMatrix" -> gMatrix,
    "Pairs" -> pairs, "Monomials" -> monomials,
    "V1Basis" -> v1Basis, "V1Dim" -> v1Dim,
    "V2Basis" -> v2Basis, "V2Dim" -> v2Dim,
    "R0" -> r0, "V1Polynomials" -> basisPolynomials,
    "ReflectedPolynomials" -> reflectedPolynomials,
    "V5Euler" -> eulerOperators, "V6Operator" -> operator,
    "V6Basis" -> v6Basis, "V6Dim" -> v6Dim,
    "PDPolynomial" -> pdPolynomial|>
];

emitInvariantCensus[data_, prefix_String] := Module[
  {v1 = data["V1Polynomials"], reflected = data["ReflectedPolynomials"], r0 = data["R0"]},
  emit[prefix <> "_MONOMIAL_ORDERING", data["Monomials"]];
  emit[prefix <> "_V1_BASIS", data["V1Basis"]];
  emit[prefix <> "_V1_DIM", Length[data["V1Basis"]]];
  emit[prefix <> "_V2_BASIS", data["V2Basis"]];
  emit[prefix <> "_V2_DIM", Length[data["V2Basis"]]];
  emit[prefix <> "_V3_DIFFERENCE", Length[data["V1Basis"]] - Length[data["V2Basis"]]];
  emit[prefix <> "_R0_MATRIX", r0];
  emit[prefix <> "_R0_ORTHOGONALITY_RESIDUAL",
    Transpose[r0].r0 - IdentityMatrix[data["D"]]];
  emit[prefix <> "_R0_DETERMINANT", Det[r0]];
  emit[prefix <> "_V4_REFLECTED", reflected];
  emit[prefix <> "_V4_RESIDUAL", MapThread[Expand[#1 - #2] &, {reflected, v1}]];
  emit[prefix <> "_V4_SUM", MapThread[Expand[#1 + #2] &, {reflected, v1}]];
  emit[prefix <> "_V5_EULER_LAGRANGE", data["V5Euler"]];
  emit[prefix <> "_V6_OPERATOR", data["V6Operator"]];
  emit[prefix <> "_V6_BASIS", data["V6Basis"]];
  emit[prefix <> "_V6_DIM", Length[data["V6Basis"]]];
  emit[prefix <> "_V7_RESIDUAL",
    Length[data["V6Basis"]] - (Length[data["V1Basis"]] - Length[data["V2Basis"]])];
];

packageTermBlueprint[package_String] := Switch[package,
  "MAIN", {{muR/2, "curl"}, {bComp/2, "div"}},
  "XFORM_CURLONLY", {{muR/2, "curl"}},
  "XFORM_DIVONLY", {{bComp/2, "div"}},
  "XFORM_TRACELESS", {{muR/2, "curl"}, {muBr, "symtl"}},
  "XFORM_EXTRA", {{muR/2, "curl"}, {bComp/2, "div"}, {beta/2, "pd"}},
  "XCOEF_BSCALE", {{muR/2, "curl"}, {sScale bComp/2, "div"}},
  "XCOEF_BSIGN", {{muR/2, "curl"}, {-bComp/2, "div"}}
];

makeTermRecords[package_String, densityJet_, densityCoordinate_] := Module[{blueprint},
  blueprint = packageTermBlueprint[package];
  Map[
    Function[item, <|"Factor" -> item[[1]], "Kind" -> item[[2]],
      "DensityJet" -> densityJet[item[[2]]],
      "DensityCoordinate" -> densityCoordinate[item[[2]]]|>],
    blueprint
  ]
];

makePremises[package_String, d_Integer, k_, amplitudes_, sDimension_, fieldModel_] := Module[
  {common, domain, entries, joint, coefficientDimensionDeclarations},
  coefficientDimensionDeclarations = fieldModel["CoefficientDimensionDeclarations"];
  common = <|
    "PREMISE_RHO" -> rhoBr > 0,
    "PREMISE_MUR" -> muR > 0,
    "PREMISE_BCOMP" -> bComp > 0,
    "PREMISE_WAVEVECTOR" -> Total[k^2] > 0,
    "PREMISE_WAVEVECTOR_REAL" -> And @@ (Element[#, Reals] & /@ k),
    "PREMISE_AMPLITUDE_REAL" -> And @@ (Element[#, Reals] & /@ amplitudes),
    "PREMISE_DIMENSION" -> (Element[dimensionSymbol, Integers] && dimensionSymbol > 0)
  |>;
  domain = Switch[package,
    "XFORM_TRACELESS", <|"PREMISE_PACKAGE_DOMAIN" -> muBr > 0|>,
    "XFORM_EXTRA", <|"PREMISE_PACKAGE_DOMAIN" -> Element[beta, Reals]|>,
    "XCOEF_BSCALE", <|"PREMISE_PACKAGE_DOMAIN" ->
      And[sScale > 0, sScale != 1,
        And @@ Thread[sDimension == coefficientDimensionDeclarations[sScale]]]|>,
    _, <|"PREMISE_PACKAGE_DOMAIN" -> True|>
  ];
  entries = Join[common, domain];
  joint = And @@ Values[entries];
  <|"Entries" -> entries, "Joint" -> joint, "AtDimension" -> (joint /. dimensionSymbol -> d)|>
];

emitPremises[premises_, prefix_String, fieldModel_, actionModel_, bulkModel_] := Module[{},
  KeyValueMap[(emit[prefix <> "_" <> #1, #2]) &, premises["Entries"]];
  emit[prefix <> "_JOINT_ASSUMPTIONS", premises["Joint"]];
  emit[prefix <> "_PREMISE_U_DIMENSION", SuppliedPremise[fieldModel["UDimension"]]];
  emit[prefix <> "_PREMISE_IN_PLANE_COMPONENTS",
    SuppliedPremise[Length[fieldModel["Amplitudes"]]]];
  emit[prefix <> "_PREMISE_BACKGROUND_STATE",
    SuppliedPremiseWithNoInCodeConsumer[actionModel["BackgroundVelocity"]]];
  emit[prefix <> "_PREMISE_NO_DISSIPATION",
    SuppliedPremiseWithNoInCodeConsumer[actionModel["DissipativeTerms"]]];
  emit[prefix <> "_PREMISE_LINEAR_RESPONSE",
    SuppliedPremise[actionModel["MaximumFieldDegree"]]];
  emit[prefix <> "_PREMISE_FROZEN_WALL_WIDTH",
    SuppliedPremiseWithNoInCodeConsumer[actionModel["WallWidthFields"]]];
  emit[prefix <> "_PREMISE_BULK_MODE", SuppliedPremise[bulkModel["ModeContent"]]];
  emit[prefix <> "_PREMISE_BULK_SPEED", bulkModel["SpeedPremise"]];
  emit[prefix <> "_PREMISE_BULK_AMPLITUDE",
    SuppliedPremiseWithNoInCodeConsumer[bulkModel["AmplitudePremise"]]];
];

extractConditionalExpressions[expr_] := DeleteDuplicates[Cases[
  expr, ConditionalExpression[value_, condition_] :> {value, condition}, Infinity]];

branchEquations[branch_] := Module[{rules, conditions},
  rules = Cases[branch, Rule[left_, right_] :> (left == right), {1}];
  conditions = Cases[branch, ConditionalExpression[_, condition_] :> condition, Infinity];
  Join[rules, conditions]
];

emitLocus[sharedBase_String, localBase_String, equations_List, variables_List, assumptions_] := Module[
  {solution, residuals, identicallySatisfied, unrestrictedReduction, inconsistent,
   realAdmissible, conditions, operands, realOperands, quantifiedVariables,
   branchReduction, simplifiedBranchEquations},
  emit[sharedBase <> "_EQUATIONS", equations];
  solution = Quiet[Solve[And @@ equations, variables]];
  emit[sharedBase <> "_SOLUTION", SolveResult[variables, solution]];
  conditions = extractConditionalExpressions[solution];
  emit[localBase <> "_CONDITIONS", conditions];
  residuals = relationResidual /@ equations;
  identicallySatisfied = And @@ (TrueQ[simpUnrestricted[# == 0]] & /@ residuals);
  emit[sharedBase <> "_IDENTICALLY_SATISFIED", identicallySatisfied];
  unrestrictedReduction = Quiet[Reduce[And @@ equations, variables, Complexes]];
  inconsistent = TrueQ[unrestrictedReduction === False];
  emit[sharedBase <> "_INCONSISTENT", inconsistent];
  realAdmissible = If[ListQ[solution],
    Map[
      Function[branch,
        operands = And[assumptions, And @@ branchEquations[branch]];
        simplifiedBranchEquations = simp[#, assumptions] & /@ branchEquations[branch];
        realOperands = And[assumptions,
          realConditionFromRelations[simplifiedBranchEquations]];
        quantifiedVariables = globalSymbolsIn[realOperands];
        branchReduction = Quiet[Reduce[realOperands, quantifiedVariables, Reals]];
        {operands, realOperands, branchReduction,
          Not[SameQ[branchReduction, False]]}
      ],
      solution
    ],
    SolverBranchesUnavailable[solution]
  ];
  emit[sharedBase <> "_REAL_ADMISSIBLE", realAdmissible];
  <|"Equations" -> equations, "Variables" -> variables, "Solution" -> solution,
    "UnrestrictedReduction" -> unrestrictedReduction,
    "RealAdmissible" -> realAdmissible|>
];

multiplicitySolutions[determinant_] := Module[{numerator, factors, pieces},
  numerator = Numerator[Together[determinant]];
  factors = Select[FactorList[numerator], ! FreeQ[First[#], omegaSquared] &];
  pieces = Map[
    Function[factorPair,
      Flatten[ConstantArray[#, factorPair[[2]]] & /@
        Quiet[Solve[factorPair[[1]] == 0, omegaSquared]], 1]
    ],
    factors
  ];
  Flatten[pieces, 1]
];

distinctUnderAssumptions[values_, assumptions_] := DeleteDuplicates[
  values, TrueQ[simp[#1 == #2, assumptions]] &
];

allMaximalMinors[matrix_, rank_Integer] := Module[{rows, columns},
  If[rank == 0, Return[{}]];
  rows = Subsets[Range[Length[matrix]], {rank}];
  columns = Subsets[Range[Length[First[matrix]]], {rank}];
  Flatten[Table[Factor[Det[matrix[[rowSet, columnSet]]]],
    {rowSet, rows}, {columnSet, columns}]]
];

computeSpectrumModes[matrix_, coefficients_, k_, assumptions_, sharedPrefix_String,
    localPrefix_String, includeRankLoci_] := Module[
  {determinant, returnedSolutions, returnedValues, distinctRoots, rootConditions,
   pairEquations, kLocus, coefficientLocus, rootObjects, root, rootParts,
   rootPrefix, rootLocalPrefix, mr, rank, nullity, stacked, stackedRank,
   transverseNullity, basis, basisDots, basisResiduals, basisCount,
   minors, minorEquations, kRankLocus, coefficientRankLocus, jointRankLocus,
   rankLoci = {}, q4Derived = {}, pairs},
  determinant = Factor[Det[matrix]];
  emit[sharedPrefix <> "_DET_M", determinant];
  returnedSolutions = multiplicitySolutions[determinant];
  emit[sharedPrefix <> "_ROOT_SOLUTION_SET", returnedSolutions];
  rootConditions = extractConditionalExpressions[returnedSolutions];
  emit[localPrefix <> "_ROOT_SOLUTION_CONDITIONS", rootConditions];
  returnedValues = (omegaSquared /. #) & /@ returnedSolutions;
  distinctRoots = distinctUnderAssumptions[returnedValues, assumptions];
  emit[sharedPrefix <> "_ROOT_DISTINCT", distinctRoots];
  emit[sharedPrefix <> "_ROOT_COUNT_ALL", Length[returnedSolutions]];
  emit[sharedPrefix <> "_ROOT_COUNT_DISTINCT", Length[distinctRoots]];
  emit[sharedPrefix <> "_ROOT_ORDERING", distinctRoots];
  Do[
    rootParts = conditionalParts[distinctRoots[[r]]];
    rootPrefix = sharedPrefix <> "_ROOT" <> ToString[r];
    rootLocalPrefix = localPrefix <> "_ROOT" <> ToString[r];
    emit[rootPrefix <> "_VALUE", rootParts[[1]]];
    emit[rootLocalPrefix <> "_CONDITION", rootParts[[2]]];
    emit[rootPrefix <> "_SIGN", {rootParts[[1]], fourWaySign[rootParts[[1]], assumptions]}],
    {r, Length[distinctRoots]}];
  pairs = Subsets[Range[Length[distinctRoots]], {2}];
  pairEquations = (simp[distinctRoots[[#[[1]]]] - distinctRoots[[#[[2]]]], assumptions] == 0) & /@ pairs;
  kLocus = emitLocus[sharedPrefix <> "_ROOT_COINCIDENCE_K",
    localPrefix <> "_ROOT_COINCIDENCE_K", pairEquations, k, assumptions];
  coefficientLocus = emitLocus[sharedPrefix <> "_ROOT_COINCIDENCE_COEFF",
    localPrefix <> "_ROOT_COINCIDENCE_COEFF", pairEquations, coefficients, assumptions];
  rootObjects = Table[
    root = distinctRoots[[r]];
    rootPrefix = sharedPrefix <> "_ROOT" <> ToString[r];
    rootLocalPrefix = localPrefix <> "_ROOT" <> ToString[r];
    mr = Map[simp[#, assumptions] &, matrix /. omegaSquared -> root, {2}];
    emit[rootPrefix <> "_N1", mr];
    rank = assumedRank[mr, assumptions];
    emit[rootPrefix <> "_N2_RANK", rank];
    nullity = Length[mr] - rank;
    emit[rootPrefix <> "_N2_NULLITY", nullity];
    stacked = Join[mr, {k}];
    stackedRank = assumedRank[stacked, assumptions];
    emit[rootPrefix <> "_N3_STACKED_RANK", stackedRank];
    transverseNullity = Length[mr] - stackedRank;
    emit[rootPrefix <> "_N3_TRANSVERSE_NULLITY", transverseNullity];
    emit[rootPrefix <> "_N4_NULLITY_DIFFERENCE", nullity - transverseNullity];
    emit[rootPrefix <> "_N5_M_DOT_K", Map[simp[#, assumptions] &, mr.k]];
    basis = assumedNullSpace[mr, assumptions];
    emit[rootPrefix <> "_N6_BASIS", basis];
    basisDots = Map[simp[#.k, assumptions] &, basis];
    emit[rootPrefix <> "_N6_DOT_K", basisDots];
    basisResiduals = MapThread[
      Function[{vector, dot}, Map[simp[#, assumptions] &,
        Total[k^2] vector - dot k]],
      {basis, basisDots}];
    emit[rootPrefix <> "_N6_RESIDUAL", basisResiduals];
    basisCount = Length[basis];
    emit[rootPrefix <> "_N7_BASIS_COUNT", basisCount];
    emit[rootPrefix <> "_N7_RESIDUAL", basisCount - nullity];
    If[TrueQ[includeRankLoci],
      minors = allMaximalMinors[mr, rank];
      emit[rootPrefix <> "_RANK_DROP_MINORS", minors];
      minorEquations = (# == 0) & /@ minors;
      kRankLocus = emitLocus[rootPrefix <> "_RANK_DROP_K",
        rootLocalPrefix <> "_RANK_DROP_K", minorEquations, k, assumptions];
      coefficientRankLocus = emitLocus[rootPrefix <> "_RANK_DROP_COEFF",
        rootLocalPrefix <> "_RANK_DROP_COEFF", minorEquations, coefficients, assumptions];
      jointRankLocus = emitLocus[rootPrefix <> "_RANK_DROP_JOINT",
        rootLocalPrefix <> "_RANK_DROP_JOINT", minorEquations,
        Join[k, coefficients], assumptions];
      AppendTo[rankLoci, <|"RootIndex" -> r, "Joint" -> jointRankLocus|>]
    ];
    AppendTo[q4Derived, <|"MDotK" -> mr.k, "BasisResiduals" -> basisResiduals|>];
    <|"Root" -> root, "Matrix" -> mr, "Rank" -> rank,
      "Nullity" -> nullity, "Basis" -> basis,
      "BasisResiduals" -> basisResiduals|>,
    {r, Length[distinctRoots]}];
  <|"Determinant" -> determinant, "Solutions" -> returnedSolutions,
    "Roots" -> distinctRoots, "RootObjects" -> rootObjects,
    "RankLoci" -> rankLoci, "Q4Derived" -> q4Derived,
    "CoincidenceK" -> kLocus, "CoincidenceCoefficient" -> coefficientLocus|>
];

dimensionVectorsEqualQ[left_, right_] := If[ListQ[left] && ListQ[right] &&
    Length[left] == Length[right],
  And @@ (TrueQ[PossibleZeroQ[Together[#]]] & /@ (left - right)), False];

dimensionOfScalar[expr_, atomDimensions_, assumptions_] := Module[{head, pieces, dimensions, unique},
  Which[
    NumberQ[expr], {0, 0, 0},
    KeyExistsQ[atomDimensions, Unevaluated[expr]], atomDimensions[Unevaluated[expr]],
    Head[expr] === ConditionalExpression,
      dimensionOfScalar[expr[[1]], atomDimensions, assumptions],
    Head[expr] === Plus,
      pieces = List @@ expr;
      dimensions = dimensionOfScalar[#, atomDimensions, assumptions] & /@ pieces;
      unique = DeleteDuplicates[dimensions, dimensionVectorsEqualQ];
      If[Length[unique] == 1, First[unique], DimensionalAlternatives @@ unique],
    Head[expr] === Times,
      dimensions = dimensionOfScalar[#, atomDimensions, assumptions] & /@ (List @@ expr);
      If[And @@ (ListQ /@ dimensions), Total[dimensions], DimensionalProduct[dimensions]],
    Head[expr] === Power && (IntegerQ[expr[[2]]] || Head[expr[[2]]] === Rational),
      expr[[2]] dimensionOfScalar[expr[[1]], atomDimensions, assumptions],
    True, UnknownDimension[expr]
  ]
];

expandedTerms[expr_] := If[TrueQ[expr === 0], {},
  If[Head[Expand[expr]] === Plus, List @@ Expand[expr], {Expand[expr]}]];

fieldDegree[term_, fieldAtoms_List] := Total[Exponent[term, #] & /@ fieldAtoms];

retainThroughFieldDegree[expr_, fieldAtoms_List, maximumDegree_Integer] :=
  Expand[Total[Select[expandedTerms[expr],
    fieldDegree[#, fieldAtoms] <= maximumDegree &]]];

dimensionInventory[expr_, atomDimensions_, assumptions_] := Module[{scalars, termDimensions},
  scalars = If[ListQ[expr], Flatten[expr], {expr}];
  Map[
    Function[scalar,
      termDimensions = dimensionOfScalar[#, atomDimensions, assumptions] & /@ expandedTerms[scalar];
      {scalar, termDimensions,
        If[termDimensions === {}, True,
          And @@ (dimensionVectorsEqualQ[#, First[termDimensions]] & /@ termDimensions)]}
    ],
    scalars
  ]
];

buildDimensions[coefficientOrdering_, actionTerms_, jetAtoms_, fieldModel_, assumptions_, prefix_String] := Module[
  {dimensionDeclarations, dimensionless, unknownCoefficients, dimensionVariables,
   coefficientDimensions,
   atomDimensions, targetDimension, equationVectors, dimensionEquations,
   flatUnknowns, solution, solutionRules, finalCoefficientDimensions,
   firstSlotVariables, firstSlotExpressions, coefficientMatrix, equationCount,
   unknownCount, unresolved, determinacy, finalAtomDimensions, actionInventory},
  dimensionDeclarations = fieldModel["CoefficientDimensionDeclarations"];
  dimensionless = Keys[Select[dimensionDeclarations, SameQ[#, {0, 0, 0}] &]];
  unknownCoefficients = Select[coefficientOrdering, ! MemberQ[dimensionless, #] &];
  dimensionVariables = Association@Table[
    coefficient -> Table[
      Symbol["dim" <> SymbolName[coefficient] <> "slot" <> ToString[slot]],
      {slot, 3}],
    {coefficient, unknownCoefficients}];
  coefficientDimensions = Association@Table[
    coefficient -> If[MemberQ[dimensionless, coefficient], {0, 0, 0},
      dimensionVariables[coefficient]],
    {coefficient, coefficientOrdering}];
  atomDimensions = Join[coefficientDimensions, jetAtoms];
  targetDimension = {2 - dimensionSymbol, -2, 1};
  equationVectors = DeleteDuplicates[Flatten[
    (dimensionOfScalar[#, atomDimensions, assumptions] & /@ expandedTerms[#]) & /@ actionTerms,
    1]];
  dimensionEquations = Flatten[(Thread[# == targetDimension] &) /@ equationVectors];
  flatUnknowns = Flatten[Values[dimensionVariables]];
  solution = Quiet[Solve[dimensionEquations, flatUnknowns]];
  emit[prefix <> "_DIM_EQUATIONS", dimensionEquations];
  emit[prefix <> "_DIM_SOLUTION", solution];
  solutionRules = If[solution === {}, {}, First[solution]];
  finalCoefficientDimensions = Table[
    {coefficient, simp[coefficientDimensions[coefficient] /. solutionRules, assumptions]},
    {coefficient, coefficientOrdering}];
  emit[prefix <> "_DIM_COEFFICIENTS", finalCoefficientDimensions];
  firstSlotVariables = If[unknownCoefficients === {}, {},
    dimensionVariables[#][[1]] & /@ unknownCoefficients];
  firstSlotExpressions = If[equationVectors === {}, {},
    (#[[1]] - targetDimension[[1]]) & /@ equationVectors];
  coefficientMatrix = If[firstSlotExpressions === {} || firstSlotVariables === {},
    ConstantArray[0, {Length[firstSlotExpressions], Length[firstSlotVariables]}],
    Normal[Last[CoefficientArrays[firstSlotExpressions, firstSlotVariables]]]];
  equationCount = If[coefficientMatrix === {} || Length[firstSlotVariables] == 0,
    0, MatrixRank[coefficientMatrix]];
  unknownCount = Length[unknownCoefficients];
  unresolved = If[solution === {}, flatUnknowns,
    Intersection[flatUnknowns,
      Cases[Values[Association[solutionRules]], z_Symbol, Infinity]]];
  determinacy = Which[
    solution === {}, OVER_DETERMINED,
    unresolved =!= {} || Length[solutionRules] < Length[flatUnknowns], UNDER_DETERMINED,
    True, EXACTLY_DETERMINED
  ];
  emit[prefix <> "_DIM_EQUATION_COUNT", equationCount];
  emit[prefix <> "_DIM_UNKNOWN_COUNT", unknownCount];
  emit[prefix <> "_DIM_COUNT_DIFFERENCE", equationCount - unknownCount];
  emit[prefix <> "_DIM_DETERMINACY", determinacy];
  finalAtomDimensions = atomDimensions /. solutionRules;
  actionInventory = dimensionInventory[#, finalAtomDimensions, assumptions] & /@ actionTerms;
  emit[prefix <> "_DIM_HOMOGENEITY_ACTION", actionInventory];
  <|"CoefficientDimensions" -> finalCoefficientDimensions,
    "AtomDimensions" -> finalAtomDimensions, "Solution" -> solution,
    "SolutionRules" -> solutionRules|>
];

strataFromRankLoci[rankLoci_, assumptions_, pointVariables_] := Module[
  {branchPairs, branches, equationSets, uniqueSets, instances},
  branchPairs = Flatten[Map[
    Function[item,
      If[ListQ[item["Joint"]["Solution"]] && ListQ[item["Joint"]["RealAdmissible"]],
        Transpose[{item["Joint"]["Solution"], item["Joint"]["RealAdmissible"]}], {}]],
    rankLoci], 1];
  branches = First /@ Select[branchPairs,
    ! TrueQ[Last[Last[#]] === False] &];
  equationSets = branchEquations /@ branches;
  uniqueSets = DeleteDuplicates[equationSets,
    Sort[ToString[#, InputForm] & /@ #1] === Sort[ToString[#, InputForm] & /@ #2] &];
  instances = Map[
    Function[equations,
      Quiet[FindInstance[And[assumptions, And @@ equations], pointVariables, Reals, 1]]],
    uniqueSets];
  Cases[MapThread[{#1, #2} &, {uniqueSets, instances}],
    {equations_, instance_} /; ListQ[instance] && Length[instance] > 0 &&
      ListQ[First[instance]] && And @@ (MatchQ[#, _Rule] & /@ First[instance]) :>
      <|"Equations" -> equations, "Point" -> First[instance]|>]
];

emitScaledRoots[spectrum_, k_, assumptions_, prefix_String] := Module[
  {root, scaled, ratio, numeratorExponent, denominatorExponent, candidate, exponent},
  Do[
    root = spectrum["Roots"][[r]];
    scaled = simp[root /. Thread[k -> lambdaScale k],
      And[assumptions, lambdaScale > 0]];
    emit[prefix <> "_ROOT" <> ToString[r] <> "_SCALED", scaled];
    emit[prefix <> "_ROOT" <> ToString[r] <> "_UNSCALED", root];
    If[TrueQ[simp[root == 0, assumptions]],
      ratio = UNDEFINED_RATIO;
      exponent = UNDEFINED_RATIO,
      ratio = simp[scaled/root, And[assumptions, lambdaScale > 0]];
      numeratorExponent = Exponent[Numerator[Together[ratio]], lambdaScale];
      denominatorExponent = Exponent[Denominator[Together[ratio]], lambdaScale];
      candidate = numeratorExponent - denominatorExponent;
      exponent = If[TrueQ[simp[ratio == lambdaScale^candidate,
        And[assumptions, lambdaScale > 0]]], candidate, NOT_A_PURE_POWER]
    ];
    emit[prefix <> "_ROOT" <> ToString[r] <> "_SCALE_RATIO", ratio];
    emit[prefix <> "_ROOT" <> ToString[r] <> "_SCALE_EXPONENT", exponent],
    {r, Length[spectrum["Roots"]]}]
];

emitQ7[packageRecords_, coordinateDensities_, prefix_String] := Module[
  {gMatrix, d, curlTerm, curlDensityValue, cVector, reference, fullW},
  gMatrix = coordinateDensities["IndependentGradient"];
  d = Length[gMatrix];
  fullW = Total[(#["Factor"] (#["DensityCoordinate"] /.
      Thread[Flatten[coordinateDensities["Gradient"]] -> Flatten[gMatrix]]) &) /@ packageRecords];
  curlTerm = Total[Map[
    Function[record,
      If[TrueQ[record["Kind"] === "curl"],
        record["Factor"] (record["DensityCoordinate"] /.
          Thread[Flatten[coordinateDensities["Gradient"]] -> Flatten[gMatrix]]), 0]],
    packageRecords]];
  curlDensityValue = coordinateDensities["Curl"] /.
    Thread[Flatten[coordinateDensities["Gradient"]] -> Flatten[gMatrix]];
  cVector = Table[Total[Flatten[Table[
    Signature[{i, j, k}] gMatrix[[j, k]], {j, d}, {k, d}]]], {i, d}];
  reference = Expand[cVector.cVector];
  emit[prefix <> "_Q7_W_FULL", Expand[fullW]];
  emit[prefix <> "_Q7_CURL_TERM", Expand[curlTerm]];
  emit[prefix <> "_Q7_CURL_DENSITY", Expand[curlDensityValue]];
  emit[prefix <> "_Q7_CURL_REFERENCE", reference];
  emit[prefix <> "_Q7_RESIDUAL", Expand[curlDensityValue - reference]];
  Expand[curlDensityValue - reference]
];

declaredSweep = <|
  "MAIN" -> {2, 3, 4, 5},
  "XFORM_CURLONLY" -> {2, 3, 4, 5},
  "XFORM_EXTRA" -> {2, 3, 4, 5},
  "XFORM_DIVONLY" -> {3, 4},
  "XFORM_TRACELESS" -> {3, 4},
  "XCOEF_BSCALE" -> {3},
  "XCOEF_BSIGN" -> {3}
|>;
declaredPairs = Flatten[KeyValueMap[Function[{package, dimensions},
  ({package, #} &) /@ dimensions], declaredSweep], 1];
runPairs = {};
invariantCache = <||>;

Do[
  invariantCache[d] = buildInvariantCensus[d];
  Do[
    If[MemberQ[declaredPairs, {package, d}],
      Module[
        {prefix, localPrefix, invariant, frame, x, t, uFields, gCoordinate,
         gVariables, gJet, vJet, amplitudes, k, phi, omegaLinear,
         fieldModel, actionModel, bulkModel, sDimension, premises,
         assumptions, bulkAssumptions, densityJet, densityCoordinate, pdJet, pdCoordinate,
         termRecords, stiffnessTermsCoordinate, stiffnessTermsJet, wCoordinate,
         kineticCoordinate, lagrangianCoordinate, kineticJet, lagrangianJet,
         jetRules, eomExpressions, eomRelations, phase, ansatzRules,
         eomPlane, strippedFactor, amplitudeEquations, matrixA,
         planeGradientRules, planeVelocityRules, planeLagrangian,
         averagedLagrangian, matrixB, routeRecords, downstreamRecord,
         downstreamMatrix, coefficientFactors, coefficientOrdering,
         dimensionlessCoefficients, jetAtomDimensions, dimensionData,
         spectrum, rootJacobian, scaledData, rankLoci, pointVariables,
         strata, stratum, stratumPrefix, stratumLocalPrefix, point,
         stratumSpectrum, restrictedJacobian,
         recomputedJacobian, q7Residual, q11Data, bulkDispersion,
         kwEquation, kwSolutionRules, kwSolutionExpression, kwParts,
         bulkModeRecords, suppliedEquationInventory, closureEquations, closureUnknowns,
         closureMatrix, closurePadded, closureRref, closureRank, q11Derived, derivedInventory,
         rootDimOverKsq, actionTermsForDimension},
        prefix = "WL_S11_" <> package <> "_D" <> ToString[d];
        localPrefix = "WL_S11_LOCAL_" <> package <> "_D" <> ToString[d];
        invariant = invariantCache[d];
        frame = coordinateFrame[d];
        x = frame["X"];
        t = frame["T"];
        uFields = frame["Fields"];
        gCoordinate = frame["Gradient"];
        gVariables = invariant["GVariables"];
        gJet = Partition[gVariables, d];
        vJet = Table[Symbol["vt" <> ToString[j]], {j, d}];
        amplitudes = Table[Symbol["a" <> ToString[j]], {j, d}];
        k = Table[Symbol["k" <> ToString[i]], {i, d}];
        phi = Symbol["phaseVariable"];
        omegaLinear = Symbol["omegaLinear"];
        sDimension = Table[Symbol["dimsslot" <> ToString[i]], {i, 3}];
        fieldModel = <|"Amplitudes" -> amplitudes, "UDimension" -> {1, 0, 0},
          "AnsatzProfile" -> ansatzProfile,
          "CoefficientDimensionDeclarations" -> If[package === "XCOEF_BSCALE",
            <|sScale -> {0, 0, 0}|>, <||>]|>;
        actionModel = <|"BackgroundVelocity" -> 0, "DissipativeTerms" -> {},
          "MaximumFieldDegree" -> 2, "WallWidthFields" -> {}|>;
        bulkModel = <|"ModeContent" -> {ScalarSoundMode},
          "SpeedPremise" -> bulkSpeedPremise,
          "AmplitudePremise" -> Element[bulkAmplitude, Reals]|>;
        premises = makePremises[package, d, k, amplitudes, sDimension, fieldModel];
        assumptions = premises["AtDimension"];
        bulkAssumptions = And[assumptions, bulkModel["SpeedPremise"]];
        emit[localPrefix <> "_SIMPLIFIER", heldMethod[engineMethods["Simplifier"]]];
        emit[localPrefix <> "_COEFFICIENT_SORT_KEY",
          heldMethod[engineMethods["CoefficientSortKey"]]];
        emitPremises[premises, prefix, fieldModel, actionModel, bulkModel];
        emitInvariantCensus[invariant, prefix];
        pdJet = invariant["PDPolynomial"];
        pdCoordinate = pdJet /. Thread[gVariables -> Flatten[gCoordinate]];
        emit[prefix <> "_PD_TERM", Expand[pdCoordinate]];
        densityJet = <|"curl" -> curlDensity[gJet], "div" -> divDensity[gJet],
          "symtl" -> symtlDensity[gJet], "pd" -> pdJet|>;
        densityCoordinate = <|"curl" -> curlDensity[gCoordinate],
          "div" -> divDensity[gCoordinate], "symtl" -> symtlDensity[gCoordinate],
          "pd" -> pdCoordinate|>;
        termRecords = Map[
          Function[record, With[
            {retainedDensity = retainThroughFieldDegree[record["DensityJet"],
                Join[vJet, gVariables], actionModel["MaximumFieldDegree"]]},
            <|"Factor" -> record["Factor"], "Kind" -> record["Kind"],
              "DensityJet" -> retainedDensity,
              "DensityCoordinate" -> retainedDensity /.
                Thread[gVariables -> Flatten[gCoordinate]]|>]],
          makeTermRecords[package, densityJet, densityCoordinate]];
        stiffnessTermsCoordinate = (#["Factor"] #["DensityCoordinate"]) & /@ termRecords;
        stiffnessTermsJet = (#["Factor"] #["DensityJet"]) & /@ termRecords;
        wCoordinate = Total[stiffnessTermsCoordinate];
        kineticJet = retainThroughFieldDegree[rhoBr/2 Total[vJet^2],
          Join[vJet, gVariables], actionModel["MaximumFieldDegree"]];
        jetRules = Join[Thread[vJet -> (D[#, t] & /@ uFields)],
          Thread[gVariables -> Flatten[gCoordinate]]];
        kineticCoordinate = kineticJet /. jetRules;
        lagrangianCoordinate = Expand[kineticCoordinate - wCoordinate];
        emit[prefix <> "_LAGRANGIAN", lagrangianCoordinate];
        emit[prefix <> "_STIFFNESS_TERMS", stiffnessTermsCoordinate];
        lagrangianJet = Expand[kineticJet - Total[stiffnessTermsJet]];
        eomExpressions = Table[Expand[
          D[(D[lagrangianJet, vJet[[j]]] /. jetRules), t] +
          Total[Table[D[(D[lagrangianJet, gJet[[i, j]]] /. jetRules), x[[i]]], {i, d}]]],
          {j, d}];
        eomRelations = (# == 0) & /@ eomExpressions;
        emit[prefix <> "_EULER_LAGRANGE_SYSTEM", eomRelations];
        phase = Total[k x] - omegaLinear t;
        ansatzRules = Table[With[{head = frame["Heads"][[j]], amp = amplitudes[[j]],
            args = frame["Arguments"], phaseValue = phase,
            profile = fieldModel["AnsatzProfile"]},
          head -> Function[Evaluate[args], Evaluate[amp profile[phaseValue]]]], {j, d}];
        eomPlane = TrigFactor[eomExpressions /. ansatzRules] /.
          omegaLinear^2 -> omegaSquared;
        strippedFactor = routeAStrippedFactor[eomPlane, amplitudes, Join[x, {t}]];
        amplitudeEquations = Map[simp[TrigFactor[#/strippedFactor], assumptions] &, eomPlane];
        matrixA = Table[Coefficient[amplitudeEquations[[i]], amplitudes[[j]]], {i, d}, {j, d}];
        emit[prefix <> "_M_ROUTE_A_STRIPPED_FACTOR", strippedFactor];
        emit[prefix <> "_M_A", matrixA];
        planeGradientRules = Thread[gVariables -> Flatten[Table[
          D[amplitudes[[j]] fieldModel["AnsatzProfile"][phi], phi] k[[i]],
          {i, d}, {j, d}]]];
        planeVelocityRules = Thread[vJet -> Table[
          D[amplitudes[[j]] fieldModel["AnsatzProfile"][phi], phi] (-omegaLinear),
          {j, d}]];
        planeLagrangian = lagrangianJet /. Join[planeGradientRules, planeVelocityRules];
        averagedLagrangian = simp[
          Integrate[Expand[planeLagrangian], {phi, 0, 2 Pi}]/(2 Pi) /.
            omegaLinear^2 -> omegaSquared, assumptions];
        matrixB = Table[simp[D[averagedLagrangian, amplitudes[[i]], amplitudes[[j]]], assumptions],
          {i, d}, {j, d}];
        emit[prefix <> "_M_B", matrixB];
        emit[prefix <> "_M_RESIDUAL", Map[simp[#, assumptions] &, matrixA - matrixB, {2}]];
        emit[prefix <> "_M_RATIO", simp[matrixA[[1, 1]]/matrixB[[1, 1]], assumptions]];
        routeRecords = {{EulerLagrangeRoute, matrixA}, {QuadraticFormRoute, matrixB}};
        emit[prefix <> "_M_RESIDUAL_SCOPE",
          CodingConsistency[First /@ routeRecords]];
        downstreamRecord = Last[routeRecords];
        downstreamMatrix = Last[downstreamRecord];
        emit[prefix <> "_M_ROUTE_USED", First[downstreamRecord]];
        coefficientFactors = Prepend[(#["Factor"] &) /@ termRecords, rhoBr/2];
        coefficientOrdering = SortBy[
          DeleteDuplicates[Flatten[globalSymbolsIn /@ coefficientFactors]],
          engineMethods["CoefficientSortKey"]];
        emit[prefix <> "_COEFFICIENT_ORDERING", coefficientOrdering];
        emit[prefix <> "_M_COEFFICIENT_JACOBIAN",
          Table[Map[simp[#, assumptions] &, D[downstreamMatrix, coefficient], {2}],
            {coefficient, coefficientOrdering}]];
        jetAtomDimensions = Association@Join[
          Thread[vJet -> ConstantArray[fieldModel["UDimension"] + {0, -1, 0}, d]],
          Thread[gVariables -> ConstantArray[fieldModel["UDimension"] + {-1, 0, 0}, d^2]],
          Thread[k -> ConstantArray[{-1, 0, 0}, d]],
          Thread[amplitudes -> ConstantArray[fieldModel["UDimension"], d]],
          {omegaSquared -> {0, -2, 0}, cs0 -> {1, -1, 0}, kwSquared -> {-2, 0, 0},
           lambdaScale -> {0, 0, 0}}
        ];
        actionTermsForDimension = Prepend[(-#) & /@ stiffnessTermsJet, kineticJet];
        dimensionData = buildDimensions[coefficientOrdering, actionTermsForDimension,
          jetAtomDimensions, fieldModel, assumptions, prefix];
        spectrum = computeSpectrumModes[downstreamMatrix, coefficientOrdering, k,
          assumptions, prefix, localPrefix, True];
        rootDimOverKsq = Table[
          dimensionOfScalar[spectrum["Roots"][[r]]/Total[k^2],
            dimensionData["AtomDimensions"], assumptions],
          {r, Length[spectrum["Roots"]]}];
        Do[emit[prefix <> "_ROOT" <> ToString[r] <> "_DIM_OVER_KSQ", rootDimOverKsq[[r]]],
          {r, Length[rootDimOverKsq]}];
        emitScaledRoots[spectrum, k, assumptions, prefix];
        rootJacobian = Table[
          simp[D[spectrum["Roots"][[r]], coefficient], assumptions],
          {r, Length[spectrum["Roots"]]}, {coefficient, coefficientOrdering}];
        emit[prefix <> "_ROOT_COEFFICIENT_JACOBIAN", rootJacobian];
        rankLoci = spectrum["RankLoci"];
        pointVariables = DeleteDuplicates@Join[globalSymbolsIn[assumptions],
          coefficientOrdering];
        strata = strataFromRankLoci[rankLoci, assumptions, pointVariables];
        strata = SortBy[strata, ToString[#["Equations"], InputForm] &];
        emit[prefix <> "_STRATUM_ORDERING", (#["Equations"] &) /@ strata];
        Do[
          stratum = strata[[stratumIndex]];
          stratumPrefix = prefix <> "_STRATUM" <> ToString[stratumIndex];
          stratumLocalPrefix = localPrefix <> "_STRATUM" <> ToString[stratumIndex];
          point = stratum["Point"];
          emit[stratumPrefix <> "_DEFINING_EQUATIONS", stratum["Equations"]];
          emit[stratumPrefix <> "_POINT", point];
          stratumSpectrum = computeSpectrumModes[downstreamMatrix /. point,
            coefficientOrdering, k, assumptions /. point, stratumPrefix,
            stratumLocalPrefix, False];
          restrictedJacobian = rootJacobian /. point;
          emit[stratumPrefix <> "_ROOT_COEFFICIENT_JACOBIAN_RESTRICTED", restrictedJacobian];
          recomputedJacobian = stratumDerivativeFailure[
            stratumSpectrum["Roots"], coefficientOrdering,
            stratum["Equations"], point];
          emit[stratumPrefix <> "_ROOT_COEFFICIENT_JACOBIAN_RECOMPUTED", recomputedJacobian],
          {stratumIndex, Length[strata]}];
        q7Residual = If[d == 3,
          emitQ7[termRecords,
            <|"Gradient" -> gCoordinate, "IndependentGradient" -> Partition[gVariables, d],
              "Curl" -> densityCoordinate["curl"]|>, prefix],
          NotApplicableQ7[d]];
        bulkModeRecords = Map[
          Function[mode, Switch[mode,
            ScalarSoundMode, <|
              "Amplitude" -> bulkAmplitude,
              "FieldEquation" -> (BulkScalarField[x, Symbol["w"], t] ==
                bulkAmplitude Cos[Total[k x] + kw Symbol["w"] - omegaLinear t]),
              "Dispersion" -> (omegaSquared == cs0^2 (Total[k^2] + kwSquared))|>,
            _, Failure["UnsupportedBulkMode", <|"Mode" -> mode|>]]],
          bulkModel["ModeContent"]];
        bulkDispersion = First[Lookup[bulkModeRecords, "Dispersion"]];
        q11Derived = Table[
          kwEquation = bulkDispersion /. omegaSquared -> spectrum["Roots"][[r]];
          emit[prefix <> "_ROOT" <> ToString[r] <> "_KW_EQUATION", kwEquation];
          kwSolutionRules = Quiet[Solve[kwEquation, kwSquared]];
          kwSolutionExpression = kwSquared /. First[kwSolutionRules];
          kwParts = conditionalParts[kwSolutionExpression];
          emit[prefix <> "_ROOT" <> ToString[r] <> "_KW_SQUARED", kwParts[[1]]];
          emit[localPrefix <> "_ROOT" <> ToString[r] <> "_KW_CONDITION", kwParts[[2]]];
          emit[prefix <> "_ROOT" <> ToString[r] <> "_KW_SIGN",
            {kwParts[[1]], fourWaySign[kwParts[[1]], bulkAssumptions]}];
          emitLocus[prefix <> "_ROOT" <> ToString[r] <> "_KW_ZERO_LOCUS",
            localPrefix <> "_ROOT" <> ToString[r] <> "_KW_ZERO_LOCUS",
            {kwParts[[1]] == 0}, Join[coefficientOrdering, {cs0}],
            bulkAssumptions];
          kwParts[[1]],
          {r, Length[spectrum["Roots"]]}];
        suppliedEquationInventory = Join[
          Table[uFields[[j]] == amplitudes[[j]] fieldModel["AnsatzProfile"][phase],
            {j, d}],
          Flatten[({#["FieldEquation"], #["Dispersion"]} &) /@ bulkModeRecords]];
        closureEquations = Select[suppliedEquationInventory,
          Function[equation, ! FreeQ[equation, bulkAmplitude] &&
            AnyTrue[amplitudes, Function[amplitude, ! FreeQ[equation, amplitude]]]]];
        closureUnknowns = Join[fieldModel["Amplitudes"],
          Lookup[bulkModeRecords, "Amplitude"]];
        closureMatrix = Table[Coefficient[relationResidual[equation], unknown],
          {equation, closureEquations}, {unknown, closureUnknowns}];
        closurePadded = Join[closureMatrix, {ConstantArray[0, Length[closureUnknowns]]}];
        closureRref = RowReduce[closurePadded];
        closureRank = Length[Select[closureRref,
          ! And @@ (zeroTest[bulkAssumptions] /@ #) &]];
        emit[prefix <> "_C1_EQUATIONS", closureEquations];
        emit[prefix <> "_C2_UNKNOWNS", closureUnknowns];
        emit[prefix <> "_C2_COUNT", Length[closureUnknowns]];
        emit[prefix <> "_C3_RANK", closureRank];
        emit[prefix <> "_C4_DIFFERENCE", Length[closureUnknowns] - closureRank];
        derivedInventory = Join[
          {dimensionInventory[spectrum["Determinant"], dimensionData["AtomDimensions"], assumptions]},
          (dimensionInventory[#, dimensionData["AtomDimensions"], assumptions] & /@ spectrum["Roots"]),
          (dimensionInventory[#MDotK, dimensionData["AtomDimensions"], assumptions] & /@
            spectrum["Q4Derived"]),
          (dimensionInventory[#BasisResiduals, dimensionData["AtomDimensions"], assumptions] & /@
            spectrum["Q4Derived"]),
          If[d == 3, {dimensionInventory[q7Residual,
            dimensionData["AtomDimensions"], assumptions]}, {}],
          (dimensionInventory[#, dimensionData["AtomDimensions"], bulkAssumptions] & /@
            q11Derived)
        ];
        emit[prefix <> "_DIM_HOMOGENEITY_DERIVED", derivedInventory];
        AppendTo[runPairs, {package, d}];
      ]
    ],
    {package, Keys[declaredSweep]}],
  {d, Sort[DeleteDuplicates[Flatten[Values[declaredSweep]]]]}];

emit["WL_S11_RUN_PAIRS", runPairs];
emit["WL_S11_SKIPPED_PAIRS", Complement[declaredPairs, runPairs]];
emit[localTagNamesTag, Append[localNames, localTagNamesTag]];
