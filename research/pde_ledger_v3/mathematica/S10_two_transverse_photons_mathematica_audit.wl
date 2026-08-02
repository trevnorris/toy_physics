(* Independent Wolfram Language audit of the antisymmetric-derivative brane
   Lagrangian.  All displayed roots, bases, ranks/nullities, orientations,
   sums, and control results are generated below from the quadratic action. *)

ClearAll["Global`*"];
$HistoryLength = 0;

truthToken[value_] := If[TrueQ[value], "TRUE", "FALSE"];
inputString[value_] := ToString[value, InputForm];

dimensions = Range[2, 5];
dThree = SelectFirst[dimensions, IntegerQ[#] && ToString[#] === "3" &];

curlSquared[gradient_] := Module[{dimension = Length[gradient]},
  Expand[
    Total[Flatten[
      Table[(gradient[[i, j]] - gradient[[j, i]])^2,
        {i, dimension}, {j, dimension}]
    ]]/2
  ]
];

fullGradientSquared[gradient_] := Expand[Total[Flatten[gradient^2]]];

stiffnessSquared[gradient_, form_] := Switch[form,
  "Antisymmetric", curlSquared[gradient],
  "FullGradient", fullGradientSquared[gradient],
  _, $Failed
];

(* Verify the relation between the dimension-independent two-form norm and
   the ordinary three-dimensional curl norm using independent gradient
   components. *)
gradientAtThree = Array[gradientComponent, {dThree, dThree}];
generalCurlNormAtThree = curlSquared[gradientAtThree];
ordinaryCurlVector = {
  gradientAtThree[[2, 3]] - gradientAtThree[[3, 2]],
  gradientAtThree[[3, 1]] - gradientAtThree[[1, 3]],
  gradientAtThree[[1, 2]] - gradientAtThree[[2, 1]]
};
ordinaryCurlNorm = Expand[ordinaryCurlVector . ordinaryCurlVector];
curlReducesAtThree = TrueQ[FullSimplify[
  generalCurlNormAtThree == ordinaryCurlNorm
]];

caseAssumptions[wavevector_] := And[
  rhoBr > 0,
  muR > 0,
  Total[wavevector^2] > 0,
  And @@ (Element[#, Reals] & /@ Join[wavevector, {rhoBr, muR}])
];

(* This is the Euler-Lagrange variation of the coordinate-space density. *)
eulerLagrangeEquations[dimension_, form_, coefficientScale_] := Module[
  {coordinates, fieldHeads, fields, gradient, kineticDensity,
   lagrangianDensity, equations},

  coordinates = Table[Symbol["x" <> ToString[i]], {i, dimension}];
  fieldHeads = Table[Symbol["u" <> ToString[i]], {i, dimension}];
  fields = (#[Sequence @@ Join[coordinates, {t}]] &) /@ fieldHeads;
  gradient = Table[D[fields[[j]], coordinates[[i]]],
    {i, dimension}, {j, dimension}];
  kineticDensity = (rhoBr/2) Total[(D[#, t] & /@ fields)^2];
  lagrangianDensity = Expand[
    kineticDensity - (coefficientScale muR/2) stiffnessSquared[gradient, form]
  ];

  equations = Table[
    FullSimplify[
      D[D[lagrangianDensity, D[fields[[component]], t]], t] +
      Sum[
        D[D[lagrangianDensity,
          D[fields[[component]], coordinates[[direction]]]],
          coordinates[[direction]]],
        {direction, dimension}
      ] - D[lagrangianDensity, fields[[component]]]
    ] == 0,
    {component, dimension}
  ];
  equations
];

(* Vary the plane-wave quadratic action with respect to every amplitude.
   The real-phase quadratic form has the same coefficient matrix as
   a Exp[I (k.x - omega t)]. *)
planeWaveMatrix[dimension_, form_, coefficientScale_] := Module[
  {wavevector, amplitude, gradientAmplitude, potentialDensity,
   kineticFrequencyDensity},

  wavevector = Array[k, dimension];
  amplitude = Array[a, dimension];
  gradientAmplitude = Outer[Times, wavevector, amplitude];
  potentialDensity = Expand[
    (coefficientScale muR/2) stiffnessSquared[gradientAmplitude, form]
  ];
  kineticFrequencyDensity = Expand[
    (rhoBr omegaSquared/2) Total[amplitude^2]
  ];
  Table[
    D[potentialDensity - kineticFrequencyDensity,
      amplitude[[i]], amplitude[[j]]],
    {i, dimension}, {j, dimension}
  ]
];

sameExpressionQ[left_, right_, assumptions_] := TrueQ[
  FullSimplify[left == right, assumptions]
];

distinctRoots[matrix_, assumptions_] := Module[{solutions, roots},
  solutions = Solve[Det[matrix] == 0, omegaSquared];
  roots = FullSimplify[omegaSquared /. solutions, assumptions];
  DeleteDuplicates[roots, sameExpressionQ[#1, #2, assumptions] &]
];

vectorOrientation[vector_, wavevector_, assumptions_] := Module[
  {dotProduct, parallelResidual, isParallel, isPerpendicular},

  dotProduct = FullSimplify[vector . wavevector, assumptions];
  parallelResidual = FullSimplify[
    Total[wavevector^2] vector - dotProduct wavevector,
    assumptions
  ];
  isParallel = TrueQ[FullSimplify[
    And @@ Thread[parallelResidual == ConstantArray[0, Length[wavevector]]],
    assumptions
  ]];
  isPerpendicular = TrueQ[FullSimplify[dotProduct == 0, assumptions]];

  Which[
    isParallel, "PARALLEL_TO_K",
    isPerpendicular, "PERPENDICULAR_TO_K",
    True, "NEITHER"
  ]
];

rootAnalysis[matrix_, root_, wavevector_, assumptions_] := Module[
  {matrixAtRoot, basis, nullity, rank, vectorOrientations,
   rootOrientation, annihilationCheck, rankNullityCheck},

  matrixAtRoot = FullSimplify[matrix /. omegaSquared -> root, assumptions];
  basis = FullSimplify[NullSpace[matrixAtRoot], assumptions];
  nullity = Length[basis];
  rank = MatrixRank[matrixAtRoot];
  vectorOrientations =
    vectorOrientation[#, wavevector, assumptions] & /@ basis;
  rootOrientation = If[
    Length[DeleteDuplicates[vectorOrientations]] == 1,
    First[vectorOrientations],
    "NEITHER"
  ];
  annihilationCheck = And @@ (
    TrueQ[FullSimplify[
      And @@ Thread[matrixAtRoot . # == ConstantArray[0, Length[wavevector]]],
      assumptions
    ]] & /@ basis
  );
  rankNullityCheck = TrueQ[rank + nullity == Length[wavevector]];

  <|
    "OmegaSquared" -> root,
    "MatrixAtRoot" -> matrixAtRoot,
    "Rank" -> rank,
    "Basis" -> basis,
    "Nullity" -> nullity,
    "VectorOrientations" -> vectorOrientations,
    "Orientation" -> rootOrientation,
    "BasisAnnihilated" -> annihilationCheck,
    "RankNullity" -> rankNullityCheck
  |>
];

analyzeCase[dimension_, form_, coefficientScale_] := Module[
  {wavevector, assumptions, equations, matrix, roots, analyses,
   nullitySum, rootChecks, determinantChecks},

  wavevector = Array[k, dimension];
  assumptions = caseAssumptions[wavevector];
  equations = eulerLagrangeEquations[dimension, form, coefficientScale];
  matrix = FullSimplify[
    planeWaveMatrix[dimension, form, coefficientScale], assumptions
  ];
  roots = distinctRoots[matrix, assumptions];
  analyses = rootAnalysis[matrix, #, wavevector, assumptions] & /@ roots;
  nullitySum = Total[Lookup[analyses, "Nullity"]];
  determinantChecks = And @@ (
    TrueQ[FullSimplify[Det[matrix /. omegaSquared -> #] == 0,
      assumptions]] & /@ roots
  );
  rootChecks = determinantChecks &&
    And @@ Lookup[analyses, "BasisAnnihilated"] &&
    And @@ Lookup[analyses, "RankNullity"];

  <|
    "Dimension" -> dimension,
    "Form" -> form,
    "CoefficientScale" -> coefficientScale,
    "Wavevector" -> wavevector,
    "Assumptions" -> assumptions,
    "EquationOfMotion" -> equations,
    "Matrix" -> matrix,
    "Roots" -> roots,
    "RootAnalyses" -> analyses,
    "NullitySum" -> nullitySum,
    "NullitySumEqualsDimension" -> TrueQ[nullitySum == dimension],
    "RootChecks" -> rootChecks
  |>
];

emitCase[prefix_, case_, includeDynamics_] := Module[
  {analyses = case["RootAnalyses"], tag},

  If[TrueQ[includeDynamics],
    Print[prefix <> "_EQUATION_OF_MOTION: " <>
      inputString[case["EquationOfMotion"]]];
    Print[prefix <> "_DYNAMICAL_MATRIX: " <>
      inputString[case["Matrix"]]];
  ];
  Print[prefix <> "_ROOTS: " <> inputString[case["Roots"]]];
  Do[
    tag = prefix <> "_ROOT_" <> ToString[rootIndex];
    Print[tag <> "_OMEGA_SQUARED: " <>
      inputString[analyses[[rootIndex, "OmegaSquared"]]]];
    Print[tag <> "_NULLITY: " <>
      inputString[analyses[[rootIndex, "Nullity"]]]];
    Print[tag <> "_BASIS: " <>
      inputString[analyses[[rootIndex, "Basis"]]]];
    Print[tag <> "_ORIENTATION: " <>
      analyses[[rootIndex, "Orientation"]]];
    , {rootIndex, Length[analyses]}];
  Print[prefix <> "_NULLITY_SUM: " <> inputString[case["NullitySum"]]];
  Print[prefix <> "_NULLITY_SUM_EQUALS_D: " <>
    truthToken[case["NullitySumEqualsDimension"]]];
];

baselineCases = analyzeCase[#, "Antisymmetric", 1] & /@ dimensions;
formControlCase = analyzeCase[dThree, "FullGradient", 1];
coefficientControlScale =
  Length[ordinaryCurlVector] - Boole[curlReducesAtThree];
coefficientControlCase =
  analyzeCase[dThree, "Antisymmetric", coefficientControlScale];

sameRootSetQ[firstCase_, secondCase_] := Module[
  {firstRoots, secondRoots, assumptions},
  firstRoots = firstCase["Roots"];
  secondRoots = secondCase["Roots"];
  assumptions = firstCase["Assumptions"];
  Length[firstRoots] == Length[secondRoots] &&
    And @@ (Function[root,
      AnyTrue[secondRoots, sameExpressionQ[root, #, assumptions] &]
    ] /@ firstRoots)
];

caseShape[case_] := ({
  #["Nullity"], #["VectorOrientations"]
} &) /@ case["RootAnalyses"];

baselineAtThree = SelectFirst[baselineCases,
  #["Dimension"] == dThree &];
formControlSensitive = !(
  sameRootSetQ[baselineAtThree, formControlCase] &&
  caseShape[baselineAtThree] === caseShape[formControlCase]
);
coefficientControlSensitive = !sameRootSetQ[
  baselineAtThree, coefficientControlCase
];

tableRows = ({
  "D" -> #["Dimension"],
  "Roots" -> (({
    "OmegaSquared" -> #["OmegaSquared"],
    "Nullity" -> #["Nullity"],
    "Orientation" -> #["Orientation"]
  } &) /@ #["RootAnalyses"])
} &) /@ baselineCases;

Print["WL_S10_CURL_REDUCES_AT_D3: " <> truthToken[curlReducesAtThree]];
Scan[
  emitCase["WL_S10_D" <> ToString[#1["Dimension"]], #1, True] &,
  baselineCases
];

Print["WL_S10_FORM_CONTROL_D3_FORM: FULL_GRADIENT_SQUARED"];
emitCase["WL_S10_FORM_CONTROL_D3", formControlCase, False];
Print["WL_S10_FORM_CONTROL_SENSITIVE: " <>
  truthToken[formControlSensitive]];

Print["WL_S10_COEFFICIENT_CONTROL_D3_SCALE: " <>
  inputString[coefficientControlScale]];
emitCase["WL_S10_COEFFICIENT_CONTROL_D3", coefficientControlCase, False];
Print["WL_S10_COEFFICIENT_CONTROL_SENSITIVE: " <>
  truthToken[coefficientControlSensitive]];

Print["WL_S10_TABLE: " <> inputString[tableRows]];

checks = Join[
  {{"CURL_REDUCTION_D3", curlReducesAtThree}},
  ({"D" <> ToString[#["Dimension"]] <> "_PIPELINE",
      #["RootChecks"] && #["NullitySumEqualsDimension"]} &) /@
    baselineCases,
  {
    {"FORM_CONTROL_PIPELINE",
      formControlCase["RootChecks"] &&
      formControlCase["NullitySumEqualsDimension"]},
    {"FORM_CONTROL_SENSITIVITY", formControlSensitive},
    {"COEFFICIENT_CONTROL_PIPELINE",
      coefficientControlCase["RootChecks"] &&
      coefficientControlCase["NullitySumEqualsDimension"]},
    {"COEFFICIENT_CONTROL_SENSITIVITY", coefficientControlSensitive}
  }
];
failedChecks = Select[checks, !TrueQ[Last[#]] &];

If[failedChecks === {},
  Print["WL_S10_VERDICT: PASS"],
  Print["WL_S10_VERDICT: FAIL " <>
    StringRiffle[First /@ failedChecks, ","]];
  Exit[1]
];
