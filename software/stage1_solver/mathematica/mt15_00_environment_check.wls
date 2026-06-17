#!/usr/bin/env wolframscript

(* M0 Mathematica 15 environment/toolchain check.
   This script is intentionally diagnostic: absent features are report rows,
   not process failures. *)

$HistoryLength = 0;

ClearAll["Global`*"];

scriptPath =
  If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    $InputFileName,
    FileNameJoin[{Directory[], "software/stage1_solver/mathematica/mt15_00_environment_check.wls"}]
  ];
scriptDir = DirectoryName[scriptPath];
If[! DirectoryQ[scriptDir],
  CreateDirectory[scriptDir, CreateIntermediateDirectories -> True]
];
reportPath = FileNameJoin[{scriptDir, "environment_report.md"}];
workingInvocation =
  "timeout 600 wolframscript -script software/stage1_solver/mathematica/mt15_00_environment_check.wls";

results = {};

ClearAll[addResult, md, shortString, heldValue, badResultQ, safeEval, statusFromBool,
  functionSymbolExistsQ, globalShadowExistsQ, safeFileDelete];

addResult[id_, routeNeed_, probe_, status_, details_: ""] :=
  AppendTo[results,
    <|
      "id" -> ToString[id],
      "routeNeed" -> ToString[routeNeed],
      "probe" -> ToString[probe],
      "status" -> ToString[status],
      "details" -> ToString[details]
    |>
  ];

md[x_] := StringReplace[ToString[x], {
    "\r\n" -> "<br>",
    "\n" -> "<br>",
    "\r" -> "<br>",
    "|" -> "\\|"
  }];

shortString[x_, depth_: 3] :=
  StringReplace[ToString[Short[x, depth], InputForm], {
    "\r\n" -> " ",
    "\n" -> " ",
    "\r" -> " "
  }];

heldValue[x_] := Apply[HoldComplete, {x}];

badResultQ[x_] := TrueQ[x === $Failed] || TrueQ[FailureQ[x]];

SetAttributes[safeEval, HoldFirst];
safeEval[expr_] := Quiet[Check[expr, $Failed]];

statusFromBool[ok_, caveat_: ""] :=
  If[TrueQ[ok],
    If[StringLength[caveat] > 0, "PRESENT (caveat: " <> caveat <> ")", "PRESENT"],
    "ABSENT"
  ];

functionSymbolExistsQ[name_String] := Length[Names["System`" <> name]] > 0;
globalShadowExistsQ[name_String] := Length[Names["Global`" <> name]] > 0;

safeFileDelete[path_String] := If[FileExistsQ[path], Quiet[DeleteFile[path]]];

(* 1. Version and license. *)
versionOK = NumericQ[$VersionNumber] && $VersionNumber >= 15.0;
addResult[
  "1.1",
  "gating prerequisite",
  "$VersionNumber >= 15.0",
  statusFromBool[versionOK],
  "$VersionNumber = " <> ToString[$VersionNumber, InputForm]
];

addResult[
  "1.2",
  "gating prerequisite",
  "$Version / $MachineName / $ProcessorCount",
  "PRESENT",
  "$Version = " <> ToString[$Version] <>
    "; $MachineName = " <> ToString[$MachineName] <>
    "; $ProcessorCount = " <> ToString[$ProcessorCount, InputForm]
];

licenseOK = TrueQ[$MaxLicenseProcesses === 2];
addResult[
  "1.3",
  "gating prerequisite",
  "$MaxLicenseProcesses == 2",
  If[licenseOK, "PRESENT",
    "PRESENT (caveat: observed cap is not the expected 2-seat value)"
  ],
  "$MaxLicenseProcesses = " <> ToString[$MaxLicenseProcesses, InputForm]
];

(* 2. FEM / NDSolve. *)
femLoad = safeEval[Needs["NDSolve`FEM`"]; True];
femLoadOK = TrueQ[femLoad];
addResult[
  "2.1",
  "stationary open-throat + linearized BVPs",
  "Needs[\"NDSolve`FEM`\"]",
  statusFromBool[femLoadOK],
  "load result = " <> shortString[femLoad]
];

sol1D = safeEval[
  NDSolveValue[
    {Derivative[2][u][x] == -1, u[0] == 0, u[1] == 0},
    u,
    {x, 0, 1},
    Method -> {"FiniteElement"}
  ]
];
val1D = If[badResultQ[sol1D], $Failed, safeEval[N[sol1D[0.5]]]];
fem1DFunctionOK = ! badResultQ[sol1D] && TrueQ[NumericQ[val1D]];
fem1DValueOK = fem1DFunctionOK && Abs[val1D - 0.125] < 10^-5;
addResult[
  "2.2",
  "stationary open-throat + linearized BVPs",
  "1-D FEM BVP u''[x] == -1 with Dirichlet endpoints",
  If[fem1DFunctionOK,
    If[fem1DValueOK, "PRESENT",
      "PRESENT (caveat: returned a function, but center value was not within the tight smoke tolerance)"
    ],
    "ABSENT"
  ],
  "Head = " <> ToString[Head[sol1D], InputForm] <>
    "; u(0.5) = " <> shortString[val1D]
];

sol2D = safeEval[
  NDSolveValue[
    {
      Laplacian[v[x, y], {x, y}] == -1,
      DirichletCondition[v[x, y] == 0, True]
    },
    v,
    Element[{x, y}, Rectangle[{0, 0}, {1, 1}]],
    Method -> {"FiniteElement"}
  ]
];
val2D = If[badResultQ[sol2D], $Failed, safeEval[N[sol2D[0.5, 0.5]]]];
fem2DOK = ! badResultQ[sol2D] && TrueQ[NumericQ[val2D]];
addResult[
  "2.3",
  "stationary open-throat + linearized BVPs",
  "2-D unit-square Poisson FEM smoke",
  statusFromBool[fem2DOK],
  "Head = " <> ToString[Head[sol2D], InputForm] <>
    "; u(0.5,0.5) = " <> shortString[val2D]
];

(* 3. Coordinate / curvilinear support. *)
chartExists = functionSymbolExistsQ["CoordinateChartData"];
chartData = If[chartExists,
  safeEval[CoordinateChartData["Spherical", "CoordinateRangeAssumptions"]],
  $Failed
];
chartOK = chartExists && ! badResultQ[chartData];
addResult[
  "3.1",
  "keep throat coordinates natural",
  "CoordinateChartData[\"Spherical\", ...]",
  statusFromBool[chartOK],
  "System symbol exists = " <> ToString[chartExists] <>
    "; sample = " <> shortString[chartData]
];

lap = safeEval[Laplacian[f[r, th, ph], {r, th, ph}, "Spherical"]];
grad = safeEval[Grad[f[r, th, ph], {r, th, ph}, "Spherical"]];
div = safeEval[
  Div[
    {ar[r, th, ph], ath[r, th, ph], aph[r, th, ph]},
    {r, th, ph},
    "Spherical"
  ]
];
lapOK = ! badResultQ[lap] && FreeQ[heldValue[lap], Laplacian, Heads -> True];
gradOK = ! badResultQ[grad] && FreeQ[heldValue[grad], Grad, Heads -> True];
divOK = ! badResultQ[div] && FreeQ[heldValue[div], Div, Heads -> True];
addResult[
  "3.2",
  "keep throat coordinates natural",
  "Laplacian/Grad/Div with \"Spherical\" coordinate spec",
  statusFromBool[lapOK && gradOK && divOK],
  "Laplacian ok = " <> ToString[lapOK] <>
    "; Grad ok = " <> ToString[gradOK] <>
    "; Div ok = " <> ToString[divOK]
];

(* 4. Eigensolvers. *)
denseEig = safeEval[Eigensystem[{{2., 1.}, {1., 3.}}]];
denseEigOK =
  ListQ[denseEig] && Length[denseEig] == 2 &&
    VectorQ[denseEig[[1]], NumericQ];
addResult[
  "4.1",
  "BdG / mixed modal extraction",
  "Eigensystem on a small dense matrix",
  statusFromBool[denseEigOK],
  "eigenvalues = " <> If[denseEigOK, shortString[denseEig[[1]]], shortString[denseEig]]
];

genEig = safeEval[Eigenvalues[{{{2., 0.}, {0., 3.}}, {{1., 0.}, {0., 2.}}}]];
genEigOK = VectorQ[genEig, NumericQ] && Length[genEig] == 2;
addResult[
  "4.2",
  "BdG / mixed modal extraction",
  "Generalized eigenproblem Eigenvalues[{A, B}]",
  statusFromBool[genEigOK],
  "eigenvalues = " <> shortString[genEig]
];

nEig = 20;
sparseEigA =
  SparseArray[
    {Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.},
    {nEig, nEig}
  ];
arnoldiEig = safeEval[Eigenvalues[sparseEigA, 3, Method -> {"Arnoldi"}]];
arnoldiEigOK = VectorQ[arnoldiEig, NumericQ] && Length[arnoldiEig] == 3;
addResult[
  "4.3",
  "BdG / mixed modal extraction",
  "Sparse Arnoldi Eigenvalues[sparseA, 3, Method -> {\"Arnoldi\"}]",
  statusFromBool[arnoldiEigOK],
  "eigenvalues = " <> shortString[arnoldiEig]
];

(* 5. Sparse LinearSolve methods. *)
ClearAll[linearSolveProbe];
nLin = 8;
linA =
  SparseArray[
    {Band[{1, 1}] -> 2., Band[{1, 2}] -> -1., Band[{2, 1}] -> -1.},
    {nLin, nLin}
  ];
linX = N[Range[nLin]/nLin];
linB = linA.linX;

linearSolveProbe[method_String] := Module[{sol, err, ok, detail},
  sol = safeEval[LinearSolve[linA, linB, Method -> method]];
  ok = VectorQ[sol, NumericQ] &&
    (err = Norm[sol - linX]/Max[Norm[linX], 1.0]; err < 10^-8);
  detail = If[TrueQ[ok],
    "relative error = " <> ToString[err, InputForm],
    "result = " <> shortString[sol]
  ];
  {ok, detail}
];

linearMethods = {"Multifrontal", "Pardiso", "Krylov", "Banded", "MUMPS",
  "DefinitelyNotAMethod"};
linearProbeData = Association[];
Do[
  Module[{probe = linearSolveProbe[m], caveat = ""},
    linearProbeData[m] = probe;
    If[m === "MUMPS" && TrueQ[probe[[1]]],
      caveat = "small sparse smoke only; no production-scale or linked-library benchmark"
    ];
    addResult[
      "5." <> m,
      "scalable linear algebra bottleneck",
      "LinearSolve[A, b, Method -> \"" <> m <> "\"]",
      If[TrueQ[probe[[1]]],
        If[StringLength[caveat] > 0, "PRESENT (caveat: " <> caveat <> ")", "PRESENT"],
        "ABSENT"
      ],
      probe[[2]]
    ];
  ],
  {m, linearMethods}
];

(* 6. Matrix decompositions. *)
ClearAll[decompositionArgs, decompositionProbe];
mat = {{4., 1.}, {2., 3.}};
spdMat = {{2., 1.}, {1., 2.}};
indefMat = {{0., 1.}, {1., 0.}};

decompositionArgs[name_String] := Switch[name,
  "CholeskyDecomposition" | "LDLDecomposition", {spdMat},
  "BunchKaufmanDecomposition", {indefMat},
  _, {mat}
];

decompositionProbe[name_String] := Module[
  {exists, shadow, sym, args, res, unevaluated, ok, detail},
  exists = functionSymbolExistsQ[name];
  shadow = globalShadowExistsQ[name];
  If[! exists,
    Return[{False, "defined = no; tiny matrix smoke = no"}]
  ];
  sym = ToExpression["System`" <> name];
  args = decompositionArgs[name];
  res = safeEval[Apply[sym, args]];
  unevaluated = ! badResultQ[res] && ! FreeQ[heldValue[res], sym, Heads -> True];
  ok = ! badResultQ[res] && ! unevaluated;
  detail =
    "defined = yes; tiny matrix smoke = " <> If[ok, "yes", "no"] <>
    "; Global shadow = " <> If[shadow, "yes", "no"] <>
    "; result head = " <> ToString[Head[res], InputForm] <>
    If[ok, "", "; result = " <> shortString[res]];
  {ok, detail}
];

decompNames = {
  "LUDecomposition",
  "CholeskyDecomposition",
  "SchurDecomposition",
  "JordanDecomposition",
  "QRDecomposition",
  "SingularValueDecomposition",
  "LDLDecomposition",
  "BunchKaufmanDecomposition",
  "OrderedSchurDecomposition",
  "PolarDecomposition",
  "PopovDecomposition",
  "JordanReduce",
  "FrobeniusReduce",
  "PartialFractions",
  "PolynomialReduction"
};
decompProbeData = Association[];
Do[
  Module[{probe = decompositionProbe[name]},
    decompProbeData[name] = probe;
    addResult[
      "6." <> name,
      "stability classification + pole extraction",
      name <> " on a tiny matrix",
      statusFromBool[probe[[1]]],
      probe[[2]]
    ];
  ],
  {name, decompNames}
];

(* 7. GPU. *)
cudaLoad = safeEval[Needs["CUDALink`"]; True];
cudaLoadOK = TrueQ[cudaLoad];
addResult[
  "7.1",
  "nice-to-have GPU acceleration",
  "Needs[\"CUDALink`\"]",
  statusFromBool[cudaLoadOK],
  "load result = " <> shortString[cudaLoad]
];

gpuArraySymbol = functionSymbolExistsQ["GPUArray"];
gpuArrayConstruct = If[gpuArraySymbol,
  safeEval[GPUArray[{1., 2.}]],
  $Failed
];
gpuArrayOK =
  gpuArraySymbol && ! badResultQ[gpuArrayConstruct] &&
    FreeQ[heldValue[gpuArrayConstruct], GPUArray, Heads -> True];
addResult[
  "7.2",
  "nice-to-have GPU acceleration",
  "GPUArray symbol and small construction smoke",
  If[gpuArrayOK, "PRESENT",
    If[gpuArraySymbol,
      "ABSENT",
      "ABSENT"
    ]
  ],
  "System symbol exists = " <> ToString[gpuArraySymbol] <>
    "; construction result = " <> shortString[gpuArrayConstruct]
];

targetDeviceLinearSolve =
  safeEval[LinearSolve[linA, linB, TargetDevice -> "GPU"]];
targetDeviceOK = VectorQ[targetDeviceLinearSolve, NumericQ] &&
  Norm[targetDeviceLinearSolve - linX]/Max[Norm[linX], 1.0] < 10^-8;
addResult[
  "7.3",
  "nice-to-have GPU acceleration",
  "LinearSolve[A, b, TargetDevice -> \"GPU\"]",
  statusFromBool[targetDeviceOK],
  "result = " <> shortString[targetDeviceLinearSolve]
];

hybridProbe = linearSolveProbe["HybridCPUGPU"];
addResult[
  "7.4",
  "nice-to-have GPU acceleration",
  "LinearSolve[A, b, Method -> \"HybridCPUGPU\"]",
  statusFromBool[hybridProbe[[1]]],
  hybridProbe[[2]]
];

(* 8. Export. *)
jsonPath = FileNameJoin[{scriptDir, "mt15_00_json_probe.json"}];
yamlPath = FileNameJoin[{scriptDir, "mt15_00_yaml_probe.yaml"}];
jsonAssoc = <|"alpha" -> 1, "beta" -> {True, "x"}|>;
safeFileDelete[jsonPath];
jsonExport = safeEval[Export[jsonPath, jsonAssoc, "JSON"]];
jsonImport = If[FileExistsQ[jsonPath], safeEval[Import[jsonPath, "RawJSON"]], $Failed];
jsonOK = AssociationQ[jsonImport] && jsonImport === jsonAssoc;
addResult[
  "8.1",
  "V2-22B packet handoff",
  "Export/Import round trip through .json",
  statusFromBool[jsonOK],
  "Export result = " <> shortString[jsonExport] <>
    "; Import result = " <> shortString[jsonImport]
];
safeFileDelete[jsonPath];

yamlListed = MemberQ[$ExportFormats, "YAML"];
safeFileDelete[yamlPath];
yamlExport = If[yamlListed, safeEval[Export[yamlPath, jsonAssoc, "YAML"]], $Failed];
yamlOK = yamlListed && ! badResultQ[yamlExport] && FileExistsQ[yamlPath];
addResult[
  "8.2",
  "human/LLM-readable adjacent exports only",
  "YAML listed in $ExportFormats and Export[..., \"YAML\"]",
  statusFromBool[yamlOK],
  "$ExportFormats contains YAML = " <> ToString[yamlListed] <>
    "; Export result = " <> shortString[yamlExport]
];
safeFileDelete[yamlPath];

(* 9. ModelFit. *)
linearFit = safeEval[LinearModelFit[{{0, 1}, {1, 3}, {2, 5}}, x, x]];
linearFitValue = If[badResultQ[linearFit], $Failed, safeEval[N[linearFit[3]]]];
linearFitOK = Head[linearFit] === FittedModel && NumericQ[linearFitValue];
addResult[
  "9.1",
  "post-freeze diagnostics only",
  "LinearModelFit smoke",
  statusFromBool[linearFitOK],
  "Head = " <> ToString[Head[linearFit], InputForm] <>
    "; fitted value at x=3 = " <> shortString[linearFitValue]
];

nonlinearFit =
  safeEval[
    NonlinearModelFit[
      {{0, 1}, {1, N[E]}, {2, N[E^2]}, {3, N[E^3]}},
      a Exp[b x],
      {{a, 0.8}, {b, 0.9}},
      x
    ]
  ];
nonlinearFitValue = If[badResultQ[nonlinearFit], $Failed, safeEval[N[nonlinearFit[1]]]];
nonlinearFitOK = Head[nonlinearFit] === FittedModel && NumericQ[nonlinearFitValue];
addResult[
  "9.2",
  "post-freeze diagnostics only",
  "NonlinearModelFit smoke",
  statusFromBool[nonlinearFitOK],
  "Head = " <> ToString[Head[nonlinearFit], InputForm] <>
    "; fitted value at x=1 = " <> shortString[nonlinearFitValue]
];

(* Summary construction. *)
ClearAll[statusById, prefixStatuses, anyAbsentQ, anyCaveatQ, compactStatus];
statusById[id_String] := Lookup[FirstCase[results, a_ /; a["id"] === id :> a, <||>], "status", "ABSENT"];
prefixStatuses[prefix_String] := Lookup[Select[results, StringStartsQ[#["id"], prefix] &], "status"];
anyAbsentQ[statuses_] := AnyTrue[statuses, # === "ABSENT" &];
anyCaveatQ[statuses_] := AnyTrue[statuses, StringStartsQ[#, "PRESENT (caveat:"] &];
compactStatus[statuses_, caveat_String: "some sub-probes have caveats"] := Which[
  Length[statuses] == 0, "ABSENT",
  AllTrue[statuses, # === "PRESENT" &], "PRESENT",
  anyAbsentQ[statuses] && AnyTrue[statuses, StringStartsQ[#, "PRESENT"] &],
    "PRESENT (caveat: " <> caveat <> ")",
  anyAbsentQ[statuses], "ABSENT",
  anyCaveatQ[statuses], "PRESENT (caveat: " <> caveat <> ")",
  True, "PRESENT"
];

linearRequiredStatuses =
  Lookup[
    Select[results, MemberQ["5." <> # & /@ {"Multifrontal", "Pardiso", "Krylov", "Banded", "MUMPS"}, #["id"]] &],
    "status"
  ];
decompCoreStatuses =
  Lookup[
    Select[results, StringStartsQ[#["id"], "6."] &],
    "status"
  ];

summaryRows = {
  {"2", "FEM / NDSolve", compactStatus[prefixStatuses["2."]]},
  {"3", "Coordinate / curvilinear support", compactStatus[prefixStatuses["3."]]},
  {"4", "Eigensolvers", compactStatus[prefixStatuses["4."]]},
  {"5", "Sparse LinearSolve methods",
    compactStatus[linearRequiredStatuses,
      "small sparse smoke only; unknown-method sanity probe is excluded from this aggregate"]},
  {"6", "Matrix decompositions",
    compactStatus[decompCoreStatuses,
      "some named symbols exist but did not pass the tiny-matrix decomposition smoke"]},
  {"7", "GPU",
    If[gpuArrayOK || targetDeviceOK || TrueQ[hybridProbe[[1]]],
      "PRESENT (caveat: partial GPU support only)",
      "ABSENT"
    ]},
  {"8", "Export", compactStatus[prefixStatuses["8."]]},
  {"9", "ModelFit", compactStatus[prefixStatuses["9."]]}
};

absentHoped =
  Select[results,
    MemberQ[
      {
        "5.MUMPS",
        "7.4",
        "6.LDLDecomposition",
        "6.BunchKaufmanDecomposition",
        "6.OrderedSchurDecomposition",
        "6.PolarDecomposition",
        "6.PopovDecomposition",
        "6.JordanReduce",
        "6.FrobeniusReduce",
        "6.PartialFractions",
        "6.PolynomialReduction"
      },
      #["id"]
    ] && #["status"] === "ABSENT" &
  ];

surprisingPresent =
  Select[results,
    MemberQ[
      {
        "5.Pardiso",
        "5.MUMPS",
        "6.LDLDecomposition",
        "6.BunchKaufmanDecomposition",
        "6.OrderedSchurDecomposition",
        "6.PolarDecomposition"
      },
      #["id"]
    ] && StringStartsQ[#["status"], "PRESENT"] &
  ];

reportLines = Join[
  {
    "# Mathematica 15 Environment Check",
    "",
    "**Generated:** " <> DateString[{"ISODate", " ", "Time", " ", "TimeZoneNameShort"}],
    "",
    "**Working invocation:** `" <> workingInvocation <> "`",
    "",
    "**Kernel:** `$VersionNumber = " <> ToString[$VersionNumber, InputForm] <> "`; `" <>
      md[$Version] <> "`",
    "",
    "**Machine/license:** `$MachineName = " <> md[$MachineName] <>
      "`; `$ProcessorCount = " <> ToString[$ProcessorCount, InputForm] <>
      "`; `$MaxLicenseProcesses = " <> ToString[$MaxLicenseProcesses, InputForm] <> "`",
    "",
    "## Compact Status",
    "",
    "| Item | Route area | Status |",
    "| --- | --- | --- |"
  },
  ("| " <> #[[1]] <> " | " <> md[#[[2]]] <> " | " <> md[#[[3]]] <> " |") & /@ summaryRows,
  {
    "",
    "## Hoped-For Capabilities Requiring Re-Plan",
    ""
  },
  If[Length[absentHoped] == 0,
    {"- None among the explicitly watched hoped-for capabilities. Notably, `Method -> \"MUMPS\"` passed the small sparse smoke probe."},
    ("- `" <> #["probe"] <> "`: " <> #["status"] <> " - " <> md[#["details"]]) & /@ absentHoped
  ],
  {
    "",
    "## Surprising / Useful Present Capabilities",
    ""
  },
  If[Length[surprisingPresent] == 0,
    {"- None beyond the baseline probes."},
    ("- `" <> #["probe"] <> "`: " <> #["status"] <> " - " <> md[#["details"]]) & /@ surprisingPresent
  ],
  {
    "",
    "## Detailed Probe Results",
    "",
    "| ID | Route need | Probe | Status | Details |",
    "| --- | --- | --- | --- | --- |"
  },
  ("| " <> md[#["id"]] <> " | " <> md[#["routeNeed"]] <> " | `" <>
      md[#["probe"]] <> "` | " <> md[#["status"]] <> " | " <> md[#["details"]] <> " |") & /@ results,
  {
    "",
    "## Notes",
    "",
    "- This is an environment capability report, not a performance benchmark.",
    "- Linear algebra `PRESENT` rows mean the method was accepted and returned a correct answer on the stated small smoke problem.",
    "- GPU rows are nice-to-have only; failures there do not block the CPU route.",
    "- The `DefinitelyNotAMethod` LinearSolve row is an internal sanity check that unknown methods are reported as absent.",
    ""
  }
];

Export[reportPath, StringRiffle[reportLines, "\n"], "Text"];

Print["Wrote " <> reportPath];
Exit[0];
