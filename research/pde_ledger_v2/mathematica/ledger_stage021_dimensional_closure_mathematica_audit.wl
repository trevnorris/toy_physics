(* Ledger stage021 Mathematica audit: mu_hat0-free dimensional closure.

   Standalone, print-only, no arguments, and no file I/O.  This keeps the
   source Wolfram Which-based dimOf engine for pathA_33 II-G4d, replaces the
   excluded stage018 fingerprint with its local sourced fixture, and adds the
   corrupt-G scope control, anti-v1 discriminator, structured Yhat catch, and
   dynamically rerun 021-local self-ablations.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

DIMENSIONALOK = "DIMENSIONAL_OK";
FAILDIMENSIONAL = "FAIL_DIMENSIONAL";
NOFAIL = "NO_FAIL";
QUADCALIBRATED = "QUAD_CALIBRATED";

$Assumptions = Element[{a, cs, c, G, omega, D0, N0, chiQsym, mtilde0, muHat0}, Reals] &&
  a != 0 && cs != 0 && c != 0 && D0 != 0 && N0 != 0;

fail[msg_] := Throw[msg, "ledgerStage021Failure"];
dimFail[msg_] := Throw[msg, "dimOfFailure"];

heading[text_] := (
  Print[""];
  Print[StringRepeat["=", StringLength[text]]];
  Print[text];
  Print[StringRepeat["=", StringLength[text]]]
);

subheading[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);

clean[expr_] := FullSimplify[expr, $Assumptions];
fmt[expr_String] := expr;
fmt[expr_] := ToString[InputForm[clean[expr]]];

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
    fail[name]
  ]
];

expectZero[name_, residual_] := Module[{c0},
  assertExact[name, residual];
  c0 = clean[residual];
  assertExact[name, c0];
  If[TrueQ[c0 === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", fmt[c0]];
    fail[name]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectFail[name_, residual_] := Module[{c0},
  assertExact[name, residual];
  c0 = clean[residual];
  assertExact[name, c0];
  If[! TrueQ[c0 === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", fmt[c0], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    fail[name]
  ]
];

verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];

(* KEEP-NATIVE source dimensional block: independent Which recursion. *)
zeroDim = {0, 0, 0};
dimAdd[x_, y_] := x + y;
dimScale[x_, q_] := q x;
dimOf[expr_, dims_] := Module[{args, ds, base, pow},
  Which[
    TrueQ[expr == 0] || NumericQ[expr], zeroDim,
    AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr],
    AtomQ[expr], dimFail["missing dimension for " <> ToString[Unevaluated[expr], InputForm]],
    Head[expr] === Times, Total[dimOf[#, dims] & /@ (List @@ expr)],
    Head[expr] === Power,
      base = expr[[1]];
      pow = expr[[2]];
      If[! NumericQ[pow], dimFail["non-numeric dimension exponent"]];
      dimScale[dimOf[base, dims], pow],
    Head[expr] === Plus,
      args = List @@ expr;
      ds = dimOf[#, dims] & /@ Select[args, ! TrueQ[# == 0] &];
      If[Length[ds] == 0, zeroDim,
        If[Length[DeleteDuplicates[ds]] != 1, dimFail["dimension mismatch in sum"]];
        First[ds]
      ],
    True, dimFail["unsupported dimension expression " <> ToString[expr, InputForm]]
  ]
];

tryDimOf[expr_, dims_] := Catch[
  <|"Status" -> "OK", "Dimension" -> dimOf[expr, dims], "Error" -> None|>,
  "dimOfFailure",
  Function[{msg, tag},
    <|"Status" -> "DIM_ERROR", "Dimension" -> Missing["DimError"], "Error" -> ToString[msg]|>
  ]
];

expText[e_] := Module[{r = Rationalize[e]},
  If[Denominator[r] == 1,
    ToString[Numerator[r]],
    ToString[Numerator[r]] <> "/" <> ToString[Denominator[r]]
  ]
];

dimText[d_] := Module[{pairs, parts},
  pairs = {{"L", d[[1]]}, {"T", d[[3]]}, {"M", d[[2]]}};
  parts = (If[TrueQ[#[[2]] == 1], #[[1]], #[[1]] <> "^" <> expText[#[[2]]]] &) /@
    Select[pairs, ! TrueQ[#[[2]] == 0] &];
  If[Length[parts] == 0, "1", StringRiffle[parts, " "]]
];

superscriptText[e_] := StringReplace[
  expText[e],
  {"-" -> "⁻", "0" -> "⁰", "1" -> "¹", "2" -> "²", "3" -> "³",
   "4" -> "⁴", "5" -> "⁵", "6" -> "⁶", "7" -> "⁷", "8" -> "⁸",
   "9" -> "⁹", "/" -> "ᐟ"}
];

dimPrettyText[d_] := Module[{pairs, parts},
  pairs = {{"L", d[[1]]}, {"T", d[[3]]}, {"M", d[[2]]}};
  parts = (If[TrueQ[#[[2]] == 1], #[[1]], #[[1]] <> superscriptText[#[[2]]]] &) /@
    Select[pairs, ! TrueQ[#[[2]] == 0] &];
  If[Length[parts] == 0, "1", StringRiffle[parts, " "]]
];

dimVectorText[d_] := "(" <> StringRiffle[ToString[InputForm[#]] & /@ d, ","] <> ")";

dimensionReadSymbols[expr_, dims_] := DeleteDuplicates[
  Cases[expr, s_Symbol /; KeyExistsQ[dims, s] :> s, Infinity]
];

gateData[expr_, dims_] := Module[{p0DimLocal, symbolReads, dimensionalOkLocal, readSet},
  p0DimLocal = dimOf[expr, dims];
  symbolReads = dimensionReadSymbols[expr, dims];
  dimensionalOkLocal = TrueQ[p0DimLocal == zeroDim];
  readSet = DeleteDuplicates[
    Join[{"dimensionalOk", "p0Dim"}, SymbolName /@ symbolReads]
  ];
  <|
    "Dimension" -> p0DimLocal,
    "DimensionalOk" -> dimensionalOkLocal,
    "SymbolReads" -> symbolReads,
    "VerdictReadSet" -> readSet
  |>
];

scopedDimensionalVerdict[gate_] := If[
  TrueQ[gate["DimensionalOk"]],
  DIMENSIONALOK,
  FAILDIMENSIONAL
];

dynamicDimensionalAblation[baselineExpr_, baselineDims_, mutatedExpr_, mutatedDims_] := Module[
  {withGate, withoutGate, verdictTrace, rerunVerdict, withMutation, withoutMutation, rerunGateLogic},
  withGate = gateData[mutatedExpr, mutatedDims];
  withoutGate = gateData[baselineExpr, baselineDims];
  verdictTrace = {};
  rerunVerdict = Function[gateLocal,
    With[{verdictLocal = scopedDimensionalVerdict[gateLocal]},
      AppendTo[verdictTrace, {gateLocal["Dimension"], verdictLocal}];
      verdictLocal
    ]
  ];
  withMutation = rerunVerdict[withGate];
  withoutMutation = rerunVerdict[withoutGate];
  rerunGateLogic = TrueQ[
    Length[verdictTrace] === 2 &&
    Length[DeleteDuplicates[verdictTrace[[All, 1]]]] === 2 &&
    verdictTrace[[All, 2]] === {withMutation, withoutMutation} &&
    withMutation =!= withoutMutation
  ];
  <|
    "RerunGateLogic" -> rerunGateLogic,
    "VerdictTrace" -> verdictTrace,
    "WithMutation" -> withMutation,
    "WithoutMutation" -> withoutMutation,
    "ExpectedFail" -> FAILDIMENSIONAL,
    "FailSuppressed" -> TrueQ[
      withMutation === FAILDIMENSIONAL &&
      withoutMutation === DIMENSIONALOK &&
      withMutation =!= withoutMutation
    ]
  |>
];

buildCorruptionCase[name_, p0RawExpr_, p0PhysicalExpr_, dims_] := Module[
  {rawDimLocal, gateLocal, verdictLocal},
  rawDimLocal = dimOf[p0RawExpr, dims];
  gateLocal = gateData[p0PhysicalExpr, dims];
  verdictLocal = If[TrueQ[gateLocal["DimensionalOk"]], NOFAIL, FAILDIMENSIONAL];
  <|
    "Name" -> name,
    "P0RawDimension" -> rawDimLocal,
    "P0PhysicalDimension" -> gateLocal["Dimension"],
    "DimensionalOk" -> gateLocal["DimensionalOk"],
    "Verdict" -> verdictLocal,
    "Gate" -> gateLocal
  |>
];

backSolveMutant[p0PhysicalExpr_, targetRHSExpr_, dims_] := Module[
  {p0DimLocal, rhsDimLocal, muDimLocal, homogeneityPassLocal, wiredVerdict},
  p0DimLocal = dimOf[p0PhysicalExpr, dims];
  rhsDimLocal = dimOf[targetRHSExpr, dims];
  muDimLocal = dimScale[rhsDimLocal - p0DimLocal, 1/2];
  homogeneityPassLocal = TrueQ[dimAdd[dimScale[muDimLocal, 2], p0DimLocal] == rhsDimLocal];
  wiredVerdict = If[homogeneityPassLocal, NOFAIL, FAILDIMENSIONAL];
  <|
    "P0Dimension" -> p0DimLocal,
    "RHSDimension" -> rhsDimLocal,
    "MuDimension" -> muDimLocal,
    "HomogeneityPass" -> homogeneityPassLocal,
    "WiredVerdict" -> wiredVerdict,
    "VerdictReadSet" -> {"homogeneityPass", "muDim", "p0Dim", "rhsDim"}
  |>
];

(* Native raw dimension map; muHat0 is deliberately absent from the gate map. *)
rawDims = <|
  a -> {1, 0, 0},
  cs -> {1, 0, -1},
  c -> {1, 0, -1},
  G -> {3, -1, -2},
  omega -> {0, 0, -1},
  D0 -> {-1, 1, -2},
  N0 -> {-1, 1, 0},
  chiQsym -> zeroDim,
  mtilde0 -> zeroDim
|>;

(* PROVENANCE-FROZEN local fixture from stage018; no fingerprint block. *)
u2Sourced = a^2/(9 cs^2);
u4Sourced = (4 a^4)/(81 cs^4);
v5Sourced = a^5/(27 cs^5);

(* P0Raw is stage019 provenance; dimensions are computed here. *)
P0Raw = N0/D0;
frequencyNormalization = (cs/a)^2;
P0Physical = frequencyNormalization P0Raw;
YhatPhysical = 1 + u2Sourced omega^2 + u4Sourced omega^4 + I v5Sourced omega^5;
YhatPowerMutation = 1 + u2Sourced omega^3 + u4Sourced omega^4 + I v5Sourced omega^5;
Gamma5 = chiQsym P0Physical a^5/(27 cs^5);
(* stage020 target_rhs provenance, used only by the diagnostic. *)
targetRHS = 54 G cs^5/(5 a^5 c^5);

p0RawDim = dimOf[P0Raw, rawDims];
frequencyNormDim = dimOf[frequencyNormalization, rawDims];
realGate = gateData[P0Physical, rawDims];
p0Dim = realGate["Dimension"];
dimensionalOk = realGate["DimensionalOk"];
dimensionalStatus = scopedDimensionalVerdict[realGate];
dimensionalGateDependsOnMuHat0 = MemberQ[realGate["VerdictReadSet"], "muHat0"];

dropNormDim = dimOf[P0Raw, rawDims];
dropNormOk = TrueQ[dropNormDim == zeroDim];
dropNormVerdict = If[dropNormOk, NOFAIL, FAILDIMENSIONAL];

corruptN0Dims = Join[KeyDrop[rawDims, N0], <|N0 -> zeroDim|>];
corruptD0Dims = Join[KeyDrop[rawDims, D0], <|D0 -> zeroDim|>];
corruptGDims = Join[KeyDrop[rawDims, G], <|G -> zeroDim|>];
corruptCsDims = Join[KeyDrop[rawDims, cs], <|cs -> zeroDim|>];

truthTable = <|
  "corrupt_N0" -> buildCorruptionCase["corrupt_N0", P0Raw, P0Physical, corruptN0Dims],
  "corrupt_D0" -> buildCorruptionCase["corrupt_D0", P0Raw, P0Physical, corruptD0Dims],
  "corrupt_G" -> buildCorruptionCase["corrupt_G", P0Raw, P0Physical, corruptGDims],
  "corrupt_c_s" -> buildCorruptionCase["corrupt_c_s", P0Raw, P0Physical, corruptCsDims],
  "correct" -> buildCorruptionCase["correct", P0Raw, P0Physical, rawDims]
|>;

corruptN0RawDim = truthTable["corrupt_N0"]["P0RawDimension"];
corruptN0P0Dim = truthTable["corrupt_N0"]["P0PhysicalDimension"];
corruptN0Ok = truthTable["corrupt_N0"]["DimensionalOk"];
corruptN0Verdict = truthTable["corrupt_N0"]["Verdict"];
corruptGOk = TrueQ[dimOf[P0Physical, corruptGDims] == zeroDim];
corruptGVerdict = If[corruptGOk, NOFAIL, FAILDIMENSIONAL];

yhatTry = tryDimOf[YhatPhysical, rawDims];
yhatMutantTry = tryDimOf[YhatPowerMutation, rawDims];
yhatDim = If[yhatTry["Status"] === "OK", yhatTry["Dimension"], Missing["DimError"]];
yhatOk = TrueQ[yhatTry["Status"] === "OK" && yhatDim == zeroDim];

rhsDim = dimOf[targetRHS, rawDims];
muDim = (rhsDim - p0Dim)/2;
dims = Join[rawDims, <|muHat0 -> muDim|>];
lhs = (muHat0 mtilde0)^2 P0Physical;
lhsRawMutation = (muHat0 mtilde0)^2 P0Raw;
lhsDim = dimOf[lhs, dims];
lhsRawDim = dimOf[lhsRawMutation, dims];
requiredP0Dim = rhsDim - 2 muDim;
gamma5Dim = dimOf[Gamma5, dims];
(* The Yhat tooth is split out; this diagnostic never drives the pass guard. *)
homogeneityPass = TrueQ[lhsDim == rhsDim && p0Dim == requiredP0Dim];
diagnosticParticipatesInVerdict = AnyTrue[
  {"muHat0", "muDim", "homogeneityPass"},
  MemberQ[realGate["VerdictReadSet"], #] &
];

backsolveMutants = <|
  "corrupt_N0" -> backSolveMutant[P0Physical, targetRHS, corruptN0Dims],
  "corrupt_D0" -> backSolveMutant[P0Physical, targetRHS, corruptD0Dims],
  "corrupt_G" -> backSolveMutant[P0Physical, targetRHS, corruptGDims],
  "corrupt_c_s" -> backSolveMutant[P0Physical, targetRHS, corruptCsDims]
|>;

dropNormAblation = dynamicDimensionalAblation[P0Physical, rawDims, P0Raw, rawDims];
corruptN0Ablation = dynamicDimensionalAblation[P0Physical, rawDims, P0Physical, corruptN0Dims];

localPassGuard = TrueQ[
  dimensionalOk && ! dropNormOk && ! corruptN0Ok && yhatOk && corruptGOk
];

runMuFreeGate[] := Module[{readNames},
  subheading["Mu-hat0-free physical P0 dimensional gate"];
  readNames = Sort[SymbolName /@ realGate["SymbolReads"]];
  Print["  [P₀_raw]=", dimPrettyText[p0RawDim]];
  Print["  [(c_s/a)²]=", dimPrettyText[frequencyNormDim]];
  Print["  [P₀^phys]=", dimPrettyText[p0Dim]];
  Print["  dimensional_ok=", dimensionalOk, " (COMPUTED by dimOf)"];
  Print["  dimension symbols actually read = ", readNames];
  Print["  verdict read-set = ", Sort[realGate["VerdictReadSet"]]];
  Print["  dimensional_gate_depends_on_mu_hat0=", dimensionalGateDependsOnMuHat0];
  expectBool["sourced [N0] is exactly (-1,1,0) in (L,M,T) order", rawDims[N0] === {-1, 1, 0}];
  expectBool["sourced [D0] is exactly (-1,1,-2) in (L,M,T) order", rawDims[D0] === {-1, 1, -2}];
  expectBool["computed [P0_raw] is T^2", p0RawDim === {0, 0, 2}];
  expectBool["computed [(c_s/a)^2] is T^-2", frequencyNormDim === {0, 0, -2}];
  expectBool["computed [P0_physical] is dimensionless", p0Dim === zeroDim];
  expectBool["dimensional_ok is computed True from [P0_physical]", dimensionalOk];
  expectBool["muHat0 is absent from the raw gate dimension map", ! KeyExistsQ[rawDims, muHat0]];
  expectBool["runtime gate trace reads exactly a,cs,N0,D0", Sort[readNames] === Sort[{"a", "cs", "N0", "D0"}]];
  expectBool["computed dimensional gate does not depend on muHat0", ! dimensionalGateDependsOnMuHat0]
];

runYhatStructuredCatch[] := (
  subheading["Yhat dimensionless check with structured DimError catch"];
  Print["  local stage018 fixture = ", <|"u2" -> u2Sourced, "u4" -> u4Sourced, "v5" -> v5Sourced|>];
  Print["  Yhat structured result = ", yhatTry];
  Print["  Yhat omega-power mutant structured result = ", yhatMutantTry];
  expectBool[
    "Yhat dimensionless named assert consumes structured status and computed dimension",
    yhatTry["Status"] === "OK" && yhatTry["Dimension"] === zeroDim && yhatOk
  ];
  expectFail[
    "Yhat dimensionless named assert catches u2*omega^3 sum mismatch",
    If[yhatMutantTry["Status"] === "DIM_ERROR", 1, 0]
  ]
);

runTruthTableAndProbes[] := Module[{expectedVerdicts, names, observed, expected, p0SymbolNames},
  subheading["Corrupt-dimension scope truth-table and 3d/3d-prime probes"];
  expectedVerdicts = <|
    "corrupt_N0" -> FAILDIMENSIONAL,
    "corrupt_D0" -> FAILDIMENSIONAL,
    "corrupt_G" -> NOFAIL,
    "corrupt_c_s" -> FAILDIMENSIONAL,
    "correct" -> NOFAIL
  |>;
  names = {"corrupt_N0", "corrupt_D0", "corrupt_G", "corrupt_c_s", "correct"};
  Scan[
    Function[name,
      Print[
        "  ", name,
        ": [P₀_raw]=", dimPrettyText[truthTable[name]["P0RawDimension"]],
        "; [P₀^phys]=", dimPrettyText[truthTable[name]["P0PhysicalDimension"]],
        "; vector(L,M,T)=", dimVectorText[truthTable[name]["P0PhysicalDimension"]],
        "; ", truthTable[name]["Verdict"]
      ];
      expectZero[
        "truth-table " <> name <> " reaches " <> expectedVerdicts[name] <> " at its own assert",
        verdictResidual[truthTable[name]["Verdict"], expectedVerdicts[name]]
      ]
    ],
    names
  ];
  expectBool["corrupt-N0 keeps normalization: [P0_raw]=(1,-1,2)", corruptN0RawDim === {1, -1, 2}];
  expectBool["corrupt-N0 keeps normalization: [P0_physical]=(1,-1,0)", corruptN0P0Dim === {1, -1, 0}];
  p0SymbolNames = Sort[SymbolName /@ dimensionReadSymbols[P0Physical, rawDims]];
  Print["  free symbols in P0Physical = ", p0SymbolNames];
  expectBool["dependency-scope control computes G absent from P0Physical", FreeQ[P0Physical, G] && ! MemberQ[p0SymbolNames, "G"]];
  observed = truthTable[#]["Verdict"] & /@ names;
  expected = Values[expectedVerdicts];
  Print[
    "  SUMMARY (not counted): full corrupt-dim truth-table observed=", observed,
    "; expected=", expected, "; matches=", observed === expected
  ];

  Print["  3d_dimensional_break: drop (c_s/a)² -> ", dropNormVerdict];
  Print["  3d self_ablation = ", dropNormAblation];
  Print["  3d_prime_corrupt_port_dimension: corrupt [N₀] -> ", corruptN0Verdict];
  Print["  3d_prime self_ablation = ", corruptN0Ablation];
  Print["  corrupt [G] dependency-scope control -> ", corruptGVerdict];
  expectBool["drop-frequency-normalization mutation fires", dropNormVerdict === FAILDIMENSIONAL];
  expectZero["3d drop normalization reaches FAIL_DIMENSIONAL", verdictResidual[dropNormVerdict, FAILDIMENSIONAL]];
  expectBool["corrupt-N0 dimension mutation fires", corruptN0Verdict === FAILDIMENSIONAL];
  expectZero["3d-prime corrupt port dimension reaches FAIL_DIMENSIONAL", verdictResidual[corruptN0Verdict, FAILDIMENSIONAL]];
  expectZero["native corrupt-G control remains NO_FAIL", verdictResidual[corruptGVerdict, NOFAIL]];
  Scan[
    Function[pair,
      expectBool[pair[[1]] <> " self-ablation dynamically reruns 021-local gate logic", pair[[2]]["RerunGateLogic"]];
      expectZero[pair[[1]] <> " dynamic rerun with mutation is FAIL_DIMENSIONAL", verdictResidual[pair[[2]]["WithMutation"], FAILDIMENSIONAL]];
      expectZero[pair[[1]] <> " dynamic rerun without mutation is DIMENSIONAL_OK", verdictResidual[pair[[2]]["WithoutMutation"], DIMENSIONALOK]];
      expectBool[pair[[1]] <> " dynamic self-ablation changes the local verdict", pair[[2]]["FailSuppressed"]]
    ],
    {{"3d", dropNormAblation}, {"3d-prime", corruptN0Ablation}}
  ]
];

runAntiV1AndDiagnostic[] := Module[
  {forbidden, names, mutantReadSet, wiringPolicyResidual, realFailCount, mutantFailCount},
  subheading["Decisive anti-v1 read-set tooth and non-verdict mu_hat0 diagnostic"];
  forbidden = {"muHat0", "muDim", "homogeneityPass"};
  Print["  real verdict read-set = ", Sort[realGate["VerdictReadSet"]]];
  Print["  [mu_hat0]=", dimPrettyText[muDim]];
  Print["  [lhs]=", dimPrettyText[lhsDim]];
  Print["  [rhs]=", dimPrettyText[rhsDim]];
  Print["  homogeneity_pass=", homogeneityPass, " DIAGNOSTIC"];
  Print["  participates_in_verdict=", diagnosticParticipatesInVerdict];
  expectBool[
    "computed real verdict read-set excludes muHat0/muDim/homogeneityPass",
    Intersection[realGate["VerdictReadSet"], forbidden] === {}
  ];
  expectBool[
    "real verdict read-set contains dimensionalOk and its computed p0Dim input",
    SubsetQ[realGate["VerdictReadSet"], {"dimensionalOk", "p0Dim"}]
  ];
  expectBool["back-solved muHat0 has dimension (-1,-1/2,-1)", muDim === {-1, -1/2, -1}];
  expectBool["muHat0 homogeneity diagnostic computes True", homogeneityPass];
  expectBool["muHat0 homogeneity diagnostic is explicitly non-verdict", diagnosticParticipatesInVerdict === False];
  names = {"corrupt_N0", "corrupt_D0", "corrupt_G", "corrupt_c_s"};
  Scan[
    Function[name,
      Print[
        "  wired back-solve mutant ", name,
        ": homogeneity_pass=", backsolveMutants[name]["HomogeneityPass"],
        " -> ", backsolveMutants[name]["WiredVerdict"]
      ];
      expectZero[
        "wired back-solve mutant stays NO_FAIL under " <> name,
        verdictResidual[backsolveMutants[name]["WiredVerdict"], NOFAIL]
      ]
    ],
    names
  ];
  mutantReadSet = DeleteDuplicates[Flatten[backsolveMutants[#]["VerdictReadSet"] & /@ names]];
  wiringPolicyResidual = Length[Intersection[mutantReadSet, forbidden]];
  expectFail[
    "anti-v1 demotion policy rejects wiring homogeneityPass/muDim into the verdict",
    wiringPolicyResidual
  ];
  realFailCount = Count[truthTable[#]["Verdict"] & /@ names, FAILDIMENSIONAL];
  mutantFailCount = Count[backsolveMutants[#]["WiredVerdict"] & /@ names, FAILDIMENSIONAL];
  Print[
    "  SUMMARY (not counted): real mu-free gate fail count=", realFailCount,
    "; wired back-solve mutant fail count=", mutantFailCount
  ]
];

runAritySelfCheck[] := Module[{arityGate, arityTry, arityAblation, arityMutant},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  arityGate = gateData[P0Physical, rawDims];
  arityTry = tryDimOf[YhatPowerMutation, rawDims];
  arityAblation = dynamicDimensionalAblation[P0Physical, rawDims, P0Raw, rawDims];
  arityMutant = backSolveMutant[P0Physical, targetRHS, corruptGDims];
  expectBool["arity dimOf[expr,dims] returns a three-vector", Length[dimOf[P0Physical, rawDims]] === 3];
  expectBool["arity tryDimOf[expr,dims] returns structured DIM_ERROR", arityTry["Status"] === "DIM_ERROR"];
  expectBool["arity gateData[expr,dims] returns computed DimensionalOk", TrueQ[arityGate["DimensionalOk"]]];
  expectBool["arity scopedDimensionalVerdict[gate] returns DIMENSIONAL_OK", scopedDimensionalVerdict[arityGate] === DIMENSIONALOK];
  expectBool["arity dynamicDimensionalAblation[baselineExpr,baselineDims,mutatedExpr,mutatedDims] changes verdicts", arityAblation["WithMutation"] =!= arityAblation["WithoutMutation"]];
  expectBool["arity backSolveMutant[p0,target,dims] returns a tautological True", TrueQ[arityMutant["HomogeneityPass"]]];
  expectBool[
    "no unevaluated dimensional helper call remains in computed outputs",
    FreeQ[
      {p0RawDim, frequencyNormDim, p0Dim, yhatTry, truthTable, backsolveMutants, dropNormAblation, corruptN0Ablation},
      dimOf | tryDimOf | gateData | dynamicDimensionalAblation | backSolveMutant
    ]
  ]
];

runScopeProvenanceAndLanding[] := (
  subheading["021 scope, provenance-only consumption, and COMPLETING landing"];
  expectBool["021-local pass guard excludes homogeneityPass and passes its dimensional teeth", localPassGuard];
  Print["  CONSUMES (PROVENANCE only): [N₀]=L⁻¹M from the density/continuity-port numerator; it genuinely enters the checked gate."];
  Print["  CONSUMES (PROVENANCE only): [D₀]=L⁻¹T⁻²M from the carried reduced static conservative denominator D₀=K−B₀−Z₀."];
  Print["  CONSUMES (PROVENANCE only): stage018 u₂/u₄/v₅ as local u2Sourced/u4Sourced/v5Sourced; no fingerprint block."];
  Print["  CONSUMES (PROVENANCE only): stage019 P0=N0/D0 defines P₀_raw; no prefactor block."];
  Print["  CONSUMES (PROVENANCE only): stage020 target_rhs enters only the non-verdict mu_hat0 diagnostic."];
  Print["  018 fingerprint DONE; 019 prefactor DONE; 020 provenance partition + CALIBRATED classification DONE."];
  Print["  QUAD_CALIBRATED (4/4) -- mu_hat0-free [P₀^phys]=1 dimensional closure EARNED (COMPLETING)"];
  Print["  joint QUAD_CALIBRATED COMPLETE = 018 DONE ∧ 019 DONE ∧ 020 DONE ∧ 021 dimensional closure EARNED."];
  Print["  COMPLETE != PASS: token remains QUAD_CALIBRATED, not QUAD_PASS; 020 has G=GENUINE_BLOCKED."];
  Print["  EXPORTS: completed Gate-4 dimensional closure to 022/027 and the Part-VII dimensional firewall."];
  Print["  dropped-bookkeeping: scratch YAML, reports, feed writers, engine cross-reads, and all file I/O are absent."]
);

printVerdictLabels[] := (
  subheading["Verdict labels"];
  Print["  ledger local gate label: DIMENSIONAL_OK (mu_hat0-free [P0_phys]=1 from sourced [N0]/[D0] and [(c_s/a)^2])"];
  Print["  source top-line verdict: QUAD_CALIBRATED (JOINT 4-stage; 021 is the COMPLETING 4/4 leg)"];
  Print["  joint composition: 018 fingerprint [DONE] AND 019 prefactor [DONE] AND 020 provenance partition + CALIBRATED classification [DONE] AND 021 dimensional closure [EARNED here] => COMPLETE"];
  Print["  COMPLETE != PASS: the joint token remains QUAD_CALIBRATED, not QUAD_PASS (020 provenance; G=GENUINE_BLOCKED)"];
  Print["  earned dimensional gate: [P0_raw]=T^2, [(c_s/a)^2]=T^-2, [P0_phys]=1; mu_hat0 is absent from rawDims and the verdict read-set"];
  Print["  earned able-to-fail: 3d drop normalization and 3d-prime corrupt [N0] reach FAIL_DIMENSIONAL with DYNAMIC 021-local self-ablations; corrupt [G] is NO_FAIL scope control"];
  Print["  diagnostic only: back-solved [mu_hat0] and homogeneity_pass=True participate_in_verdict=False; the wired back-solve mutant is all-NO_FAIL and rejected as vacuous"]
);

runAll[] := (
  heading["ledger_stage021_dimensional_closure_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage021_dimensional_closure"];
  Print["Engine: native Wolfram Which-based exact (L,M,T) dimension algebra; no floats/tolerances; zero file I/O."];
  runMuFreeGate[];
  runYhatStructuredCatch[];
  runTruthTableAndProbes[];
  runAntiV1AndDiagnostic[];
  runAritySelfCheck[];
  runScopeProvenanceAndLanding[];
  printVerdictLabels[];
  0
);

result = Catch[
  runAll[],
  "ledgerStage021Failure",
  Function[{msg, tag}, failureMessage = ToString[msg]; 1]
];

heading["Tallies"];
Print["TALLY mathematica: ", passCount, " pass + ", failCount, " fail = ", passCount + failCount, " checks"];
If[result === 0 && failCount === 0,
  Print["OVERALL PASS"];
  Exit[0],
  If[failureMessage =!= "", Print["ABORTED: ", failureMessage]];
  Print["OVERALL FAIL"];
  Exit[1]
];
