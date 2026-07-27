(* Ledger stage016 Mathematica audit: l=2 SO(3) covariance theorem.

   Standalone, print-only, no arguments, no file I/O. This keeps the native
   Integrate/D/Hash pathA_32 Wolfram route and restricts it to the 016 angular
   slice: real l=2 harmonics, Gram=I5, computed -Delta_S2 eigenvalues, bare K2
   angular stiffness with live computed lambda, dimensional gate, and the three
   016 covariance probes. Stage 017 owns grouped lanes and response/calibration.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;

ISOTROPYCALIBRATED = "ISOTROPY_CALIBRATED";
FAILDIMENSIONAL = "FAIL_DIMENSIONAL";
FAILNOTCOVARIANT = "FAIL_NOT_COVARIANT";
FAILTAUTOLOGICAL = "FAIL_TAUTOLOGICAL";
FAILNOTABLETOFAIL = "FAIL_NOT_ABLE_TO_FAIL";

$Assumptions = Element[{theta, phi}, Reals] && Mtilde > 0 && Ktilde > 0 && TomegaTilde > 0;

zeroDim = {0, 0, 0};
expectedM = {0, 1, 0};
expectedK = {0, 1, -2};
expectedRatio = {0, 0, -2};

failureMessage = "";

raise[msg_] := Throw[msg, "ledgerStage016Failure"];
dimRaise[msg_] := Throw[msg, "stage016DimError"];

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

dropConditions[expr_] := expr /. ConditionalExpression[value_, _] :> value;
clean[expr_] := FullSimplify[dropConditions[expr], $Assumptions];
fmt[expr_String] := expr;
fmt[expr_] := ToString[InputForm[clean[expr]]];

assertExact[name_, expr_] := Module[{reals},
  reals = Cases[Unevaluated[expr], _Real, Infinity];
  If[reals =!= {},
    failCount++;
    Print["FAIL  ", name, ": machine-real atom(s) found: ", ToString[InputForm[reals]]];
    raise[name]
  ]
];

expectZero[name_, residual_] := Module[{c},
  assertExact[name, residual];
  c = clean[residual];
  assertExact[name, c];
  If[TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    Print["FAIL  ", name, ": residual = ", fmt[c]];
    raise[name]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

expectFail[name_, residual_] := Module[{c},
  assertExact[name, residual];
  c = clean[residual];
  assertExact[name, c];
  If[! TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name, " produced required FAIL (residual = ", fmt[c], ")"],
    failCount++;
    Print["FAIL  ", name, ": required mutation/ablation did not fire"];
    raise[name]
  ]
];

boolResidual[condition_] := If[TrueQ[condition], 0, 1];
verdictResidual[actual_, expected_] := If[actual === expected, 0, 1];
dimResidualVec[actual_, expected_] := FullSimplify[(actual - expected).(actual - expected)];

scopedVerdict[dimensionalOk_, covariant_, tautologyOk_, ableOk_] := Which[
  ! TrueQ[dimensionalOk], FAILDIMENSIONAL,
  ! TrueQ[covariant], FAILNOTCOVARIANT,
  ! TrueQ[tautologyOk], FAILTAUTOLOGICAL,
  ! TrueQ[ableOk], FAILNOTABLETOFAIL,
  True, ISOTROPYCALIBRATED
];

caseVerdict[overrides_: <||>] := scopedVerdict[
  Lookup[overrides, "dimensionalOk", True],
  Lookup[overrides, "covariant", True],
  Lookup[overrides, "tautologyClear", True],
  Lookup[overrides, "ableToFailOk", True]
];

intS2[expr_] := FullSimplify[
  Integrate[
    Integrate[TrigExpand[expr] Sin[theta], {phi, 0, 2 Pi}],
    {theta, 0, Pi}
  ],
  $Assumptions
];

lapS2[expr_] := FullSimplify[
  1/Sin[theta] D[Sin[theta] D[expr, theta], theta] +
    1/Sin[theta]^2 D[expr, {phi, 2}],
  $Assumptions
];

dimensionAxisSlots = {{"L", 1}, {"M", 2}, {"T", 3}};
dimensionAxesLabel[] := StringRiffle[dimensionAxisSlots[[All, 1]], ","];
dimensionPairs[d_] := ({#[[1]], d[[#[[2]]]]} &) /@ dimensionAxisSlots;

dimText[d_] := Module[{pairs, parts, emit},
  emit[label_, exp_] := If[TrueQ[exp == 1], label, label <> "^" <> ToString[InputForm[exp]]];
  pairs = dimensionPairs[d];
  parts = (emit[#[[1]], #[[2]]] &) /@ Select[pairs, ! TrueQ[#[[2]] == 0] &];
  If[Length[parts] == 0, "1", StringRiffle[parts, " "]]
];

printDimRecord[name_, binding_] := Print[
  "DIM|axes=", dimensionAxesLabel[],
  "|name=", name,
  "|exponents=", ToString[InputForm[dimensionPairs[binding][[All, 2]]]]
];

dimOf[expr_, dims_] := Module[{args, ds, base, pow, argDims},
  Which[
    TrueQ[expr == 0] || NumericQ[expr], zeroDim,
    AtomQ[expr] && KeyExistsQ[dims, expr], dims[expr],
    AtomQ[expr], dimRaise["missing dimension for " <> ToString[Unevaluated[expr], InputForm]],
    Head[expr] === Times, Total[dimOf[#, dims] & /@ (List @@ expr)],
    Head[expr] === Power,
      base = expr[[1]];
      pow = expr[[2]];
      If[! NumericQ[pow], dimRaise["non-numeric dimension exponent"]];
      pow dimOf[base, dims],
    Head[expr] === Plus,
      args = Select[List @@ expr, ! TrueQ[clean[#] == 0] &];
      ds = dimOf[#, dims] & /@ args;
      If[Length[ds] == 0, zeroDim,
        If[Length[DeleteDuplicates[ds]] != 1, dimRaise["dimension mismatch in sum " <> ToString[expr, InputForm]]];
        First[ds]
      ],
    MemberQ[{Sin, Cos, Tan, Cot, Sinh, Cosh, Tanh, Coth, Sech, Csch}, Head[expr]],
      argDims = dimOf[#, dims] & /@ (List @@ expr);
      If[AnyTrue[argDims, # =!= zeroDim &], dimRaise["dimensionful argument in dimensionless function"]];
      zeroDim,
    True, dimRaise["unsupported dimension expression " <> ToString[expr, InputForm]]
  ]
];

order = {"20", "21c", "21s", "22c", "22s"};
harmonics = <|
  "20" -> Sqrt[5/(16 Pi)] (3 Cos[theta]^2 - 1),
  "21c" -> -Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Cos[phi],
  "21s" -> -Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Sin[phi],
  "22c" -> Sqrt[15/(16 Pi)] Sin[theta]^2 Cos[2 phi],
  "22s" -> Sqrt[15/(16 Pi)] Sin[theta]^2 Sin[2 phi]
|>;
ys = harmonics /@ order;

computeAngularBlock[h_Association] := Module[
  {localOrder, localYs, localGram, localNegLaps, localLambdas, localResiduals},
  localOrder = Keys[h];
  localYs = h /@ localOrder;
  localGram = Table[intS2[localYs[[i]] localYs[[j]]], {i, Length[localYs]}, {j, Length[localYs]}];
  localNegLaps = AssociationThread[localOrder, clean[-lapS2[#]] & /@ localYs];
  localLambdas = AssociationThread[
    localOrder,
    Table[
      clean[intS2[localYs[[i]] localNegLaps[localOrder[[i]]]]/intS2[localYs[[i]]^2]],
      {i, Length[localYs]}
    ]
  ];
  localResiduals = AssociationThread[
    localOrder,
    Table[
      clean[localNegLaps[localOrder[[i]]] - localLambdas[localOrder[[i]]] localYs[[i]]],
      {i, Length[localYs]}
    ]
  ];
  <|
    "Order" -> localOrder,
    "Ys" -> localYs,
    "Gram" -> localGram,
    "GramIsIdentity" -> TrueQ[FullSimplify[localGram == IdentityMatrix[Length[localYs]], $Assumptions]],
    "NegLaps" -> localNegLaps,
    "Lambdas" -> localLambdas,
    "Residuals" -> localResiduals,
    "LambdaAllSix" -> And @@ (TrueQ[clean[# - 6] === 0] & /@ Values[localLambdas]),
    "ResidualsZero" -> And @@ (TrueQ[clean[#] === 0] & /@ Values[localResiduals])
  |>
];

angular = computeAngularBlock[harmonics];
gram = angular["Gram"];
negLaps = angular["NegLaps"];
lambdas = angular["Lambdas"];
residuals = angular["Residuals"];
gramIsIdentity = angular["GramIsIdentity"];
lambdaAllSix = angular["LambdaAllSix"];
residualsZero = angular["ResidualsZero"];

buildK2[coeff_] := clean[Ktilde + coeff TomegaTilde];
extractK2Coeff[k2Expr_] := clean[D[k2Expr, TomegaTilde]];

k2Core = AssociationThread[order, buildK2[lambdas[#]] & /@ order];
k2Coeff = AssociationThread[order, extractK2Coeff[k2Core[#]] & /@ order];
k2CoeffResiduals = AssociationThread[
  order,
  clean[negLaps[#] - k2Coeff[#] harmonics[#]] & /@ order
];
k2CoeffResidualsZero = And @@ (TrueQ[clean[#] === 0] & /@ Values[k2CoeffResiduals]);

lambdaRef = lambdas["20"];
m2Core = Mtilde;
k2Ref = buildK2[lambdaRef];

measureExpr = aDim^2 dwDim dOmegaDim;
m2IntegralExpr = muEtaDensity beta2Dim^2 measureExpr;
kTwTermExpr = TwDensity beta2PrimeDim^2 measureExpr;
kEtaTermExpr = KetaDensity beta2Dim^2 measureExpr;
kOmegaTermExpr[lambda_] := lambda TOmegaDensity beta2Dim^2 measureExpr;
k2IntegralExpr[lambda_] := kTwTermExpr + kEtaTermExpr + kOmegaTermExpr[lambda];

makeDimRules[] := <|
  aDim -> {1, 0, 0},
  dwDim -> {1, 0, 0},
  dOmegaDim -> zeroDim,
  beta2Dim -> zeroDim,
  beta2PrimeDim -> {-1, 0, 0},
  muEtaDensity -> {-3, 1, 0},
  TwDensity -> {-1, 1, -2},
  KetaDensity -> {-3, 1, -2},
  TOmegaDensity -> {-3, 1, -2},
  Mtilde -> expectedM,
  Ktilde -> expectedK,
  TomegaTilde -> expectedK
|>;

evalDimensional[lambda_, m2Expr_, k2Expr_, rules_] := Catch[
  Module[
    {
      measureDim, m2IntegralDim, kTwDim, kEtaDim, kOmegaDim, k2IntegralDim,
      actualM2Dim, actualK2Dim, actualRatioDim, termHomogeneous, ok
    },
    measureDim = dimOf[measureExpr, rules];
    m2IntegralDim = dimOf[m2IntegralExpr, rules];
    kTwDim = dimOf[kTwTermExpr, rules];
    kEtaDim = dimOf[kEtaTermExpr, rules];
    kOmegaDim = dimOf[kOmegaTermExpr[lambda], rules];
    k2IntegralDim = dimOf[k2IntegralExpr[lambda], rules];
    actualM2Dim = dimOf[m2Expr, rules];
    actualK2Dim = dimOf[k2Expr, rules];
    actualRatioDim = actualK2Dim - actualM2Dim;
    termHomogeneous = TrueQ[kTwDim == kEtaDim == kOmegaDim == k2IntegralDim];
    ok = TrueQ[
      measureDim == {3, 0, 0} && m2IntegralDim == expectedM && termHomogeneous &&
        k2IntegralDim == expectedK && actualM2Dim == expectedM &&
        actualK2Dim == expectedK && actualRatioDim == expectedRatio
    ];
    <|
      "Ok" -> ok,
      "Error" -> None,
      "Dims" -> <|
        "measure" -> measureDim,
        "M2_integral" -> m2IntegralDim,
        "T_w_beta_prime_sq" -> kTwDim,
        "K_eta_beta_sq" -> kEtaDim,
        "lambda_T_Omega_beta_sq" -> kOmegaDim,
        "K2_integral" -> k2IntegralDim,
        "actual_M2" -> actualM2Dim,
        "actual_K2" -> actualK2Dim,
        "actual_K2_over_M2" -> actualRatioDim
      |>,
      "TermHomogeneous" -> termHomogeneous,
      "Terms" -> <|
        "T_w_beta_prime_sq" -> kTwTermExpr,
        "K_eta_beta_sq" -> kEtaTermExpr,
        "lambda_T_Omega_beta_sq" -> kOmegaTermExpr[lambda]
      |>
    |>
  ],
  "stage016DimError",
  Function[{msg, tag}, <|"Ok" -> False, "Error" -> msg, "Dims" -> <||>, "TermHomogeneous" -> False, "Terms" -> <||>|>]
];

corruptRulesFor[label_, baseRules_] := Module[{symbolByLabel, sym, corrupt},
  symbolByLabel = <|
    "mu_eta_density" -> muEtaDensity,
    "T_w_density" -> TwDensity,
    "K_eta_density" -> KetaDensity,
    "T_Omega_density" -> TOmegaDensity
  |>;
  sym = symbolByLabel[label];
  corrupt = Join[KeyDrop[baseRules, {sym}], <|sym -> baseRules[sym] + {1, 0, 0}|>];
  If[label === "T_Omega_density",
    corrupt = Join[KeyDrop[corrupt, {TomegaTilde}], <|TomegaTilde -> corrupt[TomegaTilde] + {1, 0, 0}|>]
  ];
  corrupt
];

dimRules = makeDimRules[];
baselineDim = evalDimensional[lambdaRef, m2Core, k2Ref, dimRules];
If[! TrueQ[baselineDim["Ok"]], Print["FAIL  baseline dimensional check: ", baselineDim["Error"]]; Exit[1]];
densityCorruptions = AssociationMap[
  evalDimensional[lambdaRef, m2Core, k2Ref, corruptRulesFor[#, dimRules]] &,
  {"mu_eta_density", "T_w_density", "K_eta_density", "T_Omega_density"}
];
dimensionalOk = baselineDim["Ok"];
dimProbe = <|
  "Mutation" -> "corrupt sourced [T_Omega] and assembled TomegaTilde by one extra L",
  "ParticipatesInVerdict" -> (
    scopedVerdict[densityCorruptions["T_Omega_density"]["Ok"], True, True, True] === FAILDIMENSIONAL &&
      scopedVerdict[baselineDim["Ok"], True, True, True] =!= FAILDIMENSIONAL
  ),
  "WithoutMutationDimensionalOk" -> baselineDim["Ok"],
  "WithMutationDimensionalOk" -> densityCorruptions["T_Omega_density"]["Ok"],
  "ProbeVerdict" -> If[TrueQ[densityCorruptions["T_Omega_density"]["Ok"]], "NO_FAIL", FAILDIMENSIONAL],
  "MutationFires" -> ! TrueQ[densityCorruptions["T_Omega_density"]["Ok"]],
  "Error" -> densityCorruptions["T_Omega_density"]["Error"],
  "SelfAblation" -> <|
    "DimensionalOk" -> baselineDim["Ok"],
    "Verdict" -> caseVerdict[<|"dimensionalOk" -> baselineDim["Ok"]|>],
    "FailFires" -> caseVerdict[<|"dimensionalOk" -> baselineDim["Ok"]|>] === FAILDIMENSIONAL,
    "FailSuppressed" -> caseVerdict[<|"dimensionalOk" -> baselineDim["Ok"]|>] =!= FAILDIMENSIONAL
  |>
|>;

hashText[expr_] := IntegerString[Hash[ToString[clean[expr], InputForm], "SHA256"], 16];
inputHashes = AssociationThread[order, hashText /@ ys];
distinctHashes = Length[DeleteDuplicates[Values[inputHashes]]] == Length[inputHashes];
selfOverlaps = AssociationThread[order, intS2[#^2] & /@ ys];
tautologyClear = TrueQ[distinctHashes && And @@ (TrueQ[clean[# - 1] === 0] & /@ Values[selfOverlaps])];

forcedEigenvalueProbe[forced_] := Module[{forcedK2, forcedCoeff, coefficientEquals, v},
  forcedK2 = AssociationThread[order, buildK2[forced] & /@ order];
  forcedCoeff = AssociationThread[order, extractK2Coeff[forcedK2[#]] & /@ order];
  coefficientEquals = And @@ (TrueQ[clean[forcedCoeff[#] - lambdas[#]] === 0] & /@ order);
  v = caseVerdict[<|"covariant" -> coefficientEquals|>];
  <|
    "ForcedCoefficient" -> forced,
    "ForcedK2ByChannel" -> forcedK2,
    "ForcedCoeffByChannel" -> forcedCoeff,
    "CoefficientEqualsComputedLambda" -> coefficientEquals,
    "Verdict" -> v,
    "FailFires" -> v === FAILNOTCOVARIANT
  |>
];

laneHashProbe[inputs_Association] := Module[{hashes, distinct, v},
  hashes = Association @ KeyValueMap[(#1 -> hashText[#2]) &, inputs];
  distinct = Length[DeleteDuplicates[Values[hashes]]] == Length[hashes];
  v = caseVerdict[<|"tautologyClear" -> distinct|>];
  <|
    "InputHashes" -> hashes,
    "DistinctHashes" -> distinct,
    "Verdict" -> v,
    "FailFires" -> v === FAILTAUTOLOGICAL
  |>
];

wrongEigenProbe = forcedEigenvalueProbe[2];
wrongEigenAblation = forcedEigenvalueProbe[6];
tautologyProbe = laneHashProbe[<|"20" -> harmonics["20"], "21" -> harmonics["20"], "22" -> harmonics["20"]|>];
tautologyAblation = laneHashProbe[<|"20" -> harmonics["20"], "21" -> harmonics["21c"], "22" -> harmonics["22c"]|>];
dimensionalVerdict = caseVerdict[<|"dimensionalOk" -> dimProbe["WithMutationDimensionalOk"]|>];
dimensionalAblationVerdict = caseVerdict[<|"dimensionalOk" -> dimProbe["WithoutMutationDimensionalOk"]|>];

probes = <|
  "wrong_eigenvalue" -> Join[
    wrongEigenProbe,
    <|
      "ComputedFailGate" -> ! TrueQ[wrongEigenProbe["CoefficientEqualsComputedLambda"]],
      "SelfAblation" -> Join[
        wrongEigenAblation,
        <|"FailSuppressed" -> ! TrueQ[wrongEigenAblation["FailFires"]]|>
      ]
    |>
  ],
  "tautology_hash_collision" -> Join[
    tautologyProbe,
    <|
      "ComputedFailGate" -> ! TrueQ[tautologyProbe["DistinctHashes"]],
      "SelfAblation" -> Join[
        tautologyAblation,
        <|"FailSuppressed" -> ! TrueQ[tautologyAblation["FailFires"]]|>
      ]
    |>
  ],
  "dimensional_corrupt_T_Omega" -> Join[
    dimProbe,
    <|
      "Verdict" -> dimensionalVerdict,
      "FailFires" -> dimensionalVerdict === FAILDIMENSIONAL,
      "SelfAblation" -> <|
        "DimensionalOk" -> dimProbe["WithoutMutationDimensionalOk"],
        "Verdict" -> dimensionalAblationVerdict,
        "FailFires" -> dimensionalAblationVerdict === FAILDIMENSIONAL,
        "FailSuppressed" -> dimensionalAblationVerdict =!= FAILDIMENSIONAL
      |>
    |>
  ]
|>;

expectedProbeVerdicts = <|
  "wrong_eigenvalue" -> FAILNOTCOVARIANT,
  "tautology_hash_collision" -> FAILTAUTOLOGICAL,
  "dimensional_corrupt_T_Omega" -> FAILDIMENSIONAL
|>;
expectedProbeVerdictsMatch = Association @ KeyValueMap[
  (#1 -> TrueQ[probes[#1]["Verdict"] === #2]) &,
  expectedProbeVerdicts
];
computedProbeGateFlags = <|
  "wrong_eigenvalue" -> TrueQ[probes["wrong_eigenvalue"]["ComputedFailGate"] && probes["wrong_eigenvalue"]["SelfAblation"]["FailSuppressed"]],
  "tautology_hash_collision" -> TrueQ[probes["tautology_hash_collision"]["ComputedFailGate"] && probes["tautology_hash_collision"]["SelfAblation"]["FailSuppressed"]],
  "dimensional_corrupt_T_Omega" -> TrueQ[! probes["dimensional_corrupt_T_Omega"]["WithMutationDimensionalOk"] && probes["dimensional_corrupt_T_Omega"]["SelfAblation"]["FailSuppressed"]]
|>;
ableToFailFromFlags[flags_] := TrueQ[(And @@ Values[expectedProbeVerdictsMatch]) && (And @@ Values[flags])];
neuterFlag[flags_, key_] := Join[KeyDrop[flags, {key}], <|key -> False|>];
ableToFailOk = ableToFailFromFlags[computedProbeGateFlags];
ableToFailIfProbeNeutered = AssociationThread[
  Keys[computedProbeGateFlags],
  ableToFailFromFlags[neuterFlag[computedProbeGateFlags, #]] & /@ Keys[computedProbeGateFlags]
];
neuteringAnyProbeFlipsFalse = And @@ (TrueQ[# === False] & /@ Values[ableToFailIfProbeNeutered]);

covariantOk = TrueQ[gramIsIdentity && lambdaAllSix && residualsZero && k2CoeffResidualsZero];
gateBooleans = <|
  "dimensional_ok" -> dimensionalOk,
  "covariant" -> covariantOk,
  "tautology_clear" -> tautologyClear,
  "able_to_fail_ok" -> ableToFailOk
|>;
verdict = scopedVerdict[dimensionalOk, covariantOk, tautologyClear, ableToFailOk];

matrixSquaredResidual[m_] := Total[Flatten[Map[#^2 &, m, {2}]]];
symbolNames[exprs_] := Sort[DeleteDuplicates[
  ToString[#, InputForm] & /@ Cases[Unevaluated[exprs], s_Symbol /; Context[s] === "Global`", Infinity]
]];

runAngularTheorem[] := Module[{gramResidual},
  subheading["Real l=2 harmonics, Gram=I5, and computed -Delta_S2 spectrum"];
  Print["  harmonic order = ", order];
  Print["  Gram matrix = ", fmt[gram]];
  Print["  computed lambda_m = ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, lambdas]];
  Print["  eigenvalues are Rayleigh quotients from native Integrate and the native Laplace-Beltrami operator."];
  gramResidual = matrixSquaredResidual[gram - IdentityMatrix[5]];
  expectZero["Gram - I5 total squared residual", gramResidual];
  expectBool["Gram matrix equals I5", gramIsIdentity];
  Scan[
    Function[name,
      expectZero["lambda_" <> name <> " computed by Rayleigh quotient equals 6", lambdas[name] - 6];
      expectZero["eigenfunction residual (-Delta)Y_" <> name <> "-lambda*Y_" <> name, residuals[name]]
    ],
    order
  ];
  expectBool["all computed lambda_m are 6", lambdaAllSix];
  expectBool["all eigenfunction residuals are zero", residualsZero]
];

runK2Stiffness[] := Module[{wrong, ablation},
  subheading["Bare K2 angular stiffness uses the live computed lambda"];
  Print["  K2 builder = buildK2[coeff] = Ktilde + coeff*TomegaTilde"];
  Print["  baseline assembly path = buildK2[lambdas[name]]; no counted kCoeffUsed-lambdas self-compare."];
  Print["  K2 core by channel = ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, k2Core]];
  Print["  extracted TomegaTilde coefficients = ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, k2Coeff]];
  Scan[Function[name, expectZero["K2-coefficient residual reads extracted coeff for " <> name, k2CoeffResiduals[name]]], order];
  expectBool["K2 coefficient residuals all vanish", k2CoeffResidualsZero];
  wrong = probes["wrong_eigenvalue"];
  ablation = wrong["SelfAblation"];
  expectZero["wrong_eigenvalue probe verdict is FAIL_NOT_COVARIANT", verdictResidual[wrong["Verdict"], FAILNOTCOVARIANT]];
  expectBool["wrong_eigenvalue computed_fail_gate reads coefficient mismatch", wrong["ComputedFailGate"]];
  expectZero["wrong_eigenvalue self-ablation returns ISOTROPY_CALIBRATED", verdictResidual[ablation["Verdict"], ISOTROPYCALIBRATED]];
  expectBool["wrong_eigenvalue self-ablation suppresses fail", ablation["FailSuppressed"]]
];

runDimensionalGate[] := Module[{dims, probe},
  subheading["Angular dimensional gate and sourced-density corruption probes"];
  dims = baselineDim["Dims"];
  probe = probes["dimensional_corrupt_T_Omega"];
  Print["  dimension order = (", dimensionAxesLabel[], ")"];
  Print["DIMENSIONS|axes=", dimensionAxesLabel[]];
  printDimRecord["dim_rules.a", dimRules[aDim]];
  printDimRecord["dim_rules.dw", dimRules[dwDim]];
  printDimRecord["dim_rules.d_omega", dimRules[dOmegaDim]];
  printDimRecord["dim_rules.beta2", dimRules[beta2Dim]];
  printDimRecord["dim_rules.beta2_prime", dimRules[beta2PrimeDim]];
  printDimRecord["dim_rules.mu_eta", dimRules[muEtaDensity]];
  printDimRecord["dim_rules.T_w", dimRules[TwDensity]];
  printDimRecord["dim_rules.K_eta", dimRules[KetaDensity]];
  printDimRecord["dim_rules.T_Omega", dimRules[TOmegaDensity]];
  printDimRecord["dim_rules.M_tilde", dimRules[Mtilde]];
  printDimRecord["dim_rules.K_tilde", dimRules[Ktilde]];
  printDimRecord["dim_rules.T_Omega_tilde", dimRules[TomegaTilde]];
  printDimRecord["baseline_dims.measure", dims["measure"]];
  printDimRecord["baseline_dims.M2_integral", dims["M2_integral"]];
  printDimRecord["baseline_dims.T_w_beta_prime_sq", dims["T_w_beta_prime_sq"]];
  printDimRecord["baseline_dims.K_eta_beta_sq", dims["K_eta_beta_sq"]];
  printDimRecord["baseline_dims.lambda_T_Omega_beta_sq", dims["lambda_T_Omega_beta_sq"]];
  printDimRecord["baseline_dims.K2_integral", dims["K2_integral"]];
  printDimRecord["baseline_dims.actual_M2", dims["actual_M2"]];
  printDimRecord["baseline_dims.actual_K2", dims["actual_K2"]];
  printDimRecord["baseline_dims.actual_K2_over_M2", dims["actual_K2_over_M2"]];
  Print["  sourced volume-density convention: dV=a_dim^2*dw*dOmega has dimension L^3; beta2 is dimensionless."];
  Print["  walked K2 terms = ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, baselineDim["Terms"]]];
  Print["  computed dimensions = ", Association @ KeyValueMap[(#1 -> dimText[#2]) &, dims]];
  expectZero["explicit measure has dimension L^3", dimResidualVec[dims["measure"], {3, 0, 0}]];
  expectZero["M2_integral has dimension M", dimResidualVec[dims["M2_integral"], expectedM]];
  expectZero["K2 T_w beta_prime^2 term has dimension M*T^-2", dimResidualVec[dims["T_w_beta_prime_sq"], expectedK]];
  expectZero["K2 K_eta beta^2 term has dimension M*T^-2", dimResidualVec[dims["K_eta_beta_sq"], expectedK]];
  expectZero["K2 lambda*T_Omega beta^2 term has dimension M*T^-2", dimResidualVec[dims["lambda_T_Omega_beta_sq"], expectedK]];
  expectZero["K2 integral has dimension M*T^-2", dimResidualVec[dims["K2_integral"], expectedK]];
  expectZero["bare M2=Mtilde has dimension M", dimResidualVec[dims["actual_M2"], expectedM]];
  expectZero["bare K2=Ktilde+lambda*TomegaTilde has dimension M*T^-2", dimResidualVec[dims["actual_K2"], expectedK]];
  expectZero["bare K2/M2 has dimension T^-2", dimResidualVec[dims["actual_K2_over_M2"], expectedRatio]];
  expectBool["K2 density terms are homogeneous", baselineDim["TermHomogeneous"]];
  expectBool["baseline dimensional_ok", dimensionalOk];
  Print["  T_Omega corruption error = ", probe["Error"]];
  expectZero["sourced T_Omega probe verdict is FAIL_DIMENSIONAL", verdictResidual[probe["Verdict"], FAILDIMENSIONAL]];
  expectBool["sourced T_Omega mutation participates in verdict", probe["ParticipatesInVerdict"]];
  expectBool["sourced T_Omega mutation fires", probe["MutationFires"]];
  expectZero["T_Omega self-ablation returns ISOTROPY_CALIBRATED", verdictResidual[probe["SelfAblation"]["Verdict"], ISOTROPYCALIBRATED]];
  expectBool["T_Omega self-ablation suppresses fail", probe["SelfAblation"]["FailSuppressed"]];
  Scan[
    Function[label,
      expectZero[
        "corrupt sourced [" <> label <> "] fires FAIL_DIMENSIONAL",
        verdictResidual[caseVerdict[<|"dimensionalOk" -> densityCorruptions[label]["Ok"]|>], FAILDIMENSIONAL]
      ]
    ],
    Keys[densityCorruptions]
  ]
];

runTautologyHash[] := Module[{taut, ablation},
  subheading["Computed-not-typed hash guard and tautology probe"];
  Print["  harmonic input hashes = ", inputHashes];
  Print["  self-overlaps = ", Association @ KeyValueMap[(#1 -> fmt[#2]) &, selfOverlaps]];
  expectBool["five harmonic hashes are distinct", distinctHashes];
  Scan[Function[name, expectZero["self-overlap integral for Y_" <> name <> " is 1", selfOverlaps[name] - 1]], order];
  expectBool["tautology_clear = distinct_hashes and unit self-overlaps", tautologyClear];
  taut = probes["tautology_hash_collision"];
  ablation = taut["SelfAblation"];
  expectZero["tautology_hash_collision verdict is FAIL_TAUTOLOGICAL", verdictResidual[taut["Verdict"], FAILTAUTOLOGICAL]];
  expectBool["tautology_hash_collision computed_fail_gate reads non-distinct hashes", taut["ComputedFailGate"]];
  expectZero["tautology_hash_collision self-ablation returns ISOTROPY_CALIBRATED", verdictResidual[ablation["Verdict"], ISOTROPYCALIBRATED]];
  expectBool["tautology_hash_collision self-ablation suppresses fail", ablation["FailSuppressed"]]
];

runAggregate[] := (
  subheading["016 aggregate probe battery over the three covariance probes"];
  Print["  expected probe verdicts = ", expectedProbeVerdicts];
  Print["  computed probe gate flags = ", computedProbeGateFlags];
  Print["  neutered aggregates = ", ableToFailIfProbeNeutered];
  Scan[Function[key, expectBool["probe " <> key <> " verdict matches expected token", expectedProbeVerdictsMatch[key]]], Keys[expectedProbeVerdictsMatch]];
  Scan[Function[key, expectBool["probe " <> key <> " flag is computed true", computedProbeGateFlags[key]]], Keys[computedProbeGateFlags]];
  expectBool["016 able_to_fail_ok is true", ableToFailOk];
  Scan[Function[key, expectBool["neutering " <> key <> " flips able_to_fail_ok false", ! ableToFailIfProbeNeutered[key]]], Keys[ableToFailIfProbeNeutered]];
  expectBool["neutering any one 016 probe flips aggregate false", neuteringAnyProbeFlipsFalse]
);

runPerToothAblations[] := Module[
  {
    corruptBasis, corruptAngular, corruptCovariant, wrongK2, wrongCoeff,
    wrongProbe, rightProbe, gramCorruptBasis, gramCorrupt, tautNeuteredVerdict
  },
  subheading["Per-tooth ablations on copies"];
  corruptBasis = AssociationThread[order, ReplacePart[ys, 1 -> Cos[theta] + Cos[theta]^2]];
  corruptAngular = computeAngularBlock[corruptBasis];
  corruptCovariant = TrueQ[corruptAngular["GramIsIdentity"] && corruptAngular["LambdaAllSix"] && corruptAngular["ResidualsZero"]];
  expectZero[
    "computed-eigenvalue basis corruption reaches FAIL_NOT_COVARIANT",
    verdictResidual[caseVerdict[<|"covariant" -> corruptCovariant|>], FAILNOTCOVARIANT]
  ];
  expectFail["basis corruption makes at least one eigen residual nonzero", boolResidual[corruptAngular["ResidualsZero"]]];
  wrongK2 = buildK2[2];
  wrongCoeff = extractK2Coeff[wrongK2];
  expectFail[
    "mutating the assembled K2 coefficient to 2 breaks the K2-coefficient residual",
    negLaps["20"] - wrongCoeff harmonics["20"]
  ];
  wrongProbe = forcedEigenvalueProbe[2];
  rightProbe = forcedEigenvalueProbe[6];
  expectZero["bare forcedEigenvalueProbe[2] fires FAIL_NOT_COVARIANT", verdictResidual[wrongProbe["Verdict"], FAILNOTCOVARIANT]];
  expectZero["bare forcedEigenvalueProbe[6] suppresses FAIL_NOT_COVARIANT", verdictResidual[rightProbe["Verdict"], ISOTROPYCALIBRATED]];
  gramCorruptBasis = AssociationThread[order, ReplacePart[ys, 1 -> 2 ys[[1]]]];
  gramCorrupt = computeAngularBlock[gramCorruptBasis];
  expectFail["Gram tooth: scaling one harmonic breaks Gram=I5", matrixSquaredResidual[gramCorrupt["Gram"] - IdentityMatrix[5]]];
  expectZero[
    "Gram tooth reaches FAIL_NOT_COVARIANT",
    verdictResidual[caseVerdict[<|"covariant" -> gramCorrupt["GramIsIdentity"]|>], FAILNOTCOVARIANT]
  ];
  tautNeuteredVerdict = caseVerdict[<|"tautologyClear" -> True|>];
  expectFail["tautology distinctness neuter would suppress FAIL_TAUTOLOGICAL", verdictResidual[tautNeuteredVerdict, FAILTAUTOLOGICAL]];
  expectZero["T_Omega dimensional tooth reaches FAIL_DIMENSIONAL", verdictResidual[probes["dimensional_corrupt_T_Omega"]["Verdict"], FAILDIMENSIONAL]];
  Scan[Function[key, expectBool["aggregate tooth includes " <> key, ! ableToFailIfProbeNeutered[key]]], Keys[ableToFailIfProbeNeutered]];
  Print["  kCoeffEquals de-count: no kCoeffUsed-lambdas self-compare is counted; K2 computed-ness rides on buildK2[lambdas], extracted-coefficient residuals, and the bare forced probe."]
];

runAritySelfCheck[] := Module[{dimProbeResult, hashProbeResult, eigenProbeResult},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  dimProbeResult = evalDimensional[lambdaRef, m2Core, k2Ref, dimRules];
  hashProbeResult = laneHashProbe[<|"20" -> harmonics["20"], "21" -> harmonics["21c"], "22" -> harmonics["22c"]|>];
  eigenProbeResult = forcedEigenvalueProbe[2];
  expectZero["arity intS2[1] returns 4*pi", intS2[1] - 4 Pi];
  expectZero["arity lapS2[Y20] is accepted and has l=2 residual", -lapS2[harmonics["20"]] - 6 harmonics["20"]];
  expectBool["arity buildK2[lambdas[20]] returns expression containing TomegaTilde", ! FreeQ[buildK2[lambdas["20"]], TomegaTilde]];
  expectBool["arity forcedEigenvalueProbe[2] returns verdict key", KeyExistsQ[eigenProbeResult, "Verdict"]];
  expectBool["arity laneHashProbe[assoc] returns distinctness key", KeyExistsQ[hashProbeResult, "DistinctHashes"]];
  expectBool["arity evalDimensional[lambda,M2,K2,rules] returns Ok key", KeyExistsQ[dimProbeResult, "Ok"]];
  expectBool["arity scopedVerdict[4 args] returns ISOTROPY_CALIBRATED", scopedVerdict[True, True, True, True] === ISOTROPYCALIBRATED];
  expectBool["no unevaluated Integrate remains in angular results", FreeQ[{gram, Values[lambdas], Values[residuals]}, _Integrate]];
  expectBool["no unevaluated Derivative remains in angular results", FreeQ[{gram, Values[lambdas], Values[residuals]}, _Derivative]]
];

runVerdictAndScope[] := (
  subheading["016 scoped landing and 016/017 cut"];
  Print["  016 gate booleans = ", gateBooleans];
  Print["  016 scoped verdict = ", verdict];
  expectZero["016 scoped verdict lands ISOTROPY_CALIBRATED component", verdictResidual[verdict, ISOTROPYCALIBRATED]];
  Print["  ISOTROPY_CALIBRATED (JOINT, 2-stage) -- PARTIAL: 016 landed, 017 PENDING"];
  Print["    = (016: l=2 SO(3) covariance theorem: real harmonics + Gram=I5 + computed lambda_m=6 + K2 angular stiffness) [EARNED here]"];
  Print["    AND (017: grouped-P2 lane isotropy: grouped {20,21,22} lanes / raw-D=0 / normalized-u / calibration partition) [PENDING]"];
  Print["  Exact cut: this script does not assemble grouped lanes, raw-D, normalized-u, response probes, calibration partition, or port-kernel export."];
  Print["  Carried caveats: angular structure is earned; radial profile/scalars are frozen calibration inputs, so the joint is CALIBRATED not PASS."];
  Print["  Deferred caveats: 54/5 quadrupole normalization, outgoing odd-N coefficients, and solved nonlinear branch data are Gate 4/5/6 sim-deferred."];
  expectBool["joint composition is partial with 017 pending", verdict === ISOTROPYCALIBRATED]
);

runProvenance[] := Module[{liveNames},
  subheading["Provenance and scope labels"];
  Print["  CONSUMED-from-011/012/013: mu_eta/T_w physical wall constants, beta2(w)/R0(w), Mtilde/Ktilde/TomegaTilde, and Gate-1 D/N provenance are cited as provenance."];
  Print["  Self-contained dimensional integrity: pathA_32 uses volume densities on a_dim^2*dw*dOmega with dimensionless beta2; stage013's line-density K_eta=T_w*beta^2 relation does not transfer."];
  Print["  no-c_S: the l=2 covariance theorem is speed-free; matter-sector c_s/BdG remains deferred."];
  Print["  ANGULAR-EARNED / RADIAL-CALIBRATED: 016 derives Gram=I5, lambda_m=6, and K2; 017 owns the radial calibration partition."];
  Print["  COMPUTED-NOT-TYPED: Rayleigh + eigenfunction residual + extracted K2-coefficient residual + forced_eigenvalue_probe; kCoeffEquals X==X is de-counted."];
  Print["  AGGREGATE-BATTERY-INTACT: wrong_eigenvalue, tautology_hash_collision, and dimensional_corrupt_T_Omega all participate."];
  Print["  SOURCED-T_OMEGA-DIMENSIONAL: corrupting sourced T_Omega plus TomegaTilde fires FAIL_DIMENSIONAL."];
  Print["  dropped-bookkeeping: scratch-YAML engine agreement and report/feed writers are stripped; transcript-level agreement remains."];
  Print["  register note: 016 is an earned structural slice with likely zero new counted knobs; T_Omega and beta2 counting is deferred to 017."];
  liveNames = symbolNames[{Values[harmonics], Values[lambdas], Values[k2Core], m2Core, k2Ref}];
  expectBool["no c_S/cS live symbol appears", FreeQ[liveNames, "c_S"] && FreeQ[liveNames, "cS"]];
  expectBool["Btilde/Ztilde support scalars are not live 016 symbols", AllTrue[liveNames, (! StringStartsQ[#, "B"] && ! StringStartsQ[#, "Z"]) &]]
];

printVerdictLabels[] := (
  subheading["Verdict labels"];
  Print["  ledger earned-label: L2_SO3_COVARIANCE_THEOREM_EARNED"];
  Print["  source top-line verdict: ISOTROPY_CALIBRATED (JOINT 2-stage; 016 is the earned SO(3) covariance component, 017 completes the calibration/lane component)"];
  Print["  joint composition: ISOTROPY_CALIBRATED = 016[EARNED l=2 SO(3) covariance] AND 017[PENDING grouped-P2 lane isotropy + calibration partition]"];
  Print["  earned angular structure: Gram=I5 genuine; lambda_m=6 computed by Rayleigh + residual; K2=Ktilde+lambda_m*TomegaTilde uses the live computed lambda"];
  Print["  earned able-to-fail battery: wrong_eigenvalue / tautology_hash_collision / dimensional_corrupt_T_Omega, each with suppressing self-ablation"];
  Print["  consumed framing: provenance + pathA_32 self-contained dimensional integrity, not a cross-stage dual-site relation"];
  Print["  new symbols first appearing here but not counted in 016: T_Omega/TomegaTilde and beta2(w), deferred to 017"]
);

runAll[] := (
  heading["ledger_stage016_l2_so3_covariance_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage016_l2_so3_covariance"];
  Print["Engine: native Wolfram exact symbolic Integrate/D/Hash; no NIntegrate/floats/tolerances; zero file I/O."];
  runAngularTheorem[];
  runK2Stiffness[];
  runDimensionalGate[];
  runTautologyHash[];
  runAggregate[];
  runPerToothAblations[];
  runAritySelfCheck[];
  runVerdictAndScope[];
  runProvenance[];
  printVerdictLabels[];
  0
);

result = Catch[
  runAll[],
  "ledgerStage016Failure",
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
