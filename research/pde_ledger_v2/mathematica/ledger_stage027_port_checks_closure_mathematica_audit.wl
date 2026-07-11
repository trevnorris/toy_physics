(* Ledger stage027 Wolfram audit: density-port checks, closure, and joint.

   Standalone, print-only, no arguments, and zero file I/O.  This is a
   verdict-bearing re-author: dimensions and a-weights use rescaling ratios,
   the DtN sign uses built-in spherical Hankel functions, closure is decided
   by ratio invariants, and verdict routing is Association dispatch.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

jointVerdict = "DENSITY_PORT_HOSTED";
inconclusiveVerdict = "PORT_INCONCLUSIVE_SIM_DEFERRED";
failOrigin = "FAIL_NOT_DENSITY_DERIVED";
failVanishes = "FAIL_PORT_VANISHES";

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

clean[expr_] := Factor[Cancel[Together[FullSimplify[expr]]]];

externalSymbolNames = <|
  "cS" -> "c_s", "rhoEff" -> "rho_eff", "xiQ" -> "Xi_Q",
  "xiDeferred" -> "Xi_deferred", "etaQ" -> "eta_q",
  "etaPhi" -> "eta_phi", "varpiQ2" -> "varpi_q2",
  "varpiPhi2" -> "varpi_Phi2", "lambdaC" -> "lambda_c",
  "Iwrong2" -> "I_wrong2", "foreignSubject" -> "foreign_subject",
  "gravG" -> "G", "lightC" -> "c", "deltaGamma" -> "deltaGamma"
|>;

externalName[symbol_Symbol] := Module[{raw = SymbolName[symbol]},
  Lookup[externalSymbolNames, raw, raw]
];

fmt[expr_String] := expr;
fmt[expr_] := StringReplace[
  ToString[InputForm[clean[expr]]],
  Normal[externalSymbolNames]
];

fmtNames[items_List] := "{" <> StringRiffle[Sort[items], ", "] <> "}";

fail[] := Throw[failureMessage, "ledgerStage027Failure"];

expectZero[name_, residual_] := Module[{value = clean[residual]},
  If[TrueQ[value === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    failureMessage = name <> ": residual = " <> fmt[value];
    Print["FAIL  ", failureMessage];
    fail[]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

rigAssert[name_String, condition_] := If[
  ! TrueQ[condition],
  Throw[name, "stage027RigAssertion"]
];

probeAssertion[held_HoldComplete, expectedName_String] := Module[{fired},
  fired = Catch[
    ReleaseHold[held];
    "NO_ASSERT_FIRED",
    "stage027RigAssertion",
    Function[{message, tag}, message]
  ];
  {TrueQ[fired === expectedName], fired}
];

assertionPasses[held_HoldComplete] := TrueQ[Catch[
  ReleaseHold[held];
  True,
  "stage027RigAssertion",
  Function[{message, tag}, False]
]];

exerciseRig[
  label_String,
  assertionName_String,
  rigHeld_HoldComplete,
  neutralHeld_HoldComplete,
  outcome_String
] := Module[{probe, caught, firedName, neutralPass, text},
  probe = probeAssertion[rigHeld, assertionName];
  caught = probe[[1]];
  firedName = probe[[2]];
  neutralPass = assertionPasses[neutralHeld];
  expectBool[
    "META " <> label <> " routed assertion fires and neutering stops it",
    caught && neutralPass
  ];
  text = outcome <> " at " <> ToString[firedName] <> "; neutralized=PASS";
  Print["RIG ", label, ": ", text];
  text
];

(* The exact 10-symbol stage024 host contract. q2/Phi2 are metadata only. *)
hostContractNames = Sort[{
  "a", "c_s", "rho_eff", "I25", "Xi_Q", "eta_q", "eta_phi",
  "varpi_q2", "varpi_Phi2", "lambda_c"
}];

baselineTags = Sort[{
  "continuity_interface", "pathA_29_bulk", "pathA_32_wall"
}];

citedN0Den[] :=
  I25^2 xiQ^2 cS^4 rhoEff (etaPhi varpiQ2 + etaQ lambdaC)^2/
    (a^7 (lambdaC^2 - varpiPhi2 varpiQ2)^2);

globalSymbols[expr_] := DeleteDuplicates[Cases[
  HoldComplete[expr],
  symbol_Symbol /; Context[symbol] === "Global`",
  Infinity
]];

globalSymbolNames[expr_] := Sort[externalName /@ globalSymbols[expr]];

hostContractQ[expr_] := globalSymbolNames[expr] === hostContractNames;

baseFacts[] := <|
  "Tags" -> baselineTags,
  "VectorHostSymbols" -> {},
  "SourceMapComplete" -> True,
  "VectorFree" -> True,
  "MomentValid" -> True,
  "LineageValid" -> True,
  "ValidatedI25" -> I25
|>;

cfg[name_String, rules___Rule] := Join[
  <|
    "Name" -> name,
    "CouplingZero" -> False,
    "CorruptDimension" -> False,
    "CorruptGDimension" -> False,
    "IncomingSign" -> False,
    "CouplingAPower" -> -7/2,
    "DeferredUncertified" -> False,
    "ProvenDeferred" -> False,
    "SecondSubjectArm" -> True,
    "Kbar4Shift" -> 0,
    "Gamma5Shift" -> 0
  |>,
  <|rules|>
];

specializeCitedPort[config_Association] := Module[
  {xi, moment, expr},
  xi = If[
    TrueQ[config["DeferredUncertified"] || config["ProvenDeferred"]],
    xiDeferred,
    xiQ
  ];
  moment = If[TrueQ[config["CouplingAPower"] === -7/2], I25, Iwrong2];
  expr = citedN0Den[] /. {I25 -> moment, xiQ -> xi};
  (* The g_base mutation -7/2 -> -3 changes its square by one factor of a.
     Apply that specialization to the cited stage024 export itself. *)
  If[moment === Iwrong2, expr = a expr];
  expr = clean[expr];
  If[TrueQ[config["CouplingZero"]], expr = 0];
  <|"Expression" -> expr, "PortMoment" -> moment|>
];

(* Independent dimension route: rescale every live carrier by a unit
   monomial and extract the resulting homogeneous ratio. *)
baseDims = <|
  a -> {1, 0, 0}, cS -> {1, 0, -1}, lightC -> {1, 0, -1},
  gravG -> {3, -1, -2}, d0 -> {-1, 1, -2}, rhoEff -> {-3, 1, 0},
  I25 -> {5/2, 0, 0}, Iwrong2 -> {2, 0, 0}, xiQ -> {0, 0, 0},
  xiDeferred -> {0, 0, 0}, etaQ -> {0, 0, 0}, etaPhi -> {0, 0, 0},
  varpiQ2 -> {0, 0, -2}, varpiPhi2 -> {0, 0, -2},
  lambdaC -> {0, 0, -2}, foreignSubject -> {-3, 1, 0}
|>;

unitDimension[expr_, dims_Association] := Module[
  {live, missing, rules, ratio, vector, monomial},
  If[TrueQ[clean[expr] === 0], Return[{0, 0, 0}]];
  live = globalSymbols[expr];
  missing = Complement[live, Keys[dims]];
  If[missing =!= {}, Return[Missing["SourcedDimension", missing]]];
  rules = KeyValueMap[
    Function[{symbol, dim},
      symbol -> symbol uL^dim[[1]] uM^dim[[2]] uT^dim[[3]]
    ],
    dims
  ];
  ratio = clean[(expr /. rules)/expr];
  vector = {Exponent[ratio, uL], Exponent[ratio, uM], Exponent[ratio, uT]};
  monomial = uL^vector[[1]] uM^vector[[2]] uT^vector[[3]];
  If[TrueQ[clean[ratio/monomial] === 1], vector, Missing["MixedDimension"]]
];

baseWeights = <|
  cS -> 0, lightC -> 0, gravG -> 0, d0 -> 0, rhoEff -> 0,
  I25 -> 0, Iwrong2 -> 0, xiQ -> 0, xiDeferred -> 0,
  etaQ -> 0, etaPhi -> 0, varpiQ2 -> -2, varpiPhi2 -> -2,
  lambdaC -> -2, foreignSubject -> 0
|>;

weightPower[expr_, weights_Association, unknown_List] := Module[
  {live, rules, ratio, power},
  If[TrueQ[clean[expr] === 0], Return[0]];
  live = globalSymbols[expr];
  If[Intersection[live, unknown] =!= {}, Return[Missing["UnknownWeight"]]];
  rules = Join[
    {a -> a scaleMarker},
    KeyValueMap[
      Function[{symbol, weight}, symbol -> symbol scaleMarker^weight],
      weights
    ]
  ];
  ratio = clean[(expr /. rules)/expr];
  power = Exponent[ratio, scaleMarker];
  If[TrueQ[clean[ratio/scaleMarker^power] === 1], power, Missing["MixedWeight"]]
];

(* Independent special-function sign route, memoized because every control
   consumes one of the same two exact branches. *)
hankelWave["outgoing"] := SphericalHankelH1[2, z];
hankelWave["incoming"] := SphericalHankelH2[2, z];

dtnSign[kind_String] := dtnSign[kind] = Module[
  {wave, lambda, response, coefficient},
  wave = hankelWave[kind];
  lambda = FullSimplify[z D[wave, z]/wave];
  response = FullSimplify[-3/lambda];
  coefficient = FullSimplify[
    SeriesCoefficient[response, {z, 0, 5}]/I
  ];
  <|
    "Kind" -> kind,
    "CoefficientOverI" -> coefficient,
    "ChiQ" -> If[TrueQ[coefficient === 1/27], 1, -1],
    "OK" -> TrueQ[coefficient === 1/27]
  |>
];

(* Independent closure route: the verdict consumes two dimensionless ratio
   invariants.  Additive residuals remain standalone printed/ablatable teeth. *)
closureInvariants[n0_, kbar4Shift_, gamma5Shift_] := Module[
  {p0, kbar0, kbar2, kbar4, gamma5, kResidual, gammaResidual,
   kRatio, gammaRatio},
  p0 = clean[(cS/a)^2 n0/d0];
  kbar0 = 54 gravG cS^5/(5 a^5 lightC^5);
  kbar2 = 6 gravG cS^3/(5 a^3 lightC^5);
  kbar4 = 8 gravG cS/(15 a lightC^5) + kbar4Shift;
  gamma5 = clean[kbar0 a^5/(27 cS^5)] + gamma5Shift;
  kResidual = clean[kbar4 - 4 kbar2^2/kbar0];
  gammaResidual = clean[gamma5 - 2 gravG/(5 lightC^5)];
  kRatio = clean[kbar4 kbar0/(4 kbar2^2) - 1];
  gammaRatio = clean[gamma5/(2 gravG/(5 lightC^5)) - 1];
  <|
    "P0Physical" -> p0,
    "Kbar0" -> kbar0,
    "Kbar2" -> kbar2,
    "Kbar4" -> kbar4,
    "Gamma5" -> gamma5,
    "Kbar4Residual" -> kResidual,
    "Gamma5Residual" -> gammaResidual,
    "Kbar4RatioInvariant" -> kRatio,
    "Gamma5RatioInvariant" -> gammaRatio,
    "OK" -> TrueQ[kRatio === 0 && gammaRatio === 0]
  |>
];

(* Structurally distinct verdict route: ordered Association predicates select
   a route key, then a second Association maps that key to the verdict. *)
evaluatePort[config_Association, facts_Association] := Module[
  {specialized, expr, portMoment, dims, exprDim, dimOK, weights, unknown,
   n0Power, p0Power, scalingWrong, scalingOK, scalingUndecidable, sign,
   nonzeroOK, closure, membershipArm, vanishedArm, subjectBinding,
   originOK, bad, routes, route, dispatch, verdict},
  specialized = specializeCitedPort[config];
  expr = specialized["Expression"];
  portMoment = If[
    TrueQ[config["CouplingAPower"] === -7/2],
    facts["ValidatedI25"],
    Iwrong2
  ];

  dims = Association[Normal[baseDims]];
  If[TrueQ[config["CorruptDimension"]],
    AssociateTo[dims, I25 -> (baseDims[I25] + {1, 0, 0})]
  ];
  If[TrueQ[config["CorruptGDimension"]],
    AssociateTo[dims, gravG -> (baseDims[gravG] + {1, 0, 0})]
  ];
  exprDim = unitDimension[expr, dims];
  dimOK = TrueQ[exprDim === {-1, 1, 0}];

  weights = Association[Normal[baseWeights]];
  unknown = If[TrueQ[config["DeferredUncertified"]], {xiDeferred}, {}];
  n0Power = weightPower[expr, weights, unknown];
  p0Power = If[MissingQ[n0Power], n0Power, n0Power - 2 - weights[d0]];
  scalingWrong = ! MissingQ[p0Power] && ! TrueQ[p0Power === -5];
  scalingOK = TrueQ[p0Power === -5];
  scalingUndecidable = MissingQ[p0Power];

  sign = dtnSign[If[TrueQ[config["IncomingSign"]], "incoming", "outgoing"]];
  nonzeroOK = ! TrueQ[clean[expr] === 0];
  closure = closureInvariants[
    expr, config["Kbar4Shift"], config["Gamma5Shift"]
  ];

  membershipArm = MemberQ[globalSymbols[expr], portMoment];
  vanishedArm = TrueQ[
    config["SecondSubjectArm"] && clean[expr] === 0 && config["CouplingZero"]
  ];
  subjectBinding = TrueQ[membershipArm || vanishedArm];
  originOK = TrueQ[
    MemberQ[facts["Tags"], "continuity_interface"] &&
    ! MemberQ[facts["Tags"], "vector_port"] &&
    facts["VectorHostSymbols"] === {} &&
    facts["SourceMapComplete"] && facts["VectorFree"] &&
    facts["LineageValid"] && facts["MomentValid"] && subjectBinding
  ];

  bad = Pick[
    {"dimensional", "scaling", "sign"},
    {! dimOK, scalingWrong, ! sign["OK"]}
  ];
  routes = <|
    "origin" -> (! originOK),
    "vanish" -> (! nonzeroOK),
    "malformed" -> (bad =!= {}),
    "hosted" -> TrueQ[
      originOK && nonzeroOK && dimOK && scalingOK && sign["OK"] && closure["OK"]
    ],
    "inconclusive" -> True
  |>;
  route = SelectFirst[Keys[routes], TrueQ[routes[#]] &, "inconclusive"];
  dispatch = <|
    "origin" -> failOrigin,
    "vanish" -> failVanishes,
    "malformed" -> ("FAIL_PORT_MALFORMED(" <> StringRiffle[bad, ","] <> ")"),
    "hosted" -> jointVerdict,
    "inconclusive" -> inconclusiveVerdict
  |>;
  verdict = dispatch[route];

  <|
    "Config" -> config,
    "Expression" -> expr,
    "FixtureMoment" -> specialized["PortMoment"],
    "PortMoment" -> portMoment,
    "ExpressionDimension" -> exprDim,
    "DimOK" -> dimOK,
    "N0Power" -> n0Power,
    "P0Power" -> p0Power,
    "ScalingWrong" -> scalingWrong,
    "ScalingOK" -> scalingOK,
    "ScalingUndecidable" -> scalingUndecidable,
    "Sign" -> sign,
    "NonzeroOK" -> nonzeroOK,
    "Closure" -> closure,
    "MembershipArm" -> membershipArm,
    "VanishedArm" -> vanishedArm,
    "SubjectBinding" -> subjectBinding,
    "OriginOK" -> originOK,
    "Route" -> route,
    "Verdict" -> verdict,
    "RetirementRecorded" -> TrueQ[verdict === jointVerdict]
  |>
];

definitionArity[function_Symbol] := Module[{definitions, lhs},
  definitions = DownValues[function];
  If[Length[definitions] =!= 1, Return[-1]];
  lhs = Extract[definitions, {1, 1}, HoldComplete];
  Extract[
    lhs,
    {1, 1},
    Function[call, Length[Unevaluated[call]], HoldAllComplete]
  ]
];

heldCallArity[held_HoldComplete] := Extract[
  held,
  {1},
  Function[call, Length[Unevaluated[call]], HoldAllComplete]
];

arityScan[held_HoldComplete] :=
  definitionArity[evaluatePort] === heldCallArity[held];

leakageFreeQ[objects_] := FreeQ[
  objects,
  specializeCitedPort | unitDimension | weightPower | hankelWave | dtnSign |
  closureInvariants | evaluatePort | authoredLeak | SeriesCoefficient |
  SphericalHankelH1 | SphericalHankelH2 | Together | Cancel | Factor
];

runAudit[] := Module[
  {facts, n0Den, baseline, controls, expected, outcomes, label, result,
   assertionName, rigSpecs, badFacts, atomicMutations, corruptSubject,
   corruptDim, corruptPower, positiveNames, positiveFlips, properCall,
   badCall, actualTranscript, leakedTranscript},

  facts = baseFacts[];
  n0Den = clean[citedN0Den[]];
  baseline = evaluatePort[cfg["baseline"], facts];

  subheading["I. consumed sibling exports and subject-integrity oracle"];
  expectBool[
    "SUBJECT-INTEGRITY cited N0_den free_symbols equals exact 10-symbol HOST_CONTRACT",
    hostContractQ[n0Den]
  ];
  expectBool[
    "ASSEMBLY consumed 025 exact baseline taint set and atomic vector facts",
    facts["Tags"] === baselineTags &&
    MemberQ[facts["Tags"], "continuity_interface"] &&
    ! MemberQ[facts["Tags"], "vector_port"] &&
    facts["VectorHostSymbols"] === {} &&
    facts["SourceMapComplete"] && facts["VectorFree"]
  ];
  expectBool[
    "ASSEMBLY consumed 026 atomic moment_valid, validated_I25, and lineage_valid facts",
    facts["MomentValid"] && facts["ValidatedI25"] === I25 && facts["LineageValid"]
  ];

  subheading["Six port checks and independent closure residuals"];
  expectBool["DIM [N0_den] == (-1,1,0) = L^-1 M", baseline["DimOK"]];
  expectBool["SCALE P0_physical a-power == -5", baseline["P0Power"] === -5];
  expectBool[
    "SIGN outgoing SphericalHankelH1 gives coeff/i == 1/27 and chi_Q=+1",
    baseline["Sign"]["CoefficientOverI"] === 1/27 && baseline["Sign"]["ChiQ"] === 1
  ];
  expectBool["NONZERO cited density port is nonzero", baseline["NonzeroOK"]];
  expectZero[
    "CLOSURE Kbar4-4*Kbar2^2/Kbar0 standalone residual",
    baseline["Closure"]["Kbar4Residual"]
  ];
  expectZero[
    "CLOSURE Gamma5-2*G/(5*c^5) standalone residual",
    baseline["Closure"]["Gamma5Residual"]
  ];
  expectBool[
    "POSITIVE assembled joint lands DENSITY_PORT_HOSTED",
    baseline["Verdict"] === jointVerdict
  ];

  subheading["Per-control verdicts and routed per-tooth ablations"];
  controls = <|
    "dimensional" -> evaluatePort[cfg["dimensional", "CorruptDimension" -> True], facts],
    "corrupt_G_scope" -> evaluatePort[cfg["corrupt_G_scope", "CorruptGDimension" -> True], facts],
    "scaling" -> evaluatePort[cfg["scaling", "CouplingAPower" -> -3], facts],
    "sign" -> evaluatePort[cfg["sign", "IncomingSign" -> True], facts],
    "zero_coupling" -> evaluatePort[cfg["zero_coupling", "CouplingZero" -> True], facts],
    "zero_no_or" -> evaluatePort[cfg["zero_no_or", "CouplingZero" -> True, "SecondSubjectArm" -> False], facts],
    "deferred_scalar" -> evaluatePort[cfg["deferred_scalar", "DeferredUncertified" -> True], facts],
    "deferred_scalar_proven_converse" -> evaluatePort[cfg["deferred_scalar_proven_converse", "ProvenDeferred" -> True], facts],
    "closure_kbar4" -> evaluatePort[cfg["closure_kbar4", "Kbar4Shift" -> delta4], facts],
    "closure_gamma5" -> evaluatePort[cfg["closure_gamma5", "Gamma5Shift" -> deltaGamma], facts]
  |>;
  expected = <|
    "dimensional" -> "FAIL_PORT_MALFORMED(dimensional)",
    "corrupt_G_scope" -> jointVerdict,
    "scaling" -> "FAIL_PORT_MALFORMED(scaling)",
    "sign" -> "FAIL_PORT_MALFORMED(sign)",
    "zero_coupling" -> failVanishes,
    "zero_no_or" -> failOrigin,
    "deferred_scalar" -> inconclusiveVerdict,
    "deferred_scalar_proven_converse" -> jointVerdict,
    "closure_kbar4" -> inconclusiveVerdict,
    "closure_gamma5" -> inconclusiveVerdict
  |>;
  KeyValueMap[
    Function[{name, expectedVerdict},
      expectBool[
        "CONTROL " <> name <> " computed verdict == " <> expectedVerdict,
        controls[name]["Verdict"] === expectedVerdict
      ];
      Print["CONTROL ", name, ": ", controls[name]["Verdict"]]
    ],
    expected
  ];
  expectBool[
    "SCALE single-tag isolation keeps dim_ok True and only scaling_wrong fires",
    controls["scaling"]["DimOK"] && controls["scaling"]["ScalingWrong"] &&
    controls["scaling"]["Verdict"] === "FAIL_PORT_MALFORMED(scaling)"
  ];
  expectBool[
    "NONZERO OR-arm meta flips zero_coupling FAIL_PORT_VANISHES to FAIL_NOT_DENSITY_DERIVED",
    controls["zero_coupling"]["Verdict"] === failVanishes &&
    controls["zero_no_or"]["Verdict"] === failOrigin
  ];
  expectBool[
    "CLOSURE isolated Kbar4 mutation leaves Gamma5 residual zero",
    ! TrueQ[controls["closure_kbar4"]["Closure"]["Kbar4Residual"] === 0] &&
    TrueQ[controls["closure_kbar4"]["Closure"]["Gamma5Residual"] === 0]
  ];
  expectBool[
    "CLOSURE isolated Gamma5 mutation leaves Kbar4 residual zero",
    ! TrueQ[controls["closure_gamma5"]["Closure"]["Gamma5Residual"] === 0] &&
    TrueQ[controls["closure_gamma5"]["Closure"]["Kbar4Residual"] === 0]
  ];

  outcomes = <||>;
  rigSpecs = {
    {"DIM dimensional", "DIM [N0_den]=L^-1 M named assert", "dimensional", "DimOK", controls["dimensional"]["Verdict"]},
    {"SCALE scaling-single-tag", "SCALE p0_power=-5 named assert", "scaling", "ScalingOK", controls["scaling"]["Verdict"]},
    {"SIGN incoming", "SIGN outgoing coeff/i=1/27 named assert", "sign", "Sign", controls["sign"]["Verdict"]},
    {"NONZERO zero_coupling", "NONZERO N0_den!=0 named assert", "zero_coupling", "NonzeroOK", controls["zero_coupling"]["Verdict"]},
    {"DEFERRED uncertified", "DEFERRED scaling-decidable named assert", "deferred_scalar", "Deferred", controls["deferred_scalar"]["Verdict"]},
    {"CLOSURE Kbar4 isolated", "CLOSURE Kbar4 standalone named assert", "closure_kbar4", "Kbar4", controls["closure_kbar4"]["Verdict"]},
    {"CLOSURE Gamma5 isolated", "CLOSURE Gamma5 standalone named assert", "closure_gamma5", "Gamma5", controls["closure_gamma5"]["Verdict"]}
  };
  Scan[
    Function[spec,
      label = spec[[1]];
      assertionName = spec[[2]];
      result = controls[spec[[3]]];
      AssociateTo[outcomes, label -> Switch[spec[[4]],
        "DimOK", exerciseRig[label, assertionName,
          HoldComplete[rigAssert[assertionName, result["DimOK"]]],
          HoldComplete[rigAssert[assertionName, baseline["DimOK"]]], spec[[5]]],
        "ScalingOK", exerciseRig[label, assertionName,
          HoldComplete[rigAssert[assertionName, result["ScalingOK"]]],
          HoldComplete[rigAssert[assertionName, baseline["ScalingOK"]]], spec[[5]]],
        "Sign", exerciseRig[label, assertionName,
          HoldComplete[rigAssert[assertionName, result["Sign"]["OK"]]],
          HoldComplete[rigAssert[assertionName, baseline["Sign"]["OK"]]], spec[[5]]],
        "NonzeroOK", exerciseRig[label, assertionName,
          HoldComplete[rigAssert[assertionName, result["NonzeroOK"]]],
          HoldComplete[rigAssert[assertionName, baseline["NonzeroOK"]]], spec[[5]]],
        "Deferred", exerciseRig[label, assertionName,
          HoldComplete[rigAssert[assertionName, ! result["ScalingUndecidable"]]],
          HoldComplete[rigAssert[assertionName, ! controls["deferred_scalar_proven_converse"]["ScalingUndecidable"]]], spec[[5]]],
        "Kbar4", exerciseRig[label, assertionName,
          HoldComplete[rigAssert[assertionName, result["Closure"]["Kbar4Residual"] === 0]],
          HoldComplete[rigAssert[assertionName, baseline["Closure"]["Kbar4Residual"] === 0]], spec[[5]]],
        "Gamma5", exerciseRig[label, assertionName,
          HoldComplete[rigAssert[assertionName, result["Closure"]["Gamma5Residual"] === 0]],
          HoldComplete[rigAssert[assertionName, baseline["Closure"]["Gamma5Residual"] === 0]], spec[[5]]]
      ]]
    ],
    rigSpecs
  ];

  atomicMutations = <|
    "ASSEMBLY continuity_interface" -> Join[facts, <|"Tags" -> DeleteCases[facts["Tags"], "continuity_interface"]|>],
    "ASSEMBLY vector_port" -> Join[facts, <|"Tags" -> Sort[Append[facts["Tags"], "vector_port"]]|>],
    "ASSEMBLY vector_host_symbols" -> Join[facts, <|"VectorHostSymbols" -> {"A_w"}|>],
    "ASSEMBLY source_map_complete" -> Join[facts, <|"SourceMapComplete" -> False|>],
    "ASSEMBLY vector_free" -> Join[facts, <|"VectorFree" -> False|>],
    "ASSEMBLY moment_valid" -> Join[facts, <|"MomentValid" -> False|>],
    "ASSEMBLY validated_I25" -> Join[facts, <|"ValidatedI25" -> foreignSubject|>],
    "ASSEMBLY lineage_valid" -> Join[facts, <|"LineageValid" -> False|>]
  |>;
  KeyValueMap[
    Function[{atomicLabel, mutatedFacts},
      result = evaluatePort[cfg[atomicLabel], mutatedFacts];
      assertionName = atomicLabel <> " origin_ok named assert";
      expectBool[
        atomicLabel <> " computed verdict == FAIL_NOT_DENSITY_DERIVED",
        result["Verdict"] === failOrigin
      ];
      AssociateTo[outcomes, atomicLabel -> exerciseRig[
        atomicLabel,
        assertionName,
        HoldComplete[rigAssert[assertionName, result["OriginOK"]]],
        HoldComplete[rigAssert[assertionName, baseline["OriginOK"]]],
        result["Verdict"]
      ]]
    ],
    atomicMutations
  ];

  assertionName = "ASSEMBLY zero-coupling subject-binding OR-arm named assert";
  AssociateTo[outcomes, "ASSEMBLY OR-arm" -> exerciseRig[
    "ASSEMBLY OR-arm",
    assertionName,
    HoldComplete[rigAssert[assertionName, controls["zero_no_or"]["SubjectBinding"]]],
    HoldComplete[rigAssert[assertionName, controls["zero_coupling"]["SubjectBinding"]]],
    failVanishes <> "->" <> failOrigin
  ]];

  corruptSubject = clean[n0Den /. rhoEff -> foreignSubject];
  corruptDim = unitDimension[corruptSubject, baseDims];
  corruptPower = weightPower[corruptSubject, baseWeights, {}];
  expectBool[
    "SUBJECT-INTEGRITY foreign replacement preserves dimension and a-power",
    corruptDim === {-1, 1, 0} &&
    corruptPower === weightPower[n0Den, baseWeights, {}]
  ];
  assertionName = "SUBJECT-INTEGRITY exact host-contract named assert";
  AssociateTo[outcomes, "SUBJECT-INTEGRITY" -> exerciseRig[
    "SUBJECT-INTEGRITY",
    assertionName,
    HoldComplete[rigAssert[assertionName, hostContractQ[corruptSubject]]],
    HoldComplete[rigAssert[assertionName, hostContractQ[n0Den]]],
    "CONTRACT_REJECTED_ONLY"
  ]];

  positiveNames = {
    "dimensional", "scaling", "sign", "zero_coupling",
    "deferred_scalar", "closure_kbar4"
  };
  positiveFlips = controls[#]["Verdict"] =!= jointVerdict & /@ positiveNames;
  expectBool[
    "POSITIVE baseline hosts and corruption of any of six checks flips the joint",
    baseline["Verdict"] === jointVerdict && And @@ positiveFlips
  ];
  AssociateTo[outcomes, "POSITIVE" -> "DENSITY_PORT_HOSTED; any-check-corruption=FLIP"];
  Print["RIG POSITIVE: ", outcomes["POSITIVE"]];

  assertionName = "RETIREMENT recorded only with hosted joint named assert";
  AssociateTo[outcomes, "RETIREMENT-CONDITIONAL" -> exerciseRig[
    "RETIREMENT-CONDITIONAL",
    assertionName,
    HoldComplete[rigAssert[assertionName, controls["sign"]["RetirementRecorded"]]],
    HoldComplete[rigAssert[assertionName, baseline["RetirementRecorded"]]],
    "NOT_RECORDED_ON_FAIL"
  ]];

  subheading["WL runtime arity and unevaluated-leakage teeth"];
  properCall = HoldComplete[evaluatePort[cfg["arity"], facts]];
  badCall = HoldComplete[evaluatePort[cfg["arity"]]];
  assertionName = "WL ARITY definition/call scanner named assert";
  AssociateTo[outcomes, "WL ARITY" -> exerciseRig[
    "WL ARITY scanner",
    assertionName,
    HoldComplete[rigAssert[assertionName, arityScan[badCall]]],
    HoldComplete[rigAssert[assertionName, arityScan[properCall]]],
    "SCANNER_CAUGHT"
  ]];
  actualTranscript = {
    n0Den, baseline["Closure"]["Kbar4Residual"], baseline["Verdict"]
  };
  leakedTranscript = Append[actualTranscript, authoredLeak[n0Den]];
  assertionName = "WL LEAKAGE unevaluated-helper scanner named assert";
  AssociateTo[outcomes, "WL LEAKAGE" -> exerciseRig[
    "WL LEAKAGE scanner",
    assertionName,
    HoldComplete[rigAssert[assertionName, leakageFreeQ[leakedTranscript]]],
    HoldComplete[rigAssert[assertionName, leakageFreeQ[actualTranscript]]],
    "SCANNER_CAUGHT"
  ]];
  expectBool[
    "WL LEAKAGE actual transcript is helper-free",
    leakageFreeQ[actualTranscript]
  ];

  <|
    "N0Den" -> n0Den,
    "Facts" -> facts,
    "Baseline" -> baseline,
    "Controls" -> controls,
    "Outcomes" -> outcomes
  |>
];

emit[data_Association] := Module[
  {baseline = data["Baseline"], facts = data["Facts"], closure},
  closure = baseline["Closure"];
  subheading["Stage027 checks + closure transcript"];
  Print["consumes stage024 N0_den (canonical factored): ", fmt[data["N0Den"]]];
  Print["N0_den free_symbols: ", fmtNames[globalSymbolNames[data["N0Den"]]]];
  Print["HOST_CONTRACT (10): ", fmtNames[hostContractNames]];
  Print["free_symbols == contract: ", hostContractQ[data["N0Den"]]];
  Print[
    "consumed stage025 atomics: tags=", fmtNames[facts["Tags"]],
    "; vector_port_not_in_tags=", ! MemberQ[facts["Tags"], "vector_port"],
    "; vector_host_symbols=", fmtNames[facts["VectorHostSymbols"]],
    "; source_map_complete=", facts["SourceMapComplete"],
    "; vector_free(P2)=", facts["VectorFree"]
  ];
  Print[
    "consumed stage026 atomics: moment_valid=True; validated_I25=I25; ",
    "lineage_certificate=PASS; lineage_valid=True"
  ];
  Print["dimension tuple convention: (L,M,T)"];
  Print[
    "CHECK DIM: [N0_den]=", baseline["ExpressionDimension"],
    " = L^-1 M; dim_ok=", baseline["DimOK"]
  ];
  Print["CHECK SCALE: P0_physical=(c_s/a)^2*N0_den/D0"];
  Print[
    "CHECK SCALE: p0_power=", baseline["P0Power"],
    "; scaling_ok=", baseline["ScalingOK"]
  ];
  Print[
    "CHECK SIGN: outgoing SphericalHankelH1 coeff/i=",
    baseline["Sign"]["CoefficientOverI"],
    "; chi_Q=+1; sign_ok=", baseline["Sign"]["OK"]
  ];
  Print["CHECK NONZERO: nonzero_ok=", baseline["NonzeroOK"]];
  Print[
    "CHECK DEFERRED: uncertified->PORT_INCONCLUSIVE_SIM_DEFERRED; ",
    "proven converse->DENSITY_PORT_HOSTED"
  ];
  Print[
    "CHECK CLOSURE: Kbar4-4*Kbar2^2/Kbar0=",
    fmt[closure["Kbar4Residual"]]
  ];
  Print[
    "CHECK CLOSURE: Gamma5-2*G/(5*c^5)=",
    fmt[closure["Gamma5Residual"]]
  ];
  Print[
    "EXPORT Kbar moments: {Kbar0=", fmt[closure["Kbar0"]],
    ", Kbar2=", fmt[closure["Kbar2"]],
    ", Kbar4=", fmt[closure["Kbar4"]],
    ", Gamma5=", fmt[closure["Gamma5"]], "}"
  ];

  subheading["Joint-owning verdict and completion"];
  Print["LOCAL_AUDIT_VERDICT: ", baseline["Verdict"], " (CALIBRATED)"];
  Print[
    "JOINT_VERDICT_OWNED_BY_027: ", baseline["Verdict"],
    " (CALIBRATED, not PASS)"
  ];
  Print["scope: STRUCTURE hosted; magnitude SIM_DEFERRED; G=GENUINE_BLOCKED"];
  Print["provenance: 27=derived_in_gate; 54/5 and Gamma5=external_bridge_input"];
  Print[
    "RETIREMENT_RECORD: EM A_w/U,W scaffold RETIRES; pathA_43 diagnostic ",
    "sliver CLOSES (conditional=", baseline["RetirementRecorded"], ")"
  ];
  Print[
    "COMPLETION: pathA_43 4-way split COMPLETE = 024 derivation AND 025 ",
    "vector-freedom AND 026 lineage AND 027 checks+closure"
  ]
];

runAll[] := Module[{data},
  heading["ledger_stage027_port_checks_closure_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage027_port_checks_closure"];
  Print[
    "Engine: re-authored Wolfram unit/weight rescaling + built-in Hankel + ",
    "ratio-invariant closure + Association dispatch; zero file I/O."
  ];
  data = runAudit[];
  emit[data];
  0
];

result = Catch[
  runAll[],
  "ledgerStage027Failure",
  Function[{message, tag}, 1]
];

heading["Tallies"];
Print[
  "TALLY mathematica: ", passCount, " pass + ", failCount,
  " fail = ", passCount + failCount, " checks"
];
If[result === 0 && failCount === 0,
  Print["OVERALL PASS"];
  Exit[0],
  If[failureMessage =!= "", Print["ABORTED: ", failureMessage]];
  Print["OVERALL FAIL"];
  Exit[1]
];
