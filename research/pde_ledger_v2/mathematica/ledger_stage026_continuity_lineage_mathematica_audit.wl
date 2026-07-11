(* Ledger stage026 Wolfram audit: continuity-lineage token check.

   Standalone, print-only, no arguments, and no file I/O.  This is a native
   re-author: exact identifier tokens feed an Association/Graph ancestry walk,
   and the terminal l2 node produces I25 only after the full walk succeeds.
   Stage024's factored N0_den is cited, never reconstructed from an inverse.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

localVerdict = "CONTINUITY_LINEAGE_EARNED";
jointPartial = "DENSITY_PORT_HOSTED (3/4, CONTINUITY LINEAGE — I25 is a genuine ∫Y₂*·S_leak ℓ=2 continuity moment descended from pathA_29's operator; earns moment_valid; 024 = derivation, 025 = vector-freedom, 027 = port checks + closure)";
failVerdict = "FAIL_NOT_DENSITY_DERIVED";

continuityOperatorID = "pathA_29_projected_continuity";
continuityL0 = "M0 = Integral(S_leak d3x)";
continuityL1 = "D1_i = Integral(x_i*S_leak d3x) + Integral(j_i d3x)";
continuityL2 = "Q2_m = Integral(Y2_m_star*S_leak d3x)";
continuityL2Kernel = "Y2_m_star*S_leak";
badContinuityL2 = "GARBAGE_NOT_A_CONTINUITY_MOMENT_AT_ALL";
stuffedL0 = "NOT_M0 = FakeIntegral(S_leakage d3xyz)";

requiredTokens = <|
  "l0" -> {"M0", "Integral", "S_leak", "d3x"},
  "l1" -> {"D1_i", "Integral", "x_i", "S_leak", "j_i", "d3x"},
  "l2" -> {"Q2_m", "Integral", "Y2_m_star", "S_leak", "d3x"},
  "l2Kernel" -> {"Y2_m_star", "S_leak"}
|>;

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

clean[expr_] := Factor[Cancel[Together[expr]]];

externalSymbolNames = <|
  "cS" -> "c_s", "rhoEff" -> "rho_eff", "xiQ" -> "Xi_Q",
  "etaQ" -> "eta_q", "etaPhi" -> "eta_phi",
  "varpiQ2" -> "varpi_q2", "varpiPhi2" -> "varpi_Phi2",
  "lambdaC" -> "lambda_c", "foreignSubject" -> "foreign_subject",
  "OmegaU" -> "Omega_U", "OmegaW" -> "Omega_W",
  "Rmix" -> "R_mix", "gU" -> "g_U", "gW" -> "g_W",
  "Iwrong2" -> "I_wrong2"
|>;

externalName[symbol_Symbol] := Module[{raw = SymbolName[symbol]},
  If[KeyExistsQ[externalSymbolNames, raw], externalSymbolNames[raw], raw]
];

fmt[expr_String] := expr;
fmt[expr_] := StringReplace[
  ToString[InputForm[clean[expr]]],
  Normal[externalSymbolNames]
];

fmtNames[items_List] := "{" <> StringRiffle[Sort[items], ", "] <> "}";

fail[] := Throw[failureMessage, "ledgerStage026Failure"];

expectZero[name_, residual_] := Module[{c},
  c = clean[residual];
  If[TrueQ[c === 0],
    passCount++;
    Print["PASS  ", name],
    failCount++;
    failureMessage = name <> ": residual = " <> fmt[c];
    Print["FAIL  ", failureMessage];
    fail[]
  ]
];

expectBool[name_, condition_] := expectZero[name, If[TrueQ[condition], 0, 1]];

hostContractNames = Sort[{
  "a", "c_s", "rho_eff", "I25", "Xi_Q", "eta_q", "eta_phi",
  "varpi_q2", "varpi_Phi2", "lambda_c"
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

(* Exact [A-Za-z0-9_]+ extraction.  WordCharacter alone excludes underscore. *)
identifierTokens[text_] := DeleteDuplicates[
  StringCases[ToString[text], (WordCharacter | "_") ..]
];

tokenSetContainsQ[text_, needed_List] :=
  TrueQ[SubsetQ[identifierTokens[text], needed]];

rawSubstringContainsQ[text_, needed_List] := And @@
  (StringContainsQ[ToString[text], #] & /@ needed);

baselineLineage[] := <|
  "operator_id" -> continuityOperatorID,
  "moments" -> <|
    "l0" -> continuityL0,
    "l1" -> continuityL1,
    "l2" -> continuityL2
  |>,
  "l2_kernel" -> continuityL2Kernel,
  "lineage" -> "l0->l1->l2 from the same projected-continuity operator"
|>;

withDecoys[lineage_Association] := Join[
  lineage,
  <|"valid" -> True, "continuity_interface" -> True|>
];

replaceMoment[lineage_Association, level_String, value_String] := Module[
  {moments},
  moments = Association[Normal[lineage["moments"]]];
  AssociateTo[moments, level -> value];
  Join[lineage, <|"moments" -> moments|>]
];

(*
   Native decisive G1 route.  The graph defines the only allowed ancestry
   path.  Each path node consumes exact-token subsets; the l2 node also
   validates its kernel.  The walk's terminal state, not a flat copied
   conjunction or a self-reported flag, is verdict-bearing.
*)
ancestryWalk[lineage_Association] := Module[
  {moments, payload, graph, paths, path, nodeChecks, states, valid, tokens},
  moments = Lookup[lineage, "moments", <||>];
  payload = <|
    "l0" -> Lookup[moments, "l0", ""],
    "l1" -> Lookup[moments, "l1", ""],
    "l2" -> <|
      "moment" -> Lookup[moments, "l2", ""],
      "kernel" -> Lookup[lineage, "l2_kernel", ""]
    |>
  |>;
  graph = Graph[{
    DirectedEdge["l0", "l1"],
    DirectedEdge["l1", "l2"]
  }];
  paths = FindPath[graph, "l0", "l2"];
  path = If[paths === {}, {}, First[paths]];
  nodeChecks = <|
    "l0" -> tokenSetContainsQ[payload["l0"], requiredTokens["l0"]],
    "l1" -> tokenSetContainsQ[payload["l1"], requiredTokens["l1"]],
    "l2" -> (
      tokenSetContainsQ[payload["l2"]["moment"], requiredTokens["l2"]] &&
      tokenSetContainsQ[
        payload["l2"]["kernel"], requiredTokens["l2Kernel"]
      ]
    )
  |>;
  states = FoldList[
    Function[{state, level}, TrueQ[state && Lookup[nodeChecks, level, False]]],
    TrueQ[Lookup[lineage, "operator_id", ""] === continuityOperatorID],
    path
  ];
  valid = TrueQ[path === {"l0", "l1", "l2"} && Last[states]];
  tokens = <|
    "l0" -> identifierTokens[payload["l0"]],
    "l1" -> identifierTokens[payload["l1"]],
    "l2" -> identifierTokens[payload["l2"]["moment"]],
    "l2_kernel" -> identifierTokens[payload["l2"]["kernel"]]
  |>;
  <|
    "Valid" -> valid,
    "Path" -> path,
    "NodeChecks" -> nodeChecks,
    "States" -> states,
    "Tokens" -> tokens
  |>
];

flagOnlyWalk[lineage_Association] := <|
  "Valid" -> TrueQ[Lookup[lineage, "valid", False]],
  "Path" -> {"l0", "l1", "l2"},
  "NodeChecks" -> <||>,
  "States" -> {},
  "Tokens" -> <||>
|>;

(* G2 is produced at the terminal l2 node reached by the selected walker. *)
ancestryCertificate[
  lineage_Association,
  aPower_,
  walker_
] := Module[{walk, terminalValid, terminalMoment, terminalNode},
  walk = walker[lineage];
  terminalValid = TrueQ[walk["Valid"]];
  terminalMoment = Which[
    ! terminalValid, Iwrong2,
    TrueQ[aPower === -7/2], I25,
    True, Iwrong2
  ];
  terminalNode = <|
    "Level" -> "l2",
    "ProducedMoment" -> terminalMoment,
    "MomentValid" -> terminalValid
  |>;
  Join[walk, <|
    "EarnedMoment" -> terminalNode["ProducedMoment"],
    "MomentValid" -> terminalNode["MomentValid"],
    "TerminalNode" -> terminalNode
  |>]
];

continuityDependencyQ[
  lineageValid_,
  momentValid_,
  earnedMoment_,
  expr_,
  couplingZero_,
  secondArmEnabled_
] := Module[{membershipArm, vanishedArm},
  membershipArm = MemberQ[globalSymbols[expr], earnedMoment];
  vanishedArm = TrueQ[
    secondArmEnabled && clean[expr] === 0 && couplingZero
  ];
  TrueQ[lineageValid && momentValid && (membershipArm || vanishedArm)]
];

evaluateLineage[lineage_Association, expr_, walker_] := Module[
  {certificate, dependency},
  certificate = ancestryCertificate[lineage, -7/2, walker];
  dependency = continuityDependencyQ[
    certificate["Valid"],
    certificate["MomentValid"],
    certificate["EarnedMoment"],
    expr,
    False,
    True
  ];
  Join[certificate, <|
    "ContinuityDependencyOK" -> dependency,
    "Accepted" -> TrueQ[
      certificate["Valid"] &&
      certificate["EarnedMoment"] === I25 &&
      certificate["MomentValid"] && dependency
    ]
  |>]
];

negativeLineages[] := Module[
  {base, fake, attack2, operatorDrop, l0Drop, l1Drop, kernelDrop, stuffed},
  base = baselineLineage[];
  fake = <|
    "operator_id" -> "mis_tagged_vector_formula",
    "moments" -> <|
      "l0" -> "Omega_U",
      "l1" -> "Omega_W",
      "l2" -> "Omega_U^2*g_W + R_mix*g_U, relabeled as continuity"
    |>,
    "l2_kernel" -> "R_mix*g_U"
  |>;
  attack2 = replaceMoment[base, "l2", badContinuityL2];
  operatorDrop = Join[
    base,
    <|"operator_id" -> "pathA_29_projected_continuity_forged"|>
  ];
  l0Drop = replaceMoment[base, "l0", "mass = Integral(S_leak d3x)"];
  l1Drop = replaceMoment[
    base,
    "l1",
    "dipole = Integral(x_i*S_leak d3x) + Integral(j_i d3x)"
  ];
  kernelDrop = Join[base, <|"l2_kernel" -> "angular_kernel*S_leak"|>];
  stuffed = replaceMoment[base, "l0", stuffedL0];
  <|
    "A fake_continuity" -> withDecoys[fake],
    "B attack2" -> withDecoys[attack2],
    "C operator_id" -> withDecoys[operatorDrop],
    "C l0_token" -> withDecoys[l0Drop],
    "C l1_token" -> withDecoys[l1Drop],
    "C l2_kernel_token" -> withDecoys[kernelDrop],
    "C token_stuffing" -> withDecoys[stuffed]
  |>
];

rigAssert[name_String, condition_] := If[
  ! TrueQ[condition],
  Throw[name, "stage026RigAssertion"]
];

probeAssertion[held_HoldComplete, expectedName_String] := Module[{fired},
  fired = Catch[
    ReleaseHold[held];
    "NO_ASSERT_FIRED",
    "stage026RigAssertion",
    Function[{message, tag}, message]
  ];
  {TrueQ[fired === expectedName], fired}
];

assertionPasses[held_HoldComplete] := TrueQ[Catch[
  ReleaseHold[held];
  True,
  "stage026RigAssertion",
  Function[{message, tag}, False]
]];

exerciseRig[
  label_String,
  assertionName_String,
  rigHeld_HoldComplete,
  neutralHeld_HoldComplete,
  verdict_String : failVerdict
] := Module[{probe, caught, firedName, neutralPass, outcome},
  probe = probeAssertion[rigHeld, assertionName];
  caught = probe[[1]];
  firedName = probe[[2]];
  neutralPass = assertionPasses[neutralHeld];
  expectBool[
    "META " <> label <> " routed assertion fires and neutering stops it",
    caught && neutralPass
  ];
  outcome = verdict <> " at " <> ToString[firedName] <> "; neutralized=PASS";
  Print["RIG ", label, ": ", outcome];
  outcome
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
  definitionArity[ancestryCertificate] === heldCallArity[held];

leakageFreeQ[objects_] := FreeQ[
  objects,
  identifierTokens | tokenSetContainsQ | ancestryWalk |
  ancestryCertificate | continuityDependencyQ | evaluateLineage |
  authoredLeak | SubsetQ | StringCases | FindPath | FullSimplify
];

lineageString[lineage_Association] := Module[{moments = lineage["moments"]},
  "{operator_id: " <> lineage["operator_id"] <>
  ", l0: " <> moments["l0"] <>
  ", l1: " <> moments["l1"] <>
  ", l2: " <> moments["l2"] <>
  ", l2_kernel: " <> lineage["l2_kernel"] <> "}"
];

runAudit[] := Module[
  {n0Den, liveNames, lineage, baseline, negatives, baselineTokens,
   stuffed, stuffedTokens, corruptedGate, outcomes, fakeVectorExpr,
   result, assertName, neuteredCertificate, actualCertificate, flagAccepts,
   exactRejects, flagAssert, flipLineage, flipResult, g1, g2, g3,
   g1Neutered, gAssert, noI25, noEtaQ, foreignRho, iMutations,
   properCall, badCall, actualTranscript, leakedTranscript,
   dependencyMembershipWitness, localOK},

  n0Den = clean[citedN0Den[]];
  liveNames = globalSymbolNames[n0Den];
  lineage = baselineLineage[];
  baseline = evaluateLineage[lineage, n0Den, ancestryWalk];
  negatives = negativeLineages[];

  subheading["I. cited stage024 consumption-integrity contract"];
  expectBool[
    "I cited N0_den free_symbols equals stage024 exact 10-symbol host contract",
    liveNames === hostContractNames
  ];

  subheading["G1/G2. exact-token ancestry and moment earning"];
  baselineTokens = baseline["Tokens"];
  expectBool[
    "G1 genuine CONTINUITY_L* identifier tokens survive intact",
    SubsetQ[baselineTokens["l0"], requiredTokens["l0"]] &&
    SubsetQ[baselineTokens["l1"], requiredTokens["l1"]] &&
    SubsetQ[baselineTokens["l2"], requiredTokens["l2"]] &&
    SubsetQ[baselineTokens["l2_kernel"], requiredTokens["l2Kernel"]]
  ];
  stuffed = negatives["C token_stuffing"];
  stuffedTokens = identifierTokens[stuffed["moments"]["l0"]];
  expectBool[
    "G1 token-stuffing passes raw substring but fails exact lexical tokens",
    rawSubstringContainsQ[stuffed["moments"]["l0"], requiredTokens["l0"]] &&
    ! tokenSetContainsQ[stuffed["moments"]["l0"], requiredTokens["l0"]]
  ];
  expectBool[
    "G1 baseline continuity_lineage_valid is computed True",
    baseline["Valid"]
  ];
  expectBool[
    "G2 baseline earns exactly (I25, moment_valid=True)",
    baseline["EarnedMoment"] === I25 && baseline["MomentValid"]
  ];
  corruptedGate = evaluateLineage[
    negatives["B attack2"], n0Den, ancestryWalk
  ];
  expectBool[
    "G2 corrupted lineage earns exactly (I_wrong2, moment_valid=False)",
    corruptedGate["EarnedMoment"] === Iwrong2 &&
    ! corruptedGate["MomentValid"]
  ];
  Print[
    "DE-COUNTED diagnostic: baseline continuity_dependency_ok is True = ",
    baseline["ContinuityDependencyOK"]
  ];

  subheading["A-D/F. lineage rigs and coupling meta-tests"];
  outcomes = <||>;
  fakeVectorExpr = clean[OmegaU^2 gW + Rmix gU];
  KeyValueMap[
    Function[{label, badLineage},
      result = evaluateLineage[
        badLineage,
        If[label === "A fake_continuity", fakeVectorExpr, n0Den],
        ancestryWalk
      ];
      assertName = label <> " G1 exact-token lineage-gate assert";
      AssociateTo[outcomes, label -> exerciseRig[
        label,
        assertName,
        HoldComplete[rigAssert[assertName, result["Valid"]]],
        HoldComplete[rigAssert[assertName, baseline["Valid"]]]
      ]];
      Print[
        "RIG DATA ", label,
        ": continuity_lineage_valid=", result["Valid"],
        "; earned_moment=", externalName[result["EarnedMoment"]],
        "; moment_valid=", result["MomentValid"],
        "; decoy_valid=True; self_asserted_continuity_interface=True"
      ]
    ],
    negatives
  ];

  neuteredCertificate = ancestryCertificate[
    negatives["B attack2"], -7/2, ancestryWalk
  ];
  neuteredCertificate = Join[
    neuteredCertificate,
    <|"EarnedMoment" -> I25|>
  ];
  actualCertificate = ancestryCertificate[
    negatives["B attack2"], -7/2, ancestryWalk
  ];
  assertName = "D earning gate distinguishes I25 from I_wrong2 assert";
  AssociateTo[outcomes, "D earning_gate" -> exerciseRig[
    "D earning_gate",
    assertName,
    HoldComplete[rigAssert[
      assertName,
      neuteredCertificate["EarnedMoment"] === Iwrong2 &&
      ! neuteredCertificate["MomentValid"]
    ]],
    HoldComplete[rigAssert[
      assertName,
      actualCertificate["EarnedMoment"] === Iwrong2 &&
      ! actualCertificate["MomentValid"]
    ]]
  ]];

  flagAccepts = And @@ (
    ancestryCertificate[#, -7/2, flagOnlyWalk]["Valid"] & /@ Values[negatives]
  );
  exactRejects = And @@ (
    ! ancestryCertificate[#, -7/2, ancestryWalk]["Valid"] & /@ Values[negatives]
  );
  flagAssert = "META flag-only validator rejects every decoy-negative assert";
  AssociateTo[outcomes, "META flag_only_validator" -> exerciseRig[
    "META flag_only_validator",
    flagAssert,
    HoldComplete[rigAssert[flagAssert, ! flagAccepts]],
    HoldComplete[rigAssert[flagAssert, exactRejects]]
  ]];

  flipLineage = withDecoys[replaceMoment[
    lineage,
    "l1",
    "D1_i = Integral(x_i*S_leak d3x)"
  ]];
  flipResult = evaluateLineage[flipLineage, n0Den, ancestryWalk];
  expectBool[
    "F baseline-valid positive and ancestry-corruption flip assert",
    baseline["Accepted"] && ! flipResult["Accepted"]
  ];
  AssociateTo[
    outcomes,
    "F baseline_positive" ->
      "PASS; corrupt-any-token=FLIP to FAIL_NOT_DENSITY_DERIVED"
  ];
  Print["RIG F baseline_positive: ", outcomes["F baseline_positive"]];

  subheading["G. isolated vanished-coupling OR-arm probes"];
  g1 = continuityDependencyQ[True, True, I25, 0, True, True];
  g2 = continuityDependencyQ[True, True, I25, 0, False, True];
  g3 = continuityDependencyQ[True, True, I25, foreignSubject, True, True];
  expectBool["G(i) expr=0 and coupling_zero=True passes via arm 2", g1];
  expectBool["G(ii) expr=0 and coupling_zero=False fails", ! g2];
  expectBool[
    "G(iii) nonzero expr missing earned moment fails despite coupling_zero=True",
    ! g3
  ];
  gAssert = "G(iv) neutered vanished-coupling arm makes probe (i) fail assert";
  g1Neutered = continuityDependencyQ[True, True, I25, 0, True, False];
  AssociateTo[outcomes, "G arm2_meta" -> exerciseRig[
    "G arm2_meta",
    gAssert,
    HoldComplete[rigAssert[gAssert, g1Neutered]],
    HoldComplete[rigAssert[gAssert, g1]]
  ]];

  subheading["I. consumption-integrity mutations"];
  noI25 = clean[n0Den/I25^2];
  noEtaQ = clean[
    I25^2 xiQ^2 cS^4 rhoEff (etaPhi varpiQ2)^2/
      (a^7 (lambdaC^2 - varpiPhi2 varpiQ2)^2)
  ];
  foreignRho = clean[n0Den /. rhoEff -> foreignSubject];
  iMutations = <|
    "I drop_external_I25_squared" -> noI25,
    "I drop_eta_q_lambda_c_term" -> noEtaQ,
    "I replace_rho_eff" -> foreignRho
  |>;
  KeyValueMap[
    Function[{label, subject},
      assertName = label <> " exact host-contract assert";
      AssociateTo[outcomes, label -> exerciseRig[
        label,
        assertName,
        HoldComplete[rigAssert[
          assertName,
          globalSymbolNames[subject] === hostContractNames
        ]],
        HoldComplete[rigAssert[
          assertName,
          liveNames === hostContractNames
        ]],
        "CONTRACT_REJECTED"
      ]]
    ],
    iMutations
  ];

  subheading["H'. runtime arity and unevaluated-leakage scanners"];
  properCall = HoldComplete[
    ancestryCertificate[lineage, -7/2, ancestryWalk]
  ];
  badCall = HoldComplete[ancestryCertificate[lineage, -7/2]];
  AssociateTo[outcomes, "H' arity_scanner" -> exerciseRig[
    "H' arity_scanner",
    "H' definition/call arity scanner assert",
    HoldComplete[rigAssert[
      "H' definition/call arity scanner assert",
      arityScan[badCall]
    ]],
    HoldComplete[rigAssert[
      "H' definition/call arity scanner assert",
      arityScan[properCall]
    ]],
    "SCANNER_CAUGHT"
  ]];
  actualTranscript = {n0Den, baseline, baselineTokens};
  leakedTranscript = Append[actualTranscript, authoredLeak[n0Den]];
  AssociateTo[outcomes, "H' leakage_scanner" -> exerciseRig[
    "H' leakage_scanner",
    "H' unevaluated-leakage transcript scanner assert",
    HoldComplete[rigAssert[
      "H' unevaluated-leakage transcript scanner assert",
      leakageFreeQ[leakedTranscript]
    ]],
    HoldComplete[rigAssert[
      "H' unevaluated-leakage transcript scanner assert",
      leakageFreeQ[actualTranscript]
    ]],
    "SCANNER_CAUGHT"
  ]];
  expectBool[
    "H' actual transcript has no unevaluated authored-helper leakage",
    leakageFreeQ[actualTranscript]
  ];

  dependencyMembershipWitness = MemberQ[
    globalSymbols[n0Den], baseline["EarnedMoment"]
  ];
  Print[
    "E DE-COUNTED witness: earned_moment in N0_den.free_symbols = ",
    dependencyMembershipWitness,
    " (printed, not tallied)"
  ];

  localOK = TrueQ[
    liveNames === hostContractNames &&
    baseline["Valid"] && baseline["EarnedMoment"] === I25 &&
    baseline["MomentValid"] && baseline["ContinuityDependencyOK"]
  ];
  expectBool["LOCAL CONTINUITY_LINEAGE_EARNED = G1 and G2", localOK];

  <|
    "N0Den" -> n0Den,
    "LiveNames" -> liveNames,
    "Lineage" -> lineage,
    "Baseline" -> baseline,
    "TokenSets" -> baselineTokens,
    "StuffedTokens" -> stuffedTokens,
    "Outcomes" -> outcomes,
    "DependencyMembershipWitness" -> dependencyMembershipWitness,
    "LocalOK" -> localOK
  |>
];

emit[data_Association] := Module[{baseline = data["Baseline"], tokens},
  tokens = data["TokenSets"];
  subheading["Computed continuity-lineage transcript"];
  Print["consumes: stage024 N0_den (cited canonical factored export)"];
  Print["N0_den (canonical factored): ", fmt[data["N0Den"]]];
  Print["N0_den free_symbols: ", fmtNames[data["LiveNames"]]];
  Print["host contract (10): ", fmtNames[hostContractNames]];
  Print["free_symbols == contract: ", data["LiveNames"] === hostContractNames];
  Print["lineage dict: ", lineageString[data["Lineage"]]];
  Print["exact tokens l0: ", fmtNames[tokens["l0"]]];
  Print["exact tokens l1: ", fmtNames[tokens["l1"]]];
  Print["exact tokens l2: ", fmtNames[tokens["l2"]]];
  Print["exact tokens l2_kernel: ", fmtNames[tokens["l2_kernel"]]];
  Print["stuffed l0 exact tokens: ", fmtNames[data["StuffedTokens"]]];
  Print["tokenizer: StringCases[s, (WordCharacter | \"_\")..]"];
  Print[
    "continuity_lineage_valid=", baseline["Valid"],
    " (computed exact tokens; decoy flags ignored)"
  ];
  Print[
    "earned (earned_moment=", externalName[baseline["EarnedMoment"]],
    ", moment_valid=", baseline["MomentValid"], ")"
  ];
  Print[
    "continuity_dependency_ok=", baseline["ContinuityDependencyOK"]
  ];
  Print[
    "E dependency membership witness (de-counted): ",
    data["DependencyMembershipWitness"]
  ];
  Print[
    "exports: moment_valid=True; validated_I25=I25; ",
    "lineage_certificate=PASS"
  ];

  subheading["Verdict labels"];
  Print["LOCAL_AUDIT_VERDICT: ", localVerdict];
  Print["JOINT_LANDING_LABEL (PARTIAL): ", jointPartial]
];

runAll[] := Module[{data},
  heading["ledger_stage026_continuity_lineage_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage026_continuity_lineage"];
  Print[
    "Engine: Wolfram exact-token Graph ancestry traversal with terminal ",
    "l2 moment production; zero file I/O."
  ];
  data = runAudit[];
  emit[data];
  0
];

result = Catch[
  runAll[],
  "ledgerStage026Failure",
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
