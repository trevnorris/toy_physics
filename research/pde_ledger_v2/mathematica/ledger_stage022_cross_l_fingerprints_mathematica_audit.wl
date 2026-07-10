(* Ledger stage022 Mathematica audit: cross-ell DtN fingerprints.

   Standalone, print-only, exact, and zero-file-I/O.  This is a genuinely
   re-authored Wolfram route: built-in SphericalHankelH1/H2 objects and
   direct SeriesCoefficient extraction, not the source script's hand-built
   trigonometric branch-record translation.  Stage023 owns the joint FAIL departure.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
failureMessage = "";

CROSSLFINGERPRINTOK = "CROSS_L_FINGERPRINT_OK";
FAILFINGERPRINT = "FAIL_FINGERPRINT";
FAILQUADREGRESSION = "FAIL_QUAD_REGRESSION";
JOINTPARTIAL = "FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)";

$Assumptions = Element[z, Reals];

fail[msg_] := Throw[msg, "ledgerStage022Failure"];

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

clean[expr_] := FullSimplify[Cancel[Together[expr]], $Assumptions];
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

ells = {0, 1, 2};
expectedRadiative = <|0 -> 1, 1 -> 1/2, 2 -> 1/27|>;
expectedOrder = Association[Table[ell -> 2 ell + 1, {ell, ells}]];
stage018Typed = <|"u2_z" -> 1/9, "u4_z" -> 4/81, "v5_z" -> 1/27|>;

(* Independent special-function construction. *)
hankelBuiltin[ell_, "outgoing"] := SphericalHankelH1[ell, z];
hankelBuiltin[ell_, "incoming"] := SphericalHankelH2[ell, z];

copyFactor[ell_, "none"] := 1;
copyFactor[ell_, "radiative_slot_copy"] := 1 + I z^(2 ell + 1);
copyFactor[ell_, "lower_imaginary_copy"] := 1 + I z;
copyFactor[2, "u2_only_derivation_copy"] := 1 + z^2 + z^4/18;
copyFactor[2, "u4_only_derivation_copy"] := 1 + z^4;
copyFactor[2, "v5_only_derivation_copy"] := 1 + I z^5;

waveAt[ell_, sense_String, mutation_String] :=
  hankelBuiltin[ell, sense] copyFactor[ell, mutation];

lambdaAt[ell_, sense_String, mutation_String] :=
  lambdaAt[ell, sense, mutation] = Module[{wave},
    wave = waveAt[ell, sense, mutation];
    FullSimplify[z D[wave, z]/wave, $Assumptions]
  ];

responseAt[ell_, sense_String, mutation_String] :=
  responseAt[ell, sense, mutation] =
    FullSimplify[-(ell + 1)/lambdaAt[ell, sense, mutation], $Assumptions];

lambdaCoefficient[ell_, sense_String, mutation_String, power_Integer] :=
  lambdaCoefficient[ell, sense, mutation, power] = FullSimplify[
    SeriesCoefficient[lambdaAt[ell, sense, mutation], {z, 0, power}],
    $Assumptions
  ];

responseCoefficient[ell_, sense_String, mutation_String, power_Integer] :=
  responseCoefficient[ell, sense, mutation, power] = FullSimplify[
    SeriesCoefficient[responseAt[ell, sense, mutation], {z, 0, power}],
    $Assumptions
  ];

imaginaryPart[expr_] := FullSimplify[ComplexExpand[Im[expr]], $Assumptions];
realPart[expr_] := FullSimplify[ComplexExpand[Re[expr]], $Assumptions];

firstImaginaryOrder[ell_, sense_String, mutation_String] := Module[
  {powers, found},
  powers = Range[0, Max[7, 2 ell + 3]];
  found = SelectFirst[
    powers,
    ! TrueQ[clean[imaginaryPart[responseCoefficient[ell, sense, mutation, #]]] === 0] &,
    Missing["NotFound"]
  ];
  found
];

lowerImaginaryCoefficientsZero[ell_, sense_String, mutation_String] := And @@ Table[
  TrueQ[clean[imaginaryPart[responseCoefficient[ell, sense, mutation, power]]] === 0],
  {power, 0, 2 ell}
];

lowerRealCoefficientsAgree[ell_, outgoingMutation_String] := And @@ Table[
  TrueQ[
    clean[
      realPart[responseCoefficient[ell, "incoming", "none", power]] -
      realPart[responseCoefficient[ell, "outgoing", outgoingMutation, power]]
    ] === 0
  ],
  {power, 0, 2 ell}
];

poleLambdaStatic[ell_] := poleLambdaStatic[ell] = Module[{mutatedWave, mutatedLambda},
  mutatedWave = z SphericalHankelH1[ell, z];
  mutatedLambda = FullSimplify[z D[mutatedWave, z]/mutatedWave, $Assumptions];
  FullSimplify[SeriesCoefficient[mutatedLambda, {z, 0, 0}], $Assumptions]
];

fingerprintAudit[mutationEll_] := Module[
  {
    mutationFor, outRadiative, incomingRadiative, lambdaStatic, staticDiagnostic,
    scannedOrders, lowerImaginaryOk, incomingSignOk, incomingRealOk, matches
  },
  mutationFor[ell_] := If[IntegerQ[mutationEll] && ell === mutationEll,
    "radiative_slot_copy", "none"
  ];
  outRadiative = Association[Table[
    ell -> clean[responseCoefficient[ell, "outgoing", mutationFor[ell], 2 ell + 1]/I],
    {ell, ells}
  ]];
  incomingRadiative = Association[Table[
    ell -> clean[responseCoefficient[ell, "incoming", "none", 2 ell + 1]/I],
    {ell, ells}
  ]];
  lambdaStatic = Association[Table[
    ell -> lambdaCoefficient[ell, "outgoing", mutationFor[ell], 0],
    {ell, ells}
  ]];
  staticDiagnostic = Association[Table[
    ell -> responseCoefficient[ell, "outgoing", mutationFor[ell], 0],
    {ell, ells}
  ]];
  scannedOrders = Association[Table[
    ell -> firstImaginaryOrder[ell, "outgoing", mutationFor[ell]],
    {ell, ells}
  ]];
  lowerImaginaryOk = Association[Table[
    ell -> lowerImaginaryCoefficientsZero[ell, "outgoing", mutationFor[ell]],
    {ell, ells}
  ]];
  incomingSignOk = Association[Table[
    ell -> TrueQ[clean[incomingRadiative[ell] + outRadiative[ell]] === 0],
    {ell, ells}
  ]];
  incomingRealOk = Association[Table[
    ell -> lowerRealCoefficientsAgree[ell, mutationFor[ell]],
    {ell, ells}
  ]];
  matches = <|
    "RadiativeCoeff" -> Association[Table[
      ell -> TrueQ[clean[outRadiative[ell] - expectedRadiative[ell]] === 0],
      {ell, ells}
    ]],
    "LambdaStatic" -> Association[Table[
      ell -> TrueQ[clean[lambdaStatic[ell] + ell + 1] === 0],
      {ell, ells}
    ]],
    "ScannedOrder" -> Association[Table[
      ell -> TrueQ[scannedOrders[ell] === expectedOrder[ell]],
      {ell, ells}
    ]],
    "AllLowerImaginaryZero" -> lowerImaginaryOk,
    "IncomingFlipsSign" -> incomingSignOk,
    "IncomingLowerRealUnchanged" -> incomingRealOk
  |>;
  <|
    "OutRadiative" -> outRadiative,
    "IncomingRadiative" -> incomingRadiative,
    "LambdaStatic" -> lambdaStatic,
    "StaticDiagnostic" -> staticDiagnostic,
    "ScannedOrders" -> scannedOrders,
    "Matches" -> matches,
    "ChiQDiagnostic" -> clean[
      responseCoefficient[2, "outgoing", mutationFor[2], 5]/I/stage018Typed["v5_z"]
    ],
    "Ok" -> And @@ Flatten[Values /@ Values[matches]]
  |>
];

quadrupoleTuple[sense_String, mutation_String] := <|
  "u2_z" -> responseCoefficient[2, sense, mutation, 2],
  "u4_z" -> responseCoefficient[2, sense, mutation, 4],
  "v5_z" -> clean[responseCoefficient[2, sense, mutation, 5]/I]
|>;

gate4Audit[breakGate4_, mutation_String] := Module[{sense, tuple, matches},
  sense = If[TrueQ[breakGate4], "incoming", "outgoing"];
  tuple = quadrupoleTuple[sense, mutation];
  matches = Association[Table[
    slot -> TrueQ[clean[tuple[slot] - stage018Typed[slot]] === 0],
    {slot, Keys[stage018Typed]}
  ]];
  <|
    "BranchUsed" -> sense,
    "Mutation" -> mutation,
    "Derived" -> tuple,
    "Matches" -> matches,
    "ChiQDiagnostic" -> clean[tuple["v5_z"]/stage018Typed["v5_z"]],
    "Ok" -> And @@ Values[matches]
  |>
];

readLocalGate[name_String, value_] := (Sow[name, "022LocalRead"]; TrueQ[value]);

localAuditVerdict[fingerprintGate_, nonRegressionGate_] := Module[{packet, verdict, reads},
  packet = Reap[
    verdict = Which[
      ! readLocalGate["cross_l_fingerprints", fingerprintGate], FAILFINGERPRINT,
      ! readLocalGate["ell2_non_regression", nonRegressionGate], FAILQUADREGRESSION,
      True, CROSSLFINGERPRINTOK
    ],
    "022LocalRead"
  ];
  reads = If[Length[packet[[2]]] === 1, DeleteDuplicates[First[packet[[2]]]], {}];
  <|"Verdict" -> packet[[1]], "VerdictReadSet" -> reads|>
];

localContext[breakGate4_, mutationEll_] := Module[{fingerprints, gate4, local},
  fingerprints = fingerprintAudit[mutationEll];
  gate4 = gate4Audit[breakGate4, "none"];
  local = localAuditVerdict[fingerprints["Ok"], gate4["Ok"]];
  <|"Fingerprints" -> fingerprints, "Gate4" -> gate4, "Local" -> local|>
];

dynamicGate4Ablation[breakGate4_] := Module[
  {withContext, withoutContext, withMutation, withoutMutation, trace},
  withContext = localContext[breakGate4, None];
  withoutContext = localContext[False, None];
  withMutation = withContext["Local"]["Verdict"];
  withoutMutation = withoutContext["Local"]["Verdict"];
  trace = {
    {withContext["Gate4"]["BranchUsed"], withMutation},
    {withoutContext["Gate4"]["BranchUsed"], withoutMutation}
  };
  <|
    "RerunGateLogic" -> TrueQ[trace[[1, 2]] =!= trace[[2, 2]]],
    "VerdictTrace" -> trace,
    "WithMutation" -> withMutation,
    "WithoutMutation" -> withoutMutation,
    "ExpectedFail" -> FAILQUADREGRESSION,
    "FailSuppressed" -> TrueQ[
      withMutation === FAILQUADREGRESSION &&
      withoutMutation === CROSSLFINGERPRINTOK &&
      withMutation =!= withoutMutation
    ]
  |>
];

baselineContext = localContext[False, None];
baselineFingerprints = baselineContext["Fingerprints"];
baselineGate4 = baselineContext["Gate4"];
baselineLocal = baselineContext["Local"];

runCrossEllDerivation[] := Module[{ell, power, outgoing, incoming},
  subheading["Cross-ell outgoing DtN fingerprints (computed in z-space)"];
  Print["  Wolfram route: built-in SphericalHankelH1/H2 with direct SeriesCoefficient extraction."];
  Print["  Lambda_ell=z*h_ell'/h_ell; Yhat_ell=-(ell+1)/Lambda_ell; no hand-built trigonometric branch record is used."];
  Do[
    power = expectedOrder[ell];
    outgoing = baselineFingerprints["OutRadiative"][ell];
    incoming = baselineFingerprints["IncomingRadiative"][ell];
    Print[
      "  ell=", ell, ": radiative_coeff = ", fmt[outgoing],
      " at omega^", power, " (z^", power, "); Lambda_static = ",
      fmt[baselineFingerprints["LambdaStatic"][ell]],
      "; static = ", fmt[baselineFingerprints["StaticDiagnostic"][ell]],
      " DIAGNOSTIC; first_nonzero_imag_order = ",
      baselineFingerprints["ScannedOrders"][ell],
      " (SCANNED); incoming-flips-sign = ",
      baselineFingerprints["Matches"]["IncomingFlipsSign"][ell]
    ];
    expectZero[
      "ell=" <> ToString[ell] <> " derived radiative_coeff matches typed cross-ell target",
      outgoing - expectedRadiative[ell]
    ];
    expectZero[
      "ell=" <> ToString[ell] <> " Lambda_static derived from Hankel log-derivative is -(ell+1)",
      baselineFingerprints["LambdaStatic"][ell] + ell + 1
    ];
    expectBool[
      "ell=" <> ToString[ell] <> " scanned first nonzero imaginary power is 2*ell+1",
      baselineFingerprints["ScannedOrders"][ell] === power
    ];
    expectBool[
      "ell=" <> ToString[ell] <> " every lower imaginary coefficient vanishes",
      baselineFingerprints["Matches"]["AllLowerImaginaryZero"][ell]
    ];
    expectZero[
      "ell=" <> ToString[ell] <> " incoming radiative coefficient flips only the sign",
      incoming + outgoing
    ];
    expectBool[
      "ell=" <> ToString[ell] <> " incoming lower real coefficients equal outgoing",
      baselineFingerprints["Matches"]["IncomingLowerRealUnchanged"][ell]
    ],
    {ell, ells}
  ];
  Print["  static=1 is a de-counted diagnostic derived from -(ell+1)/Lambda_static; it is not a verdict tooth."];
  Print[
    "  chi_Q = 27*v5_z = ", fmt[baselineFingerprints["ChiQDiagnostic"]],
    " DIAGNOSTIC (de-counted because it is algebraically subsumed by v5_z=1/27)."
  ]
];

runGate4NonRegression[] := Module[{slot},
  subheading["Gate-4 ell=2 non-regression against stage018 typed fingerprint"];
  Print["  CHECKABLE consumption: stage018 independently earned the typed z-space tuple."];
  Print[
    "  ell=2 derived tuple: u2_z=", fmt[baselineGate4["Derived"]["u2_z"]],
    ", u4_z=", fmt[baselineGate4["Derived"]["u4_z"]],
    ", v5_z=", fmt[baselineGate4["Derived"]["v5_z"]]
  ];
  Do[
    expectZero[
      "ell=2 derived " <> slot <> " matches stage018 earned typed " <> fmt[stage018Typed[slot]],
      baselineGate4["Derived"][slot] - stage018Typed[slot]
    ],
    {slot, Keys[stage018Typed]}
  ];
  Print[
    "  chi_Q = ", fmt[baselineGate4["ChiQDiagnostic"]],
    " DIAGNOSTIC only; the v5_z equality owns this content."
  ]
];

runPerToothAblations[] := Module[
  {ell, mutatedCoefficient, claimedOutgoing, lowerOrder, slotMutations, slot, mutation, mutant, other},
  subheading["Per-tooth derivation-copy ablations"];
  Do[
    mutatedCoefficient = clean[
      responseCoefficient[ell, "outgoing", "radiative_slot_copy", 2 ell + 1]/I
    ];
    expectFail[
      "ell=" <> ToString[ell] <> " radiative derivation copy fires its own coefficient assert",
      mutatedCoefficient - expectedRadiative[ell]
    ];
    expectFail[
      "ell=" <> ToString[ell] <> " pole-order h_mut=z*h fires the derived Lambda_static assert",
      poleLambdaStatic[ell] + ell + 1
    ];
    expectZero[
      "ell=" <> ToString[ell] <> " H1-to-H2 is inert for Lambda_static (not its mutant)",
      lambdaCoefficient[ell, "incoming", "none", 0] -
      lambdaCoefficient[ell, "outgoing", "none", 0]
    ];
    claimedOutgoing = baselineFingerprints["IncomingRadiative"][ell];
    expectFail[
      "ell=" <> ToString[ell] <> " outgoing-to-incoming mutation fires sign-flip tooth",
      claimedOutgoing + baselineFingerprints["IncomingRadiative"][ell]
    ],
    {ell, ells}
  ];
  lowerOrder = firstImaginaryOrder[2, "outgoing", "lower_imaginary_copy"];
  Print["  +I*z lower-order Hankel-copy mutation: first_nonzero_imag_order=", lowerOrder];
  expectFail[
    "ell=2 +I*z corruption fires the SCANNED radiative-order assert",
    lowerOrder - expectedOrder[2]
  ];
  expectBool[
    "ell=2 +I*z corruption makes a lower imaginary coefficient nonzero",
    ! lowerImaginaryCoefficientsZero[2, "outgoing", "lower_imaginary_copy"]
  ];
  slotMutations = <|
    "u2_z" -> "u2_only_derivation_copy",
    "u4_z" -> "u4_only_derivation_copy",
    "v5_z" -> "v5_only_derivation_copy"
  |>;
  Do[
    mutation = slotMutations[slot];
    mutant = gate4Audit[False, mutation];
    Print[
      "  ", slot, "-only copy tuple = {u2_z:", fmt[mutant["Derived"]["u2_z"]],
      ", u4_z:", fmt[mutant["Derived"]["u4_z"]],
      ", v5_z:", fmt[mutant["Derived"]["v5_z"]], "}"
    ];
    expectFail[
      "ell=2 " <> slot <> "-only derivation copy fires its own non-regression assert",
      mutant["Derived"][slot] - stage018Typed[slot]
    ];
    Do[
      If[other =!= slot,
        expectZero[
          "ell=2 " <> slot <> "-only copy leaves " <> other <> " unchanged",
          mutant["Derived"][other] - stage018Typed[other]
        ]
      ],
      {other, Keys[stage018Typed]}
    ],
    {slot, Keys[slotMutations]}
  ]
];

runLocalVerdictAnd3e[] := Module[
  {fingerprintMutant, broken, ablation, neutered, expectedReads, forbiddenReads},
  subheading["022-local verdict, read-set, and 3e dynamic self-ablation"];
  fingerprintMutant = localContext[False, 1];
  broken = localContext[True, None];
  ablation = dynamicGate4Ablation[True];
  neutered = dynamicGate4Ablation[False];
  expectedReads = {"cross_l_fingerprints", "ell2_non_regression"};
  forbiddenReads = {"nullspace", "return_admittance", "base_verdict", "residuals", "selector"};
  Print[
    "  baseline gate values = ",
    <|"cross_l_fingerprints" -> baselineFingerprints["Ok"],
      "ell2_non_regression" -> baselineGate4["Ok"]|>
  ];
  Print["  computed LOCAL_AUDIT_VERDICT read-set = ", Sort[baselineLocal["VerdictReadSet"]]];
  expectBool[
    "LOCAL_AUDIT_VERDICT reads exactly fingerprints plus ell=2 non-regression",
    Sort[baselineLocal["VerdictReadSet"]] === Sort[expectedReads]
  ];
  Print[
    "  LOCAL_AUDIT_VERDICT read-set excludes nullspace/return/base_verdict objects = ",
    Intersection[baselineLocal["VerdictReadSet"], forbiddenReads] === {},
    " DIAGNOSTIC (de-counted; the exact read-set equality tooth owns this content)."
  ];
  expectZero[
    "baseline 022-local verdict is CROSS_L_FINGERPRINT_OK",
    verdictResidual[baselineLocal["Verdict"], CROSSLFINGERPRINTOK]
  ];
  expectZero[
    "corrupting an ell=1 fingerprint derivation reaches FAIL_FINGERPRINT",
    verdictResidual[fingerprintMutant["Local"]["Verdict"], FAILFINGERPRINT]
  ];
  Print[
    "  3e_break_gate4: incoming ell=2 gives v5_z=",
    fmt[broken["Gate4"]["Derived"]["v5_z"]],
    " -> ", broken["Local"]["Verdict"]
  ];
  expectZero[
    "3e incoming ell=2 reaches FAIL_QUAD_REGRESSION",
    verdictResidual[broken["Local"]["Verdict"], FAILQUADREGRESSION]
  ];
  Print["  3e self_ablation = ", ablation];
  expectBool["3e self-ablation dynamically reruns 022-local gate logic", ablation["RerunGateLogic"]];
  expectZero[
    "3e dynamic rerun with mutation is FAIL_QUAD_REGRESSION",
    verdictResidual[ablation["WithMutation"], FAILQUADREGRESSION]
  ];
  expectZero[
    "3e dynamic rerun without mutation is CROSS_L_FINGERPRINT_OK",
    verdictResidual[ablation["WithoutMutation"], CROSSLFINGERPRINTOK]
  ];
  expectBool["3e dynamic self-ablation changes the local verdict", ablation["FailSuppressed"]];
  expectBool[
    "neutering 3e is detected as not able to fail",
    ! neutered["FailSuppressed"] && neutered["WithMutation"] === neutered["WithoutMutation"]
  ]
];

globalSymbolNames[expr_] := Sort[DeleteDuplicates[
  SymbolName /@ Cases[
    Unevaluated[expr],
    symbol_Symbol /; Context[symbol] === "Global`",
    Infinity
  ]
]];

runAritySelfCheck[] := Module[
  {arityFingerprint, arityGate4, arityLocal, arityAblation, localDefinitions, localLHS, localArity, coreResults},
  subheading["Wolfram arity self-check and unevaluated-leakage scan"];
  arityFingerprint = fingerprintAudit[None];
  arityGate4 = gate4Audit[False, "none"];
  arityLocal = localAuditVerdict[arityFingerprint["Ok"], arityGate4["Ok"]];
  arityAblation = dynamicGate4Ablation[True];
  expectZero[
    "arity responseCoefficient[ell,sense,mutation,power] returns v5 numerator I/27",
    responseCoefficient[2, "outgoing", "none", 5] - I/27
  ];
  expectZero[
    "arity lambdaCoefficient[ell,sense,mutation,power] returns Lambda_static -3",
    lambdaCoefficient[2, "outgoing", "none", 0] + 3
  ];
  expectBool[
    "arity fingerprintAudit[mutationEll] returns Matches/Ok",
    And @@ (KeyExistsQ[arityFingerprint, #] & /@ {"Matches", "Ok"})
  ];
  expectBool[
    "arity gate4Audit[break,mutation] returns Derived/Matches/Ok",
    And @@ (KeyExistsQ[arityGate4, #] & /@ {"Derived", "Matches", "Ok"})
  ];
  expectZero[
    "arity localAuditVerdict[fingerprint,nonRegression] returns local OK",
    verdictResidual[arityLocal["Verdict"], CROSSLFINGERPRINTOK]
  ];
  expectBool[
    "arity dynamicGate4Ablation[break] changes verdicts",
    arityAblation["WithMutation"] =!= arityAblation["WithoutMutation"]
  ];
  localDefinitions = DownValues[localAuditVerdict];
  localLHS = If[
    Length[localDefinitions] === 1,
    Extract[localDefinitions, {1, 1}, HoldComplete],
    HoldComplete[]
  ];
  localArity = If[
    Length[localDefinitions] === 1,
    Extract[
      localLHS,
      {1, 1},
      Function[call, Length[Unevaluated[call]], HoldAllComplete]
    ],
    -1
  ];
  expectBool[
    "definition/call arity scan finds localAuditVerdict has exactly two arguments",
    Length[localDefinitions] === 1 && localArity === 2
  ];
  coreResults = {
    baselineContext, arityFingerprint, arityGate4, arityLocal, arityAblation,
    poleLambdaStatic /@ ells,
    firstImaginaryOrder[2, "outgoing", "lower_imaginary_copy"]
  };
  expectBool[
    "no unevaluated authored helper or SeriesCoefficient call remains in computed outputs",
    FreeQ[
      coreResults,
      hankelBuiltin | copyFactor | waveAt | lambdaAt | responseAt |
      lambdaCoefficient | responseCoefficient | imaginaryPart | realPart |
      firstImaginaryOrder | lowerImaginaryCoefficientsZero |
      lowerRealCoefficientsAgree | poleLambdaStatic | fingerprintAudit |
      quadrupoleTuple | gate4Audit | readLocalGate | localAuditVerdict |
      localContext | dynamicGate4Ablation | SeriesCoefficient | Derivative |
      D | Together | Cancel | FullSimplify | Simplify
    ]
  ]
];

runScopeAndProvenance[] := Module[{expressions, liveNames, authoredNames, forbiddenHelpers},
  subheading["022 scope and PROVENANCE consumption"];
  expressions = Flatten[Table[
    {
      waveAt[ell, "outgoing", "none"],
      lambdaAt[ell, "outgoing", "none"],
      responseAt[ell, "outgoing", "none"]
    },
    {ell, ells}
  ]];
  liveNames = globalSymbolNames[expressions];
  Print["  live symbolic names in computed fingerprint expressions = ", liveNames];
  expectBool["computed fingerprint algebra is dimensionless z-space only", liveNames === {"z"}];
  authoredNames = Names["Global`*"];
  forbiddenHelpers = {
    "Global`baseVerdict", "Global`buildTransfers", "Global`buildResiduals",
    "Global`buildRankAudit", "Global`buildPortKernel"
  };
  expectBool[
    "no source branch record, 019 prefactor, or 023 transfer/residual/nullspace helper is built",
    Intersection[authoredNames, forbiddenHelpers] === {}
  ];
  expectBool[
    "independent Wolfram route owns built-in H1/H2 definitions",
    ! FreeQ[DownValues[hankelBuiltin], SphericalHankelH1] &&
    ! FreeQ[DownValues[hankelBuiltin], SphericalHankelH2]
  ];
  Print["  z=a*omega/c_s is a provenance dictionary only; no units-restored dimensional leg is built in 022."];
  Print["  CHECKABLE: stage018 typed {u2_z=1/9,u4_z=4/81,v5_z=1/27} is the non-regression reference."];
  Print["  PROVENANCE: stage019/stage020/stage021 and the completed pathA_33 QUAD_CALIBRATED joint are cited DONE, not re-derived."];
  Print["  PROVENANCE: stage008 raw ell=0/1/2 amplitudes; stage009/stage010 bulk Helmholtz mode."];
  Print["  PROVENANCE: stage005 R1 supplies the c_s units symbol (distinct from frozen-wall c_S); no c_s value is consumed."];
  Print["  EXCLUDED/DEFERRED TO 023: residual form/sign/order, return/nullspace/selector departure, and magnitude/nonzero prediction."];
  Print["  dropped-bookkeeping: YAML/report writers, cross-engine scratch files, and all file I/O are absent."]
];

printVerdictLabels[] := (
  subheading["Verdict labels:"];
  Print["  ledger local exit-0 verdict: CROSS_L_FINGERPRINT_OK"];
  Print["  joint landing label (EARNED-first PARTIAL; FAIL not evaluated by 022): FAIL_UNDERDETERMINED_NOT_PREDICTIVE (1/2)"];
  Print["  earned here: computed cross-ell radiative coefficients/orders/Lambda_static plus the stage018 Gate-4 tuple non-regression."];
  Print["  earned able-to-fail: per-ell coefficient, pole-order, scanned-order, sign-flip, per-slot Gate-4, and dynamic 3e teeth."];
  Print["  de-counted diagnostics: static=1 and chi_Q=1; neither participates in the 022-local verdict."];
  Print["  remaining earned + FAIL-delivering/deferred content: residual form/sign/order and native nullspace/magnitude work belong to 023."];
  Print[""];
  Print["LOCAL_AUDIT_VERDICT: ", baselineLocal["Verdict"]];
  Print["JOINT_LANDING_LABEL (PARTIAL; FAIL not evaluated by 022):"];
  Print["  ", JOINTPARTIAL]
);

runAll[] := (
  heading["ledger_stage022_cross_l_fingerprints_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage022_cross_l_fingerprints"];
  Print["Engine: genuinely re-authored Wolfram built-in SphericalHankelH1/H2 + SeriesCoefficient route; no floats/tolerances; zero file I/O."];
  runCrossEllDerivation[];
  runGate4NonRegression[];
  runPerToothAblations[];
  runLocalVerdictAnd3e[];
  runAritySelfCheck[];
  runScopeAndProvenance[];
  printVerdictLabels[];
  0
);

result = Catch[
  runAll[],
  "ledgerStage022Failure",
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
