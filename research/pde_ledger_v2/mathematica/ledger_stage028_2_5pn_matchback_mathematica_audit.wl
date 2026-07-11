(* Ledger stage028 Wolfram audit: calibrated 2.5PN match-back consistency.

   Standalone, print-only, no arguments, and zero file I/O.  This is a native
   re-author: mutations replace complete moment expressions, INV1 is decided
   by a Groebner polynomial identity, and the remaining zero tests use
   PossibleZeroQ over positive-domain refined expressions.  It does not use
   the source's coefficient-config/FullSimplify reduction pipeline.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
diagnosticPassCount = 0;
diagnosticFailCount = 0;
failureMessage = "";

localVerdict = "MATCHBACK_CONSISTENT";

$Assumptions = Element[{G, cs, a, c}, Reals] && G > 0 && cs > 0 && a > 0 && c > 0;

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

fail[] := Throw[failureMessage, "ledgerStage028Failure"];

recordPass[name_String] := (
  passCount++;
  Print["PASS  ", name]
);

recordFail[name_String, detail_String] := (
  failCount++;
  failureMessage = name <> ": " <> detail;
  Print["FAIL  ", failureMessage]
);

diagnosticBool[name_String, condition_] := If[
  TrueQ[condition],
  diagnosticPassCount++;
  Print["PASS  ", name],
  diagnosticFailCount++;
  failureMessage = name;
  Print["FAIL  ", name];
  fail[]
];

rigAssert[name_String, condition_] := If[
  ! TrueQ[condition],
  Throw[name, "stage028RigAssertion"]
];

exerciseRig[mutationName_String, actual_List, expected_List] := Module[
  {assertionName, fired},
  assertionName = "ROW " <> mutationName <> " actual == EXPECTED_CAUGHT_BY";
  fired = Catch[
    rigAssert[assertionName, actual === expected];
    "NO_ASSERT_FIRED",
    "stage028RigAssertion",
    Function[{message, tag}, message]
  ];
  If[fired === "NO_ASSERT_FIRED",
    recordPass[assertionName],
    recordFail[
      assertionName,
      "actual=" <> ToString[actual, InputForm] <>
      "; expected=" <> ToString[expected, InputForm]
    ];
    fail[]
  ]
];

residualOrder = {
  "INV1_moment_invariant",
  "INV2_pathA43_form_to_BT",
  "INV3_corpus_form_to_BT",
  "INV4_redundant_cross_form",
  "INV5_Kbar0_coeff_54_5",
  "INV5_Kbar2_coeff_6_5",
  "INV5_Kbar4_coeff_8_15",
  "INV5_BT_coeff_2_5",
  "INV5_pathA43_denominator_27",
  "INV5_corpus_factor_9",
  "INV5_exp_Kbar2_5_2",
  "INV5_exp_Kbar0_3_2"
};

(* Stage027 closure moments, locally restated as PROVENANCE. *)
kbar0Base = 54 G cs^5/(5 a^5 c^5);
kbar2Base = 6 G cs^3/(5 a^3 c^5);
kbar4Base = 8 G cs/(15 a c^5);
btBase = 2 G/(5 c^5);

(* Stage028-owned anti-rig literals: never computed from the moments. *)
anchorK0 = 54/5;
anchorK2 = 6/5;
anchorK4 = 8/15;
anchorBT = 2/5;
anchorDenominator = 27;
anchorCorpusFactor = 9;
anchorExpK2 = 5/2;
anchorExpK0 = 3/2;

pointSymbols = {moment0, moment2, moment4, btTarget, pathDen, corpusFactor, exponent2, exponent0};
basePointRules = {
  moment0 -> kbar0Base,
  moment2 -> kbar2Base,
  moment4 -> kbar4Base,
  btTarget -> btBase,
  pathDen -> 27,
  corpusFactor -> 9,
  exponent2 -> 5/2,
  exponent0 -> 3/2
};

(* Direct complete-expression perturbations; no coefficient config exists. *)
mutationSpecs = {
  <|
    "Name" -> "coherent_scale_all_moments_and_BT_x2",
    "Updates" -> {
      moment0 -> 2 kbar0Base, moment2 -> 2 kbar2Base,
      moment4 -> 2 kbar4Base, btTarget -> 2 btBase
    }
  |>,
  <|
    "Name" -> "coupled_moments_x2_BT_fixed",
    "Updates" -> {
      moment0 -> 2 kbar0Base, moment2 -> 2 kbar2Base,
      moment4 -> 2 kbar4Base
    }
  |>,
  <|
    "Name" -> "Kbar4_coeff_8_15_to_8_14",
    "Updates" -> {moment4 -> 8 G cs/(14 a c^5)}
  |>,
  <|
    "Name" -> "Kbar4_sign_flip",
    "Updates" -> {moment4 -> -kbar4Base}
  |>,
  <|
    "Name" -> "Kbar2_coeff_6_5_to_7_5",
    "Updates" -> {moment2 -> 7 G cs^3/(5 a^3 c^5)}
  |>,
  <|
    "Name" -> "Kbar0_coeff_54_5_to_55_5",
    "Updates" -> {moment0 -> 55 G cs^5/(5 a^5 c^5)}
  |>,
  <|
    "Name" -> "pathA43_denominator_27_to_26",
    "Updates" -> {pathDen -> 26}
  |>,
  <|
    "Name" -> "corpus_factor_9_to_8",
    "Updates" -> {corpusFactor -> 8}
  |>,
  <|
    "Name" -> "exp_Kbar2_5_2_to_3_2",
    "Updates" -> {exponent2 -> 3/2}
  |>,
  <|
    "Name" -> "exp_Kbar0_3_2_to_1",
    "Updates" -> {exponent0 -> 1}
  |>,
  <|
    "Name" -> "BT_coeff_2_5_to_3_5",
    "Updates" -> {btTarget -> 3 G/(5 c^5)}
  |>
};

expectedCaughtBy = <|
  "coherent_scale_all_moments_and_BT_x2" -> {
    "INV5_Kbar0_coeff_54_5", "INV5_Kbar2_coeff_6_5",
    "INV5_Kbar4_coeff_8_15", "INV5_BT_coeff_2_5"
  },
  "coupled_moments_x2_BT_fixed" -> {
    "INV2_pathA43_form_to_BT", "INV3_corpus_form_to_BT",
    "INV5_Kbar0_coeff_54_5", "INV5_Kbar2_coeff_6_5",
    "INV5_Kbar4_coeff_8_15"
  },
  "Kbar4_coeff_8_15_to_8_14" -> {
    "INV1_moment_invariant", "INV5_Kbar4_coeff_8_15"
  },
  "Kbar4_sign_flip" -> {
    "INV1_moment_invariant", "INV5_Kbar4_coeff_8_15"
  },
  "Kbar2_coeff_6_5_to_7_5" -> {
    "INV1_moment_invariant", "INV3_corpus_form_to_BT",
    "INV4_redundant_cross_form", "INV5_Kbar2_coeff_6_5"
  },
  "Kbar0_coeff_54_5_to_55_5" -> {
    "INV1_moment_invariant", "INV2_pathA43_form_to_BT",
    "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
    "INV5_Kbar0_coeff_54_5"
  },
  "pathA43_denominator_27_to_26" -> {
    "INV2_pathA43_form_to_BT", "INV4_redundant_cross_form",
    "INV5_pathA43_denominator_27"
  },
  "corpus_factor_9_to_8" -> {
    "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
    "INV5_corpus_factor_9"
  },
  "exp_Kbar2_5_2_to_3_2" -> {
    "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
    "INV5_exp_Kbar2_5_2"
  },
  "exp_Kbar0_3_2_to_1" -> {
    "INV3_corpus_form_to_BT", "INV4_redundant_cross_form",
    "INV5_exp_Kbar0_3_2"
  },
  "BT_coeff_2_5_to_3_5" -> {
    "INV2_pathA43_form_to_BT", "INV3_corpus_form_to_BT",
    "INV5_BT_coeff_2_5"
  }
|>;

pointFromUpdates[updates_List] := pointSymbols /. Dispatch[Join[updates, basePointRules]];

rawResiduals[values_List] := Module[
  {k0, k2, k4, bt, denominator, factor9, exp2, exp0, pathForm, corpusForm},
  {k0, k2, k4, bt, denominator, factor9, exp2, exp0} = values;
  pathForm = k0 a^5/(denominator cs^5);
  corpusForm = factor9 k2^exp2/k0^exp0;
  {
    k4 k0 - 4 k2^2,
    pathForm - bt,
    corpusForm - bt,
    pathForm - corpusForm,
    k0 a^5 c^5/(G cs^5) - anchorK0,
    k2 a^3 c^5/(G cs^3) - anchorK2,
    k4 a c^5/(G cs) - anchorK4,
    bt c^5/G - anchorBT,
    denominator - anchorDenominator,
    factor9 - anchorCorpusFactor,
    exp2 - anchorExpK2,
    exp0 - anchorExpK0
  }
];

(* INV1's verdict route is a polynomial-ideal identity. *)
groebnerZeroQ[expr_] := Module[{numerator, basis},
  numerator = Numerator[Together[PowerExpand[Refine[expr, $Assumptions]]]];
  basis = GroebnerBasis[{numerator}, {G, cs, a, c}];
  TrueQ[basis === {}]
];

(* INV2-INV5 use a native zero predicate, not reduced-expression equality. *)
nativeZeroQ[expr_] := TrueQ[
  PossibleZeroQ[PowerExpand[Refine[expr, $Assumptions]]]
];

zeroFlags[values_List] := Module[{raw = rawResiduals[values]},
  Join[{groebnerZeroQ[First[raw]]}, nativeZeroQ /@ Rest[raw]]
];

displayResiduals[values_List] := (
  Factor[Cancel[Together[PowerExpand[Refine[#, $Assumptions]]]]] & /@
    rawResiduals[values]
);

externalText[expr_] := StringReplace[
  ToString[expr, InputForm],
  {"cs" -> "c_s"}
];

vectorText[residuals_List] := "[" <> StringRiffle[
  MapThread[#1 <> "=" <> externalText[#2] &, {residualOrder, residuals}],
  ", "
] <> "]";

zeroTestText[flags_List] := "[" <> StringRiffle[
  MapThread[#1 <> "=" <> ToString[#2, InputForm] &, {residualOrder, flags}],
  ", "
] <> "]";

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

arityScan[function_Symbol, held_HoldComplete] :=
  definitionArity[function] === heldCallArity[held];

leakageScan[objects_] := FreeQ[
  objects,
  pointFromUpdates | rawResiduals | groebnerZeroQ | nativeZeroQ | zeroFlags |
  displayResiduals | PossibleZeroQ | GroebnerBasis | Refine | PowerExpand
];

noFloatQ[objects_] := FreeQ[HoldComplete[objects], _Real];

printProvenance[] := (
  subheading["PROVENANCE and A3 boundary"];
  Print["Kbar0,Kbar2,Kbar4 are CALIBRATED closure inputs restated from upstream stage027; moments are not derived here."];
  Print["External calibrated literals: 2/5, 54/5, and G; G=GENUINE_BLOCKED; Gamma5/G is not derived here."];
  Print["The pathA_43 denominator 27 is the upstream-EARNED density-Hankel fingerprint imported from stage018."];
  Print["Corpus form 9*Kbar2^(5/2)/Kbar0^(3/2) is imported from 4d_2_5pn.tex:469, not re-derived."];
  Print["Runtime isolation: no peer Get/Run, no report/source/note/_scratch reads, no writes; zero file I/O."];
  Print["A3 BOUNDARY: INV1/INV2 are SHARED re-expressions of stage027 closure (R46)."];
  Print["A3 BOUNDARY: INV3/INV4/INV5 and the coherent-rescale matrix are stage028 content."];
  Print["A3 BOUNDARY: INV4 == INV2-INV3 identically; retained as a redundant diagnostic, not an independent tooth."]
);

runAudit[] := Module[
  {baselinePoint, baselineRaw, baselineDisplay, baselineFlags, allTranscript,
   caughtMatrix, name, point, raw, display, flags, caught, expected,
   pointCall, rawCall, groebnerCall, nativeCall, zeroCall, displayCall,
   leakageCall, floatCall, rigCall, localOK},

  printProvenance[];
  baselinePoint = pointFromUpdates[{}];
  Print["RESTATED stage027 closure_overlay moments:"];
  Print["  Kbar0=", externalText[kbar0Base]];
  Print["  Kbar2=", externalText[kbar2Base]];
  Print["  Kbar4=", externalText[kbar4Base]];
  Print["  Gamma5=Kbar0*a^5/(27*c_s^5)=", externalText[2 G/(5 c^5)]];
  Print["  A3_FIDELITY: restated literals match stage027 closure_overlay (authoring-time comparison)."];

  subheading["Baseline exact-rational residuals and zero tests"];
  baselineRaw = rawResiduals[baselinePoint];
  baselineDisplay = displayResiduals[baselinePoint];
  baselineFlags = zeroFlags[baselinePoint];
  Print["BASELINE RESIDUAL VECTOR:"];
  Print["  ", vectorText[baselineDisplay]];
  Print["BASELINE ZERO TESTS:"];
  Print["  ", zeroTestText[baselineFlags]];
  Do[
    diagnosticBool[
      "ZERO-TEST " <> residualOrder[[index]] <> " == 0",
      baselineFlags[[index]]
    ],
    {index, Length[residualOrder]}
  ];
  Print["BASELINE ALL-ZERO: PASS (aggregate positive; DE-COUNTED)"];

  subheading["Runtime-computed mutation residuals"];
  caughtMatrix = <||>;
  allTranscript = Join[baselineRaw, baselineDisplay, baselineFlags];
  Do[
    name = spec["Name"];
    point = pointFromUpdates[spec["Updates"]];
    raw = rawResiduals[point];
    display = displayResiduals[point];
    flags = zeroFlags[point];
    caught = Pick[residualOrder, flags, False];
    AssociateTo[caughtMatrix, name -> caught];
    allTranscript = Join[allTranscript, raw, display, flags];
    Print["MUTATION ", name, ":"];
    Print["  residuals: ", vectorText[display]];
    Print["  zero_tests: ", zeroTestText[flags]];
    Print["  caught_by: [", StringRiffle[caught, ", "], "]"],
    {spec, mutationSpecs}
  ];

  subheading["Immutable EXPECTED caught-by matrix"];
  Do[
    name = spec["Name"];
    expected = expectedCaughtBy[name];
    caught = caughtMatrix[name];
    Print[
      "ROW ", name,
      ": expected=[", StringRiffle[expected, ", "],
      "] actual=[", StringRiffle[caught, ", "], "]"
    ];
    exerciseRig[name, caught, expected],
    {spec, mutationSpecs}
  ];
  Print["MUTATION PROBE: PASS (11 runtime-computed rows)"];

  pointCall = HoldComplete[pointFromUpdates[{}]];
  rawCall = HoldComplete[rawResiduals[baselinePoint]];
  groebnerCall = HoldComplete[groebnerZeroQ[First[baselineRaw]]];
  nativeCall = HoldComplete[nativeZeroQ[baselineRaw[[2]]]];
  zeroCall = HoldComplete[zeroFlags[baselinePoint]];
  displayCall = HoldComplete[displayResiduals[baselinePoint]];
  leakageCall = HoldComplete[leakageScan[allTranscript]];
  floatCall = HoldComplete[noFloatQ[allTranscript]];
  rigCall = HoldComplete[exerciseRig[name, caught, expected]];
  diagnosticBool[
    "ARITY every verdict-bearing re-authored helper matches its held call",
    arityScan[pointFromUpdates, pointCall] &&
    arityScan[rawResiduals, rawCall] &&
    arityScan[groebnerZeroQ, groebnerCall] &&
    arityScan[nativeZeroQ, nativeCall] &&
    arityScan[zeroFlags, zeroCall] &&
    arityScan[displayResiduals, displayCall] &&
    arityScan[leakageScan, leakageCall] &&
    arityScan[noFloatQ, floatCall] &&
    arityScan[exerciseRig, rigCall]
  ];
  diagnosticBool[
    "LEAKAGE transcript has no unevaluated authored helper or zero-test kernel",
    leakageScan[allTranscript]
  ];

  If[noFloatQ[allTranscript],
    recordPass["NO-FLOAT exact-rational guard over raw and reduced residuals"],
    recordFail[
      "NO-FLOAT exact-rational guard over raw and reduced residuals",
      "machine Real found"
    ];
    fail[]
  ];
  Print["NO-FLOAT GUARD: PASS"];

  localOK = TrueQ[
    And @@ baselineFlags &&
    And @@ KeyValueMap[#2 === expectedCaughtBy[#1] &, caughtMatrix] &&
    noFloatQ[allTranscript]
  ];
  If[! localOK,
    failureMessage = "LOCAL MATCHBACK_CONSISTENT gate";
    fail[]
  ];

  subheading["LOCAL calibrated-consistency verdict"];
  Print["LOCAL_AUDIT_VERDICT: ", localVerdict];
  Print["scope: CALIBRATED consistency over closure moments; G=GENUINE_BLOCKED."];
  Print["scope: NOT a first-principles Gamma5/G derivation; full 1PN->4PN from-throat is SIM-DEFERRED Gate 6."];
  0
];

runAll[] := (
  heading["ledger_stage028_2_5pn_matchback_mathematica_audit"];
  Print["Target stem confirmed: ledger_stage028_2_5pn_matchback"];
  Print["Engine: Wolfram Groebner/PossibleZeroQ direct-expression route; no floats; zero file I/O."];
  runAudit[]
);

result = Catch[
  runAll[],
  "ledgerStage028Failure",
  Function[{message, tag}, 1]
];

heading["Tallies"];
Print[
  "TALLY mathematica: ", passCount, " pass + ", failCount,
  " fail = ", passCount + failCount,
  " independent EXIT-1 teeth; diagnostics=", diagnosticPassCount,
  " pass + ", diagnosticFailCount, " fail"
];
If[result === 0 && failCount === 0 && diagnosticFailCount === 0,
  Print["OVERALL PASS"];
  Exit[0],
  If[failureMessage =!= "", Print["ABORTED: ", failureMessage]];
  Print["OVERALL FAIL"];
  Exit[1]
];
