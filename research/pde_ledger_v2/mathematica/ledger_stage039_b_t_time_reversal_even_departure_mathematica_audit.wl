(* Ledger stage039 Mathematica audit: the b_T time-reversal-even departure.

   Standalone, print-only, assert-zero, exact-integer, and cross-engine-file-
   I/O free.

   Genuinely independent fixed-parity route:

   - Census records are Associations keyed by attribute names, never tuples
     read by positional indices.
   - Curl inheritance is a named-key ReplaceAll operator.  Its rotation type
     is classified independently by transforming Cross[k,u] through an
     improper rotation and testing the Det[Q] Q law.
   - The active-drain chain is a FoldList composition of parity multipliers.

   Tooth-local runtime ablation uses LEDGER_STAGE039_MUTATION.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE039_MUTATION";
activeMutation = Quiet@Check[Environment[mutationEnvironment], ""];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

verdictToken = "B_TIME_REVERSAL_EVEN";
counterfactualToken =
  "COUNTERFACTUAL_B_TIME_REVERSAL_ODD_MAXWELL_CONSISTENT";
unclassifiedToken = "UNCLASSIFIED_PARITY_STATE";

toothOrder = {
  "B_T_AXIAL_T_EVEN",
  "MAXWELL_B_T_ODD_DEPARTURE",
  "DEPARTURE_LOCALIZED_TO_T",
  "ACTIVE_DRAIN_TIME_ARROW_REQUIRED",
  "DEPARTURE_SELF_CONSISTENT",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS",
  "UNITS_RESTORED",
  "VERDICT_REDERIVATION",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "B_T_AXIAL_T_EVEN" ->
    "corrupt the named-key curl rules so derived b_T is T-odd and polar",
  "MAXWELL_B_T_ODD_DEPARTURE" ->
    "separately flip the derived-b_T side and Maxwell-benchmark side",
  "DEPARTURE_LOCALIZED_TO_T" ->
    "make the locally compared derived b_T polar",
  "ACTIVE_DRAIN_TIME_ARROW_REQUIRED" ->
    "replace active T-odd tau_d by a passive T-even FoldList root",
  "DEPARTURE_SELF_CONSISTENT" ->
    "separately break path agreement and inject time_reversal into R72",
  "TARGET_BLINDNESS" ->
    "inject a barred pathA_39 symbol into the held parity graph",
  "DUAL_ENGINE_TERMS" ->
    "drop the computed active-drain chain from the local inventory",
  "UNITS_RESTORED" ->
    "replace the cited displacement dimension L by dimensionless",
  "VERDICT_REDERIVATION" ->
    "flip computed u_T T-parity and re-run the curl/Maxwell verdict table",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "remove one scoped-out V-1 row from the live 35-tooth partition"
|>;

raise[message_] := Throw[message, "ledgerStage039Failure"];

expectBool[name_, condition_, evidence_: None] := If[
  TrueQ[condition],
  passCount += 1;
  Print["PASS  ", name],
  failCount += 1;
  Print["FIRST_FAILURE=", name];
  If[activeMutation === name, Print["FIRED_AT_OWN_ASSERT=", name]];
  Print["FAIL  ", name, ": residual = 1"];
  If[evidence =!= None, Print["      evidence = ", evidence]];
  raise[name]
];

section[text_] := (
  Print[""];
  Print[text];
  Print[StringRepeat["-", StringLength[text]]]
);


(* ---------------------------------------------------------------------- *)
(* Named-key census and native curl-parity construction.                   *)
(* ---------------------------------------------------------------------- *)

census = <|
  "s" -> <|
    "symbol" -> HoldComplete[s],
    "R_w" -> -1,
    "P_w" -> -1,
    "time_reversal" -> 1,
    "rotation" -> "scalar"
  |>,
  "V" -> <|
    "symbol" -> HoldComplete[V],
    "R_w" -> 1,
    "P_w" -> 1,
    "time_reversal" -> -1,
    "rotation" -> "polar_vector"
  |>,
  "tau_d" -> <|
    "symbol" -> HoldComplete[tauD],
    "R_w" -> 1,
    "P_w" -> 1,
    "time_reversal" -> -1,
    "rotation" -> "scalar"
  |>,
  "q_T" -> <|
    "symbol" -> HoldComplete[qT],
    "R_w" -> 1,
    "P_w" -> 1,
    "time_reversal" -> -1,
    "rotation" -> "scalar"
  |>,
  "J_T" -> <|
    "symbol" -> HoldComplete[JT],
    "R_w" -> -1,
    "P_w" -> -1,
    "time_reversal" -> 1,
    "rotation" -> "polar_vector"
  |>,
  "u_T" -> <|
    "symbol" -> HoldComplete[uT],
    "R_w" -> -1,
    "P_w" -> -1,
    "time_reversal" -> 1,
    "rotation" -> "polar_vector"
  |>,
  "b_T" -> <|
    "symbol" -> HoldComplete[bT],
    "R_w" -> -1,
    "P_w" -> -1,
    "time_reversal" -> 1,
    "rotation" -> "axial_vector"
  |>
|>;

maxwellB = <|
  "time_reversal" -> -1,
  "rotation" -> "axial_vector"
|>;

comparableKeys = {"time_reversal", "rotation"};

zeroVectorQ[vector_List] :=
  And @@ (TrueQ[PossibleZeroQ[#]] & /@ vector);

curlRotationClass[] := Module[
  {reflection, wave, displacement, reflectedCurl, axialImage, polarImage},
  reflection = DiagonalMatrix[{-1, 1, 1}];
  wave = {kx, ky, kz};
  displacement = {ux, uy, uz};
  reflectedCurl = Expand[
    Cross[reflection.wave, reflection.displacement]
  ];
  axialImage = Expand[
    Det[reflection] reflection.Cross[wave, displacement]
  ];
  polarImage = Expand[reflection.Cross[wave, displacement]];
  Which[
    zeroVectorQ[reflectedCurl - axialImage], "axial_vector",
    zeroVectorQ[reflectedCurl - polarImage], "polar_vector",
    True, "unclassified_rotation"
  ]
];

deriveCurlParity[input_Association, corrupt_: False] := Module[
  {operatorRules, numeric, rotation},
  operatorRules = If[
    TrueQ[corrupt],
    {
      CurlParity["R_w", value_] :> value,
      CurlParity["P_w", value_] :> value,
      CurlParity["time_reversal", value_] :> -value
    },
    {
      CurlParity["R_w", value_] :> value,
      CurlParity["P_w", value_] :> value,
      CurlParity["time_reversal", value_] :> value
    }
  ];
  numeric = AssociationMap[
    CurlParity[#, Lookup[input, #]] /. operatorRules &,
    {"R_w", "P_w", "time_reversal"}
  ];
  rotation = If[
    TrueQ[corrupt],
    "polar_vector",
    curlRotationClass[]
  ];
  Join[
    <|"symbol" -> HoldComplete[bT]|>,
    numeric,
    <|"rotation" -> rotation|>
  ]
];

flipTimeRecord[record_Association] := Join[
  record,
  <|"time_reversal" -> -Lookup[record, "time_reversal"]|>
];

disagreementAxes[
  model_Association,
  benchmark_Association
] := Select[
  comparableKeys,
  UnsameQ[Lookup[model, #], Lookup[benchmark, #]] &
];

departureHolds[
  model_Association,
  benchmark_Association
] := UnsameQ[
  Lookup[model, "time_reversal"],
  Lookup[benchmark, "time_reversal"]
];

localizationState[
  model_Association,
  benchmark_Association
] := {
  If[
    SameQ[Lookup[model, "rotation"], Lookup[benchmark, "rotation"]],
    Lookup[model, "rotation"],
    NoSharedRotation
  ],
  disagreementAxes[model, benchmark]
};


(* ---------------------------------------------------------------------- *)
(* FoldList active-drain chain and parity-key verdict dispatch.            *)
(* ---------------------------------------------------------------------- *)

propagateActiveDrain[tauParity_Integer] := Module[
  {multipliers, states},
  multipliers = {
    1,
    Lookup[census["s"], "time_reversal"] *
      Lookup[census["V"], "time_reversal"],
    1,
    1
  };
  states = FoldList[Times, tauParity, multipliers];
  Join[
    AssociationThread[
      {"tau_d", "q_T", "J_T", "u_T", "b_T"},
      states
    ],
    <|"moving_row_available" -> TrueQ[tauParity === -1]|>
  ]
];

citedChainState = <|
  "tau_d" -> Lookup[census["tau_d"], "time_reversal"],
  "q_T" -> Lookup[census["q_T"], "time_reversal"],
  "J_T" -> Lookup[census["J_T"], "time_reversal"],
  "u_T" -> Lookup[census["u_T"], "time_reversal"],
  "b_T" -> Lookup[census["b_T"], "time_reversal"],
  "moving_row_available" -> True
|>;

verdictFromComputed[
  model_Association,
  benchmark_Association
] := Module[{key},
  key = {
    Lookup[model, "time_reversal"],
    Lookup[benchmark, "time_reversal"],
    Lookup[model, "rotation"],
    Lookup[benchmark, "rotation"],
    disagreementAxes[model, benchmark]
  };
  Replace[
    key,
    {
      {
        1,
        -1,
        "axial_vector",
        "axial_vector",
        {"time_reversal"}
      } :> verdictToken,
      {
        -1,
        -1,
        "axial_vector",
        "axial_vector",
        {}
      } :> counterfactualToken,
      _ :> unclassifiedToken
    }
  ]
];

verdictWitnesses[mutateProduction_: False] := Module[
  {
    productionU,
    productionB,
    flippedU,
    uFlipB,
    tauFlipChain,
    tauFlipB
  },
  productionU = census["u_T"];
  If[TrueQ[mutateProduction], productionU = flipTimeRecord[productionU]];
  productionB = deriveCurlParity[productionU];

  flippedU = flipTimeRecord[census["u_T"]];
  uFlipB = deriveCurlParity[flippedU];

  tauFlipChain = propagateActiveDrain[
    -Lookup[census["tau_d"], "time_reversal"]
  ];
  tauFlipB = Join[
    deriveCurlParity[census["u_T"]],
    <|"time_reversal" -> Lookup[tauFlipChain, "b_T"]|>
  ];
  {
    verdictFromComputed[productionB, maxwellB],
    verdictFromComputed[uFlipB, maxwellB],
    verdictFromComputed[tauFlipB, maxwellB]
  }
];


(* ---------------------------------------------------------------------- *)
(* Held-symbol/DownValues blindness, computed inventory, and dimensions.   *)
(* ---------------------------------------------------------------------- *)

barredSourceMarkers = {"Nu", "aT", "aTp", "aL", "q_A_T", "q_L"};

structuralSymbolInventory[objects_List] := Module[
  {heldDefinitions, expressionGraph, names},
  heldDefinitions = DownValues /@ {
    deriveCurlParity,
    propagateActiveDrain,
    disagreementAxes,
    verdictFromComputed
  };
  expressionGraph = {objects, heldDefinitions};
  names = Cases[
    expressionGraph,
    symbol_Symbol /; Context[Unevaluated[symbol]] === "Global`" :>
      SymbolName[Unevaluated[symbol]],
    Infinity,
    Heads -> True
  ];
  DeleteDuplicates[names]
];

targetBlindnessState[
  derivedB_Association,
  chain_Association,
  mutate_: False
] := Module[{objects},
  objects = {
    census,
    derivedB,
    maxwellB,
    chain,
    disagreementAxes[derivedB, maxwellB]
  };
  If[TrueQ[mutate], AppendTo[objects, HoldComplete[Nu]]];
  Intersection[
    structuralSymbolInventory[objects],
    barredSourceMarkers
  ]
];

computedInventory[
  derivedB_Association,
  chain_Association,
  mutate_: False
] := Module[{inventory},
  inventory = <|
    "derived_b_T" -> derivedB,
    "maxwell_B" -> maxwellB,
    "departure_holds" -> departureHolds[derivedB, maxwellB],
    "disagreement_axes" -> disagreementAxes[derivedB, maxwellB],
    "active_drain_chain" -> chain,
    "verdict" -> verdictFromComputed[derivedB, maxwellB]
  |>;
  If[TrueQ[mutate], inventory = KeyDrop[inventory, "active_drain_chain"]];
  inventory
];

expectedInventory = <|
  "derived_b_T" -> census["b_T"],
  "maxwell_B" -> <|
    "time_reversal" -> -1,
    "rotation" -> "axial_vector"
  |>,
  "departure_holds" -> True,
  "disagreement_axes" -> {"time_reversal"},
  "active_drain_chain" -> citedChainState,
  "verdict" -> verdictToken
|>;

dimensionless = {0, 0, 0};
lengthDimension = {0, 1, 0};
inverseLengthDimension = {0, -1, 0};

multiplyDimensions[left_List, right_List] :=
  MapThread[Plus, {left, right}];

unitState[corrupt_: False] := Module[
  {uDimension, curlUDimension, bDimension},
  uDimension = If[TrueQ[corrupt], dimensionless, lengthDimension];
  curlUDimension = multiplyDimensions[
    inverseLengthDimension,
    uDimension
  ];
  bDimension = dimensionless;
  {curlUDimension, bDimension}
];


(* ---------------------------------------------------------------------- *)
(* Complete 35-source-tooth manifest and five stage-native authored teeth. *)
(* ---------------------------------------------------------------------- *)

sourceToothIDs = {
  "SOURCE_TRANSLATION_CONTINUITY",
  "SOURCE_NOT_IMPORTED",
  "SOURCE_BASIS",
  "PARITY_RW",
  "PARITY_PW",
  "PARITY_ROTATION",
  "PARITY_TIME_REVERSAL",
  "FIELD_IDENTITY_UNITS",
  "ACTION_KINETIC",
  "ACTION_COUPLING",
  "ACTION_STABILITY",
  "G0_DAMAGE",
  "ROUTE_INDEPENDENCE",
  "BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER",
  "BOOST_COMMON_VELOCITY",
  "DIRECT_SOURCE",
  "DIRECT_PROJECTOR",
  "DIRECT_EXCHANGE_SIGN",
  "DIRECT_FALLOFF",
  "DIRECT_VELOCITY_ORDER",
  "COMPARE_COMPUTED",
  "DELTA_RATIO",
  "CONE_RATIO",
  "QMAG_R1",
  "UNITS_RESTORED",
  "ACTIVE_FLUX_CAVEAT",
  "HOOK_LORENTZ",
  "LEDGER_READY_ROW",
  "TRUTH_TOTALITY",
  "TRUTH_PRECEDENCE",
  "LANDING_OWNERSHIP",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS"
};

citedSource = {
  "PARITY_ROTATION",
  "PARITY_TIME_REVERSAL"
};

scopedOutSource = {
  "FIELD_IDENTITY_UNITS",
  "ACTION_KINETIC",
  "ACTION_COUPLING",
  "ACTION_STABILITY",
  "G0_DAMAGE",
  "LEDGER_READY_ROW",
  "SOURCE_TRANSLATION_CONTINUITY",
  "SOURCE_NOT_IMPORTED",
  "SOURCE_BASIS",
  "PARITY_RW",
  "PARITY_PW",
  "BOOST_PROJECTOR",
  "BOOST_GENERAL_VELOCITIES",
  "BOOST_NEXT_ORDER",
  "DIRECT_SOURCE",
  "DIRECT_PROJECTOR",
  "DIRECT_EXCHANGE_SIGN",
  "DIRECT_FALLOFF",
  "DIRECT_VELOCITY_ORDER",
  "ROUTE_INDEPENDENCE",
  "BOOST_COMMON_VELOCITY",
  "COMPARE_COMPUTED",
  "DELTA_RATIO",
  "CONE_RATIO",
  "QMAG_R1",
  "TRUTH_TOTALITY",
  "TRUTH_PRECEDENCE",
  "LANDING_OWNERSHIP",
  "ACTIVE_FLUX_CAVEAT",
  "HOOK_LORENTZ"
};

buildGlobalSource = {
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS",
  "UNITS_RESTORED"
};

stageNativeAuthored = Take[toothOrder, 5];

sourceDisposition[identifier_String] := Which[
  MemberQ[citedSource, identifier],
    {"CITED", "STAGE035_V2_DONE"},
  MemberQ[buildGlobalSource, identifier],
    {"REPLACED_BY_STRONGER", "STAGE039_BUILD_GLOBAL"},
  MemberQ[Take[scopedOutSource, 6], identifier],
    {"SCOPED_OUT", "STAGE034_V1_DONE"},
  MemberQ[Take[scopedOutSource, {7, 11}], identifier],
    {"SCOPED_OUT", "STAGE035_V2_DONE"},
  MemberQ[Take[scopedOutSource, {12, 14}], identifier],
    {"SCOPED_OUT", "STAGE036_V3_DONE"},
  MemberQ[Take[scopedOutSource, {15, 25}], identifier],
    {"SCOPED_OUT", "STAGE037_V4_DONE"},
  MemberQ[Take[scopedOutSource, {26, 30}], identifier],
    {"SCOPED_OUT", "STAGE038_V5_DONE"},
  True,
    raise["UNPARTITIONED_SOURCE_TOOTH"]
];

sourceManifest = Map[
  Function[
    identifier,
    Join[{identifier}, sourceDisposition[identifier]]
  ],
  sourceToothIDs
];

manifestState[manifest_List] := Module[
  {dispositions},
  dispositions = Sort[Normal@Counts[#[[2]] & /@ manifest]];
  {
    First /@ manifest,
    Sort@Cases[manifest, {identifier_, "CITED", _} :> identifier],
    Sort@Cases[manifest, {identifier_, "SCOPED_OUT", _} :> identifier],
    Sort@Cases[
      manifest,
      {identifier_, "REPLACED_BY_STRONGER", _} :> identifier
    ],
    dispositions,
    manifest,
    Length[manifest]
  }
];

expectedManifestState = {
  sourceToothIDs,
  Sort[citedSource],
  Sort[scopedOutSource],
  Sort[buildGlobalSource],
  Sort[{
    "CITED" -> 2,
    "REPLACED_BY_STRONGER" -> 3,
    "SCOPED_OUT" -> 30
  }],
  sourceManifest,
  35
};


(* ---------------------------------------------------------------------- *)
(* Ten executable teeth.                                                   *)
(* ---------------------------------------------------------------------- *)

ok = Catch[
  If[
    activeMutation =!= "" &&
      !MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print[
    "ledger_stage039_b_t_time_reversal_even_departure Mathematica audit"
  ];
  Print[
    "ROUTE=named-key ReplaceAll curl + Cross/Det classifier + FoldList parity chain"
  ];
  Print[
    "INDEPENDENCE=no parity tuple indices; no source-.wl tuple checks; no source text scan"
  ];
  Print["FILE_IO=none; CROSS_ENGINE_COMPARE=none"];
  If[
    activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print[
      "MUTATED_PRIMITIVE=",
      Lookup[ablationDescriptions, activeMutation]
    ]
  ];

  section["Cited census -> rule-based curl inheritance"];
  derivedB = deriveCurlParity[
    census["u_T"],
    activeMutation === "B_T_AXIAL_T_EVEN"
  ];
  expectBool[
    "B_T_AXIAL_T_EVEN",
    SameQ[derivedB, census["b_T"]],
    <|
      "Derived" -> derivedB,
      "Cited035" -> census["b_T"],
      "RotationClassifier" -> "Cross/Det improper-rotation law"
    |>
  ];
  Print[
    "      b_T=(-1,-1,+1,axial_vector), derived from u_T=(-1,-1,+1,polar_vector)"
  ];

  section["Typed Maxwell-B comparison over exactly {T, rotation}"];
  productionDeparture = departureHolds[derivedB, maxwellB];
  If[
    activeMutation === "MAXWELL_B_T_ODD_DEPARTURE",
    bSideDeparture = departureHolds[
      deriveCurlParity[flipTimeRecord[census["u_T"]]],
      maxwellB
    ];
    maxwellSideDeparture = departureHolds[
      derivedB,
      Join[maxwellB, <|"time_reversal" -> 1|>]
    ];
    departureCheck = And[bSideDeparture, maxwellSideDeparture],
    bSideDeparture = None;
    maxwellSideDeparture = None;
    departureCheck = productionDeparture
  ];
  expectBool[
    "MAXWELL_B_T_ODD_DEPARTURE",
    departureCheck,
    <|
      "DerivedBT" -> Lookup[derivedB, "time_reversal"],
      "MaxwellBT" -> Lookup[maxwellB, "time_reversal"],
      "BSideAblation" -> bSideDeparture,
      "MaxwellSideAblation" -> maxwellSideDeparture
    |>
  ];
  Print["      b_T is T-even; Maxwell B is T-odd; departure=True"];

  localizedB = If[
    activeMutation === "DEPARTURE_LOCALIZED_TO_T",
    Join[derivedB, <|"rotation" -> "polar_vector"|>],
    derivedB
  ];
  liveLocalization = localizationState[localizedB, maxwellB];
  expectedLocalization = {
    "axial_vector",
    {"time_reversal"}
  };
  expectBool[
    "DEPARTURE_LOCALIZED_TO_T",
    SameQ[liveLocalization, expectedLocalization],
    <|
      "SharedRotation" -> First[liveLocalization],
      "Disagreements" -> Last[liveLocalization],
      "ComparableDomain" -> comparableKeys
    |>
  ];
  decidedAxes = disagreementAxes[derivedB, maxwellB];
  Print["      rotation=axial agrees; disagreement={time_reversal}"];

  section["Active-drain tau_d time-arrow propagation"];
  tauRoot = Lookup[census["tau_d"], "time_reversal"];
  If[
    activeMutation === "ACTIVE_DRAIN_TIME_ARROW_REQUIRED",
    tauRoot = -tauRoot
  ];
  activeChain = propagateActiveDrain[tauRoot];
  expectBool[
    "ACTIVE_DRAIN_TIME_ARROW_REQUIRED",
    SameQ[activeChain, citedChainState],
    <|
      "ComputedChain" -> activeChain,
      "CitedChain" -> citedChainState
    |>
  ];
  Print["      tau_d(-1)->q_T(-1)->J_T(+1)->u_T(+1)->b_T(+1)"];
  Print["      passive T-even throat: no O(V) moving row"];

  productionSelfState = {
    SameQ[
      Lookup[derivedB, "time_reversal"],
      Lookup[activeChain, "b_T"]
    ],
    Intersection[decidedAxes, {"sign", "magnitude"}] === {}
  };
  pathAblationState = {
    SameQ[
      -Lookup[derivedB, "time_reversal"],
      Lookup[activeChain, "b_T"]
    ],
    Last[productionSelfState]
  };
  axisAblationState = {
    First[productionSelfState],
    Intersection[
      decidedAxes,
      {"sign", "magnitude", "time_reversal"}
    ] === {}
  };
  liveSelfState = If[
    activeMutation === "DEPARTURE_SELF_CONSISTENT",
    {First[pathAblationState], Last[axisAblationState]},
    productionSelfState
  ];
  expectBool[
    "DEPARTURE_SELF_CONSISTENT",
    SameQ[liveSelfState, {True, True}],
    <|
      "Production" -> productionSelfState,
      "PathOnlyAblation" -> pathAblationState,
      "AxisOnlyAblation" -> axisAblationState,
      "R72" -> "R1_REQUIRED(electric_bc_selection)",
      "R72Axes" -> {"sign", "magnitude"}
    |>
  ];
  Print[
    "      curl-path T == tau-chain T; {time_reversal} disjoint {sign,magnitude}"
  ];
  Print[
    "      separate ablations: path={False,True}; R72-axis injection={True,False}"
  ];

  section["Build-global structural, local-inventory, and unit firewalls"];
  blindness = targetBlindnessState[
    derivedB,
    activeChain,
    activeMutation === "TARGET_BLINDNESS"
  ];
  expectBool[
    "TARGET_BLINDNESS",
    SameQ[blindness, {}],
    <|
      "BarredIntersection" -> blindness,
      "Barred" -> barredSourceMarkers
    |>
  ];
  Print[
    "      held-symbol plus core-DownValues graph; no source text scan"
  ];

  inventory = computedInventory[
    derivedB,
    activeChain,
    activeMutation === "DUAL_ENGINE_TERMS"
  ];
  expectBool[
    "DUAL_ENGINE_TERMS",
    SameQ[inventory, expectedInventory],
    <|
      "Terms" -> Keys[inventory],
      "DerivedBT" -> Lookup[inventory, "derived_b_T"],
      "MaxwellB" -> Lookup[inventory, "maxwell_B"],
      "Departure" -> Lookup[inventory, "departure_holds"],
      "Disagreements" -> Lookup[inventory, "disagreement_axes"],
      "Verdict" -> Lookup[inventory, "verdict"]
    |>
  ];
  Print[
    "      local terms=derived b_T,Maxwell B,departure,T-axis,active-drain chain,verdict"
  ];

  dimensions = unitState[
    activeMutation === "UNITS_RESTORED"
  ];
  expectBool[
    "UNITS_RESTORED",
    SameQ[First[dimensions], Last[dimensions]],
    <|
      "CurlUT" -> First[dimensions],
      "BT" -> Last[dimensions]
    |>
  ];
  Print["      [curl u_T]=L^-1*L=1=[b_T]"];

  section["Computed verdict re-derivation table"];
  actualVerdicts = verdictWitnesses[
    activeMutation === "VERDICT_REDERIVATION"
  ];
  expectedVerdicts = {
    verdictToken,
    counterfactualToken,
    counterfactualToken
  };
  If[
    activeMutation === "VERDICT_REDERIVATION",
    Print[
      "      MUTATION_REDERIVED=",
      counterfactualToken,
      " via u_T flip and tau_d flip"
    ]
  ];
  expectBool[
    "VERDICT_REDERIVATION",
    SameQ[actualVerdicts, expectedVerdicts],
    <|
      "Production" -> First[actualVerdicts],
      "UTFlip" -> actualVerdicts[[2]],
      "TauDFlip" -> Last[actualVerdicts]
    |>
  ];
  Print[
    "      production verdict re-derived from b_T/Maxwell parity objects"
  ];

  section["Canonical source-to-stage predicate manifest"];
  liveManifest = sourceManifest;
  If[
    activeMutation === "SOURCE_TO_STAGE_MANIFEST",
    liveManifest = Select[
      sourceManifest,
      First[#] =!= First[scopedOutSource] &
    ]
  ];
  liveManifestState = manifestState[liveManifest];
  expectBool[
    "SOURCE_TO_STAGE_MANIFEST",
    SameQ[liveManifestState, expectedManifestState],
    <|
      "SourceTotal" -> Last[liveManifestState],
      "Partition" -> liveManifestState[[5]],
      "Cited" -> liveManifestState[[2]],
      "ScopedOut" -> liveManifestState[[3]],
      "BuildGlobal" -> liveManifestState[[4]]
    |>
  ];
  Print[
    "      source_manifest=35; cited=2; scoped_out=30; build_global_replaced_by_stronger=3"
  ];
  Print["      CITED=", StringRiffle[citedSource, ","]];
  Print["      SCOPED_OUT=", StringRiffle[scopedOutSource, ","]];
  Print[
    "      STAGE_NATIVE_AUTHORED=",
    StringRiffle[stageNativeAuthored, ","]
  ];

  Print[""];
  Print["VERDICT_TOKEN=", verdictToken];
  Print["SCOPE_CLASS=CHARACTERIZED-DEPARTURE"];
  Print["FRAMING=FIRST_CLASS_LOAD_BEARING_NOT_EXACT_MAXWELL"];
  Print["NOT_A_BUG=TRUE; NOT_R1=TRUE; NEVER_SOFTENED=TRUE"];
  Print["MAXWELL_COMPARISON_DOMAIN=time_reversal,rotation"];
  Print["R72_CITED=R1_REQUIRED(electric_bc_selection)"];
  Print["R72_ORTHOGONAL_AXES=sign,magnitude"];
  Print["STAGE038_POST_RESOLUTION_TOKEN_REEMITTED=FALSE"];
  Print[
    "SIBLING_DEPARTURES=NATIVE_P_NO_EMERGENT_GAUSS,FAIL_CAUCHY_STRAY_LONGITUDINAL"
  ];

  If[
    activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage039Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[
  TrueQ[ok],
  Print[
    "OVERALL PASS: Mathematica independently reached ",
    verdictToken
  ];
  Exit[0],
  Print[
    "OVERALL FAIL: Mathematica stage039 audit did not close"
  ];
  Exit[1]
]
