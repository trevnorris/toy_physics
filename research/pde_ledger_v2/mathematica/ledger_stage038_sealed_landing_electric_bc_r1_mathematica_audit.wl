(* Ledger stage038 Mathematica audit: SEALED section-4 terminal landing.

   Standalone, print-only, assert-zero, exact, and cross-engine-file-I/O-free.

   Genuinely independent landing-stage route:

   - The Cartesian domain is enumerated by ordinal mixed-radix decoding, with
     the inconsistency digit fastest.  It does not use Tuples.
   - The adjudicator is a declarative ordered predicate table resolved by
     FirstCase.  It is not the source build's nested If cascade.
   - The oracle separately classifies complete cases by a composite lookup key
     and incomplete cases by the first active open-reason index.
   - Canonical rows are encoded to UTF-8 integer chunks, interleaved with byte
     10, flattened into ByteArray, and hashed at runtime.

   Tooth-local runtime ablation uses LEDGER_STAGE038_MUTATION.
*)

ClearAll["Global`*"];

passCount = 0;
failCount = 0;
mutationEnvironment = "LEDGER_STAGE038_MUTATION";
activeMutation = Quiet@Check[Environment[mutationEnvironment], ""];
If[!StringQ[activeMutation], activeMutation = ""];
activeMutation = StringTrim[activeMutation];

toothOrder = {
  "TRUTH_TOTALITY",
  "TRUTH_PRECEDENCE",
  "LANDING_OWNERSHIP",
  "ACTIVE_FLUX_CAVEAT",
  "HOOK_LORENTZ",
  "TARGET_BLINDNESS",
  "DUAL_ENGINE_TERMS",
  "UNITS_RESTORED",
  "VERDICT_REDERIVATION",
  "SOURCE_TO_STAGE_MANIFEST"
};

ablationDescriptions = <|
  "TRUTH_TOTALITY" ->
    "replace one computed cell by the defensive fall-through and drop a distinct row",
  "TRUTH_PRECEDENCE" ->
    "move the consistency arm ahead of the electric-anchor arm for the production witness",
  "LANDING_OWNERSHIP" ->
    "resolve the live electric anchor while retaining the emitted production carriage",
  "ACTIVE_FLUX_CAVEAT" ->
    "absorb F_flux into the conservative exchange and mark its value/integrability closed",
  "HOOK_LORENTZ" ->
    "close delta_BA, r_cone, higher orders, and the active-flux match",
  "TARGET_BLINDNESS" ->
    "install a nonsealed landing-token DownValue and inject every barred path marker",
  "DUAL_ENGINE_TERMS" ->
    "corrupt every local computed inventory lane",
  "UNITS_RESTORED" ->
    "corrupt every cited base dimension before deriving the hook dimensions",
  "VERDICT_REDERIVATION" ->
    "resolve the production witness anchor before running the native landing table",
  "SOURCE_TO_STAGE_MANIFEST" ->
    "remove a scoped-out source row from the live manifest"
|>;

publishedDigest =
  "983556935e50f12670fef24f17a23c048e295ddf0b6952aa0bd1618e9f179619";

raise[message_] := Throw[message, "ledgerStage038Failure"];

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
(* Neutral fact domain and ordinal mixed-radix enumeration.                *)
(* ---------------------------------------------------------------------- *)

currentDomain = {
  "CONVECTION_LIKE_CONDITIONAL",
  "CHARACTERIZED_SOURCE_DEPARTURE",
  "NULL_SOURCE",
  "R1_SOURCE_BASIS"
};
comparisonDomain = {
  "routes_agree",
  "routes_differ",
  "route_B_R1"
};
relativeDomain = {
  "relative_sign_match",
  "relative_sign_opposite",
  "leading_tensor_conflict",
  "relative_sign_anchor_conditional"
};
magnitudeDomain = {
  "magnitude_forced_by_electric",
  "magnitude_free_factor",
  "R1(magnitude)"
};
tierDomain = {
  "tier_A_gaps_closed",
  "tier_A_conditional"
};
anchorDomain = {
  "electric_anchor_closed",
  "R1_REQUIRED(bc_selection)"
};
inconsistencyDomain = {False, True};

domainAxes = {
  currentDomain,
  comparisonDomain,
  relativeDomain,
  magnitudeDomain,
  tierDomain,
  anchorDomain,
  inconsistencyDomain
};
domainKeys = {
  "current",
  "comparison",
  "relative_sign",
  "magnitude",
  "tier",
  "anchor",
  "internal"
};

decodeOrdinal[ordinal_Integer] := Module[
  {quotient = ordinal, digits, axis},
  digits = ConstantArray[1, Length[domainAxes]];
  Do[
    digits[[axis]] =
      Mod[quotient, Length[domainAxes[[axis]]]] + 1;
    quotient =
      Quotient[quotient, Length[domainAxes[[axis]]]],
    {axis, Length[domainAxes], 1, -1}
  ];
  MapThread[Part, {domainAxes, digits}]
];


(* ---------------------------------------------------------------------- *)
(* SEALED section 4: declarative adjudicator and independent oracle.       *)
(* ---------------------------------------------------------------------- *)

section4Adjudicate[case_Association, precedenceMutation_: False] := Module[
  {complete, precedenceTable},
  complete =
    case["current"] =!= "R1_SOURCE_BASIS" &&
    case["comparison"] =!= "route_B_R1" &&
    case["relative_sign"] =!= "relative_sign_anchor_conditional" &&
    case["tier"] === "tier_A_gaps_closed" &&
    case["anchor"] === "electric_anchor_closed" &&
    case["magnitude"] =!= "R1(magnitude)";

  precedenceTable = {
    {TrueQ[case["internal"]], "NO_GO(sector)"},
    {
      TrueQ[precedenceMutation] &&
        case["anchor"] === "R1_REQUIRED(bc_selection)",
      "R1_REQUIRED(consistency)"
    },
    {
      complete &&
        case["relative_sign"] === "relative_sign_opposite",
      "AMENDMENT_EXCLUDED(wrong_relative_sign)"
    },
    {
      complete &&
        case["relative_sign"] === "leading_tensor_conflict",
      "AMENDMENT_EXCLUDED(routes_leading_conflict)"
    },
    {
      complete &&
        case["comparison"] === "routes_agree" &&
        case["magnitude"] === "magnitude_forced_by_electric",
      "MAGNETISM_LORENTZ_CONSISTENT"
    },
    {
      complete && case["comparison"] === "routes_agree",
      "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)"
    },
    {
      complete && case["comparison"] === "routes_differ",
      "MAGNETISM_DEPARTURE_CHARACTERIZED"
    },
    {
      case["anchor"] === "R1_REQUIRED(bc_selection)",
      "R1_REQUIRED(electric_bc_selection)"
    },
    {
      case["current"] === "R1_SOURCE_BASIS" ||
        case["comparison"] === "route_B_R1",
      "R1_REQUIRED(direct_moving_throat)"
    },
    {
      case["magnitude"] === "R1(magnitude)",
      "R1_REQUIRED(magnitude)"
    },
    {
      case["tier"] === "tier_A_conditional" ||
        case["relative_sign"] ===
          "relative_sign_anchor_conditional",
      "R1_REQUIRED(consistency)"
    },
    {True, "R1_REQUIRED(unclassified)"}
  };

  FirstCase[
    precedenceTable,
    {True, landing_} :> landing
  ]
];

section4Oracle[case_Association] := Module[
  {
    completeFlags,
    complete,
    resolvedKey,
    resolvedMap,
    openFlags,
    openIndex,
    openMap
  },
  If[TrueQ[case["internal"]], Return["NO_GO(sector)"]];

  completeFlags = {
    case["current"] =!= "R1_SOURCE_BASIS",
    case["comparison"] =!= "route_B_R1",
    case["relative_sign"] =!=
      "relative_sign_anchor_conditional",
    case["tier"] === "tier_A_gaps_closed",
    case["anchor"] === "electric_anchor_closed",
    case["magnitude"] =!= "R1(magnitude)"
  };
  complete = And @@ completeFlags;

  If[
    complete,
    resolvedKey = StringRiffle[
      {
        case["relative_sign"],
        case["comparison"],
        case["magnitude"]
      },
      "|"
    ];
    resolvedMap = <|
      "relative_sign_opposite|routes_agree|magnitude_forced_by_electric" ->
        "AMENDMENT_EXCLUDED(wrong_relative_sign)",
      "relative_sign_opposite|routes_agree|magnitude_free_factor" ->
        "AMENDMENT_EXCLUDED(wrong_relative_sign)",
      "relative_sign_opposite|routes_differ|magnitude_forced_by_electric" ->
        "AMENDMENT_EXCLUDED(wrong_relative_sign)",
      "relative_sign_opposite|routes_differ|magnitude_free_factor" ->
        "AMENDMENT_EXCLUDED(wrong_relative_sign)",
      "leading_tensor_conflict|routes_agree|magnitude_forced_by_electric" ->
        "AMENDMENT_EXCLUDED(routes_leading_conflict)",
      "leading_tensor_conflict|routes_agree|magnitude_free_factor" ->
        "AMENDMENT_EXCLUDED(routes_leading_conflict)",
      "leading_tensor_conflict|routes_differ|magnitude_forced_by_electric" ->
        "AMENDMENT_EXCLUDED(routes_leading_conflict)",
      "leading_tensor_conflict|routes_differ|magnitude_free_factor" ->
        "AMENDMENT_EXCLUDED(routes_leading_conflict)",
      "relative_sign_match|routes_agree|magnitude_forced_by_electric" ->
        "MAGNETISM_LORENTZ_CONSISTENT",
      "relative_sign_match|routes_agree|magnitude_free_factor" ->
        "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)",
      "relative_sign_match|routes_differ|magnitude_forced_by_electric" ->
        "MAGNETISM_DEPARTURE_CHARACTERIZED",
      "relative_sign_match|routes_differ|magnitude_free_factor" ->
        "MAGNETISM_DEPARTURE_CHARACTERIZED"
    |>;
    Return[Lookup[resolvedMap, resolvedKey]]
  ];

  openFlags = {
    case["anchor"] === "R1_REQUIRED(bc_selection)",
    case["current"] === "R1_SOURCE_BASIS" ||
      case["comparison"] === "route_B_R1",
    case["magnitude"] === "R1(magnitude)",
    case["tier"] === "tier_A_conditional" ||
      case["relative_sign"] ===
        "relative_sign_anchor_conditional"
  };
  openIndex = FirstPosition[openFlags, True, Missing["NoOpenReason"]];
  openMap = <|
    1 -> "R1_REQUIRED(electric_bc_selection)",
    2 -> "R1_REQUIRED(direct_moving_throat)",
    3 -> "R1_REQUIRED(magnitude)",
    4 -> "R1_REQUIRED(consistency)"
  |>;
  If[
    MissingQ[openIndex],
    "R1_REQUIRED(unclassified)",
    Lookup[openMap, First[openIndex]]
  ]
];

truthTable[totalityMutation_: False] := Module[
  {
    expectedTotal,
    valuesByOrdinal,
    records,
    case,
    got,
    want,
    row,
    rows,
    counts,
    expectedCounts,
    byteChunks,
    streamBytes,
    digest,
    unclassified,
    namedTargets
  },
  expectedTotal = Times @@ (Length /@ domainAxes);
  valuesByOrdinal =
    decodeOrdinal /@ Range[0, expectedTotal - 1];

  records = Reap[
    Do[
      case = AssociationThread[domainKeys, valuesByOrdinal[[ordinal + 1]]];
      got = section4Adjudicate[case];
      If[TrueQ[totalityMutation] && ordinal === 0,
        got = "R1_REQUIRED(unclassified)"
      ];
      want = section4Oracle[case];
      row = StringRiffle[
        (
          If[StringQ[#], #, ToString[#, InputForm]]
        ) & /@ valuesByOrdinal[[ordinal + 1]],
        "|"
      ] <> "|" <> got;
      If[
        !(TrueQ[totalityMutation] && ordinal === expectedTotal - 1),
        Sow[{row, got, want}]
      ],
      {ordinal, 0, expectedTotal - 1}
    ]
  ][[2, 1]];

  rows = records[[All, 1]];
  counts = KeySort@Counts[records[[All, 2]]];
  byteChunks = ToCharacterCode[#, "UTF-8"] & /@ rows;
  streamBytes = Flatten[Riffle[byteChunks, {10}]];
  digest = Hash[
    ByteArray[streamBytes],
    "SHA256",
    "HexString"
  ];
  unclassified = "R1_REQUIRED(unclassified)";

  expectedCounts = KeySort@<|
    "NO_GO(sector)" -> 576,
    "MAGNETISM_LORENTZ_CONSISTENT" -> 3,
    "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)" -> 3,
    "MAGNETISM_DEPARTURE_CHARACTERIZED" -> 6,
    "AMENDMENT_EXCLUDED(wrong_relative_sign)" -> 12,
    "AMENDMENT_EXCLUDED(routes_leading_conflict)" -> 12,
    "R1_REQUIRED(electric_bc_selection)" -> 288,
    "R1_REQUIRED(direct_moving_throat)" -> 144,
    "R1_REQUIRED(magnitude)" -> 48,
    "R1_REQUIRED(consistency)" -> 60
  |>;
  namedTargets = <|
    "internal" -> "NO_GO(sector)",
    "consistent_free" ->
      "MAGNETISM_CONSISTENT_FREE_MAGNITUDE(R1)",
    "lorentz_consistent" -> "MAGNETISM_LORENTZ_CONSISTENT"
  |>;

  <|
    "agreement" ->
      And @@ MapThread[SameQ, {records[[All, 2]], records[[All, 3]]}],
    "noUnclassified" -> Lookup[counts, unclassified, 0] === 0,
    "cellCount" -> Length[records],
    "expectedTotal" -> expectedTotal,
    "counts" -> counts,
    "expectedCounts" -> expectedCounts,
    "digest" -> digest,
    "namedTargets" -> namedTargets
  |>
];


(* ---------------------------------------------------------------------- *)
(* Cited production facts, blocker carriage, and honest hook objects.      *)
(* ---------------------------------------------------------------------- *)

liveFacts = <|
  "current" -> "CONVECTION_LIKE_CONDITIONAL",
  "comparison" -> "route_B_R1",
  "relative_sign" -> "relative_sign_anchor_conditional",
  "magnitude" -> "R1(magnitude)",
  "tier" -> "tier_A_conditional",
  "anchor" -> "R1_REQUIRED(bc_selection)",
  "internal_sectors" -> {}
|>;

caseFromFacts[facts_Association] := <|
  "current" -> facts["current"],
  "comparison" -> facts["comparison"],
  "relative_sign" -> facts["relative_sign"],
  "magnitude" -> facts["magnitude"],
  "tier" -> facts["tier"],
  "anchor" -> facts["anchor"],
  "internal" -> (facts["internal_sectors"] =!= {})
|>;

emittedLanding = section4Adjudicate[caseFromFacts[liveFacts]];

neutralTokens[facts_Association] := {
  NeutralToken[CurrentIdentity, facts["current"]],
  NeutralToken[ComparisonFact, facts["comparison"]],
  NeutralToken[RelativeSignFact, facts["relative_sign"]],
  NeutralToken[MagnitudeFact, facts["magnitude"]],
  NeutralToken[TierFact, facts["tier"]],
  NeutralToken[AnchorFact, facts["anchor"]],
  NeutralToken[InternalSectors, facts["internal_sectors"]]
};

carriedNeutralTokens = {
  NeutralToken[
    CurrentIdentity,
    "CONVECTION_LIKE_CONDITIONAL"
  ],
  NeutralToken[ComparisonFact, "route_B_R1"],
  NeutralToken[
    RelativeSignFact,
    "relative_sign_anchor_conditional"
  ],
  NeutralToken[MagnitudeFact, "R1(magnitude)"],
  NeutralToken[TierFact, "tier_A_conditional"],
  NeutralToken[AnchorFact, "R1_REQUIRED(bc_selection)"],
  NeutralToken[InternalSectors, {}]
};

expectedBlockers = {
  "R1_REQUIRED(electric_bc_selection)",
  "R1_REQUIRED(direct_moving_throat)",
  "R1_REQUIRED(magnitude)",
  "R1_REQUIRED(consistency)"
};

r1Blockers[facts_Association] := Cases[
  {
    {
      facts["anchor"] === "R1_REQUIRED(bc_selection)",
      "R1_REQUIRED(electric_bc_selection)"
    },
    {
      facts["current"] === "R1_SOURCE_BASIS" ||
        facts["comparison"] === "route_B_R1",
      "R1_REQUIRED(direct_moving_throat)"
    },
    {
      facts["magnitude"] === "R1(magnitude)",
      "R1_REQUIRED(magnitude)"
    },
    {
      facts["tier"] === "tier_A_conditional",
      "R1_REQUIRED(consistency)"
    }
  },
  {True, token_} :> token
];

landingOwnershipGuard[facts_Association] := {
  r1Blockers[facts] === expectedBlockers,
  carriedNeutralTokens === neutralTokens[facts],
  emittedLanding === section4Adjudicate[caseFromFacts[facts]]
};

conservativeExchange = coefficientB velocityPair;
deltaBA = Expand[rBA - 1];

activeFluxState[mutate_: False] := If[
  TrueQ[mutate],
  <|
    "fact" -> ActiveFluxAbsorbedDecidedZero,
    "integrability" -> IntegrabilityClosed,
    "conservative" -> conservativeExchange + Fflux
  |>,
  <|
    "fact" -> ActiveFluxR1UnmatchedRemainder,
    "integrability" -> IntegrabilityR1Open,
    "conservative" -> conservativeExchange
  |>
];

activeFluxGuard[state_Association] := {
  state["fact"] === ActiveFluxR1UnmatchedRemainder,
  state["integrability"] === IntegrabilityR1Open,
  FreeQ[state["conservative"], Fflux]
};

lorentzLockState[closeAll_: False] := Module[
  {
    substitutions,
    deltaZero,
    coneOne,
    locks,
    determination
  },
  substitutions = If[TrueQ[closeAll], {rBA -> 1, rCone -> 1}, {}];
  deltaZero = TrueQ[Simplify[deltaBA /. substitutions] === 0];
  coneOne = TrueQ[Simplify[(rCone - 1) /. substitutions] === 0];
  locks = {
    deltaZero,
    coneOne,
    TrueQ[closeAll],
    TrueQ[closeAll]
  };
  determination = If[
    And @@ locks,
    LorentzDetermined,
    LorentzUndetermined
  ];
  {locks, determination}
];


(* ---------------------------------------------------------------------- *)
(* Held-DownValues seal, local term inventory, and dimensional firewall.   *)
(* ---------------------------------------------------------------------- *)

targetBlindnessGuard[mutate_: False] := Module[
  {
    allowed,
    policed,
    symbols,
    functionName,
    strings,
    violations,
    citedSymbols,
    barred,
    barredIntersection
  },
  allowed = {
    "section4Adjudicate",
    "section4Oracle",
    "truthTable",
    "targetBlindnessGuard"
  };
  policed = {
    "magnetism_lorentz_consistent",
    "amendment_excluded",
    "magnetism_departure_characterized",
    "no_go(sector)"
  };
  If[
    TrueQ[mutate],
    nonsealedLeak[] := "MAGNETISM_LORENTZ_CONSISTENT"
  ];

  symbols = Select[
    Names["Global`*"],
    DownValues[#] =!= {} &
  ];
  violations = Reap[
    Do[
      functionName = Last[StringSplit[symbol, "`"]];
      strings = Cases[
        DownValues[symbol],
        item_String :> ToLowerCase[item],
        Infinity
      ];
      Do[
        If[
          StringContainsQ[string, token] &&
            !MemberQ[allowed, functionName],
          Sow[{functionName, token}]
        ],
        {string, strings},
        {token, policed}
      ],
      {symbol, symbols}
    ]
  ][[2]];
  violations = If[violations === {}, {}, Sort@DeleteDuplicates@First[violations]];

  citedSymbols = {
    "A_E",
    "q_T",
    "r_BA",
    "delta_BA",
    "r_cone",
    "c_gamma",
    "c_E",
    "mu_R",
    "rho_br",
    "F_flux"
  };
  barred = {"N_u", "a_T", "a'_T", "a_L", "q_A^T", "q_L"};
  If[TrueQ[mutate], citedSymbols = Union[citedSymbols, barred]];
  barredIntersection = Intersection[citedSymbols, barred];
  {violations, barredIntersection}
];

dualInventory[truth_Association, mutate_: False] := Module[
  {
    counts,
    digest,
    landing,
    blockers,
    delta,
    cone,
    hook,
    locks,
    determination
  },
  counts = truth["counts"];
  digest = truth["digest"];
  landing = emittedLanding;
  blockers = r1Blockers[liveFacts];
  delta = deltaBA;
  cone = rCone;
  hook = activeFluxState[];
  {locks, determination} = lorentzLockState[];

  If[
    TrueQ[mutate],
    counts = KeyDrop[counts, First@Keys[counts]];
    digest = Hash[
      ByteArray[ToCharacterCode[digest, "UTF-8"]],
      "SHA256",
      "HexString"
    ];
    landing = section4Adjudicate[
      caseFromFacts[
        Join[liveFacts, <|"anchor" -> "electric_anchor_closed"|>]
      ]
    ];
    blockers = Most[blockers];
    delta = Expand[deltaBA + 1];
    cone = 2 rCone;
    {locks, determination} = lorentzLockState[True];
    hook = activeFluxState[True]
  ];
  {
    counts,
    digest,
    landing,
    blockers,
    delta,
    cone,
    locks,
    determination,
    hook
  }
];

dimensionless = {0, 0, 0, 0};
dimAdd[left_List, right_List] := MapThread[Plus, {left, right}];
dimSubtract[left_List, right_List] := MapThread[Subtract, {left, right}];
dimScale[power_Integer, value_List] := power value;

unitState[corrupt_: False] := Module[
  {
    aeDim,
    qtDim,
    cg2Dim,
    ce2Dim,
    muDim,
    ratioDim,
    derivedDeltaDim,
    coneDim
  },
  If[
    TrueQ[corrupt],
    aeDim = {1, 0, 1, 1};
    qtDim = {0, 1, 0, 1};
    cg2Dim = {1, 1, -1, 0};
    ce2Dim = {0, 3, -1, 1};
    muDim = {0, 0, 0, 0},
    aeDim = {0, 1, 0, 1};
    qtDim = {1, 0, -1, 0};
    cg2Dim = {0, 2, -2, 0};
    ce2Dim = {0, 2, -2, 0};
    muDim = {2, 1, -4, -1}
  ];
  ratioDim = dimSubtract[
    dimAdd[dimScale[2, qtDim], cg2Dim],
    dimAdd[muDim, aeDim]
  ];
  derivedDeltaDim = If[
    ratioDim === dimensionless,
    dimensionless,
    Missing["Inhomogeneous"]
  ];
  coneDim = dimSubtract[ce2Dim, cg2Dim];
  {
    aeDim,
    qtDim,
    cg2Dim,
    ce2Dim,
    muDim,
    ratioDim,
    derivedDeltaDim,
    coneDim
  }
];

expectedUnitState = {
  {0, 1, 0, 1},
  {1, 0, -1, 0},
  {0, 2, -2, 0},
  {0, 2, -2, 0},
  {2, 1, -4, -1},
  dimensionless,
  dimensionless,
  dimensionless
};


(* ---------------------------------------------------------------------- *)
(* Complete source-to-stage predicate manifest.                            *)
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

inScopeSource = {
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
  "PARITY_ROTATION",
  "PARITY_TIME_REVERSAL",
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
  "QMAG_R1"
};

manifestDisposition[identifier_String] := Which[
  MemberQ[inScopeSource, identifier],
    {"REPLACED_BY_STRONGER", "STAGE038_SEALED_RECONSTRUCTION"},
  MemberQ[buildGlobalSource, identifier],
    {"REPLACED_BY_STRONGER", "STAGE038_BUILD_GLOBAL"},
  MemberQ[Take[scopedOutSource, 6], identifier],
    {"SCOPED_OUT", "STAGE034_V1_DONE"},
  MemberQ[Take[scopedOutSource, {7, 13}], identifier],
    {"SCOPED_OUT", "STAGE035_V2_DONE"},
  MemberQ[Take[scopedOutSource, {14, 16}], identifier],
    {"SCOPED_OUT", "STAGE036_V3_DONE"},
  MemberQ[Take[scopedOutSource, {17, 27}], identifier],
    {"SCOPED_OUT", "STAGE037_V4_DONE"},
  True,
    raise["MANIFEST_IDENTIFIER_UNACCOUNTED"]
];

sourceManifest = (
  Join[{#}, manifestDisposition[#]] &
) /@ sourceToothIDs;

manifestDigest[manifest_List] := Module[{stream, bytes},
  stream = StringRiffle[
    Sort[StringRiffle[#, "|"] & /@ manifest],
    "\n"
  ];
  bytes = ToCharacterCode[stream, "UTF-8"];
  Hash[ByteArray[bytes], "SHA256", "HexString"]
];

manifestState[manifest_List] := {
  manifest[[All, 1]],
  KeySort@Counts[manifest[[All, 2]]],
  Cases[
    manifest,
    {identifier_, disposition_, _} /;
      disposition =!= "SCOPED_OUT" :> identifier
  ],
  Cases[
    manifest,
    {identifier_, "SCOPED_OUT", _} :> identifier
  ],
  manifestDigest[manifest]
};

expectedManifestState = manifestState[sourceManifest];


(* ---------------------------------------------------------------------- *)
(* Build-native verdict table and executable assertions.                   *)
(* ---------------------------------------------------------------------- *)

verdictTable[truth_Association, mutate_: False] := Module[
  {
    production,
    internal,
    anchorResolved,
    consistentFree,
    lorentzConsistent,
    actual,
    expected
  },
  production = If[
    TrueQ[mutate],
    Join[liveFacts, <|"anchor" -> "electric_anchor_closed"|>],
    liveFacts
  ];
  internal = Join[
    liveFacts,
    <|"internal_sectors" -> {"injected"}|>
  ];
  anchorResolved = Join[
    liveFacts,
    <|"anchor" -> "electric_anchor_closed"|>
  ];
  consistentFree = <|
    "current" -> "CONVECTION_LIKE_CONDITIONAL",
    "comparison" -> "routes_agree",
    "relative_sign" -> "relative_sign_match",
    "magnitude" -> "magnitude_free_factor",
    "tier" -> "tier_A_gaps_closed",
    "anchor" -> "electric_anchor_closed",
    "internal_sectors" -> {}
  |>;
  lorentzConsistent = Join[
    consistentFree,
    <|"magnitude" -> "magnitude_forced_by_electric"|>
  ];

  actual = {
    section4Adjudicate[caseFromFacts[production]],
    section4Adjudicate[caseFromFacts[internal]],
    section4Adjudicate[caseFromFacts[anchorResolved]],
    section4Adjudicate[caseFromFacts[consistentFree]],
    section4Adjudicate[caseFromFacts[lorentzConsistent]],
    section4Adjudicate[caseFromFacts[liveFacts], True]
  };
  expected = {
    emittedLanding,
    truth["namedTargets"]["internal"],
    "R1_REQUIRED(direct_moving_throat)",
    truth["namedTargets"]["consistent_free"],
    truth["namedTargets"]["lorentz_consistent"],
    "R1_REQUIRED(consistency)"
  };
  {actual, expected}
];

ok = Catch[
  If[
    activeMutation =!= "" &&
      !MemberQ[toothOrder, activeMutation],
    Print["FIRST_FAILURE=UNKNOWN_MUTATION"];
    Print["FAIL  UNKNOWN_MUTATION: ", activeMutation];
    raise["UNKNOWN_MUTATION"]
  ];

  Print[
    "ledger_stage038_sealed_landing_electric_bc_r1 Mathematica audit"
  ];
  Print[
    "ROUTE=mixed-radix ordinal stream + declarative FirstCase precedence + key/open-reason oracle + UTF-8 ByteArray SHA256"
  ];
  Print[
    "INDEPENDENCE=no Tuples; no nested-If adjudicator; no source-file text scan"
  ];
  Print["FILE_IO=none; CROSS_ENGINE_COMPARE=none"];
  Print["PROGRESS MATHEMATICA STAGE038 ENUMERATION START"];
  If[
    activeMutation =!= "",
    Print["ACTIVE_MUTATION=", activeMutation];
    Print[
      "MUTATED_PRIMITIVE=",
      Lookup[ablationDescriptions, activeMutation]
    ]
  ];

  section["SEALED section-4 exhaustive truth table"];
  truth = truthTable[activeMutation === "TRUTH_TOTALITY"];
  truthComponents = {
    truth["agreement"],
    truth["noUnclassified"],
    truth["cellCount"] === truth["expectedTotal"],
    truth["digest"] === publishedDigest,
    truth["counts"] === truth["expectedCounts"]
  };
  expectBool[
    "TRUTH_TOTALITY",
    And @@ truthComponents,
    <|
      "Components" -> truthComponents,
      "Cells" -> truth["cellCount"],
      "ExpectedProduct" -> truth["expectedTotal"],
      "Digest" -> truth["digest"],
      "Counts" -> truth["counts"]
    |>
  ];
  Print[
    "      cells=", truth["cellCount"],
    "; expected_product=", truth["expectedTotal"],
    "; digest=", truth["digest"]
  ];
  Print["PROGRESS MATHEMATICA STAGE038 ENUMERATION COMPLETE"];

  precedenceLanding = section4Adjudicate[
    caseFromFacts[liveFacts],
    activeMutation === "TRUTH_PRECEDENCE"
  ];
  expectBool[
    "TRUTH_PRECEDENCE",
    precedenceLanding === emittedLanding,
    <|
      "Production" -> emittedLanding,
      "Rederived" -> precedenceLanding
    |>
  ];
  Print["      production_first_match=", precedenceLanding];

  ownershipFacts = If[
    activeMutation === "LANDING_OWNERSHIP",
    Join[liveFacts, <|"anchor" -> "electric_anchor_closed"|>],
    liveFacts
  ];
  ownership = landingOwnershipGuard[ownershipFacts];
  expectBool[
    "LANDING_OWNERSHIP",
    And @@ ownership,
    <|
      "Lanes" -> ownership,
      "Blockers" -> r1Blockers[ownershipFacts],
      "Carried" -> carriedNeutralTokens,
      "Live" -> neutralTokens[ownershipFacts],
      "Emitted" -> emittedLanding,
      "LiveLanding" ->
        section4Adjudicate[caseFromFacts[ownershipFacts]]
    |>
  ];
  blockers = r1Blockers[liveFacts];
  Print["      blockers=", StringRiffle[blockers, " | "]];

  section["Honest active-flux and Lorentz hooks"];
  fluxState = activeFluxState[
    activeMutation === "ACTIVE_FLUX_CAVEAT"
  ];
  fluxGuard = activeFluxGuard[fluxState];
  expectBool[
    "ACTIVE_FLUX_CAVEAT",
    And @@ fluxGuard,
    <|"Guard" -> fluxGuard, "State" -> fluxState|>
  ];
  Print[
    "      F_flux=separate R1 O(V1*V2) remainder; full integrability=R1"
  ];

  {locks, determination} = lorentzLockState[
    activeMutation === "HOOK_LORENTZ"
  ];
  lorentzGuard = {
    determination === LorentzUndetermined,
    locks === {False, False, False, False}
  };
  expectBool[
    "HOOK_LORENTZ",
    And @@ lorentzGuard,
    <|
      "deltaBA" -> deltaBA,
      "rCone" -> rCone,
      "Locks" -> locks,
      "Determination" -> determination
    |>
  ];
  Print[
    "      emergent_lorentz=UNDETERMINED; locks=(delta_BA,r_cone,higher_orders,active_flux)"
  ];

  section["Build-global structural, local-inventory, and unit firewalls"];
  blindness = targetBlindnessGuard[
    activeMutation === "TARGET_BLINDNESS"
  ];
  expectBool[
    "TARGET_BLINDNESS",
    blindness === {{}, {}},
    <|
      "TokenViolations" -> blindness[[1]],
      "BarredMarkers" -> blindness[[2]]
    |>
  ];
  Print[
    "      landing assignment lives only in SEALED section 4; the held-DownValues visitor restricts exactly four non-R1 substrings"
  ];

  referenceInventory = dualInventory[truth];
  liveInventory = dualInventory[
    truth,
    activeMutation === "DUAL_ENGINE_TERMS"
  ];
  expectBool[
    "DUAL_ENGINE_TERMS",
    liveInventory === referenceInventory,
    <|
      "CountTerms" -> Length[liveInventory[[1]]],
      "Digest" -> liveInventory[[2]],
      "Landing" -> liveInventory[[3]],
      "Blockers" -> liveInventory[[4]],
      "Delta" -> liveInventory[[5]],
      "Cone" -> liveInventory[[6]],
      "Locks" -> liveInventory[[7]]
    |>
  ];
  Print[
    "      local_terms=landing-count map,digest,landing,blockers,delta_BA,r_cone,hook conditions"
  ];

  dimensions = unitState[
    activeMutation === "UNITS_RESTORED"
  ];
  expectBool[
    "UNITS_RESTORED",
    dimensions === expectedUnitState,
    <|
      "LiveDimensions" -> dimensions,
      "Expected" -> expectedUnitState
    |>
  ];
  Print[
    "      [A_E]=E*L; [q_T]=M*T^-1; [r_BA]=[delta_BA]=[r_cone]=1; [c_gamma^2]=[c_E^2]=L^2*T^-2"
  ];

  section["Build-native verdict re-derivation table"];
  {actualVerdicts, expectedVerdicts} = verdictTable[
    truth,
    activeMutation === "VERDICT_REDERIVATION"
  ];
  Print[
    "      REDERIVED_LANDINGS=",
    StringRiffle[actualVerdicts, " | "]
  ];
  expectBool[
    "VERDICT_REDERIVATION",
    actualVerdicts === expectedVerdicts,
    <|
      "Actual" -> actualVerdicts,
      "Expected" -> expectedVerdicts
    |>
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
    liveManifestState === expectedManifestState,
    <|
      "Partition" -> liveManifestState[[2]],
      "InScope" -> liveManifestState[[3]],
      "ScopedOut" -> liveManifestState[[4]],
      "Digest" -> liveManifestState[[5]]
    |>
  ];
  sourceTotal =
    Length[inScopeSource] +
    Length[scopedOutSource] +
    Length[buildGlobalSource];
  Print[
    "      source_manifest=", sourceTotal,
    "; in_scope=", Length[inScopeSource],
    "; scoped_out=", Length[scopedOutSource],
    "; build_global=", Length[buildGlobalSource]
  ];
  Print[
    "      SCOPED_OUT=",
    StringRiffle[scopedOutSource, ","]
  ];

  Print[""];
  Print["LANDING=", emittedLanding];
  Print["BLOCKERS=", StringRiffle[blockers, " | "]];
  Print["TRUTH_CELLS=", truth["cellCount"]];
  Print["TRUTH_DIGEST=", truth["digest"]];
  KeyValueMap[
    Print["LANDING_COUNT[", #1, "]=", #2] &,
    truth["counts"]
  ];
  Print["DEFENSIVE_UNCLASSIFIED=0"];
  Print["HOOK_ACTIVE_FLUX=R1_UNMATCHED_O(V1*V2)"];
  Print["HOOK_LORENTZ=UNDETERMINED"];
  Print[
    "STAGE039_BOUNDARY=departure characterization is not this terminal verdict"
  ];
  Print["PROGRESS MATHEMATICA STAGE038 COMPLETE"];

  If[
    activeMutation =!= "",
    Print["FIRST_FAILURE=MUTATION_DID_NOT_FIRE"];
    raise["MUTATION_DID_NOT_FIRE"]
  ];
  True,
  "ledgerStage038Failure",
  Function[{message, tag}, False]
];

Print[""];
Print["TOOTH_COUNT=", Length[toothOrder]];
Print["PASS tally: ", passCount, "; FAIL tally: ", failCount];
If[
  TrueQ[ok],
  Print[
    "OVERALL PASS: Mathematica independently reached ",
    emittedLanding
  ];
  Exit[0],
  Print[
    "OVERALL FAIL: Mathematica stage038 audit did not close"
  ];
  Exit[1]
]
